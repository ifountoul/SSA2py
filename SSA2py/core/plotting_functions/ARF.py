#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright (C) 2022 Ioannis Fountoulakis, Christos Evangelidis

#    This file is part of SSA2py.

#    SSA2py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, 
#    or any later version.

#    SSA2py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with SSA2py.  If not, see <https://www.gnu.org/licenses/>.


import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Rectangle
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.geoaxes


from scipy.stats import binned_statistic_2d
from sklearn import preprocessing as p
from scipy.stats import circmean
from scipy.special import i0

from obspy import read_inventory
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import kilometer2degrees
from obspy.geodetics.base import degrees2kilometers
from obspy import read


def _A1inv(x):
    # Approximation for _A1inv(x) according R Package 'CircStats'
    # See http://www.scienceasia.org/2012.38.n1/scias38_118.pdf, equation (4)
    if 0 <= x < 0.53:
        return 2.0 * x + x * x * x + (5.0 * x**5) / 6.0
    elif x < 0.85:
        return -0.4 + 1.39 * x + 0.43 / (1.0 - x)
    else:
        return 1.0 / (x * x * x - 4.0 * x * x + 3.0 * x)

def ARFplots(data, inv, evla, evlo, evdepth, st, origin, tt,\
             filename='ARF', outpath='.', fileformat='png', dpi=400):
    """
    The ARF results will be presented in four plots.
    --> Brightness map with the lobe at the maximum brightness time.
    --> Maximum brightness in time. 
    --> Seismic array positions with the ARF direction.
    --> Vertical Cross-Sections
    --> Synthetics plot

    Arguments:
    ----------
    data: str
         ARF results path.
    inv: Obpsy inventory object
         Inventory with stations info.
    evla, evlo, evdepth: float
         Hypocenter position
    st: Obspy Stream Object
         Streams
    origin: UTCDatime Object
         Origin Time
    tt: numpy-array
         Travel-Times
    filename: str
         Name of the Figures
    folder: str
         Folder to save figures.
    fileformat: str
         Format to save figure. 
    dpi: int
         Dpi of the figure.

    Returns:
    --------
    Saved figure
 
    """

    data_path = os.path.split(data[0])[0]

    # Read maximum brightness
    m = np.load(os.path.join(data_path, 'out_Max.npy'), allow_pickle=True)
    xmax = m[:, -1]; ymax = m[:,0];
    _lon_ = m[:, 1]; _lat_ = m[:, 2]; _depth_ = kilometer2degrees(m[:, 3]);

    _dist_ = round(max(abs(max(_lon_)-min(_lon_)), abs(max(_lat_)-min(_lat_))),2)

    max_index = np.argmax(ymax)

    #Open the Figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

    # Plot 1.
    #########
    ax1 = plt.subplot(131, adjustable='box')

    # Get the Max. Brightness time
    d_ = np.load(data[np.argmax(ymax)])

    # Normalize data
    min_max_scaler = p.MinMaxScaler()
    d_[:,0] = min_max_scaler.fit_transform(d_[:,0][:,np.newaxis]).ravel()


    #Find the unique Lat, Lon
    lat = np.unique(d_[:,-3]); lon = np.unique(d_[:,-4]);

    ret = binned_statistic_2d(d_[:,-4], d_[:,-3], d_[:,0], statistic='max', bins=len(lon))

    extent = np.min(lon), np.max(lon), np.min(lat), np.max(lat)
    sc = ax1.imshow(ret.statistic.T, origin='lower', extent=extent,\
                    cmap='jet', interpolation='spline16', vmin=0, vmax=1, alpha=1.0)

    ax1.set_xlabel('Longitude ($^\circ$)', fontsize=12)
    ax1.set_ylabel('Latitude ($^\circ$)', fontsize=12)
    ax1.grid(which='both', alpha=0.5, ls='--')

    ax1.set_xlim([_lon_[max_index]-_dist_, _lon_[max_index]+_dist_])
    ax1.set_ylim([_lat_[max_index]-_dist_, _lat_[max_index]+_dist_])
    ax1.axvline(x=_lon_[max_index], color='k', ls='--', lw=1)
    ax1.axhline(y=_lat_[max_index], color='k', ls='--', lw=1)
    ax1.set_title('ARF shape (Time: ' + str(xmax[max_index]) + ' s)', fontweight='bold', fontsize=11)


    cbar = fig.colorbar(sc, ax=ax1, shrink=0.7)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Norm. Brightness', rotation=270)

    # Plot 2.
    #########
    ax2 = plt.subplot(132, projection=ccrs.PlateCarree(), adjustable='box')

    # Get the coastline
    GSHHS = cfeature.GSHHSFeature(scale='full')

    # Add coastline
    ax2.add_feature(GSHHS, linewidth=1.0, facecolor = 'none', zorder=1)
    ax2.set_extent([_lon_[max_index]-_dist_, _lon_[max_index]+_dist_, _lat_[max_index]-_dist_, _lat_[max_index]+_dist_], crs=ccrs.PlateCarree())

    #Colorbar
    cm = plt.cm.get_cmap('jet')

    # Time
    ymaxnorm = min_max_scaler.fit_transform(ymax[:,np.newaxis]).ravel()

    # Plot scatter
    sc = ax2.scatter(m[:,1], m[:,2], s=(ymaxnorm+1.5)**6, c=xmax, cmap=cm, edgecolor='black',\
                     linewidth=1,  alpha=0.8, transform=ccrs.PlateCarree())

    ax2.set_title('ARF Maximum Brightness Positions', fontweight='bold', fontsize=11)

    cbar = fig.colorbar(sc, ax=ax2, shrink=0.7)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Time (s)', rotation=270)

    gl = ax2.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', alpha=0.2, linestyle='-', draw_labels=True)
    gl.top_labels = False
    gl.left_labels = True
    gl.right_labels = False
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator(ax1.get_xticks())
    gl.ylocator = mticker.FixedLocator(ax1.get_yticks())
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    braxes = inset_axes(ax2, '50%', '20%', loc = 'lower right')
    braxes.plot(xmax, ymaxnorm, color='r', lw=2)
    braxes.set_xlabel('Time (s)', fontweight='bold')
    braxes.set_ylabel('Max. Norm. Br.', rotation=270, labelpad=10, fontweight='bold')
    braxes.yaxis.set_label_position("left")
    braxes.xaxis.set_label_position("top")
    braxes.yaxis.set_ticks_position("left")
    braxes.xaxis.set_ticks_position("top")
    braxes.autoscale(axis='x', tight=True)

    # Plot 3
    # New sub-plot
    ax3 = plt.subplot(133, projection='polar')

    # Position North
    ax3.set_theta_zero_location('N')
    # Measure clockwise
    ax3.set_theta_direction(-1)

    # Fix limits
    ax3.set_ylim([0,10])
    # Fix grid
    ax3.grid(color='gray', ls='-', linewidth=1, alpha=0.2, zorder=0)

    # Add a star to the middle
    hypo = ax3.scatter(0,-0.5,400, marker='*', color='red', edgecolor='k', zorder=5, label = 'Hypocenter')

    # Add stations in the circumference
    azim = []
    for net in inv:
        for sta in net:
             azim.append(np.radians(gps2dist_azimuth(evla, evlo, sta.latitude, sta.longitude)[1]))
    azim = np.array(azim)

    sta = ax3.scatter(azim, np.repeat(8.9, azim.shape[0]), 130, marker = "^", color='springgreen',\
                      edgecolor='k', clip_on=False, zorder=5, alpha=0.8, label='Stations')

    # Connect with lines
    for an in azim:
        ax3.plot([0,an], [-0.5, 8.9], ls='-', color='k', lw=1, alpha=0.2)

    # Add the histogram
    # radius of inner circle
    bottom = 0
    nbins = 16
    # theta grid
    theta = np.linspace(0.0, 2*np.pi, nbins, endpoint=False)
    # histogram bin heights and edges
    bin_heights, bin_edges = np.histogram(azim, bins=nbins, range=(0,2*np.pi), density=False)
    # height of radial bins
    min_max_scaler = p.MinMaxScaler(feature_range=(0, 6))
    radii = min_max_scaler.fit_transform(bin_heights[:,np.newaxis]).ravel()
    # bin width
    width = bin_edges[1] - bin_edges[0]

    # bar plot on polar axes
    bars = ax3.bar(theta, radii, width=width, bottom=bottom, color='b', edgecolor='k',alpha=0.2, label='Histogram')

    # Von Mises Distribution
    mu = circmean(azim)
    kappa = _A1inv(np.mean(np.cos(azim - mu)))

    # Plot the Distribution
    x = np.linspace(0, 2*np.pi, num=501)
    y = np.exp(kappa*np.cos(x-mu))/(2*np.pi*i0(kappa))
    radii = min_max_scaler.fit_transform(y[:,np.newaxis]).ravel()

    VM = ax3.plot(x, radii, linewidth=2, color='firebrick', zorder=3, label='VM Distr.' + '\n' + 'Mu=' + str(round(np.degrees(mu),1)) + '$^\circ$' + 'Kappa=' + str(round(kappa,1)))

    # Add ticklabels
    ax3.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    ax3.set_yticklabels([])
    ax3.set_ylim(-0.5, 9);

    ax3.legend(fontsize='small')
    ax3.set_title('Stations Distribution', fontweight='bold', fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, filename+'_1'+'.'+fileformat), dpi=dpi)

    #########################################################
    # Plot ARF shape in-depth
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

    # Plot 1
    depth = (np.unique(d_[:,-2]))/111.11

    ret = binned_statistic_2d(d_[:,-4], d_[:,-2]/111.11,\
                              d_[:,0], statistic='max', bins=len(depth))

    extent = np.min(lon), np.max(lon), np.min(depth), np.max(depth)
    sc = ax1.imshow(ret.statistic.T, origin='lower', extent=extent,\
                    cmap='jet', interpolation='spline16', vmin=0, vmax=1, alpha=1.0)

    ax1.axvline(x=_lon_[max_index], color='k', ls='--', lw=1)
    ax1.axhline(y=_depth_[max_index], color='k', ls='--', lw=1)

    ax1.set_title('Lon. Vertical Cross Section (ARF-Time: ' + str(xmax[max_index]) + ' s)', fontweight='bold', fontsize=11)
    ax1.set_xlabel('Longitude ($^\circ$)', fontsize=12)
    ax1.set_ylabel('Depth (km)', fontsize=12)
    ax1.set_xlim([_lon_[max_index]-_dist_, _lon_[max_index]+_dist_])
    ax1.grid(which='both', alpha=0.5, ls='--')
    ax1.invert_yaxis()
    ax1.set_yticks(ax1.get_yticks())
    ax1.set_yticklabels(np.round(ax1.get_yticks()*111.12))
    ax1.autoscale()

    # Plot 2
    ret = binned_statistic_2d(d_[:,-3], d_[:,-2]/111.11,\
                              d_[:,0], statistic='max', bins=len(depth))

    extent = np.min(lat), np.max(lat), np.min(depth), np.max(depth)
    sc = ax2.imshow(ret.statistic.T, origin='lower', extent=extent,\
                    cmap='jet', interpolation='spline16', vmin=0, vmax=1, alpha=1.0)
    ax2.axvline(x=_lat_[max_index], color='k', ls='--', lw=1)
    ax2.axhline(y=_depth_[max_index], color='k', ls='--', lw=1)

    ax2.set_title('Lat. Vertical Cross Section (ARF-Time: ' + str(xmax[max_index]) + ' s)', fontweight='bold', fontsize=11)
    ax2.set_xlabel('Latitude ($^\circ$)', fontsize=12)
    ax2.set_ylabel('Depth (km)', fontsize=12)
    ax2.set_xlim([_lat_[max_index]-_dist_, _lat_[max_index]+_dist_])
    ax2.grid(which='both', alpha=0.5, ls='--')
    ax2.invert_yaxis()
    ax2.set_yticks(ax2.get_yticks())
    ax2.set_yticklabels(np.round(ax2.get_yticks()*111.12))
    ax2.autoscale()

    # Traces plot
    st.normalize(global_max=True)
    st.sort(keys=['distance'])

    scale = 1
    max_offset = st.count()
    min_offset = 1

    scale = (max_offset - min_offset) * 1/(st.count() * scale)

    count = 0
    for tr in st:
        data = (tr.data * scale) + count
        time = np.arange(tr.stats.npts) * tr.stats.delta\
                          + (tr.stats.starttime - origin)
        ax3.plot(time, data, color ='k', lw = 1.0, alpha = 0.6)
        ax3.text(math.floor(np.min(tt))-3, count+0.1, tr.stats.station, zorder=1,\
                 fontsize=8, fontstyle='italic', fontweight='demi')
        count+=1

    ax3.set_title('Synthetic Pulses', fontweight='bold', fontsize=11)
    ax3.set_xlabel('Time relative to Origin (s)', fontsize=12)
    ax3.set_yticklabels([])
    ax3.set_yticks([])
    ax3.set_xlim([math.floor(np.min(tt))-3, math.ceil(np.max(tt))+2])
    ax3.grid(which='both', alpha=0.5, ls='--')
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, filename+'_2'+'.'+fileformat), dpi=dpi)

    return


