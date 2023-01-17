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

import pyproj, multiprocessing, warnings, shapefile
import cartopy, cartopy.mpl.geoaxes, math, os, re, scipy, itertools
from scipy.interpolate import griddata, RBFInterpolator, interp1d
import numpy as np

from obspy.geodetics import kilometer2degrees
from obspy.core.inventory.inventory import read_inventory
from obspy.core import read
from obspy.core.stream import Stream
from obspy.geodetics.base import gps2dist_azimuth, degrees2kilometers, kilometer2degrees

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, BoundaryNorm, LightSource
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.geoaxes
from shapely.geometry.polygon import LinearRing
from shapely.ops import transform
from shapely.geometry import Point

from SSA2py.core import config
from SSA2py.core.basic_f.other import createDir, colorBarFix, FixPointNormalize

def NormalizeData(data):
    """
    Normalize data between 0-1.

    Arguments:
    ----------
    data: array-like
   
    Returns:
    -------
    Normalized array-like 

    From https://www.stackvidhya.com/how-to-normalize-data-between-0-and-1-range/
    """
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def compositeStack(data, evlo, evla, evdepth, hypo=True, timeStartStack=None, timeEndStack=None,\
                   minDepthStack=None, maxDepthStack=None, minLatStack=None, maxLatStack=None,\
                   minLonStack=None, maxLonStack=None, extentMinLat=None, extentMaxLat=None,\
                   extentMinLon=None, extentMaxLon=None, extentMinDepth=None, extentMaxDepth=None,\
                   maxgrid=100, colormap='jet', meridians=True, plotCross=True, interpolate=True,\
                   filename='CompositeStack', outpath='.', fileformat='pdf', dpi=400):

    """
    Composite stack of all time steps
 
    Arguments:
    ----------

    data: list
        List with .npy
    evlo: float
        Event Longitude
    evla: float
        Event Latitude
    evdepth: float
        Event Depth
    hypo: bool
        Plot hypo?
    timeStartStack: float
        Time to start stack.
    timeEndStack: float
        Time to end stack.
    minDepthStack: float
        Minimum depth to stack in x-y plane.
    maxDepthStack: float
        Maximum depth to stack in x-y plane
    minLatStack: float
        Minimum Latitude to stack in x-z plane.
    maxLatStack: float
        Maximum Latitude to stack in x-z plane.
    minLonStack: float
        Minimum Longitude to stack in y-z plane.
    maxLonStack: float
        Maximum Longitude to stak in y-z plane.
    extentMinLat: float
        Extent minimum Latitude.
    extentMaxLat: float
        Extent maximum Latitude.
    extentMinLon: float
        Extent minimum Longitude.
    extentMaxLon: float
        Extent maximum Longitude.
    extentMinDepth: float
        Extent minimum Depth.
    extentMaxDepth: float
        Extent maximum Depth.
    maximumgrid: float
        Maximum distance between grid-hypo.
    colormap: str
        Choose colormap.
    meridian: bool
        Plot meridians.
    plotCross: bool
        Plot Cross section lines in x-y plot?
    interpolate: bool
        Interpolate the planes?
    filename: str
        Filename.
    outpath: str
        Path to save.
    fileformat: str
        Format of the file.
    dpi: int
        Dpi.

    Returns:
    --------

    """

    # but let's get the data path
    data_path = os.path.split(data[0])[0]

    # read the max
    m = np.load(os.path.join(data_path, 'out_Max.npy'), allow_pickle=True)

    if (timeStartStack is None) and (timeEndStack is None):
        for i in range(len(m[:,-1])):
            if m[i,-1]>=0:
                time_S = m[i,-1]
                break
        time_E = m[-1,-1]
    if (timeStartStack is None) and (timeEndStack is not None):
        for i in range(len(m[:,-1])):
            if m[i,-1]>=0:
                time_S = m[i,-1]
                break
        time_E = timeEndStack
    if (timeStartStack is not None) and (timeEndStack is None):
        time_S = timeStartStack
        time_E = m[-1,-1]
    if (timeStartStack is not None) and (timeEndStack is not None):
        time_S = timeStartStack
        time_E = timeEndStack

    # pre-allocate the array
    resStack = np.zeros((np.load(data[0]).shape[0], 1))

    # stack the timesteps of all the grid (multiply)
    for i in range(len(data)):
       res = np.load(data[i])
       res[np.isnan(res[:,0]),0] = 0 # replace all nan with 1

       if  (res[0,-1]>= time_S) and (res[0,-1]<= time_E):
           resStack[:,0] = res[:,0]+resStack[:,0]

    # add the coordinates
    resStack = np.column_stack((resStack[:,0], res[:,1], res[:,2], res[:,3]))

    # now we have the time-stack grids move on to slices
    # 1) horizontal stack

    # get the unique lat-lon values
    lat = np.unique(resStack[:,2]); lon = np.unique(resStack[:,1]);

    #Get combinations
    x, y = np.meshgrid(lon, lat, sparse=True)

    # pre-allocate the array
    horStack = np.zeros((len(x.T), len(y)))

    if (minDepthStack is None) and (maxDepthStack is None):
        minDepthStack = np.min(resStack[:,3])
        maxDepthStack = np.max(resStack[:,3])
    if (minDepthStack is None) and (maxDepthStack is not None):
        minDepthStack = np.min(resStack[:,3])
    if (minDepthStack is not None) and (maxDepthStack is None):
        maxDepthStack = np.max(resStack[:,3])
    if (minDepthStack is not None) and (maxDepthStack is not None):
        pass

    for i in range(len(x.T)):
        for j in range(len(y)):
            pos = np.where(((resStack[:,1]==x.T[i]) & (resStack[:,2]==y[j]) &\
                  (resStack[:,3]<=maxDepthStack) & (resStack[:,3]>=minDepthStack)))[0]

            horStack[i, j] = np.sum(resStack[pos,0])

    #####################################

    horStack = NormalizeData(horStack)

    #Open the Figure
    fig = plt.figure(constrained_layout=False, figsize=(15,7))

    gs=GridSpec(7,15)

    #Plot
    ax1 = fig.add_subplot(gs[0:7,0:7], projection=ccrs.PlateCarree(), aspect="auto")

    GSHHS = cfeature.GSHHSFeature(scale='full')

    x, y = np.meshgrid(lon, lat, sparse=False)

    if interpolate==True:
        c = ax1.pcolormesh(x, y, horStack.T, cmap = colormap, shading='gouraud')
    else:
        c = ax1.pcolormesh(x, y, horStack.T, cmap = colormap)

    if (extentMinLat is not None) and (extentMaxLat is not None)\
        and (extentMinLon is not None) and (extentMaxLon is not None):
        ax1.set_extent([extentMinLon, extentMaxLon, extentMinLat, extentMaxLat])
    else:
        ax1.set_extent([evlo-kilometer2degrees(maxgrid)/2, evlo+kilometer2degrees(maxgrid)/2,\
                        evla-kilometer2degrees(maxgrid)/2, evla+kilometer2degrees(maxgrid)/2], crs=ccrs.PlateCarree())

    ax1.add_feature(GSHHS, linewidth=1.0, facecolor = 'none', zorder=1)

    if meridians==True:
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', alpha=0.2, linestyle='-', draw_labels=True)
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels =True
        gl.xlines = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

    if hypo==True:
        ax1.plot(evlo, evla, '*', color='k', linewidth=5, markersize=25, markeredgecolor='k',\
                 markeredgewidth=1.0, alpha=0.7, transform=ccrs.Geodetic())

    # add also the max brigthness in time and point out the time range we are working

    # create white box in the left

    box = inset_axes(ax1, '40%', '18%', loc = 'lower left')
    box.set_facecolor([1,1,1,0.7])
    box.set_xticks([])
    box.set_yticks([])

    braxes = inset_axes(box, '80%', '55%', loc = 'upper center')
    braxes.plot(m[:,-1], m[:,0], color='black')

    braxes.axvline(x=time_S, color = 'k', ls='--', linewidth=1, alpha=0.5)
    braxes.axvline(x=time_E, color = 'k', ls='--', linewidth=1, alpha=0.5)

    # fill the curve
    # take the slice
    tempSlice = m[(m[:,-1]>=time_S)&(m[:,-1]<=time_E), 0]

    braxes.fill_between(np.linspace(time_S,time_E,len(tempSlice)),\
                        np.zeros((1,len(tempSlice))).ravel(),\
                        tempSlice, alpha=0.5)

    braxes.set_ylabel('Max. Br.', fontsize=8, fontweight='bold')
    braxes.set_xlabel('Time (s)', fontsize=8, fontweight='bold')
    braxes.spines['right'].set_visible(False)
    braxes.spines['top'].set_visible(False)
    braxes.yaxis.set_ticks_position('none')
    braxes.set_yticks([])
    braxes.xaxis.set_ticks_position('bottom')
    braxes.patch.set_alpha(0.5)
    braxes.tick_params(axis='x', labelsize= 7)
    braxes.tick_params(axis='y', labelsize= 7)
    braxes.set_xlim([np.min(m[:,-1]), np.max(m[:,-1])])
    braxes.set_ylim([np.min(m[:,0]), np.max(m[:,0])])

    # North-south cross section
    ax2 = fig.add_subplot(gs[0:3,8:15])

    # get the unique lon-depth
    depth = np.unique(resStack[:,3])
    lon = np.unique(resStack[:,1])

    # get combinations
    x, y = np.meshgrid(lon, depth, sparse=True)

    #Choose depths
    if depth.max()<=40:
         step_ = 5
    elif depth.max()>=40 and depth.max()<=100:
         step_ = 10
    else:
         step_ = 100

    # pre-allocate the vertical stack
    verStack = np.zeros(((len(x.T),len(y))))

    if (minLatStack is None) and (maxLatStack is None):
         minLatStack = np.min(resStack[:,2])
         maxLatStack = np.max(resStack[:,2])
    if (minLatStack is None) and (maxLatStack is not None):
         minLatStack = np.min(resStack[:,2])
    if (minLatStack is not None) and (maxLatStack is None):
         maxLatStack = np.max(resStack[:,2])
    if (minLatStack is not None) and (maxLatStack is not None):
        pass

    if plotCross==True:
        ax1.axhline(y = minLatStack,\
                    color = 'k', ls='--', linewidth=2, alpha=0.5)
        ax1.axhline(y = maxLatStack,\
                    color = 'k', ls='--', linewidth=2, alpha=0.5)

    for x_ in range(len(x.T)):
        for y_ in range(len(y)):
            pos = np.where(((resStack[:,1]==x.T[x_]) & (resStack[:,3]==y[y_]) &\
                          (resStack[:,2]<=maxLatStack) & (resStack[:,2]>=minLatStack)))[0]
            verStack[x_, y_] = np.sum(resStack[pos,0])

    # normalize
    verStack = NormalizeData(verStack)

    x, y = np.meshgrid(lon, depth)

    y = y/111.1

    if interpolate==True:
        sc = ax2.pcolormesh(x, y, verStack.T, cmap =colormap, shading='gouraud')
    else:
        sc = ax2.pcolormesh(x, y, verStack.T, cmap =colormap)

    if hypo==True:
        ax2.plot(evlo, evdepth/111.1, '*', color='k', linewidth=5, markersize=20, markeredgecolor='k',\
                 markeredgewidth=1.0, alpha=0.5)

    ax2.invert_yaxis()
    ax2.set_yticks(np.arange(y.min(), y.max(), round(step_/111.1,4)))
    ax2.set_yticklabels(np.round(np.arange(y.min(), y.max(), round(step_/111.1,4)) * 111.1))
    ax2.yaxis.set_ticks_position('right')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.tick_params(axis='x', labelsize= 10)
    ax2.tick_params(axis='y', labelsize= 10)
    ax2.set_xlabel('Longitude ($^\circ$) (E-W)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Depth (km)', fontsize=10, fontweight='bold')

    if (extentMinLon is not None) and (extentMaxLon is not None):
        ax2.set_xlim([extentMinLon, extentMaxLon])
    if (extentMinDepth is not None) and (extentMaxDepth is not None):
        ax2.set_ylim([extentMinDepth, extentMaxDepth])

    # add here the colorbar
    box = inset_axes(ax2, '33%', '27%', loc = 'lower left')
    box.set_facecolor([1,1,1,0.7])
    box.set_xticks([])
    box.set_yticks([])

    cbaxes = inset_axes(box, '80%', '20%', loc = 'lower center')

    cb = plt.colorbar(sc, cax=cbaxes ,orientation="horizontal",\
         ticks=np.linspace(np.min(verStack.T), np.max(verStack.T), 3, endpoint=True))
    cb.set_label(label='Normalized Brightness', size='small', weight='bold')
    cb.ax.tick_params(labelsize='medium')
    cbaxes.xaxis.set_ticks_position("top")
    cbaxes.xaxis.set_label_position("top")

    ax3 = fig.add_subplot(gs[4:7,8:15], aspect="auto")

    # get the unique lon-depth
    depth = np.unique(resStack[:,3])
    lat = np.unique(resStack[:,2])

    # get combinations
    x, y = np.meshgrid(lat, depth, sparse=True)

    # pre-allocate the vertical stack
    verStack = np.zeros(((len(x.T),len(y))))

    if (minLonStack is None) and (maxLonStack is None):
         minLonStack = np.min(resStack[:,1])
         maxLonStack = np.max(resStack[:,1])
    if (minLonStack is None) and (maxLonStack is not None):
         minLonStack = np.min(resStack[:,1])
    if (minLonStack is not None) and (maxLonStack is None):
         maxLonStack = np.max(resStack[:,1])
    if (minLonStack is not None) and (maxLonStack is not None):
        pass

    if plotCross==True:
        ax1.axvline(x = minLonStack,\
                    color = 'k', ls='--', linewidth=2, alpha=0.5)
        ax1.axvline(x = maxLonStack,\
                    color = 'k', ls='--', linewidth=2, alpha=0.5)

    for x_ in range(len(x.T)):
        for y_ in range(len(y)):
            pos = np.where(((resStack[:,2]==x.T[x_]) & (resStack[:,3]==y[y_]) &\
                          (resStack[:,1]<=maxLonStack) & (resStack[:,1]>=minLonStack)))[0]
            verStack[x_, y_] = np.sum(resStack[pos,0])

    # normalize
    verStack = NormalizeData(verStack)

    x, y = np.meshgrid(lat, depth)

    y = y/111.1

    if interpolate==True:
        sc = ax3.pcolormesh(x, y, verStack.T, cmap =colormap, shading='gouraud')
    else:
        sc = ax3.pcolormesh(x, y, verStack.T, cmap =colormap)

    if hypo==True:
        ax3.plot(evla, evdepth/111.1, '*', color='k', linewidth=5, markersize=20, markeredgecolor='k',\
                 markeredgewidth=1.0, alpha=0.5)
    ax3.invert_yaxis()
    ax3.set_yticks(np.arange(y.min(), y.max(), round(step_/111.1,4)))
    ax3.set_yticklabels(np.round(np.arange(y.min(), y.max(), round(step_/111.1,4)) * 111.1))
    ax3.yaxis.set_ticks_position('right')
    ax3.xaxis.set_ticks_position('bottom')
    ax3.tick_params(axis='x', labelsize= 10)
    ax3.tick_params(axis='y', labelsize= 10)
    ax3.set_xlabel('Latitude ($^\circ$) (N-S)', fontsize=10, fontweight='bold')
    ax3.set_ylabel('Depth (km)', fontsize=10, fontweight='bold')

    if (extentMinLat is not None) and (extentMaxLat is not None):
        ax3.set_ylim([extentMinLat, extentMaxLat])
    if (extentMinDepth is not None) and (extentMaxDepth is not None):
        ax3.set_xlim([extentMinDepth, extentMaxDepth])

    ax1.set_title('Horizontal Cross Section ({} - {} km depth)'.format(str(minDepthStack), str(maxDepthStack)))
    ax2.set_title('East - West Vertical Cross Section')
    ax3.set_title('North - South Vertical Cross Section')

    fig.suptitle('Composite Stacks (time range {} - {} s)'.format(str(time_S),\
                 str(time_E)), weight='bold', size=15)

    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)
    plt.close('all')

    return








