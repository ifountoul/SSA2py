#/usr/bin/env python3
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

from netCDF4 import Dataset


# local functions
from SSA2py.core import config
from SSA2py.core.basic_f.other import createDir, colorBarFix, FixPointNormalize
from SSA2py.core.plotting_functions.other import circular_hist


def MaxBrightTimeStep_(brpath, brpathboot, evla, evlo, evdepth, time, inv, stations_used,\
                       startTime=None, endTime=None, minBrig=None, maxBrig=None,\
                       min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                       points_size=10, maxgrid = 100, faults=True, grid=True, hypo=True, colormap='plasma', topo=True,\
                       meridian=True, info_box=True, Test='MAIN', autoselect=True,\
                       filename='MaximumBrightness', outpath='.', fileformat='pdf', dpi=400):

    """
    Plot Maximum Brightness Per Time step

    Arguments:
    ------
    brpath: str
        Path to brightness results.
    brpathboot: str
        Path to brightness results for bootstrap.
    evla: float
        Event latitude.
    evlo: float
        Event longitude.
    evdepth: float
        Event depth.
    time: UTCDatime object
        Event Origin Time
    inv: Obspy Object
        Obspy Inventory
    stations_used: list
        Stations used
    startTime: float
        Start time to plot brightness.
    endTime: float
        End time to plot brightness.
    minBrig: float
        Minimum brightness to plot.
    maxBrig: float
        Maximum brightness to plot.
    min_lon: float
        Minimum Longitude.
    min_lat: float
        Minimim Latitude. 
    max_lon: float
        Maximum Longitude.
    max_lat: float
        Maximum Latitude.
    min_depth: float
        Minimum depth to plot.
    max_depth: float
        Maximum depth to plot.
    points_size: float
        Power to size the brightness cirlces.
    maxgrid: float
        Maximum distance between hypo and grid.
    faults: bool
        Plot faults?
    grid: bool
        Plot grid?
    hypo: bool
        Plot hypo?
    colormap: str
        Colormap name, supported from matplotlib
    topo: bool
        Plot topo?
    meridian: bool
        Plot meridians
    info_box: bool
        Plot info box
    Test: str
        MAIN = Main SSA results
    autoselect: bool
        If True the program automaticaly choose the end of the plot (when energy<0).
        If False plots all the points no mater the brightness value.
    ARF= Array Response Function
        JACK = Jacknife Test
        BOOT = Bootstrap Test 
    error_type: str
        CI=Confidence Interval
        SE=Standard Error
        STD=Standard Deviation
    filename: str
        Name of the file.
    fileformat: str
        File format such as .png, .pdf.
    outpath: str
        Path to save file.
    dpi: float
        The resolution in dots per inch.

    """
 
    # open figure
    fig = plt.figure(constrained_layout=False, figsize=(11,12))

    # break to parts
    gs = fig.add_gridspec(ncols=100, nrows=100)

    # map
    ax1 = fig.add_subplot(gs[0:65, 0:65], projection=ccrs.PlateCarree())

    # do you have info from user about the extent of the map?
    if (min_lon is None) and (min_lat is None) and (max_lon is None) and (max_lat is None):

        min_lon = evlo-kilometer2degrees(maxgrid)/2 - 0.08
        min_lat = evla-kilometer2degrees(maxgrid)/2 - 0.08
        max_lon = evlo+kilometer2degrees(maxgrid)/2 + 0.08
        max_lat = evla+kilometer2degrees(maxgrid)/2 + 0.08

    #Plot topography
    if config.cfg['Plotting']['Topography/Bathymetry'][0]==True and topo==True:

        #Get the data from dataset (This is tested for GEBCO)
        dset=Dataset(config.cfg['Plotting']['Topography/Bathymetry'][1])
        elev=dset.variables['elevation'][:].squeeze().data
        x=dset.variables['lon'][:].squeeze().data
        y=dset.variables['lat'][:].squeeze().data

        #Meshgrid
        [x,y] =np.meshgrid(x, y)
        #Flat the arrays
        x = x.flatten(); y = y.flatten(); elev = elev.flatten();
        #Filter the nc file
        indexX = np.where((x>= min_lon) & (x<= max_lon))[0]
        indexY = np.where((y>= min_lat) & (y<= max_lat))[0]
        #find common points
        horCommon = np.intersect1d(indexX, indexY)
        x = x[horCommon]; y = y[horCommon]; elev = elev[horCommon]

        ######Interpolate
        pts = 10000000
        [x_,y_] = np.meshgrid(np.linspace(np.min(x),np.max(x),int(np.sqrt(pts))),\
                            np.linspace(np.min(y),np.max(y), int(np.sqrt(pts))))
        z = griddata((x, y), elev, (x_, y_), method='cubic')
        ################# 
        # Hillshade
        ls = LightSource(azdeg=315, altdeg=45)

        dx, dy = abs(x_[0,0]-x_[0,1]), abs(y_[0,0]-y_[1,0])
        dx = dx*111.11
        dy = dy*111.11

        rgb = ls.shade(z, cmap=plt.get_cmap('Greys'), blend_mode='soft', vert_exag=1.0, dx=dx, dy=dy)
        #
        ax1.set_extent([min_lon, max_lon,\
                        min_lat, max_lat], crs=ccrs.PlateCarree())

        ax1.imshow(rgb, origin='lower', extent=[min_lon,\
                                                max_lon,\
                                                min_lat,\
                                                max_lat],\
                                                alpha=1, cmap = 'Greys', transform= ccrs.PlateCarree())
        #Fix colorbar
        colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 256))
        colors_land = plt.cm.gray(np.linspace(0.25, 1, 256))
        all_colors = np.vstack((colors_undersea, colors_land))
        terrain_map = colors.LinearSegmentedColormap.from_list(
                     'terrain_map', all_colors)
        divnorm = colors.TwoSlopeNorm(vmin=-7000., vcenter=0, vmax=7000)

        ax1.imshow(z, origin='lower', extent=[min_lon,\
                                              max_lon,\
                                              min_lat,\
                                              max_lat],\
                                              alpha=0.7, norm=divnorm, cmap = terrain_map, transform= ccrs.PlateCarree())

    else:
        ax1.set_extent([min_lon, max_lon,\
                        min_lat, max_lat], crs=ccrs.PlateCarree())

        GSHHS = cfeature.GSHHSFeature(scale='high')
        ax1.add_feature(GSHHS, linewidth=0.4, facecolor = 'lemonchiffon', zorder=1)

    
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='--', linewidth=0.5, color = 'peru', zorder=0)
    ax1.add_feature(cartopy.feature.LAKES, alpha=0.5, zorder=0)
    ax1.add_feature(cartopy.feature.RIVERS, zorder=0)

    if faults==True:
        try:
            faults = shapefile.Reader(os.path.join(config.cfg['Plotting']['Save Layers'],'FAULTS',\
                     'shapefile','gem_active_faults.shp'), encoding='latin-1')
            for tmp in faults.shapeRecords():
                 x = [i[0] for i in tmp.shape.points[:]]
                 y = [i[1] for i in tmp.shape.points[:]]
                 ax1.plot(x, y, color='k', lw=1.0, transform= ccrs.PlateCarree())
        except Exception as e:
            config.logger.warning(e)
            config.logger.warning('Faults Shapefile error')

    # plot grid
    # from the brightness path get the grid

    grid_ = np.load(os.path.join(brpath, 'grid.npy'))

    #Keep unique values in lat, lon
    lon = np.unique(np.array(grid_[:,0], dtype=float))
    lat = np.unique(np.array(grid_[:,1], dtype=float))
    depth = np.unique(grid_[:,2])

    if Test=='MAIN' or Test=='ARF':
        if grid==True:
            # plot the backround grid points
            if config.gridRules[0][0]=='box':
                lon,lat = np.meshgrid(lon, lat)
            else:
                lon = np.array(grid_[:,0])
                lat = np.array(grid_[:,1])

            if config.cfg['Plotting']['Topography/Bathymetry'][0]==True:
                ax1.scatter(lon, lat, color = 'lime', marker = 'o', s=2, alpha=0.1, transform= ccrs.PlateCarree()) 
            else:
                ax1.scatter(lon, lat, color = 'gray', marker = 'o', s=2, alpha=0.2, transform= ccrs.PlateCarree())

    # plot the brightness results
    if Test=='MAIN' or Test=='ARF':
        BR = np.load(os.path.join(brpath, 'out_Max.npy')) #1st path

        #In case of Inf replace with 0 
        BR[:,0][BR[:,0] == -np.inf] = 0

    if Test=='MAIN' or Test=='ARF':
        
        indexStart = 0; indexEnd = 0;

        # get brightness from time
        for i in range(len(BR)):
            if (startTime is None):
                if BR[i,-1]>=0:
                    indexStart = i
                    break
            else:
                if BR[i,-1]>=startTime:
                    indexStart = i
                    break

        for i in range(len(BR)):
            if (endTime is None):
                if BR[i,-1]==BR[-1,-1]:
                    indexEnd = i
                    break
            else:
                if BR[i,-1]>=endTime:
                    indexEnd = i
                    break
                else:
                    indexEnd = i

        if (minBrig is not None) and (maxBrig is not None):
            indicesBrig = np.where(((BR[:,0]>=minBrig) & (BR[:,0]<=maxBrig)))[0]
        if (minBrig is None) and (maxBrig is None):
            #indicesBrig = np.array(range(indexStart,(indexEnd+1),1))
            indicesBrig = np.where((BR[:,0]>=BR[indexStart,0])) 

        if autoselect==True:
            # common indices 
            indicesBrig = np.array(range(indexStart,(indexEnd+1),1))[np.nonzero(np.in1d(np.array(range(indexStart,(indexEnd+1),1)), indicesBrig))[0]]
        else:
            indicesBrig = np.array(range(indexStart,(indexEnd+1),1))

        #endBrig = np.where((BR[:,0]>=BR[indicesBrig[0],0]))[0][-1]
        #indicesBrig = indicesBrig[:np.where((indicesBrig==endBrig))[0][0]+1] 

        timeSl = BR[indicesBrig,-1]  #Time slice
        brSl = BR[indicesBrig,0] #Brightness slice

        #Normalize the brightness points
        brSl = (brSl - np.nanmin(brSl)) / (np.nanmax(brSl) - np.nanmin(brSl))
        #Area
        s = [(n+1)**points_size  for n in brSl]
        #Colorbar
        cm = plt.cm.get_cmap(colormap)

        # colors
        norm = plt.Normalize()
        c_ = cm(norm(timeSl))

        if Test=='MAIN' or Test=='ARF':
            #sc = ax1.scatter(BR[indicesBrig,1], BR[indicesBrig,2], s=s, c=timeSl, marker="$\u25EF$", cmap=cm, linewidth=2,  transform=ccrs.PlateCarree())
            sc = ax1.scatter(BR[indicesBrig,1], BR[indicesBrig,2], s=s, c=timeSl, cmap=cm, edgecolor='black', linewidth=1,  alpha=0.8, transform=ccrs.PlateCarree())
  
    #Epicenter in Map
    if hypo==True:
        ax1.plot(evlo, evla, '*', color='red', linewidth=5, markersize=20, markeredgecolor='k', markeredgewidth=1.0, alpha=0.6, transform=ccrs.Geodetic())

    #Time evolution
    ax2 = fig.add_subplot(gs[52:65, 67:])

    if Test=='MAIN' or  Test=='ARF':
        f_out = interp1d(BR[indicesBrig[0]:indicesBrig[-1]+1,-1], BR[indicesBrig[0]:indicesBrig[-1]+1,0])
        t_ = np.linspace(BR[indicesBrig[0],-1], BR[indicesBrig[-1],-1], 5000)
        t_out = f_out(t_)

        # colors of the cmap
        points = np.array([t_, t_out]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm_ = plt.Normalize(t_.min(), t_.max())

        lc = LineCollection(segments, cmap=cm, norm=norm_)
        lc.set_array(t_)
        lc.set_linewidth(3)
        ax2.add_collection(lc)

        if (minBrig is not None) and (maxBrig is not None):
            ylim1 = minBrig
            ylim2 = maxBrig
        if (minBrig is not None) and (maxBrig is None):
            ylim1 = minBrig
            ylim2 = np.nanmax(BR[indicesBrig,0]) + (np.nanmax(BR[indicesBrig,0])*5/100)
        if (minBrig is None) and (maxBrig is not None):
            ylim1 = np.nanmin(BR[indicesBrig,0])
            ylim2 = maxBrig
        if (minBrig is None) and (maxBrig is None):
            ylim1 = np.nanmin(BR[indicesBrig,0])
            ylim2 = np.nanmax(BR[indicesBrig,0]) + (np.nanmax(BR[indicesBrig,0])*5/100)

    ax2.plot(0, BR[indicesBrig[0],0], '*', color='red', linewidth=5, markersize=18, markeredgecolor='k', markeredgewidth=1.5, clip_on=False)

    #ax2.set_xlim([BR[indicesBrig[0],-1], BR[:,-1][BR[:,0]>=BR[indicesBrig[0],0]][-1]])
    ax2.set_xlim([BR[indicesBrig[0],-1], BR[indicesBrig[-1],-1]])

    ax2.set_ylim([ylim1, ylim2])
    ax2.set_xlabel('Relative to Origin Time (s)', fontsize = 12, labelpad=8)
    ax2.set_ylabel('Norm. Brightness', fontsize = 12, labelpad=8)
    ax2.yaxis.set_label_position("right")
    ax2.xaxis.set_label_position("top")
    ax2.tick_params(left = False, right = True , labelleft = False,\
                    labelbottom = False, bottom = False, top = True, labeltop = True, labelright = True)
    ax2.grid(True)


    #Cross 1
    ax3 = fig.add_subplot(gs[68:83, 0:49])

    if Test=='MAIN' or Test=='ARF':
        sc = ax3.scatter(BR[indicesBrig,1], BR[indicesBrig,3], s=s, c=timeSl, cmap=cm, linewidth=1, edgecolor='black', alpha=0.8)

    if hypo==True:
        ax3.plot(evlo, evdepth, '*', color='red', linewidth=5, markersize=20, markeredgecolor='k', markeredgewidth=1.5)

    ax3.set_xlim([min_lon, max_lon])

    if (min_depth is None) and (max_depth is None):
        ax3.set_ylim([math.floor(min(depth)), math.ceil(max(depth))])
    if (min_depth is None) and (max_depth is not None):
        ax3.set_ylim([math.floor(min(depth)), math.ceil(max_depth)])
    if (min_depth is not None) and (max_depth is None):
        ax3.set_ylim([math.floor(min_depth), math.ceil(max(depth))])
    if (min_depth is not None) and (max_depth is not None):
        ax3.set_ylim([math.floor(min_depth), math.ceil(max_depth)])

    ax3.set_xlabel('Longitude (°)', fontsize = 12, labelpad=8)
    ax3.set_ylabel('Depth (km)', fontsize = 12, labelpad=8)
    ax3.invert_yaxis()
    ax3.grid(True)
 
    #Cross 2
    ax4 = fig.add_subplot(gs[68:83, 51:])

    if Test=='MAIN' or Test=='ARF':
         sc = ax4.scatter(BR[indicesBrig,2], BR[indicesBrig,3], s=s, c=timeSl, cmap=cm, linewidth=1, edgecolor='black', alpha=0.8)
    ax4.set_xlabel('Latitude (°)', fontsize = 12, labelpad=8)
    ax4.yaxis.set_label_position("right")
    #ax4.set_ylabel('Depth (km)', fontsize = 12, labelpad=8)
    ax4.tick_params(left = False, right = True , labelleft = False,\
                    labelbottom = True, bottom = True, top = False, labeltop = False, labelright = True)

    ax4.set_xlim([min_lat, max_lat])

    if (min_depth is None) and (max_depth is None):
        ax4.set_ylim([math.floor(min(depth)), math.ceil(max(depth))])
    if (min_depth is None) and (max_depth is not None):
        ax4.set_ylim([math.floor(min(depth)), math.ceil(max_depth)])
    if (min_depth is not None) and (max_depth is None):
        ax4.set_ylim([math.floor(min_depth), math.ceil(max(depth))])
    if (min_depth is not None) and (max_depth is not None):
        ax4.set_ylim([math.floor(min_depth), math.ceil(max_depth)])

    if hypo==True:
        ax4.plot(evla, evdepth, '*', color='red', linewidth=5, markersize=20, markeredgecolor='k', markeredgewidth=1.5)

    ax4.invert_yaxis()
    ax4.grid(True)

    #Stations azimuths
    ##################
    ax5 = fig.add_subplot(gs[20:40, 75:95], projection='polar')
    #ax5.set_ylim([0, 100])

    ax5.set_title('Stations azimuthal distribution', fontweight='bold')

    angles = []; radii = []
    #Plot stations used
    for sta in stations_used:
        dist_theta = gps2dist_azimuth(inv.select(station=sta.split('.')[1])[0][0].latitude, inv.select(station=sta.split('.')[1])[0][0].longitude,\
                                      evla, evlo)
        #ax5.scatter(dist_theta[0]/1000, dist_theta[2])
        angles.append(dist_theta[2]); radii.append(dist_theta[0]/1000);

    circular_hist(ax5, np.radians(angles), bins=16, density=True, offset=0, gaps=False)

    # Colorbar
    ##################### 
    ax6 = fig.add_subplot(gs[0:8, 67:])

    ax6.set_yticklabels([])
    ax6.set_xticklabels([])
    ax6.set_yticks([])
    ax6.set_xticks([])

    cbaxes = inset_axes(ax6, '70%', '23%', loc = 'lower center')

    cb = plt.colorbar(sc, cax=cbaxes, orientation="horizontal")
    cb.set_label(label='Max Peaks - Time (s)', size='large', weight='bold')
    cb.ax.tick_params(labelsize='medium')
    cbaxes.xaxis.set_ticks_position("top")
    cbaxes.xaxis.set_label_position("top")

    ylabels = np.arange(math.floor(min_lat), math.ceil(max_lat + 0.20), 0.20)
    for yl in range(len(ylabels)):
        ylabels[yl] = str(round(ylabels[yl],1))

    xlabels = np.arange(math.floor(min_lon), math.ceil(max_lon + 0.20), 0.20)
    for xl in range(len(xlabels)):
        xlabels[xl] = str(round(xlabels[xl],1))

    if meridian==True:
        #Plot meridians and parallels
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=0.1, color='gray', alpha=0.2, linestyle='--', draw_labels=True)
        gl.top_labels = False
        gl.left_labels = True
        gl.right_labels =False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator(list(xlabels))
        gl.ylocator = mticker.FixedLocator(list(ylabels))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

    ### INFO BOX
    ############
    if info_box==True:
    ## Add the legend

        if Test=='MAIN' or Test=='ARF':
            from matplotlib.offsetbox import (AnchoredOffsetbox, DrawingArea, HPacker,
                                  TextArea, VPacker)
            box1 = TextArea("Information area", textprops=dict(color="k", fontsize='x-large', fontweight='bold'))

            if Test=='MAIN':
                box2 = TextArea("Description: SSA Results", textprops=dict(color="k", fontsize='large'))
            if Test=='ARF':
                box2 = TextArea("Description: Array Response Function", textprops=dict(color="k", fontsize='large'))

            box3 = TextArea("Event Details", textprops=dict(color="k", fontsize='large', fontweight='bold'))
            box4 = TextArea("Origin: " + str(time), textprops=dict(color="k", fontsize='large'))
            box5 = TextArea("Hypo: " + str(round(evlo,2)) + '/' + str(round(evla,2)) + '/' + str(round(evdepth,2)), textprops=dict(color="k", fontsize='large'))
            box6 = TextArea("Maximum Brightness Details", textprops=dict(color="k", fontsize='large', fontweight='bold'))
            box7 = TextArea("Value: " + str(np.around(BR[np.nanargmax(BR[:,0]),0],2)), textprops=dict(color="k", fontsize='large'))
            box8 = TextArea("Lon./Lat./Depth: " + str(np.around(BR[np.nanargmax(BR[:,0]),1],2)) +'/'+str(np.around(BR[np.nanargmax(BR[:,0]),2],2))\
                            + '/' + str(np.around(BR[np.nanargmax(BR[:,0]),3],2)), textprops=dict(color="k", fontsize='large'))
    
            box = VPacker(children=[box1, box2, box3, box4, box5, box6, box7, box8],
                  align="center",
                  pad=0, sep=5)
            anchored_box = AnchoredOffsetbox(loc='lower left',
                                 child=box, pad=0.,
                                 frameon=True,
                                 bbox_to_anchor=(-0.40, -2.0),
                                 bbox_transform=ax4.transAxes,
                                 borderpad=0.,
                                 )

            ax4.add_artist(anchored_box)
    fig.suptitle('Maximum Brightness Per Time Step', fontsize=15, fontweight='bold')
    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)

    return
