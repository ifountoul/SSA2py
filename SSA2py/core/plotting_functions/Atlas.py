#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright (C) 2023 Ioannis Fountoulakis, Christos Evangelidis

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

import pyproj, shapefile
import cartopy, cartopy.mpl.geoaxes, math, os
import numpy as np

from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.font_manager as font_manager

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from shapely.geometry.polygon import LinearRing
from shapely.ops import transform
from shapely.geometry import Point
from functools import partial

from SSA2py.core import config

def atlas(inv, stations_used, stations_notused, evla, evlo, evdepth,\
          extent_lon=1.1, extent_lat=1.1, extent_lon_inset=10.0, extent_lat_inset=10.0, towns=False,\
          rings_min=0, rings_max=500, rings_step=50, hypo_lines=True,\
          meridians=True, plates=True, filepath='.', filename='atlas', fileformat='png', dpi=400):
    """
    Plots Epicenter and Stations Positions

    Arguments:
    ----------
    inv: Obspy inventory Object
        Metadata
    stations_used: list
        List with stations used (net.station in list).
    stations_notused: list
        List with stations not used (net.station in list).
    evla: float
        Event latitude.
    evlo: float
        Event longitude.
    evdepth: float
        Event depth.
    extent_lon: float
        Extent of the map in degrees (from the more distant station) in longitude.
    extent_lat: float
        Extent of the map in degrees (from the more distant station) in latitude.
    extent_lon_inset: float
        Extent of the inset map in degrees (from the more distant station) in longitude.
    extent_lat_inset: float
        Extent of the inset map in degrees (from the more distant station) in latitude.
    towns: bool 
        Plot towns 
    rings_min:float
        Minimum extent of rings in km.
    rings_max: float
        Maximum extent of rings in km.
    rings_step: float
        Step of rings in km.
    hypo_lines: bool
        Plot lines crossing the hypo?.
    meridians: bool
        Plot meridians and parallels
    plates: bool
        Plot plates in inset map?
    filepath: str
        Path of the output figure (dir).
    filename: str
        Name of the file.
    fileformat: str
        File format such as .png, .pdf.
    dpi: float
        The resolution in dots per inch.
    Returns:
    --------
  
    Figure

    """

    #Calculate the max distance between station and epicenter
    maxdist = [(gps2dist_azimuth(i[0].latitude, i[0].longitude, evla, evlo)[0]/1000) for i in inv]
    maxdist = sorted(maxdist)[-2]

    #Figure
    fig = plt.figure(constrained_layout=False, figsize=(10,10))

    #Map
    ax1 = fig.add_subplot(projection=ccrs.PlateCarree(), aspect="auto")

    #OCEAN
    ocean_10m = cfeature.NaturalEarthFeature('physical', 'ocean', '10m', edgecolor='face', color = 'lightcyan')

    ax1.add_feature(ocean_10m, zorder=0)

    #COASTLINE
    GSHHS = cfeature.GSHHSFeature(scale='high')

    ax1.add_feature(GSHHS, linewidth=0.4, facecolor = 'lemonchiffon', zorder=0)

    #####
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='--', linewidth=0.5, color = 'peru', zorder=0)
    ax1.add_feature(cartopy.feature.LAKES, alpha=0.5, zorder=0)
    ax1.add_feature(cartopy.feature.RIVERS, zorder=0)

    #Add Towns

    if towns==True:
        pass

    #Plot Sectors
    proj_wgs84 = pyproj.Proj('+proj=longlat +datum=WGS84')

    def geodesic_point_buffer(lon, lat, km):
        # Azimuthal equidistant projection
        aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
        project = partial(
            pyproj.transform,
            pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
            proj_wgs84)
        buf = Point(0, 0).buffer(km * 1000)  # distance in metres
        return transform(project, buf).exterior.coords[:]


    for dist in range(rings_min, rings_max, rings_step):
        b = geodesic_point_buffer(evlo, evla, dist)

        #Get them to numpy array
        x = np.zeros((1, len(b)))
        y = np.zeros((1, len(b)))

        for i in range(len(b)):
            x[0,i] = b[i][0]
            y[0,i] = b[i][1]

        ax1.plot(x[0], y[0], '--', linewidth=1, zorder=1, color = 'k', transform=ccrs.Geodetic())
        ax1.annotate(text=str(dist) + ' ' + str('km'), xy =(max(x[0]), evla), xytext= (max(x[0]), evla),\
                     fontsize=10, weight='bold', color='k', rotation=-85, transform=ccrs.Geodetic())

    #if hypo_lines==True:
    #    ax1.axvline(x=evlo, color='k', linestyle='--', linewidth=1)
    #    ax1.axhline(y=evla, color='k', linestyle='--', linewidth=1)

    #Plot stations used
    for sta in stations_used:
        ax1.plot(inv.select(station=sta.split('.')[1])[0][0].longitude,\
                 inv.select(station=sta.split('.')[1])[0][0].latitude,\
                 "^", color='springgreen', markersize=10, markeredgecolor='k', markeredgewidth=0.5, transform=ccrs.Geodetic(), label='Active Seismic Station')

    #Plot unsed stations
    for sta in stations_notused:
        ax1.plot(inv.select(station=sta.split('.')[1])[0][0].longitude,\
                 inv.select(station=sta.split('.')[1])[0][0].latitude,\
                 "^", color='silver', markersize=10, markeredgecolor='k', markeredgewidth=0.5, transform=ccrs.Geodetic(), label='Inactive Seismic Station')
    #Plot meridians and parallels
    if meridians==True:
        gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', alpha=0.2, linestyle='-', draw_labels=True)
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels =True
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator(np.round(np.linspace(math.floor(evlo-kilometer2degrees(maxdist)-0.5), math.ceil(evlo+kilometer2degrees(maxdist)+0.5), 5),2))
        gl.ylocator = mticker.FixedLocator(np.round(np.linspace(math.floor(evla-kilometer2degrees(maxdist)-0.5), math.ceil(evla+kilometer2degrees(maxdist)+0.5), 5),2))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

    #Epicenter in Map
    ax1.plot(evlo, evla, '*', color='red', linewidth=5, markersize=20, markeredgecolor='k', markeredgewidth=1.0, alpha=0.6, transform=ccrs.Geodetic(), label='Epicenter')

    # Add legend
    handles, labels = ax1.get_legend_handles_labels()

    # Keep unique ones
    unique_pairs = {}

    # Iterate through handles and labels
    for handle, label in zip(handles, labels):
        unique_pairs[label] = handle

    # Create lists for unique handles and labels
    unique_handles = list(unique_pairs.values())
    unique_labels = list(unique_pairs.keys())

    font = font_manager.FontProperties(weight='bold', style='italic', size=12)
    ax1.legend(unique_handles, unique_labels, loc='lower left', title="Legend", title_fontproperties=font, prop={'size': 12})

    #Inset map
    ax_inset = fig.add_axes([0.712, 0.684, 0.188, 0.203], projection=ccrs.PlateCarree())
    ax_inset.set_extent([evlo-kilometer2degrees(maxdist)-extent_lon_inset, evlo+kilometer2degrees(maxdist)+extent_lon_inset,\
                   evla-kilometer2degrees(maxdist)-extent_lat_inset, evla+kilometer2degrees(maxdist)+extent_lat_inset], crs=ccrs.PlateCarree())

    #COASTLINE
    GSHHS = cfeature.GSHHSFeature(scale='full')
    ax_inset.add_feature(GSHHS, linewidth=0.3, facecolor = 'lightgray', alpha=0.6, zorder=0)

    #PLATES
    if plates==True:
        try:
            plates = shapefile.Reader(os.path.join(config.cfg['Plotting']['Save Layers'],'PLATES', 'PB2002_boundaries.shp'))
            for tmp in plates.shapeRecords():
                x = [i[0] for i in tmp.shape.points[:]]
                y = [i[1] for i in tmp.shape.points[:]]
                ax_inset.plot(x, y, color='k', lw=0.6, transform= ccrs.PlateCarree())
        except Exception as e:
            config.logger.warning(e)
            config.logger.warning('Plates Shapefile error')

    lonmin = evlo-kilometer2degrees(maxdist)-extent_lon
    lonmax = evlo+kilometer2degrees(maxdist)+extent_lon
    latmin = evla-kilometer2degrees(maxdist)-extent_lat
    latmax = evla+kilometer2degrees(maxdist)+extent_lat

    #Add red borders
    nvert = 100
    lons = np.r_[np.linspace(lonmin, lonmin, nvert),
                 np.linspace(lonmin, lonmax, nvert),
                 np.linspace(lonmax, lonmax, nvert)].tolist()
    lats = np.r_[np.linspace(latmin, latmax, nvert),
                 np.linspace(latmax, latmax, nvert),
                 np.linspace(latmax, latmin, nvert)].tolist()

    ring = LinearRing(list(zip(lons, lats)))
    ax_inset.add_geometries([ring], ccrs.PlateCarree(),
                   facecolor='none', edgecolor='red', linewidth=1)
    #ax_inset.outline_patch.set_linewidth(3.0)

    proj = ccrs.Orthographic(central_longitude=evlo, central_latitude=evla)
    ax_inset_globe = fig.add_axes([0.112, 0.784, 0.185, 0.205], projection= proj) 
    #ax_inset_globe.add_feature(cartopy.feature.LAND, facecolor = 'lemonchiffon')
    #ax_inset_globe.add_feature(cartopy.feature.OCEAN, facecolor = 'lightcyan')
    ax_inset_globe.add_feature(cartopy.feature.COASTLINE)
    ax_inset_globe.gridlines()
    ax_inset_globe.plot(evlo, evla, '*', color='red', linewidth=5, markersize=20, markeredgecolor='k', markeredgewidth=1.0, alpha=0.6, transform=ccrs.Geodetic())
    ax_inset_globe.set_global()

    plt.savefig(os.path.join(filepath, filename+'.'+fileformat), dpi=dpi)

    return
