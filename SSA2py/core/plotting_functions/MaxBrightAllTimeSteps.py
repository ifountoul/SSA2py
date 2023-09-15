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

# Imports
#########

import numpy as np
import matplotlib.pyplot as plt
import cartopy, os, glob
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from scipy.interpolate import griddata
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Local Imports
###############

from SSA2py.core import config

# Functions
###########

def find_max(pairs, values, max_values):
    """
    Find Max Values.

    """
    for p, v in zip(pairs, values): 
        key = tuple(p)
        if key in max_values:
            max_values[key] = max(max_values[key], v)
        else:
            max_values[key] = v
    return max_values

def oneD2IMSHOW(x, y, z):
    """
    Transform arrays to imshow format.

    """
    # Generate Grid Points
    X = np.linspace(np.min(x), np.max(x), 2000)
    Y = np.linspace(np.min(y), np.max(y), 2000)
    X, Y = np.meshgrid(X, Y)
    
    # Interpolate the scattered data onto the grid
    Z = griddata((x, y), z, (X, Y), method='cubic')
    return Z, X, Y


def mask(values, mask_value):
    """
    Mask arrays.
 
    """
    threshold =  mask_value * np.max(values)
    _mask_ = values > threshold
    masked_data = np.where(_mask_, np.nan, values)
    return masked_data, _mask_

# Main Functions
################

def plotMaxBrightAllTimeSteps(Data, evla, evlo, evdp,  min_lon=None, min_lat=None, max_lon=None, max_lat=None,\
                              min_depth=None, max_depth=None,\
                              hypo=True, colormap='viridis', mincolor=None, maxcolor=None,\
                              filename='.', outpath='.', fileformat='png', dpi=400, _mask_=True, mask_value = 0.9):
    """

    Plot Maximum Brightness appear at each gridnode at all timesteps.

    Input:
    ------
    Data: string
        Path with data.
    evla: float
        Event Latitude.
    evlo: float
        Event Longitude.
    evdp: float
        Event Depth.
    min_lon: float
        Minimum Lognitude.
    min_lat: float
        Minimum Latitude.
    max_lon: float
        Maximum Longitude.
    max_lat: float
        Maximum Latitude.
    min_depth: float
        Minimum Depth.
    max_depth: float
        Maximum Depth.
    hypo: bool
        Plot hypocenter?
    colormap: string
        What colormap to use?
    mincolor: float
        Minimum colormap value.
    maxcolor: float
        Maximum colormap value.
    filename: string
        Filename of the plot to save.
    outpath: string
        Where to save?
    fileformat: string
        Format to save.
    dpi: float
        Dpi
    _mask_: bool
        Mask the plot.
    mask_value: float
        Percentage?

    Output:
    -------
        Saved Figure.

    """

    # Read Maximum Brightness
    maxBr = np.load(os.path.join(Data, "out_Max.npy"))

    if (mincolor is None) and (maxcolor is None):
        vmin = np.min(maxBr[:,0]) - (0.2*np.min(maxBr[:,0]))
        vmax = np.max(maxBr[:,0])
    else:
        vmin = mincolor
        vmax = maxcolor

    # Use glob to get a list of all npy files in the directory
    npy_files = [file for file in glob.glob(os.path.join(Data, "out*.npy")) if "out_Max.npy" not in file]    

    max_valuesxy = {}
    max_valuesxz = {}
    max_valuesyz = {}

    for npy in npy_files:
        data = np.load(npy)
  
        # Lat - Lon
        lat_lon_p = data[:, 1:3]
        values = data[:,0]
        max_valuesxy = find_max(lat_lon_p, values, max_valuesxy)    
        # Lon - Depth
        lon_depth_p = data[:, [1,3]]
        max_valuesxz = find_max(lon_depth_p, values, max_valuesxz)
        # Lat - Depth     
        lat_depth_p = data[:, [2,3]]
        max_valuesyz = find_max(lat_depth_p, values, max_valuesyz)

    # Export the values from the dictionaries
    keysList_xy = np.array(list(max_valuesxy.keys()))
    values_xy = np.array(list(max_valuesxy.values()))

    keysList_xz = np.array(list(max_valuesxz.keys()))
    values_xz = np.array(list(max_valuesxz.values()))

    keysList_yz = np.array(list(max_valuesyz.keys()))
    values_yz = np.array(list(max_valuesyz.values()))

    # Start the plot
    #################

    fig, axs = plt.subplots(2, 2, figsize=(10, 10),\
               gridspec_kw={'width_ratios': [2.5, 1], 'height_ratios': [2.5, 1]})

    axs = axs.flatten()

    # Plot 1
    ########

    Z_xy, x, y = oneD2IMSHOW(keysList_xy[:,0], keysList_xy[:,1], values_xy)

    axs[0] = plt.subplot(2,2,1, projection=ccrs.PlateCarree())

    coast = cfeature.GSHHSFeature(scale='h', alpha=0.5)
    axs[0].add_feature(coast)

    im = axs[0].imshow(Z_xy, origin='lower', extent = (np.min(x), np.max(x), np.min(y), np.max(y)),
                       aspect='equal', cmap=colormap, interpolation='bicubic', vmin=vmin, vmax=vmax)

    if _mask_ == True:
        masked_data, mask_ = mask(Z_xy, mask_value)

        axs[0].pcolormesh(x, y, masked_data, cmap='gray', alpha=0.5, shading='auto')
        axs[0].contour(mask_, levels=[0.5], colors='black', linewidths=0.5,
                       extent = (np.min(x), np.max(x), np.min(y), np.max(y)))
    if hypo==True:
        axs[0].plot(evlo, evla, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6, transform=ccrs.Geodetic())

    # Plot 2
    ########

    Z_xz, x, y = oneD2IMSHOW(keysList_xz[:,0], keysList_xz[:,1], values_xz)

    axs[2].imshow(Z_xz, extent = (np.min(x), np.max(x), np.min(y), np.max(y)),\
                  origin='lower', aspect='auto', cmap=colormap, interpolation='bicubic',\
                  vmin=vmin, vmax=vmax)

    if _mask_ == True:
        masked_data, mask_ = mask(Z_xz, mask_value)

        axs[2].pcolormesh(x, y, masked_data, cmap='gray', alpha=0.5, shading='auto')
        axs[2].contour(mask_, levels=[0.5], colors='black', linewidths=0.5,
                      extent = (np.min(x), np.max(x), np.min(y), np.max(y)))

    if hypo==True:
        axs[2].plot(evlo, evdp, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6)


    axs[2].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[2].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)

    # Plot 3
    ########

    Z_yz, x, y = oneD2IMSHOW(keysList_yz[:,1], keysList_yz[:,0], values_yz)

    axs[1].imshow(Z_yz, extent = (np.min(x), np.max(x), np.min(y), np.max(y)),\
                  origin='lower', aspect='auto', cmap=colormap, interpolation='bicubic',\
                  vmin=vmin, vmax=vmax)

    if _mask_ == True:
        masked_data, mask_ = mask(Z_yz, mask_value)

        axs[1].pcolormesh(x, y, masked_data, cmap='gray', alpha=0.5, shading='auto')
        axs[1].contour(mask_, levels=[0.5], colors='black', linewidths=0.5,
                       extent = (np.min(x), np.max(x), np.min(y), np.max(y)))

    if hypo==True:
        axs[1].plot(evdp, evla, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6)


    axs[1].set_xlabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[1].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12, rotation=270)
    axs[1].tick_params(axis='y', labelright=True, labelleft=False, right=True, left=False)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.set_label_coords(1.40, 0.5)

    # Colorbar
    # Move the colorbar outside the subplot
    cbar_ax = fig.add_axes([0.70, 0.30, 0.20, 0.02])  # [x, y, width, height]
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Maximum Brightness', fontweight='bold', fontsize=10)

    # Add text
    text = 'Maximum brightness in the grid at all timesteps.\n Masked image below ' + str(mask_value*100) + '% of the values.' 
    plt.text(0.70, -5, text, size=8, rotation=0.,
             ha="center", va="center",
             bbox=dict(boxstyle="square",
                   ec='k',
                   fc='w',
                   )
            )

    fig.suptitle("Maximum Brightness at All Timesteps", fontsize=20, fontweight='bold')
    # Remove the plot
    axs[3].remove()

    # Change the limits of the plot
    if (min_lon is None) and (min_lat is None) and (max_lon is None) and (max_lat is None)\
        and (min_depth is None) and (max_depth is None):

        pass
    else:
        axs[0].set_xlim([min_lon, max_lon])
        axs[0].set_ylim([min_lat, max_lat])

        axs[1].set_xlim([min_depth, max_depth])
        axs[1].set_ylim([min_lat, max_lat])

        axs[2].set_xlim([min_lon, max_lon])
        axs[2].set_ylim([min_depth, max_depth])

    axs[2].invert_yaxis()

    # Save the figure.
    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)

    return







