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
##########

import glob, os, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation 
import cartopy.feature as cfeature
from matplotlib import rc
from scipy.stats import binned_statistic_2d
import cartopy.crs as ccrs

# Main Function
###############

def _animation_(Data, evla, evlo, evdp, min_lon=None, min_lat=None, max_lon=None, max_lat=None,\
                min_depth=None, max_depth=None, hypo=True, colormap='viridis', mincolor=None,\
                maxcolor=None,outpath='.', filename='animation', fileformat='mp4'):
    """
    Input:
    -------

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

    Output:
    -------
    Saved Animation.

    """   

 
    # Use glob to get a list of all npy files in the directory
    npy_files = [file for file in glob.glob(os.path.join(Data, "out*.npy")) if "out_Max.npy" not in file]

    # Sort this list 
    def extract_number_from_path(path):
        match = re.search(r'out_(-?\d+\.\d+)\.npy', path)
    
        if match:
            return float(match.group(1))
        return -1  # Return a default value if no match is found

    # Sorted files
    npy_files = sorted(npy_files, key=extract_number_from_path)

    # Extract the timestep values
    timesteps = list(map(extract_number_from_path, npy_files))

    # Calculate the intervals between consecutive timesteps
    intervals = np.round([timesteps[i + 1] - timesteps[i] for i in range(len(timesteps) - 1)],2)

    # Number of frames
    frames = len(npy_files)

    interval_ = int(intervals[0] * 1000) # In milliseconds
 
    # Start the Plot
    ################

    fig, axs = plt.subplots(2, 2, figsize=(10, 10), gridspec_kw={'width_ratios': [2.5, 1], 'height_ratios': [2.5, 1]})
    axs = axs.flatten()

    # Read Maximum Brightness
    maxBr = np.load(os.path.join(Data, "out_Max.npy"))

    d_ = np.load(npy_files[0])

    # Values
    lat_ = d_[:,2]; lon_ = d_[:,1]; depth_ = d_[:,3]; br_ = d_[:,0]

    if (mincolor is None) and (maxcolor is None):
        vmin = np.min(maxBr[:,0]) - (0.1*np.min(maxBr[:,0]))
        vmax = np.max(maxBr[:,0])
    else:
        vmin = mincolor
        vmax = maxcolor


    # Plot 1 #######
    axs[0] = plt.subplot(2,2,1, projection=ccrs.PlateCarree())

    ret = binned_statistic_2d(lon_, lat_, br_, statistic='max', bins=len(np.unique(lon_)))

    extent = np.min(np.unique(lon_)), np.max(np.unique(lon_)), np.min(np.unique(lat_)), np.max(np.unique(lat_))

    im1 = axs[0].imshow(ret.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='equal',\
                        alpha=1.0, animated=True, vmin=vmin,\
                        vmax=vmax)

    coast = cfeature.GSHHSFeature(scale='h', alpha=0.5)
    axs[0].add_feature(coast)

    # Add the time
    props = dict(boxstyle='square', facecolor='white', alpha=0.9)

    # place a text box in upper left in axes coords
    timetext = axs[0].text(0.03, 0.97, str(timesteps[0]), transform=axs[0].transAxes, fontsize=18,
                           verticalalignment='top', bbox=props)

    # Plot 2 #######
    ret = binned_statistic_2d(lon_, depth_, br_, statistic='max', bins=(len(np.unique(lon_)), len(np.unique(depth_))))

    extent = np.min(np.unique(lon_)), np.max(np.unique(lon_)), np.min(np.unique(depth_)), np.max(np.unique(depth_))

    im2 = axs[2].imshow(ret.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic',\
                        aspect='auto', alpha=1.0, animated=True, vmin=vmin,\
                        vmax=vmax)

    axs[2].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[2].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)

    # Plot 3 ########

    ret = binned_statistic_2d(depth_, lat_, br_, statistic='max', bins=(len(np.unique(depth_)), len(np.unique(lat_))))

    extent = np.min(np.unique(depth_)), np.max(np.unique(depth_)), np.min(np.unique(lat_)), np.max(np.unique(lat_))

    im3 = axs[1].imshow(ret.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic',\
                        aspect='auto', alpha=1.0, animated=True, vmin=vmin,\
                        vmax=vmax)

    axs[1].set_xlabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[1].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12, rotation=270)
    axs[1].tick_params(axis='y', labelright=True, labelleft=False, right=True, left=False)
    axs[1].yaxis.set_label_position("right")
    axs[1].yaxis.set_label_coords(1.40, 0.5)

    # Maximum Brightness plot
    axs[3].plot(maxBr[:,4], maxBr[:,0], color='k')
    axs[3].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
    axs[3].set_ylabel('Max. Brightness', fontweight='bold', fontsize=12, rotation=270)
    axs[3].yaxis.set_label_position("right")
    axs[3].tick_params(axis='y', labelright=True, labelleft=False, right=True, left=False)
    axs[3].yaxis.set_label_coords(1.40, 0.5)
    axs[3].set_xlim(np.min(maxBr[:,4]), np.max(maxBr[:,4]))
    axs[3].set_ylim(np.min(maxBr[:,0]), np.max(maxBr[:,0]))

    # Fill area
    axs[3].fill_between(maxBr[:,4], maxBr[:,0], color='lightgray', alpha=0.5)
    axs[3].grid()

    vline = axs[3].axvline(x=timesteps[0], color='r')

    fig.suptitle("Maximum Brightness Per Timestep Animation", fontsize=20, fontweight='bold')

    # Colorbar 
    cbar_ax = fig.add_axes([0.70, 0.93, 0.20, 0.02])  # [x, y, width, height]
    cbar = plt.colorbar(im1, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('Maximum Brightness', fontweight='bold', fontsize=10)

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

    if hypo==True:
        axs[0].plot(evlo, evla, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6, transform=ccrs.Geodetic())
        
        axs[2].plot(evlo, evdp, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6)

        axs[1].plot(evdp, evla, '*', color='red', linewidth=5, markersize=20,\
                    markeredgecolor='k', markeredgewidth=1.0, alpha=0.6)


    def updatefig(n):
    
        d_ = np.load(npy_files[n])
    
        lat_ = d_[:,2]; lon_ = d_[:,1]; depth_ = d_[:,3]; br_ = d_[:,0]
    
        # Plot 1
        ret = binned_statistic_2d(lon_, lat_, br_, statistic='max', bins=len(np.unique(lon_)))
   
        im1.set_array(ret.statistic.T)
    
        timetext.set_text(str(timesteps[n]))
    
        # Plot 2
    
        ret = binned_statistic_2d(lon_, depth_, br_, statistic='max',\
                                  bins=(len(np.unique(lon_)), len(np.unique(depth_))))
     
        im2.set_array(ret.statistic.T)
    
        # Plot 3
    
        ret = binned_statistic_2d(depth_, lat_, br_, statistic='max',\
                                  bins=(len(np.unique(depth_)), len(np.unique(lat_))))
     
        im3.set_array(ret.statistic.T)
    
        vline.set_xdata(timesteps[n])
    
        return im1, im2, im3, vline, timetext,

    ani = animation.FuncAnimation(fig, updatefig, frames = len(npy_files), interval= interval_)

    ani.save(os.path.join(outpath, filename+'.'+fileformat), writer='ffmpeg', dpi=300)




    return
