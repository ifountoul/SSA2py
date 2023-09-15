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
import os, glob, re, math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Functions
###########

def find_max(pairs, values, max_values):
    for p, v in zip(pairs, values): 
        key = tuple(p)
        if key in max_values:
            max_values[key] = max(max_values[key], v)
        else:
            max_values[key] = v
    return max_values 

def oneD2IMSHOW(x, y, z):
 
    # Generate Grid Points
    X = np.linspace(np.min(x), np.max(x), 2000)
    Y = np.linspace(np.min(y), np.max(y), 2000)
    X, Y = np.meshgrid(X, Y)
    
    # Interpolate the scattered data onto the grid
    Z = griddata((x, y), z, (X, Y), method='cubic')
    return Z, X, Y

def mask(values, mask_value):
    threshold =  mask_value * np.max(values)
    _mask_ = values > threshold
    masked_data = np.where(_mask_, np.nan, values)
    return masked_data, _mask_

# Main Function
###############

def MaxBrightTimeStep_2_(Data,  evla, evlo, evdp, min_lon=None, min_lat=None, max_lon=None,\
                         max_lat=None, usemask='True', rolling_window=0.5, mask_value=0.9,\
                         hypo=True, colormap='viridis', mincolor=None, maxcolor=None,filename='.',\
                         outpath='.', fileformat='png', dpi=400):
    """
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
    usemask: bool
        Use mask?
    rolling_window: float
        Window to find maximum.
    mask_value: float
        Percentage?
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

    Output
    ------
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

    timeranges = np.round(np.arange(timesteps[0], timesteps[-1], rolling_window), 2)

    new_timeranges = []

    for i in range(len(timeranges)):
        start = timeranges[i]

        if i == len(timeranges) - 1:
            if timeranges[i] + rolling_window > timesteps[-1]:
                end = timesteps[-1]
            else:
                end = timeranges[i] + rolling_window
        else:
            end = timeranges[i+1]
        
        new_timeranges.append((start, end))

    # Built the subplots
    max_subs = 9
    num_figs = math.ceil(len(new_timeranges)/max_subs)

    rest = len(new_timeranges)

    for f in range(num_figs):
        fig, axs = plt.subplots(3, 3, figsize=(15, 15))
        axs = axs.flatten()
    
        if f == num_figs - 1:
            loops = rest
        else:
            loops = max_subs
            rest -= max_subs
     
        for l in range(loops):
            # Keep the right files
            s = f*9 + l 
            start = new_timeranges[s][0]
            end = new_timeranges[s][1]
            new_npy = []

            for npy in npy_files:
                fnum = float(npy.split('/')[-1].split('.npy')[0].split('_')[1])
                if fnum>=start and fnum<=end:
                     new_npy.append(npy)
            # Find max in this timestep
            max_values = {}
        
            for _npy in new_npy:
                data = np.load(_npy)
                lat_lon_p = data[:, 1:3]
                values = data[:,0]
                max_values = find_max(lat_lon_p, values, max_values) 
    
            # Export the values from the dictionaries
            keysList = np.array(list(max_values.keys()))
            values = np.array(list(max_values.values())) 
    
            Z, x, y = oneD2IMSHOW(keysList[:,0], keysList[:,1], values)
        
            im = axs[l].imshow(Z, origin='lower', extent = (np.min(x), np.max(x), np.min(y), np.max(y)),
                               aspect='equal', cmap=colormap, interpolation='bicubic', vmin=vmin, vmax=vmax)
       
            if usemask==True: 
                masked_data, _mask_ = mask(Z, mask_value)

                axs[l].pcolormesh(x, y, masked_data, cmap='gray', alpha=0.5, shading='auto')
                axs[l].contour(_mask_, levels=[0.5], colors='black', linewidths=0.5,
                               extent = (np.min(x), np.max(x), np.min(y), np.max(y)))
        
            # Add the time
            props = dict(boxstyle='square', facecolor='white', alpha=0.9)

            # place a text box in upper left in axes coords
            timetext = axs[l].text(0.03, 0.97, str(start) + " / " + str(end) + ' s',\
                                   transform=axs[l].transAxes, fontsize=10,\
                                   verticalalignment='top', bbox=props)

            axs[l].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
            axs[l].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12)

            # Change the limits of the plot
            if (min_lon is None) and (min_lat is None) and (max_lon is None) and (max_lat is None):
                pass
            else:
                axs[l].set_xlim([min_lon, max_lon])
                axs[l].set_ylim([min_lat, max_lat])

            if hypo==True:
                axs[l].plot(evlo, evla, '*', color='red', linewidth=5, markersize=20,\
                            markeredgecolor='k', markeredgewidth=1.0, alpha=0.6)

        
        # Colorbar 
        cbar_ax = fig.add_axes([0.685, 0.93, 0.20, 0.02])  # [x, y, width, height]
        cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.set_label('Maximum Brightness', fontweight='bold', fontsize=10)
    
        fig.suptitle("Maximum Brightness Per Timerange (Fig {}/{})".format(f+1,num_figs), fontsize=20, fontweight='bold')
    
        if l<8:
            for j in range(l+1,max_subs):
                axs[j].remove()

        plt.subplots_adjust(wspace=0.4)
        plt.savefig(os.path.join(outpath, filename+'_' + str(f) +'.'+fileformat), dpi=dpi)










    
