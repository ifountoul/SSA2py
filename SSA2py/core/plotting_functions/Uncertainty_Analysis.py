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

import numpy as np
import os, math
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.stats import binned_statistic_2d
from obspy.geodetics.base import gps2dist_azimuth
from sklearn import preprocessing as pre
from scipy.signal import savgol_filter

def _heatmap_statistics(values, grid, type_, colormap, outpath, filename, fileformat, dpi):
    """
    Plot HeatMaps for the Unceratainty Analysis Results.

    Inputs:
    -------
    values: array
        Values to plot
    grid: array
        Grid of BP
    type_: str
        Type Of Uncertainty Analysis
    colormap: string
        Colomap to use
    outpath: string
        Path to Save
    filename: string
        Filename to use
    fileformat: string
        Fileformat
    dpi: int

    Output:
    -------
        Figure

    """     
  

    # Mean
    mean_values = np.mean(values, axis=1)
    
    ret_mean = binned_statistic_2d(grid[:,0], grid[:,1], mean_values, statistic='mean',\
                                   bins=len(np.unique(grid[:,0])))
    
    # Median
    median_values = np.median(values, axis=1)
    
    ret_median = binned_statistic_2d(grid[:,0], grid[:,1], median_values, statistic='median',\
                                     bins=len(np.unique(grid[:,0])))
    
    # Std
    std_values = np.std(values, axis=1)
    
    ret_std = binned_statistic_2d(grid[:,0], grid[:,1], std_values, statistic='std',\
                                  bins=len(np.unique(grid[:,0])))  
    
    # Max
    max_values = np.max(values, axis=1)
    
    ret_max = binned_statistic_2d(grid[:,0], grid[:,1], max_values, statistic='max',\
                                  bins=len(np.unique(grid[:,0])))

    # Plot
    extent = np.min(np.unique(grid[:,0])), np.max(np.unique(grid[:,0])),\
             np.min(np.unique(grid[:,1])), np.max(np.unique(grid[:,1]))

    fig, axs = plt.subplots(2,2, figsize=(10,10))
    axs = axs.flatten()
    
    im1 = axs[0].imshow(ret_mean.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='equal',\
                        alpha=1.0)
    
    im2 = axs[1].imshow(ret_median.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='equal',\
                        alpha=1.0)
    
    im3 = axs[2].imshow(ret_std.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='equal',\
                        alpha=1.0)
    
    im4 = axs[3].imshow(ret_max.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='equal',\
                        alpha=1.0)
    
    axs[0].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[0].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12)
    
    axs[1].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[1].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12)
    
    axs[2].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[2].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12)
    
    axs[3].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[3].set_ylabel('Latitude (°)', fontweight='bold', fontsize=12)
    
    axs[0].set_title('Mean {}'.format(type_), fontweight='bold')
    axs[1].set_title('Median {}'.format(type_), fontweight='bold')
    axs[2].set_title('Std {}'.format(type_), fontweight='bold')
    axs[3].set_title('Max. {}'.format(type_), fontweight='bold')
    
    plt.colorbar(im1, ax=axs[0], shrink=0.5)
    plt.colorbar(im2, ax=axs[1], shrink=0.5)
    plt.colorbar(im3, ax=axs[2], shrink=0.5)
    plt.colorbar(im4, ax=axs[3], shrink=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)
 
    return

def _heatmap_statistics_ver(values, grid, type_, colormap, outpath, filename, fileformat, dpi):
    """
    Plot HeatMaps for the Unceratainty Analysis Results.

    Inputs:
    -------
    values: array
        Values to plot
    grid: array
        Grid of BP
    type_: str
        Type Of Uncertainty Analysis
    colormap: string
        Colomap to use
    outpath: string
        Path to save
    filename: string
        Filename to use
    fileformat: string
        Fileformat
    dpi: int

    Output:
    -------
        Figure

    """

    # Mean
    mean_values = np.mean(values, axis=1)
    
    ret_mean = binned_statistic_2d(grid[:,0], grid[:,2], mean_values, statistic='mean',\
                                   bins=(len(np.unique(grid[:,0])), len(np.unique(grid[:,2]))))
    
    
    
    # Median
    median_values = np.median(values, axis=1)
    
    ret_median = binned_statistic_2d(grid[:,0], grid[:,2], median_values, statistic='median',\
                                     bins=(len(np.unique(grid[:,0])), len(np.unique(grid[:,2]))))
    
    # Std
    std_values = np.std(values, axis=1)
    
    ret_std = binned_statistic_2d(grid[:,0], grid[:,2], std_values, statistic='std',\
                                  bins=(len(np.unique(grid[:,0])), len(np.unique(grid[:,2]))))  
    
    # Max
    max_values = np.max(values, axis=1)
    
    ret_max = binned_statistic_2d(grid[:,0], grid[:,2], max_values, statistic='max',\
                                  bins=(len(np.unique(grid[:,0])), len(np.unique(grid[:,2]))))
    
    
    extent = np.min(np.unique(grid[:,0])), np.max(np.unique(grid[:,0])),\
             np.min(np.unique(grid[:,2])), np.max(np.unique(grid[:,2]))
    
    
    fig, axs = plt.subplots(2,2, figsize=(10,8))
    axs = axs.flatten()
    
    im1 = axs[0].imshow(ret_mean.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='auto',\
                        alpha=1.0)
    
    im2 = axs[1].imshow(ret_median.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='auto',\
                        alpha=1.0)
    im3 = axs[2].imshow(ret_std.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='auto',\
                        alpha=1.0)
    im4 = axs[3].imshow(ret_max.statistic.T, origin='lower', extent=extent,\
                        cmap=colormap, interpolation='bicubic', aspect='auto',\
                        alpha=1.0)
    
    
    axs[0].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[0].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[0].invert_yaxis()
    
    axs[1].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[1].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[1].invert_yaxis()
    
    axs[2].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[2].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[2].invert_yaxis()
    
    axs[3].set_ylabel('Depth (km)', fontweight='bold', fontsize=12)
    axs[3].set_xlabel('Longitude (°)', fontweight='bold', fontsize=12)
    axs[3].invert_yaxis()
    
    axs[0].set_title('Mean {}'.format(type_), fontweight='bold')
    axs[1].set_title('Median {}'.format(type_), fontweight='bold')
    axs[2].set_title('Std {}'.format(type_), fontweight='bold')
    axs[3].set_title('Max. {}'.format(type_), fontweight='bold')
    
    plt.colorbar(im1, ax=axs[0], shrink=0.5)
    plt.colorbar(im2, ax=axs[1], shrink=0.5)
    plt.colorbar(im3, ax=axs[2], shrink=0.5)
    plt.colorbar(im4, ax=axs[3], shrink=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)

    return


def max_br_analysis(paths, outpath, filename, fileformat, dpi):
    """
    Plot Maximum Brightness Uncertainty Analysis.

    Inputs:
    -------
    paths: list
        List with Paths
    outpath: string
        Path to save
    filename: string
        Filename to use
    fileformat: string
        Fileformat
    dpi: int

    Output:
    -------
    Figure

    """

    fig, axs = plt.subplots(2,2, figsize=(15,10))
    axs = axs.flatten()

    # Main Maximum Brightness
    max_br_ = np.load(os.path.join(paths[1], 'out_Max.npy'))

    mean_max = []
    distance = []
    depth_dist = []

    for i in range(len(os.listdir(os.path.join(paths[0], 'Detailed_Solutions')))):
         max_br = np.load(os.path.join(paths[0], 'Detailed_Solutions', str(i), 'out_Max.npy'))
         mean_max.append(max_br[:,0])
         # Distance in earth
         distance.append([gps2dist_azimuth(max_br_[j,2], max_br_[j,1], max_br[j,2], max_br[j,1])[0]/1000\
                          for j in range(max_br.shape[0])])
         # Depth Distance
         depth_dist.append([np.abs(max_br_[j,3]-max_br[j,3]) for j in range(max_br.shape[0])])
    
         axs[0].plot(max_br[:,-1], max_br[:,0], color='gray', alpha=0.3)
    
    # Plot 1    
    # dummy
    axs[0].plot(max_br[:,-1], max_br[:,0], color='gray', alpha=0.3, label='Resampling Solutions')
    #Main     
    axs[0].plot(max_br_[:,-1], max_br_[:,0], color='k', label='Main Solution')    
    #Mean
    mean_max = np.mean(np.vstack(mean_max).T, axis=1)
    axs[0].plot(max_br_[:,-1], mean_max, color='r', label='Mean Solution', ls='--')

    # Plot 2
    distance = np.vstack(distance).T
    smooth_points = math.ceil(len(distance)*0.05)
    # Mean
    mean = np.mean(distance, axis=1)
    # Max
    _max_ = np.max(distance, axis=1)

    axs[1].plot(max_br_[:,-1], _max_, color='red', label='Maximum Distance')
    axs[1].plot(max_br_[:,-1], mean, color='k', label='Mean Distance', ls='--')

    # Plot 3
    depth_dist = np.vstack(depth_dist).T

    # Mean
    mean = np.mean(depth_dist, axis=1)
    # Max
    _max_ = np.max(depth_dist, axis=1)

    axs[2].plot(max_br_[:,-1], _max_, color='red', label='Maximum Distance')
    axs[2].plot(max_br_[:,-1], mean, color='k', label='Mean Distance', ls='--')

    # Plot 4
    axs[3].plot(max_br_[:,-1], pre.MinMaxScaler().fit_transform(np.var(depth_dist, axis=1).reshape(-1, 1)), color='green', label='Depth Distance')
    axs[3].plot(max_br_[:,-1], pre.MinMaxScaler().fit_transform(np.var(distance, axis=1).reshape(-1, 1)), color='blue', label='Great Circle Distance')

    axs[0].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
    axs[0].set_ylabel('Maximum Brightness', fontweight='bold', fontsize=12)

    axs[1].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
    axs[1].set_ylabel('Great Circle Distance (km)', fontweight='bold', fontsize=12)

    axs[2].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
    axs[2].set_ylabel('Depth Distance (km)', fontweight='bold', fontsize=12)

    axs[3].set_xlabel('Time (s)', fontweight='bold', fontsize=12)
    axs[3].set_ylabel('Normalized Variance', fontweight='bold', fontsize=12)

    axs[0].autoscale(enable=True, axis='both', tight=True)
    axs[1].autoscale(enable=True, axis='both', tight=True)
    axs[2].autoscale(enable=True, axis='both', tight=True)
    axs[3].autoscale(enable=True, axis='both', tight=True)
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    axs[3].legend()
    axs[0].grid()
    axs[1].grid()
    axs[2].grid()
    axs[3].grid()

    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)
    return


