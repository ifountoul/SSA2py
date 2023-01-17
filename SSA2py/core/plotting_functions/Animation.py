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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.mpl.geoaxes

from scipy.stats import binned_statistic_2d

def brFiles(data_path):
    """
    Brightness Snapshots per time Step
    
    Arguments: 
    ------
    data_path: Path with the SSA.npy results 

    """

    #Just check that we have .npy results in the directory
    files_list = os.listdir(data_path)
    files_data = []
    for f in files_list:
        if f.endswith('.npy') and f!='out_Max.npy' and f!='tt.npy' and f!='grid.npy':
            files_data.append(f)
    if len(files_data)>0:
        pass
    else:
        print('No .npy files in the directory. Continue....')
        return

    #Sort the list based on time
    files_data.sort(key=lambda x: float('.'.join(x.replace('_', '.').split('.')[1:3])))

    #Add also the path
    files_data = [os.path.join(data_path, files_data[i]) for i in range(len(files_data))]

    return files_data

def brightAnimation(data, filepath):
    """
    Brightness Animation
    Arguments:
    ------
    data: List with .npy
    filepath: Path to save snapshots and animation

    Returns:
    ------
    Snapshots and animation

    """

    data_path = os.path.split(data[0])[0]

    m = np.load(os.path.join(data_path, 'out_Max.npy'), allow_pickle=True)

    xmax = m[:, -1]; ymax = m[:,0];

    #Open the Figure
    fig = plt.figure(constrained_layout=False, figsize=(8,8))

    def animate(i):
        #Clear the Figure
        fig.clf()

        #Load data (Rows=Grid)
        d_ = np.load(data[i])

        #Gather Data for Maximum Map View

        #Find the unique Lat, Lon
        lat = np.unique(d_[:,-3]); lon = np.unique(d_[:,-4]);

        ret = binned_statistic_2d(d_[:,-4], d_[:,-3], d_[:,0], statistic='max', bins=len(lon))

        # Add the map
        ax = plt.axes(projection=ccrs.PlateCarree(), aspect="auto")

        extent = np.min(lon), np.max(lon), np.min(lat), np.max(lat)
        sc = ax.imshow(ret.statistic.T, origin='lower', extent=extent,\
                   cmap='viridis', interpolation='spline16', vmin=0, vmax=ymax.max(), alpha=1.0)

        # Get the coastline
        GSHHS = cfeature.GSHHSFeature(scale='full')

        # Add coastline
        ax.add_feature(GSHHS, linewidth=1.0, facecolor = 'none', zorder=1)

        # Plot the point
        ax.plot(m[i][1], m[i][2], 'o', color='red', linewidth=5, markersize=5, markeredgecolor='k',\
                markeredgewidth=1.0, alpha=0.5, transform=ccrs.Geodetic())

        gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', alpha=0.2, linestyle='-', draw_labels=True)
        gl.top_labels = False
        gl.left_labels = False
        gl.right_labels =True
        gl.xlines = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Add the inset
        box = inset_axes(ax, '40%', '18%', loc = 'lower left')
        box.set_facecolor([1,1,1,0.7])
        box.set_xticks([])
        box.set_yticks([])

        braxes = inset_axes(box, '80%', '55%', loc = 'upper center')
        braxes.plot(xmax,ymax, color='black')

        braxes.axvline(x=d_[1,-1], ymin=0, ymax=10**100, color='red')
        braxes.set_xlim(min(xmax),max(xmax))

        braxes.set_ylabel('Max. Br.', fontsize=8, fontweight='bold')
        braxes.set_xlabel('Time (s)', fontsize=8, fontweight='bold')
        braxes.spines['right'].set_visible(False)
        braxes.spines['top'].set_visible(False)
        braxes.yaxis.set_ticks_position('none')
        braxes.set_yticks([])
        braxes.xaxis.set_ticks_position('bottom')
        braxes.patch.set_alpha(0.8)
        braxes.tick_params(axis='x', labelsize= 7)
        braxes.tick_params(axis='y', labelsize= 7)

        #Snapshot text
        ax.text(0.1, 1.03, 'Time (s): ', horizontalalignment='center', verticalalignment='center',\
                fontsize=12, fontweight='bold', transform=ax.transAxes)
        ax.text(0.24, 1.03, str(d_[1,-1]), horizontalalignment='center', verticalalignment='center',\
                fontsize=12, fontweight='bold', transform=ax.transAxes)

        # add here the colorbar
        box = inset_axes(ax, '40%', '12%', loc = 'lower right')
        box.set_facecolor([1,1,1,0.7])
        box.set_xticks([])
        box.set_yticks([])

        cbaxes = inset_axes(box, '80%', '20%', loc = 'lower center')

        cb = plt.colorbar(sc, cax=cbaxes ,orientation="horizontal")
        cb.set_label(label='Maximum Brightness', size='small', weight='bold')
        cb.ax.tick_params(labelsize='medium')
        cbaxes.xaxis.set_ticks_position("top")
        cbaxes.xaxis.set_label_position("top")

        fig.suptitle('Maximum Brightness Animation', weight='bold', size=15)
        return

    #How many frames per second?
    intervals = round(abs(float('.'.join(os.path.split(data[1])[-1].replace('_', '.').split('.')[1:3])) -\
                float('.'.join(os.path.split(data[2])[-1].replace('_', '.').split('.')[1:3]))),1)


    anim = FuncAnimation(fig, animate, frames= len(data), interval = int(intervals* 1000))
    anim.save(os.path.join(filepath))

    return







