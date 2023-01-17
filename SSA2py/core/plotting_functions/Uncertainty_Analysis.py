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

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from obspy.geodetics.base import kilometer2degrees
from obspy.geodetics.base import gps2dist_azimuth


import cartopy, cartopy.mpl.geoaxes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def Max_Bright_Uncertainty(path, out_path, filename, fileformat='pdf', dpi=400):
    """
    Plot the Uncertainty Analysis Results

    Arguments:
    ----------
    path: str
        Input data paths
    out_path: str
        Output path to put plots


    """

    SSA = np.load(os.path.join(path, 'Main_Solution', 'out_Max.npy'))

    fig, axs = plt.subplots(3, 3, figsize=(13, 8), constrained_layout=True)

    axs[0,0].plot(SSA[:,13], SSA[:,6]*111.1, label='Horizontal')
    axs[0,0].plot(SSA[:,13], SSA[:,10], label='Vertical')

    axs[0,0].set_xlabel('Time (s)')
    axs[0,0].set_ylabel('C.I. (km)')
    axs[0,0].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))
    axs[0,0].legend()

    axs[0,1].plot(SSA[:,13], SSA[:,8]*111.1, label='Horizontal')
    axs[0,1].plot(SSA[:,13], SSA[:,12], label='Vertical')

    axs[0,1].set_xlabel('Time (s)')
    axs[0,1].set_ylabel('S.E. (km)')
    axs[0,1].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))

    axs[0,2].plot(SSA[:,13], SSA[:,7]*111.1, label='Horizontal')
    axs[0,2].plot(SSA[:,13], SSA[:,11], label='Vertical')

    axs[0,2].set_xlabel('Time (s)')
    axs[0,2].set_ylabel('ST.D. (km)')
    axs[0,2].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))

    axs[1,0].hist(SSA[:,6]*111.1, bins='sqrt', label='Horizontal')
    axs[1,0].hist(SSA[:,10], bins='sqrt', alpha=0.7, label='Vertical')

    axs[1,0].axvline(np.mean(SSA[:,6]*111.1), ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()], label='Mean Hor.')
    axs[1,0].axvline(np.mean(SSA[:,10]), color='orange', ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()], label='Mean Ver.')

    axs[1,0].set_xlabel('C.I. (km)')
    axs[1,0].set_ylabel('Measurements')
    axs[1,0].autoscale(enable=True, axis='x', tight=True)
    axs[1,0].legend()

    axs[1,1].hist(SSA[:,8]*111.1, bins='sqrt', label='Horizontal')
    axs[1,1].hist(SSA[:,12], bins='sqrt', alpha=0.7, label='Vertical')

    axs[1,1].axvline(np.mean(SSA[:,8]*111.1), ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])
    axs[1,1].axvline(np.mean(SSA[:,12]), color='orange', ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])

    axs[1,1].set_xlabel('S.E. (km)')
    axs[1,1].autoscale(enable=True, axis='x', tight=True)

    axs[1,2].hist(SSA[:,7]*111.1, bins='sqrt', label='Horizontal')
    axs[1,2].hist(SSA[:,11], bins='sqrt', alpha=0.7, label='Vertical')

    axs[1,2].axvline(np.mean(SSA[:,7]*111.1), ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])
    axs[1,2].axvline(np.mean(SSA[:,11]), color='orange', ls='--', path_effects=[pe.Stroke(linewidth=3, foreground='k'), pe.Normal()])

    axs[1,2].set_xlabel('ST.D. (km)')
    axs[1,2].autoscale(enable=True, axis='x', tight=True)

    #########################
    #########################

    axs[2,0].plot(SSA[:,13], SSA[:,1], label='Brightness')
    axs[2,0].set_xlabel('Time (s)')
    axs[2,0].set_ylabel('C.I.')
    axs[2,0].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))
    axs[2,0].legend()

    axs[2,1].plot(SSA[:,13], SSA[:,2])
    axs[2,1].set_xlabel('Time (s)')
    axs[2,1].set_ylabel('S.E.')
    axs[2,1].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))

    axs[2,2].plot(SSA[:,13], SSA[:,3])
    axs[2,2].set_xlabel('Time (s)')
    axs[2,2].set_ylabel('ST.D.')
    axs[2,2].set_xlim(np.min(SSA[:,13]), np.max(SSA[:,13]))

    fig.suptitle('Uncertainty Analysis For Max. Brightness', fontsize=18, weight='bold')

    plt.savefig(os.path.join(out_path, filename+'.'+fileformat), dpi=dpi)

    return


def Max_Bright_Map(path, out_path, evla, evlo, evdepth, filename, fileformat='pdf', dpi=400):
    """
    Plot the Uncertainty Analysis Map
    Only for the Maximum max brightness!

    Arguments:
    ----------
    path: str
        Input data paths
    out_path: str
        Output path to put plots


    """

    # open figure
    fig, ax = plt.subplots(1,3, constrained_layout=False, figsize=(10,3))

    ax[0].scatter(evla, evlo)
    ax[1].scatter(evlo, evdepth)
    ax[2].scatter(evla, evdepth)
 
    #brMax = np.zeros((len(os.listdir(os.path.join(path, 'Detailed_Solutions'))), 5)) # list with max br 

    # get the maximum positions
    #count = 0
    #for dir_ in os.listdir(os.path.join(path, 'Detailed_Solutions')):
    #    br = np.load(os.path.join(path, 'Detailed_Solutions', dir_, 'out_Max.npy'))
    #    brMax[count, :] = br[np.argmax(br[:,0])]
    #    count+=1

    # max dist
    #dist_l = []
    #for brspot in brMax:
    #    dist = gps2dist_azimuth(evla, evlo, brspot[2], brspot[1])[0]/1000
    #    dist_l.append(dist)

    plt.savefig(os.path.join(out_path, filename+'.'+fileformat), dpi=dpi)

    return
