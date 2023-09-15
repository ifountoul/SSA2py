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
import os
import matplotlib.pyplot as plt

def recordSection(st, time_, time_min_=None, time_max_=None, dist_min=None,\
                  dist_max=None, scale=1.0, labels=True, grid=True,\
                  filename='RecordSection', outpath='.', fileformat='pdf', dpi=400):
    """
    Plot Records section with distance

    Arguments:
    ----------
    st: Obspy Stream Object
        Streams
    time_: UTCDatime Object
        Origin time
    time_min_: float
         Minimum Time to plot.
    time_max_: float
         Maximum Time to plot.
    dist_min: float
         Minimum distance.
    dist_max: float
         Maximum distance.
    scale: float
         Scale the traces
    labels: bool
         Plot labels?
    grid: bool
         Plot grid?
    filename: str
         Filename
    outpath: str
         Path to save
    fileformat: str
         Format of file
    dpi: int
         Dpi
    
    filepath: str
        Path and name of the plot to save
    """

    st.normalize()

    dists = [tr.stats.distance for tr in st]

    max_offset = max(dists)
    min_offset = min(dists)

    scale = (max_offset - min_offset) * 1/(st.count() * scale)

    fig,ax = plt.subplots(1,1,figsize=(5,5))
    time_min = []; time_max = []
    for tr in st:
        data = (tr.data * scale) + tr.stats.distance
        time = np.arange(tr.stats.npts) * tr.stats.delta\
                          + (tr.stats.starttime - time_)
        ax.plot(time, data, color ='k', lw = 0.3, alpha = 0.6)
        time_min.append(min(time))
        time_max.append(max(time))

    if (time_min_ is None) and (time_max_ is None):
        ax.set_xlim([min(time_min), max(time_max)])
    if (time_min_ is not None) and (time_max_ is not None):
        ax.set_xlim([time_min_, time_max_])
    if (time_min_ is not None) and (time_max_ is None):
        ax.set_xlim([time_min_, max(time_max)])
    if (time_min_ is None) and (time_max_ is not None):
        ax.set_xlim([min(time_min), time_max_])

    if (dist_min is not None) and (dist_max is not None):
        ax.set_ylim([dist_min, dist_max])
    if (dist_min is not None) and (dist_max is None):
        ax.set_ylim([dist_min, max(dists)])
    if (dist_min is None) and (dist_max is not None):
        ax.set_ylim([min(dists), dist_max])

    ax.set_xlabel('Time Relative to Origin (s)',\
                  fontweight = 'bold',\
                  fontsize = 10)
    ax.set_ylabel('Epicentral Distance (km)',\
                  fontweight = 'bold',\
                  fontsize = 10)
    if grid==True:
        ax.grid()

    #Add station labels to offset axis
    if labels==True:
        for tr in st:
            ax.text(max(time_max)+1, tr.stats.distance, tr.stats.station, zorder=1,\
                    fontsize=5, fontstyle='italic', fontweight='demi')

    plt.savefig(os.path.join(outpath, filename+'.'+fileformat), dpi=dpi)
    plt.close('all')

    return
