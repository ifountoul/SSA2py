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

import matplotlib.pyplot as plt
from matplotlib import lines
import numpy as np
import obspy, os, math
from obspy.core.stream import Stream
from scipy import interpolate
from scipy.signal import savgol_filter


# local functions
from SSA2py.core import config
from SSA2py.core.modules.trigger import smooth1D

def wf_(brpath, st, filename='Waveforms',\
        outpath='.', fileformat='pdf', dpi=400):
    """
    Plot used waveforms with distance

    Arguments:
    ---------
    brpath: str
        Brightness results.
    st: Obspy stream Object
        Traces
    filename: str
        Filename
    outpath: str
        Output path
    fileformat: str
        Format of the file.
    dpi: int
       Dpi
   
    Returns:
    -------- 
    """
    # read maximum Bright
    br = np.load(os.path.join(brpath, 'out_Max.npy'))
    #read tt
    tt = np.load(os.path.join(brpath, 'tt.npy'))
    # scanning time
    time = np.around(np.arange(config.scanningRules[0][0],\
                     config.scanningRules[0][1]+config.cfg['Backprojection']['Settings']['TimeShift'],\
                     config.cfg['Backprojection']['Settings']['TimeShift']), 2)
    # make sure time is float
    time = time.astype("float")

    #Find the time the zero or the first positive numbers index
    index = 0
    for i in range(len(time)):
        if time[i]>=0: #0
            index = i
            break
    #time = time[i:] #Time slice

    # station names
    st_names = [tr.stats.station for tr in st]

    # normalize the streams
    st.normalize()

    #Sort based on distance
    st = st.sort(['distance'])

    #Count the number of traces (maximum number per plot 50)
    if len(st)<=20:
        num_p = 1
    else:
        num_p = math.ceil(len(st)/20)

    for p in range(num_p):
        #Get the traces
        try:
            st_ = Stream(st[p*20:(p*20)+20])
        except:
            st_ = Stream(st[p*20:])

        #Columns per plot
        if len(st_)<=5:
            col_p = 1
        else:
            col_p = math.ceil(len(st_)/5)

        # Subplots are organized in a Rows x Cols Grid
        # Tot and Cols are known

        Tot = st_.count()
        Cols = col_p

        # Compute Rows required
        Rows = Tot // Cols
        Rows += Tot % Cols

        # Create a Position index
        Position = range(1,Tot + 1)

        # Create main figure
        plt.close('all')
        fig = plt.figure(1, figsize=(3.3*col_p, 2.0*Rows))
        for k in range(Tot):
            # add every single subplot to the figure with a for loop
            ax = fig.add_subplot(Rows,Cols,Position[k])
        axes = fig.axes

        # share x axis
        for ax in range(Tot-(Cols+1), -1, -1):
            axes[ax].set_xticklabels([])
            axes[ax].xaxis.set_ticks_position('none')
        for ax in range(Tot-1, Tot-(Cols+1), -1):
            axes[ax].set_xlabel("Time (s)", fontsize=10, fontweight="bold")

        x_min = []
        x_max = []
        for k in range(len(axes)):
            tr = st_[k].copy()

            # ymin and ymax of the plot
            ymin = np.min(tr.data) - 0.1
            ymax = np.max(tr.data) + 0.1

            # duration of the trace
            dur_ = np.linspace(float(config.cfg['Streams']['Duration'][0]),\
                               float(config.cfg['Streams']['Duration'][1]),\
                               len(tr.data))
            # plot the trace
            axes[k].plot(dur_, tr.data, c="k", lw=1.5, zorder=1)
            axes[k].set_ylim(ymin=ymin, ymax=ymax)
            axes[k].set_yticks([])
            axes[k].set_yticklabels([])
            axes[k].yaxis.tick_right()

            axes[k].text(0.02, 0.90, '{}.{} \n{} {} {} \n{} {} {}'.format(tr.stats.network,tr.stats.station,\
                        'Dist:', str(np.around(tr.stats.distance,1)), 'km', 'Azim:', str(np.around(tr.stats.azim,1)), '$^\circ$'),\
                        alpha=0.8, transform=axes[k].transAxes,
                        bbox=dict(boxstyle="round", fc="w", alpha=0.3), va="top",
                        ha="left", fontsize=8, zorder=2)

            # plot origin time
            c = axes[k].axvline(x=0, ymin=ymin, ymax=ymax, ls='--', lw = 0.5, color ='r', label = 'Orig. Time')

            #################################
            #Add in traces the amplitude positions
            tt_ = []; br_values = [];
            for br_ in br:
                # find in grid the max br posision
                idx = np.where((config.grid[:,0]==br_[1]))[0]
                idy = np.where((config.grid[:,1]==br_[2]))[0]
                idz = np.where((config.grid[:,2]==br_[3]))[0]
                #
                idxy = np.intersect1d(idx,idy)
                idxyz = np.intersect1d(idxy, idz)
                tt_.append(tt[idxyz, st_names.index(tr.stats.station)][0])
                br_values.append(br_[0])
            # relatively to origin
            tt_ = np.array(tt_)+br[:,-1]

            # normalize the br values 0-1
            br_values  = (br_values - np.min(br_values)) / (np.max(br_values) - np.min(br_values))

            # smooth data
            #add in the start and to the end zero
            br_values[0] = 0
            br_values[-1] = 0
 
            f = interpolate.interp1d(tt_, br_values, kind='nearest')
            xnew = np.linspace(tt_.min(), tt_.max(), len(tt_)) 
            ynew = f(xnew)

            #ynew = savgol_filter(ynew, int(len(ynew)/20), 3)
            ynew = smooth1D(ynew, window_len=11 ,window='hanning')

            #Normalize again
            ynew  = (ynew - np.min(ynew)) / (np.max(ynew) - np.min(ynew))

            p_ = axes[k].plot(xnew, ynew, color='r', alpha=0.6)
            axes[k].fill_between(xnew, -1, ynew, color='gray', alpha=0.1)

            x_min.append(tt_[0] - 3)
            x_max.append(tt_[-1] + 1)

        for k in range(len(axes)):
            axes[k].set_xlim(xmin = min(x_min)-1, xmax = max(x_max)+1)

        fig.suptitle('Processed Waveforms (Figure '+ str(p+1) + '/' +str(num_p) +')',\
                     fontsize=12, fontweight="bold")
        plt.subplots_adjust(hspace=0.0)
        plt.figlegend([lines.Line2D([0], [0], ls='-', c='k'),
           lines.Line2D([0], [0], ls='-', c='r', alpha=0.6)],\
           ['Trace', 'Backtraced Maximum Brightness',], loc="lower left",\
           mode="expand", borderaxespad=0, ncol=2)
        plt.tight_layout(rect=[0,0.02,1,1])
        plt.savefig(os.path.join(outpath, filename+str(p)+'.'+fileformat), dpi=dpi)
        plt.close('all')

    return
            

