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

# Imports
##########

import math, os, gc, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines


# Obspy Imports
##############
from obspy.core.stream import Stream
from obspy.geodetics.base import gps2dist_azimuth

# Local Imports
################

from SSA2py.core import config

matplotlib.use('agg')

def SNRPlot(path, streams, traces, res):
    """
    Plot SNR results

    """
   
    # normalize the stream
    traces.normalize()

    #Count the number of traces (maximum number per plot 20)
    if len(streams)<=24:
        num_p = 1
    else:
        num_p = math.ceil(len(streams)/24)

    for p in range(num_p):
    
        #Get the traces
        try:

           streams_ = streams[p*24:(p*24)+24]
           res_ = res[p*24:(p*24)+24]

        except:

           streams_ = streams[p*24:]
           res_ = res[p*24:]

        #Columns per plot
        if len(streams_)<=6:
            col_p = len(streams_)
        else:
            col_p = 6

        # Subplots are organized in a Rows x Cols Grid
        # Tot and Cols are known


        Tot = len(streams_)
        Cols = col_p

        # Compute Rows required
        Rows = math.ceil(Tot / Cols)

        # Create a Position index
        Position = range(1,Tot + 1)

        # Create main figure
        plt.close('all')

        fig = plt.figure(1, figsize=(4*col_p, 2.5*Rows))
        for k in range(Tot):
            # add every single subplot to the figure with a for loop
            ax = fig.add_subplot(Rows,Cols,Position[k])
        axes = fig.axes

        # share x axis
        for ax in range(Tot-1, Tot-(Cols+1), -1):
            axes[ax].set_xlabel("Time (s)", fontsize=10, fontweight="bold")
 
        for k in range(len(axes)):
            tr = traces.select(network= streams_[k][0],station= streams_[k][1],\
                               channel= streams_[k][3], location=streams_[k][2])[0]

            # ymin and ymax of the plot
            ymin = np.min(tr.data) - 0.2
            ymax = np.max(tr.data) + 0.2

            dur = np.linspace(0, tr.stats.npts/tr.stats.sampling_rate, tr.stats.npts)

            if len(res_[k]['SNR'])==1: 
                # plot the trace
                axes[k].plot(dur, tr.data, c="k", lw=1.5, zorder=1, alpha=0.5)
                axes[k].text(0.25, 0.5, 'Removed', weight='bold', transform=axes[k].transAxes,\
                                 fontsize=15, zorder=1)
            else:
                if res_[k]['SNR'][0]==True:
                    # plot the trace
                    axes[k].plot(dur, tr.data, c="k", lw=1.5, zorder=1)
                    # plot CF
                    axes[k].plot(dur, res_[k]['SNR'][1]/abs(np.max(res_[k]['SNR'][1])), c="r", lw=1.5, zorder=1, ls ='--')
                    # plot arrival
                    axes[k].axvline(x=res_[k]['SNR'][2]/tr.stats.sampling_rate, c='b', lw=0.8, zorder=1, ls ='--', alpha=0.7)
                    # plot span (Change this if you change the SNR time!)
                    axes[k].axvspan((res_[k]['SNR'][2]/tr.stats.sampling_rate)-5, (res_[k]['SNR'][2]/tr.stats.sampling_rate)+5,\
                                    alpha=0.3, color='gray')

                if res_[k]['SNR'][0]==False:
                    axes[k].plot(dur, tr.data, c="k", lw=1.5, zorder=1, alpha=0.5)
                    axes[k].plot(dur, res_[k]['SNR'][1]/abs(np.max(res_[k]['SNR'][1])), c="r", lw=1.5, zorder=1, ls ='--', alpha=0.5)
                    
                    axes[k].text(0.25, 0.5, 'Removed', weight='bold', transform=axes[k].transAxes,\
                                 fontsize=15, zorder=1)

            axes[k].grid(which='both', linestyle='-', linewidth=1, alpha=0.5)
            axes[k].set_xlim(xmin=dur[0], xmax=dur[-1])
            axes[k].set_ylim(ymin=ymin, ymax=ymax)
            axes[k].set_yticks([])
            axes[k].set_yticklabels([])
            axes[k].yaxis.tick_right()
            axes[k].text(0.02, 0.11, '{}.{}.{}'.format(tr.stats.network,tr.stats.station,tr.stats.channel), weight='bold', transform=axes[k].transAxes,
                         bbox=dict(boxstyle="round", fc="w", alpha=0.3), va="top",
                         ha="left", fontsize=9, zorder=2)

        fig.suptitle('SNR Test (Figure '+ str(p+1) + '/' +str(num_p) +')',\
                     fontsize=20, fontweight="bold")
        plt.figlegend([lines.Line2D([0], [0], ls='-', c='k'),
           lines.Line2D([0], [0], ls='--', c='r'), lines.Line2D([0], [0], ls='-', c='k', alpha=0.5),\
           lines.Line2D([0], [0], ls='--', c='b', alpha=0.5)],\
           ['Trace', 'CF', 'Deactivated Trace', 'Triggered Arrival'], loc="lower left",\
           mode="expand", borderaxespad=0, ncol=4, fontsize=14)
        plt.tight_layout(rect=[0,0.15,1,1])
        plt.subplots_adjust(hspace=0.5, wspace=0.2)
        plt.savefig(os.path.join(path, 'SNR'+'_'+str(p)+'.png'), dpi=400) 

        #release the memory
        plt.clf()
        plt.close()
        gc.collect()
        plt.pause(.01)

    return


def CLIPPlot(path, streams, traces, res):
    """
    Plot Clip test results

    """

    # normalize the stream
    traces.normalize()

    #from streams keep only the broadband stations
    indexes = [i for i in range(len(streams)) if streams[i][3][1]!='N']
    streams = [streams[i] for i in indexes]
    res = [res[i] for i in indexes]
    
    #Count the number of traces (maximum number per plot 50)
    if len(streams)<=24:
        num_p = 1
    else:
        num_p = math.ceil(len(streams)/24)

    for p in range(num_p):
        #Get the traces
        try:
           streams_ = streams[p*24:(p*24)+24]
           res_ = res[p*24:(p*24)+24]
        except:
           streams_ = streams[p*24:]
           res_ = res[p*24:]
        #Columns per plot
        if len(streams_)<=6:
            col_p = len(streams_)
        else:
            col_p = 6

        # Subplots are organized in a Rows x Cols Grid
        # Tot and Cols are known

        Tot = len(streams_)
        Cols = col_p

        # Compute Rows required
        Rows = math.ceil(Tot / Cols)

        # Create a Position index
        Position = range(1,Tot + 1)

        # Create main figure
        plt.close('all')

        fig = plt.figure(1, figsize=(4*col_p, 2.5*Rows))
        for k in range(Tot):
            # add every single subplot to the figure with a for loop
            ax = fig.add_subplot(Rows,Cols,Position[k])
        axes = fig.axes

        # share x axis
        for ax in range(Tot-1, Tot-(Cols+1), -1):
            axes[ax].set_xlabel("Time (s)", fontsize=10, fontweight="bold")

        for k in range(len(axes)):
            tr = traces.select(network= streams_[k][0],station= streams_[k][1],\
                               channel= streams_[k][3], location=streams_[k][2])[0]

            # ymin and ymax of the plot
            ymin = np.min(tr.data) - 0.2
            ymax = np.max(tr.data) + 0.2

            dur = np.linspace(0, tr.stats.npts/tr.stats.sampling_rate, tr.stats.npts)

            if len(res_[k]['CLIP'])==1:
                # plot the trace
                axes[k].plot(tr.data, c="k", lw=1.5, zorder=1, alpha=0.5)
            else:
                if res_[k]['CLIP'][0]==True:
                    # plot the trace
                    axes[k].plot(dur, tr.data, c="k", lw=1.5, zorder=1)
                if res_[k]['CLIP'][0]==False:
                    axes[k].plot(dur, tr.data, c="k", lw=1.5, zorder=1, alpha=0.5)
                    axes[k].plot(dur[:-1], res_[k]['CLIP'][1]/abs(np.max(res_[k]['CLIP'][1])), c="r", lw=1.5, zorder=1, ls ='--', alpha=0.5)

                    axes[k].text(0.25, 0.5, 'Removed', weight='bold', transform=axes[k].transAxes,\
                                 fontsize=15, zorder=1)

            axes[k].grid(which='both', linestyle='-', linewidth=1, alpha=0.5)
            axes[k].set_xlim(xmin=dur[0], xmax=dur[-1])
            axes[k].set_ylim(ymin=ymin, ymax=ymax)
            axes[k].set_yticks([])
            axes[k].set_yticklabels([])
            axes[k].yaxis.tick_right()
            axes[k].text(0.02, 0.13, '{}.{}.{}'.format(tr.stats.network,tr.stats.station,tr.stats.channel), weight='bold', transform=axes[k].transAxes,
                         bbox=dict(boxstyle="round", fc="w", alpha=0.3), va="top",
                         ha="left", fontsize=9, zorder=2)

        fig.suptitle('CLIP Test (Figure '+ str(p+1) + '/' +str(num_p) +')',\
                     fontsize=20, fontweight="bold")
      
        plt.figlegend([lines.Line2D([0], [0], ls='-', c='k'),
           lines.Line2D([0], [0], ls='--', c='r'), lines.Line2D([0], [0], ls='-', c='k', alpha=0.5)],
           ['Trace', 'Clipped Sample', 'Deactivated Trace'], loc="lower left",\
           mode="expand", borderaxespad=0, ncol=3, fontsize=14) 
        plt.tight_layout(rect=[0,0.15,1,1])
        plt.subplots_adjust(hspace=0.5, wspace=0.2)
        plt.savefig(os.path.join(path, 'CLIP'+'_'+str(p)+'.png'), dpi=400)

        plt.clf()
        plt.close()
        gc.collect()
        plt.pause(.01)

    del traces

    return

def TIMEPlot(path, streams, traces, res):
    """
    Time Plot. This will be presented in Records Section form

    """

    #Compute the epicentral distance to all stations
    offsets = []
    for tr in traces:
        dist = gps2dist_azimuth(config.org.latitude,config.org.longitude, config.inv.select(station=tr.stats.station,\
                                            channel=tr.stats.channel).copy()[0][0].latitude,\
                                            config.inv.select(station=tr.stats.station,\
                                            channel=tr.stats.channel).copy()[0][0].longitude)[0]/1000
        tr.stats.distance = dist 
        offsets.append(dist)
    max_offset = max(offsets)
    min_offset = min(offsets)

    #normalize trace
    traces.normalize()
    #scale traces
    scale = (max_offset - min_offset) * 1/(traces.select(component='Z').count() * 2.0) 

    for comp in ['Z', 'E', 'N']:
        st_ = Stream()

        fig,ax = plt.subplots(1,1,figsize=(5,7))
        time_min = []; time_max = []; theo_arr = []
        win = 0
   
        for i in range(len(streams)):
            if streams[i][-1][-1] == comp:
                if len(res[i]['TIME'])>1:
                    win = res[i]['TIME'][-1]
                    st_ = traces.select(station=streams[i][1], channel =streams[i][-1])
                    data = (st_[0].data * scale) + st_[0].stats.distance
                    time = np.arange(st_[0].stats.npts) * st_[0].stats.delta\
                          + (st_[0].stats.starttime - config.org.time)
                    if res[i]['TIME'][0]==False:
                        ax.plot(time, data, color ='peru', ls ='--', lw = 0.3, alpha=0.8)
                    if res[i]['TIME'][0]==True:
                        ax.plot(time, data, color ='k', lw = 0.3, alpha = 0.6)

                    Obs_arr = time[0]+(res[i]['TIME'][1] - st_[0].stats.starttime)
                    theo_arr.append([time[0]+(res[i]['TIME'][2] - st_[0].stats.starttime), st_[0].stats.distance]) 

                    ax.plot([Obs_arr, Obs_arr], [st_[0].stats.distance-(1.0*scale),\
                            st_[0].stats.distance+(1.0*scale)], color='b', ls='--', lw = 1.0)

                    time_min.append(min(time))
                    time_max.append(max(time))
                else:
                    ax.plot(time, data, color ='peru', ls ='--', lw = 1)

        #Add station labels to offset axis
        for i in range(len(streams)):
            if streams[i][-1][-1] == comp:
                tr = traces.select(station=streams[i][1], channel =streams[i][-1])[0]
                ax.text(time[-1]+1, tr.stats.distance, tr.stats.station, zorder=1,\
                        fontsize=5, fontstyle='italic', fontweight='demi')

        #Plot the theo arrivals
        theo_arr = np.array(theo_arr) 
        #sort the array
        theo_arr = theo_arr[theo_arr[:, 1].argsort()]
        ax.plot(theo_arr[:,0], theo_arr[:,1], color='maroon', ls='-', lw = 1.5, alpha=0.5)
        ax.plot(theo_arr[:,0]-win, theo_arr[:,1], color='maroon', ls='--', lw = 0.8, alpha=0.5)
        ax.plot(theo_arr[:,0]+win, theo_arr[:,1], color='maroon', ls='--', lw = 0.8, alpha=0.5)
        #################

        ax.set_xlim([min(time_min), max(time_max)])
        ax.set_xlabel('Time Relative to Origin (s)',\
                      fontweight = 'bold',\
                      fontsize = 10)
        ax.set_ylabel('Epicentral Distance (km)',\
                      fontweight = 'bold',\
                      fontsize = 10)
        ax.set_title('Component: ' + comp,\
                     fontweight = 'bold',\
                     fontsize = 12)
        ax.grid()        

        plt.figlegend([lines.Line2D([0], [0], ls='-', c='k', lw = 0.3, alpha = 0.6),
                       lines.Line2D([0], [0], color ='peru', ls ='--', lw = 0.3, alpha = 0.8),\
                       lines.Line2D([0], [0], color='b', ls='--', lw = 1.0), lines.Line2D([0], [0],\
                       color='maroon', ls='-', lw = 1.5, alpha=0.5)],\
                       ['Accepted Trace', 'Deactivated Trace', 'P Trigger', 'P Theoretical Arrival'], loc="lower left",\
                       mode="expand", borderaxespad=0, ncol=4, fontsize='x-small')
        plt.savefig(os.path.join(path, 'TIME'+'_'+ comp +'.png'), dpi=400)
        
        plt.clf()
        plt.close()
        gc.collect()
        plt.pause(.01)

    del traces
    return

