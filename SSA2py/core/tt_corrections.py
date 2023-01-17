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
from obspy.signal.cross_correlation import xcorr_max, correlate, xcorr_pick_correction
from obspy.core.stream import Stream
from obspy.core import UTCDateTime


# import local libraries
from SSA2py.core import config
from SSA2py.core.modules.time import arrival
from SSA2py.core.modules.trigger import trigSTALTA, trigKurtosis
from SSA2py.core.modules.snr import SNRampl

import matplotlib.pyplot as plt

def corr_picks(st, f, phase):
    """
    Correct traces using picked arrivals.

    Arguments:
    ---------
    st: Obspy stream Object
        Traces
    f: str
        Path for the file with picks
    phase: str
        P or S (phase to correct)

    Returns:
    -------
    st: Obspy stream Object
        Corrected Traces

    """
    # check that the path exists and is file
    if os.path.isfile(f) and os.path.exists(f):
        # read the file
        try:
            with open(f) as f_:
                lines = f_.readlines()
        except IOError:
            config.logger.warning("File not accessible")
            return st
        for line in lines:
            net = line.split()[0].split('.')[0]
            sta = line.split()[0].split('.')[1]
            if phase=='P':
                arr = UTCDateTime(line.split()[1])
            if phase=='S':
                arr = UTCDateTime(line.split()[2])
 
            # get the trace
            try:
                trace = st.select(network=net, station=sta)[0]
                # remove from stream
                st.remove(trace)
                # get the synthetic arrival     
                synth_arr = arrival(config.org.latitude,\
                                    config.org.longitude,\
                                    config.org.depth,\
                                    config.inv.select(network = trace.stats.network, station = trace.stats.station).copy()[-1][0].latitude,\
                                    config.inv.select(network = trace.stats.network, station = trace.stats.station).copy()[-1][0].longitude,\
                                    config.inv.select(network = trace.stats.network, station = trace.stats.station).copy()[-1][0].elevation, phase) 
                # to UTC
                synth_arr = config.org.time + synth_arr
                # correction (the observed arrival time of the phase minus the calculated traveltime from the hypocentre)
                corr = synth_arr - arr

                config.logger.info("{} {} {} {}".format("Station", trace.stats.network + '.' + trace.stats.station + '.' + trace.stats.channel,\
                                   "correction:", str(corr)+' sec'))

                # apply the correction to the trace
                trace.stats.starttime = trace.stats.starttime + corr


                # trim it again
                trace.trim(UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][0]),\
                           trace.stats.endtime, pad=True, nearest_sample=True, fill_value=trace.data[0])
                trace.trim(trace.stats.starttime, UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][1]),\
                          pad=True, nearest_sample=True, fill_value=trace.data[-1])
                # add back to stream
                st.append(trace)
            except:
                pass


        return st
    else:
        config.logger.warning('Error in Correction procedure')
        config.logger.warning('Problem with the input file')
        return st   

def corr_(st, f):
    """
    Correct traces just using a simple shift

    Arguments:
    ----------
    st: Obspy stream Object
        Traces
    f: str
        Path for the file

    Returns:
    -------- 
    st: Obspy stream object
        Corrected traces

    """
 
    # check that the path exists and is file
    if os.path.isfile(f) and os.path.exists(f):
        # read the file
        try:
            with open(f) as f_:
                lines = f_.readlines()
        except IOError:
            config.logger.warning("File not accessible")
            return st
        for line in lines:
            net = line.split()[0].split('.')[0]
            sta = line.split()[0].split('.')[1]
            shift = float(line.split()[1])
 
            try: 
                trace = st.select(network=net, station=sta)[0]

                # remove from stream
                st.remove(trace)

                config.logger.info("{} {} {}".format(trace.stats.network + '.' + trace.stats.station + '.' + trace.stats.channel,\
                                   "shift:", str(shift)+' sec'))

                # correct
                trace.stats.starttime = trace.stats.starttime + shift
              
                # trim it again
                trace.trim(UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][0]),\
                           trace.stats.endtime, pad=True, nearest_sample=True, fill_value=trace.data[0])
                trace.trim(trace.stats.starttime, UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][1]),\
                           pad=True, nearest_sample=True, fill_value=trace.data[-1])

                # add back to stream
                st.append(trace)

            except:
                pass
    else:
        config.logger.warning('Error in Correction procedure')
        config.logger.warning('Problem with the input file')
    return st    

def corr_cc(st, win):
    """
    
    Arguments:
    ----------
    st: Obspy stream Object
        Traces
    win: float
        Kurtosis window
    
    Returns:
    --------
    st: Obspy stream Object
        Corrected traces


    """

    st.detrend('demean')
    st.detrend('linear')
    st.resample(sampling_rate=100)

    # get triggers for these traces
    Ptrig_arrs = [trigKurtosis(tr.copy(), 1-(tr.stats.delta/win), threshold=50)[0] for tr in st]

    # for these triggers get SNR
    SNR = np.zeros((1, len(Ptrig_arrs)))
    for i in range(len(Ptrig_arrs)):
        # no trigger
        if len(Ptrig_arrs[i]) == 0:
            SNR[0, i] = 0
        else:
            SNR[0, i] = SNRampl(st[i].copy(), Ptrig_arrs[i][0], a_time=5)
    # get P theoretical arrival
    Ptheo_arrs = [arrival(config.org.latitude,\
                  config.org.longitude,\
                  config.org.depth,\
                  config.inv.select(network = tr.stats.network, station = tr.stats.station)[-1][0].latitude,\
                  config.inv.select(network = tr.stats.network, station = tr.stats.station)[-1][0].longitude,\
                  config.inv.select(network = tr.stats.network, station = tr.stats.station)[-1][0].elevation, 'P') for tr in st]        
    # get P theo - P trigger arrival
    Pdiff = np.zeros((1, len(Ptrig_arrs)))
    for i in range(len(Ptrig_arrs)):
        # no trigger
        if len(Ptrig_arrs[i]) == 0:
            Pdiff[0,i] = 10**8
        else:
            Pdiff[0,i] = abs((config.org.time+Ptheo_arrs[i]) - \
                              (st[i].stats.starttime + (Ptrig_arrs[i][0]/st[i].stats.sampling_rate))) 
    # decide the master station
    score = (SNR/np.max(SNR)) + (-Pdiff/np.max(Pdiff)) 
    index_master = int(np.argmax(score))    

    # cc with all the other triggered stations
    master = st[index_master].copy()

    correction_time = np.zeros((1,len(Ptrig_arrs)))
    for i in range(len(Ptrig_arrs)):
        if len(Ptrig_arrs[i]) != 0:
           ct, cc = xcorr_pick_correction((master.stats.starttime + (Ptrig_arrs[index_master][0]/master.stats.sampling_rate)), master,\
                                          (st[i].stats.starttime + (Ptrig_arrs[i][0]/st[i].stats.sampling_rate)), st[i],\
                                          t_before=0.3, t_after=0.3, cc_maxlag=0.3, filter="bandpass", filter_options={'freqmin': 1, 'freqmax': 10})
           correction_time[0,i] = ct

    # find the overall correction
    tCorr = np.zeros((1, len(Ptrig_arrs)))

    # for the reference station (observed arrival - theoretical)
    #tCorrREF = (master.stats.starttime + (Ptrig_arrs[index_master][0]/master.stats.sampling_rate)) -\
    #            (config.org.time+Ptheo_arrs[index_master])

    tCorrREF = ((master.stats.starttime + (Ptrig_arrs[index_master][0]/master.stats.sampling_rate)) - config.org.time) - Ptheo_arrs[index_master]
 
    for i in range(len(Ptrig_arrs)):
        if len(Ptrig_arrs[i]) != 0:
            tCorr = Ptheo_arrs[index_master] + tCorrREF + correction_time[0,i] - Ptheo_arrs[i]
            st[i].stats.starttime = st[i].stats.starttime + tCorr  
            
    st.normalize().plot()

    return


def tt_exec(st):
    """
    Execute the traveltime corrections.

    Arguments:
    ---------
    st: Obspy stream Object
        Traces

    Returns:
    --------
    st: Obspy stream Object
        Traces

    """

    st_ = Stream()

    if config.cfg['Streams']['Corrections'][1]==1:
        # do it per component
        for comp in ['Z', 'N', 'E', 'R', 'T', 'H']:
            # do we have traces?
            if st.select(component=comp).count()>0:
                st_ += corr_picks(st.select(component=comp).copy(), config.cfg['Streams']['Corrections'][2], 'P')
    if config.cfg['Streams']['Corrections'][1]==2:
        for comp in ['Z', 'N', 'E', 'R', 'T', 'H']:
            if st.select(component=comp).count()>0:
                st_ += corr_picks(st.select(component=comp).copy(), config.cfg['Streams']['Corrections'][2], 'S')      
    if config.cfg['Streams']['Corrections'][1]==3:
        for comp in ['Z', 'N', 'E', 'R', 'T', 'H']:
            if st.select(component=comp).count()>0:
                st_ += corr_(st.select(component=comp).copy(), config.cfg['Streams']['Corrections'][2])
    if config.cfg['Streams']['Corrections'][1]==0:
        for comp in ['Z', 'N', 'E', 'R', 'T', 'H']:
            if st.select(component=comp).count()>0:
                st_ += corr_cc(st.select(component=comp).copy(), 5)
    return st_


















