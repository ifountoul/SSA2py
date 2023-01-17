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
from obspy.core import UTCDateTime

from SSA2py.core import config
from SSA2py.core.modules.time import arrival
from SSA2py.core.modules.normalize import normalize

def mute_picks(st, f, phase):
    """
    Mute traces from manual picks

    Arguments:
    ----------
    st: Obspy stream object
        Traces
    f: str
        Path of the file with picks
    phase: str
        Phase to BP

    Returns:
    -------- 
    st: Obspy stream object
        Muted traces

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
            # get only the S phase
            S = UTCDateTime(line.split()[2])
            try:
                # select the trace
                trace = st.select(network=net, station=sta)[0]
                # remove from stream
                st.remove(trace)
                # arrival as sample in trace?
                sample = int((S-trace.stats.starttime) * trace.stats.sampling_rate)
                if phase=='P':
                    trace.data[sample:]=0
                if phase=='S':
                    trace.data[:sample+1]=0
                # add back to stream
                st.append(trace)
            except:
                pass
        return st 
    else:
        config.logger.warning('Error in Mute procedure')
        config.logger.warning('Problem with the input file')    
        return st

def mute_traveltimes(st, stations, phase):
    """
    Mute traces from traveltimes

    Arguments:
    ----------
    st: Obspy stream object
        Traces
    stations: dict
        Stations information
    phase: str
        Phase to BP
 
    Returns:
    -------- 
    st: Obspy stream object
        Muted traces

    """

    for i in range(st.count()):
        arr = arrival(config.org.latitude, config.org.longitude, config.org.depth,\
                      stations[st[i].stats.station][3], stations[st[i].stats.station][2],\
                      stations[st[i].stats.station][4], 'S')
        sample = int(((config.org.time + arr) - st[i].stats.starttime) * st[i].stats.sampling_rate)
        if phase=='P':
            st[i].data[sample:]=0
        if phase=='S':
            st[i].data[:sample+1]=0
    return st


def mute_traces(st, stations):
    """
    Mute traces based on traveltimes or picked arrivals

    Arguments:
    ----------
    st: Obspy stream object
        Traces
    stations: dict
        Stations information

    Returns:
    --------
    st: Obspy stream object
        Muted traces
 
    """

    if int(config.cfg['Backprojection']['Settings']['Mute'][1])==1:
        config.logger.info('Mute waveforms')
        # mute based on given arrivals
        st = mute_picks(st, config.cfg['Backprojection']['Settings']['Mute'][2],\
                        config.cfg['Backprojection']['Settings']['Phase'][0]) 
    if int(config.cfg['Backprojection']['Settings']['Mute'][1])==0:
        config.logger.info('Mute waveforms')
        # mute based on traveltimes
        st = mute_traveltimes(st, stations, config.cfg['Backprojection']['Settings']['Phase'][0]) 

    # re-normalize the traces
    for tr in st:
        tr = normalize(tr, n_root_factor = config.cfg['Streams']['Normalize'][1], \
                       norm_type = config.cfg['Streams']['Normalize'][2])
    return st
