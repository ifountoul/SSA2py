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

import os
import numpy as np

# Obspy Imports
###############

from obspy.core import UTCDateTime


# Local Imports
###############

from SSA2py.core import config
from SSA2py.core.modules._time_ import _arrival_

def exec_corr(stations, f):
    """
    Get TT Corrections from files.

    Input:
    ------
    stations: dict
        Stations dictionary

    Output:
    ------
    cor: array
         Array with corrections

    """
    
    # No Corrections yet
    cor = np.zeros((len(stations),))

    # Stations Dictionary
    keys_list = list(stations.keys())

    # The file exists?
    if not os.path.isfile(f) or not os.path.exists(f):
        config.logger.warning("Problem with Correctios file. Continue with no time shifts.")
        return cor

    # The file exists. Open the file
    try:
        with open(f) as f_:
            lines = f_.readlines()
    except IOError:
            config.logger.warning("File not accessible")
            return cor
  
    for line in lines:
        try:
            net = line.split()[0].split('.')[0]
            sta = line.split()[0].split('.')[1]
        except Exception as e:
            config.logger.warning("Cannot read properly the correction file.")
            config.logger.warning(e)
            continue

        # Simple shifts?
        if config.cfg['Streams']['Corrections'][1]==3:
            try:
                shift = float(line.split()[1])
                index = keys_list.index(sta)
                cor[index] = shift

                config.logger.info("{} {} {} {}".format("Station", sta,\
                                   "correction:", str(cor[index])+' sec'))

            except Exception as e:
                continue

        if config.cfg['Streams']['Corrections'][1]==1:
            try:
                arr = UTCDateTime(line.split()[1])
 
                # get the synthetic arrival     
                synth_arr = _arrival_(config.org.latitude,\
                                      config.org.longitude,\
                                      config.org.depth,\
                                      stations[sta][3],\
                                      stations[sta][2],\
                                      stations[sta][4], 'P') 
                # to UTC
                synth_arr = config.org.time + synth_arr
                # correction (the observed arrival time of the phase minus the calculated traveltime from the hypocentre)
                index = keys_list.index(sta)
                cor[index] = arr - synth_arr

                config.logger.info("{} {} {} {}".format("Station", sta,\
                                   "correction:", str(round(cor[index],2))+' sec'))

            except Exception as e:
                continue

        if config.cfg['Streams']['Corrections'][1]==2:
            try:
                arr = UTCDateTime(line.split()[2])

                # get the synthetic arrival     
                synth_arr = _arrival_(config.org.latitude,\
                                      config.org.longitude,\
                                      config.org.depth,\
                                      stations[sta][3],\
                                      stations[sta][2],\
                                      stations[sta][4], 'S')                                                                                       

                # to UTC
                synth_arr = config.org.time + synth_arr

                # correction (the observed arrival time of the phase minus the calculated traveltime from the hypocentre)
                index = keys_list.index(sta)
                cor[index] = arr - synth_arr

                config.logger.info("{} {} {} {}".format("Station", sta,\
                                   "correction:", str(round(cor[index],2))+' sec'))
                
            except Exception as e:
                continue 
    return cor
