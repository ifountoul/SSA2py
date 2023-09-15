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

# Local Imports
###############

from SSA2py.core import config
from SSA2py.core.modules.snr import SNRampl
from SSA2py.core.modules.trigger import trigSTALTA, trigKurtosis
from SSA2py.core.modules.clip import clipAutoDetect
from SSA2py.core.modules._time_ import _arrival_


def checkSNR(trace, Qdec, triggerType, SNR_Window, SNR_thres):
    """
    Check SNR function. Apply SNR
    
    Arguments:
    ---------
    trace: Obspy trace Object
    Qdec: dict
        Dictionary to save results
    triggerType: str
        Trigger Method KURT (Kurtosis) or STALTA
    SNR_Window: float
        Window in seconds to calculate SNR
    SNR_thres: float
        SNR threshold (usually 5)

    Returns:
    -------
    Qdec: dict 
        Updated dictionary

    """

    if triggerType == 'STALTA':
        trigs = trigSTALTA(trace.copy(), STA = 1, LTA = 20)
    else:
        Kurtosis_window = 1  
        trigs = trigKurtosis(trace, 1-(trace.stats.delta/Kurtosis_window))    

    try:
        SNR = SNRampl(trace.copy(), trigs[0], a_time=SNR_Window)

        if SNR >= SNR_thres: #SNR threshold
              Qdec['SNR'] = [True, trigs[1], trigs[0], SNR] #[True or False, CF, Arrival, SNR value]
        else:
              Qdec['SNR'] = [False, trigs[1], trigs[0], SNR]
    except:
              Qdec['SNR'] = [False]

    return Qdec
    
def checkCLIP(trace, Qdec):
   """
   Check CLIP function.

   Arguments:
   ----------
   trace: Obspy trace Object
       Trace
   Qdec: dict
        Dictionary to save results

   Returns:
   -------
   Qdec: dict
        Dictionary to save results

   """
   try:
       r = clipAutoDetect(trace.copy())    
       if r[0] == False:
            Qdec['CLIP'] = [False, r[1]]
       else:
            Qdec['CLIP'] = [True, r[1]]
   except:
       Qdec['CLIP'] = [False]
   return Qdec

def checkTIME(trace, Qdec, triggerType, P_win):
    """
    Check TIME function
   
    Arguments:
    ---------
    trace: Obspy trace Object
    Qdec: dict
        Dictionary to save results
    triggerType: str
        Trigger Method KURT (Kurtosis) or STALTA
    P_win: float
        P deviation window
 
    Returns:
    --------    
    Qdec: dict
        Dictionary to save results

    """ 

    #Get P triggers
    try:
        if triggerType =='STALTA':
            trigs = trigSTALTA(trace.copy(), STA = 1, LTA = 20)
        else: 
            Kurtosis_window = 2
            trigs = trigKurtosis(trace, 1-(trace.stats.delta/Kurtosis_window))

        #No triggers
        if len(trigs)==0:
            Qdec['TIME'] = [False]
        else: #We have trigger
            #Get the P theoretical arrival
            arr = _arrival_(config.org.latitude,\
                            config.org.longitude,\
                            config.org.depth,\
                            config.inv.select(station=trace.stats.station, channel=trace.stats.channel).copy()[0][0].latitude,\
                            config.inv.select(station=trace.stats.station, channel=trace.stats.channel).copy()[0][0].longitude,\
                            config.inv.select(station=trace.stats.station, channel=trace.stats.channel).copy()[0][0].elevation,\
                            'P')
            theo_arr = config.org.time + arr
            obs_arr = trace.stats.starttime +\
                      (trigs[0]/trace.stats.sampling_rate) 
            if abs(theo_arr-obs_arr)<=P_win:
                Qdec['TIME'] = [True, obs_arr, theo_arr, P_win]
            else:
                Qdec['TIME'] = [False, obs_arr, theo_arr, P_win]
    except:
        Qdec['TIME'] = [False]
    return Qdec
