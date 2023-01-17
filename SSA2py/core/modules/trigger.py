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


import obspy, math, sys, os, inspect
import numpy as np
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.filter import envelope


"""
Collection of python functions for P waves picking.
Simple methods to identify the primary wave arrival

So far
--> STA/LTA
--> Kurtosis

"""

def trigSTALTA(trace, STA, LTA, threshold):
    """
    Simple STA/LTA trigger Algorithm

    Arguments:
    ------
    trace: Obspy trace
        Trace
    STA: float
        Short time window in seconds
    LTA: float
        Long time window in seconds
    threshold: float
        threshold for trigger
    Returns:
    ------
    start_spot: array-like
        Trigger position
    cft: array-like
        STALTA Function

    """
    #Detrend-demean the trace
    trace.detrend('demean')
    trace.detrend('linear')

    #Sampling rate
    df = trace.stats.sampling_rate

    #Filter the signal (Bandpass 2-15 Hz
    trace.filter('bandpass', freqmin=2.0, freqmax=15.0)

    #Do the STA/LTA (trace!!)
    cft = classic_sta_lta(trace.data, int(STA*df), int(LTA*df))

    #Check where in the cft the values exceed the STA_threshold
    start_spot = np.where(cft>=threshold)

    return (start_spot[0], cft)

def trigSTALTAEnv(trace, STA, LTA, threshold):
    """
    Simple STA/LTA trigger Algorithm, with envelope and smoothing (Earle and Shearer 1994)

    Arguments:
    ------
    trace: Obspy trace
        Trace
    STA: float
        Short time window in seconds
    LTA: float
        Long time window in seconds
    threshold: float
        threshold for trigger
    Returns:
    ------
    start_spot: array-like
        Trigger position
    cft: array-like
        STALTA Function
    """

    #Detrend-demean the trace
    trace.detrend('demean')
    trace.detrend('linear')

    #Sampling rate
    df = trace.stats.sampling_rate

    #Filter the signal (Bandpass 2-15 Hz)
    trace.filter('bandpass', freqmin=2.0, freqmax=15.0)

    trace = envelope(trace.data)

    #Do the STA/LTA (trace!!)
    cft = classic_sta_lta(trace, int(STA*df), int(LTA*df))

    #Smooth the STA/LTA
    cft = smooth1D(cft, window_len=200, window='hanning')

    start_spot = np.where(cft>=threshold)

    return (start_spot[0], cft)


def KurtosisRec(data, C):
    """
    Recursive pseudo-kurtosis 

    Arguments:
    ----------
    data: array-like
        Data of the trace
    C: float
        ratio of the time-step of the data and a chosen window length
    Returns:
    ------
    cft : array-like
        Characteristic Function

    Function Recursive pseudo-kurtosis from Waveloc
    used by Langet et al. (2014).

    """
    var_trace = np.std(data)
    mean = 0; var = 0; kurt = 0;
    cft = np.empty(len(data))

    for i in range(len(data)):
         mean = mean * C + (1-C)*data[i]
         var = C*var+(1-C)*(data[i]-mean)**2
         if var > var_trace:
             kurt = C*kurt+(1-C) * (data[i]-mean)**4/var**2
         else:
             kurt = C*kurt+(1-C)*(data[i]-mean)**4/var_trace**2
         cft[i] = kurt - 3

    cft = cft+abs(np.min(cft))
    cft = smooth1D(cft, window_len=20, window='hanning')
    return cft

def positive_derivative(cft, delta):
    """
    Takes first time derivative of a stream (iterates over all available
    traces) and keep only positive values (set negative values to zero).

    Arguments:
    ----------
    cft: array-like
        Kurtosis Characteristic Function (CF)
    delta: float
        Delta of the trace
    Returns:
    -------
    cft: array-like
        First time derivative of the CF

    Based in Waveloc

    """
    xs = cft
    dt = delta
    try:
         xtemp = np.gradient(xs, dt)
         for i in range(len(xtemp)):
                if xtemp[i] < 0:
                    xtemp[i] = 0
         xs = xtemp
    except:
         pass
    return cft


def trigKurtosis(trace, C, threshold=10):
    """
    Recursive pseudo-kurtosis trigger algorithm

    Arguments:
    ----------
    tr: obspy.trace object
        trace
    C: float
        ratio of the time-step of the data and a chosen window length
    threshold: float
        Treshold for Trigger
    Returns:
    --------    
    start_spot: array-like
        Trigger sample
    cft: array-like
        Kurtosis Function
    """

    #Detrend-demean the trace
    trace.detrend('linear')
    trace.detrend('demean')
    #Taper the trace
    trace.taper(max_percentage=0.15, type='hann')
    #Filter the trace (bandpass between 1-45 Hz)
    trace.filter('bandpass', freqmin=1, freqmax=45)

    cft = KurtosisRec(trace.data, C)

    #Calculate the Derivative
    cft = positive_derivative(cft, trace.stats.delta)

    start_spot = np.where(cft>=threshold)

    return (start_spot[0], cft)

def smooth1D(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.

    from https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    """

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[int(round((window_len/2-1))):-int(round((window_len/2)))]
