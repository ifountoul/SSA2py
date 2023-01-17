#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

###Spatial weights

def NExWeight(stations):
    """
    Negative exponential weights.

    Based on:
    ---------
    Balancing unevenly distributed data in seismic tomography: a global
    adjoint tomography example(2019). Ruan et al. (doi: 10.1093/gji/ggz356)

    Input:
    ------
    station: dict
        Dicitionary with basic stations info

    Output:
    ------
    Wexpdist: array-like
        Array with weights per station

    """

    #Weight matrix
    Wexpdist = np.zeros((len(stations),\
               len(stations)))   

    for i, (k, v) in enumerate(stations.items()):
        for j, (k_, v_) in enumerate(stations.items()):
            if i==j: 
                Wexpdist[i,j] = 1
            else:
                dist = gps2dist_azimuth(stations[k][3], stations[k][2], stations[k_][3], stations[k_][2])[0]/1000
                Wexpdist[i,j] = np.exp(-(dist))
    Wexpdist = 1/np.sum(Wexpdist,axis=0)
    #Normalize
    Wexpdist = Wexpdist/np.max(Wexpdist) 
    return Wexpdist

def SPWeights(stations, typeofweight='NEW'):
    """
    Spatial weights management

    """

    if typeofweight=='NEW':
        W = NExWeight(stations)



    return W

