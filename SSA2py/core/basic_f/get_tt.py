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

from numba import jit, prange
import math, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

@jit(nopython=True, parallel=True)
def get_tt1d(grid, tables, tablpos, dist, gran, Sel):
    """
    Get tt 1D times given the above info

    Arguments:
    ---------
    grid: numpy array
        Array with numpy grid
    tables: list 
        List with numpy arrays (tt tables)
    tablpos: array-like
        Position of array fro each station
    dist: array-like
        Distance of grids from stations
    gran: float
        Granularity of the tables
    Sel: array-like
        Stations elevation

    Returns:
    --------
    tt: array-like
        Array with tt from each grid-station

    """


    #Pre-allocate the tt array
    tt = np.empty((len(grid), len(tablpos)), np.float32)

    for g in prange(len(grid)):
        for s in prange(len(tablpos)):

            #round the depth and the dist to the nearest multiple
            # of traveltimes granularity
            dist[g,s] = gran * np.around((dist[g,s]/gran))
            depth = gran * np.around((grid[g,2]+ (Sel[s]/1000))/gran)

            dist_p = int(dist[g,s]/gran)
            depth_p = int(depth/gran)
            tt[g,s] = tables[tablpos[s]][dist_p][depth_p]
    return tt

def get_tt3d(grid, stations, tables, gran, posx, posy, posz1, posz2):
    """
    Get the travel time info from the 3D tables

    Arguments:
    ----------
    grid: numpy array
        Array with numpy grid
    stations: dict
        Dictionary with stations info
    tables: list
        List with tables path
    gran: float
        Granularity of the grid
    posx: numpy array
        Positions to get from tables in x
    posy: numpy array
        Positions to get from tables in y
    posz1: float
        Minimum depth position from tables in z
    posz2: float
        Maximum depth position from tables in z

    Returns:
    -------
    tt: array-like
        Array with tt from each grid-station

    """   

 
    #Pre-allocate the tt array
    tt = np.empty((len(grid), len(tables)), np.float32)

    for s in range(len(tables)):
        #Read the table
        tab_ = np.load(tables[s])
        #Index in Z
        elev_ = gran * round(abs(stations[tables[s].split('/')[-1].split('_')[0]][4]/1000)/gran)

        idxZ1 = int(elev_/gran) + int((posz1/gran)) #Start Z
        idxZ2 = int(posz2/gran) + int(elev_/gran) #End Z

        tab_ = tab_[int(posx[0]):int(posx[-1]+1), int(posy[0]):int(posy[-1]+1), idxZ1:idxZ2+1].ravel()

        for g in range(len(grid)):
            tt[g,s] = tab_[g]
    return tt
