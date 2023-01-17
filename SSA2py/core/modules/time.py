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


import yaml, os
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

from SSA2py.core import config, backprojection
from SSA2py.core.basic_f.get_tt import get_tt1d

def arrival(lat, lon, depth, lat_, lon_, elev, phase):
    """
    Get theoretical arrival (P or S) for specific coordinates.
    Mainly for use in the TIME tool. Isolated from the rest backprojection program.

    Arguments:
    ------
    lat: float
        Latitude of the point
    lon: float
        Longitude of the point
    depth: float
        Depth of the point
    lat_: float
        Latitude of the station
    lon_: float
        Longitude of the station
    elev: float
        Elevation of the station
    phase: str
        P or S
    Returns:
    ------
    arr: float
        Arrival time after the origin time
 
    """

    #Get the tables for the specific phase (be careful with that)
    config.phase = phase
    backprojection.getTTtables()

    
    if config.model['type']=='1D':
        dist = gps2dist_azimuth(lat, lon, lat_,\
                                lon_)[0]/1000

        elevationArray = np.arange(0, float(config.model['Elevation']) +\
                                   config.model['Granularity'],\
                                   config.model['Granularity'])

        tablpos = np.array([np.array(np.where(abs(elevationArray-(elev/1000))\
                  ==np.min(abs(elevationArray-(elev/1000))))[0][0])])
        pos = np.array([[lon, lat, depth/1000]])
        #get tt
        tt = get_tt1d(pos, config.tables, tablpos,\
                      np.array([[dist]]), config.model['Granularity'], np.array([elev/1000]))[0][0]
        return tt

    if config.model['type']=='3D':
        pass      
        return
