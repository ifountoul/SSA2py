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


import os
import numpy as np
from scipy import signal

# Obspy Imports
################

from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.trace import Trace
from obspy.core.stream import Stream

#Local Imports
##############

from SSA2py.core import config
from SSA2py.core.basic_f.get_tt import get_tt1d, get_tt3d
from SSA2py.core.basic_f.other import createDir
from SSA2py.core import backprojection
from SSA2py.core.grid import grid_box3D
from SSA2py.core.plotting_functions.plot_res import plot_res_ARF

# End of Imports
################



def ARF():
    """
    - Calculate the Array Response Function for the given station geometry.
    - This Method is based on one second Ricker pulse, calculated for each station. The sytnthetic traces are backprojected.
    
    --> Main goal of this function to identify the swimming artifact from the station geometry.

    """    
 
    config.logger.info('-----------------------------------')
    config.logger.info('Array Response Function calculation')
    config.logger.info('-----------------------------------')

    # calculate the Synthetic Pulses
    # calculate distances between stations and hypocenter and finally the tt for the 1D case

    if config.model['type']=='1D':
        # identify the closer grid point to the source
        idx = np.where(abs(config.grid[:,0]-config.org.longitude)==np.min(abs(config.grid[:,0]-config.org.longitude)))[0]
        idy = np.where(abs(config.grid[:,1]-config.org.latitude)==np.min(abs(config.grid[:,1]-config.org.latitude)))[0]
        idz = np.where(abs(config.grid[:,2]-config.org.depth/1000)==np.min(abs(config.grid[:,2]-config.org.depth/1000)))[0]

        # source point
        sou_ = np.intersect1d(np.intersect1d(idx, idy), idz)[0]

        dist = [gps2dist_azimuth(config.stations[sta][3], config.stations[sta][2],\
                config.org.latitude, config.org.longitude)[0]/1000 for sta in config.stations]

        elevationArray = np.arange(0, float(config.model['Elevation']) +\
                                   config.model['Granularity'],\
                                   config.model['Granularity'])

        Scoords = np.array([v[2:5] for k,v in config.stations.items()])
  
        # for each station identify the 1D table position
        tablpos = np.array([np.where(abs(elevationArray-(sel/1000))==np.min(abs(elevationArray-(sel/1000))))[0][0]\
                            for sel in Scoords[:,2]])

        tt = get_tt1d(np.array([[config.org.longitude, config.org.latitude, config.org.depth/1000]]), config.tables, tablpos,\
                      np.array([dist]), config.model['Granularity'], Scoords[:,2])[0]


    if config.model['type']=='3D':
        # identify the closer grid point to the source
        idx = np.where(abs(config.grid[:,0]-config.org.longitude)==np.min(abs(config.grid[:,0]-config.org.longitude)))[0]
        idy = np.where(abs(config.grid[:,1]-config.org.latitude)==np.min(abs(config.grid[:,1]-config.org.latitude)))[0]
        idz = np.where(abs(config.grid[:,2]-config.org.depth/1000)==np.min(abs(config.grid[:,2]-config.org.depth/1000)))[0]

        # source point
        sou_ = np.intersect1d(np.intersect1d(idx, idy), idz)[0]

        # get the paths of the tables
        tables = [i for i in config.tables if i.split('/')[-1].split('_')[0] in [*config.stations]]

        # get the position of this point in the tt tables
        grid_ = grid_box3D(1, 1, config.org.depth/1000, config.org.depth/1000, config.model['Granularity'],\
                [config.org.latitude, config.org.longitude], 200)
 
        posx = grid_[3]; posy = grid_[4];
 
        # calculate tt
        tt = get_tt3d(np.reshape(config.grid[sou_], (1,3)), config.stations, tables,\
                      config.model['Granularity'], posx, posy, np.array([config.org.depth/1000]),\
                      np.array([config.org.depth/1000]))[0]
    
    config.st = calculatePulses(config.stations, tt)

    # ARF dir
    ARF_plots = createDir(os.path.join(config.eventdir, 'Results',\
                         'ARF', 'pulses_' + config.comp, 'Plots'))
    ARF_sol = createDir(os.path.join(config.eventdir, 'Results',\
                        'ARF', 'pulses_' + config.comp, 'Detailed_Solution'))

    config.st.write(os.path.join(config.eventdir, 'Results',\
                   'ARF', 'pulses_' + config.comp + '.mseed'), format='MSEED')

    # BP
    ba = backprojection.backprojection(config.st, config.stations, ARF_sol)

    if ba==True:
        # plots
        plot_res_ARF([ARF_sol], ARF_plots, tt)
        config.logger.info('End of Array Response Function calculation')


def calculatePulses(stations, tt):
    """
    Calculate synthetic Ricker pulses for the Array function test
    
    Arguments:
    ------
    station: dict
        Dictionary with stations information
    tt: array-like
        Traveltimes for each station
    
    Returns:
    ------
    st: Obspy stream object with pulses

    """

    length = 1; dt = 0.01

    t = np.linspace(-length/2, (length-dt)/2, int(length/dt))
    i, q, e = signal.gausspulse(t, fc=5, retquad=True, retenv=True) #1 second pulse

    st = Stream()

    for sta in stations:
        tr = Trace(data=e)
        tr.stats.sampling_rate = 100
        tr.stats.network = stations[sta][0]
        tr.stats.station = stations[sta][1]
        tr.stats.channel = config.comp
        tr.stats.stalon = stations[sta][2]
        tr.stats.stalat = stations[sta][3]
        tr.stats.staelev = stations[sta][4]
        out = gps2dist_azimuth(stations[sta][3], stations[sta][2], config.org.latitude, config.org.longitude)
        tr.stats.distance = out[0]/1000
        tr.stats.azim = out[2]
        tr.stats.starttime = config.org.time + tt[list(stations.keys()).index(sta)]
        tr.trim(starttime=config.org.time + float(config.cfg['Streams']['Duration'][0]), endtime= config.org.time +\
                float(config.cfg['Streams']['Duration'][1]), pad=True, nearest_sample=True, fill_value=0)
        st.append(tr)
    return st
