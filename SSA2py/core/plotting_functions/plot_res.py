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

import obspy, os

# local functions
from SSA2py.core import config
from SSA2py.core.basic_f.other import createDir
from SSA2py.core.plotting_functions.Atlas import atlas
from SSA2py.core.plotting_functions.MaxBrightTimeStep import MaxBrightTimeStep_
from SSA2py.core.plotting_functions.RecordsSection import recordSection
from SSA2py.core.plotting_functions.Animation import brFiles, brightAnimation
from SSA2py.core.plotting_functions.Uncertainty_Analysis import Max_Bright_Uncertainty, Max_Bright_Map
from SSA2py.core.plotting_functions.RecordswithBr import wf_
from SSA2py.core.plotting_functions.ARF import ARFplots

def plot_res_(paths, out_path, Test='MAIN', error_type=None):
    """
    Plot Results

    Arguments:
    ----------
    paths: str
        Input data paths
    out_path: str
        Output path to put plots

    """

    #Create the plot directory if does not exists
    createDir(out_path)

    if config.cfg['Plotting']['Plots'] is True:

        config.logger.info('Building Stations Map...')

        # stations
        stations_used = [tr.stats.network+'.'+tr.stats.station for tr in config.st]
        stations_notused = [i.code+'.'+i[0].code for i in config.inv if i.code+'.'+i[0].code not in stations_used]
        maxdist = max([tr.stats.distance for tr in config.st])

        atlas(config.inv, stations_used, stations_notused,\
              config.org.latitude, config.org.longitude, config.org.depth/1000,\
              extent_lon=0.2, extent_lat=0.2, extent_lon_inset=10.0, extent_lat_inset=10.0, towns=False,\
              rings_min=50, rings_max=int(maxdist)+50, rings_step=50, hypo_lines=True,\
              meridians=True, plates=True, filepath=out_path, filename='atlas', fileformat='png', dpi=400)
        
        config.logger.info('Building Maximum Brightness Per Time Step Map...')

        MaxBrightTimeStep_(paths[0], [], config.org.latitude, config.org.longitude, config.org.depth/1000, config.org.time,\
                          config.inv, stations_used, startTime=0, endTime=25, minBrig=None, maxBrig=None,\
                          min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                          points_size=9, maxgrid=max(config.gridRules[0][1], config.gridRules[0][2]),\
                          faults=True, grid=True, hypo=True, colormap='rainbow', topo=True,\
                          meridian=True, info_box=True, Test='MAIN',autoselect=False,\
                          filename='MaximumBrightness', outpath=out_path, fileformat='png', dpi=400)

        config.logger.info('Building Records Section...')

        recordSection(config.st.copy(), config.org.time, time_min_=None, time_max_=None, dist_min=None,\
                  dist_max=None, scale=1.0, labels=True, grid=True,\
                  filename='RecordSection', outpath=out_path, fileformat='png', dpi=400)
    return

def plot_res_ARF(paths, out_path, tt):
    """
    Plot Results

    Arguments:
    ----------
    paths: str
        Input data paths
    out_path: str
        Output path to put plots
    tt: numpy array
        Synthetic arrivals
    """

    #Create the plot directory if does not exists
    createDir(out_path)

    if config.cfg['Plotting']['Plots'] is True:

        config.logger.info('Building Array Response Function Plots...')

        files_data = brFiles(paths[0])
        ARFplots(files_data, config.inv, config.org.latitude, config.org.longitude,\
                 config.org.depth/1000, config.st.copy(), config.org.time, tt, filename='ARF',\
                 outpath=out_path, fileformat='png', dpi=600)  

    if config.cfg['Plotting']['Animation'] is True:
        config.logger.info('Building Maximum Brightness Animation...')
        files_data = brFiles(paths[0])
        brightAnimation(files_data, os.path.join(out_path,'animation.mp4'))

        return

def plot_Boot(paths, out_path, Test='MAIN', error_type=None):
    """
    Plot Bootstrap or Jackknife Results

    """

    #Create the plot directory if does not exists
    createDir(out_path)

    # maximum brightness map
    if config.cfg['Plotting']['Plots'] is True:
        #config.logger.info('Building Maximum Brightness Per Time Step Map for the resampling results...')

        #MaxBrightTimeStep_(paths[1], paths[0], config.org.latitude, config.org.longitude, config.org.depth/1000, config.org.time,\
        #                   startTime=None, endTime=None, minBrig=None, maxBrig=None,\
        #                   min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
        #                   maxgrid=max(config.gridRules[0][1], config.gridRules[0][2]),\
        #                   points_size=10, faults=True, grid=True, hypo=True, colormap='plasma', topo=True,\
        #                   meridian=True, info_box=True, Test=Test, error_type=error_type,\
        #                   filename='MaximumBrightness', outpath=out_path, fileformat='pdf', dpi=400) 

        config.logger.info('Building Maximum Brightness Uncertainty Analysis...')
        Max_Bright_Uncertainty(paths[0], out_path, 'MaximumBrightness_Uncertainty_Analysis', fileformat='png', dpi=400) 
         
        config.logger.info('Building Maximum Brightness Uncertainty Map...') 
        Max_Bright_Map(paths[0], out_path, config.org.latitude, config.org.longitude, config.org.depth/1000, 'MaximumBrightness_Uncertainty_Map', fileformat='png', dpi=400)


    return
 

