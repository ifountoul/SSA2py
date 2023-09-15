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

import numpy as np

# Obspy Imports
###############

import os
from obspy.geodetics.base import gps2dist_azimuth

# Local Imports
###############

from SSA2py.core import config
from SSA2py.core.basic_f.other import createDir, brFiles
from SSA2py.core.plotting_functions.Atlas import atlas
from SSA2py.core.plotting_functions.Animation import _animation_
from SSA2py.core.plotting_functions.MaxBrightTimeStep import MaxBrightTimeStep_
from SSA2py.core.plotting_functions.MaxBrightTimeStep_2 import MaxBrightTimeStep_2_
from SSA2py.core.plotting_functions.RecordsSection import recordSection
from SSA2py.core.plotting_functions.Uncertainty_Analysis import _heatmap_statistics, _heatmap_statistics_ver, max_br_analysis
from SSA2py.core.plotting_functions.ARF import ARFplots
from SSA2py.core.plotting_functions.weights_plots import plot_final_weights
from SSA2py.core.plotting_functions.MaxBrightAllTimeSteps import plotMaxBrightAllTimeSteps

def plot_res_(paths, out_path):
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

        try:

             if maxdist<=60:
                 rings = 25
                 rings_min = 25
             else:
                 rings = 50
                 rings_min = 50

             atlas(config.inv, stations_used, stations_notused,\
                   config.org.latitude, config.org.longitude, config.org.depth/1000,\
                   extent_lon=0.2, extent_lat=0.2, extent_lon_inset=8.0, extent_lat_inset=8.0, towns=False,\
                   rings_min=rings_min, rings_max=int(maxdist)+50, rings_step=rings, hypo_lines=True,\
                   meridians=True, plates=True, filepath=out_path, filename='atlas', fileformat='png', dpi=400)

        except Exception as e:
             config.logger.error(e)
             config.logger.info('Cannot Plot Atlas Map.')        

        config.logger.info('Building Maximum Brightness Per Timestep Map...')

        if config.gridRules[0][0]=='box':
            maxgrid=max(config.gridRules[0][1], config.gridRules[0][2])
        else:
            maxgrid=max([gps2dist_azimuth(config.org.latitude, config.org.longitude, g[1], g[0])[0]/1000 for g in config.grid])


        try:
            MaxBrightTimeStep_(paths[0], [], config.org.latitude, config.org.longitude, config.org.depth/1000, config.org.time,\
                               config.inv, stations_used, startTime=None, endTime=None, minBrig=None, maxBrig=None,\
                               min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                               points_size=8, maxgrid=maxgrid,\
                               faults=True, grid=True, hypo=True, colormap='rainbow', topo=True,\
                               meridian=True, info_box=True, Test='MAIN',autoselect=False,\
                               filename='MaximumBrightness', outpath=out_path, fileformat='png', dpi=400)
        except Exception as e:
             config.logger.error(e)
             config.logger.info('Cannot Plot Maximum Brightness.')

        config.logger.info('Building Maximum Brightness at All Timesteps...')
        try:
            plotMaxBrightAllTimeSteps(paths[0], config.org.latitude, config.org.longitude, config.org.depth/1000,\
                                      min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                                      hypo=True, colormap='viridis', mincolor=None, maxcolor=None,\
                                      filename='MaximumBrightnessAllTimesteps', outpath=out_path, fileformat='png', dpi=400) 
        except Exception as e:
            config.logger.error(e)
            config.logger.info('Cannot Plot Maximum Brightness at All Timesteps.')

        config.logger.info('Building Maximum Brightness Per Timerange...')

        try:
        
            MaxBrightTimeStep_2_(paths[0],  config.org.latitude, config.org.longitude, config.org.depth/1000,\
                                 min_lon=None, min_lat=None, max_lon=None,\
                                 max_lat=None, usemask=False, rolling_window=0.5, mask_value=0.9,\
                                 hypo=True, colormap='viridis', mincolor=None, maxcolor=None,filename='MaximumBrightnessPerTimeRange',\
                                 outpath=out_path, fileformat='png', dpi=400)
        except:
            config.logger.error(e)
            config.logger.info('Cannot Plot Maximum Brightness Per Timerange.')

        config.logger.info('Building Records Section...')

        recordSection(config.st.copy(), config.org.time, time_min_=None, time_max_=None, dist_min=None,\
                  dist_max=None, scale=1.0, labels=True, grid=True,\
                  filename='RecordSection', outpath=out_path, fileformat='png', dpi=400)
    
    if config.cfg['Plotting']['Animation'] == True:
        config.logger.info('Building Animation...')
        try:
            _animation_(paths[0], config.org.latitude, config.org.longitude, config.org.depth/1000,\
                        min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                        hypo=True, colormap='viridis', mincolor=None, maxcolor=None, outpath=out_path, filename='animation', fileformat='mp4')
        except:
            config.logger.error(e)
            config.logger.info('Cannot Plot Animation.')

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

    if config.cfg['Plotting']['Plots'] == True:

        config.logger.info('Building Array Response Function Plots...')

        try:
            files_data = brFiles(paths[0])
            ARFplots(files_data, config.inv, config.org.latitude, config.org.longitude,\
                     config.org.depth/1000, config.st.copy(), config.org.time, tt, filename='ARF',\
                     outpath=out_path, fileformat='png', dpi=600)  
        except Exception as e:
            config.logger.error(e)
            config.logger.info('Cannot Plot Array Response Function Plots.')

    if config.cfg['Plotting']['Animation'] == True:
        
        config.logger.info('Building Array Response Function Animation...')
   
        try:
            _animation_(paths[0], config.org.latitude, config.org.longitude, config.org.depth/1000,\
                        min_lon=None, min_lat=None, max_lon=None, max_lat=None, min_depth=None, max_depth=None,\
                        hypo=True, colormap='viridis', mincolor=None, maxcolor=None, outpath=out_path, filename='animation', fileformat='mp4')

        except Exception as e:
            config.logger.error(e)
            config.logger.info('Cannot Built Array Response Function Animation.')


    return

def plot_Un(paths, out_path):
    """
    Plot Uncertainty Analysis Results

    """

    #Create the plot directory if does not exists
    createDir(out_path)

    # maximum brightness map
    if config.cfg['Plotting']['Plots'] is True:
        config.logger.info('Building Plots for Uncertainty Analysis Results.')

        try:
            # Variance

            variance = np.load(os.path.join(paths[0], 'Main_Solution', 'variance.npy'))
            grid = np.load(os.path.join(paths[1], 'grid.npy'))

        
            _heatmap_statistics(variance, grid, type_='Variance', colormap='viridis', outpath=out_path,\
                                filename='Heatmap_Horizontal_Variance', fileformat='png', dpi=600)
        
            _heatmap_statistics_ver(variance, grid, type_='Variance', colormap='viridis', outpath=out_path,\
                                    filename='Heatmap_Vertical_Variance', fileformat='png', dpi=600)

            # Confidence
            cf = np.load(os.path.join(paths[0], 'Main_Solution', 'confidence.npy'))

            _heatmap_statistics(cf, grid, type_='Confidence Margins', colormap='viridis', outpath=out_path,\
                                filename='Heatmap_Horizontal_Confidence', fileformat='png', dpi=600)

            _heatmap_statistics_ver(cf, grid, type_='Confidence Margins', colormap='viridis', outpath=out_path,\
                                    filename='Heatmap_Vertical_Confidence', fileformat='png', dpi=600)

            # SE
            se = np.load(os.path.join(paths[0], 'Main_Solution', 'standard_error.npy'))

            _heatmap_statistics(se, grid, type_='Standard Error', colormap='viridis', outpath=out_path,\
                                filename='Heatmap_Horizontal_Standard_Error', fileformat='png', dpi=600)

            _heatmap_statistics_ver(variance, grid, type_='Standard Error', colormap='viridis', outpath=out_path,\
                                    filename='Heatmap_Vertical_Standard_Error', fileformat='png', dpi=600)

            # Plot maximum brightness plot
            max_br_analysis(paths, outpath=out_path, filename='Max_Brightness_Analysis', fileformat='png', dpi=600)

 
        except Exception as e:
            config.logger.error(e)
            config.logger.info('Cannot Built Uncertainty Analysis Plots.')

    return
