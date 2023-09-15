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

import os, copy, random, glob, re
from scipy.special import erfinv
import numpy as np
from sklearn.utils import resample
from scipy import stats

# Obspy Imports
###############

from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import kilometer2degrees
from obspy.core.stream import Stream

# Local Imports
###############

from SSA2py.core import backprojection, config
from SSA2py.core.basic_f.other import createDir

def extract_number_from_path(path):
    match = re.search(r'out_(-?\d+\.\d+)\.npy', path)

    if match:
        return float(match.group(1))
    return -1  # Return a default value if no match is found

def do_statistics(runs_br, runs_num, main_res):
    """
    Calculate Variance, Standard Error, Bias and Confidence interval for the resampling test. 

    """

    # Variance estimation
    mean_estimation = np.mean(runs_br, axis=1)
    sub = runs_br - mean_estimation[:, np.newaxis]
    
    variance_estimate = np.sum(sub**2, axis=1) / (runs_num - 1)
    
    # Standard error estimation
    standard_error_estimate = np.sqrt(variance_estimate)

    # Bias
    bias = (runs_num - 1) * (mean_estimation - main_res)
    
    # Confidence Interval
    confidence_level = 0.95
    t_score = stats.t.ppf((1 + confidence_level) / 2, df=runs_num - 1)
    
    # Calculate the margin of error
    margin_of_error = t_score * standard_error_estimate

    return variance_estimate, standard_error_estimate, bias, margin_of_error


def Jackknife_stats(st, stations, paths):
    """
    Perform the Jackknife test to estimate various statistics


    Arguments:
    ------
    st: Obspy Object 
        Obspy Stream Object
    stations: dict
        Dictionary with stations info
    paths: list
        List with various paths
    
    Returns:
    --------
    - 

    """
   
    config.logger.info('--------------')
    config.logger.info('Jackknife Test')
    config.logger.info('--------------')

    #create dirs
    createDir(paths[0])

    # detailed solutions
    path_JN_de = createDir(os.path.join(paths[0], 'Detailed_Solutions'))

    # main solution
    path_JN_main = createDir(os.path.join(paths[0], 'Main_Solution'))

    #Delete-one resampling
    for i_ in range(st.count()):
        #Copy dict and traces
        st_ = st.copy()
        stations_ = copy.deepcopy(stations)

        #Remove station
        del stations_[st_[i_].stats.station]
        st_.pop(i_)

        #Perform BP with this Stream and stations

        #Save dir
        config.job = 'Jackknife ' + '{}/{}'.format(str(i_+1), str(st.count()))        

        Jack_path_save = os.path.join(path_JN_de, str(i_))
        backprojection.backprojection(st_, stations_, Jack_path_save)

    # Calculate the statistics
    # Get the results of the first run
    runs = sorted([int(d) for d in os.listdir(path_JN_de)])

    npy_files = [file.split('/')[-1] for file in glob.glob(os.path.join(path_JN_de,str(runs[0]), "out*.npy")) if "out_Max.npy" not in file]

    # Sorted files
    npy_files = sorted(npy_files, key=extract_number_from_path)

    # How many runs ?
    runs_num = len(os.listdir(path_JN_de))
 
    # Keep to lists
    se_res = []
    var_res = []
    bias_res = []
    cf_res = []

    # Loop for the timestep files:
    for npy_i in range(len(npy_files)):
        file_ = npy_files[npy_i]

        # In this array keep the results
        runs_br = []
    
        # Loop runs
        for d in runs:
           # Load files
           new_path = os.path.join(path_JN_de, str(d), file_)

           res_load = np.load(new_path)
           runs_br.append(res_load[:,0])
        
        runs_br = np.vstack(runs_br).T

        # Main results
        main_res = np.load(os.path.join(paths[1], file_))[:,0]

        var, se, bias, cf = do_statistics(runs_br, runs_num, main_res)
  
        se_res.append(se)
        var_res.append(var)    
        bias_res.append(bias)
        cf_res.append(cf)

    se_res = np.vstack(se_res).T
    var_res = np.vstack(var_res).T
    bias_res = np.vstack(bias_res).T
    cf_res = np.vstack(cf_res).T    

    # Save the results
    config.logger.info('Save Jackknife test results.')

    np.save(os.path.join(path_JN_main, 'variance.npy'), var_res)
    np.save(os.path.join(path_JN_main, 'standard_error.npy'), se_res)
    np.save(os.path.join(path_JN_main, 'bias.npy'), bias_res)
    np.save(os.path.join(path_JN_main, 'confidence.npy'), cf_res)

    config.logger.info('End of Jackknife test')

    return

def Bootstrap_stats(st, stations, paths, repeats, perce):
   """
   Perform the Bootstrap test to determine the 95% C.I. for the Bright Spots
   Also caclulates, Standard Deviation and Standard error

   Arguments:
    ------
    st: Obspy Object 
        Obspy Stream Object
    stations: dict
        Dictionary with stations info
    paths: list
        List with various paths
    repeats: int
        Number of resampling repeats
    perce:
        Percentage of stations to resample from each sector
    
    Returns:
    --------
    - 
   """ 

   config.logger.info('--------------')
   config.logger.info('Bootstrap Test')
   config.logger.info('--------------')

   createDir(paths[0])
   # detailed solutions
   path_BOOT_de = createDir(os.path.join(paths[0], 'Detailed_Solutions'))
   # main solution
   path_BOOT_main = createDir(os.path.join(paths[0], 'Main_Solution'))

   #Get the stations from each sector
   sc = backprojection.sectorsCheck(stations, [config.org.latitude, config.org.longitude],\
                                    config.cfg['Backprojection']['Sectors'])

   #Check that the given percentage is satisfied from the number of stations 
   #in each sector
   keys = sc[0].keys()

   accepted_sectors = [] 
   for key in keys:
       if len(sc[0][key])*perce/100>=1:
           accepted_sectors.append(key)
   if len(accepted_sectors)==0:
       config.logger.info('No sector that satisfies the percentage input. Raise the percetange parameter from the configuration file')
       return
 
   #Get randomly from each sector the percentage of stations
   stations_ = []
   for i in accepted_sectors:
       number_of_station_in_sec = len(sc[0][i])
       number_of_get_stations = int(np.floor(number_of_station_in_sec * perce/100))
       stations_.append(random.sample(sc[0][i], number_of_get_stations))
   #Flatten the list
   stations_ = [x for l in stations_ for x in l]

   if len(stations_)<=3:
       config.logger.info('Small number of stations. Raise the percetange parameter from the configuration file')
       return

   ####
   # Here we will keep the station that don't participate in the 
   # bootstrap game
   stClear = Stream()
   stationClear = {}

   for tr in st:
       if (tr.stats.network+'.'+tr.stats.station) not in stations_:
           stClear.append(tr.copy())

   keys = stations.keys()
   for key in keys:
       if stations[key][0]+'.'+stations[key][1] not in stations_:
           stationClear[key] = stations[key]
   #############################################################
   #############################################################

   #Now lets start the resampling circles
   for i in range(repeats):

       # Stream and dictionary to backproject
       stBoot = Stream()
       stationsBoot = {}

       #Add the stable stations
       for tr in stClear:
           stBoot.append(tr)
       stationsBoot.update(stationClear)

       #Choose a random number of stations to get
       num_sta = random.randrange(3, len(stations_)+1)
       #resample the stations
       sta_resa = resample(stations_, replace=False, n_samples=num_sta)
       #Append these stations also to stream and dictionary
       for sta in sta_resa:
           stBoot.append(st.select(network =sta.split('.')[0],\
                         station=sta.split('.')[1])[0].copy())
           stationsBoot[sta.split('.')[1]] = stations[sta.split('.')[1]]       
       Boot_path_save = os.path.join(path_BOOT_de, str(i))
       backprojection.backprojection(stBoot, stationsBoot, Boot_path_save)

   # Calculate the statistics

   # Get the results of the first run
   runs = sorted([int(d) for d in os.listdir(path_BOOT_de)])

   npy_files = [file.split('/')[-1] for file in glob.glob(os.path.join(path_BOOT_de,str(runs[0]), "out*.npy")) if "out_Max.npy" not in file]

   # Sorted files
   npy_files = sorted(npy_files, key=extract_number_from_path)

   # How many runs ?
   runs_num = len(os.listdir(path_BOOT_de))

   # Keep to lists
   se_res = []
   var_res = []
   bias_res = []
   cf_res = []

   # Loop for the timestep files:
   for npy_i in range(len(npy_files)):
       file_ = npy_files[npy_i]

       # In this array keep the results
       runs_br = []

       # Loop runs
       for d in runs:
           # Load files
           new_path = os.path.join(path_BOOT_de, str(d), file_)

           res_load = np.load(new_path)
           runs_br.append(res_load[:,0])

       runs_br = np.vstack(runs_br).T

       # Main results
       main_res = np.load(os.path.join(paths[1], file_))[:,0]

       var, se, bias, cf = do_statistics(runs_br, runs_num, main_res)

       se_res.append(se)
       var_res.append(var)
       bias_res.append(bias)
       cf_res.append(cf)

   # Save the results
   config.logger.info('Save Bootstrap test results.')

   se_res = np.vstack(se_res).T
   var_res = np.vstack(var_res).T
   bias_res = np.vstack(bias_res).T
   cf_res = np.vstack(cf_res).T

   np.save(os.path.join(path_BOOT_main, 'variance.npy'), var_res)
   np.save(os.path.join(path_BOOT_main, 'standard_error.npy'), se_res)
   np.save(os.path.join(path_BOOT_main, 'bias.npy'), bias_res)
   np.save(os.path.join(path_BOOT_main, 'confidence.npy'), cf_res)

   config.logger.info('End of Bootstrap test')

   return
