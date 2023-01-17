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

import os, copy, scipy, random
from scipy.special import erfinv
import numpy as np
from sklearn.utils import resample

from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import kilometer2degrees
from obspy.core.stream import Stream

import matplotlib.pyplot as plt
# local
from SSA2py.core import backprojection, config
from SSA2py.core.basic_f.other import createDir, write_txt_Tests, delete_npy
from SSA2py.core.basic_f.epicDistance import _1Ddist

def bootstrap_basic_statistics(boot, confidence_level=0.95):
    """
    Calculate basic statistic parameters using the results from Bootstrap
    of Jackknife test.

    Arguments:
    ----------
    boot: array-like [grid points number, Bootstrap Efforts]
        Results from Boot or Jack

    Returns:
    --------
    SE: array-like
        Standard Error
    C.I.: array-like
        Confidence Interval
    STD: array-like
        Standard Deviation
    """


    #Calculations
    STD = np.nanstd(boot, axis=1)
    SE = STD / np.sqrt(np.count_nonzero(~np.isnan(boot), axis=1))
    
    z_score = np.sqrt(2.0)*erfinv(confidence_level)
    CI = z_score*np.array((SE))

    return  (STD, SE, CI)


def Jackknife_stats(st, stations, Jack_Paths, confidence_level=0.95):
    """
    Perform the Jackknife test to determine the 95% C.I. for the Bright Spots
    Also caclulates, Standard Deviation and Standard error


    Arguments:
    ------
    st: Obspy Object 
        Obspy Stream Object
    stations: dict
        Dictionary with stations info
    Jack_Paths: list
        List with various paths
    confidence_level: float
        Confidence level for the confidence interval of the Jackknife estimate
        Real value between (0,1)
    
    Returns:
    --------
    - 

    """
   
    config.logger.info('--------------')
    config.logger.info('Jackknife Test')
    config.logger.info('--------------')

    #create dirs
    path_JN=createDir(Jack_Paths[0])

    # detailed solutions
    path_JN_de = createDir(os.path.join(Jack_Paths[0], 'Detailed_Solutions'))

    # main solution
    path_JN_main = createDir(os.path.join(Jack_Paths[0], 'Main_Solution'))

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

        if config.cfg['Tests']['Delete'] is True:
            delete_npy(Jack_path_save, ['out_Max.npy'])                 

    # Calculate the 4-D statistics
    SSA = np.load(os.path.join(Jack_Paths[1], 'out_Max.npy')) # SSA MAX results

    time_steps = SSA[:,-1]
   
    #Pre-allocate array for the maximum brightness
    max_array = np.zeros((len(time_steps), 14), dtype = np.float64)

    if config.cfg['Tests']['Delete'] is False:
        for s in range(len(time_steps)): #Loop throught time steps, calculate C.I. for BR
            # read the SSA results for this step
            SSA_s = np.load(os.path.join(Jack_Paths[1], 'out_'+ str(float(time_steps[s])) +'.npy'))

            jack_stat_BR = np.zeros((SSA_s[:,0].shape[0], st.count())) 

            # now the jackknife results
            for i_ in range(st.count()):
                jack_s = np.load(os.path.join(path_JN_de, str(i_),\
                                'out_'+ str(float(time_steps[s])) +'.npy'))
                jack_stat_BR[:,i_] = jack_s[:,0]
            # calculate the statistics
            JackBR = bootstrap_basic_statistics(jack_stat_BR, confidence_level=0.95)
            # merge all these arrays and save it  
            all_ = np.column_stack((JackBR[2].T, JackBR[0].T, JackBR[1].T,\
                                    SSA_s[:,1], SSA_s[:,2], SSA_s[:,3])) 
            np.save(os.path.join(path_JN_main, 'out_'+ str(float(time_steps[s])) +'.npy'), all_)

    #####################################################
    #####################################################
    # do this for the max BR x and y 
    jack_stat_XY = np.zeros((len(time_steps), st.count()))
    jack_stat_Z = np.zeros((len(time_steps), st.count()))
    jack_stat_BR = np.zeros((len(time_steps), st.count()))

    for i_ in range(st.count()):
        jack_max = np.load(os.path.join(path_JN_de, str(i_),\
                           'out_Max.npy'), allow_pickle=True)
        # calculate distance with Main SSA result
        for d in range(len(jack_max[:,0])):
            jack_stat_XY[d, i_] = kilometer2degrees(gps2dist_azimuth(jack_max[d,2],\
                                  jack_max[d,1], SSA[d,2], SSA[d,1])[0]/1000)        
        jack_stat_BR[:,i_] = jack_max[:,0]
        jack_stat_Z[:,i_] = jack_max[:,3]

    # jackknife results
    JackXY = bootstrap_basic_statistics(jack_stat_XY, confidence_level=0.95)

    #print(np.max(JackXY[2][0]))
    JackZ = bootstrap_basic_statistics(jack_stat_Z, confidence_level=0.95)
    JackBR = bootstrap_basic_statistics(jack_stat_BR, confidence_level=0.95)

    # keep maximum array
    max_array[:,0] = SSA[:,0]
    max_array[:,1] = JackBR[2] #CI
    max_array[:,2] = JackBR[0] #STD
    max_array[:,3] = JackBR[1] #SE

    max_array[:,4] = SSA[:,1]
    max_array[:,5] = SSA[:,2]
    max_array[:,6] = JackXY[2] #CI
    max_array[:,7] = JackXY[0] #STD
    max_array[:,8] = JackXY[1] #SE

    max_array[:,9] = SSA[:,3]
    max_array[:,10] = JackZ[2] #CI
    max_array[:,11] = JackZ[0] #STD
    max_array[:,12] = JackZ[1] #SE

    max_array[:,13] = SSA[:,-1]

    np.save(os.path.join(path_JN_main, 'out_Max.npy'), max_array)
    config.logger.info('End of Jackknife test')

    write_txt_Tests(os.path.join(path_JN_main, 'out_Max.npy'),\
                    os.path.join(path_JN_main, 'out_Max_file.txt'))     

    return

def Bootstrap_stats(st, stations, Boot_Paths, repeats, perce, confidence_level=0.95):
   """
   Perform the Bootstrap test to determine the 95% C.I. for the Bright Spots
   Also caclulates, Standard Deviation and Standard error

   Arguments:
    ------
    st: Obspy Object 
        Obspy Stream Object
    stations: dict
        Dictionary with stations info
    Jack_Paths: list
        List with various paths
    repeats: int
        Number of resampling repeats
    perce:
        Percentage of stations to resample from each sector
    confidence_level: float
        Confidence level for the confidence interval of the Jackknife estimate
        Real value between (0,1)
    
    Returns:
    --------
    - 
   """ 

   config.logger.info('###################################')
   config.logger.info('Bootstrap Test')
   config.logger.info('###################################')

   path_BOOT=createDir(Boot_Paths[0])
   # detailed solutions
   path_BOOT_de = createDir(os.path.join(Boot_Paths[0], 'Detailed_Solutions'))
   # main solution
   path_BOOT_main = createDir(os.path.join(Boot_Paths[0], 'Main_Solution'))

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

       if config.cfg['Tests']['Delete'] is True:
            delete_npy(Boot_path_save, ['out_Max.npy'])

   #Calculate the 4-D statistics
   SSA = np.load(os.path.join(Boot_Paths[1], 'out_Max.npy')) # SSA MAX results

   time_steps = SSA[:,-1]

   #Pre-allocate array for the maximum brightness
   max_array = np.zeros((len(time_steps), 14), dtype = np.float64)

   if config.cfg['Tests']['Delete'] is False:
        for s in range(len(time_steps)): #Loop throught time steps, calculate C.I. for BR
             # read the SSA results for this step
             SSA_s = np.load(os.path.join(Boot_Paths[1], 'out_'+ str(float(time_steps[s])) +'.npy'))

             boot_stat_BR = np.zeros((SSA_s[:,0].shape[0], repeats))

             # now the boot results
             for i_ in range(repeats):
                  boot_s = np.load(os.path.join(path_BOOT_de, str(i_),\
                                  'out_'+ str(float(time_steps[s])) +'.npy'))
                  boot_stat_BR[:,i_] = boot_s[:,0]
             # calculate the statistics
             BootBR = bootstrap_basic_statistics(boot_stat_BR, confidence_level=0.95)
             # merge all these arrays and save it  
             all_ = np.column_stack((BootBR[2].T, BootBR[0].T, BootBR[1].T,\
                                     SSA_s[:,1], SSA_s[:,2], SSA_s[:,3]))
             np.save(os.path.join(path_BOOT_main, 'out_'+ str(float(time_steps[s])) +'.npy'), all_)
   #####################################################
   #####################################################
   # do this for the max BR x and y 
   boot_stat_XY = np.zeros((len(time_steps), repeats))
   boot_stat_Z = np.zeros((len(time_steps), repeats))
   boot_stat_BR = np.zeros((len(time_steps), repeats))

   for i_ in range(repeats):
       boot_max = np.load(os.path.join(path_BOOT_de, str(i_),\
                          'out_Max.npy'), allow_pickle=True)
       # calculate distance with Main SSA result
       for d in range(len(boot_max[:,0])):
           boot_stat_XY[d, i_] = kilometer2degrees(gps2dist_azimuth(boot_max[d,2],\
                                 boot_max[d,1], SSA[d,2], SSA[d,1])[0]/1000)
       boot_stat_BR[:,i_] = boot_max[:,0]
       boot_stat_Z[:,i_] = boot_max[:,3]

   # jackknife results
   BootXY = bootstrap_basic_statistics(boot_stat_XY, confidence_level=0.95)

   BootZ = bootstrap_basic_statistics(boot_stat_Z, confidence_level=0.95)
   BootBR = bootstrap_basic_statistics(boot_stat_BR, confidence_level=0.95)

   # keep maximum array
   max_array[:,0] = SSA[:,0]
   max_array[:,1] = BootBR[2] #CI
   max_array[:,2] = BootBR[0] #STD
   max_array[:,3] = BootBR[1] #SE

   max_array[:,4] = SSA[:,1]
   max_array[:,5] = SSA[:,2]
   max_array[:,6] = BootXY[2] #CI
   max_array[:,7] = BootXY[0] #STD
   max_array[:,8] = BootXY[1] #SE

   max_array[:,9] = SSA[:,3]
   max_array[:,10] = BootZ[2] #CI
   max_array[:,11] = BootZ[0] #STD
   max_array[:,12] = BootZ[1] #SE

   max_array[:,13] = SSA[:,-1]

   np.save(os.path.join(path_BOOT_main, 'out_Max.npy'), max_array)

   write_txt_Tests(os.path.join(path_BOOT_main, 'out_Max.npy'),\
                    os.path.join(path_BOOT_main, 'out_Max_file.txt'))

   config.logger.info('End of Bootstrap test') 
   return



