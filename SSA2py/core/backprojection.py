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

import os, yaml, multiprocessing, math, re, subprocess, operator
from time import perf_counter
import numpy as np

# Obspy Imports
################

from obspy.geodetics.base import gps2dist_azimuth
from obspy.core.utcdatetime import UTCDateTime

# Other Imports
###############


import numba
from numba import njit, cuda, float32, float64

numba.config.THREADING_LAYER = 'default'

# Local Imports
################

from SSA2py.core import config
from SSA2py.core.grid import grid_box, grid_box3D, grid_points
from SSA2py.core.basic_f.get_tt import get_tt1d, get_tt3d
from SSA2py.core.basic_f.other import createDir, write_txt_maxBright
from SSA2py.core.basic_f.epicDistance import _1Ddist
from SSA2py.core.weights import dist_weights, read_external_weights, azimuth_weights
from SSA2py.core.tt_corrections import exec_corr

# End of Imports
################


def backprojection(stream, stations, savedir):
    """
    Main Backprojection Function

    Arguments:
    ------
    stream: Obpsy Object
        Stream Obspy Object with traces
    stations: dict
        Dictionary with stations info
    savedir: string
        Directory to save results

    Returns:
    -------
    -

    """

    # sort the dictionary based on the station name
    stations_ = sorted(stations)
    stations = {key:stations[key] for key in stations_}

    # sort stream based on station name
    stream.sort(keys=['station'])

    # check the sectors criteria. Can we backproject?
    secCheck = sectorsCheck(stations, [config.org.latitude, config.org.longitude],\
                           config.cfg['Backprojection']['Sectors'])

    config.logger.info('-------Sectors-------')
    config.logger.info(secCheck[0])    

    if secCheck[1]==False:
        config.logger.warning('The sectors criteria are not satisfied. End execution.')
        return False   
    else:
        config.logger.info('The sectors criteria are satisfied.')

    config.logger.info('------------------------')
    config.logger.info(config.job) 
    config.logger.info('------------------------')
    config.logger.info('Waveform Type: ' + config.cfg['Streams']['Type'])
    config.logger.info('Filter: ' + str(config.fi))


    # Create the time array for BP
  
    ScanningTime = config.scanningRules[0]
  
    TimeShift = config.cfg['Backprojection']['Settings']['TimeShift']

    # Create the time array

    time = np.array([round(ScanningTime[0] + i * TimeShift, 2) for i in range(int((ScanningTime[1] - ScanningTime[0]) / TimeShift) + 1)])

    config.logger.info('Scanning Time = [' + str(ScanningTime[0]) + ',' + str(ScanningTime[1]) + '] and TimeShift = ' + str(TimeShift))

    # Main Backprojection
    # re-arrange the grid (List with x, y, depth)

    # Directory to save
    path_br=createDir(savedir)

    # Calculate the SSA Window Weights

    side_1 = int(config.cfg['Backprojection']['Settings']['MovingWindow'][0] * stream[0].stats.sampling_rate)
    side_2 = int(config.cfg['Backprojection']['Settings']['MovingWindow'][1] * stream[0].stats.sampling_rate)   

    w = weight(config.cfg['Backprojection']['Settings']['Weight'], side_1 + side_2)

    # Amplitude domination?
    stream, stations = ampliDom(stream, stations)

    config.logger.info('Get traveltimes between points of the grid - stations')

    # Convert to array for numba
    Scoords = np.array([[tr.stats.stalon, tr.stats.stalat, tr.stats.staelev] for tr in stream], dtype='float32')

    # Get possible tt corrections
    if config.cfg['Streams']['Corrections'][0]==True and config.job != 'Array Response Function':
        cor = exec_corr(stations, config.cfg['Streams']['Corrections'][2]) 

        # Corrections from seconds to samples
        cor = np.round(np.array(cor, dtype='float32'),2)
    else:
        cor = np.zeros((len(stations),), dtype='float32')

    # 1D or 3D model?

    if config.model['type']=='1D':
        
        # Elevation

        elevationArray = np.arange(0, float(config.model['Elevation']) +\
                                   config.model['Granularity'],\
                                   config.model['Granularity'])

        # Distance (round to keep 3 decimals ~ meters)

        dist = np.round(np.array(_1Ddist(config.grid, Scoords[:,0], Scoords[:,1]), dtype='float32'), decimals=3)


        # For each station identify the 1D table position
        
        tablpos = np.array([np.where(abs(elevationArray-(sel/1000))==np.min(abs(elevationArray-(sel/1000))))[0][0]\
                            for sel in Scoords[:,2]]) 

        # In traveltimes keep also three decimals
 
        tt = np.round(np.array(get_tt1d(config.grid, config.tables, tablpos,\
                      dist, np.array(config.model['Granularity']), Scoords[:,2]), dtype='float32'), decimals=3)

    if config.model['type']=='3D':
         #Remove paths missing from config.stations
         tables = [i for i in config.tables if i.split('/')[-1].split('_')[0] in [*stations]]
         tt = np.array(get_tt3d(config.grid, stations, tables,\
                       config.model['Granularity'],\
                       config.posx, config.posy,\
                       config.gridRules[0][3],\
                       config.gridRules[0][4]), dtype='float32')

    # convert st trace data to list of numpy arrays

    # Keep the traces data !

    data = np.round(np.array([tr.data for tr in stream], dtype='float32'),4)

    # Get rates !
    
    rates = np.array([tr.stats.sampling_rate for tr in stream], dtype='float32')

    # Fix look times in traces

    time_ = np.array((UTCDateTime(config.org.time) - stream[0].stats.starttime) + time, dtype='float32')
   
    # Basic parameters

    # Window

    win = np.array([config.cfg['Backprojection']['Settings']['MovingWindow'][0],\
                    config.cfg['Backprojection']['Settings']['MovingWindow'][1]], dtype='float32')

    # Max traveltime

    maxTT = np.float32(config.cfg['Backprojection']['Settings']['TTmax'])

    # Stations threshold

    StaThre = np.float32(config.cfg['Backprojection']['Settings']['StaThre'])
    
    # N-Root Power

    power = np.float32(config.cfg['Backprojection']['Settings']['Npower'])
    
    # Brightness Type

    BrType = np.float32(config.cfg['Backprojection']['Settings']['BrType'])

    # Start SSA execution 

    config.logger.info('SSA execution...')

    t1 = perf_counter()

    # GPU execution?
    if config.cfg['Backprojection']['GPU']==True:
        config.cfg['Backprojection']['GPU'] = checkGPU()
        # Still GPU True?
        if config.cfg['Backprojection']['GPU']==True:
            #decide threads and blocks
            threadsperblock = (16, 16) #threads number

            #Get grid size and blocks
            blockspergrid_x = math.ceil(tt.shape[0]/threadsperblock[0])
            blockspergrid_y = math.ceil(time_.shape[0]/threadsperblock[1])
            blockspergrid = (blockspergrid_x, blockspergrid_y)

            # Move data to GPU
            data_gpu = cuda.to_device(data)
            rates_gpu = cuda.to_device(rates)

            bright_gpu = cuda.device_array_like(np.empty((tt.shape[0], time_.shape[0])))

            # For the phase
            tt_gpu = cuda.to_device(tt)
            cor_gpu = cuda.to_device(cor)
            time_gpu = cuda.to_device(time_)

            t1 = perf_counter()

            #window?
            if win[0]==0 and win[1]==0:
                SSA_cuda[blockspergrid, threadsperblock](tt_gpu, time_gpu, data_gpu,\
                         rates_gpu, bright_gpu, maxTT, StaThre, power, BrType, cor_gpu)
            else:
                w_gpu = cuda.to_device(w)
                win_gpu = cuda.to_device(np.array([win[0],win[1]]))
                SSA_cudaW[blockspergrid, threadsperblock](tt_gpu, time_gpu, data_gpu,\
                          rates_gpu, bright_gpu, maxTT, StaThre, power, BrType, w_gpu, np.sum(w), win_gpu, cor_gpu)

            cuda.synchronize() 
            br = bright_gpu.copy_to_host()  
       
    # CPU execution
    if config.cfg['Backprojection']['GPU']==False:
        if win[0]==0 and win[1]==0:
            with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as pool:
                br = pool.starmap(SSA, zip(tt, len(tt) * [time_], len(tt) * [data], len(tt) * [rates],\
                                  len(tt)*[maxTT], len(tt)*[StaThre], len(tt)*[power], len(tt) * [BrType], len(tt) * [cor]))
        else:
            with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as pool:
                br = pool.starmap(SSAWindow, zip(tt, len(tt) * [time_], len(tt) * [data], len(tt) * [rates],\
                                  len(tt)*[w], len(tt)*[np.sum(w)], len(tt)*[win],len(tt)*[maxTT], len(tt)*[StaThre], len(tt)*[power], len(tt) * [BrType], len(tt) * [cor]))

    # Fix the output to a proper array
    br = np.array(br, dtype='float32')

    config.logger.info('Save results...')   

    # Pre-allocate array for the maximum brightness
    max_array = np.zeros((len(time), 5))

    for step in range(time.shape[0]):
        all_ = np.column_stack((br[:,step], config.grid, np.full((br[:,step].shape), time[step])))

        # Add to the max array
        max_index = np.argmax(br[:,step], axis=0)
        max_array[step,0] = br[max_index, step] #Brightness
        max_array[step,1] = config.grid[max_index,0] #Lon.
        max_array[step,2] = config.grid[max_index,1] #Lat.
        max_array[step,3] = config.grid[max_index,2] #Depth
        max_array[step,4] = time[step] #Time
        
        # Brightness <= bthre*bmax as NaN
        all_[all_[:,0]<max_array[step,0]*config.cfg['Backprojection']['Settings']['bthre']] = 0

        # Save the time step array
        np.save(os.path.join(path_br, 'out_{}.npy'.format(float(time[step]))), all_)
        
    # Save the Max array
    np.save(os.path.join(path_br, 'out_Max.npy'), max_array)

    # save the grid
    np.save(os.path.join(path_br, 'grid.npy'), config.grid)

    t2 = perf_counter()

    # normalize if product
    if int(config.cfg['Backprojection']['Settings']['BrType'])==1 or config.cfg['Backprojection']['Settings']['Normalize Results']==True:
        normalizeProd(path_br)

    config.logger.info('------------------------')
    config.logger.info('Output Saved --> ' + os.path.join(path_br))
    config.logger.info('SSA finished...')
    config.logger.info('Time: ' + str(np.around(t2-t1,1)) + ' sec') 
    config.logger.info('------------------------')
    
    # save the travel times and grid array
    np.save(os.path.join(path_br, 'tt.npy'), tt)
    
    # save txt file
    write_txt_maxBright(os.path.join(path_br, 'out_Max.npy'),\
                        os.path.join(path_br, 'out_Max_file.txt'))          

    return True     

@njit(float32[:](float32[:], float32[:], float32[:,:], float32[:], float32, float32, float32, float32, float32[:]))
def SSA(tt, time, data, rate, maxTT, StaThre, power, type_, cor):
    """
    SSA for CPUS

    Arguments:
    ----------
    tt: array-like
        Traveltimes per grid-station.
    time: array-like
        Timesteps.
    data: array-like
        Waveforms.
    rate: array-like
        Sampling rate of waveforms.
    maxTT: float
        Maximum traveltime to take into account.
    StaThre: float
        Percentage of station to get grid point.
    power: float
        Power for the brightness
    type_: int
        Type of SSA (0 sum, 1 multi)
    cor: array-like
        Time Corrections.

    Returns:
    --------
    br: array-like
        Brightness values (gridpoints x timesteps)

    """

    MAX_INT = 2**31 - 1 # Mind possible overflow

    # Pre-allocate arrays

    # Brightness

    br = np.zeros((time.shape[0]), dtype='float32')

    # Theoretical Arrivals

    _tt_ = np.zeros((time.shape[0]), dtype='float32')

    for t in range(time.shape[0]):

        _tt_ = time[t] + tt + cor # theoretical arrival time

        count_stations = 0

        for d in range(data.shape[0]): 
            if tt[d]<maxTT:
                count_stations += 1

                # Do corrections
                
                if type_ == 0 or (type_ == 1 and d==0):
                    br[t] += abs(data[d][int(np.around(_tt_[d]*rate[d]))])

                if type_ == 1 and d>0:
                    br[t] = br[t] * abs(data[d][int(np.around(_tt_[d]*rate[d]))])

        if count_stations>StaThre*data.shape[0]:
            br[t] = ((1/count_stations) * br[t]) ** power
            if br[t]>=MAX_INT:
                br[t] = 0
        else:
            br[t] = 0

    return br

@njit(float32[:](float32[:], float32[:], float32[:,:], float32[:], float32[:], float32, float32[:], float32, float32, float32, float32, float32[:]))
def SSAWindow(tt, time, data, rate, w, sumW, win, maxTT, StaThre, power, type_, cor):
    """
    SSA with window

    w: array-like
        Weights array for the window.
    sumW: float
        Sum weight.
    win: array-like
        Window in seconds for each side.

    """

    MAX_INT = 2**31 - 1 # Mind possible overflow

    # Pre-allocate the amplitudes array

    # Brightness

    br = np.zeros((time.shape[0]), dtype='float32')

    # Theoretical Arrivals

    _tt_ = np.zeros((time.shape[0]), dtype='float32')

    for t in range(time.shape[0]):

        _tt_ = tt + time[t] + cor # theoretical arrival

        count_stations = 0
        for d in range(data.shape[0]):

            if tt[d]<maxTT:
                count_stations += 1

                dt = np.abs(data[d][int(np.around(((_tt_[d])*rate[d]) - (win[0]*rate[d]))):\
                            int(np.around(((_tt_[d])*rate[d]) + (win[1]*rate[d])))])
                the_sum = 0

                for s in range(len(dt)):
                    the_sum = the_sum + (np.float32(dt[s]) * w[s])

                if type_ == 0 or (type_ == 1 and d==0):
                    br[t] += (the_sum/sumW)
                if type_ == 1 and d>0:
                    br[t] = br[t] * (the_sum/sumW)
        if count_stations>StaThre*data.shape[0]:
            br[t] = ((1/count_stations) * br[t]) ** power
            if br[t]>=MAX_INT:
                br[t] = 0
        else:
            br[t] = 0
    return br


@cuda.jit
def SSA_cuda(tt, time, data, rate, bright, maxTT, StaThre, power, type_, cor):
    """
    SSA execution for one phase (P or S) for GPU.
    NO Window case.

    Arguments:
    ----------
    tt: array-like
        Traveltimes per grid-station.
    time: array-like
        Timesteps.
    data: array-like
        Waveforms.
    rate: array-like
        Sampling rate of waveforms.
    bright: array-like
        Pre-allocated brightness array.
    maxTT: float
        Maximum traveltime to take into account.
    StaThre: float
        Percentage of station to get grid point.
    power: float
        Power for the brightness
    type_: int
        Type of SSA (0 sum, 1 multi)
    cor: array-like
        Time Corrections.

    Returns:
    --------
    br: array-like
        Brightness values (gridpoints x timesteps)

    """
    i, j = cuda.grid(2) #grid point, time step
    di, dj = cuda.gridsize(2)

    # Grid point
    for i1 in range(i, tt.shape[0], di):
        # Time        
        for j1 in range(j, time.shape[0], dj):
            count_stations = 0
            for d_ in range(data.shape[0]):
                if maxTT>=tt[i1,d_]:
                    count_stations += 1
                    if type_ == 0 or (type_ == 1 and d_ == 0):
                        bright[i1,j1] += abs(data[d_][int(round((tt[i1,d_]+time[j1]+cor[d_])*rate[d_]))])
                    if type_ == 1 and d_>0:    
                        bright[i1,j1] = bright[i1,j1] * abs(data[d_][int(round((tt[i1,d_]+time[j1]+cor[d_])*rate[d_]))])
            if count_stations>data.shape[0]*StaThre:
                bright[i1,j1] = (1/count_stations) * bright[i1,j1] ** power
            else:
                bright[i1,j1] = 0



@cuda.jit
def SSA_cudaW(tt, time, data, rate, bright, maxTT, StaThre, power, type_, w, sumW, win, cor):
    """
    SSA with window

    w: array-like
        Weights array for the window.
    sumW: float
        Sum weight.
    win: array-like
        Window in seconds for each side.

    """
    i, j = cuda.grid(2) #grid point, time step
    di, dj = cuda.gridsize(2)

    for i1 in range(i, tt.shape[0], di):
        for j1 in range(j, time.shape[0], dj):
            count_stations = 0
            for d_ in range(data.shape[0]):
                if maxTT>=tt[i1,d_]:
                    count_stations += 1
                    dt_ = data[d_][int(round(((tt[i1,d_]+time[j1]+cor[d_])*rate[d_]) - (win[0]*rate[d_]))):int(round(((tt[i1,d_]+time[j1]+cor[d_])*rate[d_]) + (win[1]*rate[d_])))]
                    the_sum = 0
                    for k in range(dt_.shape[0]):
                        the_sum += (abs(dt_[k])*w[k])
                    if type_ == 0 or (type_ == 1 and d_ == 0):
                        bright[i1,j1]+= operator.truediv(the_sum, sumW)
                    if type_ == 1 and d_>0:    
                        bright[i1,j1] = bright[i1,j1] * operator.truediv(the_sum, sumW)
            if count_stations>data.shape[0]*StaThre:
                bright[i1,j1] = (1/count_stations) * bright[i1,j1] ** power
            else:
                bright[i1,j1] = 0

def checkGPU():
    """
    Check if we have NVIDIA GPU.

    """

    try:
        subprocess.check_output('nvidia-smi')
        config.logger.info('Nvidia GPU detected!')
        return True
    except:
        config.logger.warning('No Nvidia GPU detected!')
        config.logger.info('Continue with CPU execution')
        return False



def getGrid():
    """
    Get grid based on the config input

    Returns:
    --------
        G: array-like
            3D Grid

    """

    # Check if the input model is 1D or 3D

    # What we have in the tables directory?

    with open(os.path.join(config.cfg['Traveltimes']['Save'], 'model.yml')) as file:
             char = yaml.load(file, Loader=yaml.FullLoader)
    config.model = char

    config.logger.info('Creating the Grid')

    grid_type = config.gridRules[0]


    if char['type']=='1D':

        #Check if the input is list 
        if isinstance(grid_type, list):

            ############
            # Box Grid #
            ############

            if len(grid_type)==6 and grid_type[0]=='box':   

                #Create the grid the coordinates (x, y, z)
                grid_ = grid_box(grid_type[1], grid_type[2], grid_type[3], grid_type[4], grid_type[5], [config.org.latitude, config.org.longitude])

                gx = grid_[0]; gy = grid_[1];

                horzG = np.repeat(np.column_stack((gx, gy)), len(grid_[2]), axis=0)
                vertG = np.ravel(np.tile(grid_[2], (len(gx), 1)))

                G = np.column_stack((horzG, vertG))
             
                config.logger.info('Box Grid!')

                return np.array(G, dtype='float32')

            if len(grid_type)==2 and grid_type[0]=='points': #Points

                grid_ = grid_points(grid_type[1])

                G = np.column_stack((grid_[0], grid_[1], grid_[2]))
  
                config.logger.info('Input Points Grid!')

                return np.array(G, dtype='float32')
                
                
 
    if char['type']=='3D':       

        #Respect the grid rules from the tt calculation (only Box Grid!)

        grid_ = grid_box3D(grid_type[1], grid_type[2], grid_type[3], grid_type[4], grid_type[5],\
                          [config.org.latitude, config.org.longitude], char['Distance'])          

        gx = grid_[0]; gy = grid_[1]; config.posx = grid_[3]; config.posy = grid_[4];

        horzG = np.repeat(np.column_stack((gx, gy)), len(grid_[2]), axis=0)
        vertG = np.ravel(np.tile(grid_[2], (len(gx), 1)))

        G = np.column_stack((horzG, vertG))

        config.logger.info('Box Grid!')

        return np.array(G, dtype='float32')  

def ampliDom(stream, stations):
    """
    Check to find if the amplitudes will be dominated (>50%) 
    by only one station.

    Arguments:
    ------
    stream: Obspy Object
         Obspy Stream Object
    stations: dict
         Stations Dictionary

    Returns:
    -------
    stream: Obspy Object
         Obspy Stream Object
    stations: dict
         Stations Dictionary

    in tuple format
 
    """

    Ampls = [max(abs(stream[i].data)) for i in range(len(stream))]
    Ampl_dom = np.where(Ampls>(sum(Ampls)/2))[0]

    #Delete the station from the dictionary and the traces 
    if Ampl_dom.size != 0: 
        try:
            stream.remove(stream[Ampl_dom[0]])
            del stations[stream[Ampl_dom[0]].stats.station]
        except:
            pass
    return (stream, stations)

def weight(typeN, nmPoints):
    """
    Weight the piece of window

    # Translated and adapted from original code by Honn Kao.
    # Personal communication.

    Arguments:
    ---------
    typeN: int
        Type of traveltime weightning 
    nmPoints: int
        Number of Weightning Points

    Returns:
    --------
        w: array-like (float32)
            Weightning array
    """

    if typeN==1: #Equal

        w = np.ones((nmPoints,), dtype=int)

    elif typeN==0: #Gaussian

        nwin0 = int(float(nmPoints) / 3.0 + 0.5)
        ngaussian = nmPoints - nwin0

        w = np.zeros((nmPoints,))

        for ix in range(1,int(nmPoints) + 1):

            w1 = -1.0 * float((ix - nwin0) * (ix - nwin0))

            w2 = 2.0 * float(ngaussian * ngaussian)

            w[ix - 1] = np.exp(w1 / w2)
    else:
        pass
    return w.astype(np.float32) #(Return float32 array)

def sectorsCheck(stations, epicenter, sectors):
    """
    Check if the stations distributions satisfy the sectors
    criteria.

    Arguments:
    ----------
    stations: dict
        Dictionary with stations info
    epicenter: list
        List with epicentral info
    sectors: list
        Dictionary with sectors
 
    Returns:
    --------
    sc: dict
        Dictionary with statiosn in sectors
    boolean

    in tuple format

    """


    #Count the stations in each sector
    sc = {0:list([]), 1:list([]), 2:list([]), 3:list([])}
    for sta in stations:
         #Azimuth
         azm = gps2dist_azimuth(epicenter[0], epicenter[1], stations[sta][3], stations[sta][2])[1]
         pos = math.floor(azm/90)

         sc[pos]+=[stations[sta][0]+'.'+stations[sta][1]]

    count = sum(map(lambda x : x >= sectors[1], [len(value) for key, value in sc.items()])) 
    return (sc, bool(count>=sectors[0]))

def normalizeProd(path):
    """
    In cases where the SSA comes from using product normalize the values

    Arguments:
    ----------
    path: str
        Path with results

    Returns:
    --------

    """

    # load maximum brightness
    maxBr = np.load(os.path.join(path,'out_Max.npy'))

    # get meximum value
    maxIndex = np.argmax(maxBr[:,0])
    maxValue = maxBr[maxIndex,0]

    # normalize the maximum
    maxBr[:,0] = maxBr[:,0]/maxValue

    # write
    np.save(os.path.join(path,'out_Max.npy'), maxBr) 

    # normalize the rest of the results
    for file_ in os.listdir(path):
        if file_.startswith('out') and file_.endswith('.npy') and file_!='out_Max.npy':
            all_ = np.load(os.path.join(path, file_))
            all_[:,0] = all_[:,0]/maxValue

            # write
            np.save(os.path.join(path, file_), all_)
    return


# Append multiple value to a key in dictionary
def add_dict(_dict, key, list_of_values):
    """Append multiple values to a key in the given dictionary"""
    if key not in _dict:
        _dict[key] = list()
    _dict[key].extend(list_of_values)
    return _dict

