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
#########

import yaml, logging, functools, os, time as tm

# Obspy Imports
###############

from obspy.core.stream import Stream
from obspy import UTCDateTime, Inventory
import copy

# Local Imports
###############

from SSA2py.core import event

def read(filepath):
    with open(filepath, 'r') as _f:
        return yaml.load(_f, Loader=yaml.FullLoader)

def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)
    console = logging.StreamHandler()

    form=logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    form.converter = tm.gmtime

    handler.setFormatter(form)
    console.setFormatter(form)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    console.setLevel(level)

    logger.addHandler(handler)
    logger.addHandler(console)

    logger.propagate = False

    return logger

def readLists(f_):
    """
    Keep or remove stations.
    """
    #Check if list or Path
    if os.path.exists(f_):
        try:
            with open(f_) as f:
                lines = f.readlines()
            lines = [line.rstrip() for line in lines]
            return lines
        except:   
            return 
    else:
        return

def rules(magnitude, domain):
        return [elem[2] for elem in list(filter(lambda x: True \
        if round(x[0],1)<=round(magnitude,1) and round(magnitude,1)<=round(x[1],1) \
        else False, domain))]


def init(evnt, _cfg):
    """
    Initialize process based on event's info and configuration
    """

    # global variables (they can be seen by other python modules)
    global cfg, org, mag, eventsdir, eventdir, inventorydir,\
           procdir, restdir, outputdir, st, vel_types, acc_types, distRules, \
           timerules, logger, inv, stations_, stations_status, grid, tables,\
           comp, fi, q, model, st_Synth, mt, data_sit, phase, scanningRules, gridRules, stations,\
           savedir, origin_p, posx, posy, raw, job, Qual, exists, eventid_ESM, dist_weights, density_weights, all_weights, stationslist 

    # make configuration global variable
    cfg=copy.deepcopy(_cfg)
 
    #Get the event origin and magnitude
    org=event.getOrigin(cfg,evnt, cfg['Monitor']['Historical'])
    mag=event.getMagnitude(evnt, org)
    mt=event.getMT(cfg,evnt)

    st=Stream() # init
    st_Synth = Stream() #init
    inv=Inventory(networks=[]) #init
    grid = [] #init
    tables = [] #init 
    comp = []; fi = []; q = [] #init
    model = []; #init
    phase = 'P' #init
    stations = {} #init
    origin_p = 0 #init
    posx = [] #init
    posy = [] #init
    raw = 0 #init
    job = []
    Qual = [] #init
    exists = False
    dist_weights = []
    density_weights = []
    all_weights = []

    origin_time = org.time
    eventid_ESM = 0
    eventsdir = _cfg['Events Dir']
    eventdir = os.path.join(eventsdir, str(origin_time))
    inventorydir = os.path.join(eventsdir, str(origin_time), 'Inventory')
    savedir = '' 

    stations_ = []
    stations_status = 0 
    #data_sit = 'real'

    #Create the directories
    if not os.path.exists(eventdir):
        os.makedirs(eventdir)
    else:
        exists = True
    if not os.path.exists(inventorydir):
        os.makedirs(inventorydir)

    logger = setup_logger(str(UTCDateTime.now()), os.path.join(eventdir, 'log'))
     
    # broadband
    vel_types=['M/S', 'M/SEC', 'NM/S', 'NM/SEC', 'CM/S', 'CM/SEC', 'MM/S', 'MM/SEC']

    # strong motion
    acc_types=['M/S**2', 'M/(S**2)', 'M/SEC**2', 'M/(SEC**2)', 'M/S/S', \
    'NM/S**2', 'NM/(S**2)', 'NM/SEC**2', 'NM/(SEC**2)', 'CM/S**2', \
    'CM/(S**2)', 'CM/SEC**2', 'CM/(SEC**2)', 'MM/S**2', 'MM/(S**2)', \
    'MM/SEC**2', 'MM/(SEC**2)']

    #Download distance rules
    distRules = rules(mag.mag, cfg['Download Rules']['Distance'])
    scanningRules = rules(mag.mag, cfg['Backprojection']['Settings']['ScanningTime'])
    gridRules = rules(mag.mag, cfg['Backprojection']['Grid'])
    timerules = cfg['Download Rules']['Time'] 

    if cfg['Download Rules']['Stationslist'][0]: 
        stationslist = readLists(cfg['Download Rules']['Stationslist'][1])
    else:
        stationslist = []

# decorator function for logging
def log(func):
    def inner():
        func()
        logger.info(st.__str__(extended=True))
    return inner

# decorator function for timing
def timer(func):
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        tic = tm.perf_counter()
        value = func(*args, **kwargs)
        toc = tm.perf_counter()
        elapsed_time = toc - tic
        logger.info(f"Elapsed time: {elapsed_time:0.4f} seconds")
        return value
    return wrapper_timer
# function suppressing deprecation warnings
def warns():
    #from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
    #from shapely.errors import ShapelyDeprecationWarning
    import warnings

    #warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
    #warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
    #warnings.filterwarnings('ignore', category=ShapelyDeprecationWarning)
    warnings.filterwarnings('ignore')
    

    return







