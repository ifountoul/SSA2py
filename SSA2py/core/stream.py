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

import os, sys, logging, fnmatch, time,\
obspy, multiprocessing, importlib, shutil, copy
import numpy as np
from operator import itemgetter

from obspy.signal.filter import envelope
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core import read
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.core.utcdatetime import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from obspy.signal.rotate import rotate_ne_rt
from obspy.signal.trigger import classic_sta_lta

# import local library
from SSA2py.core import config, backprojection
from SSA2py.core.config import timer
from SSA2py.core.qcplots import CLIPPlot, SNRPlot, TIMEPlot
from SSA2py.core.tt_corrections import tt_exec
from SSA2py.core.modules.normalize import normalize
from SSA2py.core.modules.qc import checkSNR, checkCLIP, checkTIME
from SSA2py.core.modules.trigger import KurtosisRec, positive_derivative
from SSA2py.core.basic_f.other import createDir

import matplotlib.pyplot as plt
import matplotlib

def stream_process():
    """
    Process the available streams based on configuration.

    - Compare the data and metadata infos. If station is missing from metadata discard it also from the Stream.
    - Remove Responce. (Optional - Parallel)
    - Quality Control. (Optional - Parallel) 
    - Rotate in Radial-Transverse. (Optional - Parallel)
    - Filtering/Change Type/Normalization/Resample of the data. (Parallel)
    - Write Streams.

    """

    #Folders
    data_dir = config.eventdir; metadata_dir = config.inventorydir;
    #Event info
    evla = config.org.latitude; evlo = config.org.longitude; origin = config.org.time;

    #Check for differences between the inventory and stream 
    config.logger.info('Check for incosistencies between the stream and the inventory')
    config.st = check(config.st, config.inv) 

    config.raw = True #The data are still raw for now keep that to know

    if config.cfg['Streams']['Response']==True:    
            config.logger.info('Remove response from traces')
            with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
                res=p.map(correct, config.st)
            #Reconstruct Stream 
            config.st = Stream([tr for tr in res])
            del res
            config.raw = False

    #Trim data to selected range, in case of larger range pad
    config.st = Stream([tr.trim(UTCDateTime(origin) + float(config.cfg['Streams']['Duration'][0]),\
                        tr.stats.endtime, pad=True, nearest_sample=True, fill_value=tr.data[0]) for tr in config.st])
    config.st = Stream([tr.trim(tr.stats.starttime, UTCDateTime(origin) + float(config.cfg['Streams']['Duration'][1]),\
                        pad=True, nearest_sample=True, fill_value=tr.data[-1]) for tr in config.st])

    #Rotate
    if config.cfg['Streams']['Rotate'] == True:
        config.logger.info('Rotate NE --> RT')
        #Get the station names from stream
        rot_list = []
        stations_name = list(dict.fromkeys([tr.stats.station for tr in config.st])) #Remove duplicates
        #For every station get only the E,N component and the lat,lon
        for station in stations_name:
            temp_st = config.st.copy().select(station=station)
            #Keep only the horizontals
            st_hori = temp_st.select(component='E') + temp_st.select(component='N')
            lat_sta = [sta.latitude for net in config.inv for sta in net if sta.code==station]
            lon_sta = [sta.longitude for net in config.inv for sta in net if sta.code==station]
            #Append all this into one list
            if st_hori.count()==2:
                rot_list.append([st_hori ,evla, evlo, lat_sta[0],lon_sta[0]])
            elif st_hori.count()==4: #Two stations with the same name (possible one broadband one accelerometer)
                st_hori.sort(keys=['channel'])
                rot_list.append([st_hori[0:2], evla, evlo, lat_sta[0],lon_sta[0]])
                rot_list.append([st_hori[2:], evla, evlo, lat_sta[0],lon_sta[0]])
            else: #We cannot handle this situation
                pass
        #Go parallel and rotate the signals
        with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
            rt=p.map(rotateRT, rot_list) 

        #Add the new RT to the already existed stream
        for s in rt:
            config.st += s
        del rt

        #Add the RT components to inventory file (take the info from the E component)
        for tr in config.st.select(component='E'):
            tempInfo = config.inv.select(station=tr.stats.station, channel=tr.stats.channel).copy()

            #Create the new channels
            for comp in ['R', 'T']:
                config.inv = add2inv(config.st, config.inv, comp)

    typ = config.cfg['Streams']['Type']

    #Add info in the inventory if we will create the horizontal component
    if config.cfg['Streams']['Combine'] is True and typ in ['ENV', 'ABS']:
        config.inv = add2inv(config.st, config.inv, 'H')    

    #Create the processed data directory
    proc_path = createDir(os.path.join(data_dir, 'Processed_Data'))

    #Filter-type-normalizarion/resample (build the list)
    config.logger.info('Start basic processing including Filtering/New Type/Normalization')
    config.logger.info(config.st.__str__(extended=True))

    #Convert the stream to the given quantities
    st_convert = commonMetric(config.st, config.raw, config.cfg['Streams']['Quantity'])

    #Clear the config
    config.st = []

    if config.cfg['Streams']['Corrections'][0] is True:
        #Apply time shifts to the traces
        st_convert = tt_exec(st_convert.copy())
    #Write the converted/corrected traces
    st_convert.write(os.path.join(data_dir, '{}_{}.mseed'.format('streams', 'corrected')), format="MSEED", encoding="FLOAT64")

    for f in config.cfg['Streams']['Filter']:
        with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
            res = p.starmap(basicProc, list(zip(st_convert, [f]*len(st_convert), [typ]*len(st_convert))))
            
        #Create a temporal stream
        temp_st = Stream([tr for tr in res])

        #Merge the horizontal components here.
        if config.cfg['Streams']['Combine'] is True and typ in ['ENV', 'ABS']:
            config.logger.info('Combine horizontal components')
            if config.cfg['Streams']['Rotate'] == True:
                temp_st = combine_horizontal(temp_st, 'RT')
            else:
                temp_st = combine_horizontal(temp_st, 'NE')

        #Write them to the right folder
        for comp in ['Z', 'N', 'E', 'R', 'T', 'H']:
            wcomp = temp_st.select(component=comp)

            #For the last time check that all waveforms have the same length
            wcomp = len_check(wcomp)

            if wcomp.count()>0:
               config.logger.info('Saving processed waveforms in one mseed file --> '+\
               os.path.join(proc_path, '{}_{}_{}.mseed'.format(typ, str(float(f[0]))+'_'+str(float(f[1])), comp,)))
               wcomp.write(os.path.join(proc_path, '{}_{}_{}.mseed'.format(typ, str(float(f[0]))+'_'+str(float(f[1])), comp))\
                          ,format="MSEED", encoding="FLOAT64") 
        
def basicProc(tr, f, typ):
    """
    Basic processing function that contains 
    -Filtering 
    -Normalize
    -Change type of trace
    -Resampling.

    Arguments:
    -----
    tr: Obspy trace object
    f: array-like
       Bandpass filter range
    typ: string
       Signal type from 'ENV': Envelope, 'OBS': Observed, 'ABS': Absolute 
    
    Returns:
    --------
    tr: Obspy trace
        Processed trace
    """

    tr.detrend('demean')
    tr.detrend('linear')
    tr.taper(max_percentage=0.05, side='both')
    #Resample the trace
    if config.cfg['Streams']['Resample'][0] == True:
        tr.resample(config.cfg['Streams']['Resample'][1], window='hanning', no_filter=True, strict_length=False)

    #Filter data
    if not (f[0] == 0 and f[1] == 0):
        tr.filter("bandpass", freqmin = f[0], freqmax = f[1])

    if typ not in ['ENV', 'ABS', 'OBS', 'STALTA', 'KURT']:
        msg = ("Method has to be one of ('ENV', 'ABS', 'OBS', 'STALTA', 'KURT')")
        raise ValueError(msg)
        pass

    if typ == 'ENV':
        tr.data = envelope(tr.data)
    elif typ == 'ABS':
        tr.data[tr.data<0] = 0
    elif typ == 'STALTA':
  
        STA = config.cfg['Streams']['Type Parameters'][0]
        LTA = config.cfg['Streams']['Type Parameters'][1]

        df = tr.stats.sampling_rate
        tr.data = classic_sta_lta(tr.data, int(STA*df), int(LTA*df))
    elif typ == 'KURT':
        Kurtosis_window = config.cfg['Streams']['Type Parameters'][0]
        C = 1-(tr.stats.delta/Kurtosis_window)
        tr.data = KurtosisRec(tr.data, C)
        tr.data = positive_derivative(tr.data, tr.stats.delta)
    else:
        pass
        
    tr.taper(max_percentage=0.05, side='both')
    #Normalize data 
    if config.cfg['Streams']['Normalize'][0] == True:
        tr = normalize(tr, n_root_factor = config.cfg['Streams']['Normalize'][1], \
                        norm_type = config.cfg['Streams']['Normalize'][2])
    return tr

def rotateRT(rot_list):
    """
    -Rotate to Radial-Transverse system.

    Arguments:
    ----------
    rot_list: list
    [st: Obspy Stream with East-West components
    lat1, lon1, lat2, lon2: source-station coordinates]

    Returns:
    --------
    st: Obspy Stream with RT
    """

    #Assign the list into variables
    st, lat1, lon1, lat2, lon2 = itemgetter(0, 1, 2, 3, 4,)(rot_list)

    #Rotate
    ba = gps2dist_azimuth(lat1, lon1, lat2, lon2)[1] 
    rt = rotate_ne_rt(st[0].data, st[1].data, ba)

    #Change the data
    st[0].data = rt[0]; st[0].data = rt[1]

    #Change the component name
    st[0].stats.channel = st[0].stats.channel[0:2] + 'R'
    st[1].stats.channel = st[1].stats.channel[0:2] + 'T'

    return st

def combine_horizontal(st, hors):
   """
   -Combine horizontal components

   Arguments:
   ----------
   st: obspy stream
       Stream with data
   hors: str
       Horizontal components to merge RT or NE
  
   Returns:
   --------
   st: Obspy stream with H component

   """

   if hors=='NE':
       st_N = st.select(component='N')
       st_E = st.select(component='E')
 
       trH = Trace()
       for trN in st_N:
           trE = st_E.select(station=trN.stats.station)[0]
           trH = trE.copy()
           trH.data = np.sqrt(trE.data**2 + trN.data**2) 
           trH.stats.channel = trE.stats.channel[0:2] + 'H'
           st.append(trH)
   else:
       st_R = st.select(component='R')
       st_T = st.select(component='T')

       trH = Trace()
       for trR in st_R:
           trT = st_T.select(station=trR.stats.station)[0]
           trH = trT.copy()
           trH.data = np.sqrt(trR.data**2 + trT.data**2)
           trH.stats.channel = trT.stats.channel[0:2] + 'H'
           st.append(trH)
   return st


def correct(tr):
    """
    -Remove Response (This can be performed only if a StationXML is available!!)

    Arguments:
    ---------
    tr: Obspy trace Object
        Trace
    Returns:
    -------- 
    tr: Obspy corrected trace Object
        Trace
    """
    #Remove responce
    tr.detrend('linear')
    try:
        if config.cfg['Streams']['Response']==True:
            tr.remove_response(output='ACC',
                               pre_filt=(0.001, 0.005, 45, 50), # bandpass frqs (Hz)
                               zero_mean=True, # detrend(demean)
                               taper=True, # cos taper from f1 to f2 and from f3 to f4
                               taper_fraction=0.05, # percentage of tapering
                               water_level=60
                               )
    except:
         config.logger.warning('Unable to remove response for trace: ' + tr.stats.network+'.'+tr.stats.station+'.'+tr.stats.channel)
    return tr

def commonMetric(st, raw, q):
    """
    -Differentiate/integrate traces for common metric Velocity/Acceleration/Displacement
  
    Arguments:
    ----------
    st: Obspy Stream Object
        Traces
    raw: bool
        True or False (With or without response)
    q: str
        Type of quantity choose between 'ACC', 'VEL', 'DISP'

    Returns:
    -------- 
    st: Obspy Stream Object 
       Traces
    """
    #Just to avoid problems with integration
    st.detrend('demean')
    st.detrend('linear');
    st.taper(max_percentage=0.05, type='hann',)
 
    # Avoid drifting
    st.filter('bandpass', freqmin=0.3, freqmax=20.0) 

    if raw==False: #Data in ACC physical quantity (without response)
        if q=='VEL':
            st.integrate(method='cumtrapz')
        elif q=='DISP':
            st.integrate(method='cumtrapz')
            st.integrate(method='cumtrapz')
        else:
            pass
    else:
        #The data here are raw (with response)
        for tr in st:
            if q=='ACC':
                if fnmatch.fnmatch(tr.stats.channel,'?N?') or fnmatch.fnmatch(tr.stats.channel,'?G?'):
                    pass
                if fnmatch.fnmatch(tr.stats.channel,'?H?') or fnmatch.fnmatch(tr.stats.channel,'?L?'):
                    tr.differentiate(method='gradient')
            elif q=='VEL':
                if fnmatch.fnmatch(tr.stats.channel,'?N?') or fnmatch.fnmatch(tr.stats.channel,'?G?'):
                    tr.integrate(method='cumtrapz')     
                if fnmatch.fnmatch(tr.stats.channel,'?H?') or fnmatch.fnmatch(tr.stats.channel,'?L?'):
                    pass
            elif q=='DISP':
                if fnmatch.fnmatch(tr.stats.channel,'?N?') or fnmatch.fnmatch(tr.stats.channel,'?G?'):
                    tr.integrate(method='cumtrapz')
                    tr.integrate(method='cumtrapz')
                if fnmatch.fnmatch(tr.stats.channel,'?H?') or fnmatch.fnmatch(tr.stats.channel,'?L?'):
                    tr.integrate(method='cumtrapz')
            else:
                if fnmatch.fnmatch(tr.stats.channel,'?N?') or fnmatch.fnmatch(tr.stats.channel,'?G?'):
                    pass
                if fnmatch.fnmatch(tr.stats.channel,'?H?') or fnmatch.fnmatch(tr.stats.channel,'?L?'):
                    tr.differentiate(method='gradient')
    return st

def len_check(st):
    """
    Recheck that all the traces in the stream have the same length.

    """

    #Maximum start and minimum end
    if st.count()>0:
         max_start = max([tr.stats.starttime for tr in st]) 
         min_end = min([tr.stats.endtime for tr in st])
         st.trim(max_start, min_end)
    return st


def check(st, inv):
    """
    -Check for incosistencies between the stream and the input inventory (XML or YAML)

    Arguments:
    ----------
    st: Obspy Stream Object
        Traces in Stream
    inv: Obspy inventory object
        Metadata file
  
    Returns:
    -------
    st: Obspy Stream Object
        Obspy clean Object
    """

    #Get the station names from stream
    station_stream = [tr.stats.station for tr in st]    

    #Get the station names from stream
    station_stream = [tr.stats.station for tr in st]

    stations_XML = [sta.code for net in inv for sta in net] #Name of the stations in XML

    diff_sta = list(set(station_stream) - set(stations_XML)) #Stations missing from the XML
    #Remove those traces
    for station_name in diff_sta:
        st_ = st.select(station=station_name)
        for tr in st_:
            st.remove(tr)
            config.logger.warning(str(tr.stats.network + '.' + tr.stats.station + '.' + tr.stats.location +'.'  + tr.stats.channel) +\
                                  ' was removed due to inconsistencies')
    return st        


def clean():
    """
    Clean waveforms by:
    - Removing waveforms with data and meta-data inconsistency
    - Removing waveforms with gaps
    - Rotate waveforms/update inventory
    - Run quality tests and remove traces
    - Prioritize stations
    """

    def verify(tr):
        """
        Verify consistency of data and meta-data
        """
        try:
            tr.verify()
            return True
        except:
            return False 

    config.logger.info('Removing waveforms with inconsistency ' +\
    'between data and meta-data')
    config.st=Stream(list(filter(verify, config.st)))

    config.logger.info('Removing waveforms with gaps')
    gaps=config.st.get_gaps()
    config.st=Stream(list(filter(lambda tr: not bool([_ for _ in gaps \
    if _[0]==tr.stats.network and _[1]==tr.stats.station and \
    _[2]==tr.stats.location and _[3]==tr.stats.channel]), config.st)))

    #In case of Gaps left do merge and interpolate the missing values
    config.st.merge(method=0, fill_value='interpolate')

    

    #Rotate waveforms
    try:
        config.logger.info('Rotate to --> ZNE system')
        config.st.rotate(method="->ZNE", inventory=config.inv)

        #Update the inventory file for the misaligned components
        for comps in config.cfg['Download Rules']['Components']:
            for cha in config.inv.select(channel= '*' + comps[0]):
                cha[0][0].azimuth, cha[0][0].dip, cha[0][0].code = 10, -90, cha[0][0].code[:2] + 'Z'
            for cha in config.inv.select(channel= '*' + comps[1]):
                cha[0][0].azimuth, cha[0][0].dip, cha[0][0].code = 0, 0, cha[0][0].code[:2] + 'N'
            for cha in config.inv.select(channel= '*' + comps[2]):
                cha[0][0].azimuth, cha[0][0].dip, cha[0][0].code = 90, 0, cha[0][0].code[:2] + 'E'
    except Exception as e:
        config.logger.warning("Unable to Rotate traces, only ZNE will be used")
        pass

    #Re-check that we have only E,N,Z components
    config.st = config.st.select(component="E") + config.st.select(component="N") + config.st.select(component="Z")

    #Keep the SM station if other stations coexist with the SM
    rem_ = Stream([config.st.select(station = tr.stats.station, channel = tr.stats.channel)[0] for tr in config.st if \
                   len(config.st.select(station = tr.stats.station, component = tr.stats.channel[2]))>1])

    config.st = Stream(list(filter(lambda tr: not bool([ tr for r in rem_ if tr.stats.station==r.stats.station \
                   and tr.stats.channel==r.stats.channel and r.stats.channel[1] != 'N']), config.st))).sort(keys=['network', 'station'])

    # delete the rem_ parameter
    del rem_

    if len(config.cfg['Streams']['Quality Control'][0])>0:
        # create a list of Stream object, one per net,sta,loc,cha
        streams=list(set([(tr.stats.network, tr.stats.station, \
        tr.stats.location, tr.stats.channel) for tr in config.st]))

        #Sort the Stream
        streams = sorted(streams, key=lambda tup: (tup[1],tup[3])) 

        config.logger.info('Starting waveform checks. Checking for disturbances (based on configuration enabling).')

        with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
           config.Qual=p.map(checkD, streams)
 
        #Write the results to a special .txt file and do some plots
        txtpath = createDir(os.path.join(config.eventdir, 'QC'))

        #deep copy results
        Qual_copy = copy.deepcopy(config.Qual)
        streams_ = copy.deepcopy(streams)

        plotStreams = config.st.copy()

        for check in ['CLIP', 'SNR', 'TIME']:
            if check in config.cfg['Streams']['Quality Control'][0]:
                if check=='CLIP' and config.cfg['Plotting']['Plots']==True:
                    try: 
                        config.logger.info('Plot CLIP results')
                        CLIPPlot(txtpath, streams_, plotStreams, Qual_copy)
                    except Exception as e:
                        config.logger.warning('Error in CLIP plot')
                        pass

                if check=='SNR' and config.cfg['Plotting']['Plots']==True:
                    try: 
                        config.logger.info('Plot SNR results')
                        SNRPlot(txtpath, streams_, plotStreams, Qual_copy)
                    except Exception as e:
                        config.logger.warning('Error in SNR plot')
                        pass

                if check=='TIME' and config.cfg['Plotting']['Plots']==True:
                    try:
                        config.logger.info('Plot TIME results')
                        TIMEPlot(txtpath, streams_, plotStreams, Qual_copy)
                    except Exception as e:
                        config.logger.warning('Error in TIME plot')
                        pass

                #remove from stream and results
                index_with_issue = []
                for i in range(len(Qual_copy)):
                    try:
                        if Qual_copy[i][check][0]==False:
                            index_with_issue.append(i)
                    except Exception as e:
                        pass
                streams_ = [streams_[i] for i in range(len(Qual_copy)) if i not in index_with_issue]
                Qual_copy = [Qual_copy[i] for i in range(len(Qual_copy)) if i not in index_with_issue]

        del plotStreams

        #Remove traces (if at least one False in the quality check)
        for i in range(config.st.count()):
            keys = config.Qual[i].keys()
            for k in keys:
                if config.Qual[i][k][0]==False:
                    config.logger.info('Trace ' + '{}.{}.{}.{}'.format(streams[i][0], streams[i][1],\
                                   streams[i][2], streams[i][3]) + ' removed from ' + k + ' test')
                    config.st.remove(config.st.select(network=streams[i][0], station=streams[i][1],\
                                     location=streams[i][2], channel=streams[i][3]).copy()[0])
                    break
         
 
    return

def commandRemoveST():
    """
    -Remove stations/components from the stream given from the command window

    """
    
    if config.stations_status == 'r': #Remove these stations/components
        for sta in config.stations_:
            try:
                if len(sta.split('.')) == 2:
                    config.st.remove(config.st.select(station=sta.split('.')[0], component=sta.split('.')[1])[0])
                else:
                    tempST = config.st.select(station=sta)
                    for tr in tempST:
                        config.st.remove(tr)
            except:
                config.logger.info('The ' + sta + ' cannot be removed')    
    if config.stations_status == 'k': #Keep only these stations/components
        pass
        st = Stream()
        for sta in config.stations_:
            try:
                if len(sta.split('.')) == 2:
                    st+=config.st.select(station=sta.split('.')[0], component=sta.split('.')[1])
                else:
                    st+=config.st.select(station=sta)
            except:
                pass
        config.st = st

def checkD(stream):
    """
    -Check disturbances per trace.
  
    Arguments:
    ----------
    stream: tuple
        Tuple with traces info

    Returns:
    --------
    Qdec: dict
        dictionary with check results

    """
    #Make upper the list
    config.cfg['Streams']['Quality Control'][0] = [each_method.upper() for each_method in config.cfg['Streams']['Quality Control'][0]]

    #Dictionary with decision characteristics
    Qdec = {}

    if 'CLIP' in config.cfg['Streams']['Quality Control'][0] and (stream[3][1]!='N' or stream[3][1]!='G'):
        Qdec = checkCLIP(config.st.select(station=stream[1], channel=stream[3])[0], Qdec)
    if 'SNR' in config.cfg['Streams']['Quality Control'][0]:
        Qdec = checkSNR(config.st.select(station=stream[1], channel=stream[3])[0],\
                        Qdec, 'KURT', 5, 5)
    if 'TIME' in config.cfg['Streams']['Quality Control'][0]:
        Qdec = checkTIME(config.st.select(station=stream[1], channel=stream[3])[0],\
                         Qdec, 'KURT', 3)
    return Qdec

def add2inv(st, inv, comp):
    """
    Add to inventory.
    In case of new conponents in the stream e.g. R,T,H.

    Arguments:
    ----------
    st: Obspy stream object
        Data
    inv: Obspy inventory object
        Metadata
    comp: str
        Component to add.

    Returns:
    --------
    inv: New inventory with new component

    """

    for tr in config.st.select(component='E'):
        tempInfo = inv.select(station=tr.stats.station, channel=tr.stats.channel).copy()

        cha = Channel(code=tr.stats.channel[:-1]+comp, location_code=tr.stats.location,\
                      latitude=tempInfo[0][0].latitude, longitude=tempInfo[0][0].longitude,\
                      elevation=tempInfo[0][0].elevation, depth=0)
        sta = Station(code=tempInfo[0][0].code, latitude=tempInfo[0][0].latitude, longitude=tempInfo[0][0].longitude,\
                      elevation=tempInfo[0][0].elevation, creation_date=tempInfo[0][0].creation_date)
        net = Network(code=tempInfo[0].code, stations=[])
        sta.channels.append(cha)
        net.stations.append(sta)
        inv.networks.append(net)
    return inv
   


def StreamReady(path):
    """
    Read processed Stream, insert stations info and prepare them for backprojection.
    Also create stations dictionary

    Arguments:
    ------
    path: str
        MSEED path

    Returns:
    --------
    st: Obspy Stream Object
        Traces

    """

    if os.path.isfile(path):
        st = read(path)
        config.logger.info('Succesfully read the stream file: ' + os.path.basename(os.path.normpath(path)))
        config.logger.info('Remove stations based on distance')
        #Filter the Stream Based on distance
        #First the inventory
        inv = config.inv.select(latitude = config.org.latitude,\
                                longitude = config.org.longitude,\
                                channel="*"+config.comp,\
                                maxradius =  kilometer2degrees(config.cfg['Backprojection']['Selection']['Distance']))
        #Stations in the main stream
        staStr_  = [tr.stats.station for tr in st]
        #New Stream
        st_ = Stream([])
        #Stations dictionary
        stations = {}
 
        for tr in st:
            try:
                inv_ = inv.copy().select(station=tr.stats.station,\
                       channel=tr.stats.channel)
                out = gps2dist_azimuth(config.org.latitude,\
                                       config.org.longitude,\
                                       inv_[0][0].latitude,
                                       inv_[0][0].longitude)
                dist = out[0]/1000
                azim = out[1]

                #Add the distance and the station position info
                tr.stats.distance = round(dist, 3)
                tr.stats.stalon = inv_[0][0].longitude
                tr.stats.stalat = inv_[0][0].latitude
                tr.stats.staelev = inv_[0][0].elevation
                tr.stats.azim = round(azim,1)
                stations[inv_[0][0].code] = [inv_[0].code, inv_[0][0].code, inv_[0][0].longitude,\
                                         inv_[0][0].latitude, inv_[0][0].elevation]
                st_.append(tr.copy())
            except Exception as e:
                pass

        #Sort the stream and the stations dict based on station name
        stations_ = sorted(stations)
        stations = {key:stations[key] for key in stations_}
        st_.sort(keys=['station'])
        return (st_, stations)
    else:
        config.logger.warning('The MSEED path is wrong. Return empty Stream')
        return (Stream([]), {})

