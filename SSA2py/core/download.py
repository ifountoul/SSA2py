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

import yaml, obspy, time, os, urllib.request, shutil, zipfile
import numpy as np

# Obspy Imports
###############

from obspy import read_inventory
from obspy.core import read
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn import RoutingClient
from obspy.geodetics.base import kilometers2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy import UTCDateTime, Stream
from obspy.core.inventory import Inventory, Network, Station, Channel

from obspy.clients.filesystem.sds import Client as SDS
from obspy.clients.seedlink.basic_client import Client as SeedLink

# Local Imports
###############

from SSA2py.core import config


def getInventory():
    """
    Get inventory from stationXML file or FDSNWS
    """

    # find min and max distance of accepted rules
    # for defining inventory geobox selection
    mindist=min([rule[0] for rule in config.distRules])
    maxdist=max([rule[1] for rule in config.distRules])

    #Time window to retreive data
    origin_time = config.org.time
    starttime = origin_time + config.timerules[0] 
    endtime= origin_time + config.timerules[1]

    config.logger.info(('Retrieving inventory with working epochs from {}' + \
    ' to {} and stations between {} and {} km from event\'s ' + \
    'location').format(starttime, endtime, mindist, maxdist))

    inv=Inventory()
    _inv=Inventory()

    # use queue priority as priority host service
    for service in config.cfg['Download Service']['Inventory']:
        # loops will finish when no service is left for requesting data
        if service[0]=='FDSNWS':
            try:
                # connect to host
                if service[1] in ["iris-federator", "eida-routing"]:
                    fdsnws = RoutingClient(service[1], timeout=500)
                    config.logger.info('Connecting to ' +service[1]) 
                else:
                    fdsnws = Client(service[1], eida_token=service[2], timeout=500)
                    config.logger.info('Connecting to FDSNWS host ' +service[1])

                if service[1]=='https://esm-db.eu':
                    _inv=fdsnws.get_stations(starttime=starttime,
                         endtime=endtime, latitude=config.org.latitude, \
                         longitude=config.org.longitude, \
                         minradius=kilometers2degrees(mindist), \
                         maxradius=kilometers2degrees(maxdist), \
                         level='channel')
                else:
                    # download inventory from host based on config rules
                    _inv=fdsnws.get_stations(starttime=starttime,
                         endtime=endtime, latitude=config.org.latitude, \
                         longitude=config.org.longitude, \
                         minradius=kilometers2degrees(mindist), \
                         maxradius=kilometers2degrees(maxdist), \
                         level='response')
            except Exception as e:
                if e.__class__.__name__=='FDSNNoDataException':
                    config.logger.info('Could not retrieve inventory from FDSNWS host '+service[1])
                else:
                    config.logger.info(e)
                _inv=Inventory()
                pass

        elif service[0]=='StationXML':
            try:
                config.logger.info('Reading Inventory from ' +service[1])
                _inv=read_inventory(service[1], format='STATIONXML').select(starttime=starttime, \
                     endtime=endtime, latitude=config.org.latitude, longitude=config.org.longitude, \
                     minradius=kilometers2degrees(mindist), \
                     maxradius=kilometers2degrees(maxdist))
            except:
                config.logger.info('Could not load Inventory file '+service[1])
                _inv=Inventory()
                pass
        
        elif service[0]=='StationYAML':
            if os.path.exists(service[1]):
                config.logger.info('Reading Inventory from ' +service[1])
                try:
                    with open(service[1]) as outfile:
                        sta_all = yaml.load(outfile, Loader=yaml.FullLoader)
                        #Build an XML file
                        for sta in sta_all:
                            net = Network(code=sta_all[sta]['network'], stations=[], start_date=obspy.UTCDateTime('1970-01-01'))
                            sta_ = Station(code=sta_all[sta]['station'], latitude=sta_all[sta]['coords'][0], longitude=sta_all[sta]['coords'][1],\
                                           elevation=sta_all[sta]['coords'][2], creation_date=obspy.UTCDateTime('1970-01-01'))
                            for comp in sta_all[sta]['channels']:
                                cha = Channel(code=comp, location_code='', latitude=sta_all[sta]['coords'][0], longitude=sta_all[sta]['coords'][1],\
                                              elevation=sta_all[sta]['coords'][2], depth=0)
                                sta_.channels.append(cha)
                            net.stations.append(sta_)

                            if gps2dist_azimuth(sta_all[sta]['coords'][0], sta_all[sta]['coords'][1], config.org.latitude, config.org.longitude)[0]/1000<=maxdist and\
                               gps2dist_azimuth(sta_all[sta]['coords'][0], sta_all[sta]['coords'][1], config.org.latitude, config.org.longitude)[0]/1000>=mindist:
                                _inv.networks.append(net)
                except:
                    config.logger.info('Could not load Inventory file '+service[1])                
                    _inv=Inventory()
                    pass
            else:
                config.logger.warning('The StationYAML path does not exists. Please change the Path') 
                _inv=Inventory()

        # remove duplicates (keep only the first occurance)
        for name in list(set([sta.code for net in inv for sta in net])):
            _inv=_inv.remove(station=name)
        # append rest
        inv+=_inv
    _inv=inv
    # get only selected streams
    chas=[]

    # for each accepted channel get all possible component combinations
    for cha in [case for rule in config.distRules for case in rule[2]]:
        chas+=[cha+orient for comp in config.cfg['Download Rules']['Components'] \
        for orient in comp]
    chas=sorted(list(set(chas)))

    config.logger.info('Accepted channels: '+ str(chas))

    inv=Inventory()
    for cha in chas:
        inv+=_inv.select(channel=cha).copy()
        
    # keep unique station names
    _inv=Inventory()
    for name in list(set([sta.code for net in inv for sta in net])):
        _inv+=inv.select(station=name)

    # remove or keep stations based on input file
    if len(config.stationslist) == 0:
        config.inv = _inv
    else:

        _inv_=Inventory()

        for sta in config.stationslist:
            try:
                sta_net = (sta.split()[0]).split('.')[0]
                sta_name = (sta.split()[0]).split('.')[1]
                sta_type = (sta.split()[0]).split('.')[2]

                des = sta.split()[1]

                if des == 'True':
                    _inv_ += _inv.select(network=sta_net, station=sta_name, channel=sta_type+'?').copy()
                else:
                    pass
            except:
                pass

        config.inv = _inv_    
  
    return
 
def do_ESM_request(_inv, eventid, _dir_, pt, dt, token=None):
    """
    Do data request from ESM.

    Input:
    ------

    _inv: Obspy inventory object
        Inventory with stations information.
    eventid: str
        Event id (ESM or USGS or EMSC)
    _dir_: path
        Temporal directory to save requests
    pt: str
        Prosessing type (CV:data converted to cm/s^2 or MP:manually processed or AP: automatically processed)
    dt: str 
        Data type (ACC: acceleration, in cm/s^2, VEL: velocity, in cm/s, DISP: displacement, in cm)
    token: path 
        Path to token .txt file

    Output:
    -------

    st: Obspy stream Object
        Stream with waveforms

    """

    # Create directory to save
    if not os.path.exists(os.path.join(_dir_, 'tmp')):
        os.makedirs(os.path.join(_dir_, 'tmp'))
    else:
        shutil.rmtree(os.path.join(_dir_, 'tmp'))
        os.makedirs(os.path.join(_dir_, 'tmp'))

    #Start to built the request
    if eventid.startswith('us'):
        cat = 'USGS'
    if eventid.startswith('INT') or eventid.startswith('EMSC'):
        cat = 'ESM'
    if eventid.startswith('20'):
        cat = 'EMSC'

    # Data format
    _format_ = 'mseed'

    # Per station
    for net in _inv:
        for sta in net:
            station = sta.code

            if token==None:
                req = 'https://esm-db.eu/esmws/eventdata/1/query?eventid={}&catalog={}&station={}&format={}&processing-type={}&data-type={}'\
                      .format(eventid, cat, station, _format_, pt, dt)

                #Do the request
                urllib.request.urlretrieve(req, os.path.join(_dir_, 'tmp', station+'_query.zip'))
            else:
                #Need testing....
                req = 'https://esm-db.eu/esmws/eventdata/1/query?eventid={}&catalog={}&station={}&format={}&processing-type={}&data-type={}'\
                      .format(eventid, cat, station, _format_, pt, dt)

                curl = 'curl -X POST -F '
                mess = """"message={}" """.format(token)
                req_ = """"{}" """.format(req)
                end = '-o ' + os.path.join(_dir_, 'tmp', station+'_query.zip')

                reqc = curl + mess + req_ + end
                os.system(reqc)
    # Unzip the files
    for _zip_ in os.listdir(os.path.join(_dir_, 'tmp')):
        try:
            with zipfile.ZipFile(os.path.join(_dir_, 'tmp', _zip_), 'r') as zip_ref:
                zip_ref.extractall(os.path.join(_dir_, 'tmp'))
        except:
            pass

    # Read the mseed files
    try:
        st = read(os.path.join(_dir_, 'tmp', '*mseed'))
    except:
        st = Stream([])

    # Remove the temporal directory
    shutil.rmtree(os.path.join(_dir_, 'tmp'))
    return st
                         
def getWaveforms():
    """
    Retrieve Waveforms from SDS or FDSNWS
    """

    starttime = config.org.time + config.timerules[0]
    endtime=config.org.time + config.timerules[1]

    # wait if endtime has not come yet
    if UTCDateTime.now() <= endtime:
        time.sleep(endtime-UTCDateTime.now()+1)

    config.logger.info(('Retrieving waveforms based on inventory from {} '+ \
    'to {}').format(starttime, endtime))

    # create a list of (net,sta,loc,cha,start,end) 
    streams=[(net.code, sta.code, cha.location_code, cha.code, \
                UTCDateTime(starttime)-1, UTCDateTime(endtime)+1) \
                for net in config.inv for sta in net for cha in sta]

    # use queue priority as priority host service
    for service in config.cfg['Download Service']['Stream']:
        # loops will finish by two ways: 
        # i) no streams left for download 
        # ii) or no service left for requesting data

        # if no streams left for downloading, exit function
        if not streams: return

        if service[0]=='SeedLink':
            try:
                # connect to SeedLink host
                sl=SeedLink(service[1].split(':')[0], port=int(service[1].split(':')[1]))
                config.logger.info('Connecting to SeedLink host ' +service[1])

                # make one request for all streams
                multiselect = ','.join(list(set(["%s_%s:%s%s" % (s[0], s[1], s[2], s[3][:2]+'?') for s in streams])))
                config.st+=sl._multiselect_request(multiselect, UTCDateTime(starttime)-1, UTCDateTime(endtime)+1)

            except:
                config.logger.info('Could not retrieve data from SeedLink host '+service[1])
                pass

        elif service[0]=='SDS':
            try:
                # connect to SDS host
                sds=SDS(service[1])
                config.logger.info('Connecting to SDS directory ' +service[1])

                # make one request per quad (net,sta,loc,cha)
                for s in streams:
                    try:
                        config.st+=sds.get_waveforms(*s)
                    except:
                        continue
            except:
                config.logger.info('Could not find SDS directory '+service[1])
                pass

        elif service[0]=='FDSNWS':
            try:
                # connect to host
                if service[1] in ["iris-federator", "eida-routing"]:
                    fdsnws = RoutingClient(service[1], timeout=500)
                    config.logger.info('Connecting to ' +service[1])
                    config.st+=fdsnws.get_waveforms_bulk(streams)
                elif service[1] == 'https://esm-db.eu':
                    ESM_st=do_ESM_request(config.inv, config.eventid_ESM, config.eventdir, 'MP', 'ACC', service[2])
                    for tr in ESM_st: #to m/s2
                        tr.data=tr.data/100
                    config.st+=ESM_st 
                else:
                    fdsnws = Client(service[1], eida_token=service[2], timeout=500)
                    config.logger.info('Connecting to FDSNWS host ' +service[1])
                    # make one request for all streams
                    config.st+=fdsnws.get_waveforms_bulk(streams)
            except Exception as e:
                if e.__class__.__name__=='FDSNNoDataException':
                    config.logger.info('Could not retrieve data from FDSNWS host '+service[1])
                else:
                    config.logger.info(e)
                pass

        elif service[0]=='MSEED':
            if os.path.exists(service[1]):
                try:
                    config.logger.info('Reading MSEED from ' +service[1])
                    config.st+= read(service[1])
                except Exception as e:
                    config.logger.info(e) 
                    config.logger.warning('Unable to read the MSEED file!')
            else:
                config.logger.warning('The MSEED path does not exists. Please change the Path')

        # filter streams by those that have been downloaded
        streams=list(filter(lambda s: not bool([tr for tr in config.st \
        if s[1]==tr.stats.station]), streams))
       
    #detrend-demean-taper before trim
    config.st.detrend('linear');config.st.detrend('demean');
    config.st.taper(max_percentage=0.05, side='both')
 
    # trim to exactly start and end time
    if service[0]!='MSEED':
        config.st.trim(starttime, endtime)
    else:
        #Trim data to selected range, in case of larger range pad
        config.st = Stream([tr.trim(UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][0]),\
                            tr.stats.endtime, pad=True, nearest_sample=True, fill_value=tr.data[0]) for tr in config.st])
        config.st = Stream([tr.trim(tr.stats.starttime, UTCDateTime(config.org.time) + float(config.cfg['Streams']['Duration'][1]),\
                            pad=True, nearest_sample=True, fill_value=tr.data[-1]) for tr in config.st])

    # get traces with capable time length and number of samples
    #config.st=Stream([tr for tr in config.st \
    #if tr.stats.endtime-tr.stats.starttime >= endtime-starttime-2])

    # save waveforms
    config.logger.info(('Saving waveforms in one mseed file in ' + \
    '{}').format(os.path.join(config.eventdir,'streams.mseed')))

    # write to float64
    for tr in config.st:
        tr.data = np.float64(tr.data)

    try:
        config.st.write(os.path.join(config.eventdir,'streams.mseed'), \
        format='MSEED', encoding="FLOAT64")
    except Exception as error:
        config.logger.warning('Warning! Possible empty Stream.')
        config.logger.warning(error)
        

def getClear():
    """
    - Compare the initial downloaded Stream with the Inventory. Remove stations that don't exists in the Stream.
    Basically remove inactive stations in the current period of time.
    """

    #empty inventory
    _inv=Inventory()
    tempST = config.st.copy().merge(method=0, fill_value='interpolate')

    #remove traces that dont have info in metadata
    traces2go = []
    for j in range(tempST.count()):
        try:
            config.inv.select(station=tempST[j].stats.station,\
            channel=tempST[j].stats.channel).copy()[0]
        except:
            traces2go.append([tempST[j].stats.station, tempST[j].stats.channel])
    for remTr in traces2go:
        trace = tempST.select(station=remTr[0], channel=remTr[1])[0].copy()
        tempST.remove(trace)      

    for tr in tempST:
        _inv += config.inv.select(station=tr.stats.station, channel=tr.stats.channel)

    config.inv = _inv
    config.st = tempST

    try:
        #Write the Stream
        config.st.write(os.path.join(config.eventdir,'streams.mseed'), \
        format='MSEED')
    except Exception as error:
        config.logger.warning('Warning! Possible Empty Stream and inconsistencies with inventory.')
        config.logger.warning(error)

    #Write the inventory
    config.logger.info(('Saving inventory in one StationXML file in ' + \
    '{}').format(os.path.join(config.inventorydir,'inventory.xml')))
    _inv.write(os.path.join(config.inventorydir,'inventory.xml'), \
    format='StationXML')
        
    #Also attach the inventory to the stream
    config.st.attach_response(config.inv)

    return


