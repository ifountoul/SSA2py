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

from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
from obspy.core.event.origin import Origin
from obspy.core.event.magnitude import Magnitude
from obspy.core.event.catalog import Catalog
from obspy.core.event.event import Event

import numpy as np
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, Polygon, Feature
import math, os

from SSA2py.core import config

def getCatalog(text):

    def getEvent(info):
        """
        Converting text from format:
        2020-06-04T18:12:03,3.0,ML,37.24,20.49,4.1
        to ObsPy Event object 
        """
        evt=Event()
        evt.event_type='earthquake'
        evt.event_type_certainty='suspected'
        org=Origin()
        org.time=UTCDateTime(info[0])
        org.latitude=float(info[3])
        org.longitude=float(info[4])
        org.depth=float(info[5])*1000
        mag=Magnitude()
        mag.mag=float(info[1])
        mag.magnitude_type=info[2]
        mag.origin_id=org.resource_id
        evt.origins.append(org)
        evt.magnitudes.append(mag)
        evt.preferred_origin_id=org.resource_id
        evt.preferred_magnitude_id=mag.resource_id
        return evt

    cat=Catalog()
    cat.append(getEvent(text))
    return cat

def getFDSNWSCatalog(cfg, log, minMag=None, starttime=None, endtime=None, eventid=None):
    """
    # Returns a Catalog of Events based on starttime-endtime or eventid from a FDSNWS-event
    # Cases of Events requests (gisola API):
    # -e eventid (FDSNWS-event required)
    # -e datetime (FDSNWS-event required)
    # -e datetime lat lon mag depth
    # -p starttime endtime (FDSNWS-event required)
    # -r (-p (now-offset) - now)
    # -f eventfile where f: datetime lat lon mag depth (multiple lines)
    """
    log.info('Connecting to FDSNWS-event host:  ' + cfg['Download Service']['Event Info']['Host'])
    fdsnws = Client(cfg['Download Service']['Event Info']['Host'])

    if starttime and endtime:
        log.info('Requesting event in timespan: {} - {}'.format(UTCDateTime(starttime), UTCDateTime(endtime)))
        log.info('Minimum Magnitude and type is set to: '+str(minMag)+ \
        ' and '+ (cfg['Monitor']['Magnitudetype'] if cfg['Monitor']['Magnitudetype'] else 'None') +' respectively')

    elif eventid:
        log.info('Requesting event with ID: ' + eventid)

    try:
        if cfg['Download Service']['Event Info']['Host'] == 'https://esm-db.eu':
                return fdsnws.get_events(eventid=eventid or None,
                                         starttime=UTCDateTime(starttime) if starttime else None,
                                         endtime=UTCDateTime(endtime) if endtime else None,
                                         minmagnitude=minMag or None if not eventid else None)
        else:
                return fdsnws.get_events(eventid=eventid or None,
                                         starttime=UTCDateTime(starttime) if starttime else None,
                                         endtime=UTCDateTime(endtime) if endtime else None,
                                         minmagnitude=minMag or None if not eventid else None,
                                         magnitudetype=cfg['Monitor']['Magnitudetype'] or None if not eventid else None,
                                         includeallorigins=True, includeallmagnitudes=True, includearrivals=True,
                                         orderby='time-asc' if not eventid else None)

    except Exception as e:
        if e.__class__.__name__=='FDSNNoDataException':
            log.info('No event is found...')
        else:
            log.info(e)

        # in many cases Nodes dont support some features for example INGV hits
        # bad request "includeallorigins" parameter is allowed only with "eventid" or "originid".

        # try something else
        # now we will ingore includeallorigins, includeallmagnitudes, includearrivals and magnitudetype
        log.info('Lets try something else...')
        log.info('Ignore magnitudetype.')
        try:
            return fdsnws.get_events(eventid=eventid or None,
                                     starttime=UTCDateTime(starttime) if starttime else None,
                                     endtime=UTCDateTime(endtime) if endtime else None,
                                     minmagnitude=minMag or None if not eventid else None,
                                     orderby='time-asc' if not eventid else None)
        except Exception as e:
            log.info('The exception insists...')
            log.info(e)


def monitor(cfg, log):

    minMag = float(cfg['Monitor']['MagnitudeTresh'])
    # get the maximum duration of SSA
    dur = [np.abs(_[2][1]-_[2][0]) for _ in cfg['Backprojection']['Settings']['ScanningTime']]

    endtime= UTCDateTime.now()- math.ceil(min(dur))- 2 - cfg['Monitor']['Playback']
    starttime= endtime- math.ceil(max(dur)) - cfg['Monitor']['Range']

    #starttime='2022-02-01T00:00:00.00000Z' #set for debugging
    #endtime='2022-02-28T00:00:00.00000Z' #set for debugging

    cat=getFDSNWSCatalog(cfg, log, minMag=minMag, starttime=starttime, endtime=endtime)

    vcat=Catalog()
    for evt in (cat if cat else []):
        try:
            org=getOrigin(cfg,evt,cfg['Monitor']['Historical'])
            if (not cfg['Monitor']['Geobox']) or boolean_point_in_polygon(\
            Feature(geometry=Point((org.longitude, org.latitude))),\
            Polygon([[eval('('+x.replace(')','').replace('(','')+')') for x in cfg['Monitor']['Geobox'].replace(', ', ',').split('),(')]])):

                # check for quality of origin
                if qualityEvent(cfg, evt, org):
                    workdir=os.path.join(cfg['Events Dir'], str(org.time))
                    if not os.path.exists(workdir): 
                        vcat.append(evt)
                    else:
                        log.info('{}\nalready calculated. It is removed from process'.format(evt))
                else:
                    log.info('{}\ndoes not meet quality thresholds. It is removed from process'.format(evt))

        except Exception as e:
            log.info(e+'\nError in event parsing. Moving to next event, if any')

    log.info(vcat)
    return vcat    


def getMagnitude(evt, org):
    # retrieve associated magnitude
    mag = []
    for m in evt.magnitudes:
        if m.origin_id==org.resource_id:
            mag.append(m)
    if len(mag)!=0:
        return mag[0]
    elif len(mag)==0:
        for m in evt.magnitudes:
            mag.append(m)
        if len(mag)!=0:
            return mag[0]
    else:
        return []


def qualityEvent(cfg, evt, org):
    try:
        inittime=min([o.creation_info.creation_time for o in evt.origins])
        magnitude=getMagnitude(evt, org)

        # if data exit for realtime scenario
        windowRules=config.rules(magnitude.mag, cfg['Backprojection']['Settings']['ScanningTime'])

        maxtl=max([abs(tl[1]) for tl in windowRules])

        if org.time+math.ceil(maxtl)+2>UTCDateTime.now()-cfg['Monitor']['Playback']:
            return False

        if org.time+math.ceil(maxtl)+2+cfg['Monitor']['Quality']['Timeout']>=UTCDateTime.now()-cfg['Monitor']['Playback'] or \
           org.evaluation_mode=='manual' or \
           (round(org.time_errors.uncertainty,1)<=cfg['Monitor']['Quality']['Time'] and \
            round(org.depth_errors.uncertainty/1000,1)<=cfg['Monitor']['Quality']['Depth'] and \
            round(org.depth_errors.uncertainty/1000,1)!=0 and \
            round(org.latitude_errors.uncertainty,1)<=cfg['Monitor']['Quality']['Latitude'] and \
            round(org.longitude_errors.uncertainty,1)<=cfg['Monitor']['Quality']['Longitude'] and \
            round(magnitude.mag_errors.uncertainty,1)<=cfg['Monitor']['Quality']['Magnitude']):

            return True

    except Exception as e:
        print(e)
        return False


def getOrigin(cfg, evt, historical):

    if not historical:
        # retrieve "best" origin info (if any), else last found
        try:
            org=evt.preferred_origin()
            if not org: raise
        except:
            try:
                org=sorted(evt.origins, key=lambda o: o.creation_info.creation_time)[-1]
            except:
                org=evt.origins[-1]
    else:
        try:
            for i, org in enumerate(sorted(evt.origins, key=lambda o: o.creation_info.creation_time)):
                # first found (mimics real-time status)
                if qualityEvent(cfg,evt,org):
                    return org
        except:
            for i, org in enumerate(evt.origins):
                # first found (mimics real-time status)
                if qualityEvent(cfg,evt,org):
                    return org 
    return org

def getMT(cfg, evt):
    #retrieve moment tensor
    try:
        mt = [ evt.focal_mechanisms[0].moment_tensor.tensor.m_rr, evt.focal_mechanisms[0].moment_tensor.tensor.m_tt,\
             evt.focal_mechanisms[0].moment_tensor.tensor.m_pp, evt.focal_mechanisms[0].moment_tensor.tensor.m_rt,\
             evt.focal_mechanisms[0].moment_tensor.tensor.m_rp, evt.focal_mechanisms[0].moment_tensor.tensor.m_tp]
    except:
        mt = []
    return mt
