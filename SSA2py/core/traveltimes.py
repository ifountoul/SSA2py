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

import subprocess, sys, tempfile, shutil, multiprocessing,\
       os, obspy, sys, yaml, skfmm, math
import numpy as np
from array import array
from scipy.interpolate import interp1d
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from obspy.taup import TauPyModel
from obspy.taup.taup_create import build_taup_model
from obspy.geodetics.base import kilometer2degrees
from obspy.geodetics.base import gps2dist_azimuth

from shapely.geometry import shape, Point, Polygon
from pyproj import Transformer

# import local library
from SSA2py.core import config

def calculateTraveltimes(_cfg):
    """
    Calculate Traveltimes based on configuration
    """

    config.cfg = _cfg # give the config file
    # desicion Variables
    create_tables = False
    model=''
    type_=''

    # construct the directory
    if os.path.exists(_cfg['Traveltimes']['Save']):
        for f in os.listdir(_cfg['Traveltimes']['Save']):
            try:
                shutil.rmtree(os.path.join(_cfg['Traveltimes']['Save'], f))
            except Exception as e:
                os.remove(os.path.join(_cfg['Traveltimes']['Save'], f))
                pass
    else:
        os.makedirs(_cfg['Traveltimes']['Save'])

    # check initialy the 1D models
    if config.cfg['Traveltimes']['Priority']=='1D':
        # find the appropriate velocity model based on the location
        for crustal in config.cfg['Traveltimes']['Crustals1D']:
            if crustal['Geobox'] is None: #Take this model
                create_tables = True
                model=crustal; type_='1D';
            else:
                if geoboxCheck(crustal['Geobox'], config.org.longitude, config.org.latitude):
                    create_tables = True
                    model=crustal; type_='1D';
        # now check the 3D models
        if create_tables is False:
            for crustal in config.cfg['Traveltimes']['Crustals3D']:
                if crustal['Geobox'] is None: #Take this model
                    create_tables = True
                    model=crustal; type_='3D';
                else:
                    if geoboxCheck(crustal['Geobox'], config.org.longitude, config.org.latitude):
                        create_tables = True
                        model=crustal; type_='3D';
    else:
        for crustal in config.cfg['Traveltimes']['Crustals3D']:
            if crustal['Geobox'] is None: #Take this model
                create_tables = True
                model=crustal; type_='3D';
            else:
                if geoboxCheck(crustal['Geobox'], config.org.longitude, config.org.latitude):
                    create_tables = True
                    model=crustal; type_='3D';
        #Now check the 1D model
        if create_tables is False:
            for crustal in config.cfg['Traveltimes']['Crustals1D']:
                if crustal['Geobox'] is None: #Take this model
                    create_tables = True
                    model=crustal; type_='1D';
                else:
                    if geoboxCheck(crustal['Geobox'], config.org.longitude, config.org.latitude):
                        create_tables = True
                        model=crustal; type_='1D';  
    if create_tables:
        if type_=='1D':
            # traveltime save directory
            config.logger.info('Selected Velocity ' + type_+ ' Model --> ' + model['Filename'])
            tt_dir = os.path.join(config.cfg['Traveltimes']['Save'], '1D')
            os.makedirs(tt_dir)
            if config.cfg['Traveltimes']['Package'].upper() == 'NNLOC': #Run NNLOC package
                config.logger.info('Run NNLOC Package')
                config.vel2grid=os.path.join(config.cfg['System']['NonNinLoc'], 'Vel2Grid')
                config.grid2time=os.path.join(config.cfg['System']['NonNinLoc'], 'Grid2Time')
                for phase in config.cfg['Traveltimes']['Phase']:
                    for elevation in np.arange(0, model['Elevation'] + \
                                                model['Granularity'],
                                                model['Granularity']):
                        # create Control file for Nonlinloc
                        text=getControl(model['Filename'], phase, elevation,
                        model['Distance'], model['Depth'],model['Granularity'], model['VPVS'])

                        # generate command
                        command=("bash -c 'echo -e \"{}\" >./Control;{} ./Control;"
                        "{} ./Control'").format(text, config.vel2grid, config.grid2time)
                        # work on the tttables dir
                        proc=subprocess.Popen(command, cwd=tt_dir, shell=True,
                        universal_newlines=True)
                        proc.wait()

                        with open (os.path.join(tt_dir, 'time_layer.'+phase+'.1.time.hdr'), 'rb') as _:
                            lines = _.readlines()
                        nx, ny, nz=[int(x) for x in lines[0].split()[:3]]
                        # read buf data
                        b = array('f')
                        with open (os.path.join(tt_dir, 'time_layer.'+phase+'.1.time.buf'), 'rb') as _:
                             b.fromfile(_, nx*ny*nz)
                        data=np.array(b).reshape(nx, ny, nz)[0]
                       #Save the data to an numpy array
                        np.save(os.path.join(tt_dir,'model' + '_' + phase + '_' + str(elevation) + '.npy') , data)
                        #Also remove the nonninloc files
                        os.remove(os.path.join(tt_dir,'Control'))
                        for f in os.listdir(tt_dir):
                            if f.startswith('layer') or f.startswith('time'):
                                os.remove(os.path.join(tt_dir,f)) 
            if config.cfg['Traveltimes']['Package'].upper() == 'TAUP': #Run TAUP package
                config.logger.info('Run TAUP Package')
                customTAUP(model['Filename'], model['VPVS']) #Custom velocity model 
                #Build and read the model
                build_taup_model(os.path.join(tt_dir, 'Model.tvel'), tt_dir)
                model_ = TauPyModel(model=os.path.join(tt_dir,'Model.npz'))
                for phase in config.cfg['Traveltimes']['Phase']:
                    if phase == 'P':
                        phase_list=('P', 'p')
                    if phase == 'S':
                        phase_list=('S', 's')
                    for elevation in np.arange(0,model['Elevation'] + \
                                               model['Granularity'],
                                               model['Granularity']):
                        #Empty array
                        data = np.zeros((model['Distance'] + 1, int(model['Depth']/model['Granularity']) +\
                                         1 + int(elevation/model['Granularity'])))
                        counter_distance = 0
                        for distance in np.arange(0, model['Distance']+\
                                                  model['Granularity'],
                                                  model['Granularity']):
                            counter_depth = 0
                            for depth in np.arange(0, model['Depth']+\
                                                   model['Granularity'] + elevation,
                                                   model['Granularity']):
                                #Calculate arrivals
                                arrivals = model_.get_travel_times(source_depth_in_km = depth,\
                                           distance_in_degree=kilometer2degrees(distance), \
                                           phase_list=phase_list, receiver_depth_in_km= 0)
                                data[counter_distance, counter_depth] = arrivals[0].time
                                counter_depth+=1
                            counter_distance+=1
                        np.save(os.path.join(tt_dir,'model' + '_' + phase + '_' + str(elevation) + '.npy') , data) 
                #Remove the TAUP files
                os.remove(os.path.join(tt_dir, 'Model.npz'))
                os.remove(os.path.join(tt_dir, 'Model.tvel'))
            if config.cfg['Traveltimes']['Package'].upper() == 'FMM': #Run FMM package
                    config.logger.info('Run FMM Package')
                    for elevation in np.arange(0, model['Elevation'] + \
                                                model['Granularity'],
                                                model['Granularity']):
                       #Calculate
                       ttP, ttS = _1DFMM(model['Filename'], model['Distance'],\
                                         model['Depth'],model['Granularity'],\
                                         elevation, model['VPVS'])    
                       for phase in config.cfg['Traveltimes']['Phase']:
                           if phase == 'P':
                               #Save the data to an numpy array
                               np.save(os.path.join(tt_dir,'model' + '_' + phase + '_' + str(elevation) + '.npy') , ttP)
                           if phase == 'S':
                               np.save(os.path.join(tt_dir,'model' + '_' + phase + '_' + str(elevation) + '.npy') , ttS)
            #Update the dictionary
            model['type'] = type_
            model['method'] = config.cfg['Traveltimes']['Package'].upper()
           
            #Write the characteristics of the velocity model into a yaml file before the 'break'
            with open(os.path.join(config.cfg['Traveltimes']['Save'],'model' + '.yml'), 'w') as outfile:
                yaml.dump(crustal, outfile, default_flow_style=False)          
        if type_=='3D':
            #Create the traveltime tables for 3D model
            tt_dir = os.path.join(config.cfg['Traveltimes']['Save'], '3D')
            os.makedirs(tt_dir)
            config.logger.info('Run FMM Package')

            if 'P' in config.cfg['Traveltimes']['Phase']: 
                config.logger.info('Getting ' + type_+ ' Velocity Model for P phase --> ' + model['Filename VP'])

                modelP = _3DModelCrop(model['Filename VP'],\
                         [config.org.longitude, config.org.latitude, config.org.depth],\
                         model['Distance'])

            if 'S' in config.cfg['Traveltimes']['Phase']:
                if model['VPVS'] is None:
                    config.logger.info('Getting ' + type_+ ' Velocity Model for S phase --> ' + model['Filename VS'])

                    modelS = _3DModelCrop(model['Filename VS'],\
                         [config.org.longitude, config.org.latitude, config.org.depth],\
                         model['Distance'])

                else:
                    config.logger.info('Ignore S phase input model and use Vp/Vs= ' + str(model['VPVS']))
                    
                    modelS = _3DModelCrop(model['Filename VP'],\
                         [config.org.longitude, config.org.latitude, config.org.depth],\
                         model['Distance'])

                    modelS[:,3] = modelS[:,3]/model['VPVS']

            #Get the utm code
            utm_code = getEPSG(config.org.longitude, config.org.latitude)

            #Keep only the one component in the inventory
            invE = config.inv.copy().select(channel='*E')
            invZ = config.inv.copy().select(channel='*Z')
            invN = config.inv.copy().select(channel='*N') 
            inv = [invE, invZ, invN]

            inv = inv[np.argmax([len(invE), len(invZ), len(invN)])]

            #Create tt tables in parallel
            if 'P' in config.cfg['Traveltimes']['Phase']:
                with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
                    res = p.starmap(_3DCalc, list(zip(inv, [modelP]*len(inv), [model]*len(inv),\
                                   [utm_code]*len(inv), ['P']*len(inv), [tt_dir]*len(inv))))
            if 'S' in config.cfg['Traveltimes']['Phase']:
                with multiprocessing.Pool(processes=int(config.cfg['Backprojection']['NumCPUS'])) as p:
                    res = p.starmap(_3DCalc, list(zip(inv, [modelS]*len(inv), [model]*len(inv),\
                                   [utm_code]*len(inv), ['S']*len(inv), [tt_dir]*len(inv))))
            config.logger.info('3D Travel Time Tables Ready!') 
            #Update the dictionary
            model['type'] = type_
            model['method'] = config.cfg['Traveltimes']['Package'].upper()

            #Write the characteristics of the velocity model into a yaml file before the 'break'
            with open(os.path.join(config.cfg['Traveltimes']['Save'],'model' + '.yml'), 'w') as outfile:
                yaml.dump(crustal, outfile, default_flow_style=False)

    if not create_tables: #We cannot find velocity model / end the process
        config.logger.error('No Available Velocity Model for this Event. Exiting....')    
        sys.exit()

#Functions used
def _3DCalc(sta, model_, model, utm_code, phase, tt_dir):
    """
    Calculate the tt tables for stations. For parallel execution.

    Arguments:
    ------
    sta: Object 
        Obspy station Inventory Object
    model_: array-like
        Selected 3D Velocity model
    model: dictionary
        Characteristics of the selected model
    utm_code: int
        Utm Code for the transformation
    phase: str
        Phase to work with
    tt_dir: str
        Output Directory

    Returns:
    ------

    -

    """

    #Round to the nearest multiple of granularity
    elev_ = model['Granularity'] * round(abs(sta[0].elevation/1000)/model['Granularity'])

    #To help the computations give a stable granularity of the 3D model
    granHor = model['Granularity']; granVer = model['Granularity']  #in km

    #UTM code
    utm_code = getEPSG(sta[0].longitude, sta[0].latitude)

    #Move to UTM the model coords
    x_utmp, y_utmp = coords2utm(model_[:,0:2], utm_code, inverse='False')

    # Stack the model
    model_UTM = np.column_stack((x_utmp, y_utmp, model_[:,2], model_[:,3]))

    #get new grid
    gx, gy, gz = newGrid(config.org.longitude, config.org.latitude, model['Distance'],\
                 elev_, model['Depth'], granHor, granVer)

    #get the grid in UTM
    x_utmG, y_utmG = coords2utm(np.column_stack((gx.ravel(), gy.ravel())), utm_code, inverse='False')

    #reshape the grid
    x_utmG = x_utmG.reshape((gx.shape))
    y_utmG = y_utmG.reshape((gy.shape))

    # interpolate the grid
    V = interpGrid(model_, x_utmG/1000, y_utmG/1000, gz) 

    #Move to UTM station
    x_sta, y_sta = coords2utm(np.column_stack((sta[0].longitude, sta[0].latitude)), utm_code, inverse='False')

    #Compute the granularities in UTM
    granx = np.mean(np.diff((x_utmG/1000)[:,0,0], axis=0)) 
    grany = np.mean(np.diff((y_utmG/1000)[0,:,0], axis=0))
    granz = granVer

    #Calculate tt tables
    tt = _3DFMM(x_utmG/1000, y_utmG/1000, gz, V, [x_sta/1000, y_sta/1000, elev_], granx, grany, granz)   

    #Save table     
    np.save(os.path.join(tt_dir, sta[0].code + '_' + phase + '.npy'), tt)
    config.logger.info('TT Table for ' + sta[0].code + ' ' + phase + ' phase is ready')
 
    return


def newGrid(cLon, cLat, distance, elev_, depth, granHor, granVer):
    """
    Get the new resampled grid

    Arguments:
    ----------
    cLon: float
        Center Longitude
    cLat: float
        Center Latitude
    distance: float
        Distance
    elev_: float
        Elevation
    depth: float
        Depth
    granHor: float
        Horizontal granularity.
    granVer: float
        Vertical granularity.
    Returns:
    --------

    -

    """

    gx = np.linspace(cLon-kilometer2degrees(round(distance)),\
                     cLon+kilometer2degrees(round(distance)),\
                     round(distance*2/granHor)+1)
    gy = np.linspace(cLat-kilometer2degrees(round(distance)),\
                     cLat+kilometer2degrees(round(distance)),\
                     round(distance*2/granHor)+1)
    gz = np.arange(-elev_, (depth)+granVer, granVer)

    gx, gy, gz = np.meshgrid(gx, gy, gz, indexing = 'ij', sparse = False)
   

    return (gx, gy, gz)

def interpGrid(model, gx, gy, gz):
    """
    Interpolate the grid

    Arguments:
    ----------
    model: array-like
        3D velocity Model.
    gx: array-like
        Grid in X.
    gy: array-like
        Grid in Y.
    gz: array-like
        Grid in Z.
    Returns:
    --------

    -

    """
    #Move model to km
    model[:,0] = model[:,0]/1000;
    model[:,1] = model[:,1]/1000;

    #interpolate
    points = np.column_stack((model[:,0], model[:,1], model[:,2]))
    NInter= NearestNDInterpolator(points, model[:,3])

    #Interpolate to the new grid
    V = NInter(gx, gy, gz)

    return (V)

def _3DFMM(x, y, z, V, source, granx, grany, granz):
    """
    Compute traveltimes for 3D velocity model using the Fast Marching Method

    Arguments:
    -------
    x,y,z: array-like
        Coordinates of the grid in UTM and in km's
    V: array-like
        Velocity
    source: array-like
        Source position in [x,y,z]
    granHor: float
        Granularity of the grid in km (Horizontal)
    granVer: float
        Granularity of the grid in km (Vertical)

    Returns:
        tt: array-like
            Traveltimes
    """

    #Identify the closer grid point to the source
    idx = np.where(abs(x-source[0])==np.min(abs(x-source[0])))[0][0]
    idy = np.where(abs(y-source[1])==np.min(abs(y-source[1])))[1][0]
    idz = np.where(abs(z-source[2])==np.min(abs(z-source[2])))[2][0]

    #Create the phi array
    phi = -np.ones((x.shape[0], y.shape[1], z.shape[2]))
    #Put the source
    phi[idx,idy,idz] = 0
    #Calculate travel times
    tt = skfmm.travel_time(phi, V, dx=[granx, grany, granz])
    return tt


def getEPSG(lon,lat):
    """
    Based on lon, lat get the best utm epsg-code

    from --> https://gis.stackexchange.com/questions/269518/auto-select-suitable-utm-zone-based-on-grid-intersection

    Arguments:
    ------
    lon: float
       Longitude
    lat: float
       Latitude

    Returns:
    -------

    epsg_code: int
        epsg code for the coordinates pair

    """
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
        return epsg_code
    epsg_code = '327' + utm_band
    return epsg_code


def coords2utm(coords, epsg, inverse='False'):
    """
    Convert coordinates to UTM
   
    Arguments:
    ------
    coords: array-like
        Array with Lon, Lat information in columns
    epsg: int 
        EPSG code for transformation (if inverse=='True' input projection)
    inverse: bool
        If True from local to coordinates again

    Returns: array-like
    -------
        Array with x,y or lon,lat
    """

    lon = coords[:,0]; lat = coords[:,1];

    if not inverse:
        transformer = Transformer.from_crs("epsg:{}".format(epsg), "epsg:4326", always_xy=True)
        x,y = transformer.transform(lon, lat)
    else:
        transformer = Transformer.from_crs("epsg:4326", "epsg:{}".format(epsg), always_xy=True)
        x,y = np.array(transformer.transform(lon, lat))

    return x,y


def _3DModelCrop(path, event, dist):
    """
    Crop from the 3D model a specific part based on the distance from the event

    Arguments:
    ------
    path: str
        path of the 3D velocity model
    event: list
        [Lon, Lat, Depth] of the seismic event
    dist: float 
        Distance for travel time calculations

    Returns:
    ------
    model: list 
        Cropped 3D model [Lon, Lat, Depth, V]

    """

    #Read the model
    model = np.loadtxt(path)
    #Filter the model
    modelF = [ [point[1],point[0],point[2],point[3]] \
               for point in model if (gps2dist_azimuth(event[1], event[0], point[1], point[0])[0]/1000)<=dist]
    return np.array(modelF)

def _1DFMM(path, distance_g, depth_g, gran, elevation, VPVS=None):
    """
    Compute traveltime tables for 1D velocity model using the Fast Marching Method

    Arguments:
    -------
    path: str
        Path of the 1D velocity Model
    distance_g: float
        Distance of the grid (km)
    depth_g: float 
        Depth of the grid (km)
    gran: float
        Granularity or step of the grid points (km)
    elevation: float
        Elevation in km to start
    VPVS: float
        Create S model using VPVS
  
    Returns:
    -------
    (ttP, ttS): tuple
        Traveltime tables (numpy arrays) for P and S waves

    """

    model = np.loadtxt(path)
    Depth = model[:,0]; VP = model[:,1];

    if VPVS is not None:
        VS = VP/VPVS
    else:
        VS = model[:,2];

    #Interpolate the velocity model in depth
    fp = interp1d(Depth,VP, bounds_error=False, fill_value=(VP[0], VP[-1]))
    fs = interp1d(Depth,VS, bounds_error=False, fill_value=(VS[0], VS[-1]))

    #Setup the used grid
    gx = np.linspace(0, distance_g+1, int(distance_g/gran)+1)
    gz = np.linspace(-elevation, depth_g+1, int((depth_g+elevation)/gran) + 1)

    iVP = fp(gz); iVS = fs(gz);

    #Convert to 2D
    iVP = np.tile(iVP, (gx.shape[0], 1))
    iVS = np.tile(iVS, (gx.shape[0], 1))

    phi = -np.ones((gx.shape[0], gz.shape[0]))

    phi[0,0] = 0

    ttP = skfmm.travel_time(phi, iVP, dx=[gran, gran])
    ttS = skfmm.travel_time(phi, iVS, dx=[gran, gran])

    return (ttP, ttS)


def customTAUP(model_path, VPVS=None):
    """
    Create custom TAUP model using the input model and the IASP91 as the rest of the model

    Arguments:
    ------
    model_path: string
        Path of the 1D velocity model

    Returns:
    ------
    -

    """

    # read crustal layers from given file
    arr=np.loadtxt(model_path)

    DEPTH = arr[:,0]
    VP = arr[:,1]

    if VPVS is not None:
        VS = VP/VPVS
    else:
        VS = arr[:,2]

    # if only depth, vp and vs given, calculate rest needed values
    if arr.shape[1]==3:
        arr=modelCalc(DEPTH, VP, VS)
    #Read the IASP velocity model from Obspy
    path = os.path.abspath(obspy.__file__).split('/')[:-1] #Obspy path
    iaspPath = os.path.join('/'.join(path), 'taup/data/iasp91.npz')
 
    #Open the velocity model
    m_ = np.load(iaspPath) 
    model_ = np.zeros((m_['v_mod.layers'].size, 12))
 
    for i in range(m_['v_mod.layers'].size):
        model_[i,:] = np.stack(m_['v_mod.layers'][i])
    #keep depth, Vp, Vs, Density
    depths= model_[:,0]
    Vp = model_[:,2]; Vs = model_[:,4]
    density = model_[:,6]

    #interpolatete the arrays
    d_ = np.arange(min(arr[:,0]), max(arr[:,0])+1, 1)
    Vp_input = np.interp(d_, arr[:,0], arr[:,1])
    Vs_input = np.interp(d_, arr[:,0], arr[:,3])
    density_input = np.interp(d_, arr[:,0], arr[:,5])

    #Maximum position to get
    keep_maxPos = max(np.where(Vp>max(Vp_input))[0][0], np.where(Vs>max(Vs_input))[0][0])

    #Stack the IASP and the local Model
    VP = np.append(Vp_input, Vp[keep_maxPos:])
    VS = np.append(Vs_input, Vs[keep_maxPos:])
    DEN = np.append(density_input, density[keep_maxPos:])
    DEPTHS = np.append(d_, depths[keep_maxPos:])

    text = 'IASPEI91 P Model ({} layers) no "discontinuity" at 120, 760 km)\nIASPEI91 S Model ({} values to cmb)\n'\
       .format(len(DEPTHS), np.where(DEPTHS == 2889.000)[0][0])

    f = open(os.path.join(config.cfg['Traveltimes']['Save'], '1D', 'Model.tvel'), 'w+')
    f.write(text)
    for i in range(len(DEPTHS)):
        if DEPTHS[i] == 2889.0:
            f.write('{:12.3f}{:9.4f}{:9.4f}{:9.4f}\n'.format(2889.000, 13.6908, 7.3015, 5.5515))
        if DEPTHS[i] == 5153.900:
            f.write('{:12.3f}{:9.4f}{:9.4f}{:9.4f}\n'.format(5153.900,10.2578,0.0000,12.1391))
        f.write('{:12.3f}{:9.4f}{:9.4f}{:9.4f}\n'.format(DEPTHS[i], VP[i], VS[i], DEN[i]))
    return

def getControl(model_path, phase, elevation, distance, depth, step, VPVS=None):
  """
  Prepare the Control file needed for the pre-calculation

  Arguments:
  ------
  model_path: string
      Path of the 1D velocity model
  phase: string
      P or S phase
  elevation: float
      Elevation of computation (km)
  distance: float
      Maximum distance (km)
  depth: float
      Maximum depth (km)
  step: float
      Granularity of computation (km)

  Returns:
  ------
  string with the cotrol file for the NONINLOC

  """
  # calculate degrees from distance
  deg=kilometer2degrees(distance)

  # read crustal layers from given file
  arr=np.loadtxt(model_path)

  DEPTH = arr[:,0] 
  VP = arr[:,1]

  # if only depth, vp and vs given, calculate rest needed values
  if VPVS is not None:
      VS = VP/VPVS
  else:
      VS = arr[:,2]

  arr=modelCalc(DEPTH, VP, VS, 
      gradient=True)

  layers=''.join([('LAYER {:.2f} {:.3f} {:.3f} {:.3f} {:.3f} '
  '{:.3f} {:.3f}\n').format(*row) for row in arr])


  # create Control file format
  return ('CONTROL 0 54321\n'
  'TRANS  LAMBERT  WGS-84 37.98 23.72 {} {} 0.0\n'
  'VGOUT ./layer\n'
  'VGTYPE {}\n'
  'VGGRID 2 {} {} 0.0 0.0 {} {} {} {} SLOW_LEN\n'
  '{}\n'
  'GTFILES ./layer ./time_layer {}\n'
  'GTMODE GRID3D ANGLES_NO\n'
  'GTSRCE 1 LATLON 37.98 23.72 0.0 {}\n'
  'GT_PLFD  1.0e-3  0\n').format(str(round(37.98-deg,1)), 
  str(round(37.98+deg,1)), phase, str(int(distance/step)+1), 
  str(int((depth+elevation)/step)+1), str(-elevation), str(step), 
  str(step), str(step), layers, phase, str(elevation))


def modelCalc(depth, vp, vs, gradient=False):
   """
   Calculate basic values for the Velocity model.

   Arguments:
   ---------
   vp: numpy array
       P wave velocity
   vs: numpy Array
       S wave velocity
   depth: numpy array
       depth values (numpy Array)
   gradient: boolean
       True/False for constant or gradient velocity layers.

   Returns:
   ------
   numpy array with the model

   """
   # density
   rho = 1.7 + 0.2 * vp

   # init values
   vp_grad = np.zeros(depth.shape[0])
   vs_grad = np.zeros(depth.shape[0])
   rho_grad = np.zeros(depth.shape[0])

   if gradient:
      for d in range(depth.shape[0]-1):
         depth_diff = depth[d+1] - depth[d] 
         vs_diff = vs[d+1] - vs[d]
         vp_diff = vp[d+1] - vp[d] 
         rho_diff = rho[d+1] - rho[d]

         vp_grad[d] = vp_diff/depth_diff
         vs_grad[d] = vs_diff/depth_diff
         rho_grad[d] = rho_diff/depth_diff

   return np.concatenate((depth[np.newaxis].T,
                         vp[np.newaxis].T,
                         vp_grad[np.newaxis].T,
                         vs[np.newaxis].T,
                         vs_grad[np.newaxis].T,
                         rho[np.newaxis].T,
                         rho_grad[np.newaxis].T),
                         axis=1)


def geoboxCheck(geobox, lon, lat):
    """
    Check if the hypocenter is inside the given geobox
    
    Arguments:
    ------
    geobox: area (Check format in config file)
    lon: float
        Longitude 
    lat: float
        Latitude

    Returns:
    ------
    boolean

    """

    area = Polygon([[p[0], p[1]] for p in eval(geobox)])
    #Check if the hypo in the area
    point = Point(lon, lat)
    if area.contains(point):
        return True
    else:
        return False


if __name__ == "__main__":

    calculateTraveltimes(_cfg)

