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

import os, sys


def check(cfg, logger):
   """
   This function will check the configuration file for inputs that can cause errors.
   In most the of the problematic cases the program will end.

   Input:
   ------
   cfg: Configuration file

   Output:
   -------
   Message about the problematic input.

   """ 

   logger.info('Config file checks....')

   # Check that we have all the keywords
   if 'Version' not in cfg:
       logger.warning("'Version' in the config file is missing.")
   if 'Events Dir' not in cfg:
       logger.warning("'Events Dir' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Traveltimes' not in cfg:
       logger.warning("'Traveltimes' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Download Service' not in cfg:
       logger.warning("'Download Service' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'System' not in cfg:
       logger.warning("'System' in the config file is missing.")
       logger.info('Ignore this message if you are not using any external code.')
   if 'Download Rules' not in cfg:
       logger.warning("'Download Rules' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Streams' not in cfg:
       logger.warning("'Streams' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Backprojection' not in cfg:
       logger.warning("'Backprojection' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Tests' not in cfg:
       logger.warning("'Tests' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Plotting' not in cfg:
       logger.warning("'Plotting' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if 'Monitor' not in cfg:
       logger.warning("'Monitor' in the config file is missing.")
       logger.info('Ignore this message if you are not using SSA2py in realtime.')

   if cfg['Version']!=1.0:
      logger.warning("'Version' in the config file should be 1.0")
   if not os.path.isdir(cfg['Events Dir']):
      logger.error("'Events Dir' should be directory!")
      logger.info('Program is terminated...')
      sys.exit()
    
   # Check if important components are missing from traveltimes
   if 'Crustals1D' not in cfg['Traveltimes']:
      logger.warning("'Crustals1D' in the config file is missing.")
      logger.info('Program is terminated...')
      sys.exit()
   if 'Crustals3D' not in cfg['Traveltimes']:
      logger.warning("'Crustals3D' in the config file is missing.")
      logger.info('Ignore this message if you are not using 3D model.')
   if 'Priority' not in cfg['Traveltimes']:
      logger.warning("'Priority' in the config file is missing.")
      logger.info('Program is terminated...')
      sys.exit()
   if 'Phase' not in cfg['Traveltimes']:
      logger.warning("'Phase' in the config file is missing.")
      logger.info('Program is terminated...')
      sys.exit()
   if 'Package' not in cfg['Traveltimes']:
      logger.warning("'Package' in the config file is missing.")
      logger.info('Program is terminated...')
      sys.exit()
   if 'Save' not in cfg['Traveltimes']:
      logger.warning("'Save' in the config file is missing.")
      logger.info('Program is terminated...')
      sys.exit()
   
   if not os.path.isdir(cfg['Traveltimes']['Save']):
      logger.error("'Traveltimes/Save' should be directory!")
      logger.info('Program is terminated...')
      sys.exit() 

   if bool(cfg['Traveltimes']['Package']):
      if cfg['Traveltimes']['Package'] not in ['NNLOC', 'TAUP', 'FMM']:
          logger.error("'Traveltimes/Package' should be NNLOC/TAUP/FMM!")
          logger.info('Taking FMM as default')
          cfg['Traveltimes']['Package'] = 'FMM'
   else:
       logger.info('Taking FMM as default')
       cfg['Traveltimes']['Package'] = 'FMM' 
 
   if bool(cfg['Traveltimes']['Priority']):
      if cfg['Traveltimes']['Priority'] not in ['1D', '3D']:
          logger.warning("'Traveltimes/Priority' should be 1D or 3D!")
          logger.info('Taking 1D as default')
          cfg['Traveltimes']['Priority']='1D'
   else:
      logger.error("No value in 'Traveltimes/Priority'")
      logger.info('Program is terminated...')
      sys.exit()  

   if bool(cfg['Traveltimes']['Phase']):
      for i in cfg['Traveltimes']['Phase']:
         if i == 'P' or i == 'S':
            pass
         else:
            logger.error("'Traveltimes/Phase' should be P,S not " + str(i))
            logger.info('Program is terminated...')
            sys.exit()
   else:
      logger.error("No value in 'Traveltimes/Phase'")
      logger.info('Program is terminated...')
      sys.exit() 


   for model in cfg['Traveltimes']['Crustals1D']:
       if 'Filename' not in model or 'Elevation' not in model or 'Depth' not in model or 'Distance' not in model\
       or 'Granularity' not in model or 'VPVS' not in model or 'Geobox' not in model:
           logger.error("The 'Crustals1D' should have the following variables: 'Filename', 'Elevation', 'Depth', 'Distance', 'Granularity', 'VPVS', 'Geobox'")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['Filename']):
           if not os.path.isfile(model['Filename']):
               logger.warning("Problem with path in Crustals1D " + model['Filename'])
               logger.info('If this model is not the selected one ignore this message')
       else:
           logger.error("Empty path in Crustals1D")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['Elevation']):
           if isinstance(model['Elevation'], str):
               logger.error('Elevation in Crustals1D models must be number not string!')
               logger.info('Program is terminated...')
               sys.exit()
           else:
               if model['Elevation']<0 or model['Elevation']>=10:
                   logger.warning('Elevation in Crustals1D model has extremely low/high value.')
                   logger.info('Use default Elevation: 2')
                   model['Elevation'] = 2
       else:
           logger.error("Empty parameter in Crustals1D/Elevation")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['Depth']):
           if isinstance(model['Depth'], str):
               logger.error('Depth in Crustals1D models must be number not string!')
               logger.info('Program is terminated...')
               sys.exit()
           else:
               if model['Depth']<0 or model['Depth']>=300:
                   logger.warning('Depth in Crustals1D model has extremely low/high value.')
                   logger.info('Use default Depth: 80')
                   model['Depth'] = 80
       else:
           logger.error("Empty parameter in Crustals1D/Depth")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['Distance']):
           if isinstance(model['Distance'], str):
               logger.error('Distance in Crustals1D models must be number not string!')
               logger.info('Program is terminated...')
               sys.exit()
           else:
               if model['Distance']<=0 or model['Distance']>=800:
                   logger.warning('Distance in Crustals1D model has extremely low/high value.')
                   logger.info('Use default Distance: 100')
                   model['Distance'] = 100
       else:
           logger.error("Empty parameter in Crustals1D/Distance")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['Granularity']):
           if isinstance(model['Granularity'], str):
               logger.error('Granularity in Crustals1D models must be number not string!')
               logger.info('Program is terminated...')
               sys.exit()
           else:
               if model['Granularity']<=0 or model['Granularity']>=100:
                   logger.warning('Granularity in Crustals1D model has extremely low/high value.')
                   logger.info('Use default Granularity: 1')
                   model['Granularity'] = 1
       else:
           logger.error("Empty parameter in Crustals1D/Granularity")
           logger.info('Program is terminated...')
           sys.exit()

       if bool(model['VPVS']):
           if isinstance(model['VPVS'], str) and model['VPVS']!='null':
               logger.error('VPVS in Crustals1D models must be number or null not string!')
               logger.info('Program is terminated...')
               sys.exit()


   if 'Crustals3D' in cfg['Traveltimes']:
       for model in cfg['Traveltimes']['Crustals3D']:
           if 'Filename VP' not in model or 'Filename VS' not in model or 'Depth' not in model or 'Distance' not in model\
           or 'Granularity' not in model or 'VPVS' not in model or 'Geobox' not in model:
              logger.error("The 'Crustals3D' should have the following variables: 'Filename VP', 'Filename VS', 'Elevation', 'Depth', 'Distance', 'Granularity', 'VPVS', 'Geobox'")
              logger.info('Program is terminated...')
              sys.exit()

           if bool(model['Filename VP']):
               if not os.path.isfile(model['Filename VP']):
                   logger.warning("Problem with path in Crustals3D " + os.path.isfile(model['Filename VP']))
           else:
               logger.warning("'Filename VP' is empty")
          
           if bool(model['Filename VS']):
               if not os.path.isfile(model['Filename VS']):
                   logger.error("Problem with path in Crustals3D " + os.path.isfile(model['Filename VS']))  
           else:
               logger.warning("'Filename VS' is empty")

           if bool(model['Depth']):
               if isinstance(model['Depth'], str):
                  logger.error('Depth in Crustals1D models must be number not string!')
                  logger.info('Program is terminated...')
                  sys.exit()
               else:
                  if model['Depth']<0 or model['Depth']>=300:
                      logger.warning('Depth in Crustals3D model has extremely low/high value.')
                      logger.info('Use default Depth: 80')
                      model['Depth'] = 80
                  
           else:
              logger.error("Empty parameter in Crustals1D/Depth")
              logger.info('Program is terminated...')
              sys.exit()

           if bool(model['Distance']):
               if isinstance(model['Distance'], str):
                  logger.error('Distance in Crustals1D models must be number not string!')
                  logger.info('Program is terminated...')
                  sys.exit()
               else:
                  if model['Distance']<=0 or model['Distance']>=800:
                     logger.warning('Distance in Crustals1D model has extremely low/high value.')
                     logger.info('Use default Distance: 100')
                     model['Distance'] = 100
           else:
              logger.error("Empty parameter in Crustals1D/Distance")
              logger.info('Program is terminated...')
              sys.exit()

           if bool(model['Granularity']):
               if isinstance(model['Granularity'], str):
                   logger.error('Granularity in Crustals1D models must be number not string!')
                   logger.info('Program is terminated...')
                   sys.exit()
               else:
                   if model['Granularity']<=0 or model['Granularity']>=100:
                      logger.warning('Granularity in Crustals1D model has extremely low/high value.')
                      logger.info('Use default Granularity: 1')
                      model['Granularity'] = 1
           else:
              logger.error("Empty parameter in Crustals1D/Granularity")
              logger.info('Program is terminated...')
              sys.exit()

           if bool(model['VPVS']):
               if isinstance(model['VPVS'], str) and model['VPVS']!='null':
                  logger.error('VPVS in Crustals1D models must be number or null not string!')
                  logger.info('Program is terminated...')
                  sys.exit()
   # > SYSTEM
   ############## 
   if bool(cfg['System']):
       if 'NonNinLoc' in cfg['System']:
          pass
       else:
          logger.warning("Only 'NonNinLoc' input is accepted in System")
          logger.info('Program is terminated...')
          sys.exit()

   # > Download Service
   ####################
   if 'Event Info' not in cfg['Download Service']:
       logger.warning("'Event Info' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Host' not in cfg['Download Service']['Event Info']:
       logger.warning("'Event Info/Host' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if not bool(cfg['Download Service']['Event Info']['Host']):
       logger.warning("'Event Info/Host' is empty.")

   if 'Inventory' not in cfg['Download Service']:
       logger.warning("'Inventory' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Download Service']['Inventory']):
       for inv in cfg['Download Service']['Inventory']:
           if len(inv)!=3:
               logger.error('Inventory lists always accepts 3 variables.')
               logger.info('Program is terminated...')
               sys.exit()
           if inv[0] not in ['FDSNWS', 'StationXML', 'StationYAML']:
               logger.error("Inventory input should be from 'FDSNWS', 'StationXML', 'StationYAML' not " + str(inv[0]))
               logger.info('Program is terminated...')
               sys.exit()
   else:
       logger.error("Empty inventory")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Stream' not in cfg['Download Service']:
       logger.warning("'Stream' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Download Service']['Stream']):
       for st in cfg['Download Service']['Stream']:
           if len(st)!=3:
               logger.error('Stream lists always accepts 3 variables.')
               logger.info('Program is terminated...')
               sys.exit()
           if st[0] not in ['SeedLink', 'SDS', 'FDSNWS', 'MSEED']:
               logger.error("Stream input should be from 'SeedLink', 'SDS', 'FDSNWS', 'MSEED' not " + str(st[0]))
               logger.info('Program is terminated...')
               sys.exit()
   else:
       logger.error("Empty Stream lists.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Distance' not in cfg['Download Rules']:
       logger.warning("'Distance' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Download Rules']['Distance']):
       for i in cfg['Download Rules']['Distance']:
           if len(i)!=3:
               logger.error('Download Rules/Distance always accepts 3 variables.')
               logger.info('Program is terminated...')
               sys.exit()
           if i[0]>=i[1]:
               logger.error('Download Rules/Distance problem with given magnitudes.')
               logger.info('Program is terminated...')
               sys.exit()
           if len(i[2])!=3:
               logger.error('Download Rules/Distance always accepts 3 variables.')
               logger.info('Program is terminated...')
               sys.exit() 
           if i[2][0]>=i[2][1]:
               logger.error('Download Rules/Distance problem with given distances.')
               logger.info('Program is terminated...')
               sys.exit()
           for j in i[2][2]:
               if j!='HH' and j!='HN' and j!='BN' and j!='BH' and j!='EH' and j!='EN':
                   logger.error('Download Rules/Distance problem with data identity. Accepts only IDS HH: Broadband, HN: Strong motion, BH: Broadband Borehole, BN: Strong motion Borehole.')
                   logger.info('Program is terminated...')
                   sys.exit()
   else:
       logger.error("Download Rules/Distance")
       logger.info('Program is terminated...')
       sys.exit()    
   
   if 'Time' not in cfg['Download Rules']:
       logger.warning("'Time' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Download Rules']['Time']):
       if len(cfg['Download Rules']['Time'])!=2:
           logger.error('Download Rules/Time always accepts 2 variables')
           logger.info('Program is terminated...')
           sys.exit()
       else:
           if cfg['Download Rules']['Time'][0]>=cfg['Download Rules']['Time'][1]:
               logger.error('Download Rules/Time problematic values')
               logger.info('Program is terminated...')
               sys.exit()
   else:
      logger.error("Download Rules/Time empty")
      logger.info('Program is terminated...')
      sys.exit()

   if 'Components' not in cfg['Download Rules']:
       logger.warning("'Components' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Download Rules']['Components']):
       if len(cfg['Download Rules']['Components'])==0:
           logger.error("Download Rules/Components empty list. Accept default values ['ZNE']")
           cfg['Download Rules']['Components'] = ['ZNE']
       else:
           for v in cfg['Download Rules']['Components']:
               if v not in ['Z23', 'Z12', '123', 'ZNE']:
                   logger.error("Download Rules/Components must be from Z23', 'Z12', '123', 'ZNE' not " + str(st[0]))
                   logger.info('Program is terminated...')
                   sys.exit()
   else:
       logger.error("Download Rules/Components empty. Accept default values ['ZNE']")
       cfg['Download Rules']['Components'] = ['ZNE']

   if 'Stationslist' not in cfg['Download Rules']:
       logger.warning("'Stationslist' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   if bool(cfg['Download Rules']['Stationslist']):
       if len(cfg['Download Rules']['Stationslist'])==0:
           logger.error("Download Rules/Stationslist is empty list. Accept default values [False, null]")
           cfg['Download Rules']['Stationslist'] = [False, 'null']
       else:
           if cfg['Download Rules']['Stationslist'][0]==True and len(cfg['Download Rules']['Stationslist'])==2:
               if not os.path.exists(cfg['Download Rules']['Stationslist'][1]):
                   logger.error('Stations file does not exist in Download Rules/Stationslist.')
                   cfg['Download Rules']['Stationslist'] = [False, 'null']
           else:
               cfg['Download Rules']['Stationslist'] = [False, 'null']
   
   # > Streams
   ############

   if 'Duration' not in cfg['Streams']:
       logger.error("'Streams/Duration' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Quality Control' not in cfg['Streams']:
       logger.error("'Streams/Quality Control' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Resample' not in cfg['Streams']:
       logger.error("'Streams/Resample' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Filter' not in cfg['Streams']:
       logger.error("'Streams/Filter' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Type' not in cfg['Streams']:
       logger.error("'Streams/Type' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Type Parameters' not in cfg['Streams']:
       logger.error("'Streams/Type Parameters' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Rotate' not in cfg['Streams']:
       logger.error("'Streams/Rotate' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Response' not in cfg['Streams']:
       logger.error("'Streams/Response' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Quantity' not in cfg['Streams']:
       logger.error("'Streams/Quantity' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Normalize' not in cfg['Streams']:
       logger.error("'Streams/Normalize' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Corrections' not in cfg['Streams']:
       logger.error("'Streams/Corrections' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Combine' not in cfg['Streams']:
       logger.error("'Streams/Combine' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Pre Filter' not in cfg['Streams']:
       logger.error("'Streams/Pre Filter' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if bool(cfg['Streams']['Duration']):
       if len(cfg['Streams']['Duration'])==2:
           if cfg['Streams']['Duration'][0]>0 or cfg['Streams']['Duration'][0]>cfg['Streams']['Duration'][1] or cfg['Streams']['Duration'][1]<=0:
               logger.warning("Streams/Duration starttime is problematic. Use as Default [-10, 100]")
               cfg['Streams']['Duration'] = [-10, 100]
       else:
           logger.warning("Streams/Duration is problematic. Use as Default [-10, 100]")
           cfg['Streams']['Duration'] = [-10, 100] 
   else:
       logger.warning("Streams/Duration is emtpy. Use as Default [-10, 100]")
       cfg['Streams']['Duration'] = [-10, 100]
 
   if bool(cfg['Streams']['Resample']):
       if len(cfg['Streams']['Resample'])==2:
           if not str(cfg['Streams']['Resample'][0])=='True' and not str(cfg['Streams']['Resample'][0])=='False':
               logger.warning("Streams/Resampling should be True/False. Use as Default True")
               cfg['Streams']['Resample'][0] = True
           if isinstance(cfg['Streams']['Resample'][1], str):
               logger.warning("Streams/Resampling position 2 should be Float. Use as Default 100")
               cfg['Streams']['Resample'][1] = 100
           if cfg['Streams']['Resample'][1]==0 or cfg['Streams']['Resample'][1]>1000:
               logger.warning("Streams/Resampling position 2 is too low/large. Use as Default 100")
               cfg['Streams']['Resample'][1] = 100
       else:
           logger.warning("Streams/Resample is problematic. Use as Default [True, 100]")
           cfg['Streams']['Resample'] = [True, 100]
   else:
       logger.warning("Streams/Resample is empty. Use as Default [True, 100]")
       cfg['Streams']['Resample'] = [True, 100]      
  
   if bool(cfg['Streams']['Type']):
       if not isinstance(cfg['Streams']['Type'], str):
           logger.warning("Streams/Type should be string 'ENV': Envelopes, 'OBS': Observed waveforms,'ABS' : Absolute part of the Observed (>1), 'STALTA': STALTA method, 'KURT': Kurtosis method")
           logger.warning("Use as Default ENV")
           cfg['Streams']['Type'] = 'ENV'   
       if cfg['Streams']['Type'] not in ['ENV', 'OBS', 'ABS', 'STALTA', 'KURT']:
           logger.warning("Streams/Type should be string 'ENV': Envelopes, 'OBS': Observed waveforms,'ABS' : Absolute part of the Observed (>1), 'STALTA': STALTA method, 'KURT': Kurtosis method")  
           logger.warning("Use as Default ENV")
           cfg['Streams']['Type'] = 'ENV'
   else: 
      logger.warning("Streams/Type is empty. Use as Default ENV")
      cfg['Streams']['Type'] = 'ENV'

   if bool(cfg['Streams']['Type Parameters']):
      if cfg['Streams']['Type'] == 'KURT':
          if len(cfg['Streams']['Type Parameters'])==1:
              if not isinstance(cfg['Streams']['Type'][0], (float, int)):
                  logger.warning("Streams/Type Parameters should be float or integer. Use as default [0.5]")
                  cfg['Streams']['Type Parameters'] = [0.5]
          else:
              logger.warning("Streams/Type Parameters should have 1 value list. Use as default [0.5]")
              cfg['Streams']['Type Parameters'] = [0.5]
      if cfg['Streams']['Type'] == 'STALTA':
          if len(cfg['Streams']['Type Parameters'])==2:
              if not isinstance(cfg['Streams']['Type'][0], (float, int)) or not isinstance(cfg['Streams']['Type'][1], (float, int)):
                  logger.warning("Streams/Type Parameters should be float aor integer.")
                  sys.exit()
          else:
              logger.warning("Streams/Type Parameters should have 2 values list.") 
              sys.exit() 
   else:
      if cfg['Streams']['Type'] == 'KURT':
          logger.warning("Streams/Type Parameters is empty. Use as default [0.5]")
          cfg['Streams']['Type Parameters'] = [0.5]

   if not cfg['Streams']['Rotate']=="":
       if not isinstance(cfg['Streams']['Rotate'], bool):
           logger.warning("Streams/Rotate should be True/False")
           logger.warning("Use as Default False")
           cfg['Streams']['Rotate'] = False
       if str(cfg['Streams']['Rotate']) not in ['True', 'False']:
           logger.warning("Streams/Rotate should be True/False")
           logger.warning("Use as Default False")
           cfg['Streams']['Rotate'] = False
   else:
      logger.warning("Streams/Rotate is empty. Use as Default False")
      cfg['Streams']['Rotate'] = False

   if not cfg['Streams']['Response']=="":
       if not isinstance(cfg['Streams']['Response'], bool):
           logger.warning("Streams/Response should be True/False")
           logger.warning("Use as Default False")
           cfg['Streams']['Response'] = False
       if str(cfg['Streams']['Response']) not in ['True', 'False']:
           logger.warning("Streams/Response should be True/False")
           logger.warning("Use as Default False")
           cfg['Streams']['Response'] = False
   else:
      logger.warning("Streams/Response is empty. Use as Default False")
      cfg['Streams']['Response'] = False
 
   if bool(cfg['Streams']['Quantity']):
      if not isinstance(cfg['Streams']['Quantity'], str):
          logger.warning("Streams/Quantity must be 'VEL', 'ACC', 'DISP'")
          logger.warning("Use as Default VEL")
          cfg['Streams']['Quantity'] = 'VEL'
      if cfg['Streams']['Quantity'] not in ['VEL', 'ACC', 'DISP']:
          logger.warning("Streams/Quantity must be 'VEL', 'ACC', 'DISP'")
          logger.warning("Use as Default VEL")
          cfg['Streams']['Quantity'] = 'VEL'
   else:
      logger.warning("Streams/Quantity is empty. Use as Default VEL")
      cfg['Streams']['Quantity'] = 'VEL'

   if bool(cfg['Streams']['Normalize']):
       if len(cfg['Streams']['Normalize'])!=3:
           logger.warning("Streams/Normalize should have 3 position list. Use as default [True, 1, 1]")
           cfg['Streams']['Normalize'] = [True, 1, 1]     
       if str(cfg['Streams']['Normalize'][0]) not in ['True', 'False']:
           logger.warning("Streams/Normalize 1st position should be True/False. Use as default True")
           cfg['Streams']['Normalize'][0] = True
       if not isinstance(cfg['Streams']['Normalize'][1], (float, int)) or not isinstance(cfg['Streams']['Normalize'][2], (float, int)):
           logger.warning("Streams/Normalize 2nd/3rd position should be integers. Use as default  1, 1")
           cfg['Streams']['Normalize'][1] = 1
           cfg['Streams']['Normalize'][2] = 1       
   else:
       logger.warning("Streams/Normalize is empty. Use as default [True, 1, 1]")
       cfg['Streams']['Normalize'] = [True, 1, 1]

   if bool(cfg['Streams']['Corrections']):
      if len(cfg['Streams']['Corrections'])!=3:
          logger.warning("Streams/Corrections should have 3 position list. Use as default [False, 0, '']")
          cfg['Streams']['Corrections'] = [False, 0, '']
      if str(cfg['Streams']['Corrections'][0])  not in ['True', 'False']:
          logger.warning("Streams/Corrections 1st position should be True/False. Use as default False")    
          cfg['Streams']['Corrections'][0] = False
      if not isinstance(cfg['Streams']['Corrections'][1], (float, int)) or cfg['Streams']['Corrections'][1] not in [0, 1, 2, 3]:
          logger.warning("Streams/Corrections should have values between 0-3. Use as default 0")
   else:
       logger.warning("Streams/Corrections is empty. Use as default [False, 0, '']")
       cfg['Streams']['Corrections'] = [False, 0, '']

   if not isinstance(cfg['Streams']['Combine'], bool):
       logger.warning("Streams/Combine should be True/False")
       logger.warning("Use as Default False")
       cfg['Streams']['Combine'] = False
   if str(cfg['Streams']['Combine']) not in ['True', 'False']:
       logger.warning("Streams/Combine should be True/False")
       logger.warning("Use as Default False")
       cfg['Streams']['Combine'] = False

   if bool(cfg['Streams']['Filter']):
      if not isinstance(cfg['Streams']['Filter'], list):
          logger.warning("Streams/Filter should be list. Use [2,8] as default")
          cfg['Streams']['Filter'] = [[2,8]]
      if not len(cfg['Streams']['Filter'])>0:
          logger.warning("Streams/Filter should be list. Use [2,8] as default")
          cfg['Streams']['Filter'] = [[2,8]] 
      for i in cfg['Streams']['Filter']:
          if len(i)!=2:
              logger.warning("Streams/Filter should have lists with two variables. Use [2,8] as default")
              i = [2,8]
          if not isinstance(i[0], (float, int)) or not isinstance(i[1], (float, int)):
              logger.warning("Streams/Filter should have lists with numbers. Use [2,8] as default")
              i = [2,8]
          if isinstance(i[0], (float, int)) and isinstance(i[1], (float, int)):
              if i[0]>=i[1] and i[0]!=0 and i[1]!=0:
                  logger.warning("Streams/Filter should have list with the first number smaller that the second. Use [2,8] as default")
                  i = [2,8]
   else:
      logger.warning("Streams/Filter is emtpy.")
      logger.info('Program is terminated...')
      sys.exit()

   if bool(cfg['Streams']['Quality Control']):
       if isinstance(cfg['Streams']['Quality Control'], list):
                for i in cfg['Streams']['Quality Control'][0]:
                   if i not in ['SNR', 'CLIP', 'TIME']:
                      logger.error("Streams/Quality Control Problem")
                      sys.exit() 
       else:
           cfg['Streams']['Quality Control'] = []
   else:
      logger.warning("Streams/Quality Control is emtpy.")

   if bool(cfg['Streams']['Pre Filter']):
       if not isinstance(cfg['Streams']['Pre Filter'], list):
           logger.warning("Streams/Pre Filter should be list. Use [0.005, 45] as default")
           cfg['Streams']['Pre Filter'] = [0.005, 45]
       if not len(cfg['Streams']['Pre Filter'])==2:
           logger.warning("Streams/Pre Filter should be list. Use [0.005, 45] as default")
           cfg['Streams']['Pre Filter'] = [0.005, 45]

   if 'GPU' not in cfg['Backprojection']:
       logger.error("'Backprojection/GPU' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'NumCPUS' not in cfg['Backprojection']:
       logger.error("'Backprojection/NumCPUS' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
 
   if 'Grid' not in cfg['Backprojection']:
       logger.error("'Backprojection/Grid' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Selection' not in cfg['Backprojection']:
       logger.error("'Backprojection/Selection' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Sectors' not in cfg['Backprojection']:
       logger.error("'Backprojection/Sectors' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
 
   if 'Settings' not in cfg['Backprojection']:
       logger.error("'Backprojection/Settings' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Components' not in cfg['Backprojection']['Selection']:
       logger.error("'Backprojection/Selection/Components' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Distance' not in cfg['Backprojection']['Selection']:
       logger.error("'Backprojection/Selection/Distance' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Phase' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/Phase' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'TTmax' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/TTmax' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'ScanningTime' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/ScanningTime' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'TimeShift' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/TimeShift' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'MovingWindow' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/MovingWindow' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'BrType' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/BrType' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()  

   if 'StaThre' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/StaThre' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()  

   if 'bthre' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/bthre' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Weight' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/Weight' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Npower' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/Npower' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Normalize Results' not in cfg['Backprojection']['Settings']:
       logger.error("'Backprojection/Settings/Normalize Results' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()
   
   if cfg['Backprojection']['GPU'] is not None:
       if not isinstance(cfg['Backprojection']['GPU'], (bool)):
           logger.warning('Backprojection/GPU should be True/False. Use as default False.')
           cfg['Backprojection']['GPU'] = False
   else:
       logger.warning('Backprojection/GPU is empty. Use as default False.')
       cfg['Backprojection']['GPU'] = False

   if cfg['Backprojection']['NumCPUS'] is not None:
      if not isinstance(cfg['Backprojection']['NumCPUS'], (int)):
          logger.warning('Backprojection/NumCPUS should be integer or null. Use as default null.')
          cfg['Backprojection']['NumCPUS'] = None
   else:
       pass

   if cfg['Backprojection']['Grid'] is not None:
       if not isinstance(cfg['Backprojection']['Grid'], (list)):
            logger.error("'Backprojection/Settings/Grid should be list.")
            logger.info('Program is terminated...')
            sys.exit()
       else:
           for i in cfg['Backprojection']['Grid']:
               if not isinstance(i, (list)):
                   logger.error("'Backprojection/Settings/Grid should be list.")
                   logger.info('Program is terminated...')
                   sys.exit()
               if len(i)!=3:
                   logger.error("'Backprojection/Settings/Grid error.")
                   logger.info('Program is terminated...')
                   sys.exit()
   else:
       logger.error("'Backprojection/Settings/Grid is empty.")
       logger.info('Program is terminated...')
       sys.exit()

   if cfg['Backprojection']['Selection']['Components'] is not None:
       if not isinstance(cfg['Backprojection']['Selection']['Components'], (list)):
           logger.warning("Backprojection/Selection/Components should be list.  Use as default ['N', 'E']")
           cfg['Backprojection']['Selection']['Components'] = ['N', 'E']
       else:
           for i in cfg['Backprojection']['Selection']['Components']:
               if i not in ['N', 'E', 'Z', 'H', 'R', 'T']:
                   logger.warning("Backprojection/Selection/Components should be 'N', 'E', 'Z', 'H', 'R', 'T'. Use as default ['N', 'E'].")
                   cfg['Backprojection']['Selection']['Components'] = ['N', 'E']
                   break
   else:
       logger.warning("Backprojection/Selection/Components is empty. Use as default ['N', 'E']")
       cfg['Backprojection']['Selection']['Components'] = ['N', 'E']

   if cfg['Backprojection']['Selection']['Distance'] is not None:
       if not isinstance(cfg['Backprojection']['Selection']['Distance'], (float, int)):
           logger.warning("Backprojection/Selection/Components should be float or integer. Use as default 100.")
           cfg['Backprojection']['Selection']['Distance'] = 100
   else:
       logger.warning("Backprojection/Selection/Components is empty. Use as default 100.")
       cfg['Backprojection']['Selection']['Distance'] = 100

   if cfg['Backprojection']['Sectors'] is not None:
       if not isinstance(cfg['Backprojection']['Sectors'], (list)):
           logger.warning("Backprojection/Sectors sould be list. Use as default [2, 2].")
           cfg['Backprojection']['Sectors'] = [2,2]
       else:
           if len(cfg['Backprojection']['Sectors'])!=2:
               logger.warning("Backprojection/Sectors be 2 positions list. Use as default [2, 2].")
               cfg['Backprojection']['Sectors'] = [2,2]
           else:
               if cfg['Backprojection']['Sectors'][0]>4:
                   logger.warning("Backprojection/Sectors should have 1-4 sectors selected. Use as default [2, 2].")
                   cfg['Backprojection']['Sectors'] = [2,2]
   else:
       logger.warning("Backprojection/Sectors is emtpy. Use as default [2, 2].")
       cfg['Backprojection']['Sectors'] = [2,2]

   if not isinstance(cfg['Backprojection']['Settings']['Phase'], (str)):
       logger.warning("Backprojection/Settings/Phase sould be string. Use as default ['S'].")
       cfg['Backprojection']['Settings']['Phase'] = 'S'
   if cfg['Backprojection']['Settings']['Phase'] not in ['P', 'S']:
       logger.warning("Backprojection/Settings/Phase should be 'P' or 'S'. Use as default 'S'.")
       cfg['Backprojection']['Settings']['Phase'] = 'S'

   if cfg['Backprojection']['Settings']['TTmax'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['TTmax'], (float, int)):
           logger.warning("Backprojection/Settings/TTmax should be float or integer. Use as default 100.")
           cfg['Backprojection']['Settings']['TTmax'] = 100
       else:
           if cfg['Backprojection']['Settings']['TTmax'] < 1: 
               logger.warning("Backprojection/Settings/TTmax should be positive integer. Use as default 100.")
               cfg['Backprojection']['Settings']['TTmax'] = 100
   else:
       logger.warning("Backprojection/Settings/TTmax is empty. Use as default 100.")
       cfg['Backprojection']['Settings']['TTmax'] = 100

   if cfg['Backprojection']['Settings']['ScanningTime'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['ScanningTime'], (list)):
           logger.error("'Backprojection/Settings/ScanningTime should be list.")
           logger.info('Program is terminated...')
           sys.exit()
       else:
           for i in cfg['Backprojection']['Settings']['ScanningTime']:
               if len(i)!=3:
                   logger.error("Backprojection/Settings/ScanningTime should have 3 list positions.")
                   logger.info('Program is terminated...')
                   sys.exit()
               else:
                  if i[2][0]>=i[2][1]:
                       logger.error("Backprojection/Settings/ScanningTime problem.")
                       logger.info('Program is terminated...')
                       sys.exit()
                  if not isinstance(i[0], (float, int)) or not isinstance(i[1], (float, int)) or not isinstance(i[2][0], (float, int)) or not isinstance(i[2][1], (float, int)):
                       logger.error("Backprojection/Settings/ScanningTime problem.")
                       logger.info('Program is terminated...')
                       sys.exit()
   else:
       logger.error("'Backprojection/Settings/ScanningTime' in the config file is emtpy.")
       logger.info('Program is terminated...')
       sys.exit()


   if cfg['Backprojection']['Settings']['TimeShift'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['TimeShift'], (float, int)):
           logger.warning("Backprojection/Settings/TimeShift should be float or integer. Use as default 1.")
           cfg['Backprojection']['Settings']['TimeShift'] = 1
       else:
           if cfg['Backprojection']['Settings']['TimeShift']<=0:
               logger.warning("Backprojection/Settings/TimeShift should be positive integer. Use as default 1.")
               cfg['Backprojection']['Settings']['TimeShift'] = 1
   else:
      logger.warning("Backprojection/Settings/TimeShift is empty. Use as default 1.")
      cfg['Backprojection']['Settings']['TimeShift'] = 1


   if cfg['Backprojection']['Settings']['MovingWindow'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['MovingWindow'], (list)):
           logger.warning("Backprojection/Settings/MovingWindow should be list. Use as default [0.15, 0.15]")
           cfg['Backprojection']['Settings']['MovingWindow'] = [0.15, 0.15]
       else:
           if not len(cfg['Backprojection']['Settings']['MovingWindow'])==2:
               logger.warning("Backprojection/Settings/MovingWindow should be 2 position list. Use as default [0.15, 0.15]")
               cfg['Backprojection']['Settings']['MovingWindow'] = [0.15, 0.15]
           else:
               if not isinstance(cfg['Backprojection']['Settings']['MovingWindow'][0], ((float, int))) or not isinstance(cfg['Backprojection']['Settings']['MovingWindow'][1], ((float, int))): 
                   logger.warning("Backprojection/Settings/MovingWindow should have 2 floats or integers. Use as default [0.15, 0.15]")
                   cfg['Backprojection']['Settings']['MovingWindow'] = [0.15, 0.15]
               if cfg['Backprojection']['Settings']['MovingWindow'][0]<=0 or cfg['Backprojection']['Settings']['MovingWindow'][1]<=0:
                   logger.warning("Backprojection/Settings/MovingWindow should be positive numbers. Use as default [0.15, 0.15]")
                   cfg['Backprojection']['Settings']['MovingWindow'] = [0.15, 0.15]
   else:
       logger.warning("Backprojection/Settings/MovingWindow is empty. Use as default [0.15, 0.15]")
       cfg['Backprojection']['Settings']['MovingWindow'] = [0.15, 0.15]
  
   if cfg['Backprojection']['Settings']['BrType'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['BrType'], (int)):
           logger.warning("Backprojection/Settings/BrType should be 0 or 1. Use as default 1.")
           cfg['Backprojection']['Settings']['BrType'] = 1
       else:
           if cfg['Backprojection']['Settings']['BrType'] not in [0, 1]:
              logger.warning("Backprojection/Settings/BrType should be 0 or 1. Use as default 1.")
              cfg['Backprojection']['Settings']['BrType'] = 1
   else:
       logger.warning("Backprojection/Settings/BrType is empty. Use as default 1")
       cfg['Backprojection']['Settings']['BrType'] = 1


   if cfg['Backprojection']['Settings']['StaThre'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['StaThre'], (float, int)):
           logger.warning("Backprojection/Settings/StaThre should be between 0-1. Use as default 0.9")
           cfg['Backprojection']['Settings']['StaThre'] = 0.9
       else:
           if cfg['Backprojection']['Settings']['StaThre']<0 and cfg['Backprojection']['Settings']['StaThre']>1.0:
                logger.warning("Backprojection/Settings/StaThre should be between 0-1. Use as default 0.9")
                cfg['Backprojection']['Settings']['StaThre'] = 0.9
   else:
       logger.warning("Backprojection/Settings/StaThre is empty. Use as default 0.9")
       cfg['Backprojection']['Settings']['StaThre'] = 0.9

   if cfg['Backprojection']['Settings']['bthre'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['bthre'], (float, int)):
           logger.warning("Backprojection/Settings/bthre should be float or integer, Use as default 0.0")
           cfg['Backprojection']['Settings']['bthre'] = 0.0
       else:
           if cfg['Backprojection']['Settings']['bthre']<0:
               logger.warning("Backprojection/Settings/bthre should be positive number. Use as default 0.0")
               cfg['Backprojection']['Settings']['bthre'] = 0.0
   else:
       logger.warning("Backprojection/Settings/bthre is emtpy. Use as default 0.0")
       cfg['Backprojection']['Settings']['bthre'] = 0.0

   if cfg['Backprojection']['Settings']['Npower'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['Npower'], (int)):
           logger.warning("Backprojection/Settings/Npower Results should be integer. Use as default 1")
           cfg['Backprojection']['Settings']['Npower'] = 1
   else:
       logger.warning("Backprojection/Settings/Npower is empty. Use as default 1")
       cfg['Backprojection']['Settings']['Npower'] = 1

   if cfg['Backprojection']['Settings']['Normalize Results'] is not None:
       if not isinstance(cfg['Backprojection']['Settings']['Normalize Results'], (bool)):
           logger.warning("Backprojection/Settings/Normalize Results should be True/False. Use as default False")
           cfg['Backprojection']['Settings']['Normalize Results'] = False
   else:
       logger.warning("Backprojection/Settings/Normalize Results is empty. Use as default False")
       cfg['Backprojection']['Settings']['Normalize Results'] = False

   if 'Save Layers' not in cfg['Plotting']:
       logger.error("'Plotting/Save Layers' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Topography/Bathymetry' not in cfg['Plotting']:
       logger.error("'Plotting/Topography/Bathymetry' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Plots' not in cfg['Plotting']:
       logger.error("'Plotting/Plots' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Animation' not in cfg['Plotting']:
       logger.error("'Plotting/Animation' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if cfg['Plotting']['Save Layers'] is not None:
       if not isinstance(cfg['Plotting']['Save Layers'], str):
           logger.error("'Plotting/Save Layers' shoul be string.")
           logger.info('Program is terminated...')
           sys.exit()
   else:
       logger.error("'Plotting/Save Layers' is empty.")
       logger.info('Program is terminated...')
       sys.exit()

   if cfg['Plotting']['Topography/Bathymetry'] is not None:
       if not isinstance(cfg['Plotting']['Topography/Bathymetry'], (list)):
           logger.warning("Plotting/Topography/Bathymetry should be a list. Use as default False")
           cfg['Plotting']['Topography/Bathymetry'] = [False, '']
       else:
           if not len(cfg['Plotting']['Topography/Bathymetry'])==2:
               logger.warning("Plotting/Topography/Bathymetry should be two psition list. Use as default False")
               cfg['Plotting']['Topography/Bathymetry'] = [False, '']
           else:
               if not isinstance(cfg['Plotting']['Topography/Bathymetry'][0], (bool)):
                    cfg['Plotting']['Topography/Bathymetry'] = [False, '']
   else:
       logger.warning("Plotting/Topography/Bathymetry is empty. Use as default False")
       cfg['Plotting']['Topography/Bathymetry'] = [False, '']

   if cfg['Plotting']['Plots'] is not None:
       if not isinstance(cfg['Plotting']['Plots'], (bool)):
           logger.warning("'Plotting/Plots' should be True/False. Use as default False")
           cfg['Plotting']['Plots'] = False
   else:
       logger.warning("'Plotting/Plots' is empty. Use as default False")
       cfg['Plotting']['Plots'] = False
     

   if cfg['Plotting']['Animation'] is not None:
       if not isinstance(cfg['Plotting']['Animation'], (bool)):
           logger.warning("'Plotting/Animation' should be True/False. Use as default False")
           cfg['Plotting']['Animation'] = False
   else:
       logger.warning("'Plotting/Animation' is empty. Use as default False")
       cfg['Plotting']['Animation'] = False

   if 'Array Response Function' not in cfg['Tests']:
       logger.error("'Tests/Array Response Function' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Jackknife' not in cfg['Tests']:
       logger.error("'Tests/Jackknife' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Bootstrap' not in cfg['Tests']:
       logger.error("'Tests/Bootstrap' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if 'Delete' not in cfg:
       logger.error("'Tests/Delete' in the config file is missing.")
       logger.info('Program is terminated...')
       sys.exit()

   if cfg['Tests']['Array Response Function'] is not None:
       if not isinstance(cfg['Tests']['Array Response Function'], (bool)):
           logger.warning("'Tests/Array Response Function' should be True/False. Use as default False.")
           cfg['Tests']['Array Response Function'] = False
   else:
       logger.warning("'Tests/Array Response Function' is empty. Use as default False.")
       cfg['Tests']['Array Response Function'] = False

   if cfg['Tests']['Jackknife'] is not None:
       if not isinstance(cfg['Tests']['Jackknife'], (bool)):
           logger.warning("'Tests/Jackknife' should be True/False. Use as default False.")
           cfg['Tests']['Jackknife'] = False
   else:
       logger.warning("'Tests/Jackknife' is empty. Use as default False.")
       cfg['Tests']['Jackknife'] = False


   if cfg['Delete'] is not None:
       if not isinstance(cfg['Delete'], (bool)):
           logger.warning("'Delete' should be True/False. Use as default True.")
           cfg['Tests']['Delete'] = True
   else:
       logger.warning("'Tests/Delete' is empty. Use as default True.")
       cfg['Tests']['Delete'] = False

   if cfg['Tests']['Bootstrap'] is not None:
       if not isinstance(cfg['Tests']['Bootstrap'], (list)):
           logger.warning("'Tests/Bootstrap should be list. Use as default [False, 5, 50]'")
           cfg['Tests']['Bootstrap'] = [False, 5, 50]
       else:
           if len(cfg['Tests']['Bootstrap'])!=3:
               logger.warning("'Tests/Bootstrap should be list with 3 positions. Use as default [False, 5, 50]'")
               cfg['Tests']['Bootstrap'] = [False, 5, 50]
           else:
               if not isinstance(cfg['Tests']['Bootstrap'][0], (bool)):
                   cfg['Tests']['Bootstrap'] = [False, 5, 50]
               if not isinstance(cfg['Tests']['Bootstrap'][1], (int, float)) or not isinstance(cfg['Tests']['Bootstrap'][2], (int, float)):
                   cfg['Tests']['Bootstrap'] = [False, 5, 50]
   else:
       cfg['Tests']['Bootstrap'] = [False, 5, 50]
 

   return cfg









   




















  









 











   return
