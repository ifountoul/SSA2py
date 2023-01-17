#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime, sys, argparse, yaml, logging, shutil, copy
from argparse import RawTextHelpFormatter
import textwrap, os, glob, json, sys, copy
from obspy import Catalog, Stream, read_events, read
import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.inventory.inventory import read_inventory
from obspy.core import UTCDateTime


#Local
from SSA2py.core.traveltimes import calculateTraveltimes
from SSA2py.core import event, download, config, stream, BP_statistics
from SSA2py.core import backprojection, ARF
from SSA2py.core.modules.mute import mute_traces
# plotting
from SSA2py.core.plotting_functions.plot_res import plot_res_, plot_Boot
from SSA2py.core.plotting_functions.FigureDown import FigDown

from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


__copyright__ = "🄯 "+str(datetime.datetime.utcnow().year)+", Institute of Geodynamics - National Observatory of Athens"
__credits__ = ["Ioannis Fountoulakis (ifountoul@noa.gr), Christos Evangelidis (cevan@noa.gr)"]
__author__= 'Ioannis Fountoulakis (ifountoul@noa.gr)'
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Ioannis Fountoulakis"
__email__ = "ifountoul@noa.gr"
__status__ = "Production"

os.environ['COLUMNS'] = "90"


def print_help(func):
    def inner():
        print("""
░██████╗░██████╗░█████╗░██████╗░██████╗░██╗░░░██╗
██╔════╝██╔════╝██╔══██╗╚════██╗██╔══██╗╚██╗░██╔╝
╚█████╗░╚█████╗░███████║░░███╔═╝██████╔╝░╚████╔╝░
░╚═══██╗░╚═══██╗██╔══██║██╔══╝░░██╔═══╝░░░╚██╔╝░░
██████╔╝██████╔╝██║░░██║███████╗██║░░░░░░░░██║░░░
╚═════╝░╚═════╝░╚═╝░░╚═╝╚══════╝╚═╝░░░░░░░░╚═╝░░░""")
        print('SSA2py: Source Scanning Algorithm in Python\n')
        print('Version: ' + __version__)
        print('License: ' + __license__)
        print('Author: ' + __author__)
        print ('Credits: '+ textwrap.TextWrapper().fill(text=''.join(__credits__)))
        print( __copyright__+'\n')
        func()
    return inner

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return textwrap.wrap(text, width)


class ExtendAction(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest) or []
        items.extend(values)
        setattr(namespace, self.dest, items)

def required_length(nmin,nmax):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest,nmin=nmin,nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength

parser = argparse.ArgumentParser(prog='SSA2py', description='SSA2py: Source Scanning Algorithm in Python\nFind more info at:',\
                                formatter_class=LineWrapRawTextHelpFormatter)

parser.register('action', 'extend', ExtendAction)

#Arguments in parser
parser.add_argument('-c', '--config', metavar=('FILEPATH'), help='override default configuration file (./config.yaml)', nargs='?', default='./config.yaml', type=str)
#####
help_event='EVENT can be in any of the following formats: \n(i) DATETIME MAGNITUDE TYPE LATITUDE LONGITUDE DEPTH\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat() + ' 3.0 ML 37.24 20.49 4.1\n(ii) DATETIME\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat()+'\n(iii) EVENTID (event identifiers are data center specific)\ne.g.: ' + __file__ + ' -e noa2020owyrp\nIn cases (ii) and (iii) the rest of the information is retrieved by the FDSNWS-event \nIn more than one results, only the first event is returned \nPassing milliseconds is optional.'

#####
parser.add_argument('-e', '--event', metavar=('EVENT'), help=help_event, action='append', type=str, nargs='+', default=None)
parser.add_argument('--event-file', metavar=('FILEPATH'), help='parse and run a file with EVENT lines', type=argparse.FileType('r'), default=None)
parser.add_argument('--event-xml', metavar=('FILEPATH'), help='parse and run a file in QuakeML', type=argparse.FileType('r'), default=None)
parser.add_argument("--real-time", action="store_true", help='invoke --datetime-range for real-time use (FDSN bounded)', default=False)
parser.add_argument('-d', '--datetime-range', metavar=('TIME'), help='invoke SSA computation for all events found in specific datetime range (FDSN bounded)', action='append', nargs=2, default=None)


parser.add_argument("-s", "--station", metavar=('STA[.NEZ]'), help='override default stations selection. Optionally,\
                    components could be also specified.\nIt can be combined with --remove for the reverse result', action="extend", nargs="+", type=str, default=[])

parser.add_argument("--repeat", help='run again the SSA. Suggested to use -s to choose or exclude stations', action="store_true", default=False)
parser.add_argument("--remove", help='it can only be used with --station. Invokes reverse result', action="store_true", default=False)

parser.add_argument("--disable-quality", help='disable quality control in Stream processing', action="store_false", default=False)

parser.add_argument("--download", help='download Faults and Tectonic Plates', action="store_true", default=False)

parser.add_argument('-l', '--log', metavar=('FILEPATH'), help='override default main log file (./log)', nargs='?', default='./log', type=str)
parser.add_argument('-v', '--version', action='version', version='%(prog)s v'+__version__)

# decorate (overload) print help message
parser.print_help=print_help(parser.print_help)

args=parser.parse_args(args=None if sys.argv[1:] else ['--help'])

# initiate SSA logger
logger = config.setup_logger('SSA_log', args.log)

# suppressing deprecation warnings
config.warns()

# read the configuration file
_cfg=config.read(args.config)


# download from github layers for the plots
if args.download:
   #Check if we have the figure props
    logger.info('Download plates boundaries, Faults etc. ---> ' + _cfg['Plotting']['Save Layers'])    
    FigDown(_cfg['Plotting']['Save Layers']) 
    sys.exit()

# disable quality control 
if args.disable_quality:
    _cfg['Streams']['Quality Control'] = args.disable-quality

cat=Catalog()

# get events from text file
if args.event_file:
    with args.event_file as _f:
        args.event=[_e.rstrip().split(' ') for _e in _f.readlines() if not _e[0]=='#']

# get events from XML
if args.event_xml:
    cat=read_events(args.event_xml)

# get events for datetime range
if args.datetime_range:
    start=args.datetime_range[0][0]
    end=args.datetime_range[0][1]

    minMag=min([_[0] for _ in _cfg['Backprojection']['Settings']['ScanningTime']])
    cat+=event.getFDSNWSCatalog(_cfg, logger, minMag=minMag, starttime=UTCDateTime(start), endtime=UTCDateTime(end))

# real time
if args.real_time:
    cat+=event.monitor(_cfg, logger)

if args.event: # get info of the event from FDSN of Local input
    for text in args.event:
        try:
            # if event is FDSN bound
            if len(text)==1:
                try:
                    cat.append(event.getFDSNWSCatalog(_cfg, logger, eventid=text[0])[0])
                except:
                    cat.append(event.getFDSNWSCatalog(_cfg, logger, starttime=UTCDateTime(text)-1, endtime=UTCDateTime(text)+1)[0])
            # create a manual object
            elif len(text)==6:
                cat.append(event.getCatalog(text)[0])
        except:
            pass

# get data for each event in catalog
for evt in cat:
    config.init(evt, _cfg)
    config.cfg = config.read(args.config)

    try:
        # get waveforms and inventory if not repeat
        if not args.repeat:
            config.logger.info('Download Data and Inventory for {}\n'.format(evt))
            config.logger.info('Save to: {}'.format(config.eventdir))
            # get inventory
            download.getInventory()       
            # get waveforms
            download.getWaveforms()
            # remove from INVENTORY stations that don't exists in the initial downloaded STREAM
            download.getClear()

            # traveltime tables calculation
            config.logger.info('Traveltime tables calculation')
            calculateTraveltimes(_cfg)

            # quality Control
            stream.clean()

            # waveforms processing
            stream.stream_process()
 
        else:
            config.logger.info('Revise the SSA solution')
            # read the STREAMS and the INVENTORY
            files = list((file for file in os.listdir(config.eventdir) if os.path.isfile(os.path.join(config.eventdir, file))))
            try:
                config.st+=read(os.path.join(config.eventdir, files[0]))
                config.inv = read_inventory(os.path.join(config.inventorydir, 'inventory.xml'))
                config.st.attach_response(config.inv)
            except Exception as e:
                config.logger.warning('No stream or inventory available. End of program.')
                pass

            # remove stations or components from the STREAM based on command line
            if len(args.station)>0:
                if args.remove: #Exclude stations from stream 
                    config.stations_status = 'r' #Remove
                    config.stations_ = args.station
                else:
                    config.stations_status = 'k' #Keep
                    config.stations_ = args.station
            stream.commandRemoveST()

            # quality Control
            stream.clean()

            # waveforms processing
            stream.stream_process()

        # backprojection  
        # get Grid based on config file
        config.grid = backprojection.getGrid()
   
        # get tttables for the model/phase
        config.phase = config.cfg['Backprojection']['Settings']['Phase'][0]
        backprojection.getTTtables()
        # loop per component
        for comp in config.cfg['Backprojection']['Selection']['Components']: #Component
            for fi in config.cfg['Streams']['Filter']: #Filter
               #MSEED path to use
                path = os.path.join(config.eventdir, 'Processed_Data',\
                                   '{:s}_{:s}_{:s}_{:s}{:s}'.format(config.cfg['Streams']['Type'],\
                                    str(float(fi[0])), str(float(fi[1])), comp, '.mseed'))
                config.comp = comp; config.fi = fi;

                # read the Stream
                config.st, config.stations = stream.StreamReady(path) 

                # mute traces
                if config.cfg['Backprojection']['Settings']['Mute'][0]==True:
                    config.st = mute_traces(config.st, config.stations)

                # SSA solutions folder
                config.job = 'Main SSA'
                SSA_output_dir = os.path.join(config.eventdir, 'Results', 'SSA',\
                                 os.path.basename(os.path.normpath(path)).split('.mseed')[0],\
                                'Detailed_Solution')
                ba = backprojection.backprojection(config.st, config.stations, SSA_output_dir)
                if ba==True:
                     # plot results
                     plot_save_dir = os.path.join(config.eventdir, 'Results', 'SSA',\
                                     os.path.basename(os.path.normpath(path)).split('.mseed')[0],\
                                    'Plots')

                     plot_res_([SSA_output_dir], plot_save_dir)    

                     if config.cfg['Tests']['Jackknife'] is True:
                         config.logger.info('Moving to Jackknife Test...')

                         # organize the paths
                         Jack_Paths = []
                         # main Jack 1)
                         Jack_Paths.append(os.path.join(config.eventdir, 'Results', 'JACKKNIFE',\
                                           os.path.basename(os.path.normpath(path)).split('.mseed')[0])) 
                         # SSA main solutions 4)
                         Jack_Paths.append(SSA_output_dir)

                         # perform the Jackknife test to determine the 95% Confidence Interval
                         BP_statistics.Jackknife_stats(config.st.copy(),\
                                                       config.stations, Jack_Paths, confidence_level=0.95)

                         # plot Jackknife results 
                         # Save plots dir
                         Jack_Paths_out = os.path.join(config.eventdir, 'Results', 'JACKKNIFE',\
                                          os.path.basename(os.path.normpath(path)).split('.mseed')[0], 'Plots')

                         plot_Boot(Jack_Paths, Jack_Paths_out, Test='JACK',\
                                   error_type=config.cfg['Tests']['Uncertainty Parameter'])


                     if config.cfg['Tests']['Bootstrap'][0] is True:  
                         config.logger.info('Moving to Bootstrap Test...')
                        
                         # organize the paths
                         Boot_Paths = []

                         # main Jack 1)
                         Boot_Paths.append(os.path.join(config.eventdir, 'Results', 'BOOTSTRAP',\
                                           os.path.basename(os.path.normpath(path)).split('.mseed')[0]))

                         Boot_Paths.append(SSA_output_dir)

                         # perform the Bootstrap test to determine the 95% Confidence Interval
                         BP_statistics.Bootstrap_stats(config.st.copy(),\
                                                       config.stations, Boot_Paths, config.cfg['Tests']['Bootstrap'][1],\
                                                       config.cfg['Tests']['Bootstrap'][2], confidence_level=0.95) 
                                      
                         Boot_Paths_out = os.path.join(config.eventdir, 'Results', 'BOOTSTRAP',\
                                          os.path.basename(os.path.normpath(path)).split('.mseed')[0], 'Plots')
                           
                         plot_Boot(Boot_Paths, Boot_Paths_out, Test='BOOT',\
                                   error_type=config.cfg['Tests']['Uncertainty Parameter'])
                             

        # for each phase do the ARF test
        if config.cfg['Tests']['Array Response Function'] is True:
            
            config.job = 'Array Response Function'

            old_fi = config.fi
            old_scan = config.scanningRules[0]
            old_shift = config.cfg['Backprojection']['Settings']['TimeShift']
            old_window = config.cfg['Backprojection']['Settings']['MovingWindow']
            old_type = config.cfg['Streams']['Type']

            config.fi = []
            config.scanningRules[0] = [-2,2]
            config.cfg['Backprojection']['Settings']['TimeShift'] = 0.1
            config.cfg['Backprojection']['Settings']['MovingWindow'] = [0, 0]

            config.cfg['Streams']['Type'] = 'Synthetic Pulses'

            config.phase = config.cfg['Backprojection']['Settings']['Phase'][0]
            backprojection.getTTtables()

            for comp in config.cfg['Backprojection']['Selection']['Components']:         

                #MSEED path to use
                path = os.path.join(config.eventdir, 'Processed_Data',\
                                  '{:s}_{:s}_{:s}_{:s}{:s}'.format(old_type,\
                                   str(float(config.cfg['Streams']['Filter'][0][0])),\
                                   str(float(config.cfg['Streams']['Filter'][0][1])), comp, '.mseed'))
                config.comp = comp
                # read the Stream
                config.st, config.stations = stream.StreamReady(path)

                # call the function
                ARF.ARF()

            config.fi = old_fi
            config.scanningRules[0] = old_scan
            config.cfg['Backprojection']['Settings']['TimeShift'] = old_shift
            config.cfg['Backprojection']['Settings']['MovingWindow'] = old_window
            config.cfg['Streams']['Type'] = old_type

    except:
         config.logger.exception('Fatal Error occurred')
         config.logger.info('Moving to next event, if any')
    


