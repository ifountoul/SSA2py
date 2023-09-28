#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#############
import datetime, sys, argparse, multiprocessing, textwrap, os

# Obspy Imports
###############
from obspy import Catalog, read_events, read
from obspy.core.inventory.inventory import read_inventory
from obspy.core import UTCDateTime


#Local Imports
##############
from SSA2py.core.traveltimes import calculateTraveltimes
from SSA2py.core import event, download, config, stream, BP_statistics, check_config
from SSA2py.core import backprojection, ARF
from SSA2py.core.basic_f.other import _compress_, delete_npy 
from SSA2py.core.basic_f.get_tt import getTTtables

# plotting
from SSA2py.core.plotting_functions.plot_res import plot_res_, plot_Un
from SSA2py.core.plotting_functions.FigureDown import FigDown

# End of Imports
#################



#######################
# Start the Main Code #
#######################


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

#######################
# Arguments in parser #
#######################

parser.add_argument('-c', '--config', metavar=('FILEPATH'), help='Override default configuration file (./config.yaml).', nargs='?', default='./config.yaml', type=str)

#####
help_event='EVENT can be in any of the following formats: \n(i) DATETIME MAGNITUDE TYPE LATITUDE LONGITUDE DEPTH\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat() + ' 3.0 ML 37.24 20.49 4.1\n(ii) DATETIME\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat()+'\n(iii) EVENTID (event identifiers are data center specific)\ne.g.: ' + __file__ + ' -e noa2020owyrp\nIn cases (ii) and (iii) the rest of the information is retrieved by the FDSNWS-event \nIn more than one results, only the first event is returned \nPassing milliseconds is optional.'
#####


parser.add_argument('-e', '--event', metavar=('EVENT'), help=help_event, action='append', type=str, nargs='+', default=None)
parser.add_argument('--event-file', metavar=('FILEPATH'), help='Parse and run a file with EVENT lines.', type=argparse.FileType('r'), default=None)
parser.add_argument('--event-xml', metavar=('FILEPATH'), help='Parse and run a file in QuakeML.', type=argparse.FileType('r'), default=None)
parser.add_argument("--real-time", action="store_true", help='Invoke --datetime-range for real-time use (FDSN bounded).', default=False)
parser.add_argument('-d', '--datetime-range', metavar=('TIME'), help='Invoke SSA computation for all events found in specific datetime range (FDSN bounded).', action='append', nargs=2, default=None)
parser.add_argument("-s", "--station", metavar=('STA[.NEZ]'), help='Override default stations selection. Optionally,\
                    components could be also specified.\nIt can be combined with --remove for the reverse result.', action="extend", nargs="+", type=str, default=[])
parser.add_argument("--repeat", help='Run again the SSA. Suggested to use -s to choose or exclude stations.', action="store_true", default=False)
parser.add_argument("--remove", help='It can only be used with --station. Invokes reverse result.', action="store_true", default=False)
parser.add_argument("--download", help='Download Faults and Tectonic Plates.', action="store_true", default=False)
parser.add_argument('-l', '--log', metavar=('FILEPATH'), help='Main log file.', nargs='?', default='./log', type=str)
parser.add_argument('-v', '--version', action='version', version='%(prog)s v'+__version__)

###########################
# End of Parser Arguments #
###########################

# decorate (overload) print help message
parser.print_help=print_help(parser.print_help)

args=parser.parse_args(args=None if sys.argv[1:] else ['--help'])


# Initiate SSA logger
if args.log is None:
    args.log = './log'

logger = config.setup_logger('SSA_log', args.log)
        

# Suppressing deprecation warnings
config.warns()


#############################
# Read the configuration file
_cfg=config.read(args.config)

# Check for input errors
_cfg=check_config.check(_cfg, logger)
#############################



# Number of cores
if _cfg['Backprojection']['NumCPUS']==None:
    _cfg['Backprojection']['NumCPUS']=multiprocessing.cpu_count()



# Download from github layers for the plots
if args.download:
   #Check if we have the figure props
    logger.info('Download plates boundaries, Faults etc. ---> ' + _cfg['Plotting']['Save Layers'])    
    FigDown(_cfg['Plotting']['Save Layers']) 
    sys.exit()

##########################
# Get events

cat=Catalog()

logger.info("Get the seismic events.")

# Get events from text file
if args.event_file:
    try:
        if os.path.exists(args.event_file):
            with args.event_file as _f:
                args.event=[_e.rstrip().split(' ') for _e in _f.readlines() if not _e[0]=='#']
        else:
            logger.info("The path '{}' does not exist.".format(args.event_file))
            logger.info("Exiting the Program.")
            sys.exit()
    except Exception as e:
         logger.error("An error occurred: {}".format(e))
         logger.info("Exiting the Program.")
         sys.exit()

# Get events from XML
if args.event_xml:
    try:
        if os.path.exists(args.event_xml):
            cat=read_events(args.event_xml)
        else:
            logger.info("The path '{}' does not exist.".format(args.event_xml))
            logger.info("Exiting the Program.")
            sys.exit()
    except Exception as e:
         logger.error("An error occurred: {}".format(e))
         logger.info("Exiting the Program.")
         sys.exit()

# Get events for datetime range
if args.datetime_range:

    try:
        start=args.datetime_range[0][0]
        end=args.datetime_range[0][1]

        minMag=min([_[0] for _ in _cfg['Backprojection']['Settings']['ScanningTime']])
        temp_event = event.getFDSNWSCatalog(_cfg, logger, minMag=minMag, starttime=UTCDateTime(start), endtime=UTCDateTime(end))

        if temp_event!=None:
            cat+=temp_event
    except Exception as e:
         logger.error("An error occurred: {}".format(e))
         logger.info("Exiting the Program.")
         sys.exit()


# Real time
if args.real_time:
    cat+=event.monitor(_cfg, logger)

# Seismic events
if args.event:
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
            logger.warning(" Problem with event {} .".format(text))
            continue

# Do we have seismic events in catalog ?
if len(cat)==0:
    logger.info("No seismic events in the catalog. Exiting the Program.")
    sys.exit()
else:
    logger.info('{} seismic event/s in the catalog.'.format(len(cat)))

##########################

# Get data for each event in catalog

for evt in cat:

    try:

        config.init(evt, _cfg)
        config.cfg = config.read(args.config)

        # Do we use ESM servise?
        for service in config.cfg['Download Service']['Stream']:
            if service[1] == 'https://esm-db.eu':
                config.eventid_ESM = str(evt.resource_id).split('event_id=')[1]

        # Check that the event was not executed before
        if config.exists==True and not args.repeat:
             # Events exists
             config.logger.info('Event has been already executed move to the next event in the catalog...')
             continue

        try:
             # Get waveforms and inventory if not repeat
            if not args.repeat:
                config.logger.info('Download Data and Inventory for {}\n'.format(evt))
                config.logger.info('Save to: {}'.format(config.eventdir))
                # get inventory
                download.getInventory()       

                # Empty Metadata raise exception!
                if len(config.inv) == 0:
                     raise Exception('Empty Metadata Variable!')
 
                # get waveforms
                download.getWaveforms()

                # Empty Waveforms raise exception!
                if len(config.st) == 0:
                    raise Exception('Empty Streams Variable!')

                # remove from INVENTORY stations that don't exists in the initial downloaded STREAM
                download.getClear()

        except Exception as e:
                config.logger.error(e)
                config.logger.info('Error trying to get data - metadata. Continue to next event.')
                continue

        # Traveltimes Calculation
        try:
            # If not repeat
            if not args.repeat:
                config.logger.info('Traveltime tables calculation')
                calculateTraveltimes(_cfg)
        except Exception as e:
                config.logger.error(e)
                config.logger.info('Error trying to calculate the traveltimes. Continue to next event.')
                continue
    
        # If we want to repeat the procedure

        if args.repeat:
            config.logger.info('Revise the SSA solution')
 
            # read the STREAMS and the INVENTORY
            files = list((file for file in os.listdir(config.eventdir) if os.path.isfile(os.path.join(config.eventdir, file))))

            try:
                    config.st+=read(os.path.join(config.eventdir, files[0]))
                    config.inv = read_inventory(os.path.join(config.inventorydir, 'inventory.xml'))
                    config.st.attach_response(config.inv)
            except Exception as e:
                    config.logger.error(e)
                    config.logger.info('No stream or inventory available. Continue to next event.')
                    continue

            # remove stations or components from the STREAM based on command line
            if len(args.station)>0:
                if args.remove: #Exclude stations from stream 
                    config.stations_status = 'r' #Remove
                    config.stations_ = args.station
                else:
                    config.stations_status = 'k' #Keep
                    config.stations_ = args.station
            stream.commandRemoveST()

        # QC
        stream.clean()

        # waveforms processing
        stream.stream_process()
 
        ###############################
        # Backprojection Starts Here! #
        ###############################
 
        # Get Grid based on config file

        config.grid = backprojection.getGrid()
   
        # Get tttables for the model/phase
        config.phase = config.cfg['Backprojection']['Settings']['Phase']

        # Get the TT tables

        getTTtables()

        # Loop per component

        # Component
        for comp in config.cfg['Backprojection']['Selection']['Components']:
            for fi in config.cfg['Streams']['Filter']: #Filter
                #MSEED path to use
                path = os.path.join(config.eventdir, 'Processed_Data',\
                                   '{:s}_{:s}_{:s}_{:s}{:s}'.format(config.cfg['Streams']['Type'],\
                                    str(float(fi[0])), str(float(fi[1])), comp, '.mseed'))
                config.comp = comp; config.fi = fi;

                # read the Stream
                config.st, config.stations = stream.StreamReady(path) 

                # SSA solutions folder
                config.job = 'Main SSA'

                SSA_output_dir = os.path.join(config.eventdir, 'Results', 'SSA',\
                                              os.path.basename(os.path.normpath(path)).split('.mseed')[0],\
                                             'Detailed_Solution')

                ba = backprojection.backprojection(config.st, config.stations, SSA_output_dir)

                if ba==True:

                    # Plot results
                    plot_save_dir = os.path.join(config.eventdir, 'Results', 'SSA',\
                                    os.path.basename(os.path.normpath(path)).split('.mseed')[0],\
                                    'Plots')

                    plot_res_([SSA_output_dir], plot_save_dir)    

                    if config.cfg['Tests']['Jackknife'] is True or config.cfg['Tests']['Bootstrap'][0] is True:
                         config.logger.info('Moving to Uncertainty Test...')

                         if config.cfg['Tests']['Jackknife'] is True:
                             Unpaths = []
                             Unpaths.append(os.path.join(config.eventdir, 'Results', 'JACKKNIFE',\
                                            os.path.basename(os.path.normpath(path)).split('.mseed')[0]))
                             Unpaths.append(SSA_output_dir)

                             # Perform the Jackknife test to determine the 95% Confidence Interval
                             BP_statistics.Jackknife_stats(config.st.copy(), config.stations, Unpaths)

                             # plot Jackknife results 
                             # Save plots dir
                             Jack_Paths_out = os.path.join(config.eventdir, 'Results', 'JACKKNIFE',\
                                              os.path.basename(os.path.normpath(path)).split('.mseed')[0], 'Plots')

                             plot_Un(Unpaths, Jack_Paths_out)

                             # Compress the results
                             if config.cfg['Delete']==True and os.path.exists(os.path.join(Unpaths[0], 'Main_Solution')):
                                 _compress_([os.path.join(Unpaths[0], 'Main_Solution')])
                                 for d in os.listdir(os.path.join(Unpaths[0], 'Detailed_Solutions')):
                                     delete_npy(os.path.join(Unpaths[0], 'Detailed_Solutions', d),[])

                         if config.cfg['Tests']['Bootstrap'][0] is True:  
                             Unpaths = []
                             Unpaths.append(os.path.join(config.eventdir, 'Results', 'BOOTSTRAP',\
                                            os.path.basename(os.path.normpath(path)).split('.mseed')[0]))
                             Unpaths.append(SSA_output_dir)
       

                             config.logger.info('Moving to Bootstrap Test...')
                       
                             BP_statistics.Bootstrap_stats(config.st.copy(), config.stations, Unpaths,\
                                                           config.cfg['Tests']['Bootstrap'][1], config.cfg['Tests']['Bootstrap'][2])  
 
                                      
                             Boot_Paths_out = os.path.join(config.eventdir, 'Results', 'BOOTSTRAP',\
                                              os.path.basename(os.path.normpath(path)).split('.mseed')[0], 'Plots')
                           
                             plot_Un(Unpaths, Boot_Paths_out)
            
                             # Compress the results
                             if config.cfg['Delete']==True and os.path.exists(os.path.join(Unpaths[0], 'Main_Solution')):
                                 _compress_([os.path.join(Unpaths[0], 'Main_Solution')])
                                 for d in os.listdir(os.path.join(Unpaths[0], 'Detailed_Solutions')):
                                     delete_npy(os.path.join(Unpaths[0], 'Detailed_Solutions', d),[])

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
             getTTtables()

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

                 if config.cfg['Delete']==True:
                     _compress_([os.path.join(config.eventdir, 'Results', 'ARF', 'pulses_'+comp, 'Detailed_Solution')])

             config.fi = old_fi
             config.scanningRules[0] = old_scan
             config.cfg['Backprojection']['Settings']['TimeShift'] = old_shift
             config.cfg['Backprojection']['Settings']['MovingWindow'] = old_window
             config.cfg['Streams']['Type'] = old_type
                

        # Compress the results
        if config.cfg['Delete']==True and os.path.exists(SSA_output_dir):
            _compress_([SSA_output_dir]) 

    except:
         config.logger.exception('Fatal Error occurred')
         config.logger.info('Moving to next event, if any')
    


