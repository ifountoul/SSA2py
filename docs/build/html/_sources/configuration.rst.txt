
.. configuration:

=============
Configuration
=============

.. raw:: html

    <div class='Configuration'>

One important factor for the proper use of SSA2py is the wise calibration of the configuration file. A prototype is given together with the code (config.yaml).
SSA2py initially tries to reach the default configuration file (config.yaml) but the user can given alternatives with the use of the **-c** prefix and a path.
The user can use different configuration files based to his needs. In that way, SSA2py becomes adjustable to fit in any circumstances.
For example the user can define various configurations for sub-regions or having different workspaces where the SSA solutions are archived, etc. 

The parameters of the configuration file will be presented here in detail:

- ``Version: 1.0`` *Version of the code.*

- ``Events Dir: SSA2py/path/to dir`` *Directory where all the calculations will be performed and archived. For every seismic event a separate directory is created, named with the event's origin time.*

- ``Traveltimes:`` *Traveltimes section.*

  ``Crustals1D:`` *1D Velocity Models subsection.*

  ``Filename: SSA2py/path/to dir/test.vz`` *Path of the Velocity Model.*

  ``Elevation: int`` *Calculate traveltime tables reaching this elevation in km. e.g. 3*

  .. note:: int: integers, float: floating-point numbers, str: strings. 

  ``Depth: float`` *Maximum depth of the traveltime tables (km) e.g. 90.0*

  ``Distance: float`` *Maximum distance of the traveltime tables (km) e.g. 250.0*

  ``Granularity: float`` *Step for the traveltime calculations (km) e.g. 1.0*

  ``VPVS: float/null`` *You can set a Vp/Vs value. In that way, the SSA2py can use this value to calculate S traveltimes. If you set this value the S wave column in the given model will be ignored if exists.*

  ``Geobox: (Lon. 1, Lat. 1), (Lon. 2, Lat. 2), (Lon. 3, Lat. 3), (Lon. 4, Lat. 4)`` *The 'Geobox' is a polygon that defines the region where a crustal model is going to be used if the initial hypocenter's location is within it. If 'Geobox' is set to null, then this crustal model will be triggered if no other model found (used as default crustal model). Only one crustal model with null 'Geobox' can be defined.*

  ``Crustals3D:`` *3D Velocity Models subsection.*

  ``Filename VP: SSA2py/path/to dir/test.vp`` *Path of P Velocity Model.*

  ``Filename VS: SSA2py/path/to dir/test.vs`` *Path of S Velocity Model.*

  ``Depth: float`` *Maximum depth of the traveltime tables (km).*

  ``Distance: float`` *Maximum distance of the traveltime tables (km).*

  ``Granularity: float`` *Step for the traveltime calculations (km).*

  ``VPVS: float/null`` *Vp/Vs value*

  ``Geobox: null`` *Geobox*

  ``Priority: '1D'`` *Available options '1D'/'3D'. The 'Priority' parameter indicates the crustal models priority. If '1D' the programs start to searching for  the best suited model (based on Geobox) on the 'Crustals1D' sub-section, and if '3D' on the 'Crustals3D' sub-section.*

  ``Phase: ['P', 'S']`` *Available options 'P'/'S'. Phases to calculate traveltimes.*

  ``Package: 'FMM'`` *Available options 'FMM'/'TAUP'/'NNLOC'. Methods to calculate traveltimes.*

  **Available Methods:** 'FMM' - Fast Marching Method | 'TAUP' - Taup | 'NNLOC' - Eikonal finite-difference.

  ``Save: SSA2py/path/to dir`` *Path to save traveltime tables*

- ``System:`` *Set system paths useful for the program.*

  ``NonNinLoc: path/to dir`` *Set NonNinLoc executables path. For traveltimes calculations.*

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

**Comments about the Travetimes section**

- The user can give multiple 1D/3D candidate velocity models.

- The structure of the given velocity file is simple. For the 1D case we need 3 columns depth(km)/Vp(km/s)/Vs(km/s). For the 3D we give seperately the Vp and Vs info. Each file must contain 3 columns Lon/Lat/Depth/Velocity. Examples are provided.

- For the 3D models only the FMM method is provided.

- The NonNinLoc is optional. You can leave empty path and ingore the NNLOC choice.

- The traveltime calculations are really fast. We suggest to always calculate traveltimes for both P and S. 

- Granularity around 1 km is satisfying for quite detailed calculations. You can set lower values but this could be computational expensive.

- Keep in mind that the Distance+Depth>SSA Grid Size.

- In 3D models SSA2py internally interpolates the velocity, also in areas outside given positions. It is the users responsibility to provide 3D model that complies with the SSA grid.

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**


- ``Download Service:`` *Section dedicated to provide SSA2py with event information, data and metadata.* 

  ``Event Info:`` *This service is used to provide SSA2py with a catalog of seismic events, either for real-time or past-time operation. Currently, only one FDSNWS-event service can only be used here.*

  ``Host: URL``

  ``Inventory:`` *The Inventory refers to the stations' characteristics such as the poles-zeros, digitizer's gain etc. Available options 'FDSNWS'/'StationXML'/'StationYAML' formats. Multiple sources can be used.*

  ``- ['FDSNWS', URL, null]``

  ``- ['StationXML', local/path, null]``

  ``- ['StationYAML', local/path, null]``

  ``Stream:`` *Stream sub-section refers to the waveforms data. Available options 'FDSNWS'/'SeedLink'/'SDS'/'MSEED'. Multiple sources can be used.*

  ``- ['FDSNWS', URL, null]``

  ``- ['MSEED', local/path, null]`` 

  ``- ['SDS', local/path, null]``

  ``- ['SeedLink', URL, null]``

- ``Download Rules:`` *Retrieve channel types based on distance and magnitude rules.*

  ``Distance:`` 
  
  ``- [minMag, maxMag, [minDist (km), maxDist (km), [channels]]]`` *Distance rules.*

  ``- [5.0, 6.0, [0, 80, ['HN']]]`` *For example, if 5.1>= Mag <=6.0, SSA2py will retrieve stations within a range of: 0 km <= loc<= 80 km where loc is the initial location (hypocenter) of the event and only for channel type of HN (Strong motion stations).*

  ``Time: [seconds before, seconds after Origin Time]`` *Traces duration.*

  ``Components: ['Z23', 'Z12', '123', 'ZNE']`` *Accepted types of components. Note that SSA2py uses only E,N,Z components, so for that reason rotation based on metadata will be performed.*

  ``Blacklist: []`` *Blacklisted Stations. Stations/channels mentioned here will be discared. Accepted format network.station (e.g. HL.KLV) or network.station.channel (e.g. HL.KLV.HHE)*

  ``Whitelist: []`` *Whitelisted Stations. Only stations mentioned here will be used.*

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

**Comments about the Download Service and Download Rules sections**

- You can use as many services as you want in the Stream and Inventory subsections and the stations' data/metadata will be merged by looking in FIFO order. The more services you provided, the more time will be consumed for retrieval. If XML/YAML file is used, you are responsible to provide updated meta-data relevant to the seismic events or real-time scenario (e.g. daily inventory update cron-job). In case you want to use federator please provide "eida-routing" or "iris-federator" as second parameter. 

- The third parameter in Inventory and Stream ([.., .., 3]) is used for cases with restricted data. Please provide the token in the third position.

- For realtime operation is recommended a SeedLink service. However, it takes more time than the other services.

- In case of metadata (e.g 'StationYAML') without information about the orientation, only E,N,Z components will kept and user is responsible for the orientation.

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**


- ``Streams:`` *The 'Streams' section is all about the processing of the data, that will be used as input to SSA2py.*

  ``Duration: [seconds before, seconds after origin time]`` *Time duration of the trace (in seconds) before and after the origin time (if smaller pad with zeros).*

  ``Quality Control:`` *The modules help to choose good quality waveforms.*
  
  ``-['TIME', 'SNR', 'CLIP']`` *The possible inputs are: SNR, TIME, CLIP. Use an empty list [] if no modules need to be set.*
  
   **More info about Quality Control** 

   SNR: Real-time signal to noise ratio (based on module's threshold)

   CLIP: Real-time check for clipped waveforms (based on module's threshold)

   TIME: Real-time check for accuracy in instruments timing (based on module's threshold)
  

  ``Resample: [True/False, freq]`` *If True resample to the given frequency - freq (in Hz) all the traces.*

  ``Filter:`` *Bandpass filter for the processed traces*

  ``-[low_freq, high_freq]`` *e.g. '2-8': Lower limit 2, high 8 Hz (>1). Use [0,0] for no filtering.*

  ``Type: 'ENV'`` *Convert wavefroms. Available options 'ENV', 'OBS', 'ABS', 'STALTA', 'KURT'.*

  **Available types:** 'ENV': Envelopes, 'OBS': Observed waveforms, 'ABS' : Absolute part of the Observed (>1), 'STALTA': STALTA characteristic function, 'KURT': Kurtosis characteristic function.

  ``Type Parameters: [input 1, input 2]`` *For 'STALTA' and 'KURT' methods. In case of 'STALTA' - [STA, LTA], 'KURT'   - [window]. If empty default values will be used.*

  ``Rotate: True/False`` *If True Rotate to the Radial-Transverse system.*

  ``Response: True/False`` *Remove Response if True. Poles/zeros info are needed. Response removal is possible only if the user provides STATIONXML file.*

  ``Quantity: 'ACC'`` *Physical quantity of all traces. Choose between 'ACC': Acceleration, 'DISP': Displacement, 'VEL': Velocity. This process is available also for the raw records by integration/differentiation. Create a uniform dataset concerning the physical quantity.*

  ``Normalize: [True/False, input 1, input 2]`` *Normalize the traces if True. Input 1 is the root factor and input 2 the normalization limit.*

   **More info about Normalize**
   
   Input 1: Root factor (usually 1)

   Input 2: 0 - means linear normalization with 2 standard deviation of absolute amplitude being 1, =1-99  means the largest amplitude is set to this number and >=101 means the average of a number of the largest amplitudes (the exact number is set as navemax in the include file) is set to be (normal_type - 100)
  
  ``Corrections: [True/False, input 1, input 2 (if applies)]`` *Correct traces if is True. Input 1 identifies the correction method.*

   **More info about Corrections:**

   Input 1: 0 means time shifts using the P phase and CC. For the identification of the phase, picking algorithm is used. Method as applied in Evangelidis and Kao (2014). 1 means time shifts using the P phase. P picked arrival from external source is needed. 2 means time shifts using the S phase. S picked arrival from external source is needed. 3 means time shifts using directly the shift time. Shift time from external source is needed.

   Input 2: path for text file with manual picks. Format of the file NET.STA UTCDatetime (P arrival) UTCDatetime (S arrival) (For Options 1,2) and NET.STA float (For Option 3).

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

**Comments about the Stream section**

- If you want to use SSA2py for location purposes 'STALTA' or 'KURT' methods can be used.

- Be careful with the corrections bacause can really distrurb the results!

- Normalize of absolute amplitude being 1 and Root factor 1 is usually applied. 

- Recommended for SSA is to use positive data e.g. Envelopes.

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

``Backprojection:`` *Section about the Source-Scanning Settings.*

``GPU: True/False`` *If True the SSA will be performed in the GPU. Otherwise the SSA procedure will be hosted by the CPUS.*

``Grid:``  *Define the SSA grid.*

``- [min. Magnitude, max. Magnitude, ['box', x side length (km), y side length (km), z side starting depth (km), z side ending depth (km),  step (km)]]``

``Selection:``

``Components: []`` *Denife components to apply SSA. Available options E (East-West), N (North-South), Z (Vertical), T (Transverse), R (Radial).*

``Distance: float`` *Maximum distance between epicenter and stations to be selected*

``Sectors: [minSectors, minStations per sector]`` *Minimum sector requirements to run SSA. The parameter has the form of: [minSectors, minStations per sector] 4 sectors of 90 degrees are defined around the given hypocenter. The first number indicates the number of sectors (minSectors), where at least minStations are found, to be considered. The aim of this operator is to prevent SSA executions with weak azimuthal distribution.*

``Settings:`` *Settings for SSA.*

``Phase: []`` *P or S phase for the BP calculation. Must be consistent with the pre-calculated tttables.*

``Mute: [True/False, Input 1, Input 2 (if applies)]`` *Mute the P or S wave part of the waveforms (depending the BP phase), to get rid of artifacts caused by intruding phases.*

 **More infos about Mute:**
 
 Input 1: 0 - Theoretical mute based on traveltimes. 1 - Mute waveforms based on given picks

 Input 2: path - path of file with P and S picks. Format of the file: NET.STA UTCDatetime (P arrival) UTCDatetime (S arrival).

``TTmax: float`` *Only traveltimes <= TTmax will be used in brightness calculation. Otherwise if the traveltime between station-grid point exceeds this value the station will be ignored.*

``ScanningTime:`` *Scanning time before and after the origin time based on the magnitude of the event.*

``- [min. Magnitude, max. Magnitude, [seconds before, seconds after origin time]]`` 

``TimeShift: float`` *Time shift for each step (in sec).*

``MovingWindow: [half window 1, half window 2]`` *Half length of the moving window (in sec).*

``BrType: 0/1`` *Available options: 0 - means average of sum between stations (Honn Kao and Shao-Ju Shan, 2007), 1 - means average of multiplication between stations (Modification based on Roessler et. al, 2010)*

``StaThre: float`` *Only points calculated from (StaThre)% of total stations are kept put 0-1.0 (0-100%). Otherwise grid will be ignored.*

``bthre: float`` *Only points with brightness >= bthre*bmax will be outputted.*

``Weight: 0/1`` *Weighting of the Half length moving window. Available options: 0 - Gaussian weighting, 1 - Equal weighting*

``Nroot: int`` *This will raise the calculated brightness to its N-th root (Check Normalize).*


**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

**Comments about the Backprojection section**

- TTmax parameter is closely connented with the max. distance between station-grid point. The main aim of this is to help the user control the 'infection' of the dataset from secondary seismic phases.

**--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------**

``Tests:``

``Array Response Function: True/False`` *If True the Array Response Function is performed.*

``Uncertainty Parameter: SE/CI/STD`` *Uncertainty parameter to plot. Choose between: CI: Confidence of Interval, SE: Standard Error, STD: Standard Deviation.*

``Jackknife: True/False`` *Perform Jackknife test.*

``Bootstrap: [True/False, int, int]`` *Bootstrap test [True/False, number of repeats, percentage of stations].*

 **More infos about Bootstrap:**

 The Bootstrap test receives as inputs the number of resampling repeats and the percentage of stations taken into account in the reasmpling procedure from each 90 degree sector. For example in case: [False, 40, 50] the program will choose randomly the 50% of stations in each sector and will repeat the test 40 times. Each time the algorithm will select a random number of the 50% excluded stations, these stations will be restored back to the initial dataset and the SSA will be performed.

``Delete: True/False`` *Due to the significant number of binary files that the tests export you have the option to delete these binary files and keep only the simplified .txt files.*


``Plotting:`` *Plotting section*

``Save Layers: path/to dir`` *Path to download basic shapefiles for the plots. Execute SSA2py.py --download.*

``Topography/Bathymetry: [True/False, path/to topo]`` *If True provide a Topography for the plots, netCDF format only!*

.. note:: A nice suggestion for global available data is GEBCO (https://download.gebco.net).

``Plots: True/False`` *If True draw all the analysis plots.*

``Animation: True/False`` *Animation for the SSA results. Keep in mind that the animation function is relatively slow!*

.. warning:: Animating the results can be a seriously slow procedure, escpecially in cases with small timeshift (e.g. 0.1 sec). For real-time cases it is suggested to turn off this function or to adopt small timesteps. 

- ``Monitor:`` *Section about the seismic events monitoring.*

``Magnitudetype: str`` *Type of magnitude (null for all).*

``MagnitudeTresh: float`` *Magnitude threshold.*

``Range: float`` *Check interval in sec.*

``Playback: float`` *Set sec for past-time run.*

``Historical: True/False`` *Use this only if Playback is not zero and you want to reproduce the exact scenario.*

``Geobox: (Lon. 1, Lat. 1), (Lon. 2, Lat. 2), (Lon. 3, Lat. 3), (Lon. 4, Lat. 4)`` *Geobox for seismic events similar to traveltimes.*

.. note:: You can use more compicated coordinate sets, escaping from the 'Box' logic and using 'Polygon' shapes instead.

``Quality:`` *Uncertainty of the seismic event solution.*

``Time: float``

``Depth: float`` 

``Latitude: float``

``Longitude: float``

``Magnitude: float``

``Timeout: float`` *Must get associated with the Range value (<=Range)*

.. raw:: html

    </div>
