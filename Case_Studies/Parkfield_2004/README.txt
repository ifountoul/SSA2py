PARKFIELD 2004 EARTHQUAKE
------------------------

This case study focuses on the 2004 Parkfield earthquake, which took place near the town of Parkfield in California, United States. 
It serves as an illustrative example to showcase the application of the SSA2py software.

Under the 'Case_Studies/Parkfield_2004' directory you can find:
--> configParkfield.yaml: Configuration file
--> correction.txt: Text file that contains P and S arrival times for each station, picked manually.
--> data.mseed: Seismic data in mseed format.
--> Parkfield.yaml: Inventory file that contains information about the seismic stations. 
--> vmodel.txt: Velocity model of the area.


Files description:
------------------
vmodel.txt: Velocity model for the area used in Custodio et al. 2005. 
            The format of the file is: 1st column depth (km), 2nd Vp (km/s) and 3rd column Vs (km/s)

Parkfield.yaml: Inventory file with format     36450:             --> Station Name
                                                  station: 36450  --> Station Name
                                                  network: GS     --> Network
                                                  channels:       --> Three channels
                                                  - HNZ           
                                                  - HNN
                                                  - HNE
                                                  coords:
                                                  - 35.77         --> Station Latitude
                                                  - -120.24       --> Station Longitude
                                                  - 0.0           --> Station Elevation
                                               ...... Another station 

correction.txt: P and S arrivals used in the traveltimes corrections. 
                The format of the file is: 1st station network. station name, 2nd UTCDateTime P arrival and 3rd column UTCDateTime S arrival.


How to use these file:
----------------------
Open the configParkfield.yaml configuration file. 
In the Traveltimes/Crustals1D/Filename add the vmodel.txt path
In the Download Service/Inventory --> - ['StationYAML', path of the Parkfield.yaml, null]
In the Download Service/Stream --> - ['MSEED', path of the data.mseed, null]
In the Streams/Corrections --> [False, 2, path of the correction.txt]

Execute in the SSA2py Directory: ./SSA2py.py -c configParkfield.yaml -e 2004-09-28T17:15:24.000000Z 6.0 Mw 35.815 -120.374 8
Check results in SSA2py/events/2004-09-28T17:15:24.000000Z directory.


References:
-----------
Custódio, S., Liu, P., and Archuleta, R. J. (2005), The 2004 Mw6.0 Parkfield, California, earthquake: Inversion of near-source 
ground motion using multiple data sets,Geophys. Res. Lett., 32, L23312, doi:10.1029/2005GL024417. 
