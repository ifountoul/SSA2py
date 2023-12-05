RIDGECREST 2019 EARTHQUAKE
--------------------------

This case study focuses on the 2019 Ridgecrest earthquake, which occurred north and northeast of the town of Ridgecrest, California, United States. 
It serves as an illustrative example to showcase the application of the SSA2py software.

Under the 'Case_Studies/Ridgcrest_2019' directory you can find:
--> configRidgecrest.yaml: Configuration file
--> vmodel.txt: Velocity model of the area.
--> stations.txt: Stations to use

Files description:
------------------
vmodel.txt: Velocity model for the area produced by combining various models (Lin et al. 2007; Zhang and Lin 2014; White et al. 2021)
            The format of the file is: 1st column depth (km), 2nd Vp (km/s) and 3rd column Vs (km/s)
stations.txt: Stations with True will be used in the analysis.


How to use these file:
----------------------
Open the configParkfield.yaml configuration file.
In the Traveltimes/Crustals1D/Filename add the vmodel.txt path
In the Download Rules/Stationslist --> [True, path/to/stations.txt]

Execute in the SSA2py Directory: ./SSA2py.py -c configRidgecrest.yaml -e 2019-07-06T03:19:52.000000Z 7.1 Mw 35.766 -117.605 8
Check results in SSA2py/events/2019-07-06T03:19:52.000000Z directory.


References:
-----------

Lin, G., Shearer, P. M., Hauksson, E., and Thurber, C. H. (2007), A three-dimensional crustal seismic velocity model for southern California from a composite event method,
J. Geophys. Res., 112, B11306, doi:10.1029/2007JB004977. 

Zhang, Q., & Lin, G. (2014). Three-dimensional Vp and Vp/Vs models in the Coso geothermal area, California: Seismic characterization of the magmatic system. 
Journal of Geophysical Research: Solid Earth, 119, 4907–4922. https://doi.org/10.1002/2014JB010992

White, M. C. A., Fang, H., Catchings, R. D., Goldman, M. R., Steidl, J. H., & Ben-Zion, Y. (2021). Detailed traveltime tomography and seismic
catalogue around the 2019 Mw7.1 Ridgecrest, California, earthquake using dense rapid-response seismic data. Geophysical Journal Interna-
tional, 227(1), 204–227. https://doi.org/10.1093/gji/ggab224
