Example Applications
********************

Minor seismic event
--------------------------------

We present an application of the SSA method to a minor seismic event (3.01 MLh) that occured in the Gulf of Corinth (Greece).
Several studies until now have demostrated that SSA can be used to identify and image rupture details even in minor and moderate magnitude seismic events. 
Basic requirement is the existence of dense seismic network around the source area. 

Suggested works for SSA applications to moderate - minor events:

- `An Atypical Shallow Mw 5.3, 2021 Earthquake in the Western Corinth Rift (Greece) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2022JB024221>`_ by Jiří Zahradník et al. (2022)

- `Backprojection Imaging of the 2020 Mw 5.5 Magna, Utah, Earthquake Using a Local Dense Strong‐Motion Network <https://pubs.geoscienceworld.org/ssa/srl/article/92/2A/640/592979/Backprojection-Imaging-of-the-2020-Mw-5-5-Magna>`_ by Maria Mesimeri et al. (2021) 

- `Rupturing of small natural earthquakes in West Bohemia investigated by source scanning <https://link.springer.com/article/10.1007/s10950-021-10043-y>`_ by Lávička and Fischer (2022)

In this case we used envelopes of high-frequency (2-8 Hz) S waves on the E-W horizontal component. 
at 13 broadband stations located at distances unto 100 km from the epicenter. 
The specific example can be repeated from the user with the command ``./SSA2py.py -c config.yaml -e noa2023phuht``. 
The used velocity model is from `Rigo et al. (1996) <https://academic.oup.com/gji/article/126/3/663/635741?login=false>`_. 
Time shift corrections have been applied.

- Stations Used

.. image:: atlas_mod.png
   :scale: 13%

*Description: Stations used in the SSA analysis.*

- Maximum Brightness Results

.. image:: MaximumBrightness_mod.png
   :scale: 13%

*Description: Maximum brightness locations in time. The time interval between two consecutive bright spots is 0.1 s. 
The colored plotted circles at the mapview and the vertical cross sections represent maximum brightness location in each time step, sized proportionally to their values.*

.. image:: MaximumBrightnessAllTimesteps.png
   :scale: 13%

*Description: Composite maximum brightness plot of all time steps. Mapview and vertical cross sections.*

.. image:: MaximumBrightnessPerTimeRange_1.png
   :scale: 11%

*Description: Composite maximum brightness plot of given time ranges. Mapview observations.*

- Waveforms (Envelopes)

.. image:: RecordSection_mod.png
   :scale: 25%

*Description: Backprojected waveforms.*

- ARF Results

.. image:: ARF_1_mod.png
   :scale: 8%
 
.. image:: ARF_2_mod.png
   :scale: 8%

*Description: Array Response Function Results.*


SSA2py Publication Case Studies
--------------------------------

You can replicate the SSA results for the Parkfield (2004) and Ridgecrest (2019) earthquakes by following the detailed tutorials in the SSA2py `repo <https://github.com/ifountoul/SSA2py>`_.
First, you have to unzip the Case_studies.zip file and after that to read carefully the README files in each case.   




