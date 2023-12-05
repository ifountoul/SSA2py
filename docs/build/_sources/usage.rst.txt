
.. usage:


==========
Usage
==========

.. raw:: html

    <div class='Usage'>

Terminal Setup
--------------

Enter the SSA2py directory:

``cd SSA2py``

Enable the execution flag to SSA2py.py (only needed once):

``chmod +x ./SSA2py.py``

Run the program:

``./SSA2py.py``

However you can run without execution flag enabled:

``python3 ./SSA2py.py``


Usage Examples
--------------

Before you do anything else make sure that the configuration file is calibrated based to your needs.
You can find more details about the configuration `here <configuration.html>`_.

Execute SSA2py for new events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this subsection you can find examples that could help you understand the terminal-based philosophy 
of SSA2py.

**Examples**

1. Calculate SSA for a single event (info about the event retrieved from FDSNWS-event).

``./SSA2py.py -c config.yaml -e noa2019gfqmt``

Command analysis:

- Analyze the event (**-e**) with id noa2020diego (from FDSNWS-event service). 

- Use the configuration file (**-c**) with the name config.yaml.

2. Calculate SSA for two events and define a new log file.

``./SSA2py.py -c config.yaml -e noa2020diego -e noa2019gfqmt -l new_log`` 


Command analysis:

- Analyze the events (**-e**) with ids noa2020diego, noa2021asfg (from FDSNWS-event service).

- Use the configuration file (**-c**) config.yaml. 

- Bypass the default log file path (**-l**) (./log) and use another one (new_log).

3. Calculate SSA for all the events found in a datetime range (info about the events retrieved from FDSNWS-event).

``./SSA2py.py -c config.yaml --datetime-range 2018-01-01T00:00:00 2019-01-01T00:00:00``

Command analysis:

- Analyze all the events in the timerange (**- - datetime-range**) between 2018-01-01T00:00:00 and 2019-01-01T00:00:00. Restrictions about the events (e.g. Magnitude) are set to the configuration file.

- Use the configuration file (**-c**) with the name config.yaml.

4. Calculate SSA for a single event (without FDSNWS-event).

``./SSA2py.py -c config_new.yaml -e 2021-04-25T15:13:39 3.0 ML 37.24 20.49 4.1``


Command analysis:

- Use the configuration file (**-c**) config_new.yaml

- The events information (**-e**) are given manualy (datetime, magnitude, magnitude type, latitude, longitude, depth)

5. Calculate SSA for a seismic catalog described in QuakeML format file.

``./SSA2py.py --event-xml ./catalog.xml``

Command analysis:

- The events information (**- - event-xml**) are taken from the catalog.xml file.


6. Calculate SSA for a seismic catalog in plain text.

``./SSA2py.py --event-file ./catalog.txt``

Command analysis:

- The events information (**- - event-file**) are taken from the catalog.txt file. Each line of the file can have an event_id or the description as in example 4. 


7. Calculate SSA for a single event (with FDSNWS-event) removing stations from the procedure. 

``./SSA2py.py -c config.yaml -e noa2020diego -s KLV ATH.Z --remove``

Command analysis:

- Execute the event (**-e**) with id noa2020diego.

- Select stations/components (**-s**) KLV (all the components) and ATH.Z (only the vertical component) and remove them (**- - remove**). 

Repeat the computation
~~~~~~~~~~~~~~~~~~~~~~

SSA2py offers the possibility to repeat the SSA computation quickly by removing stations/components 
or even changing the configuration in order to improve the quality of the results. The user with the repeat mode (**- - repeat**)
can get away with the time consuming download part of wavefroms/metadata and the calculation of traveltime tables.

**Example**

``./SSA2py.py -c config.yaml -e noa2020diego --repeat -s KLV ATH.Z --remove``

Command analysis:

- Execute the event (**-e**) with id noa2020diego.

- Repeat the procedure (**- - repeat**) selecting stations/components (**-s**) KLV (all the components) and ATH.Z (only the vertical component) and remove them (**- - remove**).

Listen for upcoming events
~~~~~~~~~~~~~~~~~~~~~~~~~~

The below command needs to be combined with a tool such as the `cron <https://en.wikipedia.org/wiki/Cron>`_ utility to matters.


``./SSA2py.py -c config.yaml --real-time``

Analyze non-cataloged events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SSA2py primarily focuses on analyzing cataloged seismic events. Nevertheless, users have the flexibility to handle seismic data of an agnostic nature through the manual 
execution mode in SSA2py and using some dummy imputs.

**Example**

``./SSA2py.py -c config.yaml -e 2021-04-25T15:13:39 3.0 ML 37.24 20.49 4.1``

``./SSA2py.py -c config.yaml -e --scanning start time-- --dummy magnitude-- ML --central grid latitude-- --central grid longitude-- --dummy depth--``

Of cource you also have to tune the duration of scanning based on your preference the duration of scanning throught the configuration file. 

Cron Job Examples
~~~~~~~~~~~~~~~~~

Every minute listen to an FDSNWS-event for new incoming events.

``SHELL=/bin/bash``

``* * * * * source /home/.bash_profile; cd /home/SSA2py; python3 ./SSA2py.py -c config.yaml --real-time;``

Help Screen
~~~~~~~~~~~

You can see the help screen with the commands:

- Simply type ``./SSA2py.py``

- ``./SSA2py.py -h``



░██████╗░██████╗░█████╗░██████╗░██████╗░██╗░░░██╗
██╔════╝██╔════╝██╔══██╗╚════██╗██╔══██╗╚██╗░██╔╝
╚█████╗░╚█████╗░███████║░░███╔═╝██████╔╝░╚████╔╝░
░╚═══██╗░╚═══██╗██╔══██║██╔══╝░░██╔═══╝░░░╚██╔╝░░
██████╔╝██████╔╝██║░░██║███████╗██║░░░░░░░░██║░░░
╚═════╝░╚═════╝░╚═╝░░╚═╝╚══════╝╚═╝░░░░░░░░╚═╝░░░
SSA2PY: Source Scanning Algorithm in Python

Version: 1.0

License: GPLv3

Author: Ioannis Fountoulakis (ifountoul@noa.gr)

Credits: Ioannis Fountoulakis (ifountoul@noa.gr), Christos Evangelidis
(cevan@noa.gr)

🄯 2023, Institute of Geodynamics - National Observatory of Athens

usage: SSA2PY [-h] [-c [FILEPATH]] [-e EVENT [EVENT ...]] [--event-file FILEPATH]
              [--event-xml FILEPATH] [--real-time] [-d TIME TIME] [-s STA[.NEZ]
              [STA[.NEZ] ...]] [--repeat] [--remove] [--disable-quality] [--download]
              [-l [FILEPATH]] [-v]

SSA2py: Source Scanning Algorithm in Python

Find more info at:

optional arguments:

``-h, --help`` show this help message and exit

``-c [FILEPATH], --config [FILEPATH]`` default configuration file (config.yaml)

``-e EVENT [EVENT ...], --event EVENT [EVENT ...]`` EVENT can be in any of the following formats: (i) DATETIME
MAGNITUDE TYPE LATITUDE LONGITUDE DEPTH e.g.: /home/john/Projects/SSA2PY/./SSA2PY.py -e 2022-04-28T09:03:22.933595 3.0 ML 37.24 20.49 4.1 (ii) DATETIME 
e.g.: /home/john/Projects/SSA2PY/./SSA2PY.py -e 2022-04-28T09:03:22.933599 (iii) EVENTID (event identifiers are data center specific) e.g.: 
/home/john/Projects/SSA2PY/./SSA2PY.py -e noa2020owyrp In cases (ii) and (iii) the rest of the information is retrieved by the FDSNWS-event In more than one results, only the first event is
returned Passing milliseconds is optional.

``--event-file FILEPATH`` parse and run a file with EVENT lines

``--event-xml FILEPATH``  parse and run a file in QuakeML

``--real-time`` invoke --datetime-range for real-time use (FDSN bounded)

``-d TIME TIME, --datetime-range TIME TIME`` invoke SSA computation for all events found in specific datetime range (FDSN bounded)

``-s STA[.NEZ] [STA[.NEZ] ...], --station STA[.NEZ] [STA[.NEZ] ...]`` override default stations selection. Optionally, components could be also specified.
 It can be combined with --remove for the reverse result

``--repeat`` run again the SSA. Suggested to use -s to choose or exclude stations

``--remove`` it can only be used with --station. Invokes reverse result

``--download`` download Faults and Tectonic Plates

``l [FILEPATH], --log [FILEPATH]`` override default main log file (./log)

``-v, --version``  show program's version number and exit


Extra Comments
~~~~~~~~~~~~~~

- Some of the commands can be triggered in long and in short mode (e.g. ``-e`` or ``--event``).

- Event commands (``--event``, ``--event-xml``, ``--event-file``, ``--datetime-range``) can be combined together.

- The Station command ``--station`` can be used multiple times and can be inverted by attaching the ``--remove``. 

- You can use multiple simultaneous cron-jobs in order to run different configuration files. This could be useful if you want to run different setups based on (i) regions,  (ii) quality (sparse/dense 4D SSA grid search), (iii) event catalog (FDSNWS-event service).

.. raw:: html

    </div>
