.. installation:

============
Installation
============

.. raw:: html

    <div class='Basic Info'>


Download and/or clone repository
--------------------------------

You can always clone the code and docs by running:

``git clone https://github.com/ifountoul/SSA2py.git``

Hardware Requirements
---------------------

- The software is developed in a manner to be able to be deployed from PC/Laptop to 
  High-Performance Servers with or without a GPU device.
- Currently, it supports only Linux OS. Tested on Ubuntu 22.04.2 LTS.
- The GPU device must support CUDA GPU programming standard.

Prerequisites
-------------

- Install **conda** (`here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_). Tested on 4.12.0.

- (*Optional*) Install **NonNinLoc** (`here <http://alomax.free.fr/nlloc/>`_). For traveltime tables calculation.

.. note:: In order to succesfully compile NonNinLoc you will need a C compiler, such as gcc.
          

Installation
------------

1. Download and extract the `source code <https://github.com/ifountoul/SSA2PY/archive/refs/heads/main.zip>`_. 

2. Setup Python3 and the associated packages on a conda environment. 
   You will need to source the environment.yml inside your conda workspace. 
   From the terminal run:

   ``cd SSA2py``

   ``conda env create -f environment.yml``

3. Activate conda environment.

   ``conda activate SSA2py``

4. Create the default events (Events Dir) and traveltimes (Traveltimes/save) directories
   as defined in the default config file.  

5. Download files necessary for plotting.

   ``python3 SSA2py.py --download``

Test Functionality
------------------
Evaluate the setup by running in terminal:

- Check the help screen

  ``python3 SSA2py.py`` 

- Run a benchmark event (usage in detail `here <usage.html>`_) 

  ``python3 SSA2py.py -c config.yaml -e noa2019gfqmt`` 

.. raw:: html

    </div>
