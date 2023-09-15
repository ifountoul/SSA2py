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
#########

import os
from git import Repo

def FigDown(path):
    """
    Download for the figures:
    1) Plate margins (https://github.com/fraxen/tectonicplates)
    2) Faults (https://github.com/GEMScienceTools)
    3) Topography

    Arguments:
    ------
    path: str
        Path to directory.

    """
    #Create the plot directory if does not exists
    try:
        os.makedirs(os.path.join(path,'FAULTS'))
        os.makedirs(os.path.join(path,'PLATES'))
    except:
        pass

    try:
        Repo.clone_from('https://github.com/fraxen/tectonicplates.git',\
                        os.path.join(path,'PLATES'))
    except:
        pass

    try:
        Repo.clone_from('https://github.com/GEMScienceTools/gem-global-active-faults.git',\
                         os.path.join(path,'FAULTS'))
    except:
        pass

    return
