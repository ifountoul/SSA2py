#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright (C) 2022 Ioannis Fountoulakis, Christos Evangelidis

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

import numpy as np
from math import ceil

# Obspy Imports
################
from obspy.geodetics.base import kilometer2degrees

"""
Collection of python functions in order to create and manipulate the Grids and the tttables for the SSA2py.

"""

def grid_box(x_side, y_side, z_side1, z_side2, incr, epic):
    """
    This function will create the 3D box grid.

    Arguments:
    ------ 
    x_side: float 
        The x side in km
    y_side: float 
        The y side in km
    z_side: float 
        Depth of the grid in km
    incr: float
        Increment of the grid nodes
    epic: list 
        List with the lat,lon of the epicenter

    Returns:
    ------
    Numpy matrices with the grid in coordinates

    in tuple format
    """

    #Leave only two decimal in epicenter
    lat = round(epic[0],2); lon = round(epic[1],2)
   
    #x-y-z sides (keep them in coordinates for now)
    x = np.linspace(lon-kilometer2degrees(round(x_side/2)),\
                    lon+kilometer2degrees(round(x_side/2)),\
                    round(x_side/incr)+1)

    y = np.linspace(lat-kilometer2degrees(round(y_side/2)),\
                    lat+kilometer2degrees(round(y_side/2)),\
                    round(y_side/incr)+1)

    z = np.linspace(z_side1, z_side2, ceil((z_side2-z_side1)/incr) + 1)

    x, y = np.meshgrid(x, y, indexing = 'ij', sparse = False)        
    
    #Flat the array
    x, y = np.vstack([x.ravel(), y.ravel()]) 

    return (x, y, z)


def grid_box3D(x_side, y_side, z_side1, z_side2, incr, epic, distance):
   """
   Box grid for the 3D implementation

   Arguments:
   ------
   x_side: float 
        The x side in km
   y_side: float 
        The y side in km
   z_side: float 
        Depth of the grid in km
   incr: float
        Increment of the grid nodes (must be the same with the granularity in tt settings)
   epic: list 
        List with the lat,lon of the epicenter
   distance: float
        Distance setting in the tt calculations

   Returns:
   ------
   Numpy matrices with the grid in coordinates
   and positions of the grids in the initial tt array

   in tuple format  

   """

   #Leave only two decimal in epicenter
   lat = round(epic[0],2); lon = round(epic[1],2)

   #Reconstruct the grid from the tt calculations
   gx = np.linspace(lon-kilometer2degrees(round(distance)),\
                    lon+kilometer2degrees(round(distance)),\
                    round(distance*2/incr)+1)

   gy = np.linspace(lat-kilometer2degrees(round(distance)),\
                    lat+kilometer2degrees(round(distance)),\
                    round(distance*2/incr)+1)

   #keep only the grid points inside the given box size
   posx = np.where((gx>=lon-kilometer2degrees(round(x_side/2)))\
          & (gx<=lon+kilometer2degrees(round(x_side/2))))[0]

   posy = np.where((gy>=lat-kilometer2degrees(round(y_side/2)))\
          & (gy<=lat+kilometer2degrees(round(y_side/2))))[0]

   x = gx[posx]
   y = gy[posy]

   z = np.linspace(z_side1, z_side2, ceil((z_side2-z_side1)/incr) + 1)

   x, y = np.meshgrid(x, y, indexing = 'ij', sparse = False)

   #Flat the array
   x, y = np.vstack([x.ravel(), y.ravel()])

   return (x, y, z, posx, posy)


def grid_points(file_path):
    """
    Read the .txt file with the grid points

    Arguments:
    ---------
    file_path: string
         Path of the .txt file
    Returns:
    --------
    x, y, z: list
         lists with values from every column.

    """
    try:
        with open(file_path, 'r') as file:
            x_vals, y_vals, z_vals = [], [], []
            for line in file:
                x, y, z = line.strip().split()
                x_vals.append(float(x))
                y_vals.append(float(y))
                z_vals.append(float(z))
            return x_vals, y_vals, z_vals
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except ValueError:
        print(f"Error: Invalid data format in file '{file_path}'.")

