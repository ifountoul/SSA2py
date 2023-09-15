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

import numpy as np
import os

# Obspy Imports
###############

from obspy.geodetics.base import gps2dist_azimuth

def normalize_value(value, min_value, max_value, min_normalized, max_normalized):
    """
    Normalize Weights

    """
    normalized_value = ((value - min_value) / (max_value - min_value)) * (max_normalized - min_normalized) + min_normalized
    return normalized_value

def density_weights(stations, threshold_distance, max_weight=5, min_weight=1, min_norm=1, max_norm=3):
    """
    Calculate weights to reduce stations density effects in the SSA solution.
    The density is calculated by identifing the number of neighbour stations that each station has.

    You have to give the minimum and maximum weights.

    Input:
    -------
    stations: List
        List with stations coordinates (Lon, Lat).
    threshold_distance: float
        Distance to count stations.
    max_weight: float
        Maximum weight.
    min_weight: float
        Maximum weight.
    min_norm: float
        Minimum weight to normalize.
    max_norm: float
        Maximum weight to normalize.

    Output:
    -------
    weights: array
         Array with station weights.

    """

    weights = []
    num_stations = len(stations)

    for i in range(num_stations):
        station = stations[i]
        num_nearby_stations = 0

        for j in range(num_stations):
            if i != j:
                distance = gps2dist_azimuth(station[1], station[0], stations[j][1], stations[j][0])[0]/1000
                if distance <= threshold_distance:
                    num_nearby_stations += 1
        # Adjust the weight based on the nearby station count
        weight = max(min(max_weight - num_nearby_stations, max_weight), min_weight)
        weights.append(weight)

    normalized_values = [normalize_value(w, np.min(weights), np.max(weights), min_norm, max_norm) for w in weights]

    return np.array(normalized_values)

def dist_weights(stations, evla, evlo, min_weight=1, max_weight=3):
    """
    Calculates weights of stations based on distance.
    Exponential weights.

    Input:
    ------
    stations: Array
        Array with stations coordinates.
    evla: float
        Event Latitude.
    evlo: float
        Event Longitude.
    min_weight: float
        Minimum Weight
    max_weigh: float
        Maximum Weight
 
    Output:
    -------
    weights: Array
        Distance weights

    """

    distances = []
    for sta in stations:
         dist, azim, _ = gps2dist_azimuth(evla, evlo, sta[1], sta[0])
         distances.append(dist)

    min_distance = np.min(distances)
    weights = np.exp(-0.2 * (distances - min_distance) / min_distance)

    normalized_values = [normalize_value(w, np.min(weights), np.max(weights), min_weight, max_weight) for w in weights]

    return np.array(normalized_values)


def combine_weights(wden, wdist, min_w = 1, max_w = 3):
   """
   Combine Weights and normalize.

   Input:
   ------
   wden: Array
       Density weights.
   wdist: Array
       Distance weights.

   Output:
   -------
   normalized_values: Final Weights

  
   """

   weights = wden*wdist
   normalized_values = [normalize_value(w, np.min(weights), np.max(weights), min_w, max_w) for w in weights]
   return np.array(normalized_values)


def read_external_weights(path, st):
    """
    Read weights from external txt file.

    Input:
    ------
    path: str
        Path to txt file
    st: Obspy Stream
        Obspy Stream Object

    Return:
    -------
    w: Array
       External weights

    """
     
    if os.path.exists(path): 
       f = open(path, 'r')
       lines = f.readlines()
       f.close()
    else:
        return

    w = []
    for tr in st.select():
        station_name_ = tr.stats.network+'.'+tr.stats.station
        # Added weight
        add=False
        for line in lines:
            station_name = line.split()[0]
            if station_name==station_name_:
                w.append(float(line.split()[1]))
                add=True
        if add==False:
            w.append(float(1)) 
    
    return np.array(w)


def azimuth_weights(stations, evla, evlo):
    """
    Correcting for non-uniform station distribution based on Jakka et al, 2010
 
    Input:
    ------
    stations: Array
        Array with stations coordinates.
    evla: float
        Event Latitude.
    evlo: float
        Event Longitude.

    Output:
    -------
    weights: Array
        Array with weights.

    """

    def calculate_angle_differences(angles, reference_angle):
        """
        Calculate the differences between an array of angles and a reference angle.
        The result is always in the range [-180, 180) degrees.
        """
        differences = np.degrees(np.arctan2(np.sin(np.radians(angles - reference_angle)), np.cos(np.radians(angles - reference_angle))))
        return differences

    # Backazimuth of stations
    backazim = np.array([gps2dist_azimuth(evla, evlo, sta[1], sta[0])[2] for sta in stations])

    # Sort the backazimuths
    sorted_indices = np.argsort(np.array(backazim))
    backazim = backazim[sorted_indices]
    stations = stations[sorted_indices]

    # Define backazimuth of bisectors   
    bisectors = [ (backazim[i] + abs(backazim[i+1] - backazim[i])/2) for i in range(len(stations)-1)]

    # Define weights
    weights = []

    for i in range(len(stations)):
        f = 0
        for j in range(len(bisectors)):
             # Bottom
            fjk = 0
            for i_ in range(len(stations)):
                fjk+=abs(calculate_angle_differences(bisectors[j], backazim[i_]))
            fji = abs(calculate_angle_differences(bisectors[j], backazim[i]))

            f+=fji/fjk
        weights.append(f)
    weights = np.array(weights)
    return weights[np.argsort(sorted_indices)] 
