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

import math, numba
import numpy as np
from numba import jit, prange

numba.config.THREADING_LAYER = 'safe'

@jit(nopython=True, parallel=True)
def _1Ddist(grid, sx, sy):
    """
    Calculate the distance between grid points and stations.

    Input:
    ------
    grid: numpy array
        grid of SSA
    sx,sy: numpy arrays
        x and y positions of stations

    Output:
    ------
    dist_: array
        Table with the distance per grid point ---> stations

    from https://www.johndcook.com/blog/2018/11/24/spheroid-distance/
    Modified for Numba 

    """
    #Pre-allocate the dist array
    dist_ = np.zeros((len(grid), len(sx)))

    for g in prange(len(grid)):
        for s_ in prange(len(sx)):
            long1 = math.radians(grid[g,0])
            lat1 = math.radians(grid[g,1])

            long2 = math.radians(sx[s_])
            lat2 = math.radians(sy[s_])
           
            a = 6378137.0 # equatorial radius in meters 
            f = 1/298.257223563 # ellipsoid flattening 
            b = (1 - f)*a
            tolerance = 1e-11 # to stop iteration

            phi1, phi2 = lat1, lat2
            U1 = np.arctan((1-f)*np.tan(phi1))
            U2 = np.arctan((1-f)*np.tan(phi2))
            L1, L2 = long1, long2
            L = L2 - L1

            lambda_old = L + 0

            while True:
                t = (np.cos(U2)*np.sin(lambda_old))**2
                t += (np.cos(U1)*np.sin(U2) - np.sin(U1)*np.cos(U2)*np.cos(lambda_old))**2
                sin_sigma = t**0.5
                cos_sigma = np.sin(U1)*np.sin(U2) + np.cos(U1)*np.cos(U2)*np.cos(lambda_old)
                sigma = np.arctan2(sin_sigma, cos_sigma)

                sin_alpha = np.cos(U1)*np.cos(U2)*np.sin(lambda_old) / sin_sigma
                cos_sq_alpha = 1 - sin_alpha**2
                cos_2sigma_m = cos_sigma - 2*np.sin(U1)*np.sin(U2)/cos_sq_alpha
                C = f*cos_sq_alpha*(4 + f*(4-3*cos_sq_alpha))/16

                t = sigma + C*sin_sigma*(cos_2sigma_m + C*cos_sigma*(-1 + 2*cos_2sigma_m**2))
                lambda_new = L + (1 - C)*f*sin_alpha*t
                if abs(lambda_new - lambda_old) <= tolerance:
                     break
                else:
                     lambda_old = lambda_new
            u2 = cos_sq_alpha*((a**2 - b**2)/b**2)
            A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320 - 175*u2)))
            B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
            t = cos_2sigma_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m**2))
            t -= (B/6)*cos_2sigma_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2sigma_m**2)
            delta_sigma = B * sin_sigma * t
            s = b*A*(sigma - delta_sigma)

            dist_[g,s_] = s/1000
    return dist_
