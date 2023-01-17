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

import os
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def createDir(path):
    """
    Manipulate directories

    Arguments:
    ---------
    path: str
        Path to the directory
    
    Returns:
    --------
    -

    """
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        #Clear the inside of the directory
        shutil.rmtree(path)
        os.makedirs(path)
    return path

def delete_npy(path, exclude):
    """
    Delete .npy files
 
    Arguments:
    ----------
    path: str
        Path to the .npy files
    exclude: list
        Files to exclude from delete.    

    Returns:
    --------

    """

    #Remove everything except the Max files
    for file_ in os.listdir(path):
        if file_.endswith('.npy') and not file_ in exclude:
            try:
                os.remove(os.path.join(path, file_))
            except:
                pass

    return


class FixPointNormalize(matplotlib.colors.Normalize):
    """ 
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint 
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level" 
    to a color in the blue/turquise range. 
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.21875, clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def colorBarFix():
    """
    Merge colorbars

    """
    colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
    colors_land = plt.cm.terrain(np.linspace(0.25, 1, 200))

    #combine them and build a new colormap
    colors = np.vstack((colors_undersea, colors_land))
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)

    return cut_terrain_map


def write_txt_maxBright(pathin, pathout):
    """
    Write to .txt file the Maximum brightness results.

    Arguments:
    ----------
    pathin: str
        Path to .npy file
    pathout: str
        Path to write the file   

    Returns:
    --------
    -

    """

    # Read the SSA results
    res = np.load(pathin)

    # dynamically build a format string using the maximum number of digits
    digits_t = len(str(max(res[:,-1])))
    digits_br = len(str(max(np.around(res[:,0],3))))
    digits_lon = len(str(max(np.around(res[:,1],3))))
    digits_lat = len(str(max(np.around(res[:,2],3))))
    digits_dep = len(str(max(np.around(res[:,3],3))))

    f_ = '{0:>%s} {1:>%s} {2:>%s} {3:>%s} {4:>%s}\n'\
         % (digits_t, digits_br, digits_lon, digits_lat, digits_dep)
  
    #Open the .txt file
    with open(pathout, 'w+') as f:
        f.write('{}|{}|{}|{}|{}\n'.format('Time (s)', 'Max. Bright', 'Lon. (°)', 'Lat. (°)', 'Dep. (km)'))
        for r in res:
            f.write(f_.format(str(np.around(r[-1],2)).ljust(digits_t,"0"),\
                            str(np.around(r[0],3)).ljust(digits_br,"0"),\
                            str(np.around(r[1],3)).ljust(digits_lon,"0"),\
                            str(np.around(r[2],3)).ljust(digits_lat,"0"),\
                            str(np.around(r[3],3)).ljust(digits_dep,"0")))
    return

def write_txt_Tests(pathin, pathout):
    """
    Write to .txt file Uncertainty results

    Arguments:
    ----------
    pathin: str
        Path to .npy file
    pathout: str
        Path to write the file   

    Returns:
    --------
    -

    """

    # Read the SSA results
    res = np.load(pathin)
    
    # dynamically build a format string using the maximum number of digits
    digits_t = len(str(max(res[:,-1])))
    digits_br = len(str(max(np.around(res[:,0],3))))
    digits_lon = len(str(max(np.around(res[:,4],3))))
    digits_lat = len(str(max(np.around(res[:,5],3))))
    digits_dep = len(str(max(np.around(res[:,9],3))))
    digits_CI = len(str(max(np.around(res[:,1],1)))) 
    digits_SE = len(str(max(np.around(res[:,2],1))))
    digits_STD = len(str(max(np.around(res[:,3],1))))

    
    f_ = '{0:>%s} {1:>%s} {2:>%s} {3:>%s} {4:>%s} {5:>%s} {6:>%s} {7:>%s} {8:>%s} {9:>%s} {10:>%s} {11:>%s} {12:>%s} {13:>%s}\n'\
         % (digits_t, digits_br, digits_CI, digits_SE, digits_STD, digits_lon, digits_lat, digits_CI, digits_SE, digits_STD, digits_dep, digits_CI, digits_SE, digits_STD)

    #Open the .txt file
    with open(pathout, 'w+') as f:
        f.write('{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|\n'.format('Time (s)', 'Max. Bright', 'CI', 'STD', 'SE', 'Lon. (°)',\
               'Lat. (°)', 'CI (°)', 'STD (°)', 'SE (°)', 'Dep. (km)', 'CI (km)', 'STD (km)', 'SE (km)'))
        for r in res:
            f.write(f_.format(str(np.around(r[-1],2)).ljust(digits_t,"0"),\
                            str(np.around(r[0],3)).ljust(digits_br,"0"),\
                            str(np.around(r[1],3)).ljust(digits_CI,"0"),\
                            str(np.around(r[2],3)).ljust(digits_SE,"0"),\
                            str(np.around(r[3],3)).ljust(digits_STD,"0"),\
                            str(np.around(r[4],3)).ljust(digits_lon,"0"),\
                            str(np.around(r[5],3)).ljust(digits_lat,"0"),\
                            str(np.around(r[6],3)).ljust(digits_CI,"0"),\
                            str(np.around(r[7],3)).ljust(digits_SE,"0"),\
                            str(np.around(r[8],3)).ljust(digits_STD,"0"),\
                            str(np.around(r[9],3)).ljust(digits_dep,"0"),\
                            str(np.around(r[10],3)).ljust(digits_CI,"0"),\
                            str(np.around(r[11],3)).ljust(digits_SE,"0"),\
                            str(np.around(r[12],3)).ljust(digits_STD,"0")))
    return



