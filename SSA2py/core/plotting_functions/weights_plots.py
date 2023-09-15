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
##########

import matplotlib.pyplot as plt
import os
from matplotlib.patches import Circle
from matplotlib.cm import viridis


def plot_distance_weights(weights, stations, path='.', dpi=400, fformat='png'):
    """
    Plot distance weights.

    Input:
    ------
    weights: Array
         Array with weights.
    station: Array
         Array with stations Lon. Lat.
    path: str
         Path to save figure.
    dpi: int
    fformat: str
       Type of file to save.
    
    Output:
    -------
    Saved figure in the give path.

    """

    x = stations[:, 0]
    y = stations[:, 1]
    z = weights

    fig, (ax) = plt.subplots(nrows=1)

    ax.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
    cntr = ax.tricontourf(x, y, z, levels=14, cmap="viridis")

    # Plot stations
    ax.scatter(x, y, color='r', alpha=0.5, label='Stations')

    # Add colorbar
    cbar = fig.colorbar(cntr, ax=ax)
    cbar.set_label('Weight Value', rotation=270, labelpad=15)

    ax.set_title('Distance Weights', fontweight='bold', fontsize=12)
    ax.set_aspect('equal')
    ax.set_xlabel('Longitude (°)', fontweight='bold')
    ax.set_ylabel('Latitude (°)', fontweight='bold')

    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(path, 'distance_weights.' + fformat), dpi=dpi)

    return


def plot_density_weights(weights, stations, radius, path='.', dpi=400, fformat='png'):
    """
    Plot the density weights, that have been identified by the neighbours check.

    Input:
    ------
    weights: array
        Array of weights.
    stations: array
        Array of stations Lon., Lat.
    radius: float
        Radius of each check.
    path: str
        Path to save.
    dpi: int
    fformat: str
       Type of file to save.
 
    Output:
    -------
    saves figure to directory.
  

    """


    fig, ax1= plt.subplots(1, 1)

    x = stations[:,0]
    y = stations[:,1]
    z = weights

    circle_radius_km = radius
    circle_radius_data = circle_radius_km / 111

    # Normalize the z values between 0 and 1
    norm = plt.Normalize(min(z), max(z))

    for xi, yi, zi in zip(x, y, z):
        circle = Circle((xi, yi), circle_radius_data, edgecolor='black', facecolor=viridis(norm(zi)), alpha=0.2, ls='--')
        ax1.add_patch(circle)
    scatter = ax1.scatter(x, y, c=z, cmap='viridis')

    # Add colorbar
    cbar = fig.colorbar(scatter, ax=ax1)
    cbar.set_label('Weight Value', rotation=270, labelpad=15)

    ax1.set_title('Density Weights', fontweight='bold', fontsize=12)
    ax1.set_aspect('equal')
    ax1.set_xlabel('Longitude (°)', fontweight='bold')
    ax1.set_ylabel('Latitude (°)', fontweight='bold')

    plt.grid()
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'density_weights.' + fformat), dpi=dpi)

    return


def plot_final_weights(weights, stations, path='.', dpi=400, fformat='png'):
    """
    Plot final weights (Distance*Density).

    Input:
    ------
    weights: Array
        Weights array.
    stations: Aray
        Stations array.
    path: str
        Path to save figure
    dpi: int
    fformat: str
       Type of file to save.

    Output:
    -------
    Saved figure.
  

    """

    x = stations[:, 0]
    y = stations[:, 1]
    z = weights

    fig, (ax) = plt.subplots(nrows=1)

    ax.tricontour(x, y, z, levels=10, linewidths=0.5, colors='k', linestyles='--')
    cntr = ax.tricontourf(x, y, z, levels=10, cmap="viridis")

    # Plot stations
    ax.scatter(x, y, color='r', alpha=0.5, label='Stations')

    # Add colorbar
    cbar = fig.colorbar(cntr, ax=ax)
    cbar.set_label('Weight Value', rotation=270, labelpad=15)

    ax.set_title('Weights', fontweight='bold', fontsize=12)
    ax.set_aspect('equal')
    ax.set_xlabel('Longitude (°)', fontweight='bold')
    ax.set_ylabel('Latitude (°)', fontweight='bold')

    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(path, 'Final_weights.' + fformat), dpi=dpi)


    return


