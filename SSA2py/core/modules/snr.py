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

import obspy
import numpy as np

def SNRampl(trace, start_spot, a_time=10):
   """
   Simple SNR based on amplitude.

   Arguments:
   -------
   trace: Obpsy trace
   start_spot: int
       trigger position
   a_time : float 
       time around the trigger for the SNR calculation (sec)

   Returns:
   ------
   SNR: float
       SNR value
   """
   df = trace.stats.sampling_rate

   #Translate the seconds in number of points
   slice_pos = int(a_time * df)
   noise_window = trace.data[(start_spot-slice_pos):start_spot]
   signal_window = trace.data[start_spot: (start_spot + slice_pos)]

   #Calculate the SNR(Signal to Ratio) by dividing the RMS amplitude values
   Srms = np.sqrt(1/len(signal_window) * np.sum(np.square(signal_window)))
   Nrms = np.sqrt(1/len(noise_window) * np.sum(np.square(noise_window)))

   SNR = Srms/Nrms

   return SNR
