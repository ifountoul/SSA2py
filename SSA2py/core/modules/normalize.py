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


np.seterr(divide='ignore', invalid='ignore')

def normalize(trace, n_root_factor = 1, norm_type = 0):
   """
   In-between the pre-processing and the processing level function.
   For each trace the Nth-root is taken and the trace is normalized.

   Arguments:
   ---------- 
   n_root_factor: float (usually 1)
   norm_type: int
       -1 = no normalization
       0 = linear normalization with 2 S.D. as 1
       >=1  linear normalization with the maximum set to this value
       * for tremor study, set to 0

   Return:
   -------
   trace: Obspy trace object
       Normalized Trace
   """

   #Average of the trace
   yav = np.average(trace.data)

   #SD of the trace
   ysd = np.std(trace.data)

   #Take the N-root
   if n_root_factor>1:
      # absolute
      trace.data = np.abs(trace.data)
      # n-root
      trace.data = trace.data ** (1/n_root_factor)
  
   #Normalize the trace 
   if norm_type == 0:
      ynorm = 2 * ysd
      trace.data = (trace.data - yav) / ynorm
   elif norm_type>=1:
      trace.data = trace.normalize().data * norm_type 
   else:
      #No normalization
      pass

   return trace
