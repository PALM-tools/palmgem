#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2018-2024 Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Martin Bures, Jaroslav Resler.
#
# This file is part of PALM-GeM.
#
# PALM-GeM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-GeM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-GeM. If not, see <https://www.gnu.org/licenses/>.

mt=  [[11100, 203, 901],
      [11210, 203, 902],
      [11220, 203, 903],
      [11230, 203, 902],
      [11240, 203, 903],
      [11300, 203, 903],
      [12100, 203, 906],
      [12210, 201, 906],
      [12220, 202, 906],
      [12230, 209, 906],
      [12300, 203, 906],
      [12400, 103, 906],
      [13100, 101, 906],
      [13300, 101, 906],
      [13400, 108, 906],
      [14100, 118, 906],
      [14200, 103, 906],
      [21000, 102, 906],
      [22000, 101, 906],
      [23000, 103, 906],
      [24000, 102, 906],
      [25000, 102, 906],
      [31000, 117, 906],
      [32000, 110, 906],
      [33000, 112, 906],
      [40000, 110, 906],
      [50000, 301, 906],
     ]

def assign_type(mi):
    if mi < 200:
        txt = 'vegetation_type'
    elif mi < 300:
        txt = 'pavement_type'
    elif mi < 400:
        txt = 'water_type'
    elif mi < 500:
        txt = 'soil_type'
    else:
        txt = 'building_type'
    return txt

for m in mt:
    m1_text = assign_type(m[1])
    m2_text = assign_type(m[2])
    print('{}, {}, {}, {}, {}'.format(m[0], m[1], m1_text, m[2], m2_text))
