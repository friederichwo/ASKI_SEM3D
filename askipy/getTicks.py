#----------------------------------------------------------------------------
#   Copyright 2023 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------------
#   Return nicely spaced tick marks for axis annotation
#
import numpy as np


def getTicks(x1,x2,dx):
    """
    Compute intermediate ticks for an axis reaching from x1 to x2
    which are an integer multiple of the spacing dx.
    x1 may be greater than x2. 

    :param x1: value at left/bottom end of axis 
    :param x2: value at right/top end of axis
    :param dx: desired tick spacing 
    :return ticks: list of tick values in specified range
    """
    descending = x1 > x2
    length = abs(x2-x1)
    if descending:
        k1 = int(np.ceil(x2/dx))                        # next integer multiple of dx >= x2
        if abs(x2-(k1-1)*dx)/length < 1.e-3: k1 = k1-1
        k2 = int(np.floor(x1/dx))                       # next integer multiple of dx <= x1
        if abs(x1-(k2+1)*dx)/length < 1.e-3: k2 = k2+1
    else:
        k1 = int(np.ceil(x1/dx))                        # next integer multiple of dx >= x1
        if abs(x1-(k1-1)*dx)/length < 1.e-3: k1 = k1-1
        k2 = int(np.floor(x2/dx))                       # next integer multiple of dx <= x2
        if abs(x2-(k2+1)*dx)/length < 1.e-3: k2 = k2+1

    return [k*dx for k in range(k1,k2+1)]




