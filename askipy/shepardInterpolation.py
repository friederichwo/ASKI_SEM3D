# ----------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------
import numpy as np


def shepard_interpolation(x, y, z, cc):
    """
    Calculate interpolation weights for anchor points
    :param x: x-coordinate of target point
    :param y: y-coordinate of target point
    :param z: z-coordinate of target point
    :param cc: coordinates of interpolation anchors (3,na)
    :return: interpolation weights for anchors
    """
    d = np.sqrt((cc[0,:]-x)**2+(cc[1,:]-y)**2+(cc[2,:]-z)**2)
    dmin = np.min(d)
    imin = np.argmin(d)
    dmax = np.max(d)
    nc = d.size
    w = np.zeros(nc)

    #  if point is very close to an anchor point, set that weight to 1 and all others zero
    if dmin/dmax < 1.e-6:
        w[imin] = 1.0
        return w

    #  choose radius of a sphere centered at target containing all anchor points
    #  interpolation weigth = 0 on this radius
    r = 1.5*dmax

    # compute values s_i
    s = np.zeros(nc)
    for i in range(nc):
        if d[i] < r/3.0:
            s[i] = 1.0/d[i]                        # near zero was checked above
        else:
            h = d[i]/r - 1.0
            s[i] = 27*h*h/(4*r)

    # define interpolation weights from s and
    # direction factors t (interpolation function f_3 in the paper)
    sum_s = np.sum(s)
    t = np.zeros(nc)
    for i in range(nc):
        t[i] = 0.0
        for j in range(nc):
            if j == i:
                continue
            h = (cc[0,i] - x)*(cc[0,j] - x) + (cc[1,i] - y)*(cc[1,j] - y) + (cc[2,i] - z)*(cc[2,j] - z)
            t[i] = t[i] + s[j]*(1.0-h/(d[i]*d[j]))
        t[i] = t[i]/sum_s
    w = s*s*(1.0+t)
    w = w/sum(w)
    return w
