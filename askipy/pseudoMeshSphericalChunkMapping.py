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


def map_spherical_chunk_to_pseudo_mesh(ccs, re):
    """
    Maps Cartesian coordinates in a spherical chunk to
    coordinates in the pseudo mesh.
    ccs: Cartesian coordinates in the spherical chunk (array)
    returns: Cartesian coordinates in pseudo mesh (array)
    """
    xs, ys, zs = ccs
    r = np.sqrt((zs + re)**2+xs**2+ys**2)
    lat = np.arcsin(ys/r)
    rp = r * np.cos(lat)
    lon = np.arcsin(xs/rp)
    x = re*lon
    y = re*lat
    z = r-re
    return np.array([x, y, z])
