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
#   Read fault lines from shapefile

import shapefile as shp
import numpy as np

def getFaultDict(shapefile):
    """
    Read shapefile, extract fault lines and fill a dictionary
    used as kwargs to plotVslice.
    :param shapefile: File with fault line specifications
    :return:
    """
    sf = shp.Reader(shapefile)
    ls = {1: 'solid', 2: 'solid', 3: 'dashed'}
    lw = {1: 1.6, 2: 0.8, 3: 1.6}
    #
    #  faultdict is a dictionary containing as key the order index (oid) with value a dictionary containing:
    #  (lon,lat)-point coordinates as Nx2-array, lineweight and linestyle
    #
    faultdict = dict()
    for el in sf.shapeRecords():
        ft = el.record['fault_type']
        faultdict[el.shape.oid] = {'lonlat': np.array(el.shape.points), 'lw': lw[ft], 'ls': ls[ft]}
    return faultdict

