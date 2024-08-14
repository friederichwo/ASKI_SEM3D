#----------------------------------------------------------------------------
#   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
#

import numpy as np
from scipy.special import erf
from matplotlib.colors import ListedColormap
from askipy.regularSphericalGrid import RegularSphericalGrid

def getHsliceColsumDict(filename, sldict, csrsg, *args):
    """
    In first call read column sums 3D file, create rsg-object
    and produce color map. Also initalize csdict.
    In further calls, modify array values in csdict.
    Colormap:
    Column sums values are normalized to 1. Given a threshold w,
    alpha = erf(3*(w-cs)/w) for cs < w and alpha = 0.0 for cs > w.
    Colors RGB are 1 (= white).
    :param filename: name of file from which 3D column sum data are read
    :param sldict: dictionary defining kwargs for hslices overlays
    :param csrsg: regular spherical grid object for 3D colum sums
    :param args: 6 arguments containg values for csbad, ncol, radius, nlat, nlon, slname
    """
    radius = args[2]; nlat = args[3]; nlon = args[4]; slname = args[5]
    if not sldict:
        csbad = args[0]; ncol = args[1]
        rsg = RegularSphericalGrid.readFromHDF(filename)
        hslice = rsg.extractHorizontalSlice(radius, nlat, nlon, slname)
        csval, vmin, vmax = hslice.getReshapedField('colsum', depth_normalize=False, normalize=False,
                                                    cliprel=-1.0, clipabs=-1.0, sym=False, percent=False,
                                                    extrema='None')
        # construct colormap
        vals = np.ones((ncol,4))
        kw = int(csbad*ncol)
        cs = np.linspace(0.0, csbad, kw)
        vals[0:kw,3] = erf(2.0*(csbad-cs)/csbad)/erf(2.0)
        vals[kw:,3] = 0.0
        cmap = ListedColormap(vals)
        sldict = {'array': csval, 'colormap': cmap, 'vmin': 0.0, 'vmax': 1.0, 'shading': 'gouraud'}
        return sldict, rsg
    else:
        hslice = csrsg.extractHorizontalSlice(radius, nlat, nlon, slname)
        csval, vmin, vmax = hslice.getReshapedField('colsum', depth_normalize=False, normalize=False,
                                                    cliprel=-1.0, clipabs=-1.0, sym=False, percent=False,
                                                    extrema='None')
        sldict['array'] = csval
        return sldict


