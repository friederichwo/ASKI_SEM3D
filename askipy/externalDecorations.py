#!/usr/bin/env python3
# ----------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------
"""
Collection of funtions called by addExternalDecoration when plotting
points, lines, or overlays into horizontal or vertical slices
"""
import numpy as np

# global module variables to avoid repeated reading files
csrsg = None
statlist = None
sf = None


def getStationDict(file, slspecs, *args):
    """
    Read coordinates from station file
    :param file: station file
    :param slspecs: specifications of current slice
    :param args: size, marker, color
    :return: dictionary as specified in RegularSphericalGrid
    """
    from askipy.readEventStationFile import stationList
    global statlist

    if not statlist: statlist = stationList(file, list_type='standard')
    nstat = statlist.nstat
    slonlat = np.array([[float(statlist.stations[i]['lon']) for i in range(0, nstat)],
                        [float(statlist.stations[i]['lat']) for i in range(0, nstat)]])
    return {'coor': slonlat, 'size': args[0], 'marker': args[1], 'color': args[2]}


def getFaultDict(file, slspecs, *args):
    """
    Read shapefile, extract fault lines and fill a dictionary
    used as kwargs to plotHslice.
    :param file: File with fault line specifications
    :param slspecs: specifications of current slice
    :param args:
    :return: dictionary used as kwargs for plotHslice
    """
    import shapefile as shp
    global sf

    if not sf: sf = shp.Reader(file)
    ls = {1: 'solid', 2: 'solid', 3: 'dashed'}
    lw = {1: 1.6, 2: 0.8, 3: 1.6}
    color = 'black'
    #
    #  faultdict is a dictionary containing as key the order index (oid) with value a dictionary containing:
    #  (lon,lat)-point coordinates as Nx2-array, lineweight and linestyle and color
    #
    faultdict = dict()
    for el in sf.shapeRecords():
        ft = el.record['fault_type']
        faultdict[el.shape.oid] = {'coor': np.array(el.shape.points), 'lw': lw[ft], 'ls': ls[ft], 'color': color}
    return faultdict


def getColsumDict(file, slspecs, *args):
    """
    :param file: name of file from which 3D column sum data are read
    :param slspecs: specifications of current slice
    :param args: 6 arguments containg values for csbad, ncol, radius, nlat, nlon, slname
    """
    from matplotlib.colors import ListedColormap
    from askipy.regularSphericalGrid import RegularSphericalGrid
    from scipy.special import erf
    global csrsg

    if not csrsg: csrsg = RegularSphericalGrid.readFromHDF(file)
    csbad = args[0]; ncol = args[1]
    # case horizontal slice
    if 'radius' in slspecs:
        radius = slspecs['radius']; nlat = slspecs['nlat']; nlon = slspecs['nlon']; slname = slspecs['name']
        hslice = csrsg.extractHorizontalSlice(radius, nlat, nlon, slname)
        csval, vmin, vmax = hslice.getReshapedField('colsum', depth_normalize=False, normalize=False,
                                                     cliprel=-1.0, clipabs=-1.0, sym=False, percent=False,
                                                     extrema=[0,1])
    # case vertical slice
    elif 'gc' in slspecs:
        nr = slspecs['nr']; ndel = slspecs['ndel']; name = slspecs['name']
        vslice, gc = csrsg.extractVerticalSlice(*slspecs['gc'], nr, ndel, name)
        csval, vmin, vmax = vslice.getReshapedField('colsum', depth_normalize=False, normalize=False,
                                                     cliprel=-1.0, clipabs=-1.0, sym=False, percent=False,
                                                     extrema=[0,1])
    else:
        csval = np.zeros((1,1))

    # construct colormap
    vals = np.ones((ncol,4))
    kw = int(csbad*ncol)
    cs = np.linspace(0.0, csbad, kw)
    vals[0:kw,3] = erf(2.0*(csbad-cs)/csbad)/erf(2.0)
    vals[kw:,3] = 0.0
    cmap = ListedColormap(vals)

    return {'array': csval, 'colormap': cmap, 'vmin': 0.0, 'vmax': 1.0, 'shading': 'gouraud'}
