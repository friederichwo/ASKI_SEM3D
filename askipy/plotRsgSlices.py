#!/usr/bin/env python3
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


import sys
import logging
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
import json
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.regularSphericalGrid import RegularSphericalGrid
from askipy import externalDecorations

logger = logging.getLogger(__name__)


def addTextDecoration(textspecs, slidx, name):
    """
    Add text decorations to a slice.
    A dictionary as specified in RegularSphericalGrid for adding text
    in slices must be provided in the configuration file
    :param textspecs: dictionary with text specifications taken from config file
    :param slidx: string index of current slice
    :param name: name of slice
    :return: a text dictionary as specified in RegularSphericalGrid
    """
    jdc = 1                                             # text decoration counter
    tdict = {}
    for idx in textspecs:                               # go through text items
        selslices = textspecs[idx]['slices']            # text should be added to selected slices
        if slidx in selslices or selslices == 'all':    # check if current slice slidx is selected
            tdict[jdc] = textspecs[idx]                 # add text to text plotting subdict
            jdc = jdc+1
            logger.info("Text decoration %s added to slice %s", textspecs[idx]["text"], name)
    return tdict


def addExternalDecoration(decospecs, slidx, slspecs):
    """
    Add externally specified decoration to a slice (points, lines, overlays).
    A function and its arguments specfied in the config is called that
    returns a dictionary as required by RegularSphericalGrid for plotting
    external decorations in vertical or horizontal slices
    :param decospecs: dictionary in config with keys <file>, <function>, <args> and <slices>
    :param slidx: str index of current slice
    :param slspecs: dict with specifications of current slice
    :return: a dictionary as required by RegularSphericalGrid for the corresponding decoration
    """
    jdc = 1                                                   # decoration counter
    decodict = {}
    for idx in decospecs:                                      # go through decorations
        selslices = decospecs[idx]['slices']                   # decoration should be added to selected slices
        if slidx in selslices or selslices == 'all':           # check if current slice is selected
            funcname = decospecs[idx]['function']              # name of function to be called
            args = decospecs[idx]['args']                      # list of extra arguments from config
            file = decospecs[idx]['file']                      # file providing overlay data
            # call function
            decodict[jdc] = getattr(externalDecorations, funcname)(file, slspecs, *args)
            jdc = jdc+1
            logger.info("External decoration using function %s added to slice %s", funcname, slspecs["name"])
    return decodict


def plotRsgSlices(clargs):
    """
    A generic routine to organise plotting of any number of horizontal and vertical slices
    through a regular spherical grid <rsgfile> which may optionally be decorated by text, lines,
    points and overlaying slices. All plot options supported by RegularSphericalGrid
    can be set through a configuration dictionary <sliceconf> read from a JSON file.

    This function can be run standalone or may be called from a main script written
    to perform a specific plotting task. The main script passes a regular spherical
    grid file and a configuration file to this function. In addition it passes a basename
    for the output png files.

    To add the dedorations, this function assumes that the configuration sepcifies a file
    containing data and a function that reads the data and returns a dictionary as required
    by RegularSphericalGrid for plotting. These functions should be put into a module called
    <externalDecorations> which is imported by this function and called using the <getattr> utility.
    """
    ap = ArgumentParser(description="Plot slices through a3D regular spherical grid",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("rsgfile", help="Name of HDF file with regular spherical grid data")
    ap.add_argument("slconffile", help="Name of configuration file")
    ap.add_argument("outbase", help="base name including directory for output files")
    ap.add_argument("--joblog", help="Logging of execution process", default='plotRsgSlices.log')
    args = ap.parse_args(clargs)
    outbase = args.outbase

    # configure logging
    logging.basicConfig(filename=args.joblog, filemode='w', level=logging.INFO,
                        format='{asctime:s}|{module:s}|{funcName:s}:    {message:s}',
                        datefmt='%Y-%m-%d %H:%M:%S', style='{', force=False)

    # read out 3D rsg file and configuration
    rsg = RegularSphericalGrid.readFromHDF(args.rsgfile)
    logger.info("3D grid read from %s ", args.rsgfile)
    with open(args.slconffile, 'r') as fp:
        sliceconf = json.load(fp)
    logger.info("Configuration read from %s ", args.slconffile)
    # ----------------------------------------------
    #  plot vertical slices into separate figures
    # ----------------------------------------------
    slconf = sliceconf['vslices']
    gcdata = []
    if len(slconf) > 0:
        # initialize plotting dictionary with common plotting options
        pldict = slconf['plot_options']

        # go through vertical slices specified in the configuration
        for i in range(slconf['number']):
            slidx = str(i+1)                                     # str slice index used as key for slice
            slspecs = slconf[slidx]                              # dict with slice specifiations
            lat1, lon1, lat2, lon2 = slspecs['gc']               # start and end point of great circle

            # extract vertical slice object from 3D grid, also return great circle object
            nr = slspecs['nr']; ndel = slspecs['ndel']; name = slspecs['name']
            vslice, gc = rsg.extractVerticalSlice(lat1, lon1, lat2, lon2, nr, ndel, name)
            logger.info("Vslice %s along gc <%f %f %f %f> extracted", name, lat1, lon1, lat2, lon2)

            # append great circle to great circle list used later for horizontal sections
            gcdata.append(gc)
            pldict['gc'] = gc                                    # set gc key in plotting dictionary

            # set extrema in plotting dictionary from slice specs
            pldict['extrema'] = slspecs['extrema']

            # add text decorations to slice if desired
            textspecs  = slconf['text']
            if len(textspecs) > 0:
                pldict['text'] = addTextDecoration(textspecs, slidx, name)

            # add external decorations
            decorations = ("points", "lines", "overlays")
            for deco in decorations:
                decospecs = slconf[deco]
                if len(decospecs) > 0:
                    pldict[deco] = addExternalDecoration(decospecs, slidx, slspecs)

            # plot slice
            fig = vslice.plotVslice(outbase, name, slspecs['field'], **pldict)
            logger.info("Vslice %s for field %s plotted", name, slspecs['field'])
            plt.close(fig)
    # -------------------------------------------------------------------------------
    #  now deal with horizontal slices
    # ------------------------------------------------------------------------------
    slconf = sliceconf['hslices']
    if len(slconf) > 0:
        # set projection
        borders = cf.BORDERS
        coastline = cf.COASTLINE
        if abs(rsg.lat_pole-0.5*np.pi) < 1.e-4:
            clon = np.rad2deg(rsg.lonmin+0.5*rsg.dlon*(rsg.nlon-1))
            clat = np.rad2deg(rsg.latmin+0.5*rsg.dlat*(rsg.nlat-1))
            stdpar = slconf['plot_options']['standard_parallels']
            proj = ccrs.LambertConformal(central_longitude=clon, central_latitude=clat,
                                         standard_parallels=tuple(stdpar))
            projdata = ccrs.PlateCarree()
            logger.info("Projection for hslices set to LambertConformal")
        else:
            proj = ccrs.RotatedPole(pole_longitude=np.rad2deg(rsg.lon_pole),
                                    pole_latitude=np.rad2deg(rsg.lat_pole))
            projdata = proj
            logger.info("Projection for hslices set to RotatedPole")

        # initialize plotting dictionary with common plotting options
        pldict = slconf['plot_options']
        pldict['gcdata'] = gcdata

        # loop over horizontal slices
        for i in range(slconf['number']):
            slidx = str(i+1)
            slspecs = slconf[slidx]
            radius = slspecs['radius']; nlat = slspecs['nlat']; nlon = slspecs['nlon']; name = slspecs['name']
            hslice = rsg.extractHorizontalSlice(radius, nlat, nlon, name)
            logger.info("Hslice %s at %f extracted", name, radius)
 
            # set extrema from slice specs
            pldict['extrema'] = slspecs['extrema']

            # add text decorations to slice if desired
            textspecs  = slconf['text']
            if len(textspecs) > 0:
                pldict['text'] = addTextDecoration(textspecs, slidx, name)

            # add external decorations
            decorations = ("points", "lines", "overlays")
            for deco in decorations:
                decospecs = slconf[deco]
                if len(decospecs) > 0:
                    pldict[deco] = addExternalDecoration(decospecs, slidx, slspecs)

            # plot slice
            fig = hslice.plotHslice(outbase, name, slspecs['field'], borders, coastline, proj, projdata, **pldict)
            logger.info("Hslice %s for field %s plotted", name, slspecs['field'])
            plt.close(fig)


#  allow module to be run as a script
if __name__ == "__main__":
    plotRsgSlices(sys.argv[1:])
