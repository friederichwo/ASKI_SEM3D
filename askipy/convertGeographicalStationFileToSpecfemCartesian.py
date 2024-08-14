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
#  Tool to convert a station file with geographical coordinates
#  to Cartesian coordinates used by SPECFEM
# ---------------------------------------------------------------
import sys
import numpy as np
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.readEventStationFile import stationList
import askipy.axesRotation as ar
from askipy.fromMainParfile import fromMainParfile


def convertGeographicalStationFileToSpecfemCartesian(clargs):
    """
    :param clargs: list of command line arguments
    """
    #  read command line arguments
    #
    ap = ArgumentParser(description="Convert station file from geographical to SPECFEM Cartesian coordinates",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("stationfile", help="Name of station file with geographical coordinates")
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--ignore_alt", help="If set, station altitude is set to zero", action='store_true')
    args = ap.parse_args(clargs)

    #  get parameters from ASKI main parameter file
    fmp = fromMainParfile(args.main_parfile, checkfile=False)
    clat = fmp['invgrid_center_lat'] * np.pi/180.
    clon = fmp['invgrid_center_lon'] * np.pi/180.
    rearth = fmp['REARTH']
    outfile = fmp['file_station_list']

    #  read geographical station file
    sl = stationList(args.stationfile, "standard")

    #  open output file
    fout = open(outfile, 'w')

    #  write C into first line
    #  walk through stations
    fout.write('C\n')
    for js in range(sl.nstat):
        statinf = sl.stations[js]
        staname = statinf['staname']
        netcode = statinf['netcode']
        latrad = float(statinf['lat']) * np.pi/180.
        lonrad = float(statinf['lon']) * np.pi/180.
        alt = float(statinf['alt'])
        if args.ignore_alt:
            alt = 0.0
        xg, yg, zg = ar.coordinatesLCfromLS(rearth+alt, 0.5*np.pi-latrad, lonrad)
        xlc, ylc, zlc = ar.coordinatesLCfromGC(0.5*np.pi-clat, clon, xg, yg, zg)
        x, y, z = ar.coordinatesRCfromLC(0.5*np.pi, xlc, ylc, zlc)
        txt = "{name:<6}{net:<3}{x:15.3f}{y:15.3f}{z:15.3f}\n"
        fout.write(txt.format(name=staname, net=netcode, x=x, y=y, z=z-rearth))
    fout.close()

    sys.stdout.write('Geographical station information taken from: ' + args.stationfile + '\n')
    sys.stdout.write('Specfem Cartesian station information written to: ' + outfile + '\n')

#  allow module to be run as a script
if __name__ == "__main__":
    convertGeographicalStationFileToSpecfemCartesian(sys.argv[1:])