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
#   Wrapper script for making diverse kernel slices
#
import sys
import subprocess
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.fromIterParfile import fromIterParfile
from askipy.plotRsgSlices import plotRsgSlices


def makeKernelSlices(clargs):
    ap = ArgumentParser(description="Make slices through selected kernels",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("event", help="Name of event (e.g. 20180114_091852.a)")
    ap.add_argument("nsname", help="Net-station name (e.g. GR.MOX)")
    ap.add_argument("--ifreq", help="kernel for given frequency index", default='12')
    ap.add_argument("--it", help="Take model from this iteration", default='current')
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile
    if args.it == 'current':
        fmp = fromMainParfile(main_parfile, checkdir=False)
        it = fmp['current_iteration_step']
    else:
        it = int(args.it)
        fmp = fromMainParfile(main_parfile, it=it, checkdir=False)

    #  run throwKernelOnSphericalPseudoMesh
    execpath = fmp['path_aski_code'] + 'bin/'
    executable = 'throwKernelOnSphericalPseudoMesh'
    jobargs = [main_parfile, args.event, args.nsname, '-ifreq', args.ifreq, '-it', str(it), '-prop', 'vp', '-comp', 'UP']
    print([execpath+executable] + jobargs)
    subprocess.check_call([execpath+executable] + jobargs)

    # set folder and file base names
    filebase = "pseudo_mesh_vp_f{:03d}_it{:03d}".format(int(args.ifreq), it)
    folder = fmp['path_output_files'] + 'kernelrsg_' + args.event.strip('.a') + '_' + args.nsname + '/'

    #  run computeRegularSphericalGridSlices
    execpath = fmp['path_aski_code'] + 'bin/'
    executable = 'computeRegularSphericalGridSlices'
    filename = filebase + '.hdf'
    jobargs = [folder + filename, '-slicedefs', 'kernel_slicedef.txt']
    print([execpath+executable] + jobargs)
    subprocess.check_call([execpath+executable] + jobargs)

    #  plot slices
    plotRsgSlices([folder+filebase,
                                  '--slicedefs', 'kernel_slicedef.txt',
                                  '--stationfile', fmp['file_station_list']+ '_S'])



#  allow module to be run as a script
if __name__ == "__main__":
    makeKernelSlices(sys.argv[1:])
