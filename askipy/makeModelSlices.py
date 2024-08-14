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
#   Wrapper script for making diverse model slices
#
import sys
import logging
import subprocess
from os import path as os_path
from os import system as os_system
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.plotRsgSlices import plotRsgSlices
from askipy.basejob import Basejob

logger = logging.getLogger(__name__)


def makeModelSlices(clargs):
    ap = ArgumentParser(description="Make slices through inverted model",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--it", help="Take model from this iteration", default='current')
    ap.add_argument("--itr", help="Take reference model from this iteration", default='0')
    ap.add_argument("--modelfile", help="Name of model file in output_files", default='total_absdev_model.kim')
    ap.add_argument("--rsgpath", help="Path to rsg-folder under output_files with slash", default='kimrsg/')
    ap.add_argument("--msldef", help="Text file defining model slices", default='kim_slicedef.txt')
    ap.add_argument("--csldef", help="Text file defining column sum slices", default='colsum_slicedef.txt')
    ap.add_argument("--csitlist", help="List of iteartion indices for column sum collection", default='Auto')
    ap.add_argument("--csbad", help="Lower bound of normalized column sums siginifying lacking resolution", default='0.05')
    ap.add_argument("--joblog", help="Logging of execution process", default='makeModelSlices.log')
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile
    modelfile = args.modelfile
    rsgpath = args.rsgpath
    slicedef = args.msldef
    csldef = args.csldef
    csitlist = args.csitlist
    csbad = args.csbad
    joblog = args.joblog

    # configure logging
    logging.basicConfig(filename=joblog, filemode='w', level=logging.INFO,
                        format='{asctime:s}|{module:s}|{funcName:s}|{message:s}',
                        datefmt='%Y-%m-%d %H:%M:%S', style='{', force=False)

    #  create job objects
    basejob = Basejob(queue='low', hostname='minos15', memory='2G')

    itref = int(args.itr)
    if args.it == 'current':
        fmp = fromMainParfile(main_parfile, itref=itref, checkdir=False)
        it = fmp['current_iteration_step']
    else:
        it = int(args.it)
        fmp = fromMainParfile(main_parfile, it=it, itref=itref, checkdir=False)

    #  run throwKimOnSphericalPseudoMesh for 3D regular grid
    program = 'throwKimOnSphericalPseudoMesh'
    jobargs = " ".join([main_parfile, '-it', str(it), '-itref', str(itref), '-prop', 'vp',
                        '-modelfile', modelfile, '-rsgpath', rsgpath])
    cp = basejob.submit(fmp['path_aski_code']+'bin/', program, jobargs, joblog)
    cp.check_returncode()

    # run collectColsums for 3D regular grid
    program = 'collectColsums'
    if csitlist == 'Auto': csitlist = '"' + ' '.join([str(i) for i in range(1,it+1)]) + '"'
    jobargs = " ".join([main_parfile, '-itlist', csitlist, '-rsgpath', rsgpath, '-prop', 'vp'])
    cp = basejob.submit(fmp['path_aski_code']+'bin/', program, jobargs, joblog)
    cp.check_returncode()

    #  run computeRegularSphericalGridSlices for model om slices
    program = 'computeRegularSphericalGridSlices'
    filename = "pseudo_mesh_vp_it_{:02d}_itr_{:02d}.hdf".format(it,itref)
    jobargs = " ".join([fmp['path_output_files'] + rsgpath + filename, '-slicedefs', slicedef])
    cp = basejob.submit(fmp['path_aski_code']+'bin/', program, jobargs, joblog)
    cp.check_returncode()

    #  run computeRegularSphericalGridSlices for collected column sums on slices
    program = 'computeRegularSphericalGridSlices'
    filename = "pseudo_mesh_vp_colsum_allf.hdf"
    jobargs = " ".join([fmp['path_output_files'] + rsgpath + filename, '-slicedefs', slicedef])
    cp = basejob.submit(fmp['path_aski_code']+'bin/', program, jobargs, joblog)
    cp.check_returncode()

    #  plot model slices (and overlay with collected column sums)
    slicebase = "pseudo_mesh_vp_it_{:02d}_itr_{:02d}".format(it,itref)
    colsumbase = "pseudo_mesh_vp_colsum_allf"
    plotRsgSlices([fmp['path_output_files']+rsgpath+slicebase,
                                  '--slicedefs', slicedef,
                                  '--stationfile', fmp['file_station_list']+ '_S',
                                  '--faultfile', 'tectonic_maps_4dmb_2020_09_17/shape_files/faults_alcapadi',
                                  '--colsumbase', fmp['path_output_files']+rsgpath+colsumbase,
                                  '--ncol', '101', '--csbad', csbad])

    #  plot collected column sum slices
    #  slice definitions should contain the same profiles, otherwise great circles do not fit
    slicebase = "pseudo_mesh_vp_colsum_allf"
    plotRsgSlices([fmp['path_output_files']+rsgpath+slicebase,
                                  '--slicedefs', csldef,
                                  '--stationfile', fmp['file_station_list']+ '_S',
                                  '--faultfile', 'None',
                                  '--colsumbase', 'None'])

    # convert model from spm to rsg
    clat = str(fmp['invgrid_center_lat'])
    clon = str(fmp['invgrid_center_lon'])
    program = 'evaluateSpmOnRsg'
    spmfile = "pseudo_mesh_vp_it_{:02d}_itr_{:02d}.hdf".format(it,itref)
    rsgfile = "vgrid_vp_it_{:02d}_itr_{:02d}.hdf".format(it,itref)
    jobargs = " ".join([fmp['path_output_files'] + rsgpath + spmfile, fmp['path_output_files'] + rsgpath + rsgfile,
               '-clat', clat, '-clon', clon, '-hwlat', '6.0', '-hwlon', '11.0', '-rmax', '6371.0', '-rmin', '5771.0',
               '-dlat', '0.1', '-dlon', '0.1', '-dr', '10.0'])
    cp = basejob.submit(fmp['path_aski_code']+'bin/', program, jobargs, joblog)
    cp.check_returncode()


#  allow module to be run as a script
if __name__ == "__main__":
    makeModelSlices(sys.argv[1:])
