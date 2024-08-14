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
#   Run ASKI iteration workflow
#
import sys
import logging
from glob import glob
from shutil import copy
from os import path as os_path
from os import system as os_system
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.fromIterParfile import fromIterParfile
from askipy.prepareSpecfemKernel import prepareSpecfemKernel
from askipy.createDmspaceFile import createDmspaceFile
from askipy.gpujob import Gpujob
from askipy.mpijob import Mpijob

logger = logging.getLogger(__name__)

def doAskiIteration(clargs):
    ap = ArgumentParser(description="Perform one ASKI iteration starting with prepareSpecfemKernel --disp",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--newfreq", help="invert with new frequency set", action='store_true')
    ap.add_argument("--gpugt", help="number of gpus for green tensor", default='4')
    ap.add_argument("--nogt", help="Stop workflow before Green tensor calculations", action='store_true')
    ap.add_argument("--forward", help="Stop workflow after transformSynthetics", action='store_true')
    ap.add_argument("--joblog", help="Logging of execution process (relative to output_files", default='doAskiIteration.log')
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile
    newfreq = args.newfreq
    gpugt = args.gpugt
    nogt = args.nogt
    forward = args.forward
    joblog = args.joblog

    #  get main and iteration paramters
    fmp = fromMainParfile(main_parfile)
    fip = fromIterParfile(fmp['parfile_iteration_step'])
    it =  fmp['current_iteration_step']

    #  configure logging
    logfile = fmp['path_output_files'] + joblog
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO,
                        format='{asctime:s}|{module:s}|{funcName:s}|{message:s}',
                        datefmt='%Y-%m-%d %H:%M:%S', style='{', force=False)

    #  create job objects
    jobmpi = Mpijob(queue='low', hostname='minos15', memory='2G')
    jobgpu = Gpujob(queue='low', hostname='minos15', memory='2G')

    # ---------------------------- SPECFEM DISP ---------------------------------------------
    #
    # run prepareSpecfemKernel.py with --disp
    prepareSpecfemKernel([main_parfile, '--disp'])

    # run xgenerate_databases
    program = 'xgenerate_databases'
    logfile = fmp['path_specfem_logs'] + 'kdisp/' + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, nproc=fip['nproc'])
    cp.check_returncode()

    # run xspecfem3D
    program = 'xspecfem3D'
    logfile = fmp['path_specfem_logs'] + 'kdisp/' + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobgpu.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, ngpu=2, nproc=fip['nproc'])
    cp.check_returncode()
    for f in glob('OUTPUT_FILES/*.txt'):
        copy(f, fmp['path_specfem_logs'] + 'kdisp/')
    logger.info("SPECFEM_DISP output files copied to %s.", fmp['path_specfem_logs'] + 'kdisp/')

    # transform new Specfem synthetics
    program = 'transformSpecfemSynthetics'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=1)
    cp.check_returncode()

    # exit here if forward is set
    if forward: exit()

    # ---------------------------- DMSPACE  ---------------------------------------------
    #
    # it it > 1, copy ASKI_dmspace_masked file from previous iteration
    # to current one with name ASKI_dmspace
    # else create a new one
    if it > 1:
        if newfreq:
            createDmspaceFile([main_parfile])
        else:
            copy(fmp['path_dmspace'].replace(fmp['iteration_step_path'], fmp['prev_iteration_step_path']) + 'ASKI_dmspace_masked',
                 fmp['path_dmspace'] + 'ASKI_dmspace')
            logger.info("Copied masked dmspace file from previous to current iteration")
    else:
        createDmspaceFile([main_parfile, '--it', '1'])

    # ---------------------------- MISFITS ---------------------------------------------
    #
    # compute misfits unmasked
    program = 'computeMisfits'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=1)
    cp.check_returncode()

    # compute misfits masked
    logfile = fmp['path_output_files'] + program + 'Masked.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-masked']), logfile, nproc=1)
    cp.check_returncode()

    # produce dmspace file for all frequency sets for current iteration
    ifreq_set = ['6-8-10-12', '9-11-13-15', '12-14-16-18']
    for ifreq in ifreq_set:
        createDmspaceFile([main_parfile, '--ifreq', ifreq, '--it', str(it)])
    createDmspaceFile([main_parfile, '--ifreq', 'allf', '--it', str(it)])

    # compute misfits for all frequency sets, normal and masked
    program = 'computeMisfits'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    logmaskfile = fmp['path_output_files'] + program + 'Masked.log'
    if os_path.exists(logmaskfile): os_system('rm -f ' + logmaskfile)
    for ifreq in ifreq_set:
        ext = ifreq
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-ext', ext]), logfile, nproc=1)
        cp.check_returncode()
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-ext', ext, '-masked']), logmaskfile, nproc=1)
        cp.check_returncode()

    # compute misfits for all frequencies
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-ext', 'allf']), logfile, nproc=1)
    cp.check_returncode()
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-ext', 'allf', '-masked']), logmaskfile, nproc=1)
    cp.check_returncode()

    # run decompose dmspace masked
    program = 'decomposeDmspace'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, ' '.join([main_parfile, '-masked']), logfile, nproc=1)
    cp.check_returncode()

    # exit here if nogt is set
    if nogt: exit()

    # ---------------------------- SPECFEM GT  ---------------------------------------------
    #
    # run prepareSpecfemKernel.py with --gt
    prepareSpecfemKernel([main_parfile, '--gt'])

    # run xspecfem3D gt
    program = 'xspecfem3D'
    logfile = fmp['path_specfem_logs'] + 'kgt/' + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobgpu.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, ngpu=int(gpugt), nproc=fip['nproc'])
    cp.check_returncode()
    for f in glob('OUTPUT_FILES/*.txt'):
        copy(f, fmp['path_specfem_logs'] + 'kgt/')
    logger.info("SPECFEM_DISP output files copied to %s.", fmp['path_specfem_logs'] + 'kgt/')

    # ---------------------------- KERNELS, SOLVER, UPDATE, NORM  ---------------------------------------------
    #
    # run kernels
    program = 'computeKernelsDmspaceParallel'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=fip['nproc'], memory='6G')
    cp.check_returncode()

    # run solver
    program = 'solveCglsKernelSystem'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=fip['nproc'], memory='8G')
    cp.check_returncode()

    # run update
    program = 'updateModelPerturbations'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=1)
    cp.check_returncode()

    # run model norm
    program = 'computeModelNorm'
    logfile = fmp['path_output_files'] + program + '.log'
    if os_path.exists(logfile): os_system('rm -f ' + logfile)
    cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=1, memory='2G')
    cp.check_returncode()

    #  run throwKimOnSphericalPseudoMesh for 3D model
    program = 'throwKimOnSphericalPseudoMesh'
    logfile = fmp['path_output_files'] + program + '.log'
    jobargs = " ".join([main_parfile, '-it', str(it), '-itref', '0', '-prop', 'vp'])
    cp = jobmpi.submit(fmp['path_aski_code']+'bin/', program, jobargs, logfile, nproc=1, memory='2G')
    cp.check_returncode()

    # run collectColsums for model
    program = 'collectColsums'
    logfile = fmp['path_output_files'] + program + '.log'
    csitlist = '"' + ' '.join([str(i) for i in range(1,it+1)]) + '"'
    jobargs = " ".join([main_parfile, '-itlist', csitlist, '-prop', 'vp'])
    cp = jobmpi.submit(fmp['path_aski_code']+'bin/', program, jobargs, logfile, nproc=1, memory='2G')
    cp.check_returncode()

    # convert model from spm to rsg
    clat = str(fmp['invgrid_center_lat'])
    clon = str(fmp['invgrid_center_lon'])
    program = 'evaluateSpmOnRsg'
    logfile = fmp['path_output_files'] + program + '.log'
    spmfile = "pseudo_mesh_vp_it_{:02d}_itr_{:02d}.hdf".format(it, 0)
    rsgfile = "vgrid_vp_it_{:02d}_itr_{:02d}.hdf".format(it, 0)
    jobargs = " ".join([fmp['path_output_files'] + 'kimrsg/' + spmfile, fmp['path_output_files'] + 'kimrsg/' + rsgfile,
               '-clat', clat, '-clon', clon, '-hwlat', '6.0', '-hwlon', '11.0', '-rmax', '6371.0', '-rmin', '5771.0',
               '-dlat', '0.1', '-dlon', '0.1', '-dr', '10.0'])
    cp = jobmpi.submit(fmp['path_aski_code']+'bin/', program, jobargs, logfile, nproc=1, memory='2G')
    cp.check_returncode()


#  allow module to be run as a script
if __name__ == "__main__":
    doAskiIteration(sys.argv[1:])


