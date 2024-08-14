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
#   Run some ASKI executables on SGE
#
import sys
import logging
from glob import glob
from shutil import copy
from os import path as os_path
from os import system as os_system
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.gpujob import Gpujob
from askipy.mpijob import Mpijob
from askipy.pyjob import Pyjob
from askipy.fromMainParfile import fromMainParfile
from askipy.fromIterParfile import fromIterParfile

logger = logging.getLogger(__name__)


def runAskiProgram(clargs):
    tasklist = "kernels solver specfem_databases specfem_disp specfem_gt transform_syn misfits decompose_masked update norm iterate".split()
    ap = ArgumentParser(description="Submit an ASKI or SPECFEM executable to SGE",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("task", help="Task to be carried out", choices=tasklist)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--gpu", help="Number of GPUS", default='0')
    ap.add_argument("--joblog", help="Logging of execution process", default='runAskiProgram.log')
    args = ap.parse_known_args(clargs)
    task = args[0].task
    main_parfile = args[0].main_parfile
    ngpu = int(args[0].gpu)
    jobargs = " " + " ".join(args[1])
    joblog = args[0].joblog

    # configure logging
    logging.basicConfig(filename=joblog, filemode='w', level=logging.INFO,
                        format='{asctime:s}|{module:s}|{funcName:s}|{message:s}',
                        datefmt='%Y-%m-%d %H:%M:%S', style='{', force=False)

    #  paramters from main parfile
    fmp = fromMainParfile(main_parfile)
    fip = fromIterParfile(fmp['parfile_iteration_step'])

    #  create job obejcts
    jobmpi = Mpijob(queue='low', hostname='minos15', memory='2G')
    jobgpu = Gpujob(queue='low', hostname='minos15', memory='2G')
    jobpy  = Pyjob(queue='low', hostname='minos15', memory='2G')

    if task == 'kernels':
        program = 'computeKernelsDmspaceParallel'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=fip['nproc'], memory='6G')
        cp.check_returncode()

    elif task == 'solver':
        program = 'solveCglsKernelSystem'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile, logfile, nproc=fip['nproc'], memory='8G')
        cp.check_returncode()

    elif task == 'specfem_databases':
        program = 'xgenerate_databases'
        logfile = fmp['path_specfem_logs'] + 'kdisp/' + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, nproc=fip['nproc'])
        cp.check_returncode()

    elif task == 'specfem_disp':
        program = 'xspecfem3D'
        logfile = fmp['path_specfem_logs'] + 'kdisp/' + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobgpu.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, ngpu=ngpu, nproc=fip['nproc'])
        cp.check_returncode()
        for f in glob('OUTPUT_FILES/*.txt'):
            copy(f, fmp['path_specfem_logs'] + 'kdisp/')
        logger.info("SPECFEM_DISP output files copied to {:s}.", fmp['path_specfem_logs'] + 'kdisp/')

    elif task == 'specfem_gt':
        program = 'xspecfem3D'
        logfile = fmp['path_specfem_logs'] + 'kgt/' + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobgpu.submit(fmp['path_specfem_code'] + 'bin/', program, '', logfile, ngpu=ngpu, nproc=fip['nproc'])
        cp.check_returncode()
        for f in glob('OUTPUT_FILES/*.txt'):
            copy(f, fmp['path_specfem_logs'] + 'kgt/')
        logger.info("SPECFEM_GT output files copied to {:s}.", fmp['path_specfem_logs'] + 'kgt/')

    elif task == 'transform_syn':
        program = 'transformSpecfemSynthetics'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile + jobargs, logfile, nproc=1)
        cp.check_returncode()

    elif task == 'misfits':
        program = 'computeMisfits'
        if jobargs.find('masked') > 0:
            logfile = fmp['path_output_files'] + program + 'Masked.log'
            if os_path.exists(logfile): os_system('rm -f ' + logfile)
        else:
            logfile = fmp['path_output_files'] + program + '.log'
            if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile + jobargs, logfile, nproc=1)
        cp.check_returncode()

    elif task == 'decompose_masked':
        program = 'decomposeDmspace'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile + ' -masked', logfile, nproc=1)
        cp.check_returncode()

    elif task == 'update':
        program = 'updateModelPerturbations'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile + jobargs, logfile, nproc=1)
        cp.check_returncode()

    elif task == 'norm':
        program = 'computeModelNorm'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobmpi.submit(fmp['path_aski_code'] + 'bin/', program, main_parfile + jobargs, logfile, nproc=1, memory='2G')
        cp.check_returncode()

    elif task == 'iterate':
        program = 'doAskiIteration.py'
        logfile = fmp['path_output_files'] + program + '.log'
        if os_path.exists(logfile): os_system('rm -f ' + logfile)
        cp = jobpy.submit(fmp['path_aski_code'] + 'askipy/', program, main_parfile + jobargs, logfile, memory='1G')
        cp.check_returncode()

    else:
        raise Exception('runAskiProgram: unknown task' + task)

#  allow module to be run as a script
if __name__ == "__main__":
    runAskiProgram(sys.argv[1:])

