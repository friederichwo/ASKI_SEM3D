#!/usr/bin/env python3
#
# ----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#   Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
# ----------------------------------------------------------------------------
#
import sys
from os import path as os_path
from os import mkdir as os_mkdir
from os import symlink as os_symlink
from os import environ as os_environ
from shutil import copy
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.helperFunctions import check_directory, check_free, get_time_string


def createIterationFolders(clargs):
    """
    Create all subfolders of main inversion folder for one iteration.
    Copy templates of parameter files to their proper places
    and fill in values taken from the main parameter file.
    Requires that main inversion folder exists and main parameter file
    contains reasonable values and was copied to main inversion folder.
    If citer > 1, only deal with iteration subdirectories and iter_parfile.
    Run this script from the main inversion folder.
    """
    
    # command line
    ap = ArgumentParser(description="Create inversion subfolders and edit parameter files",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of main ASKI parameter file")
    ap.add_argument("--template_dir", help="Name of ASKI template folder reltive to home directory",
                    default='work_git/ASKI_SEM3D_CAR/templates/')
    args = ap.parse_args(clargs)
    template_dir = os_environ['HOME'] + '/' + args.template_dir

    # start logging
    sys.stdout.write("createIterationFolders started at " + get_time_string() + "\n")
    sys.stdout.write("Arguments: " + str(args) + "\n")
    sys.stdout.write("Create inversion subfolders and copy parameter files" + "\n")
    check_directory(template_dir)

    # get ASKI main parameters
    fmp = fromMainParfile(args.main_parfile, checkdir=False, checkfile=False)
    
    #  create all inversion subfolders only for first iteration
    if fmp['current_iteration_step'] == 1:
        # DATA folder + meshfem3d_files
        if not os_path.exists(fmp['path_specfem_input']):
            os_mkdir(fmp['path_specfem_input'])
            sys.stdout.write("Make directory " + fmp['path_specfem_input'] + "\n")
        meshfem = os_path.join(fmp['path_specfem_input'], 'meshfem3D_files')
        if not os_path.exists(meshfem):
            os_mkdir(meshfem)
            sys.stdout.write("Make directory " + meshfem + "\n")

        # DATABSES_MPI folder
        if not os_path.exists(fmp['path_specfem_databases']):
            os_mkdir(fmp['path_specfem_databases'])
            sys.stdout.write("Make directory " + fmp['path_specfem_databases'] + "\n")

        # PHASE_END_TIMES folder
        path_pet = fmp['bigssd_path_inversion'] + fmp['path_phase_end_times']
        if not os_path.exists(path_pet):
            os_mkdir(path_pet)
            sys.stdout.write("Make directory " + path_pet + "\n")
        if fmp['bigssd_path_inversion'] != fmp['main_path_inversion']:
            link = path_pet.replace(fmp['bigssd_path_inversion'],fmp['main_path_inversion']).rstrip('/')
            if not os_path.exists(link):
                os_symlink(path_pet, link)

        # MESH folder
        specfem_mesh = 'MESH/'
        if not os_path.exists(specfem_mesh):
            os_mkdir(specfem_mesh)
            sys.stdout.write("Make directory " + specfem_mesh + "\n")

        # OUTPUT_FILES folder
        specfem_of = 'OUTPUT_FILES/'
        if not os_path.exists(specfem_of):
            os_mkdir(specfem_of)
            sys.stdout.write("Make directory " + specfem_of + "\n")

        # copy Par_file and Par_file_ASKI and CMTSOLUTION and FORCESOLUTION from template to DATA
        copy(os_path.join(template_dir, 'Par_file'), fmp['path_specfem_input'])
        copy(os_path.join(template_dir, 'Par_file_ASKI'), fmp['path_specfem_input'])
        copy(os_path.join(template_dir, 'CMTSOLUTION'), fmp['path_specfem_input'])
        copy(os_path.join(template_dir, 'FORCESOLUTION'), fmp['path_specfem_input'])
        sys.stdout.write("Copied Par_file, Par_file_ASKI, CMTSOLUTION and FORCESOLUTION to  " + fmp['path_specfem_input'] + "\n")

        # copy Mesh_Par_file, interfaces.dat and interfaces_01_zcoor.dat to DATA/meshfem3d_files
        copy(os_path.join(template_dir, 'Mesh_Par_file'), meshfem)
        copy(os_path.join(template_dir, 'interfaces.dat'), meshfem)
        copy(os_path.join(template_dir, 'interface_01_zcoor.dat'), meshfem)
        sys.stdout.write("Copied Mesh_Par_file and interfaces files to  " + meshfem + "\n")

        # Gemini subfolder
        gemini = fmp['path_gemini']
        if not os_path.exists(gemini):
            os_mkdir(gemini)
            sys.stdout.write("Make directory " + gemini + "\n")

        # copy parfile_gemini to GEMINI subpath
        copy(os_path.join(template_dir, 'parfile_gemini'), gemini)
        sys.stdout.write("Copied parfile_gemini to " + gemini + "\n")

        # path to measured data
        if not os_path.exists(fmp['path_measured_data']):
            os_mkdir(fmp['path_measured_data'])
            sys.stdout.write("Make directory " + fmp['path_measured_data'] + "\n")

        # path to measured seismograms
        if not os_path.exists(fmp['path_measured_seis']):
            os_mkdir(fmp['path_measured_seis'])
            sys.stdout.write("Make directory " + fmp['path_measured_seis'] + "\n")

        # path to injection seismograms
        # put on BIGSSD with symlink if available
        # if not available assume BIGSSD == MAIN_PATH_INVERSION
        path_injection_seis = fmp['bigssd_path_inversion'] + fmp['path_injection_seis']
        if not os_path.exists(path_injection_seis):
            os_mkdir(path_injection_seis)
            sys.stdout.write("Make directory " + path_injection_seis + "\n")
        if fmp['bigssd_path_inversion'] != fmp['main_path_inversion']:
            link = path_injection_seis.replace(fmp['bigssd_path_inversion'], fmp['main_path_inversion']).rstrip('/')
            if not os_path.exists(link):
                os_symlink(path_injection_seis, link)

        # copy property correlation file to inversion folder
        copy(os_path.join(template_dir, fmp['property_correlation_file']), './')
        sys.stdout.write("Copied property correlation to  " + './' + "\n")

    #  create iteration folder
    if not os_path.exists(fmp['iteration_step_path']):
        os_mkdir(fmp['iteration_step_path'])
        sys.stdout.write("Make directory " + fmp['iteration_step_path'] + "\n")
    citer = fmp['current_iteration_step']

    #  create iteration folder on bigssd if available (kdisp, kgt and sensitivities)
    bigssd_iter_path = fmp['bigssd_path_inversion'] +  fmp['iteration_step_path']
    if fmp['bigssd_path_inversion'] != fmp['main_path_inversion']:
        if not os_path.exists(bigssd_iter_path):
            os_mkdir(bigssd_iter_path)
            sys.stdout.write("Make directory " + bigssd_iter_path + "\n")

    #  copy iter_parfile template to iteration folder
    if citer == 1:
        copy(template_dir + 'iter_parfile', fmp['iteration_step_path'])
        sys.stdout.write("Copied iter_parfile from templates to " + fmp['iteration_step_path'] + "\n")
    else:
        copy(fmp['parfile_iteration_step'].replace(fmp['iteration_step_path'],fmp['prev_iteration_step_path']),
             fmp['iteration_step_path'])
        sys.stdout.write("Copied iter_parfile from previous iteration to " + fmp['iteration_step_path'] + "\n")

    #  create iteration subfolders in main path inversion
    for key in ['path_aski_main_files', 'path_synthetic_data', 'path_output_files', 'path_dmspace',
                'path_vtk_files', 'path_specfem_logs']:
        if not os_path.exists(fmp[key]):
            os_mkdir(fmp[key])
            sys.stdout.write("Make directory " + fmp[key] + "\n")

    if not os_path.exists(fmp['path_specfem_logs'] + 'kdisp/'): os_mkdir(fmp['path_specfem_logs'] + 'kdisp/')
    if not os_path.exists(fmp['path_specfem_logs'] + 'kgt/'):   os_mkdir(fmp['path_specfem_logs'] + 'kgt/')

    #  subfolders of bigssd iteration folder if available
    for key in ['path_kernel_displacements', 'path_kernel_green_tensors', 'path_sensitivity_kernels']:
        subfolder = fmp['bigssd_path_inversion'] + fmp[key]
        if not os_path.exists(subfolder):
            os_mkdir(subfolder)
            sys.stdout.write("Make directory " + subfolder + "\n")
        if fmp['bigssd_path_inversion'] != fmp['main_path_inversion']:
            if not os_path.exists(fmp['main_path_inversion'] + fmp[key].strip('/')):
                os_symlink(fmp['bigssd_path_inversion'] + fmp[key], fmp['main_path_inversion'] + fmp[key].rstrip('/'))
                sys.stdout.write("Make symlink from " + fmp['bigssd_path_inversion'] + fmp[key] + ' to' +
                                 fmp['main_path_inversion'] + fmp[key].strip('/') + "\n")
#
    sys.stdout.close()


#  allow module to be run as a script
if __name__ == "__main__":
    createIterationFolders(sys.argv[1:])


