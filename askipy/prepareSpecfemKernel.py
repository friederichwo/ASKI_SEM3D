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
from os import path as os_path
from os import mkdir as os_mkdir
from os import system as os_system
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.fromIterParfile import fromIterParfile
from askipy.helperFunctions import check_file, get_time_string
from askipy.prepareSpecfemHelpers import set_parfile, read_timing_injected_wavefields, get_sequential_force, \
    writeSpecfemStations


def prepareSpecfemKernel(clargs):
    ap = ArgumentParser(description="Prepare SPECFEM kernel displacement job for ASKI",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--disp", action="store_true", help="prepare kernel displacement computation")
    group.add_argument("--gt",  action="store_true", help="prepare kernel Green tensor computation")
    ap.add_argument("--no_gpu", action="store_true", help="SPECFEM does not use GPU")
    args = ap.parse_args(clargs)
    if args.no_gpu:
        specfem_uses_gpu = '.false.'
    else:
        specfem_uses_gpu = '.true.'

    print("Start log at " + get_time_string())
    print("Arguments: " + str(args))

    #  extract main ASKI parameters
    fmp = fromMainParfile(args.main_parfile)
    df = fmp['measured_data_frequency_step']
    print("Extracted settings from ASKI main parameter file")

    #  read and extract ASKI iteration parameters
    fip = fromIterParfile(fmp['parfile_iteration_step'])
    nf = fip['iteration_step_number_of_freq']
    ifreq = fip['iteration_step_index_of_freq']
    output_prev = fmp['path_output_files'].replace(fmp['iteration_step_path'], fmp['prev_iteration_step_path'])
    file_kim = output_prev + fip['file_aski_inverted_model']
    if fip['impose_aski_inverted_model'] == '.true.':
        check_file(file_kim)
    print("Extracted settings from ASKI iteration parameter file")

    #  create STATIONS file for SPECFEM
    #  if fmp['statlist'] contains x,y,z instead of lat,lon,alt, write STATIONS file in
    #  the order y,x,depth because in get_elevation they will be exchanged.
    if args.disp:
        writeSpecfemStations(os_path.join(fmp['path_specfem_input'], 'STATIONS'), fmp['statlist'])
    else:
        writeSpecfemStations(os_path.join(fmp['path_specfem_input'], 'STATIONS'), fmp['gtreclist'])
    print("STATIONS file for SPECFEM generated")

    #  get and check names of SPECFEM Par_file and SPECFEM ASKI parfile
    Par_file = fmp['path_specfem_input'] + 'Par_file'
    check_file(Par_file)
    Par_file_ASKI = fmp['path_specfem_input'] + 'Par_file_ASKI'
    check_file(Par_file_ASKI)

    #  set some parameters that depend on --disp or --gt options
    if args.disp:                             # kernel displacement
        use_force_point_source = '.false.'
        couple_with_injection = '.true.'
        sequential_sources_file = "kernel_disp_specfem_injection"
        aski_decon = '.false.'
        print_stf = '.false.'
        aski_apply_phase_taper = fip['aski_dft_apply_taper']
        taper_length = fip['aski_dft_taper_length_kd']
        tt_file = fmp['aski_ttime_table_file_kd']
        print("Specfic setting for kernel displacement: ")
        print("use force point source = " + use_force_point_source)
        print("couple with injection = " + couple_with_injection)
        print("sequential sources file = " + fmp['path_specfem_input'] + sequential_sources_file)
        print("ASKI deconvolution = " + aski_decon)
        print("print source time function = " + print_stf)
    else:                                     # kernel green tensor
        use_force_point_source = '.true.'
        couple_with_injection = '.false.'
        sequential_sources_file = "kernel_gt_specfem_forces"
        print_stf = '.true.'
        aski_decon = '.true.'
        aski_apply_phase_taper = fip['aski_dft_apply_taper']
        taper_length = fip['aski_dft_taper_length_gt']
        tt_file = fmp['aski_ttime_table_file_gt']
        print("Specfic setting for kernel Green tensor: ")
        print("use force point source = " + use_force_point_source)
        print("couple with injection = " + couple_with_injection)
        print("sequential sources file = " + fmp['path_specfem_input'] + sequential_sources_file)
        print("ASKI deconvolution = " + aski_decon)
        print("print source time function = " + print_stf)

    #  write Par_file
    set_parfile(Par_file, 
                [('USE_FORCE_POINT_SOURCE', use_force_point_source),
                 ('USE_SOURCES_RECEIVERS_Z', '.true.'),
                 ('PRINT_SOURCE_TIME_FUNCTION', print_stf),
                 ('USE_RICKER_TIME_FUNCTION', '.false.'), ('GPU_MODE', specfem_uses_gpu),
                 ('COUPLE_WITH_INJECTION_TECHNIQUE', couple_with_injection),
                 ('INJECTION_TECHNIQUE_TYPE', '4'), ('TRACTION_PATH', fmp['path_injection_seis']),
                 ('NSTEP', str(fip['nstep'])), ('DT', str(fip['dt'])),
                 ('SEQUENTIAL_SOURCES_DESCRIPTION_FILE', sequential_sources_file),
                 ('MULTIPLE_SEQUENTIAL_SOURCES', '.true.'),
                 ('MOVIE_SURFACE', fip['movie_surface']), ('MOVIE_VOLUME', fip['movie_volume'])])
    print("SPECFEM Par_file written: ")

    #  write Par_file_ASKI
    set_parfile(Par_file_ASKI, 
        [('USE_ASKI_BACKGROUND_MODEL', fmp['use_aski_background_model']), 
         ('FILE_ASKI_BACKGROUND_MODEL', fmp['file_aski_background_model']),
         ('IMPOSE_ASKI_INVERTED_MODEL', fip['impose_aski_inverted_model']),
         ('FILE_ASKI_INVERTED_MODEL', file_kim),
         ('COMPUTE_ASKI_OUTPUT', '.true.'),
         ('ASKI_USES_GPU', fmp['aski_uses_gpu']),
         ('ASKI_MAIN_FILE_ONLY', fip['aski_main_file_only']),
         ('ASKI_MAIN_FILE_WRITE', fip['aski_main_file_write']),
         ('ASKI_MAIN_PATH', fmp['path_aski_main_files']),
         ('ASKI_outfile', 'will be overwritten'),
         ('ASKI_output_ID', 'will be overwritten'),
         ('ASKI_DECONVOLVE_STF', aski_decon),
         ('COMPUTE_PHASE_END_TIMES', fip['compute_phase_end_times']),
         ('PATH_PHASE_END_TIMES', fmp['path_phase_end_times']),
         ('ASKI_df', str(df)), ('ASKI_nf', str(nf)), ('ASKI_jf', " ".join(map(str,ifreq))),
         ('ASKI_DFT_method', fip['aski_dft_method']), ('ASKI_DFT_double', fip['aski_dft_double']),
         ('ASKI_DFT_apply_taper', aski_apply_phase_taper), ('ASKI_DFT_taper_length', str(taper_length)),
         ('ASKI_type_inversion_grid', fmp['aski_type_inversion_grid']),
         ('ASKI_wx', str(fmp['aski_wx'])), ('ASKI_wy', str(fmp['aski_wy'])), ('ASKI_wz', str(fmp['aski_wz'])),
         ('ASKI_rot_X', str(fmp['aski_rot_x'])), ('ASKI_rot_Y', str(fmp['aski_rot_y'])), ('ASKI_rot_Z', str(fmp['aski_rot_z'])),
         ('ASKI_cx', str(fmp['aski_cx'])), ('ASKI_cy', str(fmp['aski_cy'])), ('ASKI_cz', str(fmp['aski_cz'])),
         ('ASKI_TTIME_TABLE_FILE', tt_file),
         ('ASKI_lat_box_center', str(fmp['invgrid_center_lat'])),
         ('ASKI_lon_box_center', str(fmp['invgrid_center_lon'])),
         ('ASKI_rearth', str(fmp['rearth']))])
    print("SPECFEM Par_file_ASKI written: ")

    #  write sequential output description file
    #  If --disp:
    #  Loop over event dictionary to create list of jobs
    #  compute GEMINI_SYNSEIS_FILE, NSTEP, aski_outfile, aski_output_id
    #  If --gt:
    #  Loop over station dictionary to create list of jobs for SPECFEM
    with open(fmp['path_specfem_input']+sequential_sources_file, "w") as fkd:
        if args.disp:
            valid_sources = [key for key in fmp['evlist'].events if type(key) is str]
            fkd.write("{:5d}\n".format(3))                                          # sequential sources mode = 3
        else:
            valid_sources = [key for key in fmp['statlist'].stations if type(key) is str]
            fkd.write("{:5d}\n".format(2))                                          # sequential sources mode = 2
        fkd.write("{:5d}\n".format(len(valid_sources)))                             # number of source lines
        for source in valid_sources:
            if args.disp:                                                           # kernel displacement settings
                gemsynfile = source + '_seis.hdf'
                nstep = read_timing_injected_wavefields(source, float(fip['dt']), fmp['path_injection_seis'])
                aski_outfile = fmp['path_kernel_displacements'] + source
                aski_output_id = source
                fkd.write("{:8d} {:s} {:s} {:s}\n".format(nstep, gemsynfile, aski_output_id, aski_outfile))
            else:                                                                   # kernel Green tensor setting
                aski_output_id = source + "_" + fip['green_tensor_component']
                aski_outfile = fmp['path_kernel_green_tensors'] + aski_output_id
                line = get_sequential_force(fmp['statlist'], source, fip['green_tensor_component'],
                                            fip['hdur_green_tensor'], fmp['rearth'],
                                            fmp['invgrid_center_lat'], fmp['invgrid_center_lon'],
                                            fip['green_tensor_depth_shift'])
                fkd.write("{:s} {:s} {:s}\n".format(line, aski_output_id, aski_outfile))
            if not os_path.exists(aski_outfile + '_OUTPUT_FILES'):
                os_mkdir(aski_outfile + '_OUTPUT_FILES')
                print("mkdir " + aski_outfile + '_OUTPUT_FILES')
            else:
                os_system('rm -rf ' + aski_outfile + '_*')
                os_system('mkdir -p ' + aski_outfile + '_OUTPUT_FILES')
                print(source + ": remove existing files:")
                print('rm -rf ' + aski_outfile + '_*')
                print('mkdir -p ' + aski_outfile + '_OUTPUT_FILES')
#
    print("sequential sources file written:")
    print("End log at " + get_time_string())


#  allow module to be run as a script
if __name__ == "__main__":
    prepareSpecfemKernel(sys.argv[1:])
