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
from askipy.mainparKeys import mpkeys
from askipy.inputParameter import inputParameter
from askipy.readEventStationFile import stationList, eventList
from askipy.helperFunctions import check_directory, check_file, check_in_mainpath, check_keys

def fromMainParfile(file, **kwargs):
    """
    return a dictionary with information taken from main parfile
    kwargs:  implemented keys:
        key = it:         value = iteration step (default=current iteration)
        key = itref:      value = iteration step for reference (only used if there)
        key = checkdir:   value = true (check existence of directories), false (do not check), default=True
        key = checkfile:  value = true (check existence of files), false (do not check), default=True
    """
    check_file(file)
    checkdir = True; checkfile = True
    if 'checkdir' in kwargs: checkdir = kwargs['checkdir']
    if 'checkfile' in kwargs: checkfile = kwargs['checkfile']
    fmp = {'iteration_step_path_ref': 'None', 'path_output_files_ref': 'None'}

    mainpar = inputParameter(file)
    noKeys = mainpar.keysNotPresent(mpkeys)
    check_keys(noKeys)

    #  choose lower case version of keyword as key for mainpar dictionary
    #  treat all values as strings
    for key in mpkeys:
        fmp[key.lower()] = mainpar.sval(key)

    #  check if script calling this function is run from main inversion folder
    check_in_mainpath(fmp['main_path_inversion'])

    #  overwrite non-string parameter values by converted ones
    for key in ['REARTH', 'INVGRID_CENTER_LAT', 'INVGRID_CENTER_LON', 'MEASURED_DATA_FREQUENCY_STEP',
                'UNIT_FACTOR_MEASURED_DATA', 'ASKI_wx', 'ASKI_wy', 'ASKI_wz', 'ASKI_rot_X', 'ASKI_rot_Y',
                'ASKI_rot_Z', 'ASKI_cx', 'ASKI_cy', 'ASKI_cz']:
        fmp[key.lower()] = mainpar.fval(key)
    fmp['current_iteration_step'] = mainpar.ival('CURRENT_ITERATION_STEP')
    fmp['invgrid_dimensions'] = mainpar.ilist('INVGRID_DIMENSIONS',3)
    fmp['pseudo_mesh_spacing'] = mainpar.flist('PSEUDO_MESH_SPACING',3)
    fmp['measured_data_number_of_freq'] = mainpar.ival('MEASURED_DATA_NUMBER_OF_FREQ')
    fmp['measured_data_index_of_freq'] = mainpar.ilist('MEASURED_DATA_INDEX_OF_FREQ', fmp['measured_data_number_of_freq'])

    # check files in main iversion folder
    if checkfile:
        for key in ['file_station_list', 'file_event_list', 'property_correlation_file', 'file_gt_receiver_locations']:
            check_file(fmp[key])

    # check directories in main inversion folder
    if checkdir:
        for key in ['path_gemini', 'path_measured_data', 'path_measured_seis', 'path_injection_seis',
                    'path_specfem_input', 'path_specfem_databases', 'path_phase_end_times']:
            check_directory(fmp[key])

    #  event and station list
    fmp['statlist'] = stationList(fmp['file_station_list'], list_type='standard')
    fmp['evlist'] = eventList(fmp['file_event_list'], list_type='standard')
    fmp['gtreclist'] = stationList(fmp['file_gt_receiver_locations'], list_type='standard')

    # put aski background model under specfem input
    fmp['file_aski_background_model'] = fmp['path_specfem_input'] + fmp['file_aski_background_model']
    if checkfile: check_file(fmp['file_aski_background_model'])

    # put ASKI time tables under gemini
    fmp['aski_ttime_table_file_kd'] = fmp['path_gemini'] + fmp['aski_ttime_table_file_kd']
    fmp['aski_ttime_table_file_gt'] = fmp['path_gemini'] + fmp['aski_ttime_table_file_gt']
    if checkfile:
        check_file(fmp['aski_ttime_table_file_kd'])
        check_file(fmp['aski_ttime_table_file_gt'])

    #  deal with iteration specific folders under iteration step path
    #  put paths under iteration step path
    iter_step_base = fmp['iteration_step_path']
    fmp['iter_step_base'] = iter_step_base
    outpath = fmp['path_output_files']
    if 'it' in kwargs:
        it = kwargs['it']
    else:
        it = fmp['current_iteration_step']

    fmp['iteration_step_path'] = iter_step_base + '{:03d}'.format(it) + '/'
    for key in ['path_aski_main_files', 'path_kernel_displacements', 'path_kernel_green_tensors',
          'path_sensitivity_kernels', 'path_synthetic_data', 'path_output_files', 'path_dmspace',
          'path_vtk_files', 'path_specfem_logs']:
        fmp[key] = fmp['iteration_step_path'] + fmp[key]
        if checkdir: check_directory(fmp[key])

    fmp['parfile_iteration_step'] = fmp['iteration_step_path'] + fmp['parfile_iteration_step']
    if checkfile: check_file(fmp['parfile_iteration_step'])

    #  previous iteration path
    if it > 1:
        fmp['prev_iteration_step_path'] = iter_step_base + '{:03d}'.format(it-1) + '/'
        if checkdir: check_directory(fmp['prev_iteration_step_path'])
    else:
        fmp['prev_iteration_step_path'] = fmp['iteration_step_path']

    #  iteration from which reference data are taken
    if 'itref' in kwargs:
        itref = kwargs['itref']
        if itref > 0:
            fmp['iteraton_step_path_ref'] = iter_step_base + '{:03d}'.format(itref) + '/'
            fmp['path_output_files_ref'] = fmp['iteraton_step_path_ref'] + outpath
            if checkdir: check_directory(fmp['iteraton_step_path_ref'])
            if checkdir: check_directory(fmp['path_output_files_ref'])

    return fmp
