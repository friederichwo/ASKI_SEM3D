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
from askipy.iterparKeys import ipkeys
from askipy.inputParameter import inputParameter
from askipy.helperFunctions import check_keys, check_file

def fromIterParfile(file):
    """return a dictionary with information taken from iter parfile
    """
    check_file(file)
    fip = {}
    iterpar = inputParameter(file)
    noKeys = iterpar.keysNotPresent(ipkeys)
    check_keys(noKeys)

    # integer parameters
    for key in ['ITERATION_STEP_NUMBER_OF_FREQ', 'NPROC', 'MAX_NUM_CG_ITERATIONS', 'NITER_WINDOW_STA',
                'NITER_WINDOW_LTA', 'NSTEP']:
        fip[key.lower()] = iterpar.ival(key)

    # float parameters
    for key in ['RESIDUAL_SCATTERING_TIME_KERNEL_TAPER', 'VP_SEDIMENT', 'CALIBRATION_ERROR_FACTOR',
                'MAX_TIME_RESIDUAL', 'VTK_COORDS_SCALING_FACTOR', 'DT', 'HDUR_GREEN_TENSOR',
                'GREEN_TENSOR_DEPTH_SHIFT', 'ASKI_DFT_taper_length_kd', 'ASKI_DFT_taper_length_gt',
                'DSMBOOST', 'BOUNDCHOKE']:
        fip[key.lower()] = iterpar.fval(key)

    # string parameters
    for key in ['MODEL_PROPERTIES_INVERTED_FOR', 'VTK_GEOMETRY_TYPE', 'GREEN_TENSOR_COMPONENT',
                'FILE_ASKI_INVERTED_MODEL', 'ASKI_DFT_method']:
        fip[key.lower()] = iterpar.sval(key)

    # logical parameters, keep them as strings
    for key in ['SCALE_VTK_COORDS', 'MOVIE_SURFACE', 'MOVIE_VOLUME', 'IMPOSE_ASKI_INVERTED_MODEL',
                'COMPUTE_PHASE_END_TIMES', 'ASKI_MAIN_FILE_ONLY', 'ASKI_MAIN_FILE_WRITE',
                'ASKI_DFT_double', 'ASKI_DFT_apply_taper']:
        fip[key.lower()] = iterpar.sval(key)

    # lists
    fip['vscal_smoothing_mantle'] = iterpar.flist('VSCAL_SMOOTHING_MANTLE', 3)
    fip['vscal_smoothing_crust'] = iterpar.flist('VSCAL_SMOOTHING_CRUST', 3)
    fip['vscal_damping_mantle'] = iterpar.flist('VSCAL_DAMPING_MANTLE', 3)
    fip['vscal_damping_crust'] = iterpar.flist('VSCAL_DAMPING_CRUST', 3)
    fip['iteration_step_index_of_freq'] = iterpar.ilist('ITERATION_STEP_INDEX_OF_FREQ', fip['iteration_step_number_of_freq'])

    return fip
