#----------------------------------------------------------------------------
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
#----------------------------------------------------------------------------
#
def readSliceSpecifications(file, inbase, typ):
    """
    Read file with specifictions of slices for later plotting.
    Put information into a nested dictionary of the following struture:
    Top level: keys =   Slice_0, Slice_1 .....
               values = dictionaries with slice data
    Slice data: keys =    Name, Parameter, Filename, Field_names
                values =  slice name, slice paramters, file name, field names
    inbase: path + basename of HDF-slice data
    typ:  either VERTICAL or HORIZONTAL (default)
    """
    slspec = dict()
    with open(file, 'r') as fsl:
        while True:
            line = fsl.readline()
            if not line:
                break
            if line.strip() == typ:
                nsl = int(fsl.readline().strip())
                for i in range(nsl):
                    key = "Slice_{:d}".format(i)
                    slspec[key] = dict()
                    pro = fsl.readline().strip()
                    slspec[key]['Name'] = pro
                    slspec[key]['Parameter'] = (fsl.readline().strip()).split()
                    if typ == 'VERTICAL':
                        slspec[key]['Filename'] = inbase + '_vslice_' + pro + '.hdf'
                    else:
                        slspec[key]['Filename'] = inbase + '_hslice_' + pro + '.hdf'
                    slspec[key]['Field_names'] = (fsl.readline().strip()).split()
    return slspec

