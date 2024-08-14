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
import h5py
import axesRotation as aR
import numpy as np


def writeSpecfemStations(file, statlist):
    nstat = statlist.nstat
    with open(file, 'w') as fout:
        for i in range(nstat):
            d = statlist.stations[i]
            txt = "{name:<6}{net:<3}{lon:>15s}{lat:>15s}{zero:15.3f}{alt:>15s}\n"
            fout.write(txt.format(name=d['staname'], net=d['netcode'], lon=d['lon'], lat=d['lat'], zero=0.0, alt=d['alt']))


def set_parfile(filename, keys_vals):
    """
    edit parameter files Par_file and Par_file_ASKI according to
    settings of parameters for this script
    """
    values = dict(keys_vals)
    keys = list(values.keys())
    #
    # iterate over all lines and modify line if necessary
    #
    lines_new = []
    with open(filename, "r") as fp:
        for line in fp:
            # ignore empty lines, comment lines and invalid lines (which do not contain any '=' in front of a comment)
            if line.strip() == '' or line.strip()[0:1] == '#' or '=' not in line.split('#')[0]:
                lines_new.append(line)
                continue
            #
            # if the key,value pair of this valid line is to be modified, do so
            #
            key_line = line.split('=')[0].strip()
            val_line = line.split('=')[1].split('#')[0].strip()
            if key_line in keys:
                # first, remove newline character from end of line and append it to modified line
                line = line.strip('\n')
                # replace old value by new value, but keep all whitespace and commentary of this line as was before
                key_part = line.split('=')[0]
                val_part = line.split('=')[1].split('#')[0]
                # comment part ends on newline character
                if '#' in line:
                    comment_part = '#' + line.split('=')[1].split('#')[1] + '\n'
                else:
                    comment_part = '\n'
                if val_line == '':
                    val_part = ' ' + values[key_line] + ' '
                else:
                    val_part = val_part.replace(val_line, values[key_line])
                line = key_part + '=' + val_part + comment_part
                # remove key of this line from keys list
                # to check (in the end) if all keys were found in the file
                keys.remove(key_line)
            # add the line (if modified or not) to list of new lines
            lines_new.append(line)
    #
    # check if there are any keys which were not found on valid lines in the file
    if len(keys) > 0:
        raise Exception("could not find the following parameters in parameter file '" + filename + "': " +
                        "'" + "', '".join(keys) + "'")
    #
    # if every key was found and the respective value was modified, write modified lines to file
    with open(filename, 'w') as fp:
        fp.writelines(lines_new)


def read_timing_injected_wavefields(event, dts, path):
    """Calculate NSTEP for SPECFEM from timing of injected wavefields.
       Read tred und tend from <timing> attribute of HDF file, take difference
       and divide by SPECFEM DT"""
    with h5py.File(path + event + '_seis.hdf', 'r') as fid:
        [_, _, _, tred, tend, _] = fid.attrs["timing"]
        [dtg, _] = fid.attrs["samplingIntervalLength"]
    lgem = dtg * int((tend - tred) / dtg)  # length of injected seismogram
    return int(lgem / float(dts))  # number of dts-intervals = nstep


def get_sequential_force(statlist, nsname, gtcomp, hdur, rearth, lat_box_center, lon_box_center, depth_shift):
    """Produce one line of the sequential source description file for the 
       kernel Green tensor case on the basis of the station file. Optionally
       move source point away from surface by depth_shift > 0 (m].
    """ 
    valid_force_components = ['CX', 'CY', 'CZ', 'UP']
    if gtcomp in valid_force_components:
        x = float(statlist.stations[nsname]['lat'])   # = X coordinate in SPECFEM box
        y = float(statlist.stations[nsname]['lon'])   # = Y coordinate in SPECFEM box
        z = float(statlist.stations[nsname]['alt'])   # = Z coordinate in SPECFEM box
        r, delta, xi = aR.calculateSphericalCoordinates(x, y, z + float(rearth), lat_box_center, lon_box_center)
        #
        #  new SPECFEM box coordinates for depth shifted source
        if depth_shift > 0.0:
            x,y,z = aR.calculateBoxCoordinates(r-depth_shift, 90.0-delta, xi, lat_box_center, lon_box_center) 
            z = z-float(rearth)
        if gtcomp == 'CX':
            fcomp = [1.0, 0.0, 0.0]
        elif gtcomp == 'CY':
            fcomp = [0.0, 1.0, 0.0]
        elif gtcomp == 'CZ':
            fcomp = [0.0, 0.0, 1.0]
        elif gtcomp == 'UP':
            # calculate FX, FY and FZ components for force vector pointing local up at station on sphere
            # transform local up pointing force at station into box cartesian frame
            fg = aR.vectorLCfromLS(np.deg2rad(delta), np.deg2rad(xi), [1., 0., 0.])
            fl = aR.vectorLCfromGC(0.5 * np.pi - np.deg2rad(lat_box_center), np.deg2rad(lon_box_center), fg)
            fr = aR.vectorRCfromLC(0.5 * np.pi, fl)
            fcomp = fr
        else:
            fcomp = [0.0, 0.0, 0.0]
        #
        #  exchange x and y here because it will be exchanged again by get_elevation in locate_sources
        line = "0.000 {:.3f} {:13.3f} {:13.3f} {:13.3f} 1.0 {:.6f} {:.6f} {:.6f}".\
            format(hdur, y, x, z, fcomp[0], fcomp[1], fcomp[2])
        return line

    raise Exception("invalid force component specified")
