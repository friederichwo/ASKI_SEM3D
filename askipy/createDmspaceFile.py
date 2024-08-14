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
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.fromIterParfile import fromIterParfile


def createDmspaceFile(clargs):
    """
    Create a data model space info file
    :param clargs: list of command line arguments
    """
    ap = ArgumentParser(description="Create a data-model-space-info file for selected events and stations",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of main ASKI parameter file")
    ap.add_argument("--ifreq", help="List of minus-separated frequency indices (overwrites those from iter_parfile", default='None')
    ap.add_argument("--it", help="Iteration step for which dmsp-file is generated", default='current')
    print(clargs)
    args = ap.parse_args(clargs)

    # get ASKI main parameters
    if args.it == 'current':
        fmp = fromMainParfile(args.main_parfile)
    else:
        fmp = fromMainParfile(args.main_parfile, it=int(args.it))

    # get ASKI iteration parameters
    fip = fromIterParfile(fmp['parfile_iteration_step'])

    # get optional string of frequencies (e.g. 12-14-16-18)
    ifreq = args.ifreq

    # set name of output dmspace file
    dmspfile = fmp['path_dmspace'] + "ASKI_dmspace"
    if ifreq != 'None':
        if ifreq == 'allf':
            dmspfile = dmspfile + '_' + 'allf'
        else:
            dmspfile = dmspfile + '_' + ifreq

    #  get number of model properties from command line
    nprops = len(fip['model_properties_inverted_for'].split())
    ncomp = len(fip['green_tensor_component'].split())

    #  open dmspace info output file and write "MODEL VALUES"
    #  plus headers to DATA SAMPLES block
    with open(dmspfile, 'w') as fout:
        fout.write("MODEL VALUES\n")
        fout.write(str(nprops) + " " + fip['model_properties_inverted_for'] + "\n")
        fout.write("DATA SAMPLES\n")
        fout.write("WEIGHTING BY_PATH\n")
        fout.write("MASKING NONE\n")
        fout.write("COMPONENTS ALL\n")
        fout.write(str(ncomp) + " " + fip['green_tensor_component'] + "\n")
        fout.write("FREQUENCIES ALL\n")
        if ifreq == 'None':
            txt = "{:d} {:s}".format(fip['iteration_step_number_of_freq'],
                                     ' '.join(map(str,fip['iteration_step_index_of_freq'])))
        else:
            if ifreq == 'allf':
                txt = "{:d} {:s}".format(fmp['measured_data_number_of_freq'],
                                         ' '.join(map(str,fmp['measured_data_index_of_freq'])))
            else:
                txt = "{:d} {:s}".format(len(ifreq.split('-')), ' '.join(ifreq.split('-')))
        fout.write(txt + "\n")

    # create a list of net.statname from station list (sorted according to occurrence)
    selected_netstatnames = [key for key in fmp['statlist'].stations if type(key) is str]

    #  Loop over selected events
    selected_events = [key for key in fmp['evlist'].events if type(key) is str]

    #  Use a dictionary to associate events and weights with a given net.staname
    npath = 0
    mapEventWeightToNetsta = dict()
    for evid in selected_events:
        dmsp_block = fmp['path_measured_seis'] + evid + "/ASKI_dmspace_block"
        with open(dmsp_block, "r") as fcs:
            nstat = fcs.readline().rstrip()
            while True:
                line = fcs.readline()
                if not line:
                    break
                ns = line.split()[1]
                weight = line.split()[2]
                if ns not in selected_netstatnames:
                    #print('Not in station file: ',ns)
                    continue
                if ns in mapEventWeightToNetsta:
                    mapEventWeightToNetsta[ns].append([evid, weight])
                else:
                    mapEventWeightToNetsta[ns] = [[evid, weight]]
                npath = npath+1

    #  sort map according to length of its values
    sortedMap = {key: val for key, val in sorted(mapEventWeightToNetsta.items(), key=lambda ele: len(ele[1]), reverse=True)}

    #  write paths to ASKI_dmspace file
    with open(dmspfile, "a") as fout:
        fout.write(str(npath) + '\n')
        for ns in sortedMap:
            for evid, weight in sortedMap[ns]:
                txt = "{evid:<20}{ns:<10}{weight:<18}\n"
                fout.write(txt.format(evid=evid, ns=ns, weight=weight))

    sys.stdout.write(str(npath) + ' ' + "paths written to file " + dmspfile + "\n")


#  allow module to be run as a script
if __name__ == "__main__":
    createDmspaceFile(sys.argv[1:])

