#!/usr/bin/env python3
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
import sys
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.helperFunctions import check_directory, get_time_string


def findAllStationsUsed(clargs):
    ap = ArgumentParser(description="Find all stations used after having computed source wavelets",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of main ASKI parameter file")
    ap.add_argument("--minev", help="Min number of events per station", default='1')
    args = ap.parse_args(clargs)
    minev = int(args.minev)

    sys.stdout.write("Start log at " + get_time_string() + "\n")
    sys.stdout.write("Arguments: " + str(args) + "\n")

    #  get info from main parfile
    fmp = fromMainParfile(args.main_parfile)

    #  Loop over selected events
    #  Collect station lines from clean station file for all events
    #  prepend netstaname to all lines for sorting
    selected_events = [key for key in fmp['evlist'].events if type(key) is str]
    all_stations = []
    for evid in selected_events:
        check_directory(fmp['path_measured_seis'] + evid)
        with open(fmp['path_measured_seis'] + evid + '/ASKI_clean_station_file', 'r') as fstat:
            stat_specs = fstat.readlines()
            for line in stat_specs[1:]:
                lspl = line.strip().split()
                ns = lspl[0]+'.'+lspl[1]
                lspl.insert(0,ns)
                all_stations.append(lspl)

    #  Sort and make unique the list of all stations
    #  and count how often it occurs
    sorted_all_stat_specs = sorted(all_stations, key=lambda item: item[0])
    unique_stat_specs = []
    specref = sorted_all_stat_specs[0]
    count = 0
    for spec in sorted_all_stat_specs:
        if spec[0] != specref[0]:
            unique_stat_specs.append([count]+specref[1:])
            specref = spec
            count = 1
        else:
            count = count + 1

    #  Sort station list according to occurrence and write to file
    ntot = 0; nkeep = 0; ndrop = 0
    with open(fmp['path_measured_seis'] + 'stations_sorted_minev_{:d}'.format(minev), 'w') as fout:
        with open(fmp['path_measured_seis'] + 'dropped_stations_sorted_maxev_{:d}'.format(minev-1), 'w') as fdrop:
            fout.write('S' + '\n')
            fdrop.write('S' + '\n')
            for spec in sorted(unique_stat_specs, key=lambda item: item[0], reverse=True):
                ntot = ntot+1
                txt1 = "{name:<6}{net:<3}{lat:15.3f}{lon:15.3f}{alt:15.3f}\n"
                txt2 = "{name:<6}{net:<3}{lat:15.3f}{lon:15.3f}{alt:15.3f}{count:6d}\n"
                if spec[0] >= minev:
                    fout.write(txt1.format(name=spec[1], net=spec[2], lat=float(spec[3]), lon=float(spec[4]), alt=float(spec[5])))
                    sys.stdout.write(txt2.format(count=spec[0], name=spec[1], net=spec[2], lat=float(spec[3]),
                                                 lon=float(spec[4]), alt=float(spec[5])))
                    nkeep = nkeep+1
                else:
                    fdrop.write(txt1.format(name=spec[1], net=spec[2], lat=float(spec[3]), lon=float(spec[4]), alt=float(spec[5])))
                    sys.stdout.write(txt2.format(count=spec[0], name=spec[1], net=spec[2], lat=float(spec[3]),
                                                 lon=float(spec[4]), alt=float(spec[5])))
                    ndrop = ndrop+1

    sys.stdout.write("Min number of events per station: " + str(minev) + "\n")
    sys.stdout.write("Total number of stations used: " + str(ntot) + "\n")
    sys.stdout.write("Total number of stations kept : " + str(nkeep) + "\n")
    sys.stdout.write("Total number of stations dropped: " + str(ndrop) + "\n")


#  allow module to be run as a script
if __name__ == "__main__":
    findAllStationsUsed(sys.argv[1:])
