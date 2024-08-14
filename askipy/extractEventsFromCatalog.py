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
from obspy.core.event import read_events
from askipy.helperFunctions import check_file


def extractEventsFromCatalog(clargs):
    ap = ArgumentParser(description="Extract ASKI events from AlpArray CMT catalog based on a selection file",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("selection", help="Input file with Marcel's event ids of selected events")
    ap.add_argument("--cat", help="name of catalog file", default='alparray_cmt_catalog.qml')
    ap.add_argument("--out", help="name of output file", default='ASKI_events_S_selected')
    args = ap.parse_args(clargs)
    selfile = args.selection
    catfile = args.cat
    outfile = args.out
    check_file(selfile)
    check_file(catfile)

    #  read event catalog and quality file
    cat = read_events(catfile, format='QUAKEML')
    with open(selfile,'r') as fq:
        sel = fq.readlines()
    evlist = [line.split()[1] for line in sel[1:]]
    qlist = [0.5*(float(line.split()[4])+float(line.split()[5])) for line in sel[1:]]

    #  Walk through events in catalog and check for occurrence in selection
    #  Construct Marcels event id from centroid origin times
    #  and get other infos needed for ASKI_events_S type output
    nev = 0
    foundlist = []
    for ev in cat.events:
        for orig in ev.origins:
            if orig.origin_type == 'centroid':
                t = orig.time
                evid = '{:4d}{:02d}{:02d}_{:02d}{:02d}{:02d}.a'.format(t.year,t.month,t.day,t.hour,t.minute,t.second)
                tstring = '{:4d}{:02d}{:02d}_{:02d}{:02d}{:02d}_{:09d}'. \
                    format(t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond * 1000)
                lat = orig.latitude
                lon = orig.longitude
                depth = orig.depth

        mt = ev.focal_mechanisms[0].moment_tensor.tensor
        for m in ev.magnitudes:
            if m.magnitude_type == 'Mwc': mag = m.mag

        #   Add to aski_event_lines if found in selection
        try:
            j = evlist.index(evid)
        except ValueError:
            continue
        nev = nev+1
        aski_event_line = '{:s} {:s} {:9.3f} {:10.3f} {:10.1f} {:6.2f}  0.00  1  {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e} {:10.2e}'.\
               format(evid, tstring, lat, lon, depth, mag, mt.m_rr, mt.m_tt, mt.m_pp, mt.m_rt, mt.m_rp, mt.m_tp)
        foundlist.append([qlist[j],aski_event_line])

    #  sort aski_event_lines according to quality
    #  and write to file
    sorted_aski_events = sorted(foundlist, key=lambda el: el[0], reverse=True)
    fout = open(outfile, 'w')
    fout.write('S \n')
    for el in sorted_aski_events:
        fout.write(el[1] + '\n')
    fout.close()
    print('Found ',nev,' events of ',len(evlist))


#  allow module to be run as a script
if __name__ == "__main__":
    extractEventsFromCatalog(sys.argv[1:])
