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
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.readEventStationFile import eventList
from askipy.helperFunctions import check_file


def eventsDepthHisto(clargs):
    ap = ArgumentParser(description="Plot depth histogram of events",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--eventfile", help="ASKI_S type event file", default='None')
    ap.add_argument("--nbin", help="Number of depth bins in histogram", default='100')
    args = ap.parse_args(clargs)
    eventfile = args.eventfile
    nbin = int(args.nbin)

    #  get info from main parfile
    fmp = fromMainParfile(args.main_parfile)

    aski_event_file = fmp['file_event_list']
    if eventfile == 'None':
        eventfile = aski_event_file
        evlist = fmp['evlist']
    else:
        check_file(eventfile)
        evlist = eventList(eventfile, list_type='standard')
        
    nev = evlist.nev
    sdepth = np.array([float(evlist.events[i]['sdepth']) for i in range(nev)])*1.e-3
    srad = fmp['rearth']*1.e-3-sdepth

    fig = plt.figure(figsize=(15,8))
    fig.suptitle('Eventfile: ' + eventfile, fontsize=14)
    ax = fig.add_subplot(1, 2, 1)
    ax.hist(sdepth, bins=nbin, alpha=1.0, color='b',range=(0,700.0))
    ax = fig.add_subplot(1, 2, 2)
    ax.hist(srad, bins=nbin, alpha=1.0, color='b',range=(5671.0,6371.0))
    fig.savefig(fmp['main_path_inversion'] + 'events_depth_histo.png', dpi=300, bbox_inches='tight')

    plt.show()

    srad_sort = np.sort(srad)
    for r in srad_sort:
        sys.stdout.write('{:10.2f}'.format(r))


#  allow module to be run as a script
if __name__ == "__main__":
    eventsDepthHisto(sys.argv[1:])


