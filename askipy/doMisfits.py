#!/usr/bin/env python3
# ----------------------------------------------------------------------------
#   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
#   Compute misfits
#
import sys
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.runAskiProgram import runAskiProgram
from askipy.createDmspaceFile import createDmspaceFile


def doMisfits(clargs):
    ap = ArgumentParser(description="create dmspace files and compute misfits for various frequency sets",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--it", help="Iteration step for which misfits are computed", default='current')
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile
    it = args.it
    itm = it
    if it == 'current':
        itm = '0'

    ifreq_set = ['6 8 10 12', '9 11 13 15', '12 14 16 18']

    # create dmspace file with corresponding extension
    for ifreq in ifreq_set:
        createDmspaceFile([main_parfile, '--ifreq', ifreq, '--it', it])

    createDmspaceFile([main_parfile, '--ifreq', 'allf', '--it', it])

    # compute misfits for all frequency sets, normal and masked
    for ifreq in ifreq_set:
        ext = '-'.join(ifreq.split())
        runAskiProgram(['misfits', main_parfile, '-ext', ext, '-iter', itm])
        runAskiProgram(['misfits', main_parfile, '-ext', ext, '-masked', '-iter', itm])

    # all frequencies
    runAskiProgram(['misfits', main_parfile, '-ext', 'allf', '-iter', itm])
    runAskiProgram(['misfits', main_parfile, '-ext', 'allf', '-masked', '-iter', itm])


#  allow module to be run as a script
if __name__ == "__main__":
    doMisfits(sys.argv[1:])

