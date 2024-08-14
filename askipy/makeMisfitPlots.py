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
#   Make misfit plots
#
import sys
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.plotMisfits import plotMisfits
from askipy.plotMisfitVersusNorm import plotMisfitVersusNorm
from askipy.fromMainParfile import fromMainParfile


def makeMisfitPlots(clargs):
    ap = ArgumentParser(description="make varisou misfit plots",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--it", help="Iteration step from which misfit values are taken", default='current')
    ap.add_argument("--opt",help="Colon separated key/value plot options (vmax:4, drange:2/4, ylim:0/2)", default = None)
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile

    if args.it == 'current':
        fmp = fromMainParfile(main_parfile, checkdir=False)
        it = str(fmp['current_iteration_step'])
    else:
        it = args.it

    # plot options
    vmax = '8.0'; drange = '-8.0 +8.0'; ylim = '2 6'
    plot_options = args.opt.split()
    print('Plot options: ',plot_options)
    kw = [f.split(':')[0] for f in plot_options]
    val = [f.split(':')[1] for f in plot_options]
    pldict = dict(zip(kw, val))
    if 'vmax' in pldict: vmax = pldict['vmax']
    if 'drange' in pldict: drange = ' '.join(pldict['drange'].split('/'))
    if 'ylim' in pldict: ylim = ' '.join(pldict['ylim'].split('/'))

    ifreq_set = ['6-8-10-12', '9-11-13-15', '12-14-16-18', 'allf']

    # plot residual statistics for all frequency sets, masked only
    for ext in ifreq_set:
        plotMisfits([main_parfile, '--all', '--ext', ext, '--masked', '--it', it, '--itr', '1',
                    '--figsize', '16 12', '--noshow', '--drange', drange])

    # plot event_frequency misfit matrix for all frequency sets, masked only
    for ext in ifreq_set:
        plotMisfits([main_parfile, '--evfreq', '--ext', ext, '--masked', '--it', it, '--itr', '1',
                    '--vmax', vmax, '--figsize', '18 16', '--noshow'])

    # plot misfit versus norm
    plotMisfitVersusNorm([main_parfile, '--it', it, '--noshow', '--ylim', ylim])


#  allow module to be run as a script
if __name__ == "__main__":
    makeMisfitPlots(sys.argv[1:])

