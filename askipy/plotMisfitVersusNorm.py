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
import sys
from os import path as os_path
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile


def getModelNorm(fmp,itmax):
    it = 1
    norm2 = [0.0]
    while it <= itmax:
        iterpath = fmp['iter_step_base'] + '{:03d}'.format(it) + '/'
        normfile = fmp['path_output_files'].replace(fmp['iteration_step_path'], iterpath) + 'modelnorm.log'
        if not os_path.exists(normfile):
            print("Warning: File " + normfile + 'does not exist')
            break

        with open(normfile, 'r') as fn:
            for line in fn:
                if line.strip().find('Quadratic / rms model norm') >= 0:
                    norm2.append( float(line.strip().split(':')[1].split()[0]) )

        it = it+1
    return norm2


def getMisfits(fmp,ext,itmax):
    it = 1
    mf = []
    while it <= itmax:
        iterpath = fmp['iter_step_base'] + '{:03d}'.format(it) + '/'
        misfitfile = fmp['path_output_files'].replace(fmp['iteration_step_path'], iterpath) + 'misfits_masked_' + ext +'.log'
        if not os_path.exists(misfitfile):
            print("Warning: File " + misfitfile + 'does not exist')
            break

        with open(misfitfile, 'r') as fm:
            for line in fm:
                if line.strip().find('Total quadratic / rms misfit') >= 0:
                    mf.append( float(line.strip().split(':')[1].split()[0]) )

        it = it+1
    return mf


def plotMisfitVersusNorm(clargs):
    """
    Take misfits and norm from output files of iterations and create a graph
    """
    ap = ArgumentParser(description="Take misfits and norm from output files of iterations and create a graph",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of main ASKI parameter file")
    ap.add_argument("--it", help="Iteration step up to which misfit values are taken", default='current')
    ap.add_argument("--noshow", help="do not show figures / only to png", action='store_true')
    ap.add_argument("--ylim", help="lower and upper misfit limit", default='2 6')
    args = ap.parse_args(clargs)
    main_parfile = args.main_parfile
    noshow = args.noshow
    mflim = [float(v) for v in args.ylim.split()]

    extlist = ["6-8-10-12", "9-11-13-15", "12-14-16-18", "allf"]
    flist = ["30-40-50-60", "45-55-65-75", "60-70-80-90", "allf"]
    finv = ["30-40-50-60", "45-55-65-75", "60-70-80-90"]
    symb = ['D', 'H', 's','o']
    col = ['k', 'b', 'g', 'm']

    if args.it == 'current':
        fmp = fromMainParfile(main_parfile, checkdir=False)
        itmax = fmp['current_iteration_step']
    else:
        itmax = int(args.it)
        fmp = fromMainParfile(main_parfile, it=itmax, checkdir=False)

    # setup plot
    fig = plt.figure(figsize=(10, 15))
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Squared model norm', fontsize=18)
    ax.set_ylabel('Misfit', fontsize=18)
    ax.grid(True, linestyle='--')
    plt.ylim(bottom=mflim[0], top=mflim[1])

    for i, ext in enumerate(extlist):
        norm2 = getModelNorm(fmp, itmax)
        mf = getMisfits(fmp, ext, itmax)
        print("Frequency set: ",flist[i])
        print("Iteration      Misfit       Norm2")
        for j in range(len(mf)):
            print("{:9d}{:12.3f}{:12.3f}".format(j, mf[j], norm2[j]))
            if i==0: ax.text(norm2[j], mf[j]*0.96, "It {:d}".format(j), ha = 'center', va='top', fontsize=18)
        ax.plot(norm2[:len(mf)], mf, marker=symb[i], color=col[i], markersize=15, label='Frequencies: ' + flist[i], linestyle='--')

    ax.legend(fontsize=18)
    fig.savefig(fmp['path_output_files'] + 'misfit_vs_norm.png', dpi=300, bbox_inches='tight')
    if not noshow: plt.show()

#  allow module to be run as a script
if __name__ == "__main__":
    plotMisfitVersusNorm(sys.argv[1:])
