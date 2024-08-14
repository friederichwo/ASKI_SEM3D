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
import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.helperFunctions import check_file


############################################################
# Functions
############################################################
def plot3curves(ax, x, p, q, r, ylab):
    ax.grid(True, linestyle='-')
    ax.set_ylabel(ylab)
    ax.plot(x, p, color='k', marker='o')
    ax.plot(x, q, color='b', marker='v')
    ax.plot(x, r, color='r', marker='^')


def plot2curves(ax, x, p, q, ylab):
    ax.grid(True, linestyle='-')
    ax.set_ylabel(ylab)
    ax.plot(x, p, color='k', marker='o')
    ax.plot(x, q, color='r', marker='v')


def plothist2(ax, p, q, n, drange, alfa, ylab):
    ax.grid(True, linestyle='-')
    ax.set_xlabel(ylab)
    ax.set_ylabel('Count')
    ax.hist(p, bins=n, alpha=alfa, color='b',range=drange)
    ax.hist(q, bins=n, alpha=alfa, color='g',range=drange)


def phasediff(phi1, phi2):
    d0 = phi1-phi2 
    d1 = np.where(d0 > +180., d0-360., d0)
    d2 = np.where(d1 < -180., d1+360., d1)
    return d2


def plotKernelLinearSystemData(clargs):
    ap = ArgumentParser(description="Plot data used in kernelLinearSystem",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--hist", help="Only plot histograms", action='store_true')
    group.add_argument("--data", help="Show individual data", default='store_true')
    ap.add_argument("--it", help="take data from this iteration", default='current')
    ap.add_argument("--range", help="restrict to data index range (e.g. 0 99)", default='all')
    ap.add_argument("--nbin", help="number of bins for histograms", type=int, default=200)
    ap.add_argument("--norm", help="normalize residuals by amplitude of data", action='store_true')
    args = ap.parse_args(clargs)

    idx_range = args.range
    nbin = int(args.nbin)
    do_hist = args.hist
    do_data = args.data
    normalize = args.norm

    #  get info from main parfile
    if args.it == 'current':
       fmp = fromMainParfile(args.main_parfile)
       it = fmp['current_iteration_step']
    else:
        it = int(args.it)
        fmp = fromMainParfile(args.main_parfile, it=it)

    hdffile_ini = fmp['path_output_files']  + 'klsdata_init' + '.hdf'
    check_file(hdffile_ini)
    hdffile_fin = fmp['path_output_files']  + 'klsdata_final' + '.hdf'
    check_file(hdffile_fin)
    df = fmp['measured_data_frequency_step']

    print('Take initial data from iteration {0:4d} and file {1}'.format(it, hdffile_ini))
    print('Take final data from iteration {0:4d} and file {1}'.format(it, hdffile_fin))

    #  read data from HDF files
    with h5py.File(hdffile_ini, 'r') as fid:
        dset_ini = np.copy(fid["klsdata"])
        ifreq_ini = np.copy(fid["ifreq"])
    with h5py.File(hdffile_fin, 'r') as fid:
        dset_fin = np.copy(fid["klsdata"])
        ifreq_fin = np.copy(fid["ifreq"])

    ndata = dset_ini.shape[1]
    if idx_range == 'all':
        j1 = 0; j2 = ndata
    else:
        j1 = int(idx_range.split()[0])
        j2 = int(idx_range.split()[1])
    x = np.arange(j1, j2, 2)/2

    #  Data samples and synthetic samples (first index)
    #  Data is ordered in (Re, Im) for each frequency, component, and path (second index)
    redat = dset_ini[0, j1:j2:2]
    imdat = dset_ini[0, j1+1:j2+1:2]
    magdat = np.hypot(redat, imdat)
    phidat = np.arctan2(imdat, redat)*180./np.pi

    resyn_ini = dset_ini[1, j1:j2:2]
    imsyn_ini = dset_ini[1, j1+1:j2+1:2]
    magsyn_ini = np.hypot(resyn_ini, imsyn_ini)
    phisyn_ini = np.arctan2(imsyn_ini, resyn_ini)*180./np.pi

    resyn_fin = dset_fin[1, j1:j2:2]
    imsyn_fin = dset_fin[1, j1+1:j2+1:2]
    magsyn_fin = np.hypot(resyn_fin, imsyn_fin)
    phisyn_fin = np.arctan2(imsyn_fin, resyn_fin)*180./np.pi

    resabs_ini = np.hypot(redat-resyn_ini, imdat-imsyn_ini)
    resabs_fin = np.hypot(redat-resyn_fin, imdat-imsyn_fin)
    timediff_ini = phasediff(phidat, phisyn_ini)/(360.*df*ifreq_ini[j1:j2:2])
    timediff_fin = phasediff(phidat, phisyn_fin)/(360.*df*ifreq_fin[j1:j2:2])

    if normalize:
        norm = magdat
    else:
        norm = np.ones((j2-j1+1)//2)

    #  plotting of data versus synthetics
    if do_data:
        figds = plt.figure(figsize=(20, 20))
        ax1 = figds.add_subplot(4, 1, 1)
        ax1.set_title('Data and Synthetics')
        plot3curves(ax1, x, redat, resyn_ini, resyn_fin, 'Real part')

        ax2 = figds.add_subplot(4, 1, 2)
        plot3curves(ax2, x, imdat, imsyn_ini, imsyn_fin, 'Imaginary part')

        ax3 = figds.add_subplot(4, 1, 3)
        plot3curves(ax3, x,  magdat, magsyn_ini, magsyn_fin, 'Magnitude')

        ax4 = figds.add_subplot(4, 1, 4)
        plot3curves(ax4, x, phidat, phisyn_ini, phisyn_fin, 'Phase')
        figds.savefig(fmp['path_output_files'] + 'kls_data_vs_syn.png', dpi=300, bbox_inches='tight')

        #  plotting of residuals
        figres = plt.figure(figsize=(20, 20))
        bx1 = figres.add_subplot(4, 1, 1)
        bx1.set_title('Residuals')
        plot2curves(bx1, x, (redat-resyn_ini)/norm, (redat-resyn_fin)/norm, 'Real part')

        bx2 = figres.add_subplot(4, 1, 2)
        plot2curves(bx2, x, (imdat-imsyn_ini)/norm, (imdat-imsyn_fin)/norm, 'Imaginary part')

        bx3 = figres.add_subplot(4, 1, 3)
        plot2curves(bx3, x, resabs_ini/norm, resabs_fin/norm, 'Magnitude')

        bx4 = figres.add_subplot(4, 1, 4)
        plot2curves(bx4, x, timediff_ini, timediff_fin, 'Time shift')
        figres.savefig(fmp['path_output_files'] + 'kls_residuals.png', dpi=300, bbox_inches='tight')

    # plotting of histograms for residuals
    if do_hist:
        fighist = plt.figure(figsize=(20, 20))
        fighist.suptitle('Blue bars: ' + hdffile_ini + '\n'+ 'Green bars: ' + hdffile_fin, fontsize=14)
        cx1 = fighist.add_subplot(2, 2, 1)
        cx1.set_title('Residuals')
        plothist2(cx1, (redat-resyn_ini)/norm, (redat-resyn_fin)/norm, nbin, (-8.,8.), 0.5, 'Real part')

        cx2 = fighist.add_subplot(2, 2, 2)
        plothist2(cx2, (imdat-imsyn_ini)/norm, (imdat-imsyn_fin)/norm, nbin, (-8.,8.), 0.5, 'Imaginary part')

        cx3 = fighist.add_subplot(2, 2, 3)
        plothist2(cx3, resabs_ini/norm, resabs_fin/norm, nbin, (0.,8.), 0.5, 'Magnitude')

        cx4 = fighist.add_subplot(2, 2, 4)
        plothist2(cx4, timediff_ini, timediff_fin, nbin, (-5.,5.), 0.5, 'Time shift')
        fighist.savefig(fmp['path_output_files'] + 'kls_histo.png', dpi=300, bbox_inches='tight')

    plt.show()


#  allow module to be run as a script
if __name__ == "__main__":
    plotKernelLinearSystemData(sys.argv[1:])
