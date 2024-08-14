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
#
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.heatmap import heatmap


def mask_below(a, val):
    return np.ma.array(a, mask=(a < val))


def readSpecificResiduals(file):
    with h5py.File(file, 'r') as fid:
        resabs = np.copy(fid["absres"])
        restime = np.copy(fid["timeres"])
        sigma = np.copy(fid["sigma"])
        ifreq = np.copy(fid["ifreq"])
    return resabs, restime, sigma, ifreq


def readFlattenedResiduals(file):
    with h5py.File(file, 'r') as fid:
        resabs = np.copy(fid["resabs"])
        restime = np.copy(fid["restime"])
        res_real = np.copy(fid["res_real"])
        res_imag = np.copy(fid["res_imag"])
    return res_real, res_imag, resabs, restime


def readEventFrequencyMatrix(file):
    evlist = []
    mf = []
    start = False
    with open(file, 'r') as fd:
        for line in fd:
            if line.strip().find('Event-frequency-wise misfit') < 0 and not start:
                continue
            else:
                start = True
                if line.strip().find('.a') >= 0:
                    h = line.strip().split()
                    evlist.append(h[1])
                    mf.append(list(map(float, h[2:])))

                if line.strip().find('Index') >= 0:
                    flist = list(map(float, line.strip().split()[2:]))
                    freqlist = ["{:3.0f}".format(f*1000.) for f in flist]

    return freqlist, evlist, mf


def plotSpecificResiduals(ax1, ax2, ifreq, x, h, vrange):
    """ Plot event or station specific residuals for one frequency
    """
    ax1.set_title('Frequency index {:3d}'.format(ifreq), loc='left')
    ax1.grid(True, linestyle='--')
    ax1.bar(x, h, align='center', color='b', alpha=0.5)
    ax2.set_title('Frequency index {:3d}'.format(ifreq), loc='left')
    ax2.grid(True, linestyle='--')
    ax2.hist(h, bins=50, alpha=1.0, color='b', range=vrange)


def addRefSpecificResiduals(ax1, ax2, x, h, vrange):
    ax1.bar(x, h, align='center', color='g', alpha=0.5)
    ax2.hist(h, bins=50, alpha=0.5, color='g', range=vrange)


def plothist2(ax, p, q, n, drange, alfa, ylab):
    ax.tick_params(labelsize=16)
    ax.grid(True, linestyle='-')
    ax.set_xlabel(ylab, fontsize=18)
    ax.set_ylabel('Count', fontsize=18)
    ax.hist(p, bins=n, alpha=alfa, color='b', range=drange)
    ax.hist(q, bins=n, alpha=alfa, color='g', range=drange)


def plotMisfits(clargs):
    """Main program: plot various kinds of misfits
       clargs: list of command line arguments
    """
    ap = ArgumentParser(description="Plot various residuals of FWI",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--evid", help="Show residuals for selected event", default='None')
    group.add_argument("--nsname", help="Show residuals for selected net-station", default='None')
    group.add_argument("--all", help="Show residuals over all data", action='store_true')
    group.add_argument("--evfreq", help="Plot event-frequency misfit matrix", action='store_true')
    ap.add_argument("--it", help="use residuals for this iteration", type=int, default=1)
    ap.add_argument("--itr", help="use residuals for this iteration as reference", type=int, default=1)
    ap.add_argument("--masked", help="use masked residuals", action='store_true')
    ap.add_argument("--ext", help="use residuals with ifreq-extension", default='None')
    ap.add_argument("--ifreq", help="show residuals for selected frequency indices", default='None')
    ap.add_argument("--nbin", help="number of bins for histograms", type=int, default=100)
    ap.add_argument("--vmax", help="Saturate color scale for event-frequency misfit matrix", type=float, default=None)
    ap.add_argument("--figsize", help="2-tuple with figure size", default='20 15')
    ap.add_argument("--noshow", help="do not show figures / only to png", action='store_true')
    ap.add_argument("--drange", help="2-tuple with data range for flattened residuals", default='-8.0 8.0')
    args = ap.parse_args(clargs)

    #  extract command line arguments
    it = args.it
    itref = args.itr
    evid = args.evid
    nsname = args.nsname
    use_all_data = args.all
    evfmat = args.evfreq
    wx, wy = args.figsize.split()
    nbin = args.nbin
    use_masked_dmspace = args.masked
    resext = args.ext
    ifreq_sel = args.ifreq
    vmax = args.vmax
    noshow = args.noshow
    drange = tuple([float(v) for v in args.drange.split()])

    #  get info from main parfile
    fmp = fromMainParfile(args.main_parfile, it = it, itref = itref, checkdir=False)

    #  flattened residuals case
    if use_all_data:
        fnbase = 'resduals_flattened'
        if use_masked_dmspace:
            fnbase = 'residuals_flattened_masked'
        if resext != 'None':
            fnbase = 'residuals_flattened_' + resext
        if resext != 'None' and use_masked_dmspace:
            fnbase = 'residuals_flattened_masked_' + resext
        resfile = fmp['path_output_files'] + fnbase + '.hdf'
        resreffile = fmp['path_output_files_ref'] + fnbase + '.hdf'
        pngfile = fmp['path_output_files'] + fnbase + '.png'
        print("residuals from :", resfile)
        print("ref residuals from :", resreffile)

        res_real, res_imag, resabs, restime = readFlattenedResiduals(resfile)
        res_real_ref, res_imag_ref, resabs_ref, restime_ref = readFlattenedResiduals(resreffile)

        #  histograms for real, imag, abs and time
        fig = plt.figure(figsize=(float(wx), float(wy)))
        fig.suptitle('Blue bars: ' + resfile + '\n' + 'Green bars: ' + resreffile, fontsize=18, y=0.95)
        cx1 = fig.add_subplot(2, 2, 1)
        plothist2(cx1, res_real, res_real_ref, nbin, drange, 0.5, 'Real part')
        cx2 = fig.add_subplot(2, 2, 2)
        plothist2(cx2, res_imag, res_imag_ref, nbin, drange, 0.5, 'Imaginary part')
        cx3 = fig.add_subplot(2, 2, 3)
        plothist2(cx3, resabs, resabs_ref, nbin, (0., drange[1]), 0.5, 'Magnitude')
        cx4 = fig.add_subplot(2, 2, 4)
        plothist2(cx4, restime, restime_ref, nbin, (0.5*drange[0], 0.5*drange[1]), 0.5, 'Time shift')
        fig.savefig(pngfile, dpi=300, bbox_inches='tight')

    #  event frequency wise misfit matrix
    elif evfmat:
        fnbase = 'misfits'
        if use_masked_dmspace:
            fnbase = 'misfits_masked'
        if resext != 'None':
            fnbase = 'misfits_' + resext
        if resext != 'None' and use_masked_dmspace:
            fnbase = 'misfits_masked_' + resext
        logfile = fmp['path_output_files'] + fnbase + '.log'
        pngfile = fmp['path_output_files'] + fnbase + '.png'
        print("event-frequency misfit matrix from :", logfile)

        freqlist, evlist, mf = readEventFrequencyMatrix(logfile)
        nevh = len(evlist)//2
        fig = plt.figure(figsize=(float(wx), float(wy)))
    #   fig.suptitle('Event frequency wise misfits: ' + logfile, fontsize=14, y=0.90)
        cbar_kw = dict(shrink=0.55, pad=0.05)
        cx1 = fig.add_subplot(1, 2, 1)
        cx1.tick_params(labelsize=16)
        im, cbar = heatmap(np.array(mf[0:nevh]), evlist[0:nevh], freqlist, ax=cx1, 
                           cbar_kw=cbar_kw, cbarlabel="Misfit", vmin=0.0, vmax=vmax)
        cx2 = fig.add_subplot(1, 2, 2)
        cx2.tick_params(labelsize=16)
        im, cbar = heatmap(np.array(mf[nevh+1:]), evlist[nevh+1:], freqlist, ax=cx2, 
                           cbar_kw=cbar_kw, cbarlabel="Misfit", vmin=0.0, vmax=vmax)
        fig.savefig(pngfile, dpi=300, bbox_inches='tight')

    #  specific residuals case
    else:
        fnbase = 'resduals_specific'
        if use_masked_dmspace:
            fnbase = 'residuals_specific_masked'
        if resext != 'None':
            fnbase = 'residuals_specific_' + resext
        if resext != 'None' and use_masked_dmspace:
            fnbase = 'residuals_specific_masked_' + resext
        resfile = fmp['path_output_files'] + fnbase + '.hdf'
        resreffile = fmp['path_output_files_ref'] + fnbase + '.hdf'
        pngbase = fmp['path_output_files'] + fnbase
        print("residuals from :",resfile)
        print("ref residuals from :", resreffile)

        resabs, restime, sigma, ifreq = readSpecificResiduals(resfile)
        resabs_ref, restime_ref, sigma_ref, ifreq_ref = readSpecificResiduals(resreffile)
        nf, nev, nstat = resabs.shape

        if ifreq_sel != 'None':
            ifreq_sel = list(map(int, ifreq_sel.split()))
            nfsel = len(ifreq_sel)
        else:
            nfsel = nf
            ifreq_sel = list(ifreq)

        if evid != 'None':
            evidx = list(fmp['evlist'].events.keys()).index(evid)
            print('Event ' + evid + ' has index {:d} in event list'.format(evidx))
            heading = evid
            statidx = -1
        else:
            statidx = list(fmp['statlist'].stations.keys()).index(nsname)
            print('Station ' + nsname + ' has index {:d} in station list'.format(statidx))
            heading = nsname
            evidx = -1

        #  residuals of absolute value
        fig = plt.figure(figsize=(float(wx), float(wy)))
        gs = gridspec.GridSpec(nfsel, 2, width_ratios=[2, 1])
        fig.suptitle('Error normalized residuals of absolute value for {:s} '.format(heading), fontsize=14, y=0.95)
        if evidx >= 0:
            x = np.linspace(1, nstat, nstat)
            pngfile = pngbase + '_' + evid + '_resabs.png'
        else:
            x = np.linspace(1, nev, nev)
            pngfile = pngbase + '_' + nsname + '_resabs.png'
        for jf in range(nfsel):
            kf = list(ifreq).index(ifreq_sel[jf]) 
            ax1 = fig.add_subplot(gs[jf, 0])
            ax2 = fig.add_subplot(gs[jf, 1])
            if evidx >= 0:
                h = mask_below(resabs[kf, evidx, :], 0.0)
                href = mask_below(resabs_ref[kf, evidx, :], 0.0)
            else:
                h = mask_below(resabs[kf, :, statidx], 0.0)
                href = mask_below(resabs_ref[kf, :, statidx], 0.0)
            plotSpecificResiduals(ax1, ax2, ifreq_sel[jf], x, h, (0., 8.))
            addRefSpecificResiduals(ax1, ax2, x, href, (0., 8.))
        fig.savefig(pngfile, dpi=300, bbox_inches='tight')

        #  residuals of time
        fig = plt.figure(figsize=(float(wx), float(wy)))
        gs = gridspec.GridSpec(nfsel, 2, width_ratios=[2, 1])
        fig.suptitle('Time residuals for {:s} '.format(heading), fontsize=14, y=0.95)
        for jf in range(nfsel):
            kf = list(ifreq).index(ifreq_sel[jf]) 
            ax1 = fig.add_subplot(gs[jf, 0])
            ax2 = fig.add_subplot(gs[jf, 1])
            if evidx >= 0:
                h = mask_below(restime[kf, evidx, :], -999.0)
                href = mask_below(restime_ref[kf, evidx, :], -999.0)
                pngfile = pngbase + '_' + evid + '_restime.png'
            else:
                h = mask_below(restime[kf, :, statidx], -999.0)
                href = mask_below(restime_ref[kf, :, statidx], -999.0)
                pngfile = pngbase + '_' + nsname + '_restime.png'
            plotSpecificResiduals(ax1, ax2, ifreq_sel[jf], x, h, (-4., 4.))
            addRefSpecificResiduals(ax1, ax2, x, href, (-4., 4.))
        fig.savefig(pngfile, dpi=300, bbox_inches='tight')

    if not noshow: plt.show()


#  allow module to be run as a script
if __name__ == "__main__":
    plotMisfits(sys.argv[1:])
