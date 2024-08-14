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
import cartopy.crs as ccrs
import cartopy.feature as feature
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.readEventStationFile import eventList
from askipy.helperFunctions import check_file
from askipy.fromMainParfile import fromMainParfile


def eventsPolarPlot(clargs):
    ap = ArgumentParser(description="Plot azimuthal earthquake distribution",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("eventfile", help="ASKI_S type event file")
    ap.add_argument("--range", help="Index range of events to be plotted", default='all')
    ap.add_argument("--xaski", help="Exclude events from ASKI_events_S", action='store_true')
    ap.add_argument("--noevid", help="Omit event ids in plot", action='store_true')
    args = ap.parse_args(clargs)

    eventfile = args.eventfile
    if args.range != 'all':
        evidx = [int(e) for e in args.range.split()]
    xaski = args.xaski
    check_file(eventfile)

    #  get info from main parfile
    fmp = fromMainParfile(args.main_parfile)

    clat = fmp['invgrid_center_lat']
    clon = fmp['invgrid_center_lon']
    if xaski:
        aski_event_list = fmp['evlist']

    evlist = eventList(eventfile, list_type='standard')
    nev = evlist.nev
    if args.range == 'all':
        evidx = [0,nev]

    #  open figure
    coastline = feature.COASTLINE
    proj = ccrs.AzimuthalEquidistant(central_longitude=clon, central_latitude=clat)
    fig = plt.figure(figsize=(15,15))
    #fig.suptitle('Eventfile: ' + eventfile, fontsize=14)
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    ax.add_feature(coastline, zorder=10)

    #  plot epicenters
    if xaski:
        slon = np.array([float(evlist.events[i]['slon']) for i in range(evidx[0],evidx[1]) if evlist.events[i]['evid'] not in aski_event_list.events])
        slat = np.array([float(evlist.events[i]['slat']) for i in range(evidx[0],evidx[1]) if evlist.events[i]['evid'] not in aski_event_list.events])
        mag = np.array([float(evlist.events[i]['mag']) for i in range(evidx[0],evidx[1]) if evlist.events[i]['evid'] not in aski_event_list.events])
        evsel = [evlist.events[i]['evid'] for i in range(evidx[0],evidx[1]) if evlist.events[i]['evid'] not in aski_event_list.events]
        qual = np.array(np.linspace(1.0, 0.0, slon.size))
        print("Additional events found: ",slon.size)
    else:
        slon = np.array([float(evlist.events[i]['slon']) for i in range(evidx[0],evidx[1])])
        slat = np.array([float(evlist.events[i]['slat']) for i in range(evidx[0],evidx[1])])
        mag = np.array([float(evlist.events[i]['mag']) for i in range(evidx[0],evidx[1])])
        evsel = [evlist.events[i]['evid'] for i in range(evidx[0],evidx[1])]
        qual = np.array(np.linspace(1.0, 0.0, slon.size))

    xyz = proj.transform_points(ccrs.PlateCarree(), slon, slat)
    cf = ax.scatter(xyz[:,0], xyz[:,1], s=150*(mag-5.0)**2, marker='o', c=qual, edgecolors='k', cmap='viridis')
    if not args.noevid:
        for j in range(slon.size):
            ax.text(xyz[j,0],xyz[j,1],evsel[j], fontsize=12)

    fig.colorbar(cf, shrink=0.8, pad=0.05)

    #  AlpArray center
    ax.plot(0.0, 0.0, 'k*' , markersize=50)

    #  plot distance circles 30, 60, 90 and 120
    for dis in [30.0, 60.0, 90.0]:
        x,y = proj.transform_point(clon, clat-dis, ccrs.PlateCarree())
        r = np.sqrt(x**2+y**2)
        x = r*np.cos(np.linspace(0.0,2.*np.pi,101))
        y = r*np.sin(np.linspace(0.0,2.*np.pi,101))
        ax.plot(x, y, linestyle = '--', color = 'grey')

    fig.savefig(fmp['main_path_inversion'] + 'events_polar_plot.png', dpi=300, bbox_inches='tight')
    plt.show()


#  allow module to be run as a script
if __name__ == "__main__":
    eventsPolarPlot(sys.argv[1:])


