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
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from askipy.fromMainParfile import fromMainParfile
from askipy.readEventStationFile import stationList
from askipy.helperFunctions import check_file


def stationMap(clargs):
    ap = ArgumentParser(description="Plot a station map",
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("main_parfile", help="Name of ASKI main parameter file")
    ap.add_argument("--stationfile", help="ASKI_S type station file", default='ASKI_stations_S')
    ap.add_argument("--idns", help="Add net.station names in plot", action='store_true')
    ap.add_argument("--idn", help="Add netcodes in plot", action='store_true')
    ap.add_argument("--ids", help="Add station names in plot", action='store_true')
    ap.add_argument("--figsize", help="2-tuple with figure size", default='16 14')
    ap.add_argument("--hwlat", help="half width of map in latitude", default='6.0')
    ap.add_argument("--hwlon", help="half width of map in longitude", default='10.0')
    args = ap.parse_args(clargs)
    stationfile = args.stationfile
    check_file(stationfile)
    wx, wy = args.figsize.split()
    hwlat = float(args.hwlat)
    hwlon = float(args.hwlon)

    #  read ASKI main parameter file
    fmp = fromMainParfile(args.main_parfile)
    clat = fmp['invgrid_center_lat']
    clon = fmp['invgrid_center_lon']

    #  get station info
    statlist = stationList(stationfile, list_type='standard')
    nstat = statlist.nstat
    slon = np.array([float(statlist.stations[i]['lon']) for i in range(0,nstat)])
    slat = np.array([float(statlist.stations[i]['lat']) for i in range(0,nstat)])

    #  map setup
    borders = cf.BORDERS
    coastline = cf.COASTLINE
    land = cf.LAND.with_scale('50m')
    ocean = cf.OCEAN.with_scale('50m')
    proj = ccrs.LambertConformal(central_longitude = clon, central_latitude = clat,
                                 standard_parallels=(clat-3.0, clat+3.0))
    font = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': 16}

    #  figure
    fig = plt.figure(figsize=(float(wx), float(wy)))
    ax = plt.axes(projection=proj)

    ax.add_feature(land, color=(0.5,0.5,0,0.5), zorder=0)
    ax.add_feature(ocean, color=(0.0,0.0,1.0,0.5), zorder=0)
    ax.add_feature(borders, zorder=10, linewidth=0.5, color='k')
    ax.add_feature(coastline, zorder=10, linewidth=0.5, color='k')

    ax.set_extent((clon-hwlon,clon+hwlon,clat-hwlat,clat+hwlat))
    ax.set_title('Station map ' + ' / ' + stationfile , fontdict=font)
    ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--',
                 x_inline=False, y_inline=False, rotate_labels=False)

    #  stations
    ax.scatter(slon, slat, s=15, marker='o', color=(0.,0.,0.), edgecolors='k', transform=ccrs.PlateCarree())
    if args.idns or args.idn or args.ids:
        for j in range(nstat):
            netcode = statlist.stations[j]['netcode']
            staname = statlist.stations[j]['staname']
            if args.idn:
                ax.text(slon[j],slat[j],netcode, fontsize=12, transform=ccrs.PlateCarree())
            if args.ids:
                ax.text(slon[j],slat[j],staname, fontsize=12, transform=ccrs.PlateCarree())
            if args.idns:
                ax.text(slon[j],slat[j],netcode+'.'+staname, fontsize=12, transform=ccrs.PlateCarree())

    fig.savefig(fmp['main_path_inversion'] + 'station_map.png', dpi=300, bbox_inches='tight')
    plt.show()


#  allow module to be run as a script
if __name__ == "__main__":
    stationMap(sys.argv[1:])
