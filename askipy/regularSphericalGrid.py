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
#  Python class for spherical regular grid
#
import logging
import numpy as np
import h5py
import copy
from os import path as os_path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from askipy.getTicks import getTicks
import askipy.ctrans as ctrans
import askipy.rsgrid as rsgrid

logger = logging.getLogger(__name__)


class RegularSphericalGrid(object):
    """
    Class specifying instances of a regular spherical grid (RSG) defined in 3D.
    Further documentation is placed in doc/sphinx/build/html/index.html.

        :key = shade:    value = 'nearest' or 'gouraud'
        :key = extrema:  value = set fixed extrema for values (default = None)
        :key = cliprel:  value = clipping level as fraction of absmax (default=-1.0)
        :key = clipabs:  value = clipping level as fraction of absmax (default=-1.0)
        :key = norm:     value = False or True
        :key = log:      value = False or True, take base-10 logarithm of values (default = 'False')
        :key = sym:      value = False or True
        :key = ticklabelsize: value = size of font (default = 14)
        :key = niso:     value = number of isolines
        :key = conlab:   value = False or True
        :key = percent:  value = False or True (multiply field values with 100)
        :key = wx:       value = width of figure window (default = 16)
        :key = wy:       value = height of figure window (default = 8)
        :key = cmap:     value = name of color map (default = coolwarm_r)
        :key = cmapiso:  value = name of color map for isolines (default = seismic_r)
        :key = isofmt:   value = format for iso lines (default = '%.2f')
        :key = colbar:   value = False or True (default = True)
        :key = cbpad:    value = padding space between graph and colorbar (default = 0.05)
        :key = cbshrink: value = fraction of axes length = height of colorbar (default = 0.75)
    """

    plot_defaults = {
        'dzline': 100000, 'ddist': 100000, 'dll': 1.0, 'dllfmt': '{:2.0f}', 'shade': 'gouraud',
        'extrema': None, 'cliprel': -1.0, 'clipabs': -1.0, 'dnorm': False, 'norm': False,
        'log': False, 'sym': False, 'ticklabelsize': 14, 'niso': -1, 'conlab': False,
        'percent': False, 'wx': 16, 'wy': 8, 'isofmt': '%.2f', 'cmap': 'coolwarm_r',
        'cmapiso': 'seismic_r', 'colbar': False, 'cbpad': 0.05, 'cbshrink': 0.75
    }

    def __init__(self, dims, spacings, minvals, pole, rsgname, field_name_string, field_values):
        """
        Initialize regular spherical grid object from data:
        :param dims: dimensions of grid (tuple: nlon, nlat, nr)
        :param spacings: spacings of grid (tuple: dlon, dlat, dr) in (rad, rad, m)
        :param minvals: minimum values of grid coordinates (tuple: lonmin, latmin, rmin) in (rad, rad, m)
        :param pole: geographical latitude and longitude of pole the grid coordinates refer to (rad, rad)
        :param rsgname: name of grid
        :param field_name_string: colon-separated string composed of names for data fields defined on grid
        :param field_values: array of dimensions (nfield,npts) containing field data
        """
        self.nlon, self.nlat, self.nr = dims
        self.dlon, self.dlat, self.dr = spacings
        self.lonmin, self.latmin, self.rmin = minvals
        self.lon_pole, self.lat_pole = pole
        self.rsgname = rsgname
        self.field_name_string = field_name_string
        self.field = field_values
        self.nfield = self.field.shape[0]
        self.npts = self.nr*self.nlat*self.nlon
        self.typ = '3D'
        self.dims = dims; self.spacings = spacings; self.minvals = minvals; self.pole = pole
        if self.nr == 1 and self.nlon > 1:  self.typ = 'HS'
        if self.nr == 1 and self.nlon == 1: self.typ = 'GC'
        if self.nlat > 1 and self.nlon == 1: self.typ = 'VS'
        if self.nlat == 1: self.typ = 'RL'

    @classmethod
    def readFromHDF(cls, filename):
        """
        Class method as second constructor to initialize instance of grid from a HDF file
        :param filename: full path name of HDF file
        :return: an instance of the class
        """
        if not os_path.exists(filename):
            raise Exception('File ' + filename + ' does not exist')
        with h5py.File(filename, 'r') as fid:
            rsgname = str(fid.attrs['rsgname'], 'utf-8')
            field_name_string = str(fid.attrs['field_names'], 'utf-8')
            pole = tuple(fid.attrs['crs_pole'])
            dims = tuple(fid.attrs['dimensions'])
            spacings = tuple(fid.attrs['grid_spacings_rad_m'])
            minvals = tuple(fid.attrs['minimum_values_rad_m'])
            field = np.copy(fid['field_values'])
        return cls(dims, spacings, minvals, pole, rsgname, field_name_string, field)


    def extractVerticalSlice(self, lat1d, lon1d, lat2d, lon2d, nr, ndel, slname):
        """
        Extract a vertical slice from a 3D regular spherical grid.
        Slice is taken along a great circle from point (lat1,lon1) to point (lat2,lon2).
        Coordinates of these points refer to the true North Pole.
        The great circle is sampled at ndel points. Vertically the grid is sampled at nr points.
        :param lat1d: Latitude of starting point of great circle (deg)
        :param lon1d: Longitude of starting point of great circle (deg)
        :param lat2d: Latitude of end point of great circle (deg)
        :param lon2d: Longitude of end point of great circle (deg)
        :param nr: Number of samples in new slice along radius
        :param ndel: Number of samples in new slice along great circle
        :param slname: Name of new slice
        :return: a RSG-object containing vertical slice and one containing the great circle
        """
        if self.typ != '3D':
            raise Exception('extractVerticalSlice: RSG object is not of type 3D')

        # convert starting and end point to rad
        [lat1, lon1, lat2, lon2] = np.deg2rad([lat1d, lon1d, lat2d, lon2d])

        # calculate coordinates of end point relative to a pole at starting point
        # giving us length and azimuth of great circle
        [delmax], [xigc] = ctrans.from_gs_to_ls(lat1, lon1, [lat2], [lon2])
        delmax = 0.5*np.pi-delmax

        # coordinates of great circle relative to pole at starting point
        delta = np.linspace(0.0, delmax, ndel)
        xi = xigc*np.ones(ndel)

        # coordinates of great circle relative to North Pole
        latgc, longc = ctrans.from_ls_to_gs(lat1, lon1, 0.5*np.pi-delta, xi)

        # create a gc-RSG object
        gc = RegularSphericalGrid((1,ndel,1), (0.0,delmax/(ndel-1),0.0), (xigc,0.0,self.rmin), (lat1,lon1),
                                  slname, 'latitude:longitude', np.array([latgc, longc]))

        # coordinates of great circle relative to pole of grid
        latp, lonp = ctrans.from_gs_to_ls(self.lat_pole, self.lon_pole, latgc, longc)

        # pole is on opposite hemisphere of points, subtract pi
        lonp = np.where(np.cos(self.lon_pole-longc) < 0.0, lonp-np.pi, lonp)
        lonp = np.where(lonp < -np.pi, lonp+2*np.pi, lonp)
        
        # radii
        rp = np.linspace(self.rmin, self.rmin+(self.nr-1)*self.dr, nr)

        # compute field values at points
        vsfield = rsgrid.field_values_at_points(self.dims, self.spacings, self.minvals, self.field, self.typ,
                                                'VS', rp, latp, lonp, nr*ndel)

        # create a VS-RSG object
        vslice = RegularSphericalGrid((1,ndel,nr), (0.0,delmax/(ndel-1),self.dr), (xigc,0.0,self.rmin), (lat1,lon1),
                                      slname, self.field_name_string, vsfield)
        return vslice, gc


    def extractHorizontalSlice(self, r, nlat, nlon, slname):
        """
        Extract horizontal slice from 3D grid at given radius.
        :param r: Radius of silice
        :param nlat: Number of points along latitude
        :param nlon: Number of points along longitude
        :param slname: Name of slice
        :return: an HS-RSG object
        """
        if self.typ != '3D':
            raise Exception('extractHorizontalSlice: RSG object is not of type 3D')

        # create points of hslice in RSG specific coordinates
        latp = np.linspace(self.latmin, self.latmin+(self.nlat-1)*self.dlat, nlat)
        lonp = np.linspace(self.lonmin, self.lonmin+(self.nlon-1)*self.dlon, nlon)
        dlat = (self.nlat-1)*self.dlat/(nlat-1)
        dlon = (self.nlon-1)*self.dlon/(nlon-1)
        rp = r*np.ones(1)
        hsfield = rsgrid.field_values_at_points(self.dims, self.spacings, self.minvals, self.field, self.typ,
                                                'HS', rp, latp, lonp, nlat*nlon)

        # create a HS-RSG object
        hslice = RegularSphericalGrid((nlon,nlat,1), (dlon,dlat,0.0), (self.lonmin,self.latmin,r), self.pole,
                                      slname, self.field_name_string, hsfield)
        return hslice


    def setupQuadmesh(self):
        """
        Generate Cartesian coordinates of cell centers and cell edges
        for vertical 2D slice.
        See matplotlib docu on Quadmesh for a description
        """
        if self.typ != 'VS':
            raise Exception('setupQuadmesh: Slice object is not of type VS (vertical slice)')
        ndel = self.nlat
        delmin = self.latmin
        delstep = self.dlat
        delmax = delmin+(ndel-1)*delstep
        delmid = 0.5*(delmin+delmax)
        delta = np.linspace(delmin, delmax, num=ndel)
        xg = np.zeros((self.nr, ndel))
        yg = np.zeros((self.nr, ndel))

        # grid points as cell centers
        for k in range(0, self.nr):
            r = self.rmin + k*self.dr
            xg[k, :] = r*np.sin(delta-delmid)
            yg[k, :] = r*np.cos(delta-delmid)

        #  cell edges
        xedge = np.zeros((self.nr+1, ndel+1))
        yedge = np.zeros((self.nr+1, ndel+1))
        for k in range(0, self.nr+1):
            r = self.rmin - 0.5*self.dr + k*self.dr
            xedge[k, 0:ndel] = r*np.sin(delta-0.5*delstep-delmid)
            yedge[k, 0:ndel] = r*np.cos(delta-0.5*delstep-delmid)
            xedge[k, ndel] = r*np.sin(delta[ndel-1]+0.5*delstep-delmid)
            yedge[k, ndel] = r*np.cos(delta[ndel-1]+0.5*delstep-delmid)
        #
        return xg,yg,xedge,yedge


    def computeIsoDepthLines(self,dzline):
        """
        Calculate Cartesian coordinates of iso-depth lines
        :param dzline: spacing between iso-depth lines in meters
        """
        if self.typ != 'VS':
            raise Exception('computeIsoDepthLines: Slice object is not of type VS (vertical slice)')
        if dzline < 0:
            xzl = np.zeros((1, 1))
            yzl = np.zeros((1, 1))
            return xzl,yzl

        ndel = self.nlat
        delmin = self.latmin
        delstep = self.dlat
        delmax = delmin + (ndel - 1) * delstep
        delmid = 0.5 * (delmin + delmax)
        delta = np.linspace(delmin, delmax, num=ndel)
        rsurf = 6371000.0
        nzl = int(np.floor((rsurf-self.rmin)/dzline))
        xzl = np.zeros((nzl, ndel))
        yzl = np.zeros((nzl, ndel))
        for i in range(nzl):
            zline = (i+1)*dzline
            xzl[i,:] = (rsurf-zline)*np.sin(delta-delmid)
            yzl[i,:] = (rsurf-zline)*np.cos(delta-delmid)
        return xzl,yzl


    def computeTopBottomBoundaryLines(self):
        """
        Calculate Cartesian coordinates of top and bottom boundary of vertical section
        """
        if self.typ != 'VS':
            raise Exception('computeTopBottomBoundaryLines: Slice object is not of type VS (vertical slice)')

        ndel = self.nlat
        delmin = self.latmin
        delstep = self.dlat
        delmax = delmin + (ndel - 1) * delstep
        delmid = 0.5 * (delmin + delmax)
        delta = np.linspace(delmin, delmax, num=ndel)
        rtop = self.rmin+(self.nr-1)*self.dr
        xzl = np.zeros((2,ndel))
        yzl = np.zeros((2,ndel))
        xzl[0,:] = rtop*np.sin(delta-delmid)
        yzl[0,:] = rtop*np.cos(delta-delmid)
        xzl[1,:] = self.rmin*np.sin(delta-delmid)
        yzl[1,:] = self.rmin*np.cos(delta-delmid)
        return xzl,yzl


    def computeIsoDistanceLines(self, ddis):
        """
        Calculate coordinates of radial iso distance lines.
        :param ddis: spacing between iso distance lines in meters
        """
        if self.typ != 'VS':
            raise Exception('computeIsoDistanceLines: Slice object is not of type VS (vertical slice)')
        if ddis < 0:
            xdl = np.zeros((1, 1))
            ydl = np.zeros((1, 1))
            return xdl, ydl

        nr = self.nr
        rmin = self.rmin
        rstep = self.dr
        rmax = rmin+(nr-1)*rstep
        dmax = rmin*self.dlat*(self.nlat-1)
        delmid = self.latmin+0.5*(self.nlat-1)*self.dlat
        r = np.linspace(rmin, rmax, num=nr)
        ndl = int(np.floor(dmax/ddis))
        xdl = np.zeros((ndl, nr))
        ydl = np.zeros((ndl, nr))
        for i in range(ndl):
            dis = (i+1)*ddis
            xdl[i,:] = r*np.sin(dis/rmin-delmid)
            ydl[i,:] = r*np.cos(dis/rmin-delmid)
        return xdl, ydl


    def computeLeftRightBoundaryLines(self):
        """
        Calculate coordinates of left and right bounary lines of vertical slice
        """
        if self.typ != 'VS':
            raise Exception('computeLeftRightBoundaryLines: Slice object is not of type VS (vertical slice)')

        nr = self.nr
        rmin = self.rmin
        rstep = self.dr
        rmax = rmin+(nr-1)*rstep
        delmax = self.latmin+self.dlat*(self.nlat-1)
        delmid = self.latmin+0.5*(self.nlat-1)*self.dlat
        r = np.linspace(rmin, rmax, num=nr)
        xdl = np.zeros((2, nr))
        ydl = np.zeros((2, nr))
        xdl[0,:] = r*np.sin(self.latmin-delmid)
        ydl[0,:] = r*np.cos(self.latmin-delmid)
        xdl[1,:] = r*np.sin(delmax-delmid)
        ydl[1,:] = r*np.cos(delmax-delmid)
        return xdl, ydl


    def getExtentHS(self):
        """
        Calculate extent of horizontal slice in degrees.
        """
        if self.typ != 'HS':
            raise Exception('getExtentHS: Slice object is not of type HS (horizontal slice)')
        return np.rad2deg( np.array([self.lonmin, self.lonmin+(self.nlon-1)*self.dlon,
                                     self.latmin, self.latmin+(self.nlat-1)*self.dlat]) )


    def getMeshGridHS(self):
        """
        Return latitude and longitude of grid points of HS in degrees as meshgrid
        """
        if self.typ != 'HS':
            raise Exception('getExtentHS: Slice object is not of type HS (horizontal slice)')
        lat = np.rad2deg( np.linspace(self.latmin, self.latmin+(self.nlat-1)*self.dlat,self.nlat) )
        lon = np.rad2deg( np.linspace(self.lonmin, self.lonmin+(self.nlon-1)*self.dlon,self.nlon) )
        long,latg = np.meshgrid(lon,lat)
        return long, latg


    def getFieldNames(self):
        """
        Return list of field names of grid
        """
        return self.field_name_string.split(':')


    def getReshapedField(self,field_name, depth_normalize=False, normalize=False, cliprel=-1.0,
                         clipabs=-1.0, sym=False, percent=False, extrema=None, log=False):
        """
        Reshapes flattened field array to more than one dimension according to type of slice for given field name.
        :param field_name: name of data field to be processed
        :param depth_normalize: flag for depth normalization of field values
        :param normalize: flag for max-normalization of field data
        :param cliprel: flag and value for clipping of data at a given fraction of the maximum
        :param clipabs: flag and value for clipping of data at agiven value
        :param sym: flag for symmetrising the field data
        :param percent: flag for conversion of field data to percent
        :param extrema: flag and values of min/max-values for color scales and isoline
        :param log: flag for taking base 10 log of values and zero waterlevel
        :return: multi-dimensional array containing modified field data
        """
        fn = self.getFieldNames()
        try:
            jf = fn.index(field_name)
            if self.typ == 'HS':
                val = np.reshape(self.field[jf,:], (self.nlat, self.nlon), order='C')
            elif self.typ == 'VS':
                val = np.reshape(self.field[jf,:], (self.nr, self.nlat), order='C')
            elif self.typ == 'VS':
                val = np.reshape(self.field[jf,:], (self.nr, self.nlat), order='C')
            elif self.typ == '3D':
                val = np.reshape(self.field[jf,:], (self.nr, self.nlat, self.nlon), order='C')
            elif self.typ == 'GC':
                val = np.reshape(self.field[jf,:], self.nlat, order='C')
            elif self.typ == 'RL':
                val = np.reshape(self.field[jf,:], self.nr, order='C')
            else:
                raise Exception('Invalid type of regular spherical grid')

            valmin = np.min(val); valmax = np.max(val)
            logger.info("True extrema of field %s: %f %f", field_name, valmin, valmax)

            # depth normalization
            if depth_normalize and self.typ == 'VS':
                for k in range(0, self.nr):
                    valrms = np.sum(np.fabs(val[k,:]))/self.nlat
                    if valrms > 0:
                        val[k, :] = val[k, :] / valrms

            # normalize to max absolute value
            if normalize:
                val = val / np.max(np.fabs(val))

            # take base 10 log of values and waterlevel at 10**-10
            if log:
                if valmin < 0:
                    logger.error('Minimum of data values is negative. Cut negative values in log-plot.')
                val = np.where(val > 1.0e-10, np.log10(val), -10.0)

            # relative clipping
            if cliprel > 0.0:
                vmax = np.max(np.fabs(val))
                val = np.maximum(val,-cliprel*vmax)
                val = np.minimum(val,+cliprel*vmax)

            # absolute clipping
            if clipabs > 0.0:
                val = np.maximum(val,-clipabs)
                val = np.minimum(val,+clipabs)

            if percent:
                val = val * 1.e2

            # current extreme values
            vmin = np.min(val)
            vmax = np.max(val)

            # extrema
            if extrema:
                vmin = extrema[0]
                vmax = extrema[1]

            # symmetrize
            if sym:
                if abs(vmax) >= abs(vmin):
                    vmin = -vmax
                else:
                    vmax = -vmin
            
            return val,vmin,vmax

        except ValueError:
            logger.error('Field name %s not avilable in slice', field_name)
            exit()


    def plotPoints(self, dp, ax):
        """
        Plot external points into slice
        :param dp: dictionary ordered by index with values <coor>, <size>, <marker>, <color>
        :param ax: axis into which ponts are plotted
        :return: None
        """
        for j in dp:
            c1 = dp[j]['coor'][0]
            c2 = dp[j]['coor'][1]
            s = dp[j]['size']
            marker = dp[j]['marker']
            color = dp[j]['color']
            if self.typ == 'VS':
                ax.scatter(c1, c2, s=s, marker=marker, c=color, edgecolors='none', zorder=5.0)
            else:
                ax.scatter(c1, c2, s=s, marker=marker, c=color, edgecolors='none', transform=ccrs.PlateCarree(),
                           zorder=5.0)



    def plotLines(self, dl, ax):
        """
        Plot external lines into slice
        :param dl: dictionary ordered by index. Value is another dict ordered by index and
                   value as dict with <coor>, <ls>, <lw>, <color>
        :param ax: axis into which lines are plotted
        :return: None
        """
        for j in dl:
            dl2 = dl[j]
            for k in dl2:                        # line dict may contain many lines also ordered by an index key
                c1 = dl2[k]['coor'][:, 0]
                c2 = dl2[k]['coor'][:, 1]
                ls = dl2[k]['ls']
                lw = dl2[k]['lw']
                color = dl2[k]['color']
                if self.typ == 'HS':
                    ax.plot(c1, c2, color=color, linestyle=ls, linewidth=lw, transform=ccrs.PlateCarree(), zorder=5.0)
                else:
                    ax.plot(c1, c2, color=color, linestyle=ls, linewidth=lw, zorder=5.0)


    def plotOverlays(self, dov, ax, c1g, c2g):
        """
        Plot external overlays into slice.
        :param dov: Dictionary ordered by index and value a dict with <array>, <colormap>, <vmin>, <vmax>
        :param ax: axis into which overlay is plotted
        :param c1g: array with x-values of grid coordinates
        :param c2g: array with y values of grid coordinates
        :return:
        """
        if self.typ == 'HS': pass
        for j in dov:
            cmap = dov[j]['colormap']
            ax.contourf(c1g, c2g, dov[j]['array'], cmap=cmap, vmin=dov[j]['vmin'],
                        vmax=dov[j]['vmax'], zorder=1.0, levels=np.linspace(0.0, 1.0, cmap.N))



    def plotText(self, dt, ax):
        """
        Plot text into slice
        :param dt: dictionary ordered by index and value a dict with <text>, <tpos>, <rpos>, <va>, <arrow>
                   <fontweight>, <size>
        :param ax:
        :return:
        """
        for j in dt:
            tc1, tc2 = dt[j]['tpos']
            rc1, rc2 = dt[j]['rpos']
            va = dt[j]['va']
            arw = dt[j]['arrow']
            fw = dt[j]['fontweight']
            if self.typ == 'VS':
                ax.text(tc1, tc2, dt[j]['text'], fontsize=dt[j]['size'], zorder=6.0,
                        fontweight=fw, ha='center', rotation=0.0, va=va)
                if arw:
                    ax.arrow(tc1, tc2, rc1-tc1, rc2-tc2, zorder=6.0,
                             width=0.05, color='k', head_width=0.25, length_includes_head=True)
            else:
                ax.text(tc2, tc1, dt[j]['text'], fontsize=dt[j]['size'], zorder=6.0,
                        fontweight=fw, ha='center', rotation=0.0, va=va, transform=ccrs.PlateCarree())
                if arw:
                    ax.arrow(tc2, tc1, rc2-tc2, rc1-tc1, transform=ccrs.PlateCarree(), zorder=6.0,
                             width=0.05, color='k', head_width=0.25, length_includes_head=True)



    def plotVslice(self, inbase, slname, fn, **kwargs):
        """
        Plot a vertical slice using pcolormesh
        :param inbase: basename of output png-file
        :param slname: name of vertical slice
        :param fn:     name of data field
        kwargs:  implemented keys:
            :key = gc:       value = regular spherical grid object defining the great circle of the vslice
        """
        # first set defaults for plotting options (use deep copy to avoid changing defaults)
        plopt = copy.deepcopy(plot_defaults)

        # overwrite plopt with values from kwargs
        for kw in plopt:
            if kw in kwargs:
                plopt[kw] = kwargs[kw]

        fig = plt.figure(figsize=(plopt['wx'], plopt['wy']))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_axis_off()
        ax.axis('equal')

        # Quadmesh, zlines, dislines, reshape
        xg, yg, xedge, yedge = self.setupQuadmesh()
        xzl, yzl = self.computeIsoDepthLines(plopt['dzline'])
        xdl, ydl = self.computeIsoDistanceLines(plopt['ddist'])
        xtb, ytb = self.computeTopBottomBoundaryLines()
        xlr, ylr = self.computeLeftRightBoundaryLines()
        val, vmin, vmax = self.getReshapedField(fn, depth_normalize=plopt['dnorm'], normalize=plopt['norm'],
                                                cliprel=plopt['cliprel'], clipabs=plopt['clipabs'], sym=plopt['sym'], 
                                                percent=plopt['percent'], extrema=plopt['extrema'], log=plopt['log'])
        logger.info('Adjusted extrema of field %s in slice %s: %f %f',fn, slname, vmin,vmax)
        cf = ax.pcolormesh(xg, yg, val, cmap=plopt['cmap'], vmin=vmin, vmax=vmax, shading=plopt['shade'], zorder=0.0)

        # iso lines
        if plopt['niso'] > 0:
            if plopt['extrema']:
                cslevels = np.linspace(vmin,vmax,plopt['niso']) 
                cs = ax.contour(xg, yg, val, cslevels, linewidths=0.5, alpha=1.0, cmap=plopt['cmapiso'], 
                                vmin=vmin, vmax=vmax, zorder=0.5)
                logger.info("Isoline levels: %s", " ".join(map(str, list(cslevels))))
            else:
                cs = ax.contour(xg, yg, val, plopt['niso'], linewidths=0.5, alpha=1.0, cmap=plopt['cmapiso'], 
                                vmin=vmin, vmax=vmax, zorder=0.5)
            if plopt['conlab']: ax.clabel(cs, inline=True, fontsize=plopt['ticklabelsize']-4, fmt=plopt['isofmt'], zorder=0.5)

        # depth iso lines
        if xzl.size > 1:
            nzl = xzl.shape[0]
            for i in range(nzl):
                ax.plot(xzl[i, :], yzl[i, :], color='grey', linestyle=':', linewidth=0.5, zorder=5.0)
            # annotate depths
            self.annotateDepth(ax, nzl, plopt['dzline'], xzl[:,0], yzl[:,0], xzl[:,-1], yzl[:,-1], 5, plopt['ticklabelsize'])

        # top and bottom boundary of slice
        for i in range(2):
            ax.plot(xtb[i, :], ytb[i, :], color='black', linestyle='-', linewidth=0.5, zorder=5.0)

        # radial iso distance lines
        if xdl.size > 1:
            ndl = xdl.shape[0]
            for i in range(ndl):
                ax.plot(xdl[i,:], ydl[i,:], color='grey', linestyle=':', linewidth=0.5, zorder=5.0)
            self.annotateDistance(ax, ndl, plopt['ddist'], xdl[:,0], ydl[:,0], 5, plopt['ticklabelsize'])

        # left and right boundary of slice
        for i in range(2):
            ax.plot(xlr[i, :], ylr[i, :], color='black', linestyle='-', linewidth=0.5, zorder=5.0)

        # title
        ax.text(np.min(xedge), np.max(yedge)+1.0*self.dr , slname, fontsize=plopt['ticklabelsize']+4, ha='left', va='bottom')

        # color bar and save
        if plopt['colbar']:
            cb = fig.colorbar(cf, pad=plopt['cbpad'], format=plopt['isofmt'], fraction=0.10, shrink=plopt['cbshrink'])
            cb.ax.tick_params(labelsize=plopt['ticklabelsize']-2)

        # plot individual points if desired
        # dictionary has keys: delta-r, size, marker, color
        if 'points' in kwargs:
            dp = kwargs['points']
            if len(dp) > 0: self.plotPoints(dp, ax)

        # plot individual lines if desired
        if 'lines' in kwargs:
            dl = kwargs['lines']
            if len(dl) > 0: self.plotLines(dl, ax)

        # overlay additional slices if desired
        if 'overlays' in kwargs:
            dsl = kwargs['overlays']
            if len(dsl) > 0: self.plotOverlays(dsl, ax, xg, yg)

        # plot text if desired
        if 'text' in kwargs:
            dtxt = kwargs['text']
            if len(dtxt) > 0: self.plotText(dtxt, ax)

        # add latitude or longitude marks at top boundary if gc-data are available
        if 'gc' in kwargs:
            self.annotateLatLon(ax, kwargs['gc'], plopt['dll'], plopt['dllfmt'], plopt['ticklabelsize'])

        fig.savefig(inbase + '_vslice_' + slname + '_' + fn + '.png', dpi=300, bbox_inches='tight')
        return fig


    def plotHslice(self, inbase, slname, fn, borders, coastline, proj, projdata, **kwargs):
        """
        Plot a horiontal slice using pcolormesh
        :param inbase:basename of file to which png-output is written
        :param slname: name of slice
        :param fn: name of data field
        :param borders: border object from cartopy
        :param coastline: coastline object from cartopy
        :param proj: projection from sphere to plane
        :param projdata: coordinate system in which the field data are defined
        :param kwargs:  implemented keys: (all values are strings)
            :key = gcdata:   value = list of regular spherical grid objects defining great circles
        """
        # first set defaults for plotting options (use deep copy to avoid changing defaults)
        plopt = copy.deepcopy(plot_defaults)

        # overwrite plopt with values from kwargs
        for kw in plopt:
            if kw in kwargs:
                plopt[kw] = kwargs[kw]

        font = {'family': 'serif', 'color': 'black', 'weight': 'normal', 'size': plopt['ticklabelsize']+4.0}
        fig = plt.figure(figsize=(plopt['wx'], plopt['wy']))
        ax = plt.axes(projection=proj)

        #  borders and coastline
        ax.add_feature(borders, zorder=10, alpha=0.5)
        ax.add_feature(coastline, zorder=10, alpha=0.5)

        # extent, title and gridlines
        extent = self.getExtentHS()
        ax.set_extent(extent)
        ax.set_title(slname.replace('_',' '), fontdict=font, pad=0.5*plopt['ticklabelsize']+2)

        # annotated grid lines
        ax.gridlines(draw_labels=True, linewidth=1.0, color='gray', alpha=0.5, linestyle='-',
                     x_inline=False, y_inline=False, rotate_labels=False, xlabel_style={'fontsize':plopt['ticklabelsize']},
                     ylabel_style={'fontsize':plopt['ticklabelsize']}, zorder=5.0)

        # unlabelled grid lines
        xbl, ybl = ccrs.PlateCarree().transform_point(extent[0], extent[2], projdata)        # bottom left
        xtr, ytr = ccrs.PlateCarree().transform_point(extent[1], extent[3], projdata)        # top right
        xbr, ybr = ccrs.PlateCarree().transform_point(extent[1], extent[2], projdata)        # bottom right
        xtl, ytl = ccrs.PlateCarree().transform_point(extent[0], extent[3], projdata)        # top left
        lonmin = min(xbl, xtl)
        lonmax = max(xtr, xbr)
        latmin = min(ybl, ybr)
        latmax = max(ytr, ytl)
        xlocs = getTicks(lonmin,lonmax,plopt['dll'])
        ylocs = getTicks(latmin,latmax,plopt['dll'])
        ax.gridlines(draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='--', zorder=4.9,
                     xlocs=xlocs, ylocs=ylocs)

        # field values
        val,vmin,vmax = self.getReshapedField(fn, normalize=plopt['norm'], cliprel=plopt['cliprel'], clipabs=plopt['clipabs'],
                                              sym=plopt['sym'], percent=plopt['percent'], extrema=plopt['extrema'], 
                                              log=plopt['log'])
        logger.info('Adjusted extrema of field %s in slice %s: %f %f ',fn, slname, vmin, vmax)

        # color mesh
        long, latg = self.getMeshGridHS()
        cf = ax.pcolormesh(long, latg, val, cmap=plopt['cmap'], vmin=vmin, vmax=vmax, shading=plopt['shade'], 
                           transform=projdata, zorder=0.0)

        # color bar
        if plopt['colbar']:
            cb = fig.colorbar(cf, shrink=plopt['cbshrink'], pad=plopt['cbpad'], format=plopt['isofmt'])
            cb.ax.tick_params(labelsize=plopt['ticklabelsize']-2)

        # iso lines
        if plopt['niso'] > 0:
            if plopt['extrema'] != 'None':
                cslevels = np.linspace(vmin,vmax,plopt['niso'])
                cs = ax.contour(long, latg, val, cslevels, vmin=vmin, vmax=vmax, cmap=plopt['cmapiso'], 
                                linewidths=0.5, alpha=1.0, transform=projdata, zorder=0.5)
            else:
                cs = ax.contour(long, latg, val, plopt['niso'], vmin=vmin, vmax=vmax, cmap=plopt['cmapiso'], linewidths=0.5,
                                alpha=1.0, transform=projdata, zorder=0.5)
            if plopt['conlab']: ax.clabel(cs, inline=True, fontsize=plopt['ticklabelsize']-4.0, fmt=plopt['isofmt'], zorder=0.5)

        #  plot great circles of vertical profiles if desired (geographical coordinates)
        if 'gcdata' in kwargs:
            for gc in kwargs['gcdata']:
                ax.plot(np.rad2deg(gc.field[1,:]), np.rad2deg(gc.field[0,:]), color = 'yellow',
                        linestyle = 'solid', transform=ccrs.PlateCarree(), zorder=5.0)
                ax.text(np.rad2deg(gc.field[1,0]), np.rad2deg(gc.field[0,0]), gc.rsgname, zorder=5.0,
                        fontsize=plopt['ticklabelsize'], ha='center', rotation=0.0, va='bottom', transform=ccrs.PlateCarree())

        # plot individual points if desired
        # dictionary has keys: lonlat, size, marker, color
        if 'points' in kwargs:
            dp = kwargs['points']
            if len(dp) > 0: self.plotPoints(dp, ax)

        # plot individual lines if desired
        if 'lines' in kwargs:
            dl = kwargs['lines']
            if len(dl) > 0: self.plotLines(dl, ax)

        # overlay additional slices if desired
        if 'overlays' in kwargs:
            dov = kwargs['overlays']
            if len(dov) > 0: self.plotOverlays(dov, ax, long, latg)

        # plot text if desired
        if 'text' in kwargs:
            dtxt = kwargs['text']
            if len(dtxt) > 0: self.plotText(dtxt, ax)

        # save
        fig.savefig(inbase + '_hslice_' + slname + '_' + fn + '.png', dpi=300, bbox_inches='tight')

        return fig


    def annotateDepth(self, ax, nzl, dzline, xl, yl, xr, yr, nmtval, ticklabelsize):
        """
        Annotate left and right border of vslice with depth values
        :param ax: axis object of vertical slice
        :param nzl:    number of depth lines
        :param dzline: distance between depth lines (m)
        :param xl:     array of x-values of depth lines at left boundary
        :param yl:     array of y-values of depth lines at left boundary
        :param xr:     array of x-values of depth lines at right boundary
        :param yr:     array of y-values of depth lines at right boundary
        :param nmtval: number of minor tick interval per dzline
        :param ticklabelsize: font size of tick labels
        """
        ds = 0.14*(xr[0]-xl[0])
        phi = np.arctan2(xr[nzl//2+1], yr[nzl//2+1])
        dsx = ds*np.cos(phi)
        dsy = ds*np.sin(phi)
        ax.text(xl[nzl//2+1]-dsx, yl[nzl//2+1]-dsy, 'Depth [km]', fontsize=ticklabelsize+2, rotation=90+np.rad2deg(phi),
                ha='left', va='bottom')
        for i in range(nzl):
            phi = np.arctan2(xr[i], yr[i])
            dsx = 0.15*ds*np.cos(phi)
            dsy = 0.15*ds*np.sin(phi)
            ax.text(xl[i]-dsx, yl[i]-dsy, "{:3.0f}".format((i+1)*dzline*1.0e-3), fontsize=ticklabelsize, ha='right',
                    rotation=+np.rad2deg(phi), va='center_baseline')
            ax.text(xr[i]+dsx, yr[i]-dsy, "{:3.0f}".format((i+1)*dzline*1.0e-3), fontsize=ticklabelsize, ha='left',
                    rotation=-np.rad2deg(phi), va='center_baseline')
            ax.plot(np.array([xl[i],xl[i]-0.8*dsx]), np.array([yl[i],yl[i]-0.8*dsy]), 'k-', linewidth=0.5)
            ax.plot(np.array([xr[i],xr[i]+0.8*dsx]), np.array([yr[i],yr[i]-0.8*dsy]), 'k-', linewidth=0.5)
        # minor tickmarks
        dr = dzline/nmtval
        rsurf = 6371000.0
        nmt = int(np.floor((rsurf-self.rmin)/dr))
        for phimt in [np.arctan2(xl[0], yl[0]), np.arctan2(xr[0], yr[0])]:
            for j in range(nmt):
                rmt = rsurf-(j+1)*dr
                xmt = rmt*np.sin(phimt)
                ymt = rmt*np.cos(phimt)
                dsx = 0.10*ds*np.cos(phimt)*np.sign(phimt)
                dsy = 0.10*ds*np.sin(phimt)
                ax.plot(np.array([xmt, xmt+0.8*dsx]), np.array([ymt, ymt-0.8*abs(dsy)]), 'k-', linewidth=0.25)


    def annotateDistance(self, ax, ndl, ddis, xb, yb, nmtval, ticklabelsize):
        """
        Annotate bottom of vslice with distance values
        :param ax: axis object of vertical slice
        :param ndl:    number of distance lines
        :param ddis:   interval between distance lines (m)
        :param xb:     array of x-values of distance lines at bottom of v-slice
        :param yb:     array of y-values of distance lines at bottom of v-slice
        :param nmtval: number of minor tick intervals per ddis
        :param ticklabelsize: font size of tick labels
        """
        phirange = self.dlat*(self.nlat-1)
        dr = -0.075*self.rmin*phirange
        ax.text(0.0, self.rmin+dr, 'Distance [km]', fontsize=ticklabelsize+2, ha='center')
        for i in range(ndl):
            phi = np.arctan2(xb[i], yb[i])
            dsx = 0.20*dr*np.sin(phi)
            dsy = 0.20*dr*np.cos(phi)
            ax.text(xb[i]+dsx, yb[i]+dsy, "{:4.0f}".format((i+1)*ddis*1.0e-3), fontsize=ticklabelsize, ha='center',
                    rotation=-np.rad2deg(phi), va='top')
            ax.plot(np.array([xb[i],xb[i]+0.8*dsx]), np.array([yb[i],yb[i]+0.8*dsy]), 'k-', linewidth=0.5)
        # minor tickmarks
        dphi = ddis/self.rmin/nmtval
        nmt = int(np.floor(phirange/dphi))+1
        for j in range(nmt):
            phimt = -0.5*phirange+j*dphi
            xmt = self.rmin*np.sin(phimt)
            ymt = self.rmin*np.cos(phimt)
            dsx = 0.10*dr*np.sin(phimt)
            dsy = 0.10*dr*np.cos(phimt)
            ax.plot(np.array([xmt, xmt+0.8*dsx]), np.array([ymt, ymt+0.8*dsy]), 'k-', linewidth=0.25)


    def annotateLatLon(self, ax, gc, dll, dllfmt, ticklabelsize):
        """
        Annotate latitude or longitude at top boundary of vslice a values which are a multiple of dll.
        :param ax: axis object of vertical slice
        :param gc:     Great cirle regular grid object
        :param dll:    spacing of lat/lon annotation in degrees
        :param dllfmt: format of lat/lon annotation in degrees
        :param ticklabelsize: font size of tick labels
        """
        # mark latitude or longitude
        marklat = abs(gc.field[0,0]-gc.field[0,-1]) > abs(gc.field[1,0]-gc.field[1,-1])
        if marklat:
            vals = np.rad2deg(gc.field[0,:])
        else:
            vals = np.rad2deg(gc.field[1,:])

        # determine values of lat/lon annotations
        ticks = getTicks(vals[0],vals[-1],dll)

        delmax = gc.dlat*(gc.nlat-1)
        delmid = 0.5*delmax
        h = (self.nr-1)*self.dr
        rmax = self.rmin+h
        dr = 0.075*self.rmin*delmax

        #  Axis label
        if marklat:
            ax.text(0.0, rmax+0.70*dr, 'Latitude (°)', fontsize=ticklabelsize+2, ha='center', va='bottom')
        else:
            ax.text(0.0, rmax+0.70*dr, 'Longitude (°)', fontsize=ticklabelsize+2, ha='center', va='bottom')

        #  add annotations at marks by searching for appropriate delta along profile
        for value in ticks:
            idx_greater = np.nonzero(vals > value)[0]         # indices of elements in vals greater than value
            if vals[0] > vals[-1]:
                if len(idx_greater) == 0:                     # value at left end of gc => no vals > value
                    idx = 0
                else:
                    idx = min(idx_greater[-1], gc.npts-2)    # index of element in vals left of value on gc
            else:
                if len(idx_greater) == 0:                     # value at right end of gc => no vals > value
                    idx = gc.npts-2
                else:
                    idx = idx_greater[0]-1                    # index of element in vals left of value on gc
            delta = (idx+(value-vals[idx])/(vals[idx+1]-vals[idx]))*gc.dlat   # interpolate delta

            # place lat/lon mark and plot tick
            dsx = 0.20*dr*np.sin(delta-delmid)
            dsy = 0.20*dr*np.cos(delta-delmid)
            xb = rmax*np.sin(delta-delmid)
            yb = rmax*np.cos(delta-delmid)
            ax.text(xb+dsx, yb+dsy, dllfmt.format(value), fontsize=ticklabelsize, ha='center',
                    rotation=-np.rad2deg(delta-delmid), va='bottom')
            ax.plot(np.array([xb,xb+0.8*dsx]), np.array([yb,yb+0.8*dsy]), 'k-', linewidth=0.5)


