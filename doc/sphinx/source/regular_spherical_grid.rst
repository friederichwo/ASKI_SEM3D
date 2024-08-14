=================================
 Regular spherical grid
=================================

Class definition
================
.. py:module:: askipy.regularSphericalGrid
.. py:class:: RegularSphericalGrid(dims, spacings, minvals, pole, rsgname, field_name_string, field_values)

    Class specifying instances of a regular spherical grid (RSG) defined in 3D
    by its dimensions, spacings  and minimum values in longitude, latitude and radius, respectively.
    Coordinates may refer to a pole located anywhere on the sphere.
    Several data fields may be defined on the grid identified by their names.

    The main capabilities of the class methods are the extraction of horizontal or vertical slices
    from the 3D grid as regular spherical grid objects, and the plotting of these slices by
    shading and contour plots with optional further external decorations such as text,
    collections of lines, point sets and slice overlays.

    :param dims: dimensions of the grid in longitude, latitude and radial direction,
        ``nlon, nlat, nr``
    :param spacings: grid spacings in longitude, latitude and radius,
        ``dlon, dlat, dr`` (units: rad, rad, m)
    :param minvals: minimum values of longitude, latitude and radius,
        ``lonmin, latmin, rmin`` (units: rad, rad, m)
    :param pole: tuple with geographical coordinates of pole the grid coordinates
        refer to, ``lon_pole, lat_pole`` (units: rad, rad)
    :param rsgname: name of the grid
    :param field_name_string: a colon-seprated string specifying the names of the data fields
        defined on the grid (e. g. ``reldev:absdev:colsum``)
    :param field_values: a 2D array specifying the values of the data fields at the grid points.
        The array is flattened with respect to the grid points with longitude running first,
        latitude second and radius last.

Lower dimensional grids
-----------------------

    Special lower-dimensional cases of a regular spherical grid are:

    :Vertical slice (VS):
        Here nlat > 1, nlon = 1, nr > 1, dlon = 0 and lonmin = azimuth of slice.
        The "latitude" values are abused as angular distance along slice.
    :Horizontal slice (HS):
        Here nlat > 1, nlon > 1, nr = 1, dr = 0 and rmin = radius of slice.
    :Great circle (GC):
        Here nlat > 1, nlon = 1, nr = 1, dlon = 0, dr = 0, lonmin = azimuth of great circle and
        rmin = radius of the spherical shell.
        The "latitude" values are abused as angular distance along the great circle.
        With latitude and longitude of the points as field data, a great circle can be specified.
    :Radial line (RL):
        Here nlat = 1, nlon = 1, nr > 1, dlat = 0, dlon = 0 and  latmin and lonmin define
        lateral position of the line.
        With radius as field data, it may describe a curve in a vertical plane.
        With other field data, it may represent a vertical profile of some quantity.

Attributes
==========

    :plot_defaults: The default values for all plotting options. This is a class variable common
        to all class objects

        :type: dictionary
        :value:
            ``{'dzline': 100000, 'ddist': 100000, 'dll': 1.0, 'dllfmt': '{:2.0f}', 'shade': 'gouraud',
            'extrema': None, 'cliprel': -1.0, 'clipabs': -1.0, 'dnorm': False, 'norm': False,
            'log': False, 'sym': False, 'ticklabelsize': 14, 'niso': -1, 'conlab': False,
            'percent': False, 'wx': 16, 'wy': 8, 'isofmt': '%.2f', 'cmap': 'coolwarm_r',
            'cmapiso': 'seismic_r', 'colbar': False, 'cbpad': 0.05, 'cbshrink': 0.75}``

    :nlon: Number of grid points along longitude
    :nlat: Number of grid points along latitude
    :nr: Number of grid points along radius
    :dlon: Spacing of grid points along longitude
    :dlat: Spacing of grid points along latitude
    :dr: Spacing of grid points along radius
    :lonmin: Number of grid points along latitude
    :latmin: Number of grid points along latitude
    :rmin: Number of grid points along radius
    :lon_pole: Longitude of pole grid coordinates refer to
    :lat_pole: Latitude of pole grid coordinates refer to
    :rsgname: Name of grid
    :field_name_string: Colon-seprated string with names of data fields
    :field: 2D-array of dimension ``nfield, npts`` containing data defined on grid. The array is flattened
        with repsect to the grid points in the order longitude, latitude and radius.
    :nfield: Number of data fields
    :npts: Total number of grid points
    :typ: Type of grid. Allowed values are ``(3D, VS, HS, GC, RL)``
    :dims: Tuple with dimensions
    :spacings: Tuple with spacings
    :minvals: Tuple with minimum values
    :pole: Tuple with coordinates of pole

Methods
=======

Here, the methods are listed that may be called with benefit by the user of the class.

.. py:classmethod:: readFromHDF(cls, filename)

    Class method as second constructor to initialize an instance of the grid from a HDF file

    :param cls: a reference to the class itself
    :param str filename: Full path name of the HDF file
    :return: An instance of the class

.. py:method:: extractVerticalSlice(self, lat1d, lon1d, lat2d, lon2d, nr, ndel, slname)

    Method to extract a vertical slice from a 3D regular spherical grid.
    The slice is taken along a great circle specified by its end points.

    :param lat1d: Latitude of starting point of great circle (deg)
    :param lon1d: Longitude of starting point of great circle (deg)
    :param lat2d: Latitude of end point of great circle (deg)
    :param lon2d: Longitude of end point of great circle (deg)
    :param nr: Number of samples in new slice along radius
    :param ndel: Number of samples in new slice along great circle
    :param slname: Name of new slice
    :return: a RSG-object containing the vertical slice and one containing the great circle

    .. warning::
        Coordinates of the great circle end points refer to the true North Pole and not to the pole the grid
        coordinates refer to.

.. py:method:: extractHorizontalSlice(self, r, nlat, nlon, slname)

    Extract a horizontal slice from a 3D grid at given radius.

    :param r: Radius of silice (m)
    :param nlat: Number of points along latitude
    :param nlon: Number of points along longitude
    :param slname: Name of slice
    :return: a RSG object containing the horizontal slice

.. py:method:: plotVslice(self, inbase, slname, fn, **kwargs)

    Plot a vertical slice

    :param inbase: basename of output png-files
    :param slname: name of vertical slice
    :param fn:     name of data field
    :param kwargs:  the *plotting* dictionary specifying options and decorations
        (see :ref:`sec_plotting`)

.. py:method:: plotHslice(self, inbase, slname, fn, borders, coastline, proj, projdata, **kwargs)

    Plot a horizontal slice

    :param inbase: basename of files to which png-output is written
    :param slname: name of slice
    :param fn: name of data field
    :param borders: border object from cartopy
    :param coastline: coastline object from cartopy
    :param proj: projection from sphere to plane (see cartopy)
    :param projdata: coordinate system in which the field data are defined (see cartopy)
    :param kwargs:  the *plotting* dictionary specifying options and decorations
        (see :ref:`sec_plotting`)

.. _sec_plotting:

Plotting of a slice
==========================================
Plot options
------------

    The class offers the methods ``plotVslice`` and ``plotHslice`` to plot vertical and horizontal
    slices through the 3D grid. The appearance of the plot is controlled by plotting options
    passed to these functions through a dictionary called the *plotting* dictionary in the following.
    The plotting functions consider the following keys and their values:

    :key = dzline: the spacing of constant depth lines in vertical slices (unit: m, default: 100000)
    :key = ddist: the spacing of constant distance lines in vertical slices (unit: m, default: 100000)
    :key = dll: the spacing of latitude or longitude tickmarks at the top boundary of vertical slices,
                or the spacing of unlabelled grid lines in horizontal slices (unit: deg, default: 1.0)
    :key = dllfmt: the format used to write latitude or longitude tickmarks in vertical slices,
                   (default: ``{:2.0f}``)
    :key = shade:    shading mode ('nearest' or 'gouraud', default: 'gouraud')
    :key = extrema:  fixed extrema for plotting slice array values (default = None)
    :key = cliprel:  relative clipping level as fraction of absolute maximum (default=-1.0)
    :key = clipabs:  absolute clipping level (default=-1.0)
    :key = norm:     normalize array values to absolute maximum (default: False)
    :key = log:      take base-10 logarithm of array values (default = False)
    :key = sym:      symmetrize array values (default: False)
    :key = ticklabelsize: size of font for tick labels (default = 14)
    :key = niso:     number of isolines (default: -1)
    :key = conlab:   label contour lines (default: False)
    :key = percent:  convert array values to percent by multiplying with 100 (default: False)
    :key = wx:       width of figure window (default = 16)
    :key = wy:       height of figure window (default = 8)
    :key = cmap:     name of color map for shading plot (default = coolwarm_r)
    :key = cmapiso:  name of color map for isolines (default = seismic_r)
    :key = isofmt:   format for iso lines (default = '%.2f')
    :key = colbar:   add color bar (default = True)
    :key = cbpad:    padding space between graph and colorbar (default = 0.05)
    :key = cbshrink: height of colorbar as fraction of axes length (default = 0.75)

External decorations
--------------------

    External decorations such as point clouds, lines, overlay slices and text may be added to the current
    slice. Specfications are passed to the plotting functions ``plotVslice`` and ``plotHslice``
    via the *plotting dictionary*. The following keys in this dictionary are considered by these
    functions and signify the plotting of external decorations.

    :key = points:
        Add point clouds

        :Value:
            A dictionary with string keys ``"1"``, ``"2"`` and so forth.
            Values to these keys are dictionaries with keywords ``coor`` defining the point locations
            as a 2D-array, either specifying ``(lon,lat)`` in a horizontal slice
            or ``(delta, r)`` in a vertical slice, ``size`` defining point size, ``marker`` defining
            marker type and ``color`` defining marker color.

    :key = lines:
        Add line collections

        :Value:
            A dictionary with string keys ``"1"``, ``"2"`` and so forth. Values to these keys are
            dictionaries with integer keys ``1``, ``2`` and so forth. Values to these keys are
            dictionaries with keys ``coor`` defining the locations of the points on the line as a 2D-array
            either specifying ``(lon,lat)`` in a horizontal slice or ``(delta, r)`` in a vertical slice,
            ``lw`` specifying lineweight, ``ls`` specifying line style and ``color`` specifying color.

    :key = overlays:
        Overlay external slices of same dimensions

        :Value:
            A dictionary with with string keys ``"1"``, ``"2"`` and so forth. Values to these keys
            are dictionaries with keys ``array`` specifying the values of the overlay as a 2D array with
            the same dimensions as the current slice, ``colormap`` specifying the colormap to be used for
            the overlay, ``vmin`` providing a minimum value and ``vmax`` providing a maximum value for the
            overlay and ``shading`` setting the shading mode.

    :key = text:
        Add text with optional arrows pointing to the target

        :Value:
            A dictionary with string keys ``"1"``, ``"2"`` and so forth.  Values to these keys
            are dictionaries with keys ``text`` defining the string to be plotted, ``tpos``, a tuple
            of slice coordinates defining the position of the text, ``rpos``, a tuple of slice coordinates
            defining the position of the target point, ``va`` setting the vertical alignment of the text
            (either ``top`` or ``bottom``), ``arrow`` telling if an arrow from text to target is desired,
            ``fontweight`` setting for example ``bold`` type and ``size`` defining text size.

    :key = gc:
        Only for vertical slices.
        If set, the upper boundary of the vertical slice will be annotated either by latitude or longitude
        markers, depending on whether the latitude or longitud range is larger.

        :Value:
            A regular spherical grid object defining the great circle of the vertical slice

    :key = gcdata:
        Only for horizontal slices.
        If set, the specified great circles are plotted into the horizontal slice.
        If both vertical and horzontal slices are plotted, this option may be used to show
        the location of the vertical slices in the horizontal ones.

        :Value:
            A list of regular spherical grid objects defining great circles

Further methods
==============================
These methods are of less benefit for the user as they are intended mainly for internal use.

.. py:method:: setupQuadmesh(self)

    Generate Cartesian coordinates of cell centers and cell edges for a vertical 2D slice.
    See matplotlib docu on Quadmesh for a description

    :return:
        **xg, yg, xedge, yedge**: Cartesian x and y-coordinates of grid points as 2D arrays of
        dimension ``(nr,nlat)`` as well as Cartesian x- and y-coordinates of cell edges
        as 2D arrays of dimension ``(nr+1,nlat+1)``

.. py:method:: computeIsoDepthLines(self,dzline)

    Calculate Cartesian coordinates of iso-depth lines in a vertical slice

    :param dzline: spacing between iso-depth lines in meters
    :return:
        **xdl, ydl**: Cartesian coordinates of iso depth lines as 2D array with dimension ``(nline, nlat)``

.. py:method:: computeTopBottomBoundaryLines(self)

    Calculate Cartesian coordinates of top and bottom boundary of a vertical slices

    :returrn:
        **xtl, ytl**: Cartesian coordinates of top and bottom boundary as 2D arrays with
        dimension ``(2, nlat)``

.. py:method:: computeIsoDistanceLines(self, ddis)

    Calculate coordinates of radial iso distance lines

    :param ddis: spacing between iso distance lines (unit: m)
    :return:
        **xdl, ydl**: Cartesian coordinates of iso distance lines as 2D arrays with
        dimension ``(nline, nlat)``

.. py:method:: computeLeftRightBoundaryLines(self)

    Calculate Cartesian coordinates of left and right boudnary lines of a vertical slice

    :return:
        **xdl, ydl**: Cartesian coordinates of left/right boundaries of a vertical slice as 2D arrays
            with dimension ``(nline, nr)``

.. py:method:: getExtentHS(self)

    Get the extent of a horizontal slice

    :return:
        An array with four elements ``(latmin, latmax, lonmin, lonmax)`` in degrees.

.. py:method:: getMeshGridHS(self)

    :return:
        **long, latg**: latitude and longitude of HS grid points as a meshgrid.

.. py:method:: getFieldNames(self)

    :return:
        a list of the field names associated with the grid

.. py:method:: getReshapedField(self, field_name, depth_normalize=False, normalize=False, cliprel=-1.0, clipabs=-1.0, sym=False, percent=False, extrema=None, log=False)

    Reshapes the flattened field array for the given field name to more than one dimension
    according to type of slice

    :param field_name: name of data field to be processed
    :param depth_normalize: flag for depth normalization of field values
    :param normalize: flag for max-normalization of field data
    :param cliprel: flag and value for clipping of data at a given fraction of the maximum
    :param clipabs: flag and value for clipping of data at agiven value
    :param sym: flag for symmetrising the field data
    :param percent: flag for conversion of field data to percent
    :param extrema: flag and values of min/max-values for color scales and isolines
    :param log: flag for taking base 10 log of values and zero waterlevel
    :return: multi-dimensional array containing modified field data

.. py:method:: plotPoints(self, dp, ax)

    Helper rotine to plot point cloud decorations into a slice

    :param dp: dictionary ordered by index with values as dictionaries having keys ``coor``,
        ``size``, ``marker``, ``color``
    :param ax: axis into which points are plotted
    :return: None

.. py:method:: plotLines(self, dl, ax)

    Plot a series of line collections into a slice

    :param dl: dictionary ordered by index. Value is another dict ordered by index and
        value as dict with keys ``coor``, ``ls``, ``lw``, ``color``
    :param ax: axis into which lines are plotted
    :return: None

.. py:method:: plotOverlays(self, dov, ax, c1g, c2g)

    Plot external overlays of same dimension into slice.

    :param dov: Dictionary ordered by index and values a dict with keys ``array``, ``colormap``, ``vmin``, ``vmax``
    :param ax: axis into which overlay is plotted
    :param c1g: array with values of first grid coordinate
    :param c2g: array with values of second grid coordinate
    :return: None

.. py:method:: plotText(self, dt, ax)

     Plot text into slice with optional arrows connecting the text location with the target point

    :param dt: dictionary ordered by index and value a dict with keys ``text``, ``tpos``, ``rpos``, ``va``, ``arrow``,
       ``fontweight``, ``size``
    :param ax: axis into which text is plotted
    :return: None

.. py:method:: annotateDepth(self, ax, nzl, dzline, xl, yl, xr, yr, nmtval, ticklabelsize)

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
    :return: None

.. py:method:: annotateDistance(self, ax, ndl, ddis, xb, yb, nmtval, ticklabelsize)

    Annotate bottom of vslice with distance values

    :param ax: axis object of vertical slice
    :param ndl:    number of distance lines
    :param ddis:   interval between distance lines (m)
    :param xb:     array of x-values of distance lines at bottom of v-slice
    :param yb:     array of y-values of distance lines at bottom of v-slice
    :param nmtval: number of minor tick intervals per ddis
    :param ticklabelsize: font size of tick labels
    :return: None

.. py:method:: annotateLatLon(self, ax, gc, dll, dllfmt, ticklabelsize)

    Annotate latitude or longitude at top boundary of vslice a values which are a multiple of dll.

    :param ax: axis object of vertical slice
    :param gc:     Great cirle regular grid object
    :param dll:    spacing of lat/lon annotation in degrees
    :param dllfmt: format of lat/lon annotation in degrees
    :param ticklabelsize: font size of tick labels
    :return: None




