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
"""
This module contains selected functions from the Fortran module "axesRotation.f90"

Detailed descriptions on how certain variables are defined and how specific rotations are derived
can bbe found in the documentation of that f90-file.

"""

import numpy as np


def testForTypeOfInput(vector):
    """
    Test what type the members of a given vector are and creates the new vector based on the input.

    :param vector:
    :return:
    """
    number = all(isinstance(x, (float, int)) for x in vector)
    array = all(isinstance(x, np.ndarray) for x in vector)
    if number:
        newVector = np.empty(len(vector))
    elif array:
        newVector = np.asarray([np.empty(len(vector[0])), np.empty(len(vector[0])), np.empty(len(vector[0]))])
    else:
        raise ValueError("unexpected combination of values in the input vector!")

    return newVector


def vectorGCfromLC(thetas, phis, ul):
    """
    Calculate cartesian GC vector components from given LC components

    :param thetas: colatitude (given in Rad)
    :param phis: longitude (given in Rad)
    :param ul: Vector in local cartesian coordinates
    :return: ug: Vector in global cartesian coordinates
    """
    ug = testForTypeOfInput(ul)

    ct = np.cos(thetas)
    st = np.sin(thetas)
    cp = np.cos(phis)
    sp = np.sin(phis)

    ug[0] = ct * cp * ul[0] - sp * ul[1] + st * cp * ul[2]
    ug[1] = ct * sp * ul[0] + cp * ul[1] + st * sp * ul[2]
    ug[2] = -st * ul[0] + ct * ul[2]
    return ug


def vectorLCfromGC(thetas, phis, ug):
    """
    Calculate cartesian LC vector components from given GC components

    :param thetas: colatitude (given in Rad)
    :param phis: longitude (given in Rad)
    :param ug: Vector in global cartesian coordinates
    :return: ul: Vector in local cartesian coordinates
    """

    ul = testForTypeOfInput(ug)
    ct = np.cos(thetas)
    st = np.sin(thetas)
    cp = np.cos(phis)
    sp = np.sin(phis)

    ul[0] = ct * cp * ug[0] + ct * sp * ug[1] - st * ug[2]
    ul[1] = -sp * ug[0] + cp * ug[1]
    ul[2] = st * cp * ug[0] + st * sp * ug[1] + ct * ug[2]
    return ul


def vectorLSfromLC(delta, xi, uc):
    """
    return spherical LS components of a vector from cartesian LC ones

    :param delta: Epicentral distance (given in rad)
    :param xi: azimuth at a point on the sphere counted from S over E (counterclockwise)
    :param uc: vector in the local cartesian coordinate system
    :return: us: vector in the local spherical system.
    """
    us = testForTypeOfInput(uc)
    cd = np.cos(delta)
    sd = np.sin(delta)
    cx = np.cos(xi)
    sx = np.sin(xi)
    us[0] = sd * cx * uc[0] + sd * sx * uc[1] + cd * uc[2]
    us[1] = cd * cx * uc[0] + cd * sx * uc[1] - sd * uc[2]
    us[2] = -sx * uc[0] + cx * uc[1]

    return us


def vectorLCfromLS(delta, xi, us):
    """
    return cartesian LC components of a vector from spherical LS ones
    :param delta: Epicentral distance (given in rad)
    :param xi: azimuth at a point on the sphere counted from S over E (counterclockwise)
    :param us: vector in the local spherical system.
    :return: uc: vector in the local cartesian coordinate system
    """
    uc = testForTypeOfInput(us)
    cd = np.cos(delta)
    sd = np.sin(delta)
    cx = np.cos(xi)
    sx = np.sin(xi)
    uc[0] = sd * cx * us[0] + cd * cx * us[1] - sx * us[2]
    uc[1] = sd * sx * us[0] + cd * sx * us[1] + cx * us[2]
    uc[2] = cd * us[0] - sd * us[1]
    return uc


def vectorLCfromRC(gamma, ur):
    """
    Rotate a vector from a rotated cartesian to a local cartesian coordinate system

    :param gamma: Rotation angle
    :param ur: vector in the rotated cartesian coordinate system
    :return: ul: vector in the local cartesian coordinate system
    """
    ul = vectorGCfromLC(0., gamma, ur)
    return ul


def vectorRCfromLC(gamma, ul):
    """
    Rotate a vector from a local cartesian to a rotated cartesian coordinate system

    :param gamma: Rotation angle
    :param ul: vector in the local cartesian coordinate system
    :return: ur: vector in the rotated cartesian coordinate system
    """
    ur = vectorLCfromGC(0., gamma, ul)
    return ur


def vectorZNEfromRLT(alfa, ur, ul, ut):
    """
    Rotate a vector from an RLT to a ZNE coordinate system

    :param alfa: Angle of rotation
    :param ur: radial component
    :param ul: lateral component
    :param ut: tangential component
    :return: z,n,e components of the vector
    """
    ca = np.cos(alfa)
    sa = np.sin(alfa)
    ue = sa * ul + ca * ut
    un = -ca * ul + sa * ut
    uz = ur
    return uz, un, ue


def vectorRLTfromZNE(alfa, uz, un, ue):
    """
    Rotate a vector from a ZNE to an RLT coordinate system

    :param alfa:
    :param uz: Z-component
    :param un: North-component
    :param ue: East-component
    :return: ur, ul, ut: Radial, Lateral and Tangential components
    """
    ca = np.cos(alfa)
    sa = np.sin(alfa)
    ur = uz
    ul = sa * ue - ca * un
    ut = ca * ue + sa * un
    return ur, ul, ut


def coordinatesLCfromLS(r, delta, xi):
    """
    Calculate local cartesian from given spherical coordinates.

    :param r: Radius
    :param delta: latitude
    :param xi: longitude
    :return: xg, yg, zg: global cartesian coordinates
    """

    cd = np.cos(delta)
    sd = np.sin(delta)
    cx = np.cos(xi)
    sx = np.sin(xi)

    xl = r * sd * cx
    yl = r * sd * sx
    zl = r * cd
    return xl, yl, zl


def coordinatesLSfromLC(x, y, z):
    """
    Calculate local spherical coordinates based on local cartesian coordinates

    :param x: local x-coordinate
    :param y: local y-coordinate
    :param z: local z-coordinate
    :return: r, delta, xi: Radius, epicentral distance, azimuth at a point on the sphere counted from S over E
    """
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    delta = np.arccos(z / r)
    xi = np.arctan2(y, x)
    return r, delta, xi


def coordinatesLCfromGC(thetas, phis, xg, yg, zg):
    """
    Calculate local cartesian coordinates from given global cartesian coordinates
    based on a Point on the sphere P (thetas, phis)

    :param thetas: colatitude (given in Rad)
    :param phis: longitude (given in Rad)
    :param xg: global cartesian x-coordinate
    :param yg: global cartesian y-coordinate
    :param zg: global cartesian z-coordinate
    :return: x,y,z: Local cartesian coordinates
    """

    ct = np.cos(thetas)
    st = np.sin(thetas)
    cp = np.cos(phis)
    sp = np.sin(phis)
    x = ct * cp * xg + ct * sp * yg - st * zg
    y = -sp * xg + cp * yg
    z = st * cp * xg + st * sp * yg + ct * zg
    return x, y, z


def coordinatesGCfromLC(thetas, phis, x, y, z):
    """
    Calculates global cartesian coordinates for given local cartesian coordinates.
    The angle of rotation is based on a given Point on a Sphere (thetas, phis)

    :param thetas: colatitude (given in Rad)
    :param phis: longitude (given in Rad)
    :param x: local x-coordinate
    :param y: local y-coordinate
    :param z: local z-coordinate
    :return: xg,yg,zg: global cartesian coordinates
    """
    ct = np.cos(thetas)
    st = np.sin(thetas)
    cp = np.cos(phis)
    sp = np.sin(phis)

    xg = ct * cp * x - sp * y + st * cp * z
    yg = ct * sp * x + cp * y + st * sp * z
    zg = -st * x + ct * z
    return xg, yg, zg


def coordinatesRCfromLC(gamma, x, y, z):
    """
    return cartesian RC coordinates of point in sphere from given LC coordinates
    :param gamma: Angle of rotation
    :param x: local x-coordinate
    :param y: local y-coordinate
    :param z: local z-coordinate
    :return: xr, yr, zr, rotated local cartesian coordinates
    """

    xr, yr, zr = coordinatesLCfromGC(0., gamma, x, y, z)
    return xr, yr, zr


def coordinatesLCfromRC(gamma, xr, yr, zr):
    """

    :param gamma: Angle of rotation
    :param xr: local rotated x-coordinate
    :param yr: local rotated y-coordinate
    :param zr: local rotated z-coordinate
    :return: x,y,z: local cartesian coordinates
    """

    x, y, z = coordinatesGCfromLC(0, gamma, xr, yr, zr)

    return x, y, z


def calculateBoxCoordinates(r, lat, lon, box_centerLat, box_centerLon):
    """
    Calculate the cartesian coordinates for the station in the specfem box.

    in order to get the true coordinate of the point in the local coorindate system of the specfem box, the user has to
    subtract the Earth' radius from the returned z value. (This might be adjusted in the future.)

    :param: lat, lon, box_centerLat, box_centerLon in degrees
    :return: x,y,z = Cartesian SPECFEM box coordinates
    """

    sta_colat = 0.5 * np.pi - np.deg2rad(lat)
    sta_lon = np.deg2rad(lon)
    xg, yg, zg = coordinatesLCfromLS(r, sta_colat, sta_lon)
    xc, yc, zc = coordinatesLCfromGC(0.5 * np.pi - np.deg2rad(box_centerLat), np.deg2rad(box_centerLon), xg, yg, zg)
    x, y, z = coordinatesRCfromLC(np.pi / 2, xc, yc, zc)

    return x, y, z


def calculateSphericalCoordinates(x, y, z, box_centerLat, box_centerLon):
    """
    Calculate the spherical coordinates of a point given in SpecfemBox coordinates.

    In order to get the correct spherical coordinates on the Earth for a given point on the specfem box
    the user has to add the Earth' radius to the z coordinate. (This might be adjusted in the future)

    If this function is used to calculate the geographical coordinates, the latitude needs to be
    calculated from the colatitude via
        lat = 90 - colat

    :param x: local x-coordinate
    :param y: local y-coordinate
    :param z: local z-coordinate
    :param box_centerLat: Latitude of the box' center (deg)
    :param box_centerLon: Longitude of the box' center (deg)
    :return: r colat lon (colat and lon are returned in degree)
    """

    center_colat = 0.5 * np.pi - np.deg2rad(box_centerLat)
    center_lon = np.deg2rad(box_centerLon)

    xc, yc, zc = coordinatesLCfromRC(np.pi / 2, x, y, z)
    xg, yg, zg = coordinatesGCfromLC(center_colat, center_lon, xc, yc, zc)
    r, colat, lon = coordinatesLSfromLC(xg, yg, zg)

    # convert from rad to deg
    colat = np.rad2deg(colat)
    lon = np.rad2deg(lon)

    return r, colat, lon
