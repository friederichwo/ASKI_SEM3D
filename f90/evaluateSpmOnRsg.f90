!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!  Evaluate a 3D earth model defined on Specfem's spherical pseudo mesh
!  (centered on the equator of its pole) on any other regular spherical grid.
!  Of course, this makes only sense if the two grids overlap.
!  Values at points of no overlap are set to zero.
!
program evaluateSpmOnRsg
   use mathConstants
   use regularSphericalGrid
   use argumentParser
   use hdfWrapper
   use errorMessage
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (regular_spherical_grid) :: spm,rsg
   double precision, dimension(:), allocatable :: field
   double precision :: clat,clon,hwlon,hwlat,rmax,rmin,dlon,dlat,dr,latmin,lonmin
   double precision :: r,lat,lon,latspm,lonspm
   integer :: nr,nlat,nlon,ip
   character(len=max_length_string) :: spmfile,rsgfile
   character(len=16) :: myname = 'evaluateSpmOnRsg'
!
   call new(errmsg,myname)
!
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
!  command line processing
!
   call init(ap,myname,"Evaluate a model defined on Specfem's spherical pseudo mesh on any other regular spherical grid")
   call addPosarg(ap,'spmfile','sval','HDF file with model on Specfems spherical pseudo mesh')
   call addPosarg(ap,'rsgfile','sval','HDF output file with model on new grid')
   call addOption(ap,'-clat',.true.,'latitude of grid center (deg)','dval','46.2')
   call addOption(ap,'-clon',.true.,'longitude of grid center (deg)','dval','10.87')
   call addOption(ap,'-hwlat',.true.,'latitude half width of new grid around grid center (deg)','dval','6.0')
   call addOption(ap,'-hwlon',.true.,'longitude half width of new grid around grid center (deg)','dval','11.0')
   call addOption(ap,'-rmax',.true.,'top radius of new grid (km)','dval','6371')
   call addOption(ap,'-rmin',.true.,'bottom radius of new grid (km)','dval','5771')
   call addOption(ap,'-dlat',.true.,'latitude spacing of new grid (deg)','dval','0.1')
   call addOption(ap,'-dlon',.true.,'longitude spacing of new grid (deg)','dval','0.1')
   call addOption(ap,'-dr',.true.,'depth spacing (km)','dval','10')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
! -----------------------------------------------------------------------------
   spmfile = ap.sval.'spmfile'
   rsgfile = ap.sval.'rsgfile'
   clat = (ap.dval.'-clat')*mc_deg2radd
   clon = (ap.dval.'-hwlon')*mc_deg2radd
   hwlat = (ap.dval.'-hwlat')*mc_deg2radd
   hwlon = (ap.dval.'-hwlon')*mc_deg2radd
   rmax = (ap.dval.'-rmax')*1.d3
   rmin = (ap.dval.'-rmin')*1.d3
   dlat = (ap.dval.'-dlat')*mc_deg2radd
   dlon = (ap.dval.'-dlon')*mc_deg2radd
   dr = (ap.dval.'-dr')*1.d3
   nr = ceiling((rmax-rmin)/dr)+1
   nlat = 2.*ceiling(hwlat/dlat)+1
   nlon = 2.*ceiling(hwlon/dlon)+1
   write(6,'(a,3i8)') 'Dimensions of new grid: ',nlon,nlat,nr
   latmin = clat-0.5*(nlat-1)*dlat
   lonmin = clon-0.5*(nlon-1)*dlon
!
!  read spherical pseudo mesh
!
   call readHDFRegularSphericalGrid(spm,spmfile,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  create new regular spherical grid
!
   call rsg%create(spm%nfield,0.5*mc_pid,mc_pid,nr,nlat,nlon,dr,dlat,dlon,&
                   rmin,latmin,lonmin,6371000.0d0,spm%field_name_string,'vgrid')
!
!  walk through new grid and get values of model there by interpolation
!
   allocate(field(spm%nfield))
   do while (rsg%nextPoint(ip,r,lat,lon))
      call spm%transformPointToPole(lat,lon,latspm,lonspm)
      call spm%getFieldsAtPoint(r,latspm,lonspm,field)
      rsg%field(ip,:) = real(field)
   enddo
   call rsg%writeHDF(rsgfile,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  clean up
!
   deallocate(field)
   call closeEnvironmentHDFWrapper(errmsg)
   call spm%dealloc()
   call rsg%dealloc()
   call dealloc(errmsg)
   call dealloc(ap)
!------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif
end program evaluateSpmOnRsg

