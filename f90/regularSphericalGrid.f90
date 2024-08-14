!----------------------------------------------------------------------------
!   Copyright 2023 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!   Module for a regular spherical grid with variable number of fields
!   defined on it. Field arrays are flattened.
!   Allows 2D spherical grids either as horizontal or vertical slice.
!   Allows 1D spherical grid either as great circle or radial line
!   Vertical slice:    nlon = 1, dlon = 0, nr > 1, nlat > 1, lonmin = azimuth of slice
!   Horizontal slice:  nr = 1, dr = 0, nlat > 1, nlon > 1, rmin = radius of slice
!   Great circle: nr = 1, dr = 0, nlat > 1, nlon = 1, latmin = 0, lonmin: azimuth, rmin: radial shell of great circle
!      With associated data of latitudes and longitudes at points, this gives the description of a great circle
!      With additional values of radius this could be a ray
!   Radial line:  nr > 0, nlat = 1, nlon = 1, first point of line: (rmin,latmin,lonmin)
!
module regularSphericalGrid
!
   use mathConstants
   use errorMessage
   use hdfWrapper
   use axesRotation
   use globalHdfInfo
   implicit none
!
   type :: regular_spherical_grid
      character(len=max_length_string) :: rsgname                    ! name of regular spherical grid
      character(len=2) :: typ                                        ! type of spherical grid (3D, VS, HS, GC, RL)
      integer :: nfield                                              ! number of fields associated with grid
      integer :: npts                                                ! total number of grid points
      integer :: nr,nlat,nlon                                        ! number of grid points along each dimension
      double precision :: lat_pole, lon_pole                         ! coordinates of grid refer to this pole (rad)
      double precision :: dr,dlat,dlon                               ! grid spacing along each dimension (m,rad,rad)
      double precision :: rmin,latmin,lonmin                         ! minimum value for each space coordinate (m,rad,rad)
      double precision :: rearth                                     ! earth radius in m
      real, dimension(:,:), allocatable :: field                     ! field values on grid (ip,field index)
      character(len=max_length_string) :: field_name_string          ! string composed of field names separated by ':'
   contains
      procedure :: create => createRegularSphericalGrid
      procedure :: dealloc => deallocRegularSphericalGrid
      procedure :: nextPoint => nextPointRegularSphericalGrid
      procedure :: getIndicesClosest => getIndicesClosestRegularSphericalGrid
      procedure :: getIndicesBelowSouthWest => getIndicesBelowSouthWestRegularSphericalGrid
      procedure :: getFieldsAtPoint => getFieldsAtPointRegularSphericalGrid
      procedure :: isOffGrid => isOffGridRegularSphericalGrid
      procedure :: transformPointToPole => transformPointToPoleRegularSphericalGrid
      procedure :: transformPointsToPole => transformPointsToPoleRegularSphericalGrid
      procedure :: extractVerticalSlice => extractVerticalSliceRegularSphericalGrid
      procedure :: extractHorizontalSlice => extractHorizontalSliceRegularSphericalGrid
      procedure :: writeHdf => writeHDFRegularSphericalGrid
      procedure :: readHdf => readHDFRegularSphericalGrid
   end type regular_spherical_grid
!
   contains
!---------------------------------------------------------------------------
!  create regular spherical grid
!
   subroutine createRegularSphericalGrid(this,nfield,plat,plon,nr,nlat,nlon,dr,dlat,dlon,&
                                         rmin,latmin,lonmin,rearth,field_name_string,rsgname)
      class (regular_spherical_grid) :: this
      integer :: nfield,nr,nlat,nlon
      double precision :: plat,plon,dr,dlat,dlon,rmin,latmin,lonmin,rearth
      character(len=*) :: field_name_string,rsgname
   !
      this%nfield = nfield
      this%lat_pole = plat
      this%lon_pole = plon
      this%nr = nr; this%nlat = nlat; this%nlon = nlon
      this%dr = dr; this%dlat = dlat; this%dlon = dlon
      this%rmin = rmin; this%latmin = latmin; this%lonmin = lonmin
      this%rearth = rearth
      this%npts = nr*nlat*nlon
      allocate(this%field(this%npts,nfield))
      this%field_name_string = field_name_string
      this%rsgname = rsgname
   !
   !  type of grid):
   !  3D                           (nlon>1, nr>1, nlat>1)
   !  VS: vertical slice           (nlon=1, nr>1, nlat>1)
   !  HS: horizontal slice         (nlon>1, nr=1, nlat>1)
   !  GC: great circle             (nlon=1, nr=1, nlat>1)
   !  RL: radial line              (nlon=1, nr>1, nlat=1)
   !
      this%typ = '3D'
      if (nr == 1 .and. nlon >  1) this%typ = 'HS'
      if (nr == 1 .and. nlon == 1) this%typ = 'GC'
      if (nlat > 1 .and. nlon == 1) this%typ = 'VS'
      if (nlat == 1) this%typ = 'RL'
   !
   end subroutine createRegularSphericalGrid
!----------------------------------------------------------------------------
!  Deallocate regular spherical grid object
!
   subroutine deallocRegularSphericalGrid(this)
      class (regular_spherical_grid) :: this
      if (allocated(this%field)) deallocate(this%field)
   end subroutine deallocRegularSphericalGrid
!----------------------------------------------------------------------------
!  Iterator over grid points
!  returns coordinates and flattened index of point
!
   function nextPointRegularSphericalGrid(this,ip,r,lat,lon) result(next)
      class (regular_spherical_grid) :: this
      integer :: ip
      double precision :: r,lat,lon
      logical :: next
      integer :: i = 0,j = 1,k = 1,call_count = 0
      save :: i,j,k,call_count
   !
      next = .true.
      call_count = call_count+1
      i = i+1
      if (i > this%nlon) then
         i = 1
         j = j+1
         if (j > this%nlat) then
            j = 1
            k = k+1
            if (k > this%nr) then
               next = .false.
               i = 0; j = 1; k = 1; call_count = 0
               return
            end if
         end if
      end if
      ip = call_count
      r = this%rmin+(k-1)*this%dr
      lat = this%latmin+(j-1)*this%dlat
      lon = this%lonmin+(i-1)*this%dlon
   end function nextPointRegularSphericalGrid
!------------------------------------------------------------------------
!  get indices of closest grid point
!  (rk-1)*dr = r-rmin => rk = (r-rmin)/dr+1 => k = nint(rk)
!  (rj-1)*dlat = lat-latmin  => rj = (lat-latmin)/dlat+1 => j = nint(rj)
!  (ri-1)*dlon = lon-lonmin  => ri = (lon-lonmin)/dlon+1 => i = nint(ri)
!
   subroutine getIndicesClosestRegularSphericalGrid(this,r,lat,lon,k,j,i)
      class (regular_spherical_grid) :: this
      double precision :: r,lat,lon
      integer :: k,j,i
      k = nint((r-this%rmin)/this%dr)+1
      j = nint((lat-this%latmin)/this%dlat)+1
      i = nint((lon-this%lonmin)/this%dlon)+1
   end subroutine getIndicesClosestRegularSphericalGrid
!------------------------------------------------------------------------
!  get indices of grid point below, south and west of given point
!  except at west, south and bottom boundaries
!  (k-1)*dr < r-rmin <= k*dr => k-1 < (r-rmin)/dr <= k
!  (j-1)*dlat < lat-latmin <= j*dlat
!  (i-1)*dlon < lon-lonmin <= i*dlon
!
   subroutine getIndicesBelowSouthWestRegularSphericalGrid(this,r,lat,lon,k,j,i)
      class (regular_spherical_grid) :: this
      double precision :: r,lat,lon
      integer :: k,j,i
   !
      k = ceiling((r-this%rmin)/this%dr)
      j = ceiling((lat-this%latmin)/this%dlat)
      i = ceiling((lon-this%lonmin)/this%dlon)
      if (k == 0) k = 1
      if (j == 0) j = 1
      if (i == 0) i = 1    
      if (k == this%nr) k = this%nr-1        ! allows extrapolation towards surface  
   end subroutine getIndicesBelowSouthWestRegularSphericalGrid
!-------------------------------------------------------------------------
!  get field values at given point by 8-point trilinear interpolation
!  r,lat,lon:  coordinates of point (m,rad,rad)
!  field:      (out) dummy array to be allocated by caller with dimension this%nfield
!
   subroutine getFieldsAtPointRegularSphericalGrid(this,r,lat,lon,field)
      class (regular_spherical_grid) :: this
      double precision :: r,lat,lon
      double precision, dimension(:) :: field
      double precision :: rg,latg,long,rn,la,lo
      double precision, dimension(8) :: w
      integer, dimension(8) :: ic
      integer :: i,j,k,l
   !
      field = 0.d0
   !
   !  return a zero value if point is outside grid
   !
      if (this%isOffGrid(r,lat,lon)) return
   !   
      call this%getIndicesBelowSouthWest(r,lat,lon,k,j,i)
   !
   !  normalized coordinates
   !
      long = this%lonmin+(i-1)*this%dlon
      latg = this%latmin+(j-1)*this%dlat
      rg = this%rmin+(k-1)*this%dr
      lo = (lon-long)/this%dlon
      la = (lat-latg)/this%dlat
      rn = (r-rg)/this%dr
   !
   !  interpolation weights
   !  w1+w4 = (1-lo)*(1-la); w2+w5 = lo*(1-la); w3+w6 = (1-lo)*la; w7+w8 = lo*la
   !  w1+w4+w2+w5 = 1-la; w3+w6+w7+w8 = la => sum(w) = 1
   !
      w(1) = (1.0-lo)*(1.0-la)*(1-rn)            ! i      j      k
      w(2) = lo*(1.0-la)*(1.0-rn)                ! i+1    j      k
      w(3) = (1.0-lo)*la*(1.0-rn)                ! i      j+1    k
      w(4) = (1.0-lo)*(1.0-la)*rn                ! i      j      k+1
      w(5) = lo*(1.0-la)*rn                      ! i+1    j      k+1
      w(6) = (1.0-lo)*la*rn                      ! i      j+1    k+1
      w(7) = lo*la*(1.0-rn)                      ! i+1    j+1    k
      w(8) = lo*la*rn                            ! i+1    j+1    k+1
   !   write(6,'(3i6,8f12.3)') k,j,i,(w(l),l=1,8) 
   !
   !  flattened indices of corner points
   !
      ic(1) = i+this%nlon*(j-1+(k-1)*this%nlat)
      ic(2) = ic(1)+1
      ic(3) = ic(1)+this%nlon
      ic(4) = ic(1)+this%nlat*this%nlon
      ic(5) = ic(4)+1
      ic(6) = ic(4)+this%nlon
      ic(7) = ic(3)+1
      ic(8) = ic(6)+1
      if (any(ic > this%npts)) then
         write(6,'(3i8,3f15.1)') i,j,k,long,latg,rg
         write(6,'(8i8)') ic
      endif
   !
   !  special treatment of 2D and 1D objects:
   !  horizontal slice         (nlon>1, nr=1, nlat>1)
   !  great circle             (nlon=1, nr=1, nlat>1)
   !  vertical slice           (nlon=1, nr>1, nlat>1)
   !  radial line              (nlon=1, nr>1, nlat=1)
   !
      if (this%typ == 'HS') then
         field = field + w(1)*this%field(ic(1),:) + w(2)*this%field(ic(2),:) &
                       + w(3)*this%field(ic(3),:) + w(7)*this%field(ic(7),:)
      else if (this%typ == 'VS') then
         field = field + w(1)*this%field(ic(1),:) + w(3)*this%field(ic(3),:) &
                       + w(4)*this%field(ic(4),:) + w(6)*this%field(ic(6),:)
      else if (this%typ == 'GC') then
         field = field + w(1)*this%field(ic(1),:) + w(3)*this%field(ic(3),:)
      else if (this%typ == 'RL') then
         field = field + w(1)*this%field(ic(1),:) + w(4)*this%field(ic(4),:)
      else
         do l = 1,8
            field = field + w(l)*this%field(ic(l),:)
         end do
      end if
   end subroutine getFieldsAtPointRegularSphericalGrid
!------------------------------------------------------------------------
!  is a point not in the grid ?
!  r,lat,lon:   coordinates of point (m,rad,rad)
!
   function isOffGridRegularSphericalGrid(this,r,lat,lon) result(res)
      class (regular_spherical_grid) :: this
      double precision :: r,lat,lon,rmax
      logical :: res
   !
      res = .false.
      if (lat < this%latmin .or. lat > this%latmin+(this%nlat-1)*this%dlat) then
         res = .true.
         return
      end if
      if (lon < this%lonmin .or. lon > this%lonmin+(this%nlon-1)*this%dlon) then
         res = .true.
         return
      end if
      rmax = this%rmin+(this%nr-1)*this%dr
      if (r < this%rmin .or. r > rmax) then
         res = .true.
         if (abs(r-rmax)/this%dr <= 0.5001) res = .false.    ! avoid setting values at surface to zero
         return
      end if
   end function isOffGridRegularSphericalGrid
!-----------------------------------------------------------------------------
!  transform geographical coordinates of a single point on the sphere to
!  spherical coordinates relative to the pole of the grid
!  lat, lon:    geographcal latitude and longitude of point in rad
!  latg, long:  latitude and longitude of point in rad relative to grid pole
!
   subroutine transformPointToPoleRegularSphericalGrid(this,lat,lon,latg,long)
      class (regular_spherical_grid) :: this
      double precision :: lat,lon
      double precision :: latg,long
      double precision :: thetap,r,theta,xi,xg,yg,zg,x,y,z
   !
      thetap = 0.5d0*mc_pid-this%lat_pole
      theta = 0.5d0*mc_pid-lat
      call coordinatesLCfromLSAxesRotation(1.d0,theta,lon,xg,yg,zg)             ! global Cartesian coordinates
      call coordinatesLCfromGCAxesRotation(thetap,this%lon_pole,xg,yg,zg,x,y,z)    ! x,y,z relative to pole
      call coordinatesLSfromLCAxesRotation(x,y,z,r,theta,xi)
      if (cos(this%lon_pole-lon) < 0.0d0) then                 ! pole is on opposite hemisphere, subtract pi
         long = xi-mc_pid
         if (long < -mc_pid) long = long+mc_two_pid
      else
         long = xi
      endif
      latg = 0.5d0*mc_pid-theta
   end subroutine transformPointToPoleRegularSphericalGrid
!-----------------------------------------------------------------------------
!  transform geographical coordinates of points on the sphere to
!  geographical coordinates relative to pole of grid
!  lat, lon:    arrays with geographical latitude and longitude of points in rad
!  latg, long:  arrays with geographical latitude and longitude of points in rad relative to grid pole
!
   subroutine transformPointsToPoleRegularSphericalGrid(this,lat,lon,latg,long)
      class (regular_spherical_grid) :: this
      double precision, dimension(:) :: lat,lon
      double precision, dimension(:), allocatable :: latg,long
      integer :: j
   !
      allocate(latg(size(lat)), long(size(lat)))
      do j = 1,size(lat)
         call this%transformPointToPole(lat(j),lon(j),latg(j),long(j))
      end do
   end subroutine transformPointsToPoleRegularSphericalGrid
!-----------------------------------------------------------------------------
!  extract a horizontal slice at constant radius
!  r:        radius of slice (m)
!  nlat:     number of samples of slice in latitude direction
!  nlon:     number of samples of slice in longitude direction
!  hslice:   interpolated values of fields on horizontal slice
!
   subroutine extractHorizontalSliceRegularSphericalGrid(this,r,nlat,nlon,hslice,slname)
      class (regular_spherical_grid) :: this
      double precision :: r
      integer :: nlon,nlat
      character (len=*) :: slname
      class (regular_spherical_grid) :: hslice
      integer :: i,j,nr,ic
      double precision :: dlon,dlat,dr,lon,lat
      double precision, dimension(:), allocatable :: field
   !
      allocate(field(this%nfield))
      nr = 1; dr = 0.d0
      dlat = (this%nlat-1)*this%dlat/(nlat-1)
      dlon = (this%nlon-1)*this%dlon/(nlon-1)
      write(6,'(a,2f12.3)') 'Min limits of horizontal slice: ',this%lonmin/mc_deg2rad,this%latmin/mc_deg2rad
      write(6,'(a,2f12.3)') 'Max limits of horizontal slice: ',&
                            (this%lonmin+(nlon-1)*dlon)/mc_deg2rad,(this%latmin+(nlat-1)*dlat)/mc_deg2rad
   !
   !  create slice as 2D RegularSphericalGrid object
   !
      call hslice%create(this%nfield,this%lat_pole,this%lon_pole,nr,nlat,nlon,dr,dlat,dlon,&
                         r,this%latmin,this%lonmin,this%rearth,this%field_name_string,slname)
   !
   !  Loop over points in slice
   !
      ic = 0
      do j = 1,nlat
         lat = this%latmin+(j-1)*dlat
         do i = 1,nlon
            ic = ic+1
            lon = this%lonmin+(i-1)*dlon
            call this%getFieldsAtPoint(r,lat,lon,field)
            hslice%field(ic,:) = real(field)
         end do
      end do
      deallocate(field)
   end subroutine extractHorizontalSliceRegularSphericalGrid
!-----------------------------------------------------------------------------
!  extract a vertical slice along a great circle
!  lat1, lon1:       geographical latitude and longitude of 1st end point (rad)
!  lat2, lon2:       geographical latitude and longitude of 2nd end point (rad)
!  nr:               number of samples in radial direction
!  ndel:             number of samples along great circle
!  vslice:           interpolated values of fields on slice
!  gc:               a GC regular spherical grid object containing geographical coordinates of the great circle
!
   subroutine extractVerticalSliceRegularSphericalGrid(this,lat1,lon1,lat2,lon2,nr,ndel,slname,vslice,gc)
      class (regular_spherical_grid) :: this
      double precision :: lat1,lon1,lat2,lon2
      integer :: nr,ndel
      character (len=*) :: slname
      class (regular_spherical_grid) :: vslice
      class (regular_spherical_grid) :: gc
      double precision, dimension(:), allocatable :: field
      double precision :: pih,xg,yg,zg,x,y,z,r,dis,xi,delstep,dr,delta,theta,lat,lon,latg,long,rd,delmin,dlon
      integer :: k,j,ic,nlon
   !
   !  calculate distance dis and azimuth xi between end points
   !  GC of 2nd end point --> LC relative to 1st end point --> LS relative to first end point
   !
      pih = 0.5d0*mc_pid
      call coordinatesLCfromLSAxesRotation(this%rmin,pih-lat2,lon2,xg,yg,zg)
      call coordinatesLCfromGCAxesRotation(pih-lat1,lon1,xg,yg,zg,x,y,z)
      call coordinatesLSfromLCAxesRotation(x,y,z,rd,dis,xi)
   !
   !  spacing
   !
      delstep = dis/(ndel-1)
      dr = (this%nr-1)*this%dr/(nr-1)
      nlon = 1; dlon = 0.d0; delmin = 0.d0
      write(6,'(a,f15.0)') 'Min radius of vertical slice: ',this%rmin
      write(6,'(a,f15.0)') 'Max radius of vertical slice: ',this%rmin+(nr-1)*dr
   !
   !  create slice as 2D RegularSphericalGrid object
   !
      call vslice%create(this%nfield,lat1,lon1,nr,ndel,nlon,dr,delstep,dlon,&
                         this%rmin,delmin,xi,this%rearth,this%field_name_string,slname)
      call gc%create(2,lat1,lon1,1,ndel,1,0.d0,delstep,0.d0,&
                     this%rmin,delmin,xi,this%rearth,'latitude:longitude',slname)
   !
   !  loop over samples of slice from end point 1 to end point 2
   !  LC of point relative to 1st end point --> GC --> LS relative north pole
   !
      allocate(field(this%nfield))
      ic = 0
      do k = 1,nr
         r = this%rmin+(k-1)*dr
         do j = 1,ndel
            ic = ic+1
            delta = delmin+(j-1)*delstep
            !  get geographical coordinates of point on profile
            call coordinatesLCfromLSAxesRotation(r,delta,xi,x,y,z)
            call coordinatesGCfromLCAxesRotation(pih-lat1,lon1,x,y,z,xg,yg,zg)
            call coordinatesLSfromLCAxesRotation(xg,yg,zg,rd,theta,lon)
            lat = pih-theta
            !  get coordinates of point relative to the pole the grid coordinates refer to
            call this%transformPointToPole(lat,lon,latg,long)
            !  extract field values
            call this%getFieldsAtPoint(r,latg,long,field)
            vslice%field(ic,:) = real(field)
            !  store geographical coordinates of gc-point
            if (k == 1) then
               gc%field(j,1) = lat
               gc%field(j,2) = lon
            end if
         end do
      end do
      deallocate(field)
   end subroutine extractVerticalSliceRegularSphericalGrid
!-------------------------------------------------------------------------
!  write grid to HDF file
!
   subroutine writeHDFRegularSphericalGrid(this,filename,errmsg)
      class (regular_spherical_grid) :: this
      character(len=*) :: filename
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer (kind=8) :: fid
      integer, dimension(:), allocatable :: id
      real, dimension(:), allocatable :: d
   !
      call createFileHDFWrapper(trim(filename),fid,errmsg)
      if (.level.errmsg == 2) return
   !
      call writeStringAttributeHDFWrapper(fid,'rsgname',trim(this%rsgname),errmsg)
      if (.level.errmsg == 2) return
   !
      call writeStringAttributeHDFWrapper(fid,'field_names',trim(this%field_name_string),errmsg)
      if (.level.errmsg == 2) return
   !
      id = [this%nlon,this%nlat,this%nr]
      call aria%assoc1d(id)
      call writeArrayAttributeHDFWrapper(fid,"dimensions",aria,errmsg)
      call aria%deassoc(); deallocate(id)
      if (.level.errmsg == 2) return
   !
      d = [real(this%lon_pole),real(this%lat_pole)]
      call arra%assoc1d(d)
      call writeArrayAttributeHDFWrapper(fid,"crs_pole",arra,errmsg)
      call arra%deassoc(); deallocate(d)
      if (.level.errmsg == 2) return
   !
      d = [real(this%dlon),real(this%dlat),real(this%dr)]
      call arra%assoc1d(d)
      call writeArrayAttributeHDFWrapper(fid,"grid_spacings_rad_m",arra,errmsg)
      call arra%deassoc(); deallocate(d)
      if (.level.errmsg == 2) return
   !
      d = [real(this%lonmin),real(this%latmin),real(this%rmin)]
      call arra%assoc1d(d)
      call writeArrayAttributeHDFWrapper(fid,"minimum_values_rad_m",arra,errmsg)
      call arra%deassoc(); deallocate(d)
      if (.level.errmsg == 2) return
   !
   !  write field array
   !
      call arra%assoc2d(this%field)
      call writeArrayHDFWrapper(fid,'field_values',arra,errmsg,xferprp = hdf_xferprp)
      call arra%deassoc()
      if (.level.errmsg == 2) return
   !
      call closeFileHDFWrapper(fid,errmsg)
      if (.level.errmsg == 2) return
   end subroutine writeHDFRegularSphericalGrid
!-------------------------------------------------------------------------
!  read grid from HDF file
!
   subroutine readHDFRegularSphericalGrid(this,filename,errmsg)
      class (regular_spherical_grid) :: this
      character(len=*) :: filename
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer (kind=8) :: fid
      integer :: slen
      integer, dimension(:), pointer :: id
      real, dimension(:), pointer :: d
      real, dimension(:,:), pointer :: field
      character(len=max_length_string) :: cval
   !
      call openFileRoHDFWrapper(trim(filename),fid,errmsg)
      if (.level.errmsg == 2) return
   !
      call readStringAttributeHDFWrapper(fid,'rsgname',cval,slen,errmsg)
      if (.level.errmsg == 2) return
      this%rsgname = cval(1:slen)
   !
      call readStringAttributeHDFWrapper(fid,'field_names',cval,slen,errmsg)
      if (.level.errmsg == 2) return
      this%field_name_string = cval(1:slen)
   !
      call readArrayAttributeHDFWrapper(fid,"dimensions",aria,errmsg)
      if (.level.errmsg == 2) return
      id => aria%get1d()
      this%nlon = id(1); this%nlat = id(2); this%nr = id(3)
      this%npts = id(1)*id(2)*id(3)
      deallocate(id)
   !
      call readArrayAttributeHDFWrapper(fid,"crs_pole",arra,errmsg)
      if (.level.errmsg == 2) return
      d => arra%get1d()
      this%lon_pole = d(1); this%lat_pole = d(2)
      deallocate(d)
   !
      call readArrayAttributeHDFWrapper(fid,"grid_spacings_rad_m",arra,errmsg)
      if (.level.errmsg == 2) return
      d => arra%get1d()
      this%dlon = d(1); this%dlat = d(2); this%dr = d(3)
      deallocate(d)
   !
      call readArrayAttributeHDFWrapper(fid,"minimum_values_rad_m",arra,errmsg)
      if (.level.errmsg == 2) return
      d => arra%get1d()
      this%lonmin = d(1); this%latmin = d(2); this%rmin = d(3)
      deallocate(d)
   !
   !  read field array
   !
      call readArrayHDFWrapper(fid,'field_values',arra,errmsg,xferprp = hdf_xferprp)
      if (.level.errmsg == 2) return
      field => arra%get2d()
      this%nfield = size(field,2)
      allocate(this%field(this%npts,this%nfield))
      this%field = field
      deallocate(field)
   !
      call closeFileHDFWrapper(fid,errmsg)
      if (.level.errmsg == 2) return
   end subroutine readHDFRegularSphericalGrid
!
end module regularSphericalGrid
