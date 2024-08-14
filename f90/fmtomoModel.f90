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
!  Implements a 3D velocity model created by FMTOMO
!
module fmtomoModel
   use errorMessage
   use string
   implicit none
   interface dealloc; module procedure deallocFmtomoModel; end interface
   type fmtomo_model
      integer :: ngrid                                         ! number of subgrids
      integer :: velocity_index                                ! (1) = vp, (2) = vs
      integer :: nr,nlat,nlon                                  ! number of grid points in radius, latitude and longitude
      double precision :: dr,dlat,dlon                         ! spacing along each direction (meters, radians)
      double precision :: rmin,latmin,lonmin                   ! minimum values of each spatial coordinate (meters, radians)
      double precision, dimension(:,:,:), allocatable :: vel   ! velocity value (lon,lat,r) (meters/second)
   end type fmtomo_model
!
   contains
!-------------------------------------------------------------------------------
!  read fmtomo model from file
!
   subroutine readFmtomoModel(this,lu,filename,errmsg)
      type (fmtomo_model) :: this
      integer :: lu
      character(len=*) :: filename
      type (error_message) :: errmsg
      integer :: ios,i,j,k
      character(len=15) :: myname = 'readFmtomoModel'
   !
      open(unit=lu,file=filename,status='unknown',form='formatted',action='read',iostat=ios)
      if (ios /= 0) goto 1
   !
      read(lu,*,iostat = ios) this%ngrid,this%velocity_index
      if (ios /= 0) goto 1
      if (this%ngrid /= 1) then
         call add(errmsg,2,'only one grid currently implemented',myname); return
      end if
   !
      read(lu,*,iostat=ios) this%nr,this%nlat,this%nlon
      if (ios /= 0) goto 1
   !
      read(lu,*,iostat=ios) this%dr,this%dlat,this%dlon
      if (ios /= 0) goto 1
      this%dr = this%dr*1.d3
   !
      read(lu,*,iostat=ios) this%rmin,this%latmin,this%lonmin
      if (ios /= 0) goto 1
      this%rmin = this%rmin*1.d3
   !
   !  The actual line in the vgrids file contains more than one value.
   !  These additional values are not relevant for the calculation of the velocity model and thus are not read in.
   !  The value read in is the wave velocity in m/s
   !
      print *,this%nr,this%nlat,this%nlon
      allocate(this%vel(this%nlon,this%nlat,this%nr))
      do i = 1,this%nr
         do j = 1,this%nlat
            do k = 1,this%nlon
               read(lu,*,iostat = ios) this%vel(k,j,i)
               if (ios /= 0) goto 1
            end do
         end do
      end do
      this%vel = this%vel*1.d3
   !
      close(lu)
   !
   !  error handling
   !
1     if (ios /= 0) then
         call add(errmsg,2,'problems reading fmtomo model file',myname)
      end if
   end subroutine readFmtomoModel
!--------------------------------------------------------------------------------
!  deallocate Fmotmo model
!
   subroutine deallocFmtomoModel(this)
      type (fmtomo_model) :: this
      if (allocated(this%vel)) deallocate(this%vel)
   end subroutine deallocFmtomoModel
!--------------------------------------------------------------------------------
!  evaluate fmtomo model at specified point doing 8-point trilinear interpolation
!  if point is outside grid, return a negative value
!
   function evalFmtomoModel(this,r,lat,lon,errmsg) result(res)
      type (fmtomo_model) :: this
      double precision :: r,lat,lon
      type (error_message) :: errmsg
      double precision :: res
      integer :: i,j,k
      double precision :: lonk,latj,ri,lo,la,rn
      double precision :: w1,w2,w3,w4,w5,w6,w7,w8
   !
   !  locate point in grid, e.g. rmin+(i-1)*dr < r <= rmin+i*dr and determine i
   !  (i-1)*dr < r-rmin <= i*dr => i-1 < (r-rmin)/dr <= i
   !
      i = ceiling((r-this%rmin)/this%dr)
      j = ceiling((lat-this%latmin)/this%dlat)
      k = ceiling((lon-this%lonmin)/this%dlon)
      if (r == this%rmin) i = 1
      if (lat == this%latmin) j = 1
      if (lon == this%lonmin) k = 1
   !
   !  return negative value if point is outside grid
   !
      if (i >= this%nr .or. j >= this%nlat .or. k >= this%nlon) then
         res = -1.d0
         return
      end if
      if (i < 1 .or. j < 1 .or. k < 1) then
         res = -1.d0
         return
      end if
   !
   !  compute interpolation weights using normalized coordinates
   !
      lonk = this%lonmin+(k-1)*this%dlon
      latj = this%latmin+(j-1)*this%dlat
      ri = this%rmin+(i-1)*this%dr
   !
      lo = (lon-lonk)/this%dlon
      la = (lat-latj)/this%dlat
      rn = (r-ri)/this%dr
      w1 = (1.0-lo)*(1.0-la)*(1-rn)
      w2 = lo*(1.0-la)*(1.0-rn)
      w3 = (1.0-lo)*la*(1.0-rn)
      w4 = (1.0-lo)*(1.0-la)*rn
      w5 = lo*(1.0-la)*rn
      w6 = (1.0-lo)*la*rn
      w7 = lo*la*(1.0-rn)
      w8 = lo*la*rn
      res = w1*this%vel(k,j,i) + w2*this%vel(k+1,j,i) + w3*this%vel(k,j+1,i) + w4*this%vel(k,j,i+1)+ &
            w5*this%vel(k+1,j,i+1) + w6*this%vel(k,j+1,i+1) + w7*this%vel(k+1,j+1,i) + w8*this%vel(k+1,j+1,i+1)
   end function evalFmtomoModel
!--------------------------------------------------------------
!  get type of model (either vp or vs)
!
   function getVeltypeFmtomoModel(this) result(res)
      type (fmtomo_model) :: this
      character(len=2) :: res
      res = 'vp'
      if (this%velocity_index == 2) res = 'vs'
   end function getVeltypeFmtomoModel
end module fmtomoModel
