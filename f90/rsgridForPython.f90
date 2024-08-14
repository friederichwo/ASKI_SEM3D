!----------------------------------------------------------------------------
!   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option)
!   any later version.
!
!   ASKI version 1.2 is distributed in the hope that it
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!------------------------------------------------------------------------------------------
!  Tools to deal with regular spherical grids
!  All angles in radians
!------------------------------------------------------------------------------------------
!  Given an array of points in the 3D grid, return array of flattened indices of grid points
!  located below and south-west of a given point. Index array is zero-based.
!  If any point is off the grid, the index is set to -1.
!
   subroutine field_values_at_points(dims,spacings,minvals,field,intyp,outtyp,r,lat,lon,val,n,npts,nfield,nrp,nlatp,nlonp)
      implicit none
      integer, dimension(3), intent(in) :: dims
      double precision, dimension(3), intent(in) :: spacings, minvals
      double precision, dimension(nfield,npts), intent(in) :: field
      character(len=*) :: intyp,outtyp
      double precision, dimension(nrp), intent(in) :: r
      double precision, dimension(nlatp), intent(in) :: lat
      double precision, dimension(nlonp), intent(in) :: lon
      double precision, dimension(nfield,n), intent(out) :: val
      integer, intent(in) :: n,nfield,npts,nrp,nlatp,nlonp
      integer :: ip,nr,nlat,nlon,i,j,k,m,ir,ilat,ilon,iph,m2,l
      integer, dimension(8) :: ic
      double precision :: rmin,latmin,lonmin,dr,dlat,dlon,rmax
      double precision :: long,latg,rg,lo,la,rn
      double precision, dimension(8) :: w
   !
      nr = dims(3); nlat = dims(2); nlon = dims(1)
      rmin = minvals(3); latmin = minvals(2); lonmin = minvals(1)
      dr = spacings(3); dlat = spacings(2); dlon = spacings(1)
      rmax = rmin+(nr-1)*dr
   !  write(42,*) 'Call to field_values_at_points'
   !  write(42,*) nfield,npts,intyp,outtyp
   !
   !  Loop over points at which field values are desired
   !
      do ip = 1,n
         val(:,ip) = 0.d0
         select case (trim(outtyp))
         case ('HS')
            m = int((ip-1)/nlonp)
            ilat = m+1
            ilon = ip-m*nlonp
            ir = 1
         case ('VS')
            m = int((ip-1)/nlonp)
            ir = m+1
            ilon = ip-m*nlonp
            ilat = ilon
         case ('3D')
            m = int((ip-1)/(nlonp*nlatp))
            ir = m+1
            iph = ip-m*(nlatp*nlonp)
            m2 = int((iph-1)/nlonp)
            ilat = m2+1
            ilon = iph-m2*nlonp            
         case default
            ir = 1; ilon = 1; ilat = 1
            write(6,'(a)') 'Slice type: ',trim(outtyp),' not implemented'
            return
         end select
   !
   ! check if point is off grid
   !
         if (lat(ilat) < latmin .or. lat(ilat) > latmin+(nlat-1)*dlat) then
            cycle
         else if (lon(ilon) < lonmin .or. lon(ilon) > lonmin+(nlon-1)*dlon) then
            cycle
         else if (r(ir) < rmin .or. r(ir) > rmax) then 
            if (abs(r(ir)-rmax)/dr < 0.005) then          ! avoid setting values at surface to zero
               continue
            else
               cycle
            endif
         end if
   !
   ! get unflattened indices of point below-south-west
   !
         k = ceiling((r(ir)-rmin)/dr)
         j = ceiling((lat(ilat)-latmin)/dlat)
         i = ceiling((lon(ilon)-lonmin)/dlon)
         if (k == 0) k = 1
         if (j == 0) j = 1
         if (i == 0) i = 1
         if (k == nr) k = nr-1        ! allows extrapolation towards surface
   !
   !  normalized coordinates of point
   !
         long = lonmin+(i-1)*dlon
         latg = latmin+(j-1)*dlat
         rg = rmin+(k-1)*dr
         lo = (lon(ilon)-long)/dlon
         la = (lat(ilat)-latg)/dlat
         rn = (r(ir)-rg)/dr
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
   !
   !  flattened indices of corner points
   !
         ic(1) = i+nlon*(j-1+(k-1)*nlat)
         ic(2) = ic(1)+1
         ic(3) = ic(1)+nlon
         ic(4) = ic(1)+nlat*nlon
         ic(5) = ic(4)+1
         ic(6) = ic(4)+nlon
         ic(7) = ic(3)+1
         ic(8) = ic(6)+1
   !
   !  special treatment of 2D and 1D objects:
   !  horizontal slice         (nlon>1, nr=1, nlat>1)
   !  great circle             (nlon=1, nr=1, nlat>1)
   !  vertical slice           (nlon=1, nr>1, nlat>1)
   !  radial line              (nlon=1, nr>1, nlat=1)
   !
         if (trim(intyp) == 'HS') then
            val(:,ip) = val(:,ip) + w(1)*field(:,ic(1)) + w(2)*field(:,ic(2)) + w(3)*field(:,ic(3)) + w(7)*field(:,ic(7))
         else if (intyp == 'VS') then
            val(:,ip) = val(:,ip) + w(1)*field(:,ic(1)) + w(3)*field(:,ic(3)) + w(4)*field(:,ic(4)) + w(6)*field(:,ic(6))
         else if (intyp == 'GC') then
            val(:,ip) = val(:,ip) + w(1)*field(:,ic(1)) + w(3)*field(:,ic(3))
         else if (intyp == 'RL') then
            val(:,ip) = val(:,ip) + w(1)*field(:,ic(1)) + w(4)*field(:,ic(4))
         else
            do l = 1,8
               val(:,ip) = val(:,ip) + w(l)*field(:,ic(l))
            end do
         end if
 !        if (mod(ip,100) == 0) write(42,*) ip,ilon,ilat,ir,lat(ilat),lon(ilon),r(ir),val(1,ip),field(1,ic(1))
      end do
   end subroutine field_values_at_points


