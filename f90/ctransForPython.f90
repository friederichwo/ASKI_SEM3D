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
!  Axes rotation Fortran code for Python
!  All angles in radians
!------------------------------------------------------------------------------------------
!  Return cartesian LC coordinates of points in sphere from given LS coordinates.
!  Both LS and LC refer to the same pole.
!
   subroutine from_ls_to_lc(r,lat,lon,x,y,z,n)
      double precision, dimension(n), intent(in) :: r,lat,lon
      double precision, dimension(n), intent(out) :: x,y,z
      integer, intent(in) :: n
      integer :: i
      double precision :: cd,sd,cx,sx
      do i = 1,n
         cd = sin(lat(i)); sd = cos(lat(i))
         cx = cos(lon(i)); sx = sin(lon(i))
         x(i) = r(i)*sd*cx
         y(i) = r(i)*sd*sx
         z(i) = r(i)*cd
      enddo
   end subroutine from_ls_to_lc
!--------------------------------------------------------------------------------------
!  Return spherical LS coordinates of points in sphere from given LC coordinates
!  Both LS and LC refer to the same pole.
!
   subroutine from_lc_to_ls(x,y,z,r,lat,lon,n)
      double precision, dimension(n), intent(in) :: x,y,z
      double precision, dimension(n), intent(out) :: r,lat,lon
      integer, intent(in) :: n
      integer :: i
      double precision :: pih
      pih = 3.141592653589793/2.d0
      do i = 1,n
         r(i) = sqrt(x(i)**2+y(i)**2+z(i)**2)
         lat(i) = pih-acos(z(i)/r(i))
         lon(i) = atan2(y(i),x(i))
      enddo
   end subroutine from_lc_to_ls
!--------------------------------------------------------------------------------------
!  Return cartesian GC coordinates of points in sphere from given LC coordinates.
!  GC refers to the true North Pole while LC refers to a pole with
!  geographical coordinates (latpole, lonpole).
!
   subroutine from_lc_to_gc(latpole,lonpole,x,y,z,xg,yg,zg,n)
      double precision :: latpole,lonpole
      double precision, dimension(n), intent(in) :: x,y,z
      double precision, dimension(n), intent(out) :: xg,yg,zg
      integer, intent(in) :: n
      integer :: i
      double precision :: ct,st,cp,sp,pih
      pih = 3.141592653589793/2.d0
      ct = cos(pih-latpole); st = sin(pih-latpole)
      cp = cos(lonpole); sp = sin(lonpole)
      do i = 1,n
         xg(i) = ct*cp*x(i) - sp*y(i) + st*cp*z(i)
         yg(i) = ct*sp*x(i) + cp*y(i) + st*sp*z(i)
         zg(i) =   -st*x(i)        +    ct*z(i)
      enddo
   end subroutine from_lc_to_gc
!--------------------------------------------------------------------------------------
!  Return cartesian LC coordinates of points in sphere from given GC coordinates.
!  GC refers to the true North Pole while LC refers to a pole at (latpole, lonpole).
!
   subroutine from_gc_to_lc(latpole,lonpole,xg,yg,zg,x,y,z,n)
      double precision :: latpole,lonpole
      double precision, dimension(n), intent(in) :: xg,yg,zg
      double precision, dimension(n), intent(out) :: x,y,z
      integer, intent(in) :: n
      integer :: i
      double precision :: ct,st,cp,sp,pih
      pih = 3.141592653589793/2.d0
      ct = cos(pih-latpole); st = sin(pih-latpole)
      cp = cos(lonpole); sp = sin(lonpole)
      do i = 1,n
         x(i) = ct*cp*xg(i) + ct*sp*yg(i) - st*zg(i)
         y(i) =   -sp*xg(i) +    cp*yg(i)
         z(i) = st*cp*xg(i) + st*sp*yg(i) + ct*zg(i)
      enddo
   end subroutine from_gc_to_lc
!---------------------------------------------------------------------------------------
!  Return cartesian RC coordinates of points in sphere from given LC coordinates.
!
   subroutine from_lc_to_rc(gamma,x,y,z,xr,yr,zr,n)
      double precision, intent(in) :: gamma
      double precision, dimension(n), intent(in) :: x,y,z
      double precision, dimension(n), intent(out) :: xr,yr,zr
      integer, intent(in) :: n
      call from_gc_to_lc(0.d0,gamma,x,y,z,xr,yr,zr,n)
   end subroutine from_lc_to_rc
!---------------------------------------------------------------------------------------
!  Return cartesian LC coordinates of points in sphere from given RC coordinates.
!
   subroutine from_rc_to_lc(gamma,xr,yr,zr,x,y,z,n)
      double precision, intent(in) :: gamma
      double precision, dimension(n), intent(out) :: x,y,z
      double precision, dimension(n), intent(in) :: xr,yr,zr
      integer, intent(in) :: n
      call from_lc_to_gc(0.d0,gamma,xr,yr,zr,x,y,z,n)
   end subroutine from_rc_to_lc
!---------------------------------------------------------------------------------------
!  Return spherical LS coordinates of point in sphere from given geographical coordinates.
!  LS cordinates refer to pole with geographical coordinates (latpole, lonpole).
!  Geographical coordinates refer to true North Pole.
!
   subroutine from_gs_to_ls(latpole,lonpole,glat,glon,llat,llon,n)
      double precision, intent(in) :: latpole,lonpole
      double precision, dimension(n), intent(in) :: glat,glon
      double precision, dimension(n), intent(out) :: llat,llon
      integer, intent(in) :: n
      integer :: i
      double precision :: ct,st,cp,sp,cd,sd,cx,sx,pih
      double precision :: xg,yg,zg,x,y,z,r
      pih = 3.141592653589793/2.d0
      ct = cos(pih-latpole); st = sin(pih-latpole)
      cp = cos(lonpole); sp = sin(lonpole)
      do i = 1,n
         cd = sin(glat(i)); sd = cos(glat(i))
         cx = cos(glon(i)); sx = sin(glon(i))
   ! GS to GC
         xg = sd*cx
         yg = sd*sx
         zg = cd
   ! GC to LC
         x = ct*cp*xg + ct*sp*yg - st*zg
         y =   -sp*xg +    cp*yg
         z = st*cp*xg + st*sp*yg + ct*zg
   ! LC to LS
         r = sqrt(x**2+y**2+z**2)
         llat(i) = pih-acos(z/r)
         llon(i) = atan2(y,x)
      enddo
   end subroutine from_gs_to_ls
!---------------------------------------------------------------------------------------
!  Return geographical coordinates of point in sphere from given LS coordinates.
!  LS cordinates refer to pole with geographical coordinates (latpole, lonpole).
!  Geographical coordinates refer to true North Pole.
!
   subroutine from_ls_to_gs(latpole,lonpole,llat,llon,glat,glon,n)
      double precision, intent(in) :: latpole,lonpole
      double precision, dimension(n), intent(out) :: glat,glon
      double precision, dimension(n), intent(in) :: llat,llon
      double precision :: xg,yg,zg,x,y,z,r
      integer, intent(in) :: n
      integer :: i
      double precision :: ct,st,cp,sp,cd,sd,cx,sx,pih
      pih = 3.141592653589793/2.d0
      ct = cos(pih-latpole); st = sin(pih-latpole)
      cp = cos(lonpole); sp = sin(lonpole)
      do i = 1,n
         cd = sin(llat(i)); sd = cos(llat(i))
         cx = cos(llon(i)); sx = sin(llon(i))
   ! LS to LC
         x = sd*cx
         y = sd*sx
         z = cd
   ! LC to GC
         xg = ct*cp*x - sp*y + st*cp*z
         yg = ct*sp*x + cp*y + st*sp*z
         zg =   -st*x        +    ct*z
   ! GC to LS
         r = sqrt(xg**2+yg**2+zg**2)
         glat(i) = pih-acos(zg/r)
         glon(i) = atan2(yg,xg)
      enddo
   end subroutine from_ls_to_gs

