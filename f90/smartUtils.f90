!----------------------------------------------------------------------------
!   Copyright 2023 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------
module smartUtils
   use mathConstants
   implicit none
!
contains
!--------------------------------------------------------------------
!  calculate linear cell index from the three indices of a grid point
!  ic =  = k+(j-1)*nk+(i-1)*nk*nj
!  nk: number of grid points in k-direction
!  nj: number of grid points in j-direction
!  i: layer index
!  j: row index
!  k: index along row
!
   function cellIndexFromGridIndicesSmartUtils(i,j,k,nk,nj) result(ic)
      integer :: i,j,k,nk,nj,ic
      ic = k+nk*(j-1+(i-1)*nj)
   end function cellIndexFromGridIndicesSmartUtils
!---------------------------------------------------------------------
!  calculate three indices of a grid point from linear cell index
!  given index of grid point i,j,k the global cell index is ic = k+(j-1)*nk+(i-1)*nk*nj
!  given index of cell ic, then (with r <= 1 !!):
!  ic/(nk*nj) = i-1+r => i = ic/(nk*nj)+1  or i = ic/(nk*nj) if r = 1
!  n = ic-(i-1)*nk*nj = k+(j-1)*nk
!  n/nk = j-1+r => j = n/nk+1  or j = n/nk if r = 1
!  n-(j-1)*nk = k
!  ic: linear cell index
!  nk: number of grid points in k-direction
!  nj: number of grid points in j-direction
!  i: layer index
!  j: row index
!  k: index along row
!
   subroutine gridIndicesFromCellIndexSmartUtils(ic,nk,nj,i,j,k)
      integer :: ic,nk,nj,i,j,k,n,r
      r = mod(ic,nj*nk)
      if (r > 0) then
         i = ic/(nk*nj)+1
      else
         i = ic/(nj*nk)
      end if
      n = ic-(i-1)*nj*nk
      r = mod(n,nk)
      if (r > 0) then
         j = n/nk+1
      else
         j = n/nk
      end if
      k = n-(j-1)*nk
   end subroutine gridIndicesFromCellIndexSmartUtils
!----------------------------------------------------------------------
!  find an even sharing of a linear array of size n among nproc processes
!  equal sharing of elements, remainder is distributed again starting from rank zero
!  regular share: nreg = n/nproc
!  for 0 < myrank < nremain = mod(n,nproc) take one more
!  return number of elements and offset in array for a given process rank
!  n:       size of array
!  nproc:   number of processes
!  rank:    rank of process (starting from zero)
!  nloc:    number of elements for specified rank
!  off::    offset in array of elements for array
!
   subroutine shareLinearArraySmartUtils(n,nproc,rank,nloc,off)
      integer :: n,nproc,rank,nloc,off,nr
      nloc = n/nproc
      nr = n-nloc*nproc
      if (rank < nr) then
        nloc = nloc+1
        off = rank*nloc
      else
        off = rank*nloc+nr
      end if
   end subroutine shareLinearArraySmartUtils
!----------------------------------------------------------------------
!  calculate sin x/x
!
   elemental function sincSmartUtils(x) result(res)
      double precision, intent(in) :: x
      double precision :: res
      if (dabs(x) < 1.d-6) then
         res = 1.d0
      else
         res = dsin(x)/x
      end if
   end function sincSmartUtils
! ---------------------------------------------------------
!  interpolate data to sampling ponts of Specfem synthetics
!  using windowed sin x/x -interpolation to preserve the spectrum
!  nsamp: number of samples of input data
!  tanf: start time in seconds after midnight of input data
!  dt: sampling interval of input data
!  u: values of input data
!  tanew: desired start time of output data in seconds after midnight (in)
!  dtnew: desired sampling rate of output data (in)
!  nzmax: index of cutoff zero crossing of sinc (values beyond are not considered for convolution)
!  nsnew: number of samples of output (out)
!  unew: values of output samples (out)
!
   subroutine sincInterpolateData(nsamp,tanf,dt,u,tanew,dtnew,nzmax,nsnew,unew)
      integer :: nsamp,nsnew,nzmax
      double precision :: tanf,dt,tanew,dtnew
      real, dimension(:) :: u
      real, dimension(:), allocatable :: unew
      integer :: j,n
      double precision :: t,ts,arg,w,h,taplen
   !
   !  new number of samples, go beyond last input sample or exactly there
   !
      nsnew = ceiling((tanf+(nsamp-1)*dt-tanew)/dtnew)+1
      allocate(unew(nsnew))
      do j = 1,nsnew
         t = tanew+(j-1)*dtnew
         unew(j) = 0.0
         do n = 1,nsamp
            ts = tanf+(n-1)*dt
            arg = mc_pid*(t-ts)/dt
            taplen = nzmax*dt
            w = dcos(0.5*arg*dt/taplen)**2         ! taper sinc with a Hann window
            if (abs(arg) <= mc_pid*nzmax) then
               h = sincSmartUtils(arg)
               unew(j) = unew(j)+u(n)*h*w
            end if
         end do
      end do
   end subroutine sincInterpolateData
!----------------------------------------------------------------------
!  calculate values for a for cos**2 front taper of given length
!  ntl: length of taper in samples
!  tapval: array of taper values between 0 and 1 (out)
!
   subroutine frontTaperSmartUtils(ntl,tapval)
      integer :: ntl
      double precision, dimension(:), allocatable :: tapval
      integer :: j
   !
      allocate(tapval(ntl))
      do j = 1,ntl
         tapval(j) = dcos(0.5*mc_pid*(dble(ntl-j)/dble(ntl)))**2
      end do
   end subroutine frontTaperSmartUtils
!----------------------------------------------------------------------
!  calculate values for a for cos**2 tail taper of given length
!  ntl: length of taper in samples
!  tapval: array of taper values between 0 and 1 (out)
!
   subroutine tailTaperSmartUtils(ntl,tapval)
      integer :: ntl
      double precision, dimension(:), allocatable :: tapval
      integer :: j
   !
      allocate(tapval(ntl))
      do j = 1,ntl
         tapval(j) = dcos(0.5*mc_pid*(dble(j)/dble(ntl)))**2
      end do
   end subroutine tailTaperSmartUtils
!---------------------------------------------------------------------------------
! x,y,z:   coordinates of GLL point (in)
! cc:      coordinates of interpolation anchors (3,nc) (in)
! w:       interpolation weigths (nc) (out)
!---------------------------------------------------------------------------------
   subroutine shepardInterpolationSmartUtils(x,y,z,xc,yc,zc,w)
      implicit none
      double precision :: x,y,z
      double precision, dimension(:) :: xc,yc,zc
      double precision, dimension(:) :: w
      ! local variables
      double precision :: r,h,sum_s
      double precision, dimension(:), allocatable :: t,s,d
      integer :: nc,i,j,imin
      double precision :: dmin,dmax
   !
   ! calculate the distances of (x,y,z) to interpolation anchors
   !
      nc = size(xc)
      w = 0.d0
      allocate(d(nc))
      d = sqrt((xc-x)**2+(yc-y)**2+(zc-z)**2)
      dmin = minval(d)
      imin = minloc(d,1)
      dmax = maxval(d)

   ! if the minimum distance to a cell is close to zero, just take the value there.
      if(dmin/dmax < 1.d-6) then
         w(imin) = 1.d0
         deallocate(d)
         return
      end if

   ! choose radius r, beyond which the influence of a cell center will be absolutely zero
   ! choose here dmax
      r = 1.5*dmax

   ! compute values s_i
      allocate(s(nc))
      do i = 1,nc
         if(d(i) < r/3.d0) then
            s(i) = 1.d0/d(i)                        ! near zero was checked above
         else
            h = d(i)/r - 1.d0
            s(i) = 27.*h*h/(4.d0*r)
         end if
      enddo

   ! define interpolation weights from s and
   ! direction factors t (interpolation function f_3 in the paper)
      sum_s = sum(s)
      allocate(t(nc))
      do i = 1,nc
         t(i) = 0.
         do j = 1,nc
            if (j == i) cycle
            h = (xc(i)-x)*(xc(j)-x) + (yc(i)-y)*(yc(j)-y) + (zc(i)-z)*(zc(j)-z)
            t(i) = t(i) + s(j)*(1.-h/(d(i)*d(j)))
         end do ! j
         t(i) = t(i)/sum_s
      end do ! i
      w = s*s*(1.d0+t)
      w = w/sum(w)
      deallocate(s,t,d)
   end subroutine shepardInterpolationSmartUtils
!
end module smartUtils