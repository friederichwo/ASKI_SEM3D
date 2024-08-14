!----------------------------------------------------------------------------
!   Copyright 2016, 2020 Florian Schumacher, Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!  Implements a 2nd version to ASKI_1D background model in model_aski.
!  Used here for implementing an intial model on the inversion grid
!
module askiBackgroundModel
   use errorMessage
   use string
   implicit none
   interface dealloc; module procedure deallocASKIBackgroundModel; end interface
   type aski_background_model
      integer :: nlay
      double precision :: rearth
      integer, dimension(:), allocatable :: nnodes                                     ! number of nodes per layer
      double precision, dimension(:,:), allocatable :: depth,rho,vp,vs,Qmu,Qkappa      ! depth and parameter arrays
      double precision, dimension(:,:), allocatable :: sprho,spvp,spvs,spQmu,spQkappa  ! spline parameters p = s''(zj)
   end type aski_background_model
!
contains
!------------------------------------------------------------------------------
!  deallocate ASKI background model
!
   subroutine deallocASKIBackgroundModel(this)
      type (aski_background_model) :: this
      if (allocated(this%nnodes)) deallocate(this%nnodes)
      if (allocated(this%depth)) deallocate(this%depth)
      if (allocated(this%rho)) deallocate(this%rho)
      if (allocated(this%vp)) deallocate(this%vp)
      if (allocated(this%vs)) deallocate(this%vs)
      if (allocated(this%Qmu)) deallocate(this%Qmu)
      if (allocated(this%Qkappa)) deallocate(this%Qkappa)
      if (allocated(this%sprho)) deallocate(this%sprho)
      if (allocated(this%spvp)) deallocate(this%spvp)
      if (allocated(this%spvs)) deallocate(this%spvs)
      if (allocated(this%spQmu)) deallocate(this%spQmu)
      if (allocated(this%spQkappa)) deallocate(this%spQkappa)
   end subroutine deallocASKIBackgroundModel
!------------------------------------------------------------------------------
!  read background model from a table in a file
!
   subroutine readASKIBackgroundModel(this,filename,lu,rearth,errmsg)
      type (aski_background_model) :: this
      character (len=*) :: filename
      integer :: lu
      double precision :: rearth
      type (error_message) :: errmsg
      real, dimension(:), allocatable :: d,u,wrho,wvp,wvs,wQmu,wQkappa
      integer :: nmax,il,in,i,ier,nlay
      character(len=max_length_string) :: line
      character (len=23) :: myname = 'readASKIBackgroundModel'
   !
      open(unit=lu,file=trim(filename),status='unknown',action='read',iostat=ier)
      if (ier .ne. 0) then
         call add(errmsg,2,"could not open model file",myname)
         return
      endif
   !
      this%rearth = rearth
      read(lu,*) line                     ! ignore first line!
      read(lu,*) line                     ! ignore line with zmax (not used here)
      read(lu,*) nlay
      write(*,*) 'nlayers = ',nlay
      this%nlay = nlay
      allocate(this%nnodes(nlay))
      read(lu,*) this%nnodes

      if(minval(this%nnodes) .le. 1) then
          call add(errmsg,2,"number of nodes per layer must be at least 2",myname)
      endif

      nmax = maxval(this%nnodes)
      allocate(this%depth(nlay,nmax),this%rho(nlay,nmax), &
            this%vp(nlay,nmax),this%vs(nlay,nmax),this%Qmu(nlay,nmax),this%Qkappa(nlay,nmax),&
            this%sprho(nlay,nmax),this%spvp(nlay,nmax),this%spvs(nlay,nmax),&
            this%spQmu(nlay,nmax),this%spQkappa(nlay,nmax))
      this%depth(:,:) = 0.
      this%rho(:,:) = 0.
      this%vp(:,:) = 0.
      this%vs(:,:) = 0.
      this%Qmu(:,:) = 0.
      this%Qkappa(:,:) = 0.
      this%sprho(:,:) = 0.
      this%spvp(:,:) = 0.
      this%spvs(:,:) = 0.
      this%spQmu(:,:) = 0.
      this%spQkappa(:,:) = 0.

      do il = 1,nlay
         do in = 1,this%nnodes(il)
            read(lu,*) this%depth(il,in),this%rho(il,in),this%vp(il,in),this%vs(il,in),&
                        this%Qmu(il,in),this%Qkappa(il,in)
         end do ! in
      enddo ! il

      close(lu)
   !
   ! now compute the splines for each layer (and each parameter)
   ! in form of values for sprho,spvp,spvs = p = s''(zj)
   ! (procedure, as well as notation of p,s,xj as written in "Algorithms" by
   ! Robert Sedgewick, ADDISON-WESLEY 2002, Chapter 38)
   !
      allocate(d(nmax),u(nmax),wrho(nmax),wvp(nmax),wvs(nmax),wQmu(nmax),wQkappa(nmax))
   !
   ! for each layer calculate the second derivative of the respective spline at all nodes
   ! in case of this%nnodes(il) == 2, the spline interpolation (in that case linear interpolation) works, as sprho,spvp,spvs = 0.
   !
      do il = 1,nlay
         if(this%nnodes(il) .ge. 3) then
            d(:) = 0.       ! initiate temporary variables and calculate their values
            u(:) = 0.
            wrho(:) = 0.
            wvp(:) = 0.
            wvs(:) = 0.
            wQmu(:) = 0.
            wQkappa(:) = 0.
            do i = 2,this%nnodes(il) - 1
               d(i) = 2.*(this%depth(il,i+1)-this%depth(il,i-1))
            enddo
            do i = 1,this%nnodes(il) - 1
               u(i) = this%depth(il,i+1)-this%depth(il,i)
            enddo
            do i = 2,this%nnodes(il) - 1
               wrho(i) = 6.*((this%rho(il,i+1)-this%rho(il,i))/u(i) - &
                             (this%rho(il,i)-this%rho(il,i-1))/u(i-1))
               wvp(i) = 6.*((this%vp(il,i+1)-this%vp(il,i))/u(i) - &
                            (this%vp(il,i)-this%vp(il,i-1))/u(i-1))
               wvs(i) = 6.*((this%vs(il,i+1)-this%vs(il,i))/u(i) - &
                            (this%vs(il,i)-this%vs(il,i-1))/u(i-1))
               wQmu(i) = 6.*((this%Qmu(il,i+1)-this%Qmu(il,i))/u(i) - &
                             (this%Qmu(il,i)-this%Qmu(il,i-1))/u(i-1))
               wQkappa(i) = 6.*((this%Qkappa(il,i+1)-this%Qkappa(il,i))/u(i) - &
                             (this%Qkappa(il,i)-this%Qkappa(il,i-1))/u(i-1))
            enddo
         !
         ! now calculate the second derivatives of the spline, assuming them being zero at the extremal nodes (natural boundary conditions)
         !
            this%sprho(il,1) = 0.; this%sprho(il,this%nnodes(il)) = 0.
            this%spvp(il,1) = 0.; this%spvp(il,this%nnodes(il)) = 0.
            this%spvs(il,1) = 0.; this%spvs(il,this%nnodes(il)) = 0.
            this%spQmu(il,1) = 0.; this%spQmu(il,this%nnodes(il)) = 0.
            this%spQkappa(il,1) = 0.; this%spQkappa(il,this%nnodes(il)) = 0.
         !
         ! then calculate the others by solving a tridiagonal system of equations
         !
            if(this%nnodes(il) > 3) then
               do i = 2,this%nnodes(il) - 2
                  wrho(i+1) = wrho(i+1) - wrho(i)*u(i)/d(i)
                  wvp(i+1) = wvp(i+1) - wvp(i)*u(i)/d(i)
                  wvs(i+1) = wvs(i+1) - wvs(i)*u(i)/d(i)
                  wQmu(i+1) = wQmu(i+1) - wQmu(i)*u(i)/d(i)
                  wQkappa(i+1) = wQkappa(i+1) - wQkappa(i)*u(i)/d(i)
                  d(i+1) = d(i+1) - (u(i)**2)/d(i)
               enddo
            endif
         !
            do i = this%nnodes(il)-1,2,-1
               this%sprho(il,i) = (wrho(i) - u(i)*this%sprho(il,i+1))/d(i)
               this%spvp(il,i) = (wvp(i) - u(i)*this%spvp(il,i+1))/d(i)
               this%spvs(il,i) = (wvs(i) - u(i)*this%spvs(il,i+1))/d(i)
               this%spQmu(il,i) = (wQmu(i) - u(i)*this%spQmu(il,i+1))/d(i)
               this%spQkappa(il,i) = (wQkappa(i) - u(i)*this%spQkappa(il,i+1))/d(i)
            enddo
         endif   ! nnodes(il) >= 3
      enddo ! il
      deallocate(d,u,wrho,wvp,wvs,wQmu,wQkappa)
      ! done calculating splines
   end subroutine readASKIBackgroundModel
!------------------------------------------------------------------------------------
!  evaluate background model at some depth
!
   subroutine evalASKIBackgroundModel(this,depth,rho,vp,vs,qmu,qkappa,errmsg)
      type (aski_background_model) :: this
      double precision :: depth,rho,vp,vs,qmu,qkappa,t
      type (error_message) :: errmsg
      logical :: values_defined
      integer :: in,il,nlay
   !
      nlay = this%nlay
      if(depth > this%depth(nlay,this%nnodes(nlay)) .or. depth < this%depth(1,1) ) then
         call add(errmsg,2,'evaluation depth outside depth range','evalASKIBackgroundModel')
         return
      endif
   !
   ! find layer of evaluation depth
   !
      do il = 1,nlay
         if(depth <= this%depth(il,this%nnodes(il))) exit
      enddo
   !
   !  after this loop, il should be the index of the layer
   !  which contains the current point (each layer contains its bottom depth,
   !  the first layer contains the surface of the earth)
   !
      values_defined = .false.
      do in = 2,this%nnodes(il)
         if(depth <= this%depth(il,in)) then
         ! interpolate values at current depth
            t = (depth - this%depth(il,in-1)) / &
                (this%depth(il,in) - this%depth(il,in-1))
            rho = t*this%rho(il,in) + (1.-t)*this%rho(il,in-1) + &
                  (this%depth(il,in)-this%depth(il,in-1))**2 * &
                  ((t**3-t)*this%sprho(il,in) + ((1-t)**3-(1-t))*this%sprho(il,in-1))/6.
            vp = t*this%vp(il,in) + (1.-t)*this%vp(il,in-1) + &
                 (this%depth(il,in)-this%depth(il,in-1))**2 * &
                 ((t**3-t)*this%spvp(il,in) + ((1-t)**3-(1-t))*this%spvp(il,in-1))/6.
            vs = t*this%vs(il,in) + (1.-t)*this%vs(il,in-1) + &
                 (this%depth(il,in)-this%depth(il,in-1))**2 * &
                 ((t**3-t)*this%spvs(il,in) + ((1-t)**3-(1-t))*this%spvs(il,in-1))/6.
            qmu = t*this%Qmu(il,in) + (1.-t)*this%Qmu(il,in-1) + &
                 (this%depth(il,in)-this%depth(il,in-1))**2 * &
                 ((t**3-t)*this%spQmu(il,in) + ((1-t)**3-(1-t))*this%spQmu(il,in-1))/6.
            qkappa = t*this%Qkappa(il,in) + (1.-t)*this%Qkappa(il,in-1) + &
                 (this%depth(il,in)-this%depth(il,in-1))**2 * &
                 ((t**3-t)*this%spQkappa(il,in) + ((1-t)**3-(1-t))*this%spQkappa(il,in-1))/6.
         ! position between the nodes is found, so leave the do loop
            values_defined = .true.
            exit
         endif
      enddo ! in
   !
   ! after this loop, there should have been found an interpolation
   ! if not, raise an error
   !
      if(.not.values_defined) then
         call add(errmsg,2,'no model values defined at evaluation depth','evalASKIBackgroundModel')
         return
      endif
   end subroutine evalASKIBackgroundModel
end module askiBackgroundModel

