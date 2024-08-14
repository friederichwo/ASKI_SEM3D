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
module rayKernelWavefield
   use kernelWavefield
   use errorMessage
   use globalMpiInfo
   implicit none
!
   type, extends(kernel_wavefield) :: ray_kernel_wavefield
      integer :: nps                                                   ! number of phases from source
      integer :: npr                                                   ! number of phases from receiver
      double complex, dimension(:,:,:), allocatable :: uamp            ! amplitudes (ncell, 3, nps)
      double complex, dimension(:,:), allocatable :: utime             ! travel times (ncell, nps)
      double precision, dimension(:,:,:), allocatable :: uslowness     ! slownesses (ncell,3,nps)
      double complex, dimension(:,:,:), allocatable :: gamp            ! amplitudes (ncell, 3, npr)
      double complex, dimension(:,:), allocatable :: gtime             ! travel times (ncell, npr)
      double precision, dimension(:,:,:), allocatable :: gslowness     ! slownesses (ncell,3,npr)
      double complex, dimension(:,:,:), pointer :: u => null()                   ! displacements (ncell, 3, nps)
      double complex, dimension(:,:,:), pointer :: ustr => null()                ! strains (ncell, 3, nps)
      double complex, dimension(:,:,:,:), pointer :: g => null()                 ! Green tensor (ncell, 3, forces, npr)
      double complex, dimension(:,:,:,:), pointer :: gstr => null()              ! Green strains (ncell, 6, forces, npr)
      double precision, dimension(:), pointer :: petu                 ! phase end time for kernel displacements
      double precision, dimension(:), pointer :: petg                 ! phase end time for kernel Green tensors
   contains
      procedure :: readDisplacement => readFrequencyDisplacementRayKernelWavefield
      procedure :: readGreenTensor => readFrequencyGreenTensorRayKernelWavefield
      procedure :: getDisplacement => getDisplacementRayKernelWavefield
      procedure :: getGreenTensor => getGreenTensorRayKernelWavefield
      procedure :: getStrain => getStrainRayKernelWavefield
      procedure :: getGreenStrain => getGreenStrainRayKernelWavefield
      procedure :: dealloc => deallocRayKernelWavefield
      procedure :: deallocDisplacement => deallocDisplacementRayKernelWavefield
      procedure :: deallocGreen => deallocGreenRayKernelWavefield
      procedure :: getPhaseStrain => getPhaseStrainRayKernelWavefield
      procedure :: getPhaseDisplacement => getPhaseDisplacementRayKernelWavefield
      procedure :: getPhaseGreenStrain => getPhaseGreenStrainRayKernelWavefield
      procedure :: getPhaseGreenTensor => getPhaseGreenTensorRayKernelWavefield
      procedure :: getNps => getNpsRayKernelWavefield
      procedure :: getNpr => getNprRayKernelWavefield
      procedure :: getWavenumberSum => getWavenumberSumRayKernelWavefield
      procedure :: readDisplacementPhaseEndTime => readDisplacementPhaseEndTimeRayKernelWavefield
      procedure :: readGreenTensorPhaseEndTime => readGreenTensorPhaseEndTimeRayKernelWavefield
      procedure :: getDisplacementPhaseEndTime => getDisplacementPhaseEndTimeRayKernelWavefield
      procedure :: getGreenTensorPhaseEndTime => getGreenTensorPhaseEndTimeRayKernelWavefield
   end type ray_kernel_wavefield
!
contains
!------------------------------------------------------------------------------
!  Read amplitude, travel time and slowness for all phases arriving from source
!  nwp:  here uses as number of phases (either output parameter or for checking)
!  Only read file once to get amplitude, travel time and slownesses which are frequency independent
!  Calculate u and ustr from them at further calls to read
!
   subroutine readFrequencyDisplacementRayKernelWavefield(this,nwp,jf,dfin,basename,evid,do_hyperslab,errmsg)
      class (ray_kernel_wavefield) :: this
      integer :: nwp,jf
      double precision :: dfin
      character(len=*) :: basename,evid
      logical :: do_hyperslab
      type (error_message) :: errmsg
   !   
      return
   end subroutine readFrequencyDisplacementRayKernelWavefield
!------------------------------------------------------------------------------
!  Read amplitude, travel time and slowness for all phases arriving from receiver
!  nwp:  here uses as number of phases (either output parameter or for checking)
!  Only read file once to get amplitude, travel time and slownesses which are frequency independent
!  Calculate g and gstr from them at further calls to read
!
   subroutine readFrequencyGreenTensorRayKernelWavefield(this,nwp,jf,dfin,basename,netstaname,comp,do_hyperslab,errmsg)
      class (ray_kernel_wavefield) :: this
      integer :: nwp,jf
      double precision :: dfin
      character(len=*) :: basename,netstaname
      character(len=char_len_comp), dimension(:) :: comp
      logical :: do_hyperslab
      type (error_message) :: errmsg
   !   
      return
   end subroutine readFrequencyGreenTensorRayKernelWavefield
!------------------------------------------------------------------------------
!  Read amplitude, travel time and slowness for all phases arriving from source
!  nwp:  here uses as number of phases (either output parameter or for checking)
!  Only read file once to get amplitude, travel time and slownesses which are frequency independent
!  Calculate u and ustr from them at further calls to read
!
   subroutine readDisplacementPhaseEndTimeRayKernelWavefield(this,nwp,basename,evid,do_hyperslab,errmsg)
      class (ray_kernel_wavefield) :: this
      integer :: nwp
      character(len=*) :: basename,evid
      logical :: do_hyperslab
      type (error_message) :: errmsg
   !
      return
   end subroutine readDisplacementPhaseEndTimeRayKernelWavefield
!------------------------------------------------------------------------------
!  Read amplitude, travel time and slowness for all phases arriving from receiver
!  nwp:  here uses as number of phases (either output parameter or for checking)
!  Only read file once to get amplitude, travel time and slownesses which are frequency independent
!  Calculate g and gstr from them at further calls to read
!
   subroutine readGreenTensorPhaseEndTimeRayKernelWavefield(this,nwp,basename,netstaname,do_hyperslab,errmsg)
      class (ray_kernel_wavefield) :: this
      integer :: nwp
      character(len=*) :: basename,netstaname
      logical :: do_hyperslab
      type (error_message) :: errmsg
   !
      return
   end subroutine readGreenTensorPhaseEndTimeRayKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate object
!
   subroutine deallocRayKernelWavefield(this)
      class (ray_kernel_wavefield) :: this
      if (associated(this%u)) deallocate(this%u)
      if (associated(this%ustr)) deallocate(this%ustr)
      if (associated(this%g)) deallocate(this%g)
      if (associated(this%gstr)) deallocate(this%gstr)
   end subroutine deallocRayKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate displacements of object
!
   subroutine deallocDisplacementRayKernelWavefield(this)
      class (ray_kernel_wavefield) :: this
      if (associated(this%u)) deallocate(this%u)
      if (associated(this%ustr)) deallocate(this%ustr)
   end subroutine deallocDisplacementRayKernelWavefield
!-------------------------------------------------------------------------
!  Associate displacements of this with those of that
!
   subroutine associateDisplacementRayKernelWavefield(this,u,ustr)
      class (ray_kernel_wavefield) :: this
      double complex, dimension(:,:,:), target :: u,ustr
   end subroutine associateDisplacementRayKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate Green tensor of object
!
   subroutine deallocGreenRayKernelWavefield(this)
      class (ray_kernel_wavefield) :: this
      if (associated(this%g)) deallocate(this%g)
      if (associated(this%gstr)) deallocate(this%gstr)
   end subroutine deallocGreenRayKernelWavefield
!-------------------------------------------------------------------------
!  Get total strains (sum over all phases)
!
   function getStrainRayKernelWavefield(this) result(ustr)
      class (ray_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:), pointer :: ustr
      integer :: ncell,i
      ncell = size(this%ustr,1)
      allocate(ustr(ncell,6))
      ustr = 0.d0
      do i = 1,this%nps
         ustr = ustr + this%ustr(:,:,i)
      enddo
   end function getStrainRayKernelWavefield
!-------------------------------------------------------------------------
!  Get total displacement
!
   function getDisplacementRayKernelWavefield(this) result(u)
      class (ray_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:), pointer :: u
      integer :: ncell,i
      ncell = size(this%u,1)
      allocate(u(ncell,3))
      u = 0.d0
      do i = 1,this%nps
         u = u + this%u(:,:,i)
      end do
   end function getDisplacementRayKernelWavefield
!-------------------------------------------------------------------------
!  Get total Green strains
!
   function getGreenStrainRayKernelWavefield(this) result(gstr)
      class (ray_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:,:), pointer :: gstr
      integer :: ncell,ncomp,i
      ncell = size(this%gstr,1)
      ncomp = size(this%gstr,3)
      allocate(gstr(ncell,6,ncomp))
      gstr = 0.d0
      do i = 1,this%npr
         gstr = gstr + this%gstr(:,:,:,i)
      enddo
   end function getGreenStrainRayKernelWavefield
!-------------------------------------------------------------------------
!  Get total Green tensor
!
   function getGreenTensorRayKernelWavefield(this) result(g)
      class (ray_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:,:), pointer :: g
      integer :: ncell,ncomp,i
      ncell = size(this%g,1)
      ncomp = size(this%g,3)
      allocate(g(ncell,6,ncomp))
      g = 0.d0
      do i = 1,this%npr
         g = g + this%g(:,:,:,i)
      enddo
   end function getGreenTensorRayKernelWavefield
!-------------------------------------------------------------------------
!  Get displacement phase end time
!
   function getDisplacementPhaseEndTimeRayKernelWavefield(this) result(petu)
      class (ray_kernel_wavefield), intent(in) :: this
      double precision, dimension(:), pointer :: petu
      petu => this%petu
   end function getDisplacementPhaseEndTimeRayKernelWavefield
!-------------------------------------------------------------------------
!  Get Green tensor phase end time
!
   function getGreenTensorPhaseEndTimeRayKernelWavefield(this) result(petg)
      class (ray_kernel_wavefield), intent(in) :: this
      double precision, dimension(:), pointer :: petg
      petg => this%petg
   end function getGreenTensorPhaseEndTimeRayKernelWavefield
!-------------------------------------------------------------------------
!  Get phase strains
!
   function getPhaseStrainRayKernelWavefield(this,ip) result(ustr)
      class (ray_kernel_wavefield), intent(in), target :: this
      double complex, dimension(:,:), pointer :: ustr
      integer :: ip
      ustr => this%ustr(:,:,ip)
   end function getPhaseStrainRayKernelWavefield
!-------------------------------------------------------------------------
!  Get phase displacement
!
   function getPhaseDisplacementRayKernelWavefield(this,ip) result(u)
      class (ray_kernel_wavefield), intent(in), target :: this
      double complex, dimension(:,:), pointer :: u
      integer :: ip
      u => this%u(:,:,ip)
   end function getPhaseDisplacementRayKernelWavefield
!-------------------------------------------------------------------------
!  Get phase Green strains
!
   function getPhaseGreenStrainRayKernelWavefield(this,ip) result(gstr)
      class (ray_kernel_wavefield), intent(in), target :: this
      double complex, dimension(:,:,:), pointer :: gstr
      integer :: ip
      gstr => this%gstr(:,:,:,ip)
   end function getPhaseGreenStrainRayKernelWavefield
!-------------------------------------------------------------------------
!  Get phase Green tensor
!
   function getPhaseGreenTensorRayKernelWavefield(this,ip) result(g)
      class (ray_kernel_wavefield), intent(in), target :: this
      double complex, dimension(:,:,:), pointer :: g
      integer :: ip
      g => this%g(:,:,:,ip)
   end function getPhaseGreenTensorRayKernelWavefield
!-------------------------------------------------------------------------
!  Get number of source phases
!
   function getNpsRayKernelWavefield(this) result(res)
      class (ray_kernel_wavefield), intent(in) :: this
      integer :: res
      res = this%nps
   end function getNpsRayKernelWavefield
!-------------------------------------------------------------------------
!  Get number of receiver phases
!
   function getNprRayKernelWavefield(this) result(res)
      class (ray_kernel_wavefield), intent(in) :: this
      integer :: res
      res = this%npr
   end function getNprRayKernelWavefield
!------------------------------------------------------------------------
!  Get difference of wavenumbers for given s and r
!
   subroutine getWavenumberSumRayKernelWavefield(this,js,jr,omega,q)
      class (ray_kernel_wavefield) :: this
      integer :: js,jr
      double precision :: omega
      double precision, dimension(:,:), pointer :: q
      integer :: ncell
   !
      ncell = size(this%u,1)
      allocate(q(ncell,3))
      q = omega*(this%uslowness(:,:,js)+this%gslowness(:,:,jr))
   end subroutine getWavenumberSumRayKernelWavefield
!
end module rayKernelWavefield
