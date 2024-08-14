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
!   Abstract base class for kernel displacements
!     - define real kd type as extension to kernel_displacement
!     - real kd type provides definitions of deferred procedures
!     - routines creating real kd type allocate kernel_displacement to real kd type
!     - for example:
!         + declare <class (kernel_displacement), allocatable :: kd> in the creating routine
!         + then for example: <allocate(real_kernel_displacement :: kd> and
!         + call kd%readReal(......)
!     - pointers and dummy arguments of kd used elsewhere should be declared
!       as <class (kernel_displacement), pointer :: kd>
!
module kernelWavefield
   use errorMessage
   use constants
   implicit none
   type, abstract :: kernel_wavefield
      integer :: nwp                             ! number of wavefield points
      integer :: ncomp                           ! number of forces for Green tensor wavefields
   contains
      procedure (readDisplacementKernelWavefield), deferred :: readDisplacement
      procedure (readGreenTensorKernelWavefield), deferred :: readGreenTensor
      procedure (readDisplacementPhaseEndTimeKernelWavefield), deferred :: readDisplacementPhaseEndTime
      procedure (readGreenTensorPhaseEndTimeKernelWavefield), deferred :: readGreenTensorPhaseEndTime
      procedure (getDisplacementKernelWavefield), deferred :: getDisplacement
      procedure (getGreenTensorKernelWavefield), deferred :: getGreenTensor
      procedure (getStrainKernelWavefield), deferred :: getStrain
      procedure (getGreenStrainKernelWavefield), deferred :: getGreenStrain
      procedure (deallocKernelWavefield), deferred :: dealloc
      procedure (deallocDisplacementKernelWavefield), deferred :: deallocDisplacement
      procedure (deallocGreenKernelWavefield), deferred :: deallocGreen
      procedure :: getNwp => getNwpKernelWavefield
      procedure :: getNcomp => getNcompKernelWavefield
   end type kernel_wavefield
   !------------------------------------------------------------------------
   !  Read kernel displacement from file
   !  nwp:            number of local wavefield points
   !  jf:             frequency index
   !  dfin:           frequency step
   !  basename:       base name of file
   !  evid:           event id
   !  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
   !  errmsg:         error message
   !
   abstract interface
      subroutine readDisplacementKernelWavefield(this,nwp,jf,dfin,basename,evid,do_hyperslab,errmsg)
         import kernel_wavefield
         import error_message
         class (kernel_wavefield) :: this
         integer :: nwp,jf
         double precision :: dfin
         character(len=*) :: basename,evid
         logical :: do_hyperslab
         type (error_message) :: errmsg
      end subroutine readDisplacementKernelWavefield
   end interface
   !------------------------------------------------------------------------
   !  Read kernel Green tensor  from file
   !  nwp:            number of local wavefield points
   !  jf:             frequency index
   !  dfin:           frequency step
   !  basename:       base name of file
   !  netstaname:     netcode.sta name of station
   !  comp:           array of single force components
   !  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
   !  errmsg:         error message
   !
   abstract interface
      subroutine readGreenTensorKernelWavefield(this,nwp,jf,dfin,basename,netstaname,comp,do_hyperslab,errmsg)
         import kernel_wavefield
         import error_message
         import char_len_comp
         class (kernel_wavefield) :: this
         integer :: nwp,jf
         double precision :: dfin
         character(len=*) :: basename,netstaname
         character(len=char_len_comp), dimension(:) :: comp
         logical :: do_hyperslab
         type (error_message) :: errmsg
      end subroutine readGreenTensorKernelWavefield
   end interface
   !------------------------------------------------------------------------
   !  Read kernel displacement phase end time from file
   !  nwp:            number of local wavefield points
   !  basename:       base name of file
   !  evid:           event id
   !  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
   !  errmsg:         error message
   !
   abstract interface
      subroutine readDisplacementPhaseEndTimeKernelWavefield(this,nwp,basename,evid,do_hyperslab,errmsg)
         import kernel_wavefield
         import error_message
         class (kernel_wavefield) :: this
         integer :: nwp
         character(len=*) :: basename,evid
         logical :: do_hyperslab
         type (error_message) :: errmsg
      end subroutine readDisplacementPhaseEndTimeKernelWavefield
   end interface
   !------------------------------------------------------------------------
   !  Read kernel Green tensor phase end time from file
   !  nwp:            number of local wavefield points
   !  basename:       base name of file
   !  netstaname:     netcode.sta name of station
   !  comp:           array of single force components
   !  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
   !  errmsg:         error message
   !
   abstract interface
      subroutine readGreenTensorPhaseEndTimeKernelWavefield(this,nwp,basename,netstaname,do_hyperslab,errmsg)
         import kernel_wavefield
         import error_message
         import char_len_comp
         class (kernel_wavefield) :: this
         integer :: nwp
         character(len=*) :: basename,netstaname
         logical :: do_hyperslab
         type (error_message) :: errmsg
      end subroutine readGreenTensorPhaseEndTimeKernelWavefield
   end interface
   !-------------------------------------------------------------------------
   !  Get displacement from source of kernel wavefield
   !
   abstract interface
      function getDisplacementKernelWavefield(this) result(u)
         import kernel_wavefield
         class (kernel_wavefield), intent(in) :: this
         double complex, dimension(:,:), pointer :: u
      end function getDisplacementKernelWavefield
   end interface
   !-------------------------------------------------------------------------
   !  Get Green tensor from receiver of kernel wavefield
   !
   abstract interface
      function getGreenTensorKernelWavefield(this) result(g)
         import kernel_wavefield
         class (kernel_wavefield), intent(in) :: this
         double complex, dimension(:,:,:), pointer :: g
      end function getGreenTensorKernelWavefield
   end interface
   !-------------------------------------------------------------------------
   !  Get strains from source of kernel wavefield
   !
   abstract interface
      function getStrainKernelWavefield(this) result(ustr)
         import kernel_wavefield
         class (kernel_wavefield), intent(in) :: this
         double complex, dimension(:,:), pointer :: ustr
      end function getStrainKernelWavefield
   end interface
   !-------------------------------------------------------------------------
   !  Get Green strains from receiver of kernel wavefield
   !
   abstract interface
      function getGreenStrainKernelWavefield(this) result(gstr)
         import kernel_wavefield
         class (kernel_wavefield), intent(in) :: this
         double complex, dimension(:,:,:), pointer :: gstr
      end function getGreenStrainKernelWavefield
   end interface
   !----------------------------------------------------------------------------------------------------
   !  Deallocate object
   !
   abstract interface
      subroutine deallocKernelWavefield(this)
         import kernel_wavefield
         class (kernel_wavefield) :: this
      end subroutine deallocKernelWavefield
   end interface
   !----------------------------------------------------------------------------------------------------
   !  Deallocate displacement
   !
   abstract interface
      subroutine deallocDisplacementKernelWavefield(this)
         import kernel_wavefield
         class (kernel_wavefield) :: this
      end subroutine deallocDisplacementKernelWavefield
   end interface
   !----------------------------------------------------------------------------------------------------
   !  Deallocate Green
   !
   abstract interface
      subroutine deallocGreenKernelWavefield(this)
         import kernel_wavefield
         class (kernel_wavefield) :: this
      end subroutine deallocGreenKernelWavefield
   end interface
!
   contains
!------------------------------------------------------------------------
!  get number of proc-specific wavefield points
!
   integer function getNwpKernelWavefield(this) result(res)
      class (kernel_wavefield), intent(in) :: this
      res = this%nwp
   end function getNwpKernelWavefield
!------------------------------------------------------------------------
!  get number of forces for Green tensor
!
   integer function getNcompKernelWavefield(this) result(res)
      class (kernel_wavefield), intent(in) :: this
      res = this%ncomp
   end function getNcompKernelWavefield
!
end module kernelWavefield