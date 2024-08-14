!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief hold basic requirements for a specific iteration step
!!
!! \details Basic requirements for a specific iteration step, such as the iteration 
!!  specific parfile, wavefield points, inversion grid, kernel reference model, 
!!  integration weights, etc. are supervised by this module
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module iterationStepBasics
!
   use inversionBasics
   use inputParameter
   use fileUnitHandler
   use errorMessage
   implicit none
!
   interface dealloc; module procedure deallocateIterationStepBasics; end interface
   interface operator (.iterpath.); module procedure getIterationStepPathIterationStepBasics; end interface
   interface operator (.inpar.); module procedure getInputParameterIterationStepBasics; end interface
   interface operator (.nf.); module procedure getNfIterationStepBasics; end interface
   interface operator (.ifreq.); module procedure getFrequencyIndicesIterationStepBasics; end interface
!
   type iteration_step_basics
      character(len=max_length_string) :: iter_path     ! absolute iteration step path (starting with inverison main path)
      type (input_parameter) :: inpar                   ! parfile content
      integer, dimension(:), allocatable :: ifreq       ! iteration step frequency indices
   end type iteration_step_basics
!
   contains
!
!------------------------------------------------------------------------
!  Initiate basic iteration step requirements
!  invbasics:     inversion basics object
!  lu:            file id for reading iter parfile
!  do_hyperslab:  read inversion grid wavefield point related arrays in hyperslab mode
!  errmsg:        error message
!
   subroutine initiateIterationStepBasics(this,invbasics,lu,errmsg)
      type (iteration_step_basics) :: this
      type (inversion_basics) :: invbasics
      integer :: lu
      type (error_message) :: errmsg
      type (input_parameter), pointer :: inpar_inv
      integer :: j,ios,nf
      character(len=max_length_string) :: iter_parfile
      character (len=80), dimension(19) :: iter_inpar_keys
      character(len=27) :: myname = 'initiateIterationStepBasics'
      data iter_inpar_keys/'ITERATION_STEP_NUMBER_OF_FREQ', 'ITERATION_STEP_INDEX_OF_FREQ', 'DO_PHASE_INVERSION', &
           'MODEL_PROPERTIES_INVERTED_FOR', 'CALIBRATION_ERROR_FACTOR',  'NPROC',&
           'MAX_NUM_CG_ITERATIONS', 'NITER_WINDOW_STA', 'NITER_WINDOW_LTA', &
           'VSCAL_SMOOTHING_MANTLE', 'VSCAL_DAMPING_MANTLE', 'VSCAL_SMOOTHING_CRUST', 'VSCAL_DAMPING_CRUST', &
           'DSMBOOST', 'BOUNDCHOKE', 'CRUSTAL_DEPTH', 'VTK_GEOMETRY_TYPE', 'SCALE_VTK_COORDS', &
           'VTK_COORDS_SCALING_FACTOR'/
   !
      this%iter_path = .iterpath.invbasics
      inpar_inv => getInputParameterInversionBasics(invbasics)
      iter_parfile = trim(this%iter_path)//trim(inpar_inv.sval.'PARFILE_ITERATION_STEP')
   !
   ! read input parameters
   !
      call createKeywordsInputParameter(this%inpar,iter_inpar_keys)
      call readSubroutineInputParameter(this%inpar,lu,trim(iter_parfile),errmsg)
      if(.level.errmsg == 2) return
   !
   ! check consistency of entries in iteration step specific parfile
   !
      nf = ival(this%inpar,'ITERATION_STEP_NUMBER_OF_FREQ',iostat=ios)
      if(nf > (inpar_inv.ival.'MEASURED_DATA_NUMBER_OF_FREQ')) then
         call add(errmsg,2,"NF for iteration > NF global",myname)
         return
      end if
   !
      allocate(this%ifreq(nf))
      this%ifreq = ivecp(this%inpar,'ITERATION_STEP_INDEX_OF_FREQ',nf,iostat=ios)
      if(ios /= 0) then
         call add(errmsg,2,"inconsistent number of iteration step frequencies",myname)
         return
      end if
   !
      do j = 1,nf
         if (any(this%ifreq(j+1:nf) == this%ifreq(j))) then
            call add(errmsg,2,"duplicated entries in iteration step frequencies",myname)
            return
         end if
         if(.not. any((.ifreq.invbasics) == this%ifreq(j))) then
            call add(errmsg,2,"Iteration step frequencies not contained in global ones",myname)
            return
         end if
      enddo
   !
   end subroutine initiateIterationStepBasics
!------------------------------------------------------------------------
!  deallocate iteration step basics object
!  this iteration step basics
!
   subroutine deallocateIterationStepBasics(this)
      type (iteration_step_basics) :: this
      call dealloc(this%inpar)
      deallocate(this%ifreq)
   end subroutine deallocateIterationStepBasics
!------------------------------------------------------------------------
!  get iteration step path
!  iterpath iteration step path
!  iteration step path contained in this
!
   function getIterationStepPathIterationStepBasics(this) result(iterpath)
      type (iteration_step_basics), intent(in) :: this
      character(len=350) :: iterpath
      iterpath = this%iter_path
   end function getIterationStepPathIterationStepBasics
!------------------------------------------------------------------------
!  get input_parameter object contained in iteration_step_basics object
!  input_parameter object
!  input_parameter object contained in this
!
   function getInputParameterIterationStepBasics(this) result(inpar)
      type (iteration_step_basics), intent(in), target :: this
      type (input_parameter), pointer :: inpar
      inpar => this%inpar
    end function getInputParameterIterationStepBasics
!------------------------------------------------------------------------
! get frequency indices contained in iteration_step_basics object
! ifreq frequency indices
! frequency indices contained in this
!
   function getFrequencyIndicesIterationStepBasics(this) result(ifreq)
      type (iteration_step_basics), intent(in), target :: this
      integer, dimension(:), pointer :: ifreq
      ifreq => this%ifreq
   end function getFrequencyIndicesIterationStepBasics
!------------------------------------------------------------------------
!  get number of frequency indices contained in iteration_step_basics object
!  this iteration_step_basics object
!  nf number of frequency indices (size(this%ifreq))
!  number of frequency indices contained in this
!
   function getNfIterationStepBasics(this) result(nf)
      type (iteration_step_basics), intent(in) :: this
      integer :: nf
      nf = size(this%ifreq)
   end function getNfIterationStepBasics
!
end module iterationStepBasics
