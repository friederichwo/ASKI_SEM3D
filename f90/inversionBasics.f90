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
!> \brief holds basic requirements for inversion process
!!
!! \details Basic requirements for an inversion procedure which are
!!  independent of a specific iteration step, such as the main parameter file,
!!  events and stations, station componenent transformation coefficients, etc.
!!  are supervised by this module. 
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module inversionBasics
!
   use inputParameter
   use string
   use readEventStationFile
   use seismicEventList
   use seismicNetwork
   use errorMessage
   use propertySet
!
   implicit none
!
   interface dealloc; module procedure deallocateInversionBasics; end interface
   interface init; module procedure initiateInversionBasics; end interface
   interface operator (.iterpath.); module procedure getCurrentIterPathInversionBasics; end interface
   interface operator (.inpar.); module procedure getInputParameterInversionBasics; end interface
   interface operator (.evlist.); module procedure getEventListInversionBasics; end interface
   interface operator (.statlist.); module procedure getStationListInversionBasics; end interface
   interface operator (.ifreq.); module procedure getMeasuredDataFrequencyIndicesInversionBasics; end interface
   interface operator (.df.); module procedure getMeasuredDataFrequencyStepInversionBasics; end interface
   interface operator (.ufmdata.); module procedure getMeasuredDataUnitFactorInversionBasics; end interface
   interface operator (.propset.); module procedure getPropertySetInversionBasics; end interface
!
   type inversion_basics
      character(len=max_length_string) :: iter_path       ! absolute path to directory of iteration step
      type (input_parameter) :: inpar                     ! content of main parameter file
      type (seismic_event_list) :: event_list             ! list of all events involved in this inversion
      type (seismic_network) :: station_list              ! list of all stations involved in this inversion
      type (property_set) :: propset                      ! property set object
      integer, dimension(:), allocatable :: ifreq         ! measured data frequency indices
   end type inversion_basics
!
   contains
!
!------------------------------------------------------------------------
!  Initiate basic inversion requirements
!  Main parfile is read. Everything needed for
!  the inversion, independent of each iteration step.
!
   subroutine initiateInversionBasics(this,parfile,lu,errmsg)
      type (inversion_basics) :: this
      character(len=*) :: parfile
      integer :: lu
      type (error_message) :: errmsg
      integer :: j,nf,ios
      character(len=max_length_string) :: main_path,errstr,propsetname,event_file,station_file
      character(len=23) :: myname = 'initiateInversionBasics'
      character (len=80), dimension(40) :: main_inpar_keys
      data main_inpar_keys/'MAIN_PATH_INVERSION', 'PROPERTY_SET_NAME', 'PROPERTY_CORRELATION_FILE',&
         'CURRENT_ITERATION_STEP',  'ITERATION_STEP_PATH', 'PARFILE_ITERATION_STEP', 'INVERSION_GRID',&
         'FORWARD_METHOD', 'REARTH', 'INVGRID_CENTER_LAT', 'INVGRID_CENTER_LON', 'INVGRID_DIMENSIONS',&
         'PATH_GEMINI', 'PATH_MEASURED_DATA', 'PATH_MEASURED_SEIS', 'PATH_INJECTION_SEIS',&
         'PATH_PHASE_END_TIMES', 'FILE_EVENT_LIST', 'FILE_STATION_LIST',&
         'MEASURED_DATA_FREQUENCY_STEP', 'MEASURED_DATA_NUMBER_OF_FREQ', 'MEASURED_DATA_INDEX_OF_FREQ',&
         'UNIT_FACTOR_MEASURED_DATA', 'DEFAULT_VTK_FILE_FORMAT', 'PATH_ASKI_MAIN_FILES', 'PATH_DMSPACE',&
         'PATH_KERNEL_DISPLACEMENTS','PATH_KERNEL_GREEN_TENSORS', 'PATH_SENSITIVITY_KERNELS',&
         'PATH_SYNTHETIC_DATA', 'PATH_OUTPUT_FILES', 'PATH_VTK_FILES', 'ASKI_wx',&
         'ASKI_wy', 'ASKI_wz', 'PSEUDO_MESH_SPACING', 'FILE_ASKI_BACKGROUND_MODEL',&
         'ASKI_TTIME_TABLE_FILE_KD', 'ASKI_TTIME_TABLE_FILE_GT', 'PATH_SPECFEM_INPUT'/
   !--------------------------------------------------------------------------------------------------
   !  read input parameters
   !
      call createKeywordsInputParameter(this%inpar,main_inpar_keys)
      call readSubroutineInputParameter(this%inpar,lu,trim(parfile),errmsg)
      if (.level.errmsg == 2) return
   !
      main_path = sval(this%inpar,'MAIN_PATH_INVERSION')
      event_file = main_path+sval(this%inpar,'FILE_EVENT_LIST')
      station_file = main_path+sval(this%inpar,'FILE_STATION_LIST')
   !
   !  some checks
   !
      select case(this%inpar.sval.'DEFAULT_VTK_FILE_FORMAT')
      case('ASCII','BINARY')
         ! valid, do nothing
      case default
         write(errstr,*) "parameter 'DEFAULT_VTK_FILE_FORMAT' = '"//trim(this%inpar.sval.'DEFAULT_VTK_FILE_FORMAT')//&
              "' of main parfile is not valid; must be either 'ASCII' or 'BINARY'"
         call add(errmsg,2,trim(errstr),myname)
         return
      end select
   !
      nf = ival(this%inpar,'MEASURED_DATA_NUMBER_OF_FREQ',iostat=ios)
      this%ifreq = ivecp(this%inpar,'MEASURED_DATA_INDEX_OF_FREQ',nf,iostat=ios)
      do j = 1,nf
         if (any(this%ifreq(j+1:nf) == this%ifreq(j))) then
            write(errstr,*) "vector 'MEASURED_DATA_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
                 'MEASURED_DATA_INDEX_OF_FREQ')//"' must not contain multiple entries"
            call add(errmsg,2,trim(errstr),myname)
            return
         end if
      enddo
   !
   ! Check whether paths given in the parameter file are empty (i.e. not set by the user).
   !
      if(trim(this%inpar.sval.'PATH_ASKI_MAIN_FILES') == '') then
         call add(errmsg,2,"PATH_ASKI_MAIN_FILES not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_KERNEL_DISPLACEMENTS') == '') then
         call add(errmsg,2,"PATH_KERNEL_DISPLACEMENTS not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_KERNEL_GREEN_TENSORS') == '') then
         call add(errmsg,2,"PATH_KERNEL_GREEN_TENSORS not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_SENSITIVITY_KERNELS') == '') then
         call add(errmsg,2,"PATH_SENSITIVITY_KERNELS not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_SYNTHETIC_DATA') == '') then
         call add(errmsg,2,"PATH_SYNTHETIC_DATA not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_OUTPUT_FILES') == '') then
         call add(errmsg,2,"PATH_OUTPUT_FILES not specified",myname)
         return
      end if
      if(trim(this%inpar.sval.'PATH_VTK_FILES') == '') then
         call add(errmsg,2,"PATH_VTK_FILES not specified",myname)
         return
      end if
   !
   !  create property set object
   !
      propsetname = this%inpar.sval.'PROPERTY_SET_NAME'
      if (equalString(propsetname,'isoVelocitySI')) then
         call createIsoVelocitySIPropertySet(this%propset)
      else
         call add(errmsg,2,'Property set named '//trim(propsetname)//' not implemented',myname)
         return
      end if
      call readCorrmatPropertySet(this%propset,1,this%inpar.sval.'PROPERTY_CORRELATION_FILE',errmsg)
      if (.level.errmsg == 2) return
   !
      write(this%iter_path,"(2a,i3.3,a)") trim(this%inpar.sval.'MAIN_PATH_INVERSION'),&
            trim(this%inpar.sval.'ITERATION_STEP_PATH'), &
            this%inpar.ival.'CURRENT_ITERATION_STEP','/'
   ! ------------------------------------------------------------------------
   !  read event file, create event list
   !
      call createEventListFromEventFile(trim(event_file),lu,'ASKI_events',this%event_list,errmsg)
      if (.level.errmsg == 2) return
   ! ------------------------------------------------------------------------
   !  read station file, create station list
   !
      call createStationListFromStationFile(station_file,lu,'ASKI_stations',this%station_list,errmsg)
      if (.level.errmsg == 2) return
   !
      if(.csys.this%event_list /= .csys.this%station_list) then
         call add(errmsg,1,"the coordinate system '"//.csys.this%event_list//"' used in event list differs from '"&
              //.csys.this%station_list//"' used in station list",myname)
      end if
   end subroutine initiateInversionBasics
!------------------------------------------------------------------------
!  deallocate inversion basics object
!
   subroutine deallocateInversionBasics(this)
      type (inversion_basics) :: this
      call dealloc(this%inpar)
      call dealloc(this%event_list)
      call dealloc(this%station_list)
      call dealloc(this%propset)
      deallocate(this%ifreq)
   end subroutine deallocateInversionBasics
!------------------------------------------------------------------------
!  get absolute path to directory of current iteration step
!
   function getCurrentIterPathInversionBasics(this) result(iter_path)
      type (inversion_basics), intent(in) :: this
      character(len=350) :: iter_path
      iter_path = this%iter_path
   end function getCurrentIterPathInversionBasics
!------------------------------------------------------------------------
!  get input_parameter object contained in inversion_basics object
!
   function getInputParameterInversionBasics(this) result(inpar)
      type (inversion_basics), intent(in), target :: this
      type (input_parameter), pointer :: inpar
      inpar => this%inpar
   end function getInputParameterInversionBasics
!------------------------------------------------------------------------
!  get event_list object contained in inversion_basics object
!
   function getEventListInversionBasics(this) result(event_list)
      type (inversion_basics), intent(in), target :: this
      type (seismic_event_list), pointer :: event_list
      event_list => this%event_list
   end function getEventListInversionBasics
!------------------------------------------------------------------------
!  get station_list object contained in inversion_basics object
!
   function getStationListInversionBasics(this) result(station_list)
      type (inversion_basics), intent(in), target :: this
      type (seismic_network), pointer :: station_list
      station_list => this%station_list
   end function getStationListInversionBasics
!------------------------------------------------------------------------
!  get measured data frequency indices contained in inversion_basics object
!
   function getMeasuredDataFrequencyIndicesInversionBasics(this) result(ifreq)
      type (inversion_basics), intent(in), target :: this
      integer, dimension(:), pointer :: ifreq
      ifreq => this%ifreq
   end function getMeasuredDataFrequencyIndicesInversionBasics
!------------------------------------------------------------------------
!  get measured data frequency step contained in inversion_basics parfile
!
   function getMeasuredDataFrequencyStepInversionBasics(this) result(df)
      type (inversion_basics), intent(in) :: this
      double precision :: df
      df = (this%inpar).dval.'MEASURED_DATA_FREQUENCY_STEP'
   end function getMeasuredDataFrequencyStepInversionBasics
!------------------------------------------------------------------------
!  get measured data unit factor contained in inversion_basics parfile
!
   function getMeasuredDataUnitFactorInversionBasics(this) result(uf)
      type (inversion_basics), intent(in) :: this
      double precision :: uf
      uf = (this%inpar).dval.'UNIT_FACTOR_MEASURED_DATA'
   end function getMeasuredDataUnitFactorInversionBasics
!------------------------------------------------------------------------
!  get property set contained in inversion_basics object
!
   function getPropertySetInversionBasics(this) result(pset)
      type (inversion_basics), intent(in), target :: this
      type (property_set), pointer :: pset
      pset => this%propset
   end function getPropertySetInversionBasics
!
end module inversionBasics
