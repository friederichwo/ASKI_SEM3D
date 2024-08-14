!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   Define and hold information about a set of data and a set of model values
!
!   This module defines a type containing all information about the data samples
!   and model values, i.e. data space and model space, which should be used for a
!   certain operation like setting up the kernel matrix, doing kernel focussing, etc.
!   A data sample is defined by a station, a station component, an event, a frequency index
!   and distinction between imaginary or real part (of the complex valued data).
!   A model value is defined by a property set (e.g. 'isoVelocity',
!   allowing for properties like 'rho', 'vp', 'vs'), the respective parameter (e.g. 'vs')
!   and an inversion grid cell.
!   Only indices refering to properties defined in a different object (iteration step info,
!   inversion grid) are hold in this module, like station/event index, inversion grid cell
!   index etc.
!   For parallelized application, these data samples and model values may be defined and
!   used locally and do not need to represent the whole data and model space.<BR>
!
!  \author Florian Schumacher, Wolfgang Friederich
!  \date July 2015, July 2022
!------------------------------------------------------------------------
module dataModelSpaceInfo
!
   use constants
   use seismicEventList
   use seismicNetwork
   use propertySet
   use realloc
   use errorMessage
!
   implicit none
!
   interface dealloc; module procedure deallocateDataModelSpaceInfo; end interface
   interface operator (.npath.); module procedure getNpathDataModelSpaceInfo; end interface
   interface operator (.ndata.); module procedure getNdataDataModelSpaceInfo; end interface
   interface operator (.nmval.); module procedure getNmvalDataModelSpaceInfo; end interface
   interface operator (.ncell.); module procedure getNcellDataModelSpaceInfo; end interface
   interface operator (.nprop.); module procedure getNpropDataModelSpaceInfo; end interface
   interface operator (.evid.); module procedure getEvidDataSampleDataModelSpaceInfo; end interface
   interface operator (.netstaname.); module procedure getNetStanameDataSampleDataModelSpaceInfo; end interface
   interface operator (.comp.); module procedure getCompDataSampleDataModelSpaceInfo; end interface
   interface operator (.ifreq.); module procedure getIfreqDataSampleDataModelSpaceInfo; end interface
   interface operator (.imre.); module procedure getImreDataSampleDataModelSpaceInfo; end interface
   interface operator (.prop.); module procedure getPropertyModelValueDataModelSpaceInfo; end interface
   interface operator (.cell.); module procedure getCellModelValueDataModelSpaceInfo; end interface
!
   integer, parameter :: char_len_path = max(char_len_evid,char_len_sta+char_len_netcode+1)
   integer, parameter :: maxcomp_dmsp = 3
!
!  meta info defining rows (data space) and columns (model space) of inversion matrix
!
   type data_model_space_info
      private
      integer :: ncomp_in_dmsp                                                 ! number of components occurring in dmspace
      character(len=char_len_comp), dimension(:), allocatable :: comp_in_dmsp  ! unique list of components occuring in dmspace
   ! data space / paths
      integer :: npath                                                      ! total number of paths defined by event and station
      character(len=char_len_path), dimension(:,:), allocatable :: paths    ! array of paths (evid, netstaname) in data space
      integer, dimension(:), allocatable :: offset_path                     ! offset of data samples belonging to path (npath)
      integer, dimension(:), allocatable :: ncomp_path                      ! number of components belonging to path (npath)
      character(len=char_len_comp), dimension(:,:), allocatable :: comp_path  ! components belonging to path
      integer, dimension(:), allocatable :: nfreq_path                      ! number of frequencies belonging to path (masked)
      integer, dimension(:,:), allocatable :: ifreq_path                    ! frequency indices belonging to path (masked)
      logical, dimension(:,:), allocatable :: maskfreq_path                 ! frequency mask belonging to path
   ! unmasked values
      integer, dimension(:), allocatable :: nfreq_path_ori                  ! number of frequencies belonging to path (unmasked)
      integer, dimension(:,:), allocatable :: ifreq_path_ori                ! frequency indices belonging to path (unmasked)
      double precision, dimension(:,:), allocatable :: wfreq_path_ori       ! frequency weights belonging to path (unmasked)
   ! data space / data samples
      integer :: ndata                                                      ! total number of data samples in data space
      integer, dimension(:), allocatable :: pathidx                         ! path index of data sample
      character(len=char_len_evid), dimension(:), allocatable :: evid       ! event ID of data sample
      character(len=char_len_sta+char_len_netcode+1), dimension(:), allocatable :: netstaname  ! netcode.station_name of data sample
      character(len=char_len_comp), dimension(:), allocatable :: comp       ! component of data sample
      integer, dimension(:), allocatable :: ifreq                           ! frequency index of data sample
      double precision, dimension(:), allocatable :: noise                  ! noise of each data sample = 1/weight
   ! model space
      integer :: nprop                                                      ! number of model properties (same for all model values)
      character(len=char_len_par), dimension(:), allocatable :: prop        ! material properties occuring in model space (e.g. vp,vs), sorted
      integer :: nmval                                                      ! total number of model values in model space (props x cells)
      integer :: ncell                                                      ! total number of inversion grid cells
      character(len=char_len_par), dimension(:), allocatable :: prop_mval   ! property associated with model value (props x cells)
      integer, dimension(:), allocatable :: cell_mval                       ! index of inversion grid cell belonging to model value (props x cells)
   end type data_model_space_info
!
contains
!-------------------------------------------------------------------------------------------
!  check if a key in the data model space info file is valid
!
   logical function invalidKeyDataModelSpaceInfo(key,expected,block,iline,errmsg) result(res)
      character (len=*) :: key,expected,block
      integer :: iline
      type (error_message) :: errmsg
      character (len=400) :: errstr
      if (trim(key) /= trim(expected)) then
         write(errstr,*) "Keyword <"//trim(key)//"> on line ",iline, &
                         " in block <"//trim(block)//"> not supported. <"//trim(expected)//"> expected."
         call add(errmsg,2,trim(errstr),"invalidKeyDataModelSpaceInfo")
         res = .true.
      else
         res = .false.
      endif
   end function invalidKeyDataModelSpaceInfo
!-------------------------------------------------------------------------------------------
!  print that value for a key is invalid
!
   subroutine writeInvalidValueDataModelSpaceInfo(key,val,block,iline,errmsg)
      character (len=*) :: key,val,block
      integer :: iline
      type (error_message) :: errmsg
      character (len=max_length_string) :: errstr
      write(errstr,*) "<"//trim(key)//"> specification <"//trim(val)//"> on line ",iline,&
                      " in block <"//trim(block)//"> not supported."
      call add(errmsg,2,trim(errstr),"writeInvalidValueDataModelSpaceInfo")
   end subroutine writeInvalidValueDataModelSpaceInfo
!---------------------------------------------------------------------------------------
!  Add data samples of one path to existing ones of this data model space info object
!  This routine ASSUMES all incoming values to be valid, there is no check done here!
!  It is called from routines within this module only.
!  Furthermore, it is assumed that all incoming values are combined to ncomp*nfreq*nimre
!  data samples where nimre = 2 and the order is assumed to real followed by imaginary part.
!  evid event of data sample added
!  netstaname net.station name of data samples added
!  comp array of components of data samples
!  ifreq array of frequency indices of data samples (packed to those with mask = .true.)
!  wfreq array of weights for new data samples (packed to those with mask = .true.)
!
   subroutine addDataSamplesDataModelSpaceInfo(this,jpath,evid,netstaname,comp,ifreq,wfreq)
      type (data_model_space_info) :: this
      integer :: jpath
      character(len=*) :: evid,netstaname
      character(len=*), dimension(:) :: comp
      integer, dimension(:) :: ifreq
      double precision, dimension(:) :: wfreq
      ! local
      integer :: ncomp,nfreq_used,n,nimre
      integer :: jcomp,jfreq,jimre
      integer :: ndata_new
   !
      ncomp = size(comp)
      nfreq_used = size(ifreq)
      nimre = 2
      ndata_new = ncomp * nfreq_used * nimre
      if(ndata_new .le. 0) return
   !
   ! reallocate
   !
      n = this%ndata+ndata_new
      this%pathidx = reallocate(this%pathidx,n)
      this%evid = reallocate(this%evid,n)
      this%netstaname = reallocate(this%netstaname,n)
      this%comp = reallocate(this%comp,n)
      this%ifreq = reallocate(this%ifreq,n)
      this%noise = reallocate(this%noise,n)
   !
   ! add samples
   !
      n = this%ndata
      do jcomp = 1,ncomp
         do jfreq = 1,nfreq_used
            do jimre = 1,nimre
               n = n + 1
               this%pathidx(n) = jpath
               this%evid(n) = evid
               this%netstaname(n) = netstaname
               this%comp(n) = comp(jcomp)
               this%ifreq(n) = ifreq(jfreq)
               this%noise(n) = 1.d0/wfreq(jfreq)
            end do
         end do
      end do
      this%ndata = this%ndata + ndata_new
   end subroutine addDataSamplesDataModelSpaceInfo
! ------------------------------------------------------------------------
!  From data model space input file, a data_model_space object is filled
!  Here, simply the two subroutines createDataSamplesFromFileDataModelSpaceInfo
!  and createModelValuesFromFileDataModelSpaceInfo are called
!  this:          data_model_space_info object to be created
!  eventlist:     list of all seismic events, to check validity of event IDs
!  stationlist:   list of all seismic stations, to check validity of station names
!  ifreq_valid:   array of valid frequency indices
!  propset:       property set object
!  ntot_invgrid:  total number of inversion grid cells to check valid range of cell indices in file
!  filename:      file name to use
!  lu:            file unit to use
!  errmsg:        error message
!
   subroutine createFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,&
                                               propset,ntot_invgrid,filename,lu,errmsg)
      type (data_model_space_info) :: this
      type (seismic_event_list) :: eventlist
      type (seismic_network) :: stationlist
      integer, dimension(:) :: ifreq_valid
      type (property_set) :: propset
      integer :: ntot_invgrid
      character(len=*) :: filename
      integer :: lu
      type (error_message) :: errmsg
   !
      call createDataSamplesFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,filename,lu,errmsg)
      if(.level.errmsg == 2) return
   !
      call createModelValuesFromFileDataModelSpaceInfo(this,propset,ntot_invgrid,filename,lu,errmsg)
      if(.level.errmsg == 2) return
   end subroutine createFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!  Create data samples from the DATA SAMPLES block of the data model space input file
!  This routine searches for the line "DATA SAMPLES" in the given file and reads in the
!  block following this line.
!  Format of 'DATA SAMPLES' block:
!  line of form: 'WEIGHTING value', where value is one of 'NONE', 'BY_PATH', 'BY_FREQUENCY', 'BY_PATH_AND_FREQUENCY'
!     - in case of 'NONE', all data weights are internally set to 1.0 (i.e. no actual weighting is performed)
!     - in case of 'BY_PATH' and 'BY_PATH_AND_FREQUENCY'
!       the pair 'evid netstaname' in a path block is expected to be followed by one number > 0.0 and <= 1.0;
!       the frequency dependent weighting values (in cases 'BY_FREQUENCY', 'BY_PATH_AND_FREQUENCY') are defined
!       in a separate line following the lines of form 'nfreq ifreq_1 ... ifreq_n' (for either case of 'FREQUENCIES ALL'
!       or 'FREQUENCIES SPECIFIC'. These separate lines have themselves the form 'nfreq w_1 ... w_n' defining nfreq weights
!       in range > 0.0 and <= 1.0.
!     - in case 'BY_PATH_AND_FREQUENCY', both weights (for path and frequency) are MULTIPLIED for each data sample
!  line of form: 'COMPONENTS value', where value is either 'ALL' (for all paths, the same components are used) or 'SPECIFIC'
!     - if 'COMPONENTS ALL', the next line is of form 'ncomp comp_1 ... comp_n' defining the components for all paths.
!  line of form: 'FREQUENCIES value', where value is either 'ALL' (for all paths, the same frequency indices are used) or 'SPECIFIC'
!     - if 'FREQUENCIES ALL', the next line is of form 'nfreq ifreq_1 ... ifreq_n' defining the frequency indices for all paths.
!  The following line must contain the number of paths npaths, followed by
!  npaths blocks of lines, each defining the path and the data samples for that path. These blocks constist of at least one
!  line containing the eventid/stationname pair 'evid netstaname'.
!  For each keyword 'COMPONENTS' and 'FREQUENCIES', if 'SPECIFIC',
!  one line is added to such a block of lines, in the same form as the line following value 'ALL' (see above, e.g. 'nfreq ifreq_1 ... ifreq_n'),
!  defining the specific components, frequencies  for that path.
!  Note, it is assumed that real and imaginary part of each datum are always provided and in that order.
!  this data model space info
!  eventlist list of all seismic events, to check validity of event IDs
!  stationlist list of all seismic stations, to check validity of station names
!  ifreq_valid array of valid frequency indices
!  filename file name to use
!  lu file unit, with which the input file is already open
!  errmsg error message
!
   subroutine createDataSamplesFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,filename,lu,errmsg)
      type (data_model_space_info) :: this
      type (seismic_event_list) :: eventlist
      type (seismic_network) :: stationlist
      integer, dimension(:) :: ifreq_valid
      character(len=*) :: filename
      integer :: lu
      type (error_message) :: errmsg
   !
      character(len=43) :: myname = 'createDataSamplesFromFileDataModelSpaceInfo'
      character(len=500) :: line
      integer :: ios,iline
      character(len=40) :: key,val
      integer :: npath,jpath,ncomp,jcomp,nfreq,jfreq,maxfreq,nfreq_used
      character(len=char_len_evid) :: evid
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
      character(len=char_len_comp), dimension(:), allocatable :: comp
      integer, dimension(:), allocatable :: ifreq,ifreq_ma
      logical, dimension(:), allocatable :: mask_freq
      double precision, dimension(:), allocatable :: weight_freq, weight_freq_ma
      double precision :: weight_path
      logical :: specific_comp,specific_freq,line_data_samples_found,&
                 weight_by_path,weight_by_freq,mask_by_freq
   !
      iline = 0
      maxfreq = size(ifreq_valid)
   !
      open(unit=lu,file=trim(filename),status='old',form='formatted',action='read',iostat=ios)
      if(ios /= 0) then
         call add(errmsg,2,"could not open file '"//trim(filename)//"'",myname)
         return
      endif
   !
   ! parse through whole file until line "DATA SAMPLES" and start reading in block from there
   !
      line_data_samples_found = .false.
      do while(ios == 0)
         read(lu,"(a)",iostat=ios) line; iline = iline+1
         if(ios /= 0) exit
         if(line == "DATA SAMPLES") then
            line_data_samples_found = .true.
            exit
         endif
      enddo
      if(.not. line_data_samples_found) then
         call add(errmsg,2,"could not find line 'DATA SAMPLES",myname)
         return
      endif
   !
   !  WEIGHTING
   !
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) key,val
      if (invalidKeyDataModelSpaceInfo(key,'WEIGHTING','DATA SAMPLES',iline,errmsg)) return
      select case(trim(val))
      case('NONE')
         weight_by_path = .false.
         weight_by_freq = .false.
      case('BY_FREQUENCY')
         weight_by_path = .false.
         weight_by_freq = .true.
      case('BY_PATH')
         weight_by_path = .true.
         weight_by_freq = .false.
      case('BY_PATH_AND_FREQUENCY')
         weight_by_path = .true.
         weight_by_freq = .true.
      case default
         call writeInvalidValueDataModelSpaceInfo("WEIGHTING",val,'DATA SAMPLES',iline,errmsg)
         return
      end select
   !
   !  MASKING
   !
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) key,val
      if (invalidKeyDataModelSpaceInfo(key,'MASKING','DATA SAMPLES',iline,errmsg)) return
      select case(trim(val))
      case('NONE')
         mask_by_freq = .false.
      case('BY_FREQUENCY')
         mask_by_freq = .true.
      case default
         call writeInvalidValueDataModelSpaceInfo("MASKING",val,'DATA SAMPLES',iline,errmsg)
         return
      end select
   !
   ! COMPONENTS
   !
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) key,val
      if (invalidKeyDataModelSpaceInfo(key,'COMPONENTS','DATA SAMPLES',iline,errmsg)) return
      select case(trim(val))
      case('ALL')
         specific_comp = .false.
         read(lu,"(a)",iostat=ios) line; iline = iline+1
         read(line,*,iostat=ios) ncomp
         if (ncomp > maxcomp_dmsp) then
            call add(errmsg,2,"more than maximum number of components ordered",myname)
            return
         end if
         allocate(comp(ncomp))
         read(line,*,iostat=ios) ncomp,comp
         do jcomp = 1,ncomp
            if (.not. any(valid_components == trim(comp(jcomp)))) then
               call add(errmsg,2,"data model space info file conains invalid components",myname)
               return
            endif
         enddo
         allocate(this%comp_in_dmsp(ncomp))
         this%comp_in_dmsp = comp
         this%ncomp_in_dmsp = ncomp
      case('SPECIFIC') ! val
         specific_comp = .true.
         allocate(this%comp_in_dmsp(maxcomp_dmsp))
         this%ncomp_in_dmsp = 0                            ! initialize component count here
      case default
         call writeInvalidValueDataModelSpaceInfo("COMPONENTS",val,'DATA SAMPLES',iline,errmsg)
         return
      end select
   !
   ! FREQUENCIES
   !
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) key,val
      if (invalidKeyDataModelSpaceInfo(key,'FREQUENCIES','DATA SAMPLES',iline,errmsg)) return
      select case(trim(val))
      case('ALL')
         specific_freq = .false.
         read(lu,"(a)",iostat=ios) line; iline = iline+1
         read(line,*,iostat=ios) nfreq
         allocate(ifreq(nfreq))
         read(line,*,iostat=ios) nfreq,ifreq
      !                                                       check if there are invalid frequencies
         do jfreq = 1,nfreq
            if (.not. any(ifreq_valid == ifreq(jfreq))) then
               call add(errmsg,2,"Invalid frequencies occur in data model space info file",myname)
               return
             endif
          enddo
      !                                                        weighting by frequency
         if(weight_by_freq) then
            read(lu,"(a)",iostat=ios) line; iline = iline+1
            allocate(weight_freq(nfreq))
            read(line,*,iostat=ios) jfreq,weight_freq
         else
            allocate(weight_freq(nfreq))
            weight_freq = 1.0
         end if
      !                                                            masking by frequency (.true. means use datum)
         if(mask_by_freq) then
            read(lu,"(a)",iostat=ios) line; iline = iline+1
            allocate(mask_freq(nfreq))
            read(line,*,iostat=ios) jfreq,mask_freq
         else
            allocate(mask_freq(nfreq))
            mask_freq = .true.
         end if
      !                                                               masked weights and frequency indices
         weight_freq_ma = weight_path*pack(weight_freq,mask_freq)
         ifreq_ma = pack(ifreq,mask_freq)
      case('SPECIFIC') ! val
         specific_freq = .true.
      case default
         call writeInvalidValueDataModelSpaceInfo("FREQUENCIES",val,'DATA SAMPLES',iline,errmsg)
         return
      end select ! val
   !
   ! Iteration over SPECIFIC PATHS
   !
      this%ndata = 0                                             ! initialize data samples count
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) npath
      this%npath = npath
      allocate(this%paths(2,npath))
      allocate(this%offset_path(npath))
      allocate(this%ncomp_path(npath))
      allocate(this%comp_path(maxcomp_dmsp,npath))
      allocate(this%nfreq_path(npath))
      allocate(this%ifreq_path(maxfreq,npath))
      allocate(this%maskfreq_path(maxfreq,npath))
   !
      allocate(this%nfreq_path_ori(npath))
      allocate(this%ifreq_path_ori(maxfreq,npath))
      allocate(this%wfreq_path_ori(maxfreq,npath))
      do jpath = 1,npath
         read(lu,"(a)",iostat=ios) line; iline = iline+1
         if(weight_by_path) then
            read(line,*,iostat=ios) evid,netstaname,weight_path
         else
            read(line,*,iostat=ios) evid,netstaname
         end if
      !                                                                       check if the event ID is valid
         if (.not. searchEventidSeismicEventList(eventlist,evid)) then
            call add(errmsg,2,"Event in data model space info file not present in event list",myname)
            return
         end if
      !                                                                       check if the station name is valid
         if (.not. searchNetcodeStationNameSeismicNetwork(stationlist,netstaname)) then
            call add(errmsg,2,"Station in data model space info file not in station list",myname)
         end if
         this%paths(1,jpath) = evid
         this%paths(2,jpath) = netstaname
      !
      ! SPECIFIC COMPONENTS
      !
         if(specific_comp) then
            read(lu,"(a)",iostat=ios) line; iline = iline+1
            read(line,*,iostat=ios) ncomp
            if (ncomp > maxcomp_dmsp) then
               call add(errmsg,2,"more than maximum number of components ordered",myname)
               return
            end if
            allocate(comp(ncomp))
            read(line,*,iostat=ios) ncomp,comp
         !                                                                   check if there are invalid components
            do jcomp = 1,ncomp
               if (.not. any(valid_components == trim(comp(jcomp)))) then
                  call add(errmsg,2,"data model space info file conains invalid components",myname)
                  return
               endif
               if ((this%ncomp_in_dmsp < maxcomp_dmsp) .and. (.not. any(this%comp_in_dmsp == trim(comp(jcomp))))) then
                  this%ncomp_in_dmsp = this%ncomp_in_dmsp+1
                  this%comp_in_dmsp(this%ncomp_in_dmsp) = comp(jcomp)
               end if
            enddo
         endif
      !
      ! SPECIFIC FREQUENCIES
      !
         if(specific_freq) then
            read(lu,"(a)",iostat=ios) line; iline = iline+1
            read(line,*,iostat=ios) nfreq
            allocate(ifreq(nfreq))
            read(line,*,iostat=ios) nfreq,ifreq
         !                                                             check if there are invalid frequencies
            do jfreq = 1,nfreq
               if (.not. any(ifreq_valid == ifreq(jfreq))) then
                  call add(errmsg,2,"Invalid frequencies occur in data model space info file",myname)
                  return
                endif
            enddo
         !                                                              frequency weights
            if(weight_by_freq) then
               read(lu,"(a)",iostat=ios) line; iline = iline+1
               allocate(weight_freq(nfreq))
               read(line,*,iostat=ios) jfreq,weight_freq
            else ! weight_by_freq
               allocate(weight_freq(nfreq))
               weight_freq = 1.0
            end if
         !                                                              masking by frequency (.true. means use datum)
            if(mask_by_freq) then
               read(lu,"(a)",iostat=ios) line; iline = iline+1
               allocate(mask_freq(nfreq))
               read(line,*,iostat=ios) jfreq,mask_freq
            else
               allocate(mask_freq(nfreq))
               mask_freq = .true.
            end if
         endif
      !
      !  masked weights and frequency indices
      !
         weight_freq_ma = weight_path*pack(weight_freq,mask_freq)
         ifreq_ma = pack(ifreq,mask_freq)
      !
      ! now, for this path all information on components, frequencies, imre's, weights and masks are gathered,
      ! so finally add data samples and define weight and mask array and set data index range for path
      !
         nfreq_used = count(mask_freq)
         this%offset_path(jpath) = this%ndata
         this%ncomp_path(jpath) = ncomp
         this%comp_path(1:ncomp,jpath) = comp
         this%nfreq_path(jpath) = nfreq_used                                    ! only count usable frequencies
         this%ifreq_path(1:nfreq_used,jpath) = ifreq_ma
         this%maskfreq_path(1:nfreq,jpath) = mask_freq
      !                                                                         ! original values, no masking
         this%nfreq_path_ori(jpath) = nfreq
         this%ifreq_path_ori(1:nfreq,jpath) = ifreq
         this%wfreq_path_ori(1:nfreq,jpath) = weight_path*weight_freq
      !
         call addDataSamplesDataModelSpaceInfo(this,jpath,evid,netstaname,comp,ifreq_ma,weight_freq_ma)
         if(specific_comp) deallocate(comp)
         if(specific_freq) deallocate(ifreq,weight_freq,mask_freq,ifreq_ma,weight_freq_ma)
      enddo ! jpath
      !
      ! if arrays comp,ifreq,imre were allocated above (before jpath-loop), deallocate here
      !
      if (.not. specific_comp) deallocate(comp)
      if (.not. specific_freq) deallocate(ifreq,weight_freq,mask_freq,ifreq_ma,weight_freq_ma)
      !
      close(lu)
   end subroutine createDataSamplesFromFileDataModelSpaceInfo
!-------------------------------------------------------------------------------------
!  create model values from the MODEL VALUES block of the data model space input file
!  This routine searches for the line "MODEL VALUES" in the given file and reads in the
!  block following this line.
!  Format of 'MODEL VALUES' block:
!    line of form 'nprop prop_1 ... prop_n',
!    defining the properties of all inversion grid cells (assumed to be properies of given property set).
!  this:      data model space info
!  propset:   all parameters are assumed to be of this property set
!  ncell:     total number of inversion grid cells to check valid range of cell indices in file
!  filename:  file name to use
!  lu:        file unit
!  errmsg:    error message
!
   subroutine createModelValuesFromFileDataModelSpaceInfo(this,propset,ncell,filename,lu,errmsg)
      type (data_model_space_info) :: this
      type (property_set) :: propset
      integer :: ncell
      character(len=*) :: filename
      integer :: lu
      type (error_message) :: errmsg
      character(len=500) :: line
      integer :: ios,iline
      integer :: j,jcell,n
      logical :: line_model_values_found
      character(len=43) :: myname = 'createModelValuesFromFileDataModelSpaceInfo'
   !
      iline = 0
      this%ncell = ncell
   !
      open(unit=lu,file=trim(filename),status='old',form='formatted',action='read',iostat=ios)
      if (ios /= 0) then
         call add(errmsg,2,"could not open file <"//trim(filename)//">",myname)
         return
      endif
   !
   ! parse through whole file until line "MODEL VALUES" and start reading in block from there
   !
      line_model_values_found = .false.
      do while (ios == 0)
         read(lu,"(a)",iostat=ios) line; iline = iline+1
         if(ios /= 0) exit
         if (trim(line) == "MODEL VALUES") then
            line_model_values_found = .true.
            exit
         endif
      enddo
      if (.not. line_model_values_found) then
         call add(errmsg,2,"No line MODEL VALUES in file",myname)
         return
      endif
   !
   !  read material properties of model space and check if valid
   !
      read(lu,"(a)",iostat=ios) line; iline = iline+1
      read(line,*,iostat=ios) this%nprop
      allocate(this%prop(this%nprop))
      read(line,*,iostat=ios) this%nprop,this%prop
      do j = 1,this%nprop
         if (.not. isValidNamePropertySet(propset,trim(this%prop(j)))) then
            call add(errmsg,2,"Invalid material property in dmsp file",myname)
         endif
      enddo
      call sortNamesPropertySet(propset,this%prop)
   !
   !  set attributes of model values
   !
      this%nmval = this%nprop*this%ncell
      allocate(this%prop_mval(this%nmval))
      allocate(this%cell_mval(this%nmval))
      n = 0
      do j = 1,this%nprop
         do jcell = 1,this%ncell
            n = n + 1
            this%prop_mval(n) = this%prop(j)
            this%cell_mval(n) = jcell
         enddo
      enddo
      close(lu)
   end subroutine createModelValuesFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!  deallocate data_model_space_info object
!
   subroutine deallocateDataModelSpaceInfo(this)
      type (data_model_space_info) :: this
      if(allocated(this%comp_in_dmsp)) deallocate(this%comp_in_dmsp)
      if(allocated(this%paths)) deallocate(this%paths)
      if(allocated(this%ncomp_path)) deallocate(this%ncomp_path)
      if(allocated(this%comp_path)) deallocate(this%comp_path)
      if(allocated(this%nfreq_path)) deallocate(this%nfreq_path)
      if(allocated(this%ifreq_path)) deallocate(this%ifreq_path)
      if(allocated(this%maskfreq_path)) deallocate(this%maskfreq_path)
      if(allocated(this%nfreq_path)) deallocate(this%nfreq_path_ori)
      if(allocated(this%ifreq_path)) deallocate(this%ifreq_path_ori)
      if(allocated(this%maskfreq_path)) deallocate(this%wfreq_path_ori)
   !
      if(allocated(this%pathidx)) deallocate(this%pathidx)
      if(allocated(this%offset_path)) deallocate(this%offset_path)
      if(allocated(this%evid)) deallocate(this%evid)
      if(allocated(this%netstaname)) deallocate(this%netstaname)
      if(allocated(this%comp)) deallocate(this%comp)
      if(allocated(this%ifreq)) deallocate(this%ifreq)
      if(allocated(this%noise)) deallocate(this%noise)
      if(allocated(this%prop)) deallocate(this%prop)
      if(allocated(this%prop_mval)) deallocate(this%prop_mval)
      if(allocated(this%cell_mval)) deallocate(this%cell_mval)
   end subroutine deallocateDataModelSpaceInfo
!------------------------------------------------------------------------
!  get total number of data samples
!
   function getNdataDataModelSpaceInfo(this) result(ndata)
      type (data_model_space_info), intent(in) :: this
      integer :: ndata
      ndata = this%ndata
   end function getNdataDataModelSpaceInfo
!------------------------------------------------------------------------
!  get total number of model parameters
!
   function getNmvalDataModelSpaceInfo(this) result(nmval)
      type (data_model_space_info), intent(in) :: this
      integer :: nmval
      nmval = this%nmval
   end function getNmvalDataModelSpaceInfo
!------------------------------------------------------------------------
!  get total number of inversion grid cells
!
   function getNcellDataModelSpaceInfo(this) result(ncell)
      type (data_model_space_info), intent(in) :: this
      integer :: ncell
      ncell = this%ncell
   end function getNcellDataModelSpaceInfo
!------------------------------------------------------------------------
!  get number of model properties
!
   function getNpropDataModelSpaceInfo(this) result(nprop)
      type (data_model_space_info), intent(in) :: this
      integer :: nprop
      nprop = this%nprop
   end function getNpropDataModelSpaceInfo
!------------------------------------------------------------------------
!  get total number of data paths
!
   function getNpathDataModelSpaceInfo(this) result(npath)
      type (data_model_space_info), intent(in) :: this
      integer :: npath
      npath = this%npath
   end function getNpathDataModelSpaceInfo
!------------------------------------------------------------------------
!  get number of different components occurring in dmsp
!
   function getNcompInDataModelSpaceInfo(this) result(res)
      type (data_model_space_info), intent(in) :: this
      integer :: res
      res = this%ncomp_in_dmsp
   end function getNcompInDataModelSpaceInfo
!------------------------------------------------------------------------
!  get different components occurring in dmsp
!
   function getCompsInDataModelSpaceInfo(this) result(res)
      type (data_model_space_info), intent(in), target :: this
      character(len=char_len_comp), dimension(:), pointer :: res
      res => this%comp_in_dmsp
   end function getCompsInDataModelSpaceInfo
!------------------------------------------------------------------------
!  get paths contained in this data_model_space_info object
!
   function getPathsDataModelSpaceInfo(this) result(paths)
      type (data_model_space_info), intent(in), target :: this
      character(len=char_len_path), dimension(:,:), pointer :: paths
      paths => this%paths
   end function getPathsDataModelSpaceInfo
!------------------------------------------------------------------------
!  get all path indices as pointer to array
!
   function getArrayPathidxDataSamplesDataModelSpaceInfo(this) result(res)
      type (data_model_space_info), intent(in), target :: this
      integer, dimension(:), pointer :: res
      res => this%pathidx
   end function getArrayPathidxDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!  get event ID of i-th data sample
!  i index of data sample
!
   function getEvidDataSampleDataModelSpaceInfo(this,i) result(evid)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      character(len=char_len_evid) :: evid
      if (i < 1 .or. i > this%ndata) then
         evid = ''
      else
         evid = this%evid(i)
      endif
    end function getEvidDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!  get all event IDs as pointer to array
!
   function getArrayEvidDataSamplesDataModelSpaceInfo(this) result(evid)
      type (data_model_space_info), intent(in), target :: this
      character(len=char_len_evid), dimension(:), pointer :: evid
      evid => this%evid
   end function getArrayEvidDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!  get net.station name of i-th data sample
!  i index of data sample
!
   function getNetStanameDataSampleDataModelSpaceInfo(this,i) result(netstaname)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
      if (i < 1 .or. i > this%ndata) then
         netstaname = ''
      else
         netstaname = this%netstaname(i)
      endif
   end function getNetStanameDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!  get all net.station names as pointer to array
!
   function getArrayNetStanameDataSamplesDataModelSpaceInfo(this) result(netstaname)
      type (data_model_space_info), intent(in), target :: this
      character(len=char_len_sta+char_len_netcode+1), dimension(:), pointer :: netstaname
      netstaname => this%netstaname
   end function getArrayNetStanameDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!  get frequency indices of all data as pointer to array
!
   function getIfreqArrayDataModelSpaceInfo(this) result(ifreq)
      type (data_model_space_info), intent(in), target :: this
      integer, dimension(:), pointer :: ifreq
      ifreq => this%ifreq
   end function getIfreqArrayDataModelSpaceInfo
!-------------------------------------------------------------------------
!  get station component of i-th data sample
!  i index of data sample
!
   function getCompDataSampleDataModelSpaceInfo(this,i) result(comp)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      character(len=char_len_comp) :: comp
      if (i < 1 .or. i > this%ndata) then
         comp = ''
      else
         comp = this%comp(i)
      endif
   end function getCompDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!  get frequency index of i-th data sample
!  i index of data sample
!
   function getIfreqDataSampleDataModelSpaceInfo(this,i) result(ifreq)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      integer :: ifreq
      if (i < 1 .or. i > this%ndata) then
         ifreq = -1
      else
         ifreq = this%ifreq(i)
      endif
   end function getIfreqDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!  get imre of i-th data sample
!  i index of data sample
!
   function getImreDataSampleDataModelSpaceInfo(this,i) result(imre)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      character(len=2) :: imre
      if (i < 1 .or. i > this%ndata) then
         imre = ''
      else
         if (mod(i,2) == 1) then
            imre = 're'
         else
            imre = 'im'
         end if
      endif
   end function getImreDataSampleDataModelSpaceInfo
!----------------------------------------------------------------------
!  get noise level for data samples
!
   subroutine getNoiseDataSamplesDataModelSpaceInfo(this,noise)
      type (data_model_space_info), intent(in), target :: this
      double precision, dimension(:), pointer :: noise
      noise => this%noise
   end subroutine getNoiseDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!  get propertes occuring in model space
!
   function getPropertiesDataModelSpaceInfo(this) result(prop)
      type (data_model_space_info), intent(in), target :: this
      character(len=char_len_par), dimension(:), pointer :: prop
      prop => this%prop
   end function getPropertiesDataModelSpaceInfo
!------------------------------------------------------------------------
!  get property of i-th model value
!  i index of model value
!
   function getPropertyModelValueDataModelSpaceInfo(this,i) result(prop)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      character(len=char_len_par) :: prop
      if (i < 1 .or. i > this%nmval) then
         prop = ''
      else
         prop = this%prop(i)
      endif
   end function getPropertyModelValueDataModelSpaceInfo
!------------------------------------------------------------------------
!  inversion grid cell index of i-th model parameter
!  i index of model parameter
!
   function getCellModelValueDataModelSpaceInfo(this,i) result(res)
      type (data_model_space_info), intent(in) :: this
      integer, intent(in) :: i
      integer :: res
      if (i < 1 .or. i > this%nmval) then
         res = -1
      else
         res = this%cell_mval(i)
      endif
   end function getCellModelValueDataModelSpaceInfo
!---------------------------------------------------------------------------------------
!  Mask data specified by path and frequency index
!
   subroutine setPathFrequencyMaskDataModelSpaceInfo(this,jpath,ifreq)
      type (data_model_space_info) :: this
      integer, intent(in) :: jpath,ifreq
      integer :: jf
      jf = findloc(this%ifreq_path_ori(:,jpath),ifreq,1)
      this%maskfreq_path(jf,jpath) = .false.
   end subroutine setPathFrequencyMaskDataModelSpaceInfo
!---------------------------------------------------------------------------------------
!  Extract a contiguous subset of paths from a dmspace object into a new one
!
   subroutine extractPathsDataModelSpaceInfo(this,first,last,that,errmsg)
      type (data_model_space_info) :: this,that
      integer :: first,last
      type (error_message):: errmsg
      integer :: npath,ndata,nprop,ncell,nmval,maxfreq,j,ida,ide
      character(len=30) :: myname = "extractPathsDataModelSpaceInfo"
   !
      if (first < 1 .or. first > this%npath) then
         call add(errmsg,2,"invalid first index of path selection",myname)
         return
      end if
      if (last < 1 .or. last > this%npath .or. last < first) then
         call add(errmsg,2,"invalid last index of path selection",myname)
         return
      end if
   !
      npath = last-first+1
      if (npath == this%npath) then
         that = this
         return
      end if
   !
   !  deal with paths
   !
      that%npath = npath
      maxfreq = size(this%ifreq_path,1)
      allocate(that%paths(2,npath))
      allocate(that%offset_path(npath))
      allocate(that%ncomp_path(npath))
      allocate(that%comp_path(maxcomp_dmsp,npath))
      allocate(that%nfreq_path(npath))
      allocate(that%ifreq_path(maxfreq,npath))
      that%paths = this%paths(:,first:last)
      that%offset_path = this%offset_path(first:last)
      that%ncomp_path = this%ncomp_path(first:last)
      that%comp_path = this%comp_path(:,first:last)
      that%nfreq_path = this%nfreq_path(first:last)
      that%ifreq_path = this%ifreq_path(:,first:last)
   !
   !  deal with data samples
   !
      ndata = 0
      do j = first,last
         ndata = ndata+2*this%ncomp_path(j)*this%nfreq_path(j)
      end do
      this%ndata = ndata
      allocate(that%evid(ndata))
      allocate(that%netstaname(ndata))
      allocate(that%comp(ndata))
      allocate(that%ifreq(ndata))
      allocate(that%noise(ndata))
      ida = this%offset_path(first)+1
      ide = ida+ndata-1
      that%evid = this%evid(ida:ide)
      that%netstaname = this%netstaname(ida:ide)
      that%comp = this%comp(ida:ide)
      that%ifreq = this%ifreq(ida:ide)
      that%noise = this%noise(ida:ide)
   !
   !  deal with model values
   !
      nprop = this%nprop
      nmval = this%nmval
      ncell = this%ncell
      allocate(that%prop(nprop))
      allocate(that%prop_mval(nprop*ncell))
      allocate(that%cell_mval(nprop*ncell))
      that%nprop = nprop
      that%nmval = nmval
      that%ncell = ncell
      that%prop = this%prop
      that%prop_mval = this%prop_mval
      that%cell_mval = this%cell_mval
   end subroutine extractPathsDataModelSpaceInfo
!--------------------------------------------------------------------------------------
!  iterate over paths
!  return evid and netstaname of next path. optionally, return array of all data sample indices of that path,
!  all different components of that path and all frequenc indices of next path
!  evid event ID of next path
!  netstaname net.station name of next path
!  indx optional pointer to array of indices of data samples of the next path
!  comp optional pointer to array of component names of the next path
!  ifreq optional pointer to array of frequency indices of the next path.
!  ipath_start optional integer indicating path index at which the loop should start (only considered at first call)
!  ipath_end optional integer indicating path index at which the loop should end (only considered at first call)
!  next logical indicating whether there is a next path or not (if not, the return values are dummies and not meaningful)
!
   function nextPathDataModelSpaceInfo(this,evid,netstaname,indx,comp,ifreq) result(next)
      type (data_model_space_info) :: this
      character(len=char_len_evid) :: evid
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
      integer, dimension(:), allocatable, optional :: indx
      character(len=char_len_comp), dimension(:), allocatable, optional :: comp
      integer, dimension(:), allocatable, optional :: ifreq
      logical :: next
      ! save
      integer, save :: jpath = 0
      integer :: nsamp,ncomp,nfreq,i,nimre
   !
      jpath = jpath+1
   !
   !  check end of loop
   !
      if (jpath > this%npath) then
         next = .false.
         jpath = 0
         return
      else
         next = .true.
      end if
   !
   !  set variables
   !
      nimre = 2
      evid = this%paths(1,jpath)
      netstaname = this%paths(2,jpath)
      if (present(indx)) then
         nsamp = this%ncomp_path(jpath)*this%nfreq_path(jpath)*nimre
         allocate(indx(nsamp))
         indx = [(i, i = this%offset_path(jpath)+1,this%offset_path(jpath)+nsamp)]
      endif
      if (present(comp)) then
         ncomp = this%ncomp_path(jpath)
         allocate(comp(ncomp))
         comp = this%comp_path(1:ncomp,jpath)
      end if
      if (present(ifreq)) then
         nfreq = this%nfreq_path(jpath)
         allocate(ifreq(nfreq))
         ifreq = this%ifreq_path(1:nfreq,jpath)
      end if
   end function nextPathDataModelSpaceInfo
!------------------------------------------------------------------------
!  write data model space info to file
!
   subroutine writeDataModelSpaceInfo(this,lu,filename,first,last,errmsg)
      type (data_model_space_info) :: this
      integer :: lu
      character(len=*) :: filename
      integer :: first,last
      type (error_message) :: errmsg
      integer :: ios,j,i,nf,nc
      double precision :: weight_path
      character(len=max_length_string) :: fmt
      character(len=23) :: myname = 'writeDataModelSpaceInfo'
   !
      open(unit=lu,file=trim(filename),status='UNKNOWN',action='WRITE',iostat=ios)
      if (ios /= 0) then
         call add(errmsg,2,'Cannot open file '//trim(filename),myname)
         return
      end if
      write(lu,'(a)') 'MODEL VALUES'
      write(fmt,'(a,i1,a)') '(i1,1x,',this%nprop,'a6)'
      write(lu,trim(fmt)) this%nprop,this%prop
      write(lu,'(a)') 'DATA SAMPLES'
      write(lu,'(a,1x,a)') 'WEIGHTING','BY_PATH_AND_FREQUENCY'
      write(lu,'(a,1x,a)') 'MASKING','BY_FREQUENCY'
      write(lu,'(a,1x,a)') 'COMPONENTS','SPECIFIC'
      write(lu,'(a,1x,a)') 'FREQUENCIES','SPECIFIC'
      write(lu,'(i6)')  last-first+1
      weight_path = 1.d0
      do j = first,last
         nf = this%nfreq_path_ori(j)
         nc = this%ncomp_path(j)
         write(lu,'(a,1x,a,e15.3)') this%paths(1,j),this%paths(2,j),weight_path
         write(fmt,'(a,i1,a)') '(i1,1x,',nc,'a2)'
         write(lu,trim(fmt)) nc,(this%comp_path(i,j),i = 1,nc)
         write(fmt,'(a,i3,a)') '(i3,1x',nf,'i3)'
         write(lu,trim(fmt)) nf,(this%ifreq_path_ori(i,j),i = 1,nf)
         write(fmt,'(a,i3,a)') '(i3,1x',nf,'e15.3)'
         write(lu,trim(fmt)) nf,(this%wfreq_path_ori(i,j),i = 1,nf)
         write(fmt,'(a,i3,a)') '(i3,1x',nf,'l5)'
         write(lu,trim(fmt)) nf,(this%maskfreq_path(i,j),i = 1,nf)
      end do
      close(lu)
   end subroutine writeDataModelSpaceInfo
!------------------------------------------------------------------------
!  print all data samples and model values
!
   subroutine printDataModelSpaceInfo(this)
      type (data_model_space_info) :: this
      integer :: i
      write(*,*) "#############################################################"
      write(*,*) "DATA-MODEL-SPACE-INFO"
      write(*,*) this%ndata," data samples:"
      write(*,*) "evid    netstaname   comp    ifreq    noise    imre"
      write(*,*) "-------------------------------------------------------------"
      do i = 1,this%ndata,2
         write(*,'(3a,i5,f8.3)') trim(this%evid(i)),trim(this%netstaname(i)),trim(this%comp(i)),&
                                 this%ifreq(i),this%noise(i),"    re"
         write(*,'(3a,i5,f8.3)') trim(this%evid(i+1)),trim(this%netstaname(i+1)),trim(this%comp(i+1)),&
                                 this%ifreq(i+1),this%noise(i+1),"    im"
      enddo
      write(*,*) ""
      write(*,*) this%nmval
      write(*,*) "prop     cell "
      write(*,*) "-------------------------------------------------------------"
      do i = 1,this%nmval
         write(*,'(a,i10)') trim(this%prop_mval(i)),this%cell_mval(i)
      end do
      write(*,*) "#############################################################"
   end subroutine printDataModelSpaceInfo
!------------------------------------------------------------------------
end module dataModelSpaceInfo
