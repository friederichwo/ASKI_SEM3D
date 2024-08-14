!----------------------------------------------------------------------------
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
!
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
program decomposeDmspace
   use inversionBasics
   use iterationStepBasics
   use dataModelSpaceInfo
   use smartUtils
   use argumentParser
   use errorMessage
   implicit none

   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (data_model_space_info) :: dmspace
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar,inpar_inv
   type (seismic_event_list), pointer :: evlist
   type (seismic_network), pointer :: statlist
   integer, dimension(:), pointer  :: ifreq,id
   integer :: nproc,npath,off,nloc,ip,ncell
   logical :: use_masked_dmspace
   character(len=max_length_string) :: main_parfile,dmspace_file,filename,path_dmsp,iterpath
   character(len=16) :: myname = 'decomposeDmspace'
!--------------------------------------------------------------------------
   call new(errmsg,myname)
! ------------------------------------------------------------------------
!  read commandline arguments
!
   call init(ap,myname,"Decompose data-model-space for parallelization of solveCglsKernelSystem")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-masked',.false.,'Use the masked dmspace file')
   call parse(ap)
   call document(ap)
   if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); goto 1; endif
   main_parfile = ap.sval.'main_parfile'
   use_masked_dmspace = ap.optset.'-masked'
   call dealloc(ap)
! ------------------------------------------------------------------------
!  setup basics
!
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
   inpar_inv => getInputParameterInversionBasics(invbasics)
   iterpath = .iterpath.invbasics
   path_dmsp = trim(iterpath)//trim(inpar_inv.sval.'PATH_DMSPACE')
   dmspace_file = trim(path_dmsp)//'ASKI_dmspace'
   if (use_masked_dmspace) dmspace_file = trim(path_dmsp)//'ASKI_dmspace_masked'
   propset => getPropertySetInversionBasics(invbasics)
   evlist => getEventListInversionBasics(invbasics)
   statlist => getStationListInversionBasics(invbasics)
   inpar => getInputParameterIterationStepBasics(iterbasics)
   ifreq => getFrequencyIndicesIterationStepBasics(iterbasics)
   nproc = inpar.ival.'NPROC'
   id => ivecp(inpar_inv,'INVGRID_DIMENSIONS',3)
   ncell = id(1)*id(2)*id(3)
   deallocate(id)
!
!  read data model space info
!
   call createFromFileDataModelSpaceInfo(dmspace,evlist,statlist,&
         ifreq,propset,ncell,trim(dmspace_file),1,errmsg)
   if (.level.errmsg == 2) goto 1
   if(.ndata.dmspace == 0 .or. .nmval.dmspace == 0) then
       call add(errmsg,2,"Data space or model space is empty",myname)
       goto 1
   endif
   npath = .npath.dmspace
   print *,'Number of paths: ',npath
!
!  distribute paths over parallel processes
!  and write individual dmsp files
!
   do ip = 0,nproc-1
      call shareLinearArraySmartUtils(npath,nproc,ip,nloc,off)
      print *,'paths for rank: ',ip,nloc
      write(filename,"(a,'.',i3.3)") trim(path_dmsp)//'ASKI_dmspace',ip
      call writeDataModelSpaceInfo(dmspace,1,trim(filename),off+1,off+nloc,errmsg)
      if (.level.errmsg == 2) goto 1
   end do
   call dealloc(dmspace)
!
!  error handling
!
 1 if (.level.errmsg == 2) then
      call print(errmsg)
   endif
end program decomposeDmspace
