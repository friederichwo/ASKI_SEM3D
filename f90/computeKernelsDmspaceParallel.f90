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
!
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
program computeKernelsDmspaceParallel
   use hdf5
   use inversionBasics
   use iterationStepBasics
   use dataModelSpaceInfo
   use specfem3dInversionGrid
   use rayInversionGrid
   use specfem3dKernelWavefield
   use rayKernelWavefield
   use spectralWaveformKernel
   use propertySet
   use argumentParser
   use inputParameter
   use string
   use mpiSupport
   use hdfWrapper
   use anyRankRealArray
   use anyRankIntegerArray
   use errorMessage
   use globalMpiInfo
   use globalHdfInfo
   implicit none

   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (data_model_space_info) :: dmspace
   type (any_rank_real_array) :: arra
   class (kernel_wavefield), dimension(:), allocatable :: kwev
   class (kernel_wavefield), allocatable :: kwsta
   type (spectral_waveform_kernel) :: skernel
   type (mpi_support) :: mpisup
   class (inversion_grid), allocatable :: invgrid
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar_inv,inpar_iter
   type (seismic_event) :: event
   type (seismic_event_list), pointer :: evlist
   type (seismic_network), pointer :: statlist
   integer :: npath,nev
   integer :: jf,jev,ncell,ncell_all,nfreq,nwp,nwp_all,offcell,offwp,prank
   integer (kind=8) :: fwk,fpet
   integer, dimension(:), pointer  :: ifreq,ifreq_valid
   integer, dimension(:), allocatable :: ifreq_dmsp
   real, dimension(:), pointer  :: d
   double complex, dimension(:,:), pointer :: ustr,u
   double precision, dimension(:), pointer :: petu
   double precision :: df,ufmdata,wtanf,wtend,wtf,wtpa,wtini,petsta,rsctap
   logical :: on_wp,do_hyperslab
   character(len=max_length_string) :: main_parfile,dmspace_file,iterpath,kd_filebase,path_dmsp,path_dft_data
   character(len=max_length_string) :: skernel_filebase,kgt_filebase,main_file_path,filename,timing_file,path_phase_end_times
   character(len=max_length_string) :: property_set_name,forward_method,main_path_inversion,invgrid_type
   character(len=char_len_evid) :: evid
   character(len=char_len_sta+char_len_netcode+1) :: netstaname,netstaname_prev
   character(len=char_len_comp), dimension(:), allocatable :: comp_path
   character(len=char_len_par), dimension(:), allocatable :: prop_out,prop_req
   character(len=char_len_par), dimension(:), pointer :: prop_dmsp
   character(len=29) :: myname = 'computeKernelsDmspaceParallel'
!-------------------------------------------------------------------------------------
!  initialise MPI
!
   call new(mpisup)
   myrank = .myrank.mpisup
   numtasks = .numtasks.mpisup
!
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
!  read commandline arguments
!
   call init(ap,myname,"Compute a set of kernels characterized by the data model space file")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-wp',.false.,"if set, then plain kernel values on WAVEFIELD POINTS are produced. Otherwise"//&
      " (if not set), pre-integrated kernels on inversion grid cells are computed")
   call addOption(ap,'-ifreq',.true.,"frequency index at which the kernel should be calculated. ",'ival','1')
   call parse(ap)
   if (myrank == prank) then
      call document(ap)
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap); call usage(ap)
         call add(errmsg,2,'Command line could not be read successfully',myname)
         goto 1
      endif
   endif
   main_parfile = ap.sval.'main_parfile'
   on_wp = ap.optset.'-wp'
! ------------------------------------------------------------------------
!  setup basics
!
   prank = 21
   wtini = MPI_WTIME()
   do_hyperslab = .true.
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
!  other stuff to do beforehand
!
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter => getInputParameterIterationStepBasics(iterbasics)
   main_path_inversion = inpar_inv.sval.'MAIN_PATH_INVERSION'
   propset => getPropertySetInversionBasics(invbasics)
   evlist => getEventListInversionBasics(invbasics)
   statlist => getStationListInversionBasics(invbasics)
   ufmdata = .ufmdata.invbasics
   iterpath = .iterpath.invbasics
   path_dmsp = trim(iterpath)//trim(inpar_inv.sval.'PATH_DMSPACE')
   path_dft_data = main_path_inversion+sval(inpar_inv,'PATH_MEASURED_DATA')
   path_phase_end_times = main_path_inversion+sval(inpar_inv,'PATH_PHASE_END_TIMES')
   dmspace_file = trim(path_dmsp)//'ASKI_dmspace'
   timing_file = trim(path_dft_data)//'timing_of_data.hdf'
   property_set_name = .setname.propset
   forward_method = inpar_inv.sval.'FORWARD_METHOD'
   rsctap = inpar_iter.dval.'RESIDUAL_SCATTERING_TIME_KERNEL_TAPER'
   nev = .nev.evlist
   ifreq_valid => getFrequencyIndicesIterationStepBasics(iterbasics)
   if (ap.optset.'-ifreq') then
      allocate(ifreq(1))
      ifreq(1) = ap.ival.'-ifreq'
      nfreq = 1
   else
      ifreq => getFrequencyIndicesIterationStepBasics(iterbasics)
      nfreq = .nf.iterbasics
   endif
   call dealloc(ap)
!
!  set the type of kernel displacement and kernel Green tensor
!  according to forward method
!
   if (equalString(forward_method,'specfem3d')) then
      allocate(specfem3d_kernel_wavefield :: kwev(nev))
      allocate(specfem3d_kernel_wavefield :: kwsta)
   else if (equalString(forward_method,'fm3d')) then
      allocate(ray_kernel_wavefield :: kwev(nev))
      allocate(ray_kernel_wavefield :: kwsta)
   else
      call add(errmsg,2,'Forward method '//trim(forward_method)//' not implemented',myname)
      goto 1
   end if
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(iterpath,do_hyperslab,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
!
   df = .df.invbasics
   nwp_all = invgrid%getNwpAll()
   nwp = invgrid%getNwpLocal()
   offwp = invgrid%getWpOffset()
   ncell_all = invgrid%getNcellAll()
   ncell = invgrid%getNcellLocal()
   offcell = invgrid%getCellOffset()
!
!  read data model space info
!
   call createFromFileDataModelSpaceInfo(dmspace,evlist,statlist,&
         ifreq_valid,propset,ncell_all,trim(dmspace_file),1,errmsg)
   if (.level.errmsg == 2) goto 1
   if(.ndata.dmspace == 0 .or. .nmval.dmspace == 0) then
       call add(errmsg,2,"Data space or model space is empty",myname)
       goto 1
   endif
   npath = .npath.dmspace
!
!  open file with timing of data and synthetics at stations
!
   call openFileParallelAccessHDFWrapper(timing_file,fpet,errmsg)
   if (.level.errmsg == 2) return
!
!  get model properties to be inverted for from dmsp info
!  and check consistency with correlaton matrix
!
   prop_dmsp => getPropertiesDataModelSpaceInfo(dmspace)
   if (.not. checkCorrmatPropertySet(propset,prop_dmsp,errmsg)) goto 1
   call getRequiredKernelsByNamePropertySet(propset,prop_req)
   call getOutputKernelsByNamePropertySet(propset,prop_out)
!
!  initialize computation of spectral waveform kernel
!
   call initSpectralWaveformKernel(skernel,ncell,ncell_all,offcell,nwp,nwp_all,offwp,prop_req,prop_out,df,rsctap,on_wp)
! ------------------------------------------------------------------------
!  base of all kernel displacement filenames
!
   kd_filebase = iterpath+(inpar_inv.sval.'PATH_KERNEL_DISPLACEMENTS')
   kgt_filebase = iterpath+(inpar_inv.sval.'PATH_KERNEL_GREEN_TENSORS')
   main_file_path = iterpath+(inpar_inv.sval.'PATH_ASKI_MAIN_FILES')
   if(on_wp) then
      skernel_filebase = trim(iterpath)//trim(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_ON-WP_'//trim(property_set_name)//'_'
   else
      skernel_filebase = trim(iterpath)//trim(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_'//trim(property_set_name)//'_'
   endif
!
!  write some informative output
!
   if (myrank == prank) then
       write(*,*) "will now compute spectral '"//trim(property_set_name)//&
            "' kernels of properties ",prop_out//','," for ",npath," path(s)"
       write(*,*)"ASKI main files will be read from: '",trim(main_file_path)//"'"
       write(*,*)"Kernel displacements will be read from files: '",trim(kd_filebase)//"EVENTID'"
       write(*,*)"Kernel Green tensors will be read from files: '",trim(kgt_filebase)//"NETSTANAME'"
       if(on_wp) then
          write(*,*) "will compute plain kernel values on wavefield points"
       else
          write(*,*) "will compute pre-integrated kernel values on inversion grid cells"
       end if
   endif
!
!  start loop over frequencies
!
   do jf = 1,nfreq
      wtf = MPI_WTIME()
      if (myrank == prank) print*,'======================= Work on frequency: ',jf,' ====================='
   !
   !  create spectral waveform kernel HDF output file
   !
      write(filename,'(a,i6.6,a)') trim(skernel_filebase)//'jf',ifreq(jf),'.hdf'
      call createFileParallelAccessHDFWrapper(filename,fwk,errmsg)
      if(.level.errmsg == 2) goto 1
      if (myrank == prank) print *,'created file : ',trim(filename)
      call writeMetaSpectralWaveformKernel(skernel,propset,fwk,ifreq(jf),npath,errmsg)
      if(.level.errmsg == 2) goto 1
   !
   !  start kernel computation
   !  first read kernel displacement spectra for all events into memory
   !  to avoid repeated Green tensor reads for stations
   !
      wtanf = MPI_WTIME()
      do while(nextEventSeismicEventList(evlist,event,icall=jev))
         evid = .evid.event
         call kwev(jev)%readDisplacement(nwp,ifreq(jf),df,kd_filebase,evid,do_hyperslab,errmsg)
         if (.level.errmsg == 2) goto 1
         call kwev(jev)%readDisplacementPhaseEndTime(nwp,path_phase_end_times,evid,do_hyperslab,errmsg)
         if (.level.errmsg == 2) goto 1
         if (myrank == prank) print *,'Reading wavefields for event: ',trim(evid)
      end do
      wtend = MPI_WTIME()
      if (myrank == prank) then
         print *,'Time for reading event displacements on rank: ',prank,wtend-wtanf
      endif
   !
   !  loop over paths
   !  assume sorting according to stations (see createDmspaceFile.py)
   !  if the frequency is not requested for this path go to next one
   !  the sensitivity kernel for this path will not be present in the file
   !
      netstaname_prev = ""
      do while(nextPathDataModelSpaceInfo(dmspace,evid,netstaname,comp=comp_path,ifreq=ifreq_dmsp))
         if(.not. any(ifreq_dmsp == ifreq(jf))) then
            if (myrank == prank) then
               write(*,*) 'Frequency ',ifreq(jf),' not requested for path ',trim(evid),'_',trim(netstaname)
            end if
            deallocate(comp_path,ifreq_dmsp)
            cycle
         end if
         wtpa = MPI_WTIME()
         if (myrank == prank) print *,'================= Work on path: ',trim(evid)," ",trim(netstaname),' ================='
      !
      ! read kernel Green tensor only
      !
         if (trim(netstaname) /= trim(netstaname_prev)) then
            wtanf = MPI_WTIME()
            call kwsta%deallocGreen()
            call kwsta%readGreenTensor(nwp,ifreq(jf),df,kgt_filebase,netstaname,comp_path,do_hyperslab,errmsg)
            if (.level.errmsg == 2) goto 1
            call kwsta%readGreenTensorPhaseEndTime(nwp,path_phase_end_times,netstaname,do_hyperslab,errmsg)
            netstaname_prev = netstaname
            wtend = MPI_WTIME()
            if (myrank == prank) then
               print *,'Time for reading kernel Green tensor on rank: ',prank,wtend-wtanf
            endif
         endif
      !
         wtanf = MPI_WTIME()
      !
      !  read phase end time of data rsp. synthetics for this path
      !
         call readArrayHDFWrapper(fpet,trim(evid)//"_"//trim(netstaname),arra,errmsg,xferprp = hdf_xferprp)
         if (.level.errmsg == 2) return
         d => arra%get1d()
         petsta = d(4); deallocate(d)
      !
      !  get index of desired event in list
      !
         if (.not. searchEventidSeismicEventList(evlist,evid,iev=jev)) then
            call add(errmsg,2,'Event '//trim(evid)//' not found in event list',myname)
            goto 1
         end if
      !
      !  get displacements from kwev(jev)
      !  for ray wavefields we need a special routine to return the full 3D-array of displacements
      !  and to activate the associate routine in rayKernelWavefield
      !
         select type (kwev)
         class is (specfem3d_kernel_wavefield)
            u => kwev(jev)%getDisplacement()
            ustr => kwev(jev)%getStrain()
            petu => kwev(jev)%getDisplacementPhaseEndTime()
         class is (ray_kernel_wavefield)
            continue
         class default
            call add(errmsg,2,'kernel_wavefield-type not implemented',myname)
            return
         end select
      !
      !  associate displacements of kwsta with that of kwev(j)
      !
         select type (kwsta)
         class is (specfem3d_kernel_wavefield)
            call kwsta%associateDisplacement(u,ustr)
            call kwsta%associateDisplacementPhaseEndTime(petu)
         class is (ray_kernel_wavefield)
            continue
         class default
            call add(errmsg,2,'kernel_wavefield-type not implemented',myname)
            return
         end select
      !
      !  compute kernels honouring property correlation
      !
         call computeSpectralWaveformKernel(skernel,kwsta,invgrid,ufmdata,propset,comp_path,ifreq(jf),petsta,errmsg)
         if (.level.errmsg == 2) goto 1
      !
         wtend = MPI_WTIME()
         if (myrank == prank) then
            print *,'Time for computing sensitivity kernel on rank: ',prank,wtend-wtanf
         endif
         wtanf = MPI_WTIME()
      !
         call writeSpectralWaveformKernel(skernel,fwk,evid,netstaname,errmsg)
         if (.level.errmsg == 2) goto 1
         call deallocateKernelSpectralWaveformKernel(skernel)
         deallocate(comp_path,ifreq_dmsp)
      !
         wtend = MPI_WTIME()
         if (myrank == prank) then
            print *,'Time for writing sensitivity kernel on rank: ',prank,wtend-wtanf
            print *,'Time for this path on rank: ',prank,wtend-wtpa
            print *,'Time since start of frequency loop on rank: ',prank,wtend-wtf
         endif
      !
         call barrier(mpisup)               ! synchronize processes after a path
      enddo                                 ! next path
   !
      do jev = 1,nev
         call kwev(jev)%dealloc()           ! deallocate event kernel wavefields
      end do
      call closeFileHDFWrapper(fwk,errmsg)  ! close sensitivity kernel file for this frequency
      if (.level.errmsg == 2) goto 1
   !
      wtend = MPI_WTIME()
      if (myrank == prank) then
         print *,'Total time for frequency loop on rank: ',prank,jf,wtend-wtf
      endif
   enddo                                    ! next frequency
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call deallocateMetaSpectralWaveformKernel(skernel)
   deallocate(kwev)
   deallocate(prop_req,prop_out)
   call dealloc(dmspace)
   call invgrid%dealloc()
   deallocate(invgrid)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(mpisup)
   call dealloc(errmsg)
!
   wtend = MPI_WTIME()
   if (myrank == prank) then
      print *,'Total time for run: ',wtend-wtini
   endif
!
!  error handling
!
 1 if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   endif
end program computeKernelsDmspaceParallel
