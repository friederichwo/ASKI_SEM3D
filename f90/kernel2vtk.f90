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
program kernel2vtk
   use inversionBasics
   use iterationStepBasics
   use propertySet
   use specfem3dInversionGrid
   use spectralWaveformKernel
   use vtkFile
   use argumentParser
   use string
   use fileUnitHandler
   use errorMessage
   use mpiSupport
   use hdfWrapper
   use globalHdfInfo
   use globalMpiInfo

   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (spectral_waveform_kernel) :: kernel
   type (vtk_info) :: vtkinfo
   type (mpi_support) :: mpisup
   class (inversion_grid), allocatable :: invgrid
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar_inv,inpar_iter
   type (seismic_event_list), pointer :: evlist
   type (seismic_network), pointer :: statlist
   integer (kind=8) :: fid
   integer :: j,k,ierr,jf,vtkmode
   integer :: nprop,nfreq,ncomp
   integer, dimension(:), pointer :: all_ifreq,ifreq
   logical :: kernels_on_wp
   logical :: normalize_kernels
   real :: df,knorm
   real, dimension(:), pointer :: pre,pim
   real, dimension(:), allocatable :: nkre,nkim
   double precision :: vtk_scale
   character(len=max_length_string) :: main_parfile,str
   character(len=max_length_string), dimension(:), pointer :: str_vec
   character(len=char_len_evid) :: evid
   character(len=char_len_sta+char_len_netcode+1) :: staname
   character(len=char_len_par), dimension(:), allocatable :: prop
   character(len=char_len_comp), dimension(:), allocatable :: comp
   character(len=max_length_string) :: kernel_file,kernel_filebase,vtk_file_base,vtk_file_data_name,vtk_path
   character(len=max_length_string) :: vtk_format,vtk_geometry_type,invgrid_type
   character(len=10) :: myname = 'kernel2vtk'
   nullify(str_vec,all_ifreq,ifreq,pre,pim)
!-----------------------------------------------------------------------------
!  initialise MPI
!
   call new(mpisup)
   myrank = .myrank.mpisup
   numtasks = .numtasks.mpisup
   call new(errmsg,myname)
!------------------------------------------------------------------------
!  command line processing
!
   call init(ap,myname,"Writes spectral sensitivity kernels as vtk files for specific paths, properties, "//&
        "components and frequencies. Can either handle pre-integrated kernel values on inversion grid cells, or "//&
        "original kernel values on wavefield points (flag -wp).")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-evid',.true.,"(mandatory) defines the event id of the one path (must belong to an event in"//&
        " main event list)",'sval','')
   call addOption(ap,'-staname',.true.,"(mandatory) defines the netcode.station name of the one path (must belong to a "//&
         "netcode.station in main station list",'sval','')
   call addOption(ap,'-comp',.true.,"(mandatory) vector of receiver components for which the kernel should be "//&
        "computed; see valid components in constants.f90",'svec','')
   call addOption(ap,'-prop',.true.,"(mandatory) vector of property names for which the kernel should be "//&
        "computed",'svec','')
   call addOption(ap,'-ifreq',.true.,"explicit vector of frequency indices at which the wavefield output should "//&
        "be extracted. Exactly one of options -ifreq , -all_ifreq must be set",'ivec','')
   call addOption(ap,'-all_ifreq',.false.,"if set, all frequency indices are used. Exactly one of options "//&
        "-ifreq , -all_ifreq must be set")
   call addOption(ap,"-wp",.false.,"if set, the original kernels on the WAVEFIELD POINTS are produced. If not "//&
        "set, pre-integrated kernel values on inversion grid cells are produced")
   call addOption(ap,"-norm",.false.,"if set, kernel absolute value is normalized to 1.0")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   if ((ap.optset.'-all_ifreq') .eqv. (ap.optset.'-ifreq')) then
      call add(errmsg,2,"exactly ONE of the options -ifreq and -all_ifreq must be set!",myname)
      call usage(ap)
      goto 1
   end if
   if(.not.( (ap.optset.'-evid') .and. (ap.optset.'-staname') .and. (ap.optset.'-comp') .and. (ap.optset.'-prop'))) then
      call add(errmsg,2,"all of the options -evid, -staname, -comp, -prop must be set!",myname)
      goto 1
   end if
!
   str_vec => getSvecArgumentParser(ap,'-comp')
   ncomp = size(str_vec)
   allocate(comp(ncomp))
   do j = 1,ncomp
      comp(j) = trim(str_vec(j))
      if (.not. any(valid_components == trim(comp(j)))) then
         call add(errmsg,2,"component "//trim(comp(j))//" is invalid",myname)
         goto 1
      endif
   enddo
   deallocate(str_vec)
!
   kernels_on_wp = ap.optset.'-wp'
   normalize_kernels = ap.optset.'-norm'
   main_parfile = ap.sval.'main_parfile'
! --------------------------------------------------------------------------------------------------
!  setup inversion basics
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  setup iteration step basics
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter =>  getInputParameterIterationStepBasics(iterbasics)
   propset => getPropertySetInversionBasics(invbasics)
   evlist => getEventListInversionBasics(invbasics)
   statlist => getStationListInversionBasics(invbasics)
   all_ifreq => getFrequencyIndicesIterationStepBasics(iterbasics)
   vtk_format = trim(inpar_inv.sval.'DEFAULT_VTK_FILE_FORMAT')
   vtk_scale = inpar_iter.dval.'VTK_COORDS_SCALING_FACTOR'
   vtk_geometry_type = inpar_iter.sval.'VTK_GEOMETRY_TYPE'
!
   str = ap.sval.'-evid'
   evid = str
   if (.not. searchEventidSeismicEventList(evlist,evid)) then
      call add(errmsg,2,"event ID '"//trim(evid)//"'is not contained in event list",myname)
      goto 1
   end if
!
   str = ap.sval.'-staname'
   staname = str
   if (.not. searchNetcodeStationNameSeismicNetwork(statlist,staname)) then
      call add(errmsg,2,"netcode.station name '"//trim(staname)//"'is not contained in station list",myname)
      goto 1
   end if
!
   str_vec => getSvecArgumentParser(ap,'-prop')
   nprop = size(str_vec)
   allocate(prop(nprop))
   do j = 1,nprop
      prop(j) = str_vec(j)
      if (.not. isValidNamePropertySet(propset,prop(j))) then
         call add(errmsg,2,"model property '"//trim(prop(j))//"' is invalid",myname)
         goto 1
      end if
   end do
   deallocate(str_vec)
!
   if (ap.optset.'-all_ifreq') then
      nfreq = size(all_ifreq)
      allocate(ifreq(nfreq))
      ifreq = all_ifreq
   else
      ifreq => getIvecArgumentParser(ap,'-ifreq')
      nfreq = size(ifreq)
      do j = 1,nfreq
         if (.not. any(all_ifreq == ifreq(j))) then
            call add(errmsg,2,"invalid frequency index ",myname)
            goto 1
         end if
      end do
   end if
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(.iterpath.invbasics,.false.,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
!
   if(kernels_on_wp) then
      kernel_filebase = trim(.iterpath.invbasics)//trim(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_ON-WP_'//trim(.setname.propset)//'_'
      vtkmode = 2
   else
      kernel_filebase = trim(.iterpath.invbasics)//trim(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_'//trim(.setname.propset)//'_'
      vtkmode = 1
   endif

   call initiateVtkFile(vtkinfo,invgrid,vtk_format,vtk_geometry_type,vtk_scale,vtkmode,errmsg)
   if (.level.errmsg == 2) goto 1
   vtk_path = trim(.iterpath.iterbasics)//trim(inpar_inv.sval.'PATH_VTK_FILES')//'sensitivity'
!
   do jf = 1,nfreq
      write(kernel_file,'(a,i6.6,a)') trim(kernel_filebase)//'jf',ifreq(jf),'.hdf'
      print *,'open kernel file: ',trim(kernel_file)
      call openFileRoHDFWrapper(kernel_file,fid,errmsg)
      if (.level.errmsg == 2) goto 1
      call readMetaSpectralWaveformKernel(kernel,fid,ifreq(jf),errmsg)
      if (.level.errmsg == 2) goto 1
      df = .df.kernel
      call readSpectralWaveformKernel(kernel,fid,evid,staname,errmsg)
      if (.level.errmsg == 2) goto 1
      do j = 1,ncomp
         do k = 1,nprop
            print *,'comp = ',comp(j),' prop = ',prop(k)
            call getValuesByComponentAndPropertySpectralWaveformKernel(kernel,comp(j),prop(k),pre,pim)
            if (.not. associated(pre)) then
               call add(errmsg,2,"no sensitivity kernel for "//trim(comp(j))//"_"//trim(prop(k))//" available",myname)
               goto 1
            end if
            write(vtk_file_base,"(a,4('_',a))") trim(vtk_path),trim(prop(k)),trim(evid),trim(staname),trim(comp(j))
            write(vtk_file_data_name,*) trim(comp(j)),'_',trim(prop(k)),'-kernel'
            if (normalize_kernels) then
               allocate(nkre(size(pre)),nkim(size(pre)))
               knorm = sqrt(maxval(pre**2+pim**2))
               if (knorm < TINY(1.0)) then
                  call add(errmsg,2,'Norm of sens kernel is zero.',myname)
                  goto 1
               endif
               nkre = pre/knorm
               nkim = pim/knorm
               call writeComplexDataVtkFile(vtkinfo,vtk_file_base,1,nkre,nkim,errmsg,&
                  data_name = trim(vtk_file_data_name),file_index = ifreq(jf))
               deallocate(nkre,nkim)
            else
               call writeComplexDataVtkFile(vtkinfo,vtk_file_base,1,pre,pim,errmsg,&
                  data_name = trim(vtk_file_data_name),file_index = ifreq(jf))
            endif
            if (.level.errmsg == 2) goto 1
         end do
      end do
      call h5fclose_f(fid,ierr)
      call dealloc(kernel)
      if (ierr < 0) then
         call add(errmsg,2,'Error closng HDF fle',myname)
         goto 1
      end if
   end do
!------------------------------------------------------------------------
!  clean up
!
   if (allocated(prop)) deallocate(prop)
   if (allocated(comp)) deallocate(comp)
   if (associated(ifreq)) deallocate(ifreq)
   if (associated(str_vec)) deallocate(str_vec)
   call dealloc(vtkinfo)
   call invgrid%dealloc()
   deallocate(invgrid)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(ap)
   call dealloc(mpisup)
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   end if
end program kernel2vtk
