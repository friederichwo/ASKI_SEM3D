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
program kdispl2vtk
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use specfem3dKernelWavefield
   use seismicEventList
   use seismicNetwork
   use vtkFile
   use argumentParser
   use string
   use fileUnitHandler
   use errorMessage
   use mpiSupport
   use globalMpiInfo

   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   class (kernel_wavefield), allocatable :: kd
   type (vtk_info) :: vtk
   type(mpi_support) ::mpisup
   class (inversion_grid), allocatable :: invgrid
   type (input_parameter), pointer :: inpar_inv,inpar_iter
   integer, dimension(:), pointer :: ifreq_iterbasics
   integer :: ifreq,nwp,ic
   double precision :: df,dnorm,vtk_scale
   double complex, dimension(:,:), pointer :: ustr,u
   real, dimension(:), allocatable :: datare,dataim
   logical :: normalize_displ,do_hyperslab
   character(len=max_length_string) :: ucomp
   character(len=max_length_string) :: kd_file,kd_filebase,forward_method
   character(len=max_length_string) :: vtk_file_base,vtk_file_data_name,vtk_path
   character(len=max_length_string) :: vtk_format,vtk_geometry_type
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: evid
   character (len=10) :: myname = 'kdispl2vtk'
!-----------------------------------------------------------------------------
!  initialise MPI and errmsg
!
   call new(mpisup)
   myrank = .myrank.mpisup
   numtasks = .numtasks.mpisup
   call new(errmsg,myname)
!------------------------------------------------------------------------
!  processing of command line
!
   call init(ap,myname,'Extract kernel displacement spectra to vtk files '//&
             'for certain wavefield and strain components and frequencies')
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-evid',.true.,"event id of the kernel displacement object. This option must be set.",'sval','')
   call addOption(ap,'-ifreq',.true.,"frequency index at which the wavefield output should be extracted. ",'ival','1')
   call addOption(ap,'-ucomp',.true.,"wavefield component which can be 'ux', 'uy', 'uz'"//&
         " and 'exx', 'eyy', 'ezz', 'eyz', 'exz', 'exy' (denoting the strain components). ",'sval','ux')
   call addOption(ap,"-norm",.false.,"if set, kernel displacement absolute value is normalized to 1.0")
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   call document(ap)
!
   main_parfile = ap.sval.'main_parfile'
   if (.not. (ap.optset.'-evid')) then
      call add(errmsg,2,'option -evid must be set',myname)
      call usage(ap)
      goto 1
   end if
   evid = ap.sval.'-evid'
   normalize_displ = ap.optset.'-norm'
   ucomp = ap.sval.'-ucomp'
   select case (ucomp)
   case('ux','uy','uz','exx','eyy','ezz','eyz','exz','exy')
   ! OK, do nothing
   case default
      call add(errmsg,2,'invalid component',myname)
      goto 1
   end select
!------------------------------------------------------------------------
!  setup basics
!
   do_hyperslab = .false.
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   ifreq_iterbasics => getFrequencyIndicesIterationStepBasics(iterbasics)
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter =>  getInputParameterIterationStepBasics(iterbasics)
   vtk_format = trim(inpar_inv.sval.'DEFAULT_VTK_FILE_FORMAT')
   vtk_scale = inpar_iter.dval.'VTK_COORDS_SCALING_FACTOR'
   vtk_geometry_type = inpar_iter.sval.'VTK_GEOMETRY_TYPE'
   forward_method = inpar_inv.sval.'FORWARD_METHOD'
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(.iterpath.invbasics,do_hyperslab,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
!
!  set the type of kernel displacement according to forward method
!
   if (equalString(forward_method,'specfem3d')) then
      allocate(specfem3d_kernel_wavefield :: kd)
   else
      call add(errmsg,2,'Forward method '//trim(forward_method)//' not implemented',myname)
      goto 1
   end if
!
   ifreq = ap.ival.'-ifreq'
   if (.not. any(ifreq == ifreq_iterbasics)) then
      call add(errmsg,2,'Invalid frequency index',myname)
      goto 1
   end if
!
   if (.not. searchEventidSeismicEventList(.evlist.invbasics,evid)) then
      call add(errmsg,2,"event ID '"//trim(evid)//"'is not contained in event list",myname)
      goto 1
   end if
!------------------------------------------------------------------------
!  prepare for below
!
   df = .df.invbasics
   nwp = invgrid%getNwpLocal()
   kd_filebase = trim(.iterpath.invbasics)//trim(inpar_inv.sval.'PATH_KERNEL_DISPLACEMENTS')
   write(kd_file,'(a,i6.6,a)') trim(kd_filebase)//trim(evid)//'_jf',ifreq,'.hdf'
   write(*,*) "KERNEL DISPLACEMENT FILE '",trim(kd_file),"' TO READ"
   allocate(datare(nwp),dataim(nwp))
   vtk_path = trim(.iterpath.iterbasics)//trim(inpar_inv.sval.'PATH_VTK_FILES')//'kernel_disp'
!
!  read kernel displ and strains
!
   call kd%readDisplacement(nwp,ifreq,df,kd_filebase,evid,do_hyperslab,errmsg)
   if (.level.errmsg == 2) goto 1
   u => kd%getDisplacement()
   if(size(u,1) /= nwp) then
      call add(errmsg,2,'inconsistent number of wavefiel points in kernel displacement',myname)
      goto 1
   end if
   ustr => kd%getStrain()
   ic = 1
   select case(ucomp)
      case ('ux'); ic = 1
      case ('uy'); ic = 2
      case ('uz'); ic = 3
      case ('exx'); ic = 4
      case ('eyy'); ic = 5
      case ('ezz'); ic = 6
      case ('eyz'); ic = 7
      case ('exz'); ic = 8
      case ('exy'); ic = 9
   end select
   if (ic <= 3) then
      datare = real(u(:,ic))
      dataim = real(aimag(u(:,ic)),4)
   else
      datare = real(ustr(:,ic-3))
      dataim = real(aimag(ustr(:,ic-3)),4)
   end if
!
!  finally write vtk file
!
   if (normalize_displ) then
      dnorm = sqrt(maxval(datare**2+dataim**2))
      datare = datare/dnorm
      dataim = dataim/dnorm
   endif
   write(vtk_file_base,"(a,2('_',a))") trim(vtk_path),trim(evid),trim(ucomp)
   call initiateVtkFile(vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,2,errmsg)
   if (.level.errmsg == 2) goto 1
   write(vtk_file_data_name,*) trim(ucomp),'-kdispl'
   call writeComplexDataVtkFile(vtk,vtk_file_base,1,datare,dataim,errmsg,&
                               data_name = trim(vtk_file_data_name),file_index = ifreq)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  clean up before terminating the program
!
   if(associated(ustr)) nullify(ustr)
   if(associated(u)) nullify(u)
   if(associated(ifreq_iterbasics)) nullify(ifreq_iterbasics)
   call kd%dealloc()
   call invgrid%dealloc()
   deallocate(invgrid)
   call dealloc(iterbasics); call dealloc(invbasics)
   call dealloc(errmsg)
   call dealloc(ap)
   call dealloc(vtk)
   if(allocated(datare)) deallocate(datare)
   if(allocated(dataim)) deallocate(dataim)
   deallocate(kd)
   call dealloc(mpisup)
!
!  error treatment
!
 1 if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   end if
end program kdispl2vtk
