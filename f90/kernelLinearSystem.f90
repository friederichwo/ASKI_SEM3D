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
!   Set up and handle kernel linear system
!
!   Given an object of type data_model_space_info, the kernel matrix as well as
!   the right-hand-side of the linear system is set up, not including the
!   regularization part.
!   This module assumes that the kernel linear system is split
!   by distributing the paths among available processes. Splitting is done by
!   working on a subset of the paths of the full data set.
!
!   Authors: Florian Schumacher, Wolfgang Friederich
!   Date July 2015, 2022
!
module kernelLinearSystem
!
   use mpi
   use dataModelSpaceInfo
   use spectralWaveformKernel
   use kernelInvertedModel
   use hdfWrapper
   use globalHdfInfo
   use globalMpiInfo
   use errorMessage
!
   implicit none
!
   interface dealloc; module procedure deallocateKernelLinearSystem; end interface
   interface operator (.sol.); module procedure getSolutionKernelLinearSystem; end interface
   interface operator (.ndata.); module procedure getNdataKernelLinearSystem; end interface
   interface operator (.nmval.); module procedure getNmvalKernelLinearSystem; end interface
   interface operator (.md.); module procedure getMeasuredDataKernelLinearSystem; end interface
   interface operator (.sd.); module procedure getSyntheticDataKernelLinearSystem; end interface
!
   type kernel_linear_system
      double precision, dimension(:,:), allocatable :: klsm  ! kernel system matrix, only data part
      integer :: ndata                                       ! (local) number of rows associated with data
      integer :: nmval                                       ! number of model values on inversion grid
      double precision, dimension(:), allocatable :: rhs         ! (local) right-hand-side vector of kernel system
      double precision, dimension(:), allocatable :: q           ! (local) vector q in CGLS-algorithm of size ndata
      double precision, dimension(:), allocatable :: sol         ! solution vector (has size nmval)
      double precision, dimension(:), allocatable :: g           ! gradient vector (has size nmval)
      double precision, dimension(:), allocatable :: p           ! vector p in CGLS-algorithm (has size nmval)
      double precision, dimension(:), allocatable :: mdata       ! (local) vector of measured data according to data space
      double precision, dimension(:), allocatable :: sdata       ! (local) vector of synthetic data according to data space
      double precision, dimension(:), allocatable :: sigma       ! (local) vector of data uncertainties derived from noise rms
      double precision, dimension(:), allocatable :: dnorm       ! (local) amplitude of data for normalisation of kls (set on data read)
      logical :: normalize                                       ! not used
      logical :: do_phase_inversion                              ! only invert the phase/time shift of the data
   end type kernel_linear_system
!
contains
!------------------------------------------------------------------------
!  Allocate system matrix, according to number of data and model values
!  dpi: flag for doing a phase inversion (default = .false.)
!
   subroutine initiateKernelLinearSystem(this,dmspace,dpi)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      integer :: nd
      logical :: dpi
   !
      this%ndata = .ndata.dmspace
      this%nmval = .nmval.dmspace
      nd = this%ndata
      if (dpi) then
         nd = this%ndata/2
      end if
   !
      allocate(this%klsm(nd,this%nmval))
      allocate(this%rhs(nd))
      allocate(this%q(nd))
      allocate(this%sol(this%nmval))
      allocate(this%p(this%nmval))
      allocate(this%g(this%nmval))
      allocate(this%mdata(this%ndata))
      allocate(this%sdata(this%ndata))
      allocate(this%sigma(this%ndata))
      allocate(this%dnorm(this%ndata))
   !
      this%do_phase_inversion = dpi
   end subroutine initiateKernelLinearSystem
! -----------------------------------------------------------------------------------
!   Initiate data part of kernel linear system only
!
   subroutine initiateDataKernelLinearSystem(this,dmspace)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
   !
      this%ndata = .ndata.dmspace
      allocate(this%mdata(this%ndata))
      allocate(this%sdata(this%ndata))
      allocate(this%sigma(this%ndata))
      allocate(this%dnorm(this%ndata))
   end subroutine initiateDataKernelLinearSystem
! -----------------------------------------------------------------------------------
!  Set up local part of kernel linear system related to data for one selected process
!  ifreq_iter: indices of frequencies for this iteration
!
   subroutine fillMatrixKernelLinearSystem(this,dmspace,path_sensitivity_kernels,propsetname,ifreq_iter,errmsg)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      character(len=*) :: path_sensitivity_kernels,propsetname
      integer, dimension(:) :: ifreq_iter
      type (error_message) :: errmsg
      type (spectral_waveform_kernel) :: skernel
      integer (kind=8) :: kfid
      integer :: ncell,nk,nimre,icomp,iprop,kf,off,id,jf,nfreq,idp
      integer, dimension(:), allocatable :: idx_data_path
      integer, dimension(:), allocatable :: ifreq
      real, dimension(:,:), pointer :: pre,pim
      double precision :: w1,w2,wp,resd,imsd,sdmag
      character(len=char_len_comp), dimension(:), allocatable :: comp
      character(len=char_len_evid) :: evid
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
      character(len=max_length_string) :: kernel_filebase,kernel_file
   !
      kernel_filebase = trim(path_sensitivity_kernels)//'spectral_kernel_'//trim(propsetname)//'_'
   !
   !  Loop over iteration specific frequency indices
   !  Open spectral waveform kernel file for this frequency
   !  The ones returned from data_model_space-info sould be a subset of them
   !
      nimre = 2
      nfreq = size(ifreq_iter)
      do jf = 1,nfreq
         write(kernel_file,'(a,i6.6,a)') trim(kernel_filebase)//'jf',ifreq_iter(jf),'.hdf'
         if (myrank == 0) print *,'reading sensitivity kernel for frequency ',ifreq_iter(jf),' from :',trim(kernel_file)
         call openFileParallelAccessHDFWrapper(kernel_file,kfid,errmsg)
         if (.level.errmsg == 2) return
         call readMetaSpectralWaveformKernel(skernel,kfid,ifreq_iter(jf),errmsg)
         if (.level.errmsg == 2) return
      !
      ! loop over my paths, skip path if frequency is not desired
      ! else read kernels for this path from file and add to kernel sytem matrix
      !
         do while(nextPathDataModelSpaceInfo(dmspace,evid,netstaname,idx_data_path,comp,ifreq))
            if (.not. any(ifreq == ifreq_iter(jf))) then
               deallocate(idx_data_path,comp,ifreq)
               cycle
            end if
            call readSpectralWaveformKernel(skernel,kfid,evid,netstaname,errmsg)
            if (.level.errmsg == 2) return
            ncell = .nk.skernel
            kf = findloc(ifreq,ifreq_iter(jf),1)           ! position of iter_ifreq(jf) in ifreq-array
         !
         !  add these kernels to kernel linear system
         !
            do icomp = 1,size(comp)
               off = (icomp-1)*size(ifreq)*nimre+(kf-1)*nimre   ! offset of real part of datum relative to idx_data_path(1)
               id = idx_data_path(1+off)
               call getValuesByCompSpectralWaveformKernel(skernel,comp(icomp),pre,pim)
               w1 = 1.d0/this%sigma(id)
               w2 = 1.d0/this%sigma(id+1)
               do iprop = 1,.nprop.skernel
                  if (this%do_phase_inversion) then
                     idp = (id+1)/nimre
                     resd = this%sdata(id)
                     imsd = this%sdata(id+1)
                     wp = w1/( abs(resd)+abs(imsd) )
                     this%klsm(idp,(iprop-1)*ncell+1:iprop*ncell) = wp*( pim(:,iprop)*resd-pre(:,iprop)*imsd )
                  else
                     this%klsm(id,(iprop-1)*ncell+1:iprop*ncell) = pre(:,iprop)*w1
                     this%klsm(id+1,(iprop-1)*ncell+1:iprop*ncell) = pim(:,iprop)*w2
                  end if
               end do
            end do
            call deallocateKernelSpectralWaveformKernel(skernel)
            if (allocated(idx_data_path)) deallocate(idx_data_path)
            if (allocated(comp)) deallocate(comp)
            if (allocated(ifreq)) deallocate(ifreq)
         end do
         call dealloc(skernel)
         call closeFileHDFWrapper(kfid,errmsg)
         if (.level.errmsg == 2) return
      end do
   end subroutine fillMatrixKernelLinearSystem
!------------------------------------------------------------------------
!  Read my share of measured data (Fourier transformed observed seismograms)
!
   subroutine readMeasuredDataKernelLinearSystem(this,dmspace,path_measured_data,ifreq_all,errmsg)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      character(len=*) :: path_measured_data
      integer, dimension(:) :: ifreq_all
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      integer (kind=8) :: fid
      integer, dimension(:), allocatable :: idx_data_path
      integer, dimension(:), allocatable :: ifreq
      integer :: nimre,icomp,jf,off,id,ic,maxcomp,kf
      real, dimension(:,:), pointer:: d
      character(len=max_length_string) :: specfile,pathname
      character(len=char_len_comp), dimension(:), allocatable :: comp
      character(len=char_len_comp), dimension(:), pointer :: comp_in_dmsp
      character(len=char_len_evid) :: evid,evid_prev
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
      double precision, external :: pythag
   !
      nimre = 2
      maxcomp = getNcompInDataModelSpaceInfo(dmspace)
      comp_in_dmsp => getCompsInDataModelSpaceInfo(dmspace)
   !
   !  loop over components occurring in dmspace outside of path loop
   !  to avoid opening and closing of file for each path
   !
      do ic = 1,maxcomp
         specfile = trim(path_measured_data)//'data_'//trim(comp_in_dmsp(ic))//'.hdf'
         if (numtasks > 1) then
            call openFileParallelAccessHDFWrapper(specfile,fid,errmsg)
            if (.level.errmsg == 2) return
         else
            call openFileRoHDFWrapper(specfile,fid,errmsg)
            if (.level.errmsg == 2) return
         end if
      !
      !  iterate through paths for this component
      !
         do while(nextPathDataModelSpaceInfo(dmspace,evid,netstaname,idx_data_path,comp,ifreq))
            icomp = findloc(comp,comp_in_dmsp(ic),1)
            if (icomp == 0) then
               deallocate(idx_data_path,comp,ifreq)
               cycle                                        ! path does not have this component, try next
            end if
            pathname = trim(evid)//"_"//trim(netstaname)
            call readArrayHDFWrapper(fid,trim(pathname),arra,errmsg,xferprp = hdf_xferprp)
            if (.level.errmsg == 2) return
            d => arra%get2d()
         !
         !  here we need to take into account that DFT of measured data is done for all frequencies
         !
            do jf = 1,size(ifreq)
               kf = findloc(ifreq_all,ifreq(jf),1)                       ! position of iter_ifreq(jf) in ifreq_all-array
               off = (icomp-1)*size(ifreq)*nimre+(jf-1)*nimre            ! offset of real part of datum relative to idx_data_path(1)
               id = idx_data_path(1+off)
            !
               this%dnorm(id) = pythag(dble(d(kf,2)),dble(d(kf,3)))           ! amplitude, returned
               this%dnorm(id+1) = this%dnorm(id)
               this%mdata(id) = d(kf,2)                                  ! real part
               this%mdata(id+1) = d(kf,3)                                ! imag part
            end do
            deallocate(d)
            if (allocated(idx_data_path)) deallocate(idx_data_path)
            if (allocated(comp)) deallocate(comp)
            if (allocated(ifreq)) deallocate(ifreq)
         end do
         call closeFileHDFWrapper(fid,errmsg)
         if (.level.errmsg == 2) return
      end do
   end subroutine readMeasuredDataKernelLinearSystem
!------------------------------------------------------------------------------
!  Read my share of synthetic data (Fourier transformed synthetic seismograms)
!  Path-wise construction. Multiply with weights from dmsp.
!
   subroutine readSyntheticDataKernelLinearSystem(this,dmspace,path_synthetic_data,ifreq_all,errmsg)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      character(len=*) :: path_synthetic_data
      integer, dimension(:) :: ifreq_all
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      integer (kind=8) :: fid
      integer, dimension(:), allocatable :: idx_data_path
      integer, dimension(:), allocatable :: ifreq
      integer :: nimre,icomp,jf,off,id,ic,maxcomp,kf
      real, dimension(:,:), pointer:: d
      character(len=max_length_string) :: specfile,pathname
      character(len=char_len_comp), dimension(:), allocatable :: comp
      character(len=char_len_comp), dimension(:), pointer :: comp_in_dmsp
      character(len=char_len_evid) :: evid
      character(len=char_len_sta+char_len_netcode+1) :: netstaname
   !
      nimre = 2
      maxcomp = getNcompInDataModelSpaceInfo(dmspace)
      comp_in_dmsp => getCompsInDataModelSpaceInfo(dmspace)
   !
   !  loop over components occurring in dmspace outside of path loop
   !  to avoid opening and closing of file for each path
   !  (what happens if other processes do not have this component?)
   !
      do ic = 1,maxcomp
         specfile = trim(path_synthetic_data)//'syn_'//trim(comp_in_dmsp(ic))//'.hdf'
         if (numtasks > 1) then
            call openFileParallelAccessHDFWrapper(specfile,fid,errmsg)
            if (.level.errmsg == 2) return
         else
            call openFileRoHDFWrapper(specfile,fid,errmsg)
            if (.level.errmsg == 2) return
         end if
      !
      !  iterate through paths for this component
      !
         do while(nextPathDataModelSpaceInfo(dmspace,evid,netstaname,idx_data_path,comp,ifreq))
            icomp = findloc(comp,comp_in_dmsp(ic),1)
            if (icomp == 0) cycle                               ! path does not have this component, try next
            pathname = trim(evid)//"_"//trim(netstaname)
            call readArrayHDFWrapper(fid,trim(pathname),arra,errmsg,xferprp = hdf_xferprp)
            if (.level.errmsg == 2) return
            d => arra%get2d()
         !
         !  here we take into account that DFT of synthetics data is done for all frequencies
         !
            do jf = 1,size(ifreq)
               kf = findloc(ifreq_all,ifreq(jf),1)                       ! position of iter_ifreq(jf) in ifreq_all-array
               off = (icomp-1)*size(ifreq)*nimre+(jf-1)*nimre            ! offset of real part of datum relative to idx_data_path(1)
               id = idx_data_path(1+off)
               this%sdata(id) = d(kf,2)                                  ! real part
               this%sdata(id+1) = d(kf,3)                                ! imag part
            end do
            deallocate(d)
            if (allocated(idx_data_path)) deallocate(idx_data_path)
            if (allocated(comp)) deallocate(comp)
            if (allocated(ifreq)) deallocate(ifreq)
         end do
         call closeFileHDFWrapper(fid,errmsg)
         if (.level.errmsg == 2) return
      end do
   end subroutine readSyntheticDataKernelLinearSystem
!------------------------------------------------------------------------------------
!  Assume that uncertainty of data is given by a basic calibration uncertainty
!  proprtional to |d| plus rms-noise: sigma = f*|d|+noise.
!  For the solver, residuals rhs and sensitivity kernels need to be divided by sigma.
!  Note: this routine must be called after readMeasuredData because dnorm is computed there
!
   subroutine computeSigmaKernelLinearSystem(this,dmspace,f)
      type (kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      double precision :: f
      double precision, dimension(:), pointer :: noise
   !
      call getNoiseDataSamplesDataModelSpaceInfo(dmspace,noise)
      this%sigma = f*this%dnorm+noise
      nullify(noise)
   end subroutine computeSigmaKernelLinearSystem
!------------------------------------------------------------------------------------
!  Write weighted data and current synthetics (= d-r) to HDF for all paths and frequencies
!  dmspace:  data model space info object
!  outpath: folder to which output is written
!  ext: extension to automatically generated file name
!  errmsg:  error message
!
   subroutine writeAllDataKernelLinearSystem(this,dmspace,ndata_local,outpath,ext,errmsg)
      type(kernel_linear_system) :: this
      type (data_model_space_info) :: dmspace
      integer :: ndata_local
      character(len=*) :: outpath,ext
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer, dimension(:), allocatable :: ndata_local_all,ndata_local_sb
      integer (kind=8) :: fid,dsp,dset,dsfreq,dspf
      integer (kind=8), dimension(2) :: dims2d,offset2d,count2d
      integer (kind=8), dimension(1) :: dims1d,offset1d,count1d
      integer :: ierr
      real, dimension(:,:), allocatable :: d
      double precision, dimension(:), pointer :: noise
      integer, dimension(:), pointer :: ifreq
      character (len=max_length_string) :: outfile
      character(len=30) :: myname = 'writeAllDataKernelLinearSystem'
   !
   !  collect ndata from other processes to compute data offsets
   !
      allocate(ndata_local_all(numtasks))
      allocate(ndata_local_sb(1))
      ndata_local_sb(1) = ndata_local
      call MPI_ALLGATHER(ndata_local_sb,1,MPI_INTEGER,ndata_local_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      if (ierr .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather ndata_local',myname)
         return
      end if
   !
   !  ifreq and wdata from data model space
   !
      call getNoiseDataSamplesDataModelSpaceInfo(dmspace,noise)
      ifreq => getIfreqArrayDataModelSpaceInfo(dmspace)
   !
   !  create HDF output file
   !
      outfile = trim(outpath)//'klsdata_'//trim(ext)//'.hdf'
      call createFileParallelAccessHDFWrapper(outfile,fid,errmsg)
      if (.level.errmsg == 2) return
   !
   !  setup data space and data set for the data of all processes
   !
      dims2d = [sum(ndata_local_all),4]
      call h5screate_simple_f(2,dims2d,dsp,ierr)
      if (ierr < 0) then
         call add(errmsg,2,'cannot create dataspace for kls data output',myname)
         return
      end if
      call h5dcreate_f(fid,'klsdata',H5T_NATIVE_REAL,dsp,dset,ierr)
      if (ierr < 0) then
         call add(errmsg,2,'cannot create dataset for kls datat output',myname)
         return
      endif
    !
      dims1d = [sum(ndata_local_all)]
      call h5screate_simple_f(1,dims1d,dspf,ierr)
      if (ierr < 0) then
         call add(errmsg,2,'cannot create dataspace for ifreq data output',myname)
         return
      end if
      call h5dcreate_f(fid,'ifreq',H5T_NATIVE_INTEGER,dspf,dsfreq,ierr)
      if (ierr < 0) then
         call add(errmsg,2,'cannot create dataset for ifreq data output',myname)
         return
      endif
   !
   !  write my share of data to HDF file
   !  write d-r which is the current value of the synthetic predictions
   !
      offset2d = [sum(ndata_local_all(1:myrank)),0]
      count2d = [ndata_local,4]
      allocate(d(ndata_local,4))
      d(:,1) = real(this%mdata,4)/this%sigma
      d(:,2) = real(d(:,1)-this%rhs,4)
      d(:,3) = real(this%dnorm,4)
      d(:,4) = real(noise,4)
      call arra%assoc2d(d)
      call writeArrayHDFWrapper(fid,'klsdata',arra,errmsg,xferprp = hdf_xferprp,&
                                ds = dset,offset = offset2d,count = count2d)
      if (.level.errmsg == 2) return
      call arra%deassoc()
   !
      offset1d = [sum(ndata_local_all(1:myrank))]
      count1d = [ndata_local]
      call aria%assoc1d(ifreq)
      call writeArrayHDFWrapper(fid,'ifreq',aria,errmsg,xferprp = hdf_xferprp,&
                                ds = dsfreq,offset = offset1d,count = count1d)
      if (.level.errmsg == 2) return
      call aria%deassoc()

      deallocate(d,ndata_local_all,ndata_local_sb)
      nullify(ifreq,noise)
      call h5dclose_f(dset,ierr)
      call h5sclose_f(dsp,ierr)
      call h5dclose_f(dsfreq,ierr)
      call h5sclose_f(dspf,ierr)
      if (ierr < 0) then
         call add(errmsg,2,'cannot create dataset for kls datat output',myname)
         return
      endif
      call closeFileHDFWrapper(fid,errmsg)
   end subroutine writeAllDataKernelLinearSystem
!------------------------------------------------------------------------------------
!  Set (local) right hand side as measured minus synthetic data
!
   subroutine setRhsAsDataResidualKernelLinearSystem(this)
      type (kernel_linear_system) :: this
      integer :: id,idp
      double precision :: resd,imsd,sdmag,phisyn,phidat,phidiff
      double precision :: w1,wp
      if (this%do_phase_inversion) then
         do id = 1,this%ndata,2
            phidat = atan2(this%mdata(id+1),this%mdata(id))
            resd = this%sdata(id)
            imsd = this%sdata(id+1)
            sdmag = hypot(resd,imsd)
            phisyn = atan2(imsd,resd)
            phidiff = phidat-phisyn
            if (phidiff > mc_pid)  phidiff = phidiff-2.*mc_pid
            if (phidiff < -mc_pid) phidiff = phidiff+2.*mc_pid
            idp = (id+1)/2
            w1 = 1.d0/this%sigma(id)
            wp = w1*sdmag**2/( abs(resd)+abs(imsd) )
            this%rhs(idp) = wp*phidiff
         end do
      else
         this%rhs = (this%mdata-this%sdata)/this%sigma
      end if
   end subroutine setRhsAsDataResidualKernelLinearSystem
!------------------------------------------------------------------------
!  deallocate kernel linear system
!
   subroutine deallocateKernelLinearSystem(this)
      type (kernel_linear_system) :: this
      if (allocated(this%klsm)) deallocate(this%klsm)
      if (allocated(this%rhs)) deallocate(this%rhs)
      if (allocated(this%q)) deallocate(this%q)
      if (allocated(this%sol)) deallocate(this%sol)
      if (allocated(this%g)) deallocate(this%g)
      if (allocated(this%p)) deallocate(this%p)
      if (allocated(this%mdata)) deallocate(this%mdata)
      if (allocated(this%sdata)) deallocate(this%sdata)
      if (allocated(this%sigma)) deallocate(this%sigma)
   end subroutine deallocateKernelLinearSystem
!--------------------------------------------------------------------------------
!  Set solution vector either to zero or to starting model
!
   subroutine setNonzeroSolutionKernelLinearSystem(this,kim)
      type (kernel_linear_system) :: this
      type (kernel_inverted_model) :: kim
      integer :: i,ic,ip
   !
      i = 1
      do ip = 1,kim%nprop
         do ic = 1,kim%ncell
            this%sol(i) = kim%model_values(ic,ip)
            i = i+1
         end do
      enddo
   end subroutine setNonzeroSolutionKernelLinearSystem
!--------------------------------------------------------------------------------
!  Set solution vector either to zero or to starting model
!
   subroutine setZeroSolutionKernelLinearSystem(this)
      type (kernel_linear_system) :: this
   !
      this%sol = 0.d0
   end subroutine setZeroSolutionKernelLinearSystem
!------------------------------------------------------------------------
!  compute in place rhs - A*sol (using my values of rhs only)
!
   subroutine rhsMinusMatrixDotSolKernelLinearSystem(this)
      type (kernel_linear_system) :: this
   !
   !  y <- alfa*A*x+beta*b
   !  call DGEMV('Trans',m,n,alfa,A,lda,x,incx,beta,y,incy)
   !
      if (this%do_phase_inversion) then
         call DGEMV('No-Transpose',this%ndata/2,this%nmval,-1.d0,this%klsm,this%ndata/2,this%sol,1,1.d0,this%rhs,1)
      else
         call DGEMV('No-Transpose',this%ndata,this%nmval,-1.d0,this%klsm,this%ndata,this%sol,1,1.d0,this%rhs,1)
      end if
   end subroutine rhsMinusMatrixDotSolKernelLinearSystem
!------------------------------------------------------------------------
!  compute norm of my share of rhs and then do allreduce
!
   subroutine computeNormSquaredRhsKernelLinearSystem(this,norm_sq,errmsg)
      type (kernel_linear_system) :: this
      double precision :: norm_sq,norm_sq_loc
      type (error_message) :: errmsg
      integer :: ios
      double precision :: DDOT
      external DDOT
   !
      if (this%do_phase_inversion) then
         norm_sq_loc = DDOT(this%ndata/2,this%rhs,1,this%rhs,1)
      else
         norm_sq_loc = DDOT(this%ndata,this%rhs,1,this%rhs,1)
      end if
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of rhs','computeNormSquaredRhsKernelLinearSystem')
         return
      end if
   end subroutine computeNormSquaredRhsKernelLinearSystem
!------------------------------------------------------------------------
!  compute g = A^T rhs using my local values of rhs
!
   subroutine computeGradIsMatrixTransRhsKernelLinearSystem(this,errmsg)
      type (kernel_linear_system) :: this
      double precision, dimension(:), allocatable :: g_loc
      type (error_message) :: errmsg
      integer :: ios
   !
   !  y <- alfa*A*x+beta*b
   !  call DGEMV('Trans',m,n,alfa,A,lda,x,incx,beta,y,incy)
   !
      allocate(g_loc(this%nmval))
      if (this%do_phase_inversion) then
         call DGEMV('Transpose',this%ndata/2,this%nmval,1.d0,this%klsm,this%ndata/2,this%rhs,1,0.d0,g_loc,1)
      else
         call DGEMV('Transpose',this%ndata,this%nmval,1.d0,this%klsm,this%ndata,this%rhs,1,0.d0,g_loc,1)
      end if
      call MPI_ALLREDUCE(g_loc,this%g,this%nmval,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather results of A^T*rhs','computeGradIsMatrixTransRhsKernelLinearSystem')
         return
      end if
      deallocate(g_loc)
   end subroutine computeGradIsMatrixTransRhsKernelLinearSystem
!------------------------------------------------------------------------
!  compute norm of gradient
!
   subroutine computeNormSquaredGradKernelLinearSystem(this,norm_sq)
      type (kernel_linear_system) :: this
      double precision :: norm_sq
      double precision :: DDOT
      external DDOT
   !
      norm_sq = DDOT(this%nmval,this%g,1,this%g,1)
   end subroutine computeNormSquaredGradKernelLinearSystem
!------------------------------------------------------------------------
!  compute norm of solution
!
   subroutine computeNormSquaredSolKernelLinearSystem(this,norm_sq)
      type (kernel_linear_system) :: this
      double precision :: norm_sq
      double precision :: DDOT
      external DDOT
   !
      norm_sq = DDOT(this%nmval,this%sol,1,this%sol,1)
   end subroutine computeNormSquaredSolKernelLinearSystem
!------------------------------------------------------------------------
!  compute column sums of my share of klsm
!  cs = sum_ndata [ abs(kre(i,:)) + abs(kim(i,:) ]
!  colsum is real because it is later written to HDF
!
   subroutine computeColumnSumsKernelLinearSystem(this,colsum,errmsg)
      type (kernel_linear_system) :: this
      real, dimension(:,:), pointer :: colsum
      type (error_message) :: errmsg
      real, dimension(:), allocatable :: colsum_loc
      integer :: ios,l
      double precision :: DASUM
      external DASUM
   !
      allocate(colsum_loc(this%nmval), colsum(this%nmval,1))
      if (this%do_phase_inversion) then
         do l = 1,this%nmval
            colsum_loc(l) = DASUM(this%ndata/2,this%klsm(:,l),1)
         end do
      else
         do l = 1,this%nmval
            colsum_loc(l) = DASUM(this%ndata,this%klsm(:,l),1)
         end do
      end if
   !
      call MPI_REDUCE(colsum_loc,colsum(:,1),this%nmval,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of rhs','computeColumnSumsKernelLinearSystem')
         return
      end if
      deallocate(colsum_loc)
   end subroutine computeColumnSumsKernelLinearSystem
!------------------------------------------------------------------------
!  compute or update p = g + (gamma/gamma_old)*p
!
   subroutine updatePKernelLinearSystem(this,gamma,gamma_old)
      type (kernel_linear_system) :: this
      double precision :: gamma,gamma_old
   !
   !  first scale p and then add to g
   !  in first iteration with gamma_old < 0, set p = g
   !
      if (gamma_old < 0) then
         call DCOPY(this%nmval,this%g,1,this%p,1)
      else
         call DSCAL(this%nmval,gamma/gamma_old,this%p,1)
         call DAXPY(this%nmval,1.d0,this%g,1,this%p,1)
      end if
   end subroutine updatePKernelLinearSystem
!------------------------------------------------------------------------
!  compute my share of q = A*p
!
   subroutine qIsMatrixDotPKernelLinearSystem(this)
      type (kernel_linear_system) :: this
      if (this%do_phase_inversion) then
         call DGEMV('No-Transpose',this%ndata/2,this%nmval,1.d0,this%klsm,this%ndata/2,this%p,1,0.d0,this%q,1)
      else
         call DGEMV('No-Transpose',this%ndata,this%nmval,1.d0,this%klsm,this%ndata,this%p,1,0.d0,this%q,1)
      end if
   end subroutine qIsMatrixDotPKernelLinearSystem
!------------------------------------------------------------------------
!  compute squared norm of my share of q and the do allreduce
!
   subroutine computeNormSquaredQKernelLinearSystem(this,norm_sq,errmsg)
      type (kernel_linear_system) :: this
      double precision :: norm_sq,norm_sq_loc
      type (error_message) :: errmsg
      integer :: ios
      double precision :: DDOT
      external DDOT
   !
      if (this%do_phase_inversion) then
         norm_sq_loc = DDOT(this%ndata/2,this%q,1,this%q,1)
      else
         norm_sq_loc = DDOT(this%ndata,this%q,1,this%q,1)
      end if
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of q','computeNormSquaredQKernelLinearSystem')
         return
      end if
   end subroutine computeNormSquaredQKernelLinearSystem
!------------------------------------------------------------------------
!  update rhs = rhs - alfa * q
!
   subroutine updateRhsKernelLinearSystem(this,alfa)
      type (kernel_linear_system) :: this
      double precision :: alfa
   !
      if (this%do_phase_inversion) then
         call DAXPY(this%ndata/2,-alfa,this%q,1,this%rhs,1)
      else
         call DAXPY(this%ndata,-alfa,this%q,1,this%rhs,1)
      end if
   end subroutine updateRhsKernelLinearSystem
!------------------------------------------------------------------------
!  update sol = sol + alfa * p
!
   subroutine updateSolKernelLinearSystem(this,alfa)
      type (kernel_linear_system) :: this
      double precision :: alfa
   !
      call DAXPY(this%nmval,alfa,this%p,1,this%sol,1)
   end subroutine updateSolKernelLinearSystem
!------------------------------------------------------------------------
!  return pointer to solution vector
!
   function getSolutionVectorKernelLinearSystem(this) result(res)
      type (kernel_linear_system), intent(in), target :: this
      double precision, dimension(:), pointer :: res
      res => this%sol
   end function getSolutionVectorKernelLinearSystem
!------------------------------------------------------------------------
!  return pointer to measured data vector
!
   function getMeasuredDataKernelLinearSystem(this) result(mdata)
      type (kernel_linear_system), intent(in), target :: this
      double precision, dimension(:), pointer :: mdata
      mdata => this%mdata
   end function getMeasuredDataKernelLinearSystem
!------------------------------------------------------------------------
!  return pointer to synthetic data vector
!
   function getSyntheticDataKernelLinearSystem(this) result(sdata)
      type (kernel_linear_system), intent(in), target :: this
      double precision, dimension(:), pointer :: sdata
      sdata => this%sdata
   end function getSyntheticDataKernelLinearSystem
!------------------------------------------------------------------------
!  return pointer to solution of kernel linear system (i.e. model uptdate)
!
   function getSolutionKernelLinearSystem(this) result(sol)
      type (kernel_linear_system), intent(in), target :: this
      double precision, dimension(:), pointer :: sol
      sol => this%sol
    end function getSolutionKernelLinearSystem
!------------------------------------------------------------------------
!  return number of model values nmval of this kernel linear system
!
   function getNmvalKernelLinearSystem(this) result(n)
      type (kernel_linear_system), intent(in) :: this
      integer :: n
      n = this%nmval
   end function getNmvalKernelLinearSystem
!------------------------------------------------------------------------
!  return number of data samples ndata (first rows) of this kernel linear system
!
   function getNdataKernelLinearSystem(this) result(n)
      type (kernel_linear_system), intent(in) :: this
      integer :: n
      if (this%do_phase_inversion) then
         n = this%ndata/2
      else
         n = this%ndata
      end if
   end function getNdataKernelLinearSystem
!
end module kernelLinearSystem
