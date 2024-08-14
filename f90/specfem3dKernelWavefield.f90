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
module specfem3dKernelWavefield
   use kernelWavefield
   use specfem3dForASKIFiles
   use errorMessage
   use fileUnitHandler
   use globalMpiInfo
   implicit none
!
   type, extends(kernel_wavefield) :: specfem3d_kernel_wavefield
      double complex, dimension(:,:), pointer:: u => null()           ! displacements (local wavefield points, 3)
      double complex, dimension(:,:), pointer :: ustr => null()       ! strains (local wavefield points, 3)
      double complex, dimension(:,:,:), pointer :: g => null()        ! Green tensor (local wavefield points, 3, forces)
      double complex, dimension(:,:,:), pointer :: gstr => null()     ! Green strains (local wavefield points, 6, forces)
      double precision, dimension(:), pointer :: petu => null()       ! phase end time for kernel displacements
      double precision, dimension(:), pointer :: petg => null()       ! phase end time for kernel Green tensors
   contains
      procedure :: readDisplacement => readFrequencyDisplacementSpecfem3dKernelWavefield
      procedure :: readGreenTensor => readFrequencyGreenTensorSpecfem3dKernelWavefield
      procedure :: getDisplacement => getDisplacementSpecfem3dKernelWavefield
      procedure :: getGreenTensor => getGreenTensorSpecfem3dKernelWavefield
      procedure :: getStrain => getStrainSpecfem3dKernelWavefield
      procedure :: getGreenStrain => getGreenStrainSpecfem3dKernelWavefield
      procedure :: dealloc => deallocSpecfem3dKernelWavefield
      procedure :: deallocDisplacement => deallocDisplacementSpecfem3dKernelWavefield
      procedure :: deallocGreen => deallocGreenSpecfem3dKernelWavefield
      procedure :: associateDisplacement => associateDisplacementSpecfem3dKernelWavefield
      procedure :: associateDisplacementPhaseEndTime => associateDisplacementPhaseEndTimeSpecfem3dKernelWavefield
      procedure :: readDisplacementPhaseEndTime => readDisplacementPhaseEndTimeSpecfem3dKernelWavefield
      procedure :: readGreenTensorPhaseEndTime => readGreenTensorPhaseEndTimeSpecfem3dKernelWavefield
      procedure :: getDisplacementPhaseEndTime => getDisplacementPhaseEndTimeSpecfem3dKernelWavefield
      procedure :: getGreenTensorPhaseEndTime => getGreenTensorPhaseEndTimeSpecfem3dKernelWavefield
   end type specfem3d_kernel_wavefield
!
contains
!-----------------------------------------------------------------------------------------------------
!  read specific file of displacement spectrum for the requested frequency (if existing)
!  jf frequency index for which kernel displacement should be read in (defines filename)
!
   subroutine readFrequencyDisplacementSpecfem3dKernelWavefield(this,nwp,jf,dfin,basename,evid,do_hyperslab,errmsg)
      class (specfem3d_kernel_wavefield) :: this
      integer :: nwp,jf
      double precision :: dfin
      character(len=*) :: basename,evid
      logical :: do_hyperslab
      type (error_message) :: errmsg
      ! local
      character(len=400) :: errstr
      character(len=49) :: myname = 'readFrequencyDisplacementSpecfem3dKernelWavefield'
      character(len=500) :: filename
      integer :: nproc,specfem_version,ntot,jf_file
      double precision :: df
      character(len=char_len_aski_output_id) :: id
   !
   ! read specific file for this frequency
   !
      write(filename,'(a,i6.6,a)') trim(basename)//trim(evid)//'_jf',jf,'.hdf'
      call readSpecfem3dForASKISpectralWavefieldFileHDF(filename,specfem_version,&
            id,nproc,ntot,df,jf_file,this%u,this%ustr,do_hyperslab,errmsg)
      if (.level.errmsg == 2) return

   ! check content of frequency file
      if (nproc /= numtasks .and. numtasks /= 1) then
         write(errstr,*) "nproc = ",nproc," contained in this file does not match numtasks = ",numtasks
         call add(errmsg,2,errstr,myname)
         return
      end if
      if (trim(id) /= trim(evid)) then
         call add(errmsg,2,"ID contained in this file '"//trim(id)//"' does not match ID '"//&
              trim(evid)//"' of main file",myname)
         return
      end if
      if (ntot /= nwp) then
         write(errstr,*) "ntot ",ntot," contained in this file does not match nwp ",nwp," of inversion grid"
         call add(errmsg,2,errstr,myname)
         return
      end if
      if (abs(df-dfin)/dfin > 1.e-3) then
         write(errstr,*) "df ",df," contained in this file does not match df ",dfin," of main parameter file"
         call add(errmsg,2,errstr,myname)
         return
      end if
      if (jf_file /= jf) then
         write(errstr,*) "requested frequency index ",jf," (also contained in filename) does not match the frequency index ",&
              jf_file,"actually contained in the file"
         call add(errmsg,2,errstr,myname)
         return
      end if
      this%nwp = nwp
   end subroutine readFrequencyDisplacementSpecfem3dKernelWavefield
!--------------------------------------------------------------------------------------
!  read specific files of kernel green tensor spectrum for the requested frequency (if existent)
!  jf frequency index for which kernel green tensor should be read in (defines filename)
!
   subroutine readFrequencyGreenTensorSpecfem3dKernelWavefield(this,nwp,jf,dfin,basename,netstaname,comp,do_hyperslab,errmsg)
      class (specfem3d_kernel_wavefield) :: this
      integer :: nwp,jf
      double precision :: dfin
      character(len=*) :: basename,netstaname
      character(len=char_len_comp), dimension(:) :: comp
      logical :: do_hyperslab
      type (error_message) :: errmsg
      ! local
      character(len=max_length_string) :: errstr
      character(len=48) :: myname = 'readFrequencyGreenTensorSpecfem3dKernelWavefield'
      character(len=max_length_string) :: filename
      integer :: icomp,nproc,specfem_version,ntot,jf_file,ncomp
      double precision :: df
      character(len=char_len_aski_output_id) :: staname_comp
      double complex, dimension(:,:), pointer :: g_comp,gstr_comp
   !
   ! read specific files for this frequency for all components
   !
      ncomp = size(comp)
      this%ncomp = ncomp
      allocate(this%g(nwp,3,ncomp))
      allocate(this%gstr(nwp,6,ncomp))
      do icomp = 1,ncomp
         write(filename,"(a,i6.6,a)") trim(basename)//trim(netstaname)//"_"//trim(comp(icomp))//"_jf",jf,".hdf"
         call readSpecfem3dForASKISpectralWavefieldFileHDF(filename,specfem_version,&
              staname_comp,nproc,ntot,df,jf_file,g_comp,gstr_comp,do_hyperslab,errmsg)
      !
      ! check content of file
      !
         if (trim(staname_comp) /= trim(netstaname)//'_'//trim(comp(icomp))) then
            write(errstr,*) "ASKI output ID read from file = '",trim(staname_comp),&
                 "' has not the expected value 'staname_comp' = '",trim(netstaname)//'_'//trim(comp(icomp)),&
                 "' referring to station name '",trim(netstaname),"' and component '",trim(comp(icomp)),&
                 "' of this Green function. By convention, ASKI output ID must be of this form."
            call add(errmsg,2,errstr,myname)
            return
         end if
         if (nproc /= numtasks .and. numtasks /= 1) then
            write(errstr,*) "nproc = ",nproc," contained in this file does not match numtasks = ",numtasks
            call add(errmsg,2,errstr,myname)
            return
         end if
         if (ntot /= nwp) then
            write(errstr,*) "ntot ",ntot," contained in this file does not match nwp ",nwp," of inversion grid"
            call add(errmsg,2,errstr,myname)
            return
         end if
         if (abs(df-dfin)/dfin > 1.e-3) then
            write(errstr,*) "df ",df," contained in this file does not match df ",dfin," of main parameter file"
            call add(errmsg,2,errstr,myname)
            return
         end if
         if (jf_file /= jf) then
            write(errstr,*) "requested frequency index ",jf," (also contained in filename) does not match the frequency index ",&
                 jf_file," which is actually contained in the file"
            call add(errmsg,2,errstr,myname)
            return
         end if
      !
         this%g(:,:,icomp) = g_comp
         this%gstr(:,:,icomp) = gstr_comp
         deallocate(g_comp,gstr_comp)
      end do ! icomp
      this%nwp = nwp
   end subroutine readFrequencyGreenTensorSpecfem3dKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate object
!
   subroutine deallocSpecfem3dKernelWavefield(this)
      class (specfem3d_kernel_wavefield) :: this
      if (associated(this%u)) deallocate(this%u)
      if (associated(this%ustr)) deallocate(this%ustr)
      if (associated(this%g)) deallocate(this%g)
      if (associated(this%gstr)) deallocate(this%gstr)
      if (associated(this%petu)) deallocate(this%petu)
      if (associated(this%petg)) deallocate(this%petg)
   end subroutine deallocSpecfem3dKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate displacements of object
!
   subroutine deallocDisplacementSpecfem3dKernelWavefield(this)
      class (specfem3d_kernel_wavefield) :: this
      if (associated(this%u)) deallocate(this%u)
      if (associated(this%ustr)) deallocate(this%ustr)
      if (associated(this%petu)) deallocate(this%petu)
   end subroutine deallocDisplacementSpecfem3dKernelWavefield
!-------------------------------------------------------------------------------
!  Deallocate Green tensor of object
!
   subroutine deallocGreenSpecfem3dKernelWavefield(this)
      class (specfem3d_kernel_wavefield) :: this
      if (associated(this%g)) deallocate(this%g)
      if (associated(this%gstr)) deallocate(this%gstr)
      if (associated(this%petg)) deallocate(this%petg)
   end subroutine deallocGreenSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Associate provided displacements with those of this
!
   subroutine associateDisplacementSpecfem3dKernelWavefield(this,u,ustr)
      class (specfem3d_kernel_wavefield) :: this
      double complex, dimension(:,:), target :: u,ustr
      this%u => u
      this%ustr => ustr
   end subroutine associateDisplacementSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Associate displacement phase end time of this with those of that
!
   subroutine associateDisplacementPhaseEndTimeSpecfem3dKernelWavefield(this,petu)
      class (specfem3d_kernel_wavefield) :: this
      double precision, dimension(:), target :: petu
      this%petu => petu
   end subroutine associateDisplacementPhaseEndTimeSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get strains
!
   function getStrainSpecfem3dKernelWavefield(this) result(ustr)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:), pointer :: ustr
      ustr => this%ustr
   end function getStrainSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get displacements
!
   function getDisplacementSpecfem3dKernelWavefield(this) result(u)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:), pointer :: u
      u => this%u
   end function getDisplacementSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get Green strains
!
   function getGreenStrainSpecfem3dKernelWavefield(this) result(gstr)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:,:), pointer :: gstr
      gstr => this%gstr
   end function getGreenStrainSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get Green tensor
!
   function getGreenTensorSpecfem3dKernelWavefield(this) result(g)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double complex, dimension(:,:,:), pointer :: g
      g => this%g
   end function getGreenTensorSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get displacement phase end time
!
   function getDisplacementPhaseEndTimeSpecfem3dKernelWavefield(this) result(petu)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double precision, dimension(:), pointer :: petu
      petu => this%petu
   end function getDisplacementPhaseEndTimeSpecfem3dKernelWavefield
!-------------------------------------------------------------------------
!  Get Green tensor phase end time
!
   function getGreenTensorPhaseEndTimeSpecfem3dKernelWavefield(this) result(petg)
      class (specfem3d_kernel_wavefield), intent(in) :: this
      double precision, dimension(:), pointer :: petg
      petg => this%petg
   end function getGreenTensorPhaseEndTimeSpecfem3dKernelWavefield
!------------------------------------------------------------------------
!  Read phase end time for kernel wavefield
!
   subroutine readDisplacementPhaseEndTimeSpecfem3dKernelWavefield(this,nwp,basename,evid,do_hyperslab,errmsg)
      class (specfem3d_kernel_wavefield) :: this
      integer :: nwp
      character(len=*) :: basename,evid
      logical :: do_hyperslab
      type (error_message) :: errmsg
      ! local
      character(len=400) :: errstr
      character(len=52) :: myname = 'readDisplacementPhaseEndTimeSpecfem3dKernelWavefield'
      character(len=500) :: filename
      integer :: nproc,ntot
   !
   ! read specific file
   !
      filename = trim(basename)//trim(evid)//'_pet.hdf'
      call readSpecfem3dForASKIPhaseEndTimeFileHDF(filename,nproc,ntot,this%petu,do_hyperslab,errmsg)
      if (.level.errmsg == 2) return
   !
      if (nproc /= numtasks .and. numtasks /= 1) then
         write(errstr,*) "nproc = ",nproc," contained in this file does not match numtasks = ",numtasks
         call add(errmsg,2,errstr,myname)
         return
      end if
      if (ntot /= nwp) then
         write(errstr,*) "ntot ",ntot," contained in this file does not match nwp ",nwp," of inversion grid"
         call add(errmsg,2,errstr,myname)
         return
      end if
   end subroutine readDisplacementPhaseEndTimeSpecfem3dKernelWavefield
!------------------------------------------------------------------------
!  Read phase end time for kernel Green tensor
!
   subroutine readGreenTensorPhaseEndTimeSpecfem3dKernelWavefield(this,nwp,basename,netstaname,do_hyperslab,errmsg)
      class (specfem3d_kernel_wavefield) :: this
      integer :: nwp
      character(len=*) :: basename,netstaname
      logical :: do_hyperslab
      type (error_message) :: errmsg
      ! local
      character(len=400) :: errstr
      character(len=51) :: myname = 'readGreenTensorPhaseEndTimeSpecfem3dKernelWavefield'
      character(len=500) :: filename
      integer :: nproc,ntot
   !
   ! read specific file
   !
      filename = trim(basename)//trim(netstaname)//'_pet.hdf'
      call readSpecfem3dForASKIPhaseEndTimeFileHDF(filename,nproc,ntot,this%petg,do_hyperslab,errmsg)
      if (.level.errmsg == 2) return
   !
      if (nproc /= numtasks .and. numtasks /= 1) then
         write(errstr,*) "nproc = ",nproc," contained in this file does not match numtasks = ",numtasks
         call add(errmsg,2,errstr,myname)
         return
      end if
      if (ntot /= nwp) then
         write(errstr,*) "ntot ",ntot," contained in this file does not match nwp ",nwp," of inversion grid"
         call add(errmsg,2,errstr,myname)
         return
      end if
   end subroutine readGreenTensorPhaseEndTimeSpecfem3dKernelWavefield
!
end module specfem3dKernelWavefield
