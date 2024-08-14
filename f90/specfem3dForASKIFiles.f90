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

module specfem3dForASKIFiles
      use hdf5
   use constants
   use errorMessage
   use hdfWrapper
   use globalMpiInfo
   use globalHdfInfo
   implicit none

   contains
!------------------------------------------------------------------------
!  read in specfem3d for ASKI main file
!  filename file name
!  errmsg error message
!  specfem_version specfem version
!  nproc number of procs
!  type_inversion_grid type of inversion grid
!  nwp number of wavefield points
!  df frequency step
!  nf number of frequencies
!  jf frequency indices
!  x,y,z coordinates of wavefield points
!  rho,vp,vs model parameters of reference model
!  ngllx,nglly,ngllz
!  jacobian jacobians on wavefield points
!  ncell number of invgrid cells
!  nb_idx invgrid cell neighbours (7,ncell)
!------------------------------------------------------------------------
!
   subroutine readSpecfem3dForASKIMainFileHDF(filename,&
        specfem_version,nproc,type_inversion_grid,nwp_all,nwp,offwp,df,nf,jf,&
        ngllx,nglly,ngllz,ncell,nb_idx,do_hyperslab,errmsg,x,y,z,rho,vp,vs,jacobian)
      ! incoming
      character(len=*) :: filename
      type (error_message) :: errmsg
      integer :: specfem_version,nproc,type_inversion_grid,nwp_all,nwp,offwp,nf,ngllx,nglly,ngllz,ncell
      double precision :: df
      integer, dimension(:), allocatable  :: jf
      double precision, dimension(:), allocatable, optional :: x,y,z,rho,vp,vs,jacobian
      integer, dimension(:,:), allocatable :: nb_idx
      logical :: do_hyperslab
      ! local
      integer, dimension(:,:), pointer :: neighbour
      real, dimension(:,:), pointer :: xyz,model
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer(kind=8) :: fid
      integer :: ierr,nwp_local,offwp_local
      integer(kind=8), dimension(:), allocatable :: dimsslab1d,dimsslab2d,offset1d,offset2d,count1d,count2d
      integer, dimension(:), pointer :: id
      real, dimension(:), pointer :: d
      character(len=31) :: myname = "readSpecfem3dForASKIMainFileHDF"
   !
      if (myrank == 0) then
         print *,"Open ASKI main file: ",trim(filename)
      endif
      if (numtasks > 1) then
         call openFileParallelAccessHDFWrapper(filename,fid,errmsg)
      else
         call openFileRoHDFWrapper(filename,fid,errmsg)
      endif
      if (.level.errmsg == 2) goto 1
   !
      call readArrayAttributeHDFWrapper(fid,"mixed_integers",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d()
      specfem_version = id(1)
      nproc = id(2)
      type_inversion_grid = id(3)
      nf = id(4)
      deallocate(id)
   !
   !  make sure that number of processes is equal to nproc if hyperslab reading is desired
   !
      if (do_hyperslab) then
         if (nproc /= numtasks) then
            call add(errmsg,2,"Number of processes > 1 and /= to number of processes used to compute wavefield spectra",myname)
            return
         end if
      end if
   !
      call readArrayAttributeHDFWrapper(fid,"aski_jf",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d()
      allocate(jf(size(id)))
      jf = id
      deallocate(id)
   !
      call readArrayAttributeHDFWrapper(fid,"aski_df",arra,errmsg)
      if (.level.errmsg == 2) goto 1
      d => arra%get1d()
      df = d(1)
      deallocate(d)
   !
      call readArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d()
      if (do_hyperslab) then
         nwp_local = id(myrank+1)
         offwp_local = sum(id(1:myrank))
      else
         nwp_local = sum(id)
         offwp_local = 0
      end if
      nwp_all = sum(id)
      nwp = nwp_local
      offwp = offwp_local
      dimsslab2d = [nwp_local,3]; offset2d = [offwp_local,0]; count2d = [nwp_local,3]
      dimsslab1d = [nwp_local]; offset1d = [offwp_local]; count1d = [nwp_local]
      deallocate(id)
!
      if (myrank == 0) then
         print *,"Total number of wavefield points = ",nwp_all
      endif
   !
      if (present(x) .or. present(y) .or. present(z)) then
         call readArrayHDFWrapper(fid,'wavefield_points',arra,errmsg,&
                 xferprp = hdf_xferprp,dimsslab = dimsslab2d,offset = offset2d,count = count2d)
         if (.level.errmsg == 2) goto 1
         xyz => arra%get2d()
         if (present(x)) allocate(x(nwp_local)); x = xyz(:,1)
         if (present(y)) allocate(y(nwp_local)); y = xyz(:,2)
         if (present(z)) allocate(z(nwp_local)); z = xyz(:,3)
         deallocate(xyz)
      end if
      if (myrank == 0) then
         print *," Successfully read coordinates of wavefield points"
      endif
   !
      if (present(rho) .or. present(vp) .or. present(vs)) then
         call readArrayHDFWrapper(fid,'model_values',arra,errmsg,&
                 xferprp = hdf_xferprp,dimsslab = dimsslab2d,offset = offset2d,count = count2d)
         if (.level.errmsg == 2) goto 1
         model => arra%get2d()
         if (present(rho)) allocate(rho(nwp_local)); rho = model(:,1)
         if (present(vp)) allocate(vp(nwp_local));  vp = model(:,2)
         if (present(vs)) allocate(vs(nwp_local));  vs = model(:,3)
         deallocate(model)
      end if
      if (myrank == 0) then
         print *," Successfully read model values of wavefield points"
      endif
   !
      if (type_inversion_grid == 4) then
         call readArrayAttributeHDFWrapper(fid,"cell_info",aria,errmsg)
         if (.level.errmsg == 2) goto 1
         id => aria%get1d()
         ngllx = id(1)
         nglly = id(2)
         ngllz = id(3)
         ncell = id(4)
         deallocate(id)
         if (myrank == 0) then
            print *," NGLLX, NGLLY, NGLLZ, NCELL = ",ngllx,nglly,ngllz,ncell
         endif
      !
         if (present(jacobian)) then
            call readArrayHDFWrapper(fid,'jacobian',arra,errmsg,&
                  xferprp = hdf_xferprp,dimsslab = dimsslab1d,offset = offset1d,count = count1d)
            if (.level.errmsg == 2) goto 1
            allocate(jacobian(nwp_local))
            d => arra%get1d()
            jacobian = d
            deallocate(d)
            if (myrank == 0) then
               print *,'Jacobian(1): ',jacobian(1)
               print *," Successfully read jacobian at wavefield points"
            endif
         end if
      !
      !  no distribution of neighbours across processes
      !  number of cells is 1/125 of number of wavefeld points
      !  so each process can easily hold info on neighbours of all cells
      !
         call readArrayHDFWrapper(fid,'neighbours',aria,errmsg)
         if (.level.errmsg == 2) goto 1
         neighbour => aria%get2d()
         allocate(nb_idx(7,ncell))
         nb_idx = neighbour
         deallocate(neighbour)
         if (myrank == 0) then
            print *,'nb_idx(:,1) = ',nb_idx(:,1)
            print *," Successfully read indices of neighbours"
         endif
      endif
1     call h5fclose_f(fid,ierr)
   end subroutine readSpecfem3dForASKIMainFileHDF
!  ------------------------------------------------------------------------
!  read in wavefield spectra computed by SPECFEM
!  filename file name
!  errmsg error message
!  specfem_version specfem version
!  file_ID file ID
!  nproc number of procs
!  nwp number of wavefield points
!  df frequency step
!  jf frequency index for which spectra are contained in the file
!  u displacement spectrum
!  ustr strain spectrum
!
   subroutine readSpecfem3dForASKISpectralWavefieldFileHDF(filename,specfem_version,file_ID,&
                                                           nproc,nwp_local,df,jfcur,u,ustr,do_hyperslab,errmsg)
      character (len=*) :: filename
      integer :: specfem_version,nproc,nwp_local,jfcur
      character(len=char_len_aski_output_id) :: file_ID       ! should be same as len_ASKI_output_id
      double precision :: df
      double complex, dimension(:,:), pointer :: u,ustr
      logical ::  do_hyperslab
      type (error_message) :: errmsg
      integer, dimension(:), pointer :: id
      real, dimension(:), pointer:: d
      real, dimension(:,:), pointer :: ure,uim
      type (any_rank_real_array) :: arra,arra2
      type (any_rank_integer_array) :: aria
      integer(kind=8) :: fid
      integer(kind=8), dimension(:), allocatable :: dimsslab2d,offset2d,count2d
      integer :: slen,ierr,offwp
      character(len=400) :: cval
      character(len=44) :: myname = 'readSpecfem3dForASKISpectralWavefieldFileHDF'

      if (numtasks > 1) then
         call openFileParallelAccessHDFWrapper(filename,fid,errmsg)
      else
         call openFileRoHDFWrapper(filename,fid,errmsg)
      endif
      if (.level.errmsg == 2) goto 1
   !
      call readStringAttributeHDFWrapper(fid,'aski_output_id',cval,slen,errmsg)
      if (.level.errmsg == 2) goto 1
      if (slen > char_len_aski_output_id) then
         call add(errmsg,2,'returned string length for ASKI_output_ID greater than assumed length',myname)
         goto 1
      end if
      file_ID = cval(1:slen)
   !
      call readArrayAttributeHDFWrapper(fid,"nproc_nwp_jf_version",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d(); nproc = id(1); jfcur = id(2); specfem_version = id(3)
      deallocate(id)
      !
      if (do_hyperslab) then
         if (nproc /= numtasks) then
            call add(errmsg,2,"Number of processes > 1 and /= to number of processes used to compute wavefield spectra",myname)
            return
         end if
      end if
   !
      call readArrayAttributeHDFWrapper(fid,"aski_df",arra,errmsg)
      if (.level.errmsg == 2) goto 1
      d => arra%get1d(); df = d(1)
      deallocate(d)
   !
      call readArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d()
      if (do_hyperslab) then
         nwp_local = id(myrank+1)
         offwp = sum(id(1:myrank))
      else
         nwp_local = sum(id)
         offwp = 0
      end if
      dimsslab2d = [nwp_local,9]; offset2d = [offwp,0]; count2d = [nwp_local,9]
      deallocate(id)
   !
      call readArrayHDFWrapper(fid,'spectra_real',arra,errmsg,&
            xferprp = hdf_xferprp,dimsslab = dimsslab2d,offset = offset2d,count = count2d)
      if (.level.errmsg == 2) goto 1
   !   if (myrank == 0) then
   !      print *," Successfully read real part of wavefield spectra at wavefield points"
   !   endif
      call readArrayHDFWrapper(fid,'spectra_imag',arra2,errmsg,&
            xferprp = hdf_xferprp,dimsslab = dimsslab2d,offset = offset2d,count = count2d)
      if (.level.errmsg == 2) goto 1
   !   if (myrank == 0) then
   !      print *," Successfully read imaginary part of wavefield spectra at wavefield points"
   !   endif
      ure => arra%get2d()
      uim => arra2%get2d()
      allocate(u(nwp_local,3))
      allocate(ustr(nwp_local,6))
      u = dcmplx(ure(:,1:3),uim(:,1:3))
      ustr = dcmplx(ure(:,4:9),uim(:,4:9))
      deallocate(ure,uim)
 1    call h5fclose_f(fid,ierr)
   end subroutine readSpecfem3dForASKISpectralWavefieldFileHDF
!  ------------------------------------------------------------------------
!  read in phase end times output by SPECFEM
!  filename file name
!  file_ID file ID
!  nproc number of procs
!  nwp number of wavefield points
!  pet phase end time
!  errmsg error message
!
   subroutine readSpecfem3dForASKIPhaseEndTimeFileHDF(filename,nproc,nwp_local,pet,do_hyperslab,errmsg)
      character (len=*) :: filename
      integer :: nproc,nwp_local
      double precision, dimension(:), pointer :: pet
      logical ::  do_hyperslab
      type (error_message) :: errmsg
      integer, dimension(:), pointer :: id
      real, dimension(:), pointer:: d
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer(kind=8) :: fid
      integer(kind=8), dimension(1) :: dimsslab1d,offset1d,count1d
      integer :: ierr,offwp
      character(len=39) :: myname = 'readSpecfem3dForASKIPhaseEndTimeFileHDF'

      if (numtasks > 1) then
         call openFileParallelAccessHDFWrapper(filename,fid,errmsg)
      else
         call openFileRoHDFWrapper(filename,fid,errmsg)
      endif
      if (.level.errmsg == 2) goto 1
   !
      call readArrayAttributeHDFWrapper(fid,"nproc",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d(); nproc = id(1)
      deallocate(id)
      !
      if (do_hyperslab) then
         if (nproc /= numtasks) then
            call add(errmsg,2,"Number of processes > 1 and /= to number of processes used to compute wavefield spectra",myname)
            return
         end if
      end if
   !
      call readArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      id => aria%get1d()
      if (do_hyperslab) then
         nwp_local = id(myrank+1)
         offwp = sum(id(1:myrank))
      else
         nwp_local = sum(id)
         offwp = 0
      end if
      dimsslab1d = [nwp_local]; offset1d = [offwp]; count1d = [nwp_local]
      deallocate(id)
   !
      call readArrayHDFWrapper(fid,'phase_end_time',arra,errmsg,&
            xferprp = hdf_xferprp,dimsslab = dimsslab1d,offset = offset1d,count = count1d)
      if (.level.errmsg == 2) goto 1
      d => arra%get1d()
      allocate(pet(nwp_local))
      pet = d
      deallocate(d)
 1    call h5fclose_f(fid,ierr)
   end subroutine readSpecfem3dForASKIPhaseEndTimeFileHDF
!
end module specfem3dForASKIFiles
