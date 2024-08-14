!----------------------------------------------------------------------------
!   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!  Collect kernel matrix column sums for different frequency bands
!
!  Input: main_parfile, list of iterations,
!----------------------------------------------------------------------------
program collectColsums
!
   use regularSphericalGrid
   use hdfWrapper
   use inversionBasics
   use inputParameter
   use argumentParser
   use errorMessage
   use string
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (input_parameter), pointer :: inpar_inv
   type (regular_spherical_grid) :: rsg,rsgcols
   integer :: i,it,nit,itref,icol
   integer, dimension(:), allocatable :: itlist
   real :: csmax
   character(len=max_length_string), dimension(:), pointer :: res
   character(len=char_len_par) :: prop
   character(len=max_length_string) :: rsgfile
   character(len=max_length_string) :: main_parfile
   character(len=max_length_string) :: main_path,iter_path,rsgpath,outpath
   character(len=14) :: myname = 'collectColsums'
!
   call new(errmsg,myname)
!
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  command line processing
!
   call init(ap,myname,"Collect column sums from different iterations")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-rsgpath',.true.,'path to column sums under output files','sval','kimrsg/')
   call addOption(ap,'-itlist',.true.,'List of iterations','ivec','1 3 5')
   call addOption(ap,'-prop',.true.,'Material property','sval','vp')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   rsgpath = ap.sval.'-rsgpath'
   itlist = ap.ivec.'-itlist'
   nit = size(itlist)
   prop = ap.sval.'-prop'
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
!
!  Loop over iterations
!
   itref = 0
   do i = 1,nit
      it = itlist(i)
      write(iter_path,"(2a,i3.3,a)") trim(main_path),trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it,'/'
      outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
      write(rsgfile,'(a,a4,i2.2,a5,i2.2,a)') &
            trim(outpath)//trim(rsgpath)//'pseudo_mesh_'//trim(prop),'_it_',it,'_itr_',itref,'.hdf'
      print *,'Read column sums from ',trim(rsgfile)
!
!  read column sums for given iteration
!
      call readHDFRegularSphericalGrid(rsg,trim(rsgfile),errmsg)
      if (.level.errmsg == 2) goto 1
      res => getWordsString(trim(rsg%field_name_string),':',errmsg)
      icol = findloc(res,'colsum',1)
      if (.level.errmsg == 2) goto 1
      if (i == 1) then
         call rsgcols%create(1,rsg%lat_pole,rsg%lon_pole,rsg%nr,rsg%nlat,rsg%nlon,rsg%dr,rsg%dlat,rsg%dlon,&
                            rsg%rmin,rsg%latmin,rsg%lonmin,rsg%rearth,'colsum','colsum_allf')
      end if
      rsgcols%field(:,1) = rsgcols%field(:,1)+rsg%field(:,icol)
      call rsg%dealloc()
      deallocate(res)
   end do
!
!  normalize colsums to maximum
!
   csmax = maxval(rsgcols%field(:,1))
   rsgcols%field(:,1) = rsgcols%field(:,1)/csmax
!
!  write output grid
!
   write(iter_path,"(2a,i3.3,a)") trim(main_path),trim(inpar_inv.sval.'ITERATION_STEP_PATH'),itlist(nit),'/'
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   write(rsgfile,'(a,a,a)') &
         trim(outpath)//trim(rsgpath)//'pseudo_mesh_'//trim(prop),'_colsum_allf','.hdf'
   print *,'Write collected column sums to ',trim(rsgfile)
   call writeHDFRegularSphericalGrid(rsgcols,rsgfile,errmsg)
   if (.level.errmsg == 2) goto 1
   call rsgcols%dealloc()
   call closeEnvironmentHDFWrapper(errmsg)
   call dealloc(invbasics)
   call dealloc(errmsg)
   call dealloc(ap)
!------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif
!
end program collectColsums
