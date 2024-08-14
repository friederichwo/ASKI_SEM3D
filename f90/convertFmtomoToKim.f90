!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!  This code may serves for:
!  1) calculating perturbations of vp and/or vs that can be used
!     in a SPECFEM simulation run or
!  2) serve as starting solution of an FWI iteration.
!
program convertFmtomoToKim
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use propertySet
   use kernelInvertedModel
   use askiBackgroundModel
   use fmtomoModel
   use axesRotation
   use argumentParser
   use inputParameter
   use string
   use errorMessage
   use mpiSupport
   use globalMpiInfo
   use globalHdfInfo

   implicit none

   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (kernel_inverted_model) :: kim
   type (aski_background_model) :: abm
   type (fmtomo_model) :: fmtomo
   type (mpi_support) :: mpisup
   type (input_parameter), pointer :: inpar_inv,inpar_iter
   class (inversion_grid), allocatable :: invgrid
   type (property_set), pointer :: propset
   double precision :: thetac,phic,rearth
   character(len=max_length_string) :: main_parfile,fmtomo_file,abm_file,outpath,iterpath
   character(len=max_length_string) :: kimfile,str,invgrid_type,path_specfem_input
   character(len=19) :: myname = 'convertFmtomoToKim'
!-------------------------------------------------------------------------------
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
!------------------------------------------------------------------------
!  command line processing
!
   call init(ap,myname,"Convert a FMTOMO 3D absolute model to a kernel_inverted_model (.kim)")
   call addPosarg(ap,"mainparfile","sval","name of main parameter file")
   call addPosarg(ap,"mfile","sval","FMTOMO input file of specified type")
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); goto 1; end if
   call document(ap)
!
   main_parfile = ap.sval.'mainparfile'
   fmtomo_file = ap.sval.'mfile'
   call dealloc(ap)
!------------------------------------------------------------------------
!  Setup inversion basics
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter => getInputParameterIterationStepBasics(iterbasics)
   iterpath = .iterpath.invbasics
   propset => getPropertySetInversionBasics(invbasics)
   thetac = 0.5*mc_pid-(inpar_inv.dval.'INVGRID_CENTER_LAT')*mc_deg2radd
   phic = (inpar_inv.dval.'INVGRID_CENTER_LON')*mc_deg2radd
   rearth = inpar_inv.dval.'REARTH'
   path_specfem_input = inpar_inv.sval.'PATH_SPECFEM_INPUT'
   abm_file = trim(path_specfem_input)//(inpar_inv.sval.'FILE_ASKI_BACKGROUND_MODEL')
   outpath = iterpath + (inpar_inv.sval.'PATH_OUTPUT_FILES')
   str = getBehindLastSeparatorString(trim(fmtomo_file),'/',errmsg)
   if (.level.errmsg == 2) goto 1
   kimfile = trim(str)//'.kim'
   print *,'Take aski background model from file: ',trim(abm_file)
   print *,'Take fmtomo model from file: ',trim(fmtomo_file)
   print *,'Write kim to file: ',trim(outpath)//trim(kimfile)
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(iterpath,.false.,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
!
!  compute kernel inverted model
!
   call readASKIBackgroundModel(abm,abm_file,1,rearth,errmsg)
   if (.level.errmsg == 2) goto 1
   call readFmtomoModel(fmtomo,1,fmtomo_file,errmsg)
   if (.level.errmsg == 2) goto 1
   call createFromFmtomoAbmKernelInvertedModel(kim,abm,fmtomo,invgrid,propset,thetac,phic,errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  write kim to HDF file
!
   call writeHDFKernelInvertedModel(kim,trim(outpath)//trim(kimfile),invgrid,errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(abm)
   call dealloc(fmtomo)
   call dealloc(kim)
   call invgrid%dealloc()
   deallocate(invgrid)
   call dealloc(invbasics)
   call dealloc(iterbasics)
   call dealloc(errmsg)
   call dealloc(mpisup)
!
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   end if
end program convertFmtomoToKim
