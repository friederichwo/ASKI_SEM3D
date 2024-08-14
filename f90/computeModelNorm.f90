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
!  Compute the model norm from N^2 = m^T S^T S m + m^T D^T D m
!------------------------------------------------------------------------------
program computeModelNorm
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use dataModelSpaceInfo
   use linearModelRegularization
   use propertySet
   use kernelInvertedModel
   use inputParameter
   use argumentParser
   use string
   use errorMessage
   use mpiSupport
   use smartUtils
   use globalMpiInfo
   use globalHdfInfo
   implicit none

   type (argument_parser) :: ap
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (data_model_space_info) :: dmspace
   type (linear_model_regularization) :: lmreg
   type (kernel_inverted_model) :: kim
   type (mpi_support) :: mpisup
   type (error_message) :: errmsg
   class (inversion_grid), allocatable :: invgrid
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar_inv
   type (input_parameter), pointer :: inpar_iter
   integer :: nprop_set,nprop
   integer :: nreg,offreg,ncell,j,idx,it,nmval
   double precision, dimension(:), allocatable :: vscal_smoothing_mantle,vscal_damping_mantle
   double precision, dimension(:), allocatable :: vscal_smoothing_crust,vscal_damping_crust
   double precision, dimension(:), allocatable :: vscal_smooth_mantle,vscal_damp_mantle
   double precision, dimension(:), allocatable :: vscal_smooth_crust,vscal_damp_crust
   double precision :: dsmboost,boundchoke,crdepth,wx,wy,mnorm2
   double precision, dimension(:), allocatable :: mval
   character(len=char_len_par), dimension(:), pointer :: prop_inv
   character(len=max_length_string) :: main_parfile,dmspfile,modelfile,outpath,path_dmsp,logfile
   character(len=max_length_string) :: iterpath
   character(len=max_length_string) :: main_path,invgrid_type,itchar
   character(len=16) :: myname = 'computeModelNorm'
! -------------------------------------------------------------------------------
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
!  command line processing
!
   call init(ap,myname,"Compute the norm of an inverted model")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-it',.true.,'Iteration from which inverted model is taken','sval','current')
   call parse(ap)
   call document(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line oculd not be read successfully',myname)
      goto 1
   endif
! -------------------------------------------------------------
!  get values of positional arguments
!
   main_parfile = ap.sval.'main_parfile'
   itchar = ap.sval.'-it'
! ------------------------------------------------------------------------
!  setup inversion basics
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter => getInputParameterIterationStepBasics(iterbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
!
   if (equalString(itchar,'current')) then
      it = inpar_inv.ival.'CURRENT_ITERATION_STEP'
      iterpath = .iterpath.invbasics
   else
      read(itchar,'(i3)') it
      write(iterpath,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it,'/'
   end if
   path_dmsp = trim(iterpath)//trim(inpar_inv.sval.'PATH_DMSPACE')
   dmspfile = trim(path_dmsp)//'ASKI_dmspace_masked'
   propset => getPropertySetInversionBasics(invbasics)
   outpath = iterpath+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   modelfile = trim(outpath)//'total_absdev_model.kim'
   wx = inpar_inv.dval.'ASKI_wx'
   wy = inpar_inv.dval.'ASKI_wy'
!
   nprop_set = .nprop.propset
   allocate(vscal_smoothing_mantle(nprop_set),vscal_damping_mantle(nprop_set))
   allocate(vscal_smoothing_crust(nprop_set),vscal_damping_crust(nprop_set))
   vscal_smoothing_mantle = dvec(inpar_iter,'VSCAL_SMOOTHING_MANTLE',nprop_set)
   vscal_smoothing_crust = dvec(inpar_iter,'VSCAL_SMOOTHING_CRUST',nprop_set)
   vscal_damping_mantle = dvec(inpar_iter,'VSCAL_DAMPING_MANTLE',nprop_set)
   vscal_damping_crust = dvec(inpar_iter,'VSCAL_DAMPING_CRUST',nprop_set)
   dsmboost = inpar_iter.dval.'DSMBOOST'
   boundchoke = inpar_iter.dval.'BOUNDCHOKE'
   crdepth = inpar_iter.dval.'CRUSTAL_DEPTH'
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
   ncell = invgrid%getNcellAll()
!
!  open output file
!
   logfile = trim(outpath)//'modelnorm.log'
   open(1,file = trim(logfile))
   write(1,'(a,a)') 'do model norm calculation for: ',trim(iterpath)
! ----------------------------------------------------------------------
!  read my share of the data model space
!
   call createModelValuesFromFileDataModelSpaceInfo(dmspace,propset,ncell,dmspfile,2,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,a)') 'Model space info read from file: ',trim(dmspfile)
   call readHDFKernelInvertedModel(kim,modelfile,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,a)') 'Total model perturbations read from file: ',trim(modelfile)
! ----------------------------------------------------------------------
!  Set my range of the model space for regularization
!  Equal sharing of invgrid cells, remainder is distributed again starting from rank zero
!  Regular share: nreg = ncell/numtasks
!  For 0 < myrank < nremain = mod(ncell,numtasks) take one more
!
   call shareLinearArraySmartUtils(ncell,1,0,nreg,offreg)
   write(1,'(i8,a,i8)')  nreg,' regularization equations at offset ',offreg
! --------------------------------------------------------------------------
!  Calculate indices and values in regularization equations
!
   nprop = .nprop.dmspace
   prop_inv => getPropertiesDataModelSpaceInfo(dmspace)
   if (.not. checkCorrmatPropertySet(propset,prop_inv,errmsg)) goto 1
   write(1,*) 'Inverted for properties: ',nprop,prop_inv
   write(1,'(a15,2a15)') 'Prop','Smoothing','Damping'
   allocate(vscal_smooth_mantle(nprop),vscal_damp_mantle(nprop))
   allocate(vscal_smooth_crust(nprop),vscal_damp_crust(nprop))
   do j = 1,nprop
      idx = getIndexFromNamePropertySet(propset,prop_inv(j))
      vscal_smooth_mantle(j) = vscal_smoothing_mantle(idx)
      vscal_damp_mantle(j) = vscal_damping_mantle(idx)
      vscal_smooth_crust(j) = vscal_smoothing_crust(idx)
      vscal_damp_crust(j) = vscal_damping_crust(idx)
      write(1,'(a9,a6,2e15.3)') 'Mantle:',trim(prop_inv(j)),vscal_smooth_mantle(j),vscal_damp_mantle(j)
      write(1,'(a9,a6,2e15.3)') 'Crust: ',trim(prop_inv(j)),vscal_smooth_crust(j),vscal_damp_crust(j)
   end do
   write(1,'(a,e15.3,a,e15.3)') 'Depth smoothing boost: ',dsmboost,' Boundary damping: ',boundchoke
   call initiateLinearModelRegularization(lmreg,nprop,vscal_smooth_mantle,vscal_damp_mantle,&
                                          vscal_smooth_crust,vscal_damp_crust,dsmboost,boundchoke,&
                                          offreg,nreg,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a)') 'model regularization initiated'
   call smoothingLinearModelRegularization(lmreg,invgrid,crdepth)
   write(1,'(a)') 'Smoothing equations done '
   call dampingLinearModelRegularization(lmreg,invgrid,crdepth,wx,wy)
   write(1,'(a)') 'Damping equations done'
!
!  flatten model values
!  reshape cannot be used because of type conversion real to double
!
   nmval = .nmval.dmspace
   mval = doubleFlattenedValuesKernelInvertedModel(kim)
!
!  compute model norm
!
   call qIsMatrixDotPLinearModelRegularization(lmreg,mval)
   call computeNormSquaredQLinearModelRegularization(lmreg,mnorm2,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,2f15.3)') 'Quadratic / rms model norm: ',mnorm2/nmval, sqrt(mnorm2/nmval)
   write(6,'(a,2f15.3)') 'Quadratic / rms model norm: ',mnorm2/nmval, sqrt(mnorm2/nmval)
   close(1)
!
!  clean up
!
   deallocate(mval)
   call dealloc(lmreg)
   deallocate(vscal_smooth_mantle,vscal_damp_mantle)
   deallocate(vscal_smooth_crust,vscal_damp_crust)
   nullify(prop_inv)
   call dealloc(dmspace)
   deallocate(vscal_smoothing_mantle,vscal_damping_mantle)
   deallocate(vscal_smoothing_crust,vscal_damping_crust)
   call invgrid%dealloc(); deallocate(invgrid)
   call dealloc(iterbasics)
   call dealloc(invbasics)
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(mpisup)
!------------------------------------------------------------------------------
!  error handling
!
 1 if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   endif
!
end program computeModelNorm
