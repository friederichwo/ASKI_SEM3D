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
!   Program to perform the solution of the least-squares
!   problem |(d|s)-(K|S)x|^2 -> min for kernel matrix
!   K and regularization matrix S using conjugate gradients.
!   The algorithm applied here is the "CGLS1" from paper:
!
!     Ake Bjoerck, Tommy Elfving and Zdenek Strakos
!     "Stability of conjugate gradient and Lanczos methods
!     for linear least squares problems"
!     SIAM J. Matrix Anal. APPL., 19, 720-736, 1998.'
!------------------------------------------------------------------------------
program solveCglsKernelSystem
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use dataModelSpaceInfo
   use linearModelRegularization
   use kernelLinearSystem
   use propertySet
   use kernelInvertedModel
   use vtkFile
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
   type (kernel_linear_system) :: kls
   type (linear_model_regularization) :: lmreg
   type (kernel_inverted_model) :: kim,kimcs
   type (vtk_info) :: vtk
   type (mpi_support) :: mpisup
   type (error_message) :: errmsg
   class (inversion_grid), allocatable :: invgrid
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar_inv
   type (input_parameter), pointer :: inpar_iter
   type (seismic_event_list), pointer :: evlist
   type (seismic_network), pointer :: statlist
   integer :: nsta,nlta,nproc,max_niter,kiter,nprop_set,nprop,ndata,ndata_all
   integer :: nreg,offreg,ncell,i,j,ios,idx,jdx
   integer, dimension(:), pointer :: ifreq_iter,ifreq_all
   double precision, dimension(:,:), pointer :: ref_model_values,corrmat
   double precision, dimension(:), allocatable :: vscal_smoothing_mantle,vscal_damping_mantle
   double precision, dimension(:), allocatable :: vscal_smoothing_crust,vscal_damping_crust
   double precision, dimension(:), allocatable :: vscal_smooth_mantle,vscal_damp_mantle
   double precision, dimension(:), allocatable :: vscal_smooth_crust,vscal_damp_crust
   double precision, dimension(:), allocatable :: norm_r_sta,norm_r_lta
   double precision :: dsmboost,boundchoke,crdepth,wx,wy
   double precision :: alfa,gamma,gamma_old,norm_sq_qd,norm_sq_qs,norm_sq_rd,norm_sq_rs,norm_sq_sol
   double precision :: sta,lta,lta_old,vtk_scale,calib_error_factor
   real, dimension(:), pointer :: rdata
   real, dimension(:,:), pointer :: absdev,colsum
   logical :: use_nonzero_starting_solution,init_only,do_phase_inversion
   character(len=char_len_par) :: prop
   character(len=char_len_par), dimension(:), allocatable :: prop_out
   character(len=char_len_par), dimension(:), pointer :: prop_inv
   character(len=max_length_string) :: main_parfile,dmsp_file_base,dmsp_file,startsol_file,outpath,path_dmsp
   character(len=max_length_string) :: path_measured_data,path_synthetic_data,path_sensitivity_kernels,iterpath
   character(len=max_length_string) :: vtk_format,vtkbasename,vtkpath,vtk_geometry_type,main_path,invgrid_type
   character(len=21) :: myname = 'solveCglsKernelSystem'
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
   call init(ap,myname,"Do inversion step by solving the kernel linear system defined by data model space info "//&
        "and regularization constraints in least squares optimization by parallelized conjugate-gradient method")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-startsol',.true.,'defines optional starting solution of linear system','sval','up.kim')
   call addOption(ap,'-init_only',.false.,'only read data and synthetics and produce kls_init data file')
   call parse(ap)
   if (myrank == 0) then
      call document(ap)
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap); call usage(ap)
         call add(errmsg,2,'Command line oculd not be read successfully',myname)
         goto 1
      endif
   endif
! -------------------------------------------------------------
!  get values of positional arguments
!
   main_parfile = ap.sval.'main_parfile'
   use_nonzero_starting_solution = ap.optset.'-startsol'
   if(use_nonzero_starting_solution) startsol_file = ap.sval.'-startsol'
   init_only = ap.optset.'-init_only'
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
   iterpath = getIterationStepPathIterationStepBasics(iterbasics)
   path_dmsp = trim(iterpath)//trim(inpar_inv.sval.'PATH_DMSPACE')
   dmsp_file_base = trim(path_dmsp)//'ASKI_dmspace'
   propset => getPropertySetInversionBasics(invbasics)
   evlist => getEventListInversionBasics(invbasics)
   statlist => getStationListInversionBasics(invbasics)
   path_measured_data = main_path+(inpar_inv.sval.'PATH_MEASURED_DATA')
   path_synthetic_data = iterpath+(inpar_inv.sval.'PATH_SYNTHETIC_DATA')
   path_sensitivity_kernels = iterpath+(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')
   outpath = iterpath+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   vtkpath = iterpath+(inpar_inv.sval.'PATH_VTK_FILES')
   ifreq_all => getMeasuredDataFrequencyIndicesInversionBasics(invbasics)
   wx = inpar_inv.dval.'ASKI_wx'
   wy = inpar_inv.dval.'ASKI_wy'
!
   do_phase_inversion = inpar_iter.lval.'DO_PHASE_INVERSION'
   calib_error_factor = inpar_iter.dval.'CALIBRATION_ERROR_FACTOR'
   vtk_scale = inpar_iter.dval.'VTK_COORDS_SCALING_FACTOR'
   vtk_geometry_type = inpar_iter.sval.'VTK_GEOMETRY_TYPE'
   ifreq_iter => getFrequencyIndicesIterationStepBasics(iterbasics)
   max_niter = inpar_iter.ival.'MAX_NUM_CG_ITERATIONS'
   nsta = inpar_iter.ival.'NITER_WINDOW_STA'
   nlta = inpar_iter.ival.'NITER_WINDOW_LTA'
   nproc = inpar_iter.ival.'NPROC'
   if (nproc /= numtasks) then
      call add(errmsg,2,'Number of actual processes conflicts with NPROC in iter parfile',myname)
      goto 1
   end if
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
   if (myrank == 0 .and. do_phase_inversion) print *,'Do a phase inversion'
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
! ----------------------------------------------------------------------
!  read my share of the data model space
!
   write(dmsp_file,"(a,'.',i3.3)") trim(dmsp_file_base),myrank
   call createFromFileDataModelSpaceInfo(dmspace,evlist,statlist,ifreq_iter,propset,ncell,trim(dmsp_file),1,errmsg)
   if (.level.errmsg == 2) goto 1
   if (myrank == 0) print *,'Data model space info read from file: ',trim(dmsp_file)
!
!  calculate total number of paths over all processes
!
   ndata = getNdataDataModelSpaceInfo(dmspace)
   call MPI_ALLREDUCE(ndata,ndata_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ios)
   if (ios .ne. MPI_SUCCESS) then
      call add(errmsg,2,'could not compute total number of paths',myname)
      goto 1
   end if
!
   if (myrank == 0) then
      print *,'Share of data of rank 0 is: ',ndata
      print *,'Total number of data is: ',ndata_all
   endif
! -------------------------------------------------------------------------------
!  Initalize kernel linear system, read my share of measured and synthetic data
!  and fill my share of kernel matrix and set rhs
!
   call initiateKernelLinearSystem(kls,dmspace,do_phase_inversion)
   if (myrank == 0) print *,'Initiated kernel system'
   call readMeasuredDataKernelLinearSystem(kls,dmspace,path_measured_data,ifreq_all,errmsg)
   if (.level.errmsg == 2) goto 1
   if (myrank == 0) print *,'Measured data read from path: ',trim(path_measured_data)
   call readSyntheticDataKernelLinearSystem(kls,dmspace,path_synthetic_data,ifreq_all,errmsg)
   if (.level.errmsg == 2) goto 1
   if (myrank == 0) print *,'Synthetic data read from path: ',trim(path_synthetic_data)
   call computeSigmaKernelLinearSystem(kls,dmspace,calib_error_factor)
   if (myrank == 0) print *,'Sigma calculated using: ',calib_error_factor
   call setRhsAsDataResidualKernelLinearSystem(kls)
   if (myrank == 0) print *,'Data residuals calculated '
!
!  if init_only option is set, write here KLS data and then exit code
!
   if (init_only) then
      call writeAllDataKernelLinearSystem(kls,dmspace,ndata,outpath,'init',errmsg)
      call abort(mpisup)
   endif
!
   call fillMatrixKernelLinearSystem(kls,dmspace,path_sensitivity_kernels,.setname.propset,ifreq_iter,errmsg)
   if (.level.errmsg == 2) goto 1
   if (myrank == 0) print *,'Kernel system matrix filled from path ',trim(path_sensitivity_kernels)
   if (use_nonzero_starting_solution) then
      call readHDFKernelInvertedModel(kim,startsol_file,errmsg)
      if (.level.errmsg == 2) goto 1
   end if
! ----------------------------------------------------------------------
!  Set my range of the model space for regularization
!  Equal sharing of invgrid cells, remainder is distributed again starting from rank zero
!  Regular share: nreg = ncell/numtasks
!  For 0 < myrank < nremain = mod(ncell,numtasks) take one more
!
   call shareLinearArraySmartUtils(ncell,numtasks,myrank,nreg,offreg)
   if (myrank == 0) print *,'Rank 0 gets ',nreg,' regularization equations at offset ',offreg
! --------------------------------------------------------------------------
!  Calculate indices and values in regularization equations
!
   nprop = .nprop.dmspace
   prop_inv => getPropertiesDataModelSpaceInfo(dmspace)
   if (.not. checkCorrmatPropertySet(propset,prop_inv,errmsg)) goto 1
   if (myrank == 0) print *,'Invert for properties: ',nprop,prop_inv
   if (myrank == 0) write(6,'(a15,2a15)') 'Prop','Smoothing','Damping'
   allocate(vscal_smooth_mantle(nprop),vscal_damp_mantle(nprop))
   allocate(vscal_smooth_crust(nprop),vscal_damp_crust(nprop))
   do j = 1,nprop
      idx = getIndexFromNamePropertySet(propset,prop_inv(j))
      vscal_smooth_mantle(j) = vscal_smoothing_mantle(idx)
      vscal_damp_mantle(j) = vscal_damping_mantle(idx)
      vscal_smooth_crust(j) = vscal_smoothing_crust(idx)
      vscal_damp_crust(j) = vscal_damping_crust(idx)
      if (myrank == 0) write(6,'(a9,a6,2e15.3)') 'Mantle:',trim(prop_inv(j)),vscal_smooth_mantle(j),vscal_damp_mantle(j)
      if (myrank == 0) write(6,'(a9,a6,2e15.3)') 'Crust: ',trim(prop_inv(j)),vscal_smooth_crust(j),vscal_damp_crust(j)
   end do
   if (myrank == 0) write(6,'(a,e15.3,a,e15.3)') 'Depth smoothing boost: ',dsmboost,' Boundary damping: ',boundchoke
   call initiateLinearModelRegularization(lmreg,nprop,vscal_smooth_mantle,vscal_damp_mantle,&
                                          vscal_smooth_crust,vscal_damp_crust,dsmboost,boundchoke,&
                                          offreg,nreg,errmsg)
   if (.level.errmsg == 2) goto 1
   if (myrank == 0) print *,'model regularization initiated'
   call smoothingLinearModelRegularization(lmreg,invgrid,crdepth)
   if (myrank == 0) print *,'Smoothing equations done '
   call dampingLinearModelRegularization(lmreg,invgrid,crdepth,wx,wy)
   if (myrank == 0) print *,'Damping equations done'
! ------------------------------------------------------------------------------
!  Solve kernel linear system using CGLS algorithm
!  Initialize iteration, gamma = -1.d0 indicates first iteration
!
   gamma = -1.d0
   kiter = 1
!
!  allocate arrays for sta/lta-monitoring of rhs-norm
!  fill array like a ring buffer restarting with index 1 if full
!
   allocate(norm_r_sta(nsta),norm_r_lta(nlta))
   norm_r_sta = 0.d0
   norm_r_lta = 0.d0
!
!  start with setting rhs = rhs - A*sol or keep rhs as is
!
   if (use_nonzero_starting_solution) then
      call setNonzeroSolutionKernelLinearSystem(kls,kim)
      call rhsMinusMatrixDotSolKernelLinearSystem(kls)
      call rhsMinusMatrixDotSolLinearModelRegularization(lmreg,kls%sol)
   else
      call setZeroSolutionKernelLinearSystem(kls)
   end if
!
!  write current data d and synthetics s = d-r to HDF file for initial state
!
   call writeAllDataKernelLinearSystem(kls,dmspace,ndata,outpath,'init',errmsg)
   if (.level.errmsg == 2) goto 1
!
   if (myrank == 0) then
      print *,'Start CGLS iterations'
      write(6,'(a8,7a15)') 'k_Iter','STA','LTA','NormSq_RD','NormSq_RS','Norm_Sol','Min_Sol','Max_Sol'
   endif
! -------------------------------------------------------------
!  Do CGLS iterations
!
   do while (kiter <= max_niter)
   !
   !  compute squared norm of current rhs (including allreduce of local norms)
   !
      call computeNormSquaredRhsKernelLinearSystem(kls,norm_sq_rd,errmsg)
      if (.level.errmsg == 2) goto 1
      call computeNormSquaredRhsLinearModelRegularization(lmreg,norm_sq_rs,errmsg)
      if (.level.errmsg == 2) goto 1
      call computeNormSquaredSolKernelLinearSystem(kls,norm_sq_sol)
      norm_r_sta(mod(kiter-1,nsta)+1) = sqrt(norm_sq_rd+norm_sq_rs)
      norm_r_lta(mod(kiter-1,nlta)+1) = sqrt(norm_sq_rd+norm_sq_rs)
      lta_old = lta
      sta = sum(norm_r_sta)/nsta
      lta = sum(norm_r_lta)/nlta
      if (myrank == 0) then
         if (do_phase_inversion) then
            write(6,'(i8,7d15.6)') kiter,sta,lta,2.d0*norm_sq_rd/ndata_all,norm_sq_rs/(2*ncell),&
                                   sqrt(norm_sq_sol/ncell),minval(kls%sol),maxval(kls%sol)
         else
            write(6,'(i8,7d15.6)') kiter,sta,lta,norm_sq_rd/ndata_all,norm_sq_rs/(2*ncell),&
                                   sqrt(norm_sq_sol/ncell),minval(kls%sol),maxval(kls%sol)
         end if
      end if
   !
   !  check convergence, i.e. sta > lta or previous lta < lta
   !
      if (kiter > max(nsta,nlta)) then
         if (sta/lta+1.d-2 > 1.d0) then
            if (myrank == 0) then
               write(6,'(a,a)') 'Iteration terminated because short-term average ',&
                                'of residual norm has reached the value of its long-term average'
            endif
            exit
         end if
         if (lta > lta_old) then
            if (myrank == 0) then
               write(6,'(a,a)') 'Iteration terminated because long-term average ',&
                                'of residual norm is not monotonically decreasing'
            endif
            exit
         end if
      end if
   !
   !  compute g = A^T rhs (including allreduce of local g vectors)
   !  and squared norm
   !
      gamma_old = gamma
      call computeGradIsMatrixTransRhsKernelLinearSystem(kls,errmsg)
      if (.level.errmsg == 2) goto 1
      call addGradIsMatrixTransRhsLinearModelRegularization(lmreg,kls%g,errmsg)
      if (.level.errmsg == 2) goto 1
      call computeNormSquaredGradKernelLinearSystem(kls,gamma)
   !
   !  compute p = g or update p = g + (gamma/gamma_old)*p
   !
      call updatePKernelLinearSystem(kls,gamma,gamma_old)
   !
   !  compute q = A p and squared norm of q (including allreduce)
   !
      call qIsMatrixDotPKernelLinearSystem(kls)
      call qIsMatrixDotPLinearModelRegularization(lmreg,kls%p)
      call computeNormSquaredQKernelLinearSystem(kls,norm_sq_qd,errmsg)
      if (.level.errmsg == 2) goto 1
      call computeNormSquaredQLinearModelRegularization(lmreg,norm_sq_qs,errmsg)
      if (.level.errmsg == 2) goto 1
   !
   !  update rhs = rhs - gamma/norm_sq_q * q
   !  update solution: sol = sol + gamma/norm_sq_q * p
   !
      alfa = gamma/(norm_sq_qd+norm_sq_qs)
      call updateRhsKernelLinearSystem(kls,alfa)
      call updateRhsLinearModelRegularization(lmreg,alfa)
      call updateSolKernelLinearSystem(kls,alfa)
      kiter = kiter+1
   end do
! -----------------------------------------------------------------------------------------
!  write current data d and synthetics s = d-r to HDF file for final state
!
   call writeAllDataKernelLinearSystem(kls,dmspace,ndata,outpath,'final',errmsg)
   if (.level.errmsg == 2) goto 1
!
!  compute column sums of kernel matrix (collected by rank 0 only)
!
   call computeColumnSumsKernelLinearSystem(kls,colsum,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  since all processes have the identical solution, the following is done on one rank only
!  calculate absolute model deviations honouring property correlation
!
   if (myrank == 0) then
   !
   !  calculate absolute model deviations honouring property correlation
   !
      corrmat => getCorrmatPropertySet(propset)
      ref_model_values => invgrid%getModelValues()
      call getRequiredKernelsByNamePropertySet(propset,prop_out)
      allocate(absdev(ncell,size(prop_out)))
      do i = 1,size(prop_out)
         idx = propset.index.prop_out(i)
         print *,'Compute model deviations for property: ',trim(prop_out(i)),' with index ',idx,' in property set'
         absdev(:,i) = 0.d0
      !
      !  Loop over properties inverted for and
      !  add here to the relative deviations of prop_out(i)
      !  the contribution due to solution for property prop_inv(j)
      !
         do j = 1,size(prop_inv)
            jdx = propset.index.prop_inv(j)
            absdev(:,i) = absdev(:,i) + ref_model_values(:,idx)*corrmat(idx,jdx) &
                        & *kls%sol((j-1)*ncell+1:j*ncell)/ref_model_values(:,jdx)
            print *,'Add contribution of property ',trim(prop_inv(j)),' with index ',jdx,' in property set'
            print *,'and correlation: ',corrmat(idx,jdx)
            print *,'Min and max of solution vector for ',trim(prop_inv(j)),': ',&
                   &  minval(kls%sol((j-1)*ncell+1:j*ncell)),maxval(kls%sol((j-1)*ncell+1:j*ncell))
         end do
      end do
      call createFromValuesKernelInvertedModel(kim,prop_out,absdev,.setname.propset,'absdev',errmsg)
      if (.level.errmsg == 2) goto 1
      nullify(absdev)
   !
   !  produce a kim file filled with column sum values
   !
      call createFromValuesKernelInvertedModel(kimcs,['cs    '],colsum,'None','colsum',errmsg)
      if (.level.errmsg == 2) goto 1
      nullify(colsum)
   !
   !  write kernel inverted model to files
   !
      call writeHDFKernelInvertedModel(kim,trim(outpath)//'inverted_absdev_model.kim',invgrid,errmsg)
      if (.level.errmsg == 2) goto 1
      call writeHDFKernelInvertedModel(kimcs,trim(outpath)//'column_sums.kim',invgrid,errmsg)
      if (.level.errmsg == 2) goto 1
   !
   !  write VTK files
   !
      vtk_format = trim(inpar_inv.sval.'DEFAULT_VTK_FILE_FORMAT')
      call initiateVtkFile(vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
      if (.level.errmsg == 2) goto 1
      do while (nextPropKernelInvertedModel(kim,prop,iprop = j))
         rdata => getValuesPropKernelInvertedModel(kim,j)
         write(vtkbasename,'(a,a,a)') trim(vtkpath),'inverted_absdev_model_',trim(prop)
         call writeRealDataVtkFile(vtk,trim(vtkbasename),1,rdata,errmsg,data_name='inv_absdev_'//trim(prop))
         if (.level.errmsg == 2) goto 1
      enddo
      call dealloc(vtk)
      call dealloc(kim)
      call dealloc(kimcs)
      deallocate(prop_out)
      print *,'Model written to HDF and VTK files'
   endif
!
!  clean up
!
   call barrier(mpisup)
   deallocate(norm_r_sta,norm_r_lta)
   call dealloc(lmreg)
   deallocate(vscal_smooth_mantle,vscal_damp_mantle)
   deallocate(vscal_smooth_crust,vscal_damp_crust)
   nullify(prop_inv)
   call dealloc(kls)
   call dealloc(dmspace)
   deallocate(vscal_smoothing_mantle,vscal_damping_mantle)
   deallocate(vscal_smoothing_crust,vscal_damping_crust)
   nullify(ifreq_iter,ifreq_all)
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
      if (myrank == 0) call print(errmsg)
      call abort(mpisup)
   endif

end program solveCglsKernelSystem
