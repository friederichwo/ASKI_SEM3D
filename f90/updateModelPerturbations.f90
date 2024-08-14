program updateModelPerturbation
   use inversionBasics
   use kernelInvertedModel
   use specfem3dInversionGrid
   use inputParameter
   use argumentParser
   use string
   use errorMessage
   use globalMpiInfo
   implicit none

   type (argument_parser) :: ap
   type (inversion_basics) :: invbasics
   type (error_message) :: errmsg
   type (input_parameter), pointer :: inpar_inv
   type (kernel_inverted_model) :: delmv,totmv,newmv
   class (inversion_grid), allocatable :: invgrid
   integer :: it,ierr
   real, dimension(:,:), pointer :: new_model_values
   character(len=max_length_string) :: main_parfile,outpath,outpath_prev,iterpath,iterpath_prev
   character(len=max_length_string) :: main_path,invgrid_type
   character(len=24) :: myname = 'updateModelPerturbations'
!  --------------------------------------------------------------------------
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
   call init(ap,myname,"Update total model perturbation of previous iteration")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-it',.true.,'iteration whose model deviations are taken for update','ival','0')
   call parse(ap)
   call document(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap); goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   it = ap.ival.'-it'
!  ------------------------------------------------------------------------
!  setup inversion basics
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
!
!  iteration specific stuff
!
   if (it == 0) then
      iterpath = .iterpath.invbasics
      it = inpar_inv.ival.'CURRENT_ITERATION_STEP'
   else
      write(iterpath,"(2a,i3.3,a)") trim(main_path),trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it,'/'
   end if
   outpath = iterpath+(inpar_inv.sval.'PATH_OUTPUT_FILES')
!
!  if it = 1, just copy the inverted_absdev_model_kim to toal_absdev_model_kim
!
   if (it == 1) then
      print *,'Copied inverted perturbations to total perturbations (done for it=1)'
      call system('cp -f '//trim(outpath)//'inverted_absdev_model.kim'//' '//trim(outpath)//'total_absdev_model.kim',ierr)
      call dealloc(invbasics)
      goto 1
   end if
!
!  path to previous iteration
!
   if (it > 1) then
      write(iterpath_prev,"(2a,i3.3,a)") trim(main_path),trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it-1,'/'
      outpath_prev = iterpath_prev+(inpar_inv.sval.'PATH_OUTPUT_FILES')
      print *,'Path to previous iteration: ',trim(iterpath_prev)
   end if
!
!  read model perturbations
!
   call readHDFKernelInvertedModel(delmv,trim(outpath)//'inverted_absdev_model.kim',errmsg)
   if(.level.errmsg == 2) goto 1
   call readHDFKernelInvertedModel(totmv,trim(outpath_prev)//'total_absdev_model.kim',errmsg)
   if(.level.errmsg == 2) goto 1
!
!  add model perturbations
!
   if (delmv%nprop .ne. totmv%nprop) then
      call add(errmsg,2,'Models have different number of properties',myname)
      goto 1
   end if
   if (delmv%ncell .ne. totmv%ncell) then
      call add(errmsg,2,'Models have different number of cell',myname)
      goto 1
   end if
!
   allocate(new_model_values(delmv%ncell,delmv%nprop))
   new_model_values = delmv%model_values + totmv%model_values
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
   call createFromValuesKernelInvertedModel(newmv,delmv%prop,new_model_values,delmv%propsetname,delmv%value_kind,errmsg)
   if(.level.errmsg == 2) goto 1
   call writeHDFKernelInvertedModel(newmv,trim(outpath)//'total_absdev_model.kim',invgrid,errmsg)
   if(.level.errmsg == 2) goto 1
   nullify(new_model_values)
!
!  clean up
!
   call dealloc(ap)
   call dealloc(invbasics)
   call dealloc(delmv)
   call dealloc(totmv)
   call dealloc(newmv)
   call invgrid%dealloc(); deallocate(invgrid)
!  ------------------------------------------------------------------------------
!  error handling
!
 1 if (.level.errmsg == 2) then
      if (myrank == 0) call print(errmsg)
   endif
end program updateModelPerturbation
