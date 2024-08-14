program computeMisfits
   use inversionBasics
   use dataModelSpaceInfo
   use kernelLinearSystem
   use propertySet
   use inputParameter
   use argumentParser
   use hdfWrapper
   use mpiSupport
   use string
   use mathConstants
   use errorMessage
   use globalHdfInfo
   use globalMpiInfo
   implicit none

   type (argument_parser) :: ap
   type (inversion_basics) :: invbasics
   type (data_model_space_info) :: dmspace
   type (kernel_linear_system) :: kls
   type (property_set), pointer :: propset
   type (error_message) :: errmsg
   type (mpi_support) :: mpisup
   type (any_rank_real_array) :: arra
   type (any_rank_integer_array) :: aria
   type (input_parameter), pointer :: inpar_inv
   type (input_parameter) :: inpar_iter
   type (seismic_event_list), pointer :: evlist
   type (seismic_network), pointer :: statlist
   integer :: i,j,ndata,nev,nstat,nf,is,ie,jf,iter,npath,ncell,badcount
   integer(kind=8) :: fid
   integer, dimension(:), allocatable :: nsevfr
   integer, dimension(:), pointer :: ifreq_all,ifreq_ds,pathidx_ds,griddim
   real, dimension(:,:,:), allocatable :: absres,timeres,sigma,dnorm
   real, dimension(:), allocatable :: res_real,res_imag,resabs,restime
   double precision, dimension(:), allocatable :: mf
   double precision :: total_mf,df,phidat,phisyn,phidiff,calib_error_factor,timeres_max
   logical :: use_masked_dmspace,badflag,use_allf,use_ext
   character(len=char_len_evid), dimension(:), pointer :: evid_ds
   character(len=char_len_netcode+char_len_sta+1), dimension(:), pointer :: netstaname_ds
   character(len=char_len_evid), dimension(:), allocatable :: eventid
   character(len=char_len_netcode+char_len_sta+1), dimension(:), allocatable :: nsname
   character(len=max_length_string) :: main_parfile,dmsp_file,outpath,path_dmsp,logfile
   character(len=max_length_string) :: path_measured_data,path_synthetic_data,iterpath,iter_parfile
   character(len=max_length_string) :: main_path,outfile,fmt,dmspext
   character(len=14) :: myname = 'computeMisfits'
   character (len=80), dimension(4) :: iter_inpar_keys
   data iter_inpar_keys/'ITERATION_STEP_NUMBER_OF_FREQ', 'ITERATION_STEP_INDEX_OF_FREQ', 'CALIBRATION_ERROR_FACTOR',&
                        'MAX_TIME_RESIDUAL'/
!  --------------------------------------------------------------------------
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
!
!  command line processing
!
   call init(ap,myname,"Compute various kinds of misfits from FT measured and synthetic data")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-iter',.true.,'iteration for which residuals are computed','ival','0')
   call addOption(ap,'-masked',.false.,'Use the masked dmspace file')
   call addOption(ap,'-allf',.false.,'Use a dmspace file with all frequencies')
   call addOption(ap,'-ext',.true.,'Use a dmspace file with given extension','sval','None')
   call parse(ap)
   call document(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line oculd not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   iter = ap.ival.'-iter'
   use_masked_dmspace = ap.optset.'-masked'
   use_allf = ap.optset.'-allf'
   use_ext = ap.optset.'-ext'
   dmspext = ap.sval.'-ext'
!  ------------------------------------------------------------------------
!  setup inversion basics
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
   evlist => getEventListInversionBasics(invbasics)
   statlist => getStationListInversionBasics(invbasics)
   path_measured_data = main_path+(inpar_inv.sval.'PATH_MEASURED_DATA')
   df = getMeasuredDataFrequencyStepInversionBasics(invbasics)
   ifreq_all => getMeasuredDataFrequencyIndicesInversionBasics(invbasics)
   propset => getPropertySetInversionBasics(invbasics)
   griddim => ivecp(inpar_inv,'INVGRID_DIMENSIONS',3)
   ncell = product(griddim)
!
!  iteration specific stuff
!
   if (iter == 0) then
      iterpath = .iterpath.invbasics
   else
      write(iterpath,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),iter,'/'
   end if
!
   path_dmsp = trim(iterpath)//trim(inpar_inv.sval.'PATH_DMSPACE')
   dmsp_file = trim(path_dmsp)//'ASKI_dmspace'
   if (use_masked_dmspace) dmsp_file = trim(path_dmsp)//'ASKI_dmspace_masked'
   if (use_allf) dmsp_file = trim(path_dmsp)//'ASKI_dmspace_allf'
   if (use_ext) dmsp_file = trim(path_dmsp)//'ASKI_dmspace_'//trim(dmspext)
   if (use_allf .and. use_masked_dmspace) dmsp_file = trim(path_dmsp)//'ASKI_dmspace_masked_allf'
   if (use_ext .and. use_masked_dmspace) dmsp_file = trim(path_dmsp)//'ASKI_dmspace_masked_'//trim(dmspext)
!
   path_synthetic_data = iterpath+(inpar_inv.sval.'PATH_SYNTHETIC_DATA')
   outpath = iterpath+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   iter_parfile = trim(iterpath)//trim(inpar_inv.sval.'PARFILE_ITERATION_STEP')
   call createKeywordsInputParameter(inpar_iter,iter_inpar_keys)
   call readSubroutineInputParameter(inpar_iter,1,trim(iter_parfile),errmsg)
   if(.level.errmsg == 2) return
   nf = ival(inpar_inv,'MEASURED_DATA_NUMBER_OF_FREQ')
   calib_error_factor = inpar_iter.dval.'CALIBRATION_ERROR_FACTOR'
   timeres_max = inpar_iter.dval.'MAX_TIME_RESIDUAL'
!
   logfile = trim(outpath)//'misfits.log'
   if (use_masked_dmspace) logfile = trim(outpath)//'misfits_masked.log'
   if (use_allf) logfile = trim(outpath)//'misfits_allf.log'
   if (use_ext) logfile = trim(outpath)//'misfits_'//trim(dmspext)//'.log'
   if (use_allf .and. use_masked_dmspace) logfile = trim(outpath)//'misfits_masked_allf.log'
   if (use_ext .and. use_masked_dmspace) logfile = trim(outpath)//'misfits_masked_'//trim(dmspext)//'.log'
   open(1,file = trim(logfile))
   write(1,'(a,a)') 'do residual calculation for: ',trim(iterpath)
!  ----------------------------------------------------------------------
!  read the data model space
!
   call createDataSamplesFromFileDataModelSpaceInfo(dmspace,evlist,statlist,ifreq_all,trim(dmsp_file),7,errmsg)
   if (.level.errmsg == 2) goto 1
   call createModelValuesFromFileDataModelSpaceInfo(dmspace,propset,ncell,trim(dmsp_file),7,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,a)') 'Data model space info read from file: ',trim(dmsp_file)
!
   ndata = getNdataDataModelSpaceInfo(dmspace)
   nev = .nev.evlist
   nstat = .nstat.statlist
   npath = getNpathDataModelSpaceInfo(dmspace)
!
   write(1,'(a,i8)') 'Total number of paths is: ',npath
   write(1,'(a,i8)') 'Total number of data is: ',ndata
   write(1,'(a,i8)') 'Number of events is: ',nev
   write(1,'(a,i8)') 'Number of stations is: ',nstat
   write(1,'(a,i8)') 'Number of valif frequencies is: ',nf
   write(1,*) 'Valid frequency indices: ',ifreq_all
!  -------------------------------------------------------------------------------
!  Initalize kernel linear system, read my share of measured and synthetic data
!
   call initiateDataKernelLinearSystem(kls,dmspace)
   write(1,'(a,a)') 'Initiated kernel system'
   call readMeasuredDataKernelLinearSystem(kls,dmspace,path_measured_data,ifreq_all,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,a)') 'Measured data read from path: ',trim(path_measured_data)
   call readSyntheticDataKernelLinearSystem(kls,dmspace,path_synthetic_data,ifreq_all,errmsg)
   if (.level.errmsg == 2) goto 1
   write(1,'(a,a)') 'Synthetic data read from path: ',trim(path_synthetic_data)
   call computeSigmaKernelLinearSystem(kls,dmspace,calib_error_factor)
!  -----------------------------------------------------------------------------
!  produce arrays of event ids and netstanames
!
   allocate(eventid(nev))
   do i = 1,nev
	   eventid(i) = .evid.(evlist.event.i)
   end do
   allocate(nsname(nstat))
   do i = 1,nstat
	   nsname(i) = .nsname.(statlist.station.i)
   end do
!  -----------------------------------------------------------------------------
!  compute absolute data-synthetic residual resolved for event, netstaname and frequency index
!  assumes only one component in data
!
   allocate(absres(nstat,nev,nf),timeres(nstat,nev,nf))
   allocate(sigma(nstat,nev,nf),dnorm(nstat,nev,nf))
   absres = -1.0; timeres = -1000.0                      ! set to impossible values if station is not available
   sigma = -1.0; dnorm = -1.0
!
!  flat arrays for residuals ordered according to data index
!  for overall histogram plotting
!
   allocate(res_real(ndata/2),res_imag(ndata/2),resabs(ndata/2),restime(ndata/2))
!
!  get arrays of event ids, netstanames, frequency and path indices in dmspace
!
   evid_ds => getArrayEvidDataSamplesDataModelSpaceInfo(dmspace)
   netstaname_ds => getArrayNetStanameDataSamplesDataModelSpaceInfo(dmspace)
   ifreq_ds => getIfreqArrayDataModelSpaceInfo(dmspace)
   pathidx_ds => getArrayPathidxDataSamplesDataModelSpaceInfo(dmspace)
!
!  Loop over data samples as ordered in dmspace
!
   total_mf = 0.d0
   badcount = 0
   do i = 1,ndata,2                                                ! data are always real and imaginary part
      j = (i-1)/2+1
      badflag = .false.
      is = findloc(nsname,netstaname_ds(i),1)                      ! index of sample nsname in station list
      ie = findloc(eventid,evid_ds(i),1)                           ! index of sample evid in event list
      jf = findloc(ifreq_all,ifreq_ds(i),1)                        ! index of sample frequency in iter frequency list
      res_real(j) = (kls%mdata(i)-kls%sdata(i))/kls%sigma(i)
      res_imag(j) = (kls%mdata(i+1)-kls%sdata(i+1))/kls%sigma(i)
      resabs(j) = hypot(res_real(j),res_imag(j))
      phidat = atan2(kls%mdata(i+1),kls%mdata(i))
      phisyn = atan2(kls%sdata(i+1),kls%sdata(i))
      phidiff = phidat-phisyn                                      ! phase difference is not the phase of res = d-s
      if (phidiff > mc_pi)  phidiff = phidiff-2.*mc_pi             ! bring phase difference into range -pi to +pi
      if (phidiff < -mc_pi) phidiff = phidiff+2.*mc_pi
      restime(j) = phidiff/(2.*mc_pi*df*ifreq_ds(i))
      total_mf = total_mf+resabs(j)**2                             ! update total error normalised quadratic misift
   !
   !  fill output arrays
   !
      absres(is,ie,jf) = resabs(j)
      timeres(is,ie,jf) = phidiff/(2.*mc_pi*df*ifreq_ds(i))
      sigma(is,ie,jf) = kls%sigma(i)
      dnorm(is,ie,jf) = kls%dnorm(i)
   !
   !  check for stations with resabs/sigma > 8 or abs(timeres) > 5 s
   !  and mask them in dmspace
   !
      if (absres(is,ie,jf) > 0.8/calib_error_factor) then
         call setPathFrequencyMaskDataModelSpaceInfo(dmspace,pathidx_ds(i),ifreq_ds(i))
         write(1,'(a25,a20,a25,i6,f8.2)') 'Residual exceeded: ',trim(netstaname_ds(i)),trim(evid_ds(i)), &
                                          ifreq_ds(i),resabs(j)
         badflag = .true.
      end if
      if (abs(timeres(is,ie,jf)) > timeres_max) then
         call setPathFrequencyMaskDataModelSpaceInfo(dmspace,pathidx_ds(i),ifreq_ds(i))
         write(1,'(a25,a20,a25,i6,f8.2)') 'Time shift exceeded: ',trim(netstaname_ds(i)),trim(evid_ds(i)), &
                                          ifreq_ds(i),abs(timeres(is,ie,jf))
         badflag = .true.
      end if
      if (badflag) badcount = badcount+1
   end do
   write(1,'(a,2f15.3)') 'Total quadratic / rms misfit: ',total_mf/ndata, sqrt(total_mf/ndata)
   write(1,'(a,2i8)') 'Total number of data / bad data ',ndata,2*badcount
!
!  compute cumulated misfit per event and frequency
!
   write(1,'(a)') 'Event-frequency-wise misfit'
   write(fmt,'(a,i3,a)') '(a8,a20,',nf,'f10.3)'
   write(1,trim(fmt)) 'Index','Event',(ifreq_all(jf)*df, jf = 1,nf)
   write(fmt,'(a,i3,a)') '(i8,a20,',nf,'f10.3)'
   allocate(mf(nf),nsevfr(nf))
   do ie = 1,nev
      do jf = 1,nf
         mf(jf) = 0.d0
         nsevfr(jf) = 0
         do is = 1,nstat
            if (absres(is,ie,jf) > 0.d0) then
               mf(jf) = mf(jf)+absres(is,ie,jf)**2
               nsevfr(jf) = nsevfr(jf)+1
            end if
         end do
         if (nsevfr(jf) == 0) nsevfr(jf) = 1
      end do
      write(1,trim(fmt)) ie,eventid(ie),(0.5*mf(jf)/nsevfr(jf), jf = 1,nf)    ! real and imag part
   end do
   close(1)
!
!  write new data model space info file with masked data
!
   outfile = trim(path_dmsp)//'ASKI_dmspace_masked'
   if (use_allf) outfile = trim(path_dmsp)//'ASKI_dmspace_masked_allf'
   if (use_ext) outfile = trim(path_dmsp)//'ASKI_dmspace_masked_'//trim(dmspext)
   call writeDataModelSpaceInfo(dmspace,2,outfile,1,npath,errmsg)
!  ------------------------------------------------------------------------------
!  write station, event, frequency specific residual information to HDF file
!
   outfile = trim(outpath)//'residuals_specific.hdf'
   if (use_masked_dmspace) outfile = trim(outpath)//'residuals_specific_masked.hdf'
   if (use_allf) outfile = trim(outpath)//'residuals_specific_allf.hdf'
   if (use_ext) outfile = trim(outpath)//'residuals_specific_'//trim(dmspext)//'.hdf'
   if (use_allf .and. use_masked_dmspace) outfile = trim(outpath)//'residuals_specific_masked_allf.hdf'
   if (use_ext .and. use_masked_dmspace) outfile = trim(outpath)//'residuals_specific_masked_'//trim(dmspext)//'.hdf'
   call createFileHDFWrapper(trim(outfile),fid,errmsg)
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc3d(absres)
   call writeArrayHDFWrapper(fid,"absres",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc3d(timeres)
   call writeArrayHDFWrapper(fid,"timeres",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc3d(sigma)
   call writeArrayHDFWrapper(fid,"sigma",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc3d(dnorm)
   call writeArrayHDFWrapper(fid,"dnorm",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call aria%assoc1d(ifreq_all)
   call writeArrayHDFWrapper(fid,"ifreq",aria,errmsg)
   call aria%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) goto 1
!  ------------------------------------------------------------------------------
!  write flattened residual information to HDF file
!
   outfile = trim(outpath)//'residuals_flattened.hdf'
   if (use_masked_dmspace) outfile = trim(outpath)//'residuals_flattened_masked.hdf'
   if (use_allf) outfile = trim(outpath)//'residuals_flattened_allf.hdf'
   if (use_ext) outfile = trim(outpath)//'residuals_flattened_'//trim(dmspext)//'.hdf'
   if (use_allf .and. use_masked_dmspace) outfile = trim(outpath)//'residuals_flattened_masked_allf.hdf'
   if (use_ext .and. use_masked_dmspace) outfile = trim(outpath)//'residuals_flattened_masked_'//trim(dmspext)//'.hdf'
   call createFileHDFWrapper(trim(outfile),fid,errmsg)
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc1d(res_real)
   call writeArrayHDFWrapper(fid,"res_real",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc1d(res_imag)
   call writeArrayHDFWrapper(fid,"res_imag",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc1d(resabs)
   call writeArrayHDFWrapper(fid,"resabs",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call arra%assoc1d(restime)
   call writeArrayHDFWrapper(fid,"restime",arra,errmsg)
   call arra%deassoc()
   if (.level.errmsg == 2) goto 1
!
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) goto 1
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
!
   deallocate(absres,timeres,sigma,dnorm,eventid,nsname,res_real,res_imag,resabs,restime)
   deallocate(mf,nsevfr)
   call dealloc(kls)
   call dealloc(dmspace)
   nullify(ifreq_all,griddim)
   call dealloc(inpar_iter)
   call dealloc(invbasics)
   call dealloc(mpisup)
!  ------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif
!
end program computeMisfits
