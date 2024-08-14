!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!             2021 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   Program for transforming time-domain SPECFEM 3D Cartesian synthetis
!   to ASKI-conform frequency-domain data files
!----------------------------------------------------------------------------
program transformSpecfemSynthetics
   use argumentParser
   use string
   use constants
   use hdfWrapper
   use inversionBasics
   use inputParameter
   use asciiDataIO
   use asciiSynseisIO
   use seismicEvent
   use seismicEventList
   use seismicNetwork
   use seismicStation
   use discreteFourierTransform
   use axesRotation
   use travelTimes
   use smartUtils
   use errorMessage

   implicit none
   type (argument_parser) :: ap
   type (input_parameter) :: inpar
   type (inversion_basics) :: invbasics
   type (seismic_station) :: station
   type (seismic_event) :: event
   type (discrete_fourier_transform) :: dft
   type (any_rank_real_array) :: arra
   type (error_message) :: errmsg
   type (input_parameter), pointer :: inpar_inv
   type (seismic_event_list), pointer :: evlist
   type (seismic_network) :: statlist
   type (raytable_travel_times) :: ttime_table
   real :: df
   real, dimension(:,:), pointer :: h_p
   real, dimension(:,:), allocatable:: d2
   real, dimension(:,:,:), allocatable :: spsta
   real, dimension(:), allocatable :: usemls
   double precision :: delta,tt,dt,dts
   double precision :: rearth,xlc,ylc,zlc,xg,yg,zg,twinlen,tred,taper_length,tv,tsbuf,ttmin,tend,therd
   double precision :: tanfdp,theta_sta,phi_sta,rd,lon_grid_center,lat_grid_center,thetas,phis,rs,xi
   double precision, dimension(3) :: ul,ug,us,uh
   double precision, dimension(:), allocatable :: usemx,usemy,usemz
   double precision, dimension(:), allocatable :: freq,tapval
   complex, dimension(:), allocatable :: spectrum
   integer (kind=8) :: specid
   integer :: j,nf,nsamp,ic,jh,ntapstart,ntl,ios,nstat,jsta
   integer, dimension(:), pointer :: ifreq
   logical :: apply_taper,do_open
   character(len=char_len_comp) :: comp,cf
   character(len=char_len_sta) :: staname
   character(len=max_length_string) :: main_parfile,path_specfem_seis,file_specfem_trace,main_path_inversion,path_gemini
   character(len=max_length_string) :: path_dft_data,sem_ext,iter_path,iter_parfile,path_kd,specfile,timing_file,ttfile
   character(len=max_length_string) :: pathname,netstaname,evid,path_injection_seis,path_measured_seis,station_file
   character(len=26) :: myname = 'transformSpecfemSynthetics'
   character(len=80), dimension(5) :: iter_inpar_keys
   data iter_inpar_keys/'ITERATION_STEP_NUMBER_OF_FREQ', 'ITERATION_STEP_INDEX_OF_FREQ',&
                        'ASKI_DFT_apply_taper', 'ASKI_DFT_taper_length_kd', 'DT'/
!----------------------------------------------------------------------------------------------------------------
   call init(ap,myname,'Fourier transform of SPECFEM time-domain synthetics into ASKI-conform frequency-domain')
   call addPosarg(ap,'main_parfile','sval','Main parameter file of ASKI inversion')
   call addOption(ap,'-comp',.true.,"Component of seismograms, as available in measurements",'sval','UP')
   call addOption(ap,'-sem',.true.,"SPECFEM synthetics file extension to be used (semd,semv,sema)",'sval','semd')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then; call usage(ap); call print(.errmsg.ap); goto 1; endif
   main_parfile = ap.sval.'main_parfile'
   comp = ap.sval.'-comp'
   if (.not. equalString(comp,'UP') .and. .not. equalString(comp,'N') .and. .not. equalString(comp,'E')) then
      call add(errmsg,2,'Only UP, N and E components currently implemented',myname)
      goto 1
   end if
   sem_ext = ap.sval.'-sem'
   call document(ap)
   call dealloc(ap)
!------------------------------------------------------------------------
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  setup basic stuff from ASKI main parameter file
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path_inversion = inpar_inv.sval.'MAIN_PATH_INVERSION'
   path_injection_seis = main_path_inversion+(inpar_inv.sval.'PATH_INJECTION_SEIS')
   path_gemini = main_path_inversion+(inpar_inv.sval.'PATH_GEMINI')
   path_measured_seis = main_path_inversion+(inpar_inv.sval.'PATH_MEASURED_SEIS')
   lat_grid_center = mc_deg2radd*(inpar_inv.dval.'INVGRID_CENTER_LAT')
   lon_grid_center = mc_deg2radd*(inpar_inv.dval.'INVGRID_CENTER_LON')
   rearth = inpar_inv.dval.'REARTH'
   evlist => getEventListInversionBasics(invbasics)
   ttfile = trim(path_gemini)+(inpar_inv.sval.'ASKI_TTIME_TABLE_FILE_KD')
!
!  setup frequency discretization to be used in Discrete Fourier Transform objects below
!  perfor DFT for all frequencies and not only for those specified in th iter_parfile
!
   df = rval(inpar_inv,'MEASURED_DATA_FREQUENCY_STEP')
   ifreq => getMeasuredDataFrequencyIndicesInversionBasics(invbasics)
   nf = ival(inpar_inv,'MEASURED_DATA_NUMBER_OF_FREQ')
   allocate(freq(nf))
   freq = [(ifreq(j)*df, j=1,nf)]
   allocate(d2(nf,3))
!
!  get info from iter_parfile
!
   write(iter_path,"(2a,i3.3,a)") trim(main_path_inversion),&
         trim(inpar_inv.sval.'ITERATION_STEP_PATH'), inpar_inv.ival.'CURRENT_ITERATION_STEP','/'
   iter_parfile = trim(iter_path)//trim(inpar_inv.sval.'PARFILE_ITERATION_STEP')
   call createKeywordsInputParameter(inpar,iter_inpar_keys)
   call readSubroutineInputParameter(inpar,1,trim(iter_parfile),errmsg)
   if (.level.errmsg == 2) goto 1
!
   apply_taper = inpar.lval.'ASKI_DFT_apply_taper'
   taper_length = inpar.dval.'ASKI_DFT_taper_length_kd'
   dt = inpar.dval.'DT'
   path_dft_data = trim(iter_path)//trim(inpar_inv.sval.'PATH_SYNTHETIC_DATA')
   path_kd = trim(iter_path)//trim(inpar_inv.sval.'PATH_KERNEL_DISPLACEMENTS')
   specfile = trim(path_dft_data)//'syn_'//trim(comp)//'.hdf'
   call dealloc(inpar)
!
!  precalculate taper values
!
   ntl = nint(taper_length/dt)
   call tailTaperSmartUtils(ntl,tapval)
!
!  read travel time table
!
   print *,trim(ttfile)
   call readTableTravelTimes(ttime_table,trim(ttfile),errmsg,parallel=.false.)
   if (.level.errmsg == 2) goto 1
!
!  open HDF output file for writing
!
   call createFileHDFWrapper(specfile,specid,errmsg)
   if (.level.errmsg == 2) goto 1
!--------------------------------------------------------------------------------------
!  Loop over events
!  Original SPECFEM seismograms have same number of samples for given event
!--------------------------------------------------------------------------------------
   do while (nextEventSeismicEventList(evlist,event))
      call printSeismicEvent(event)
      evid = .evid.event
      therd = .tanfdp.(.otime.event)
      path_specfem_seis = trim(path_kd)//trim(evid)//'_OUTPUT_FILES/'
      print *, trim(path_specfem_seis)
      station_file = trim(path_measured_seis)//trim(evid)//'/ASKI_clean_station_file'
      print *,trim(station_file)
      call createStationListFromStationFile(station_file,1,'ASKI_stations',statlist,errmsg)
      if (.level.errmsg == 2) goto 1
   !
   !  read reduction time of SPECFEM synthetics from injection timing files
   !
      timing_file = trim(path_injection_seis)//"timing_"//trim(evid)
      open(1,file = trim(timing_file),form='FORMATTED',status='OLD',action='READ',iostat=ios)
      if (ios /= 0) then
         call add(errmsg,2,'cannot open timing file',myname)
         goto 1
      endif
      read(1,*) twinlen,tsbuf,ttmin,tred,tend
      close(1)
   !
   !  source part of travel time calculation
   !
      thetas = 0.5*mc_pid-(.slatrad.event)
      phis = .slonrad.event
      rs = rearth-(.sdepth.event)
      call setSourceInfoTravelTimes(ttime_table,rs,errmsg)
      if (.level.errmsg == 2) goto 1
   !
   !  Loop over stations
   !
      nstat = .nstat.statlist
      allocate(spsta(nstat,nf,4))
      do_open = .true.
      do while (nextStationSeismicNetwork(statlist,station,c=jsta))
         staname = .staname.station
         netstaname = trim(.netcode.station)//'.'//trim(staname)
      !
      !  calculate delta for station
      !
         theta_sta = 0.5*mc_pid-(.latrad.station)
         phi_sta = .lonrad.station
         rd = rearth                       ! ignore altitude of station here
      !
      !  global Cartesian coordinates of station
         call coordinatesLCfromLSAxesRotation(rd,theta_sta,phi_sta,xg,yg,zg)
      !
      !  local Cartesian coordinates of station with origin at source
         call coordinatesLCfromGCAxesRotation(thetas,phis,xg,yg,zg,xlc,ylc,zlc)
      !
      !  spherical coordinates of station with origin at source
         call coordinatesLSfromLCAxesRotation(xlc,ylc,zlc,rd,delta,xi)
      !
         call getReceiverTravelTimes(ttime_table,rd,delta,tt,errmsg)
         if (.level.errmsg == 2) goto 1
         write(6,'(a13,a2,4f15.4,f12.2)') trim(netstaname),': ',(0.5*mc_pid-theta_sta)/mc_deg2rad,phi_sta/mc_deg2rad,&
                                          delta/mc_deg2rad,xi/mc_deg2rad,tt
      !
      !  go through SPECFEM components X,Y,Z
      !  jump to next station if current one does not exist
      !
         file_specfem_trace = trim(path_specfem_seis)//trim(netstaname)//'.BXX.'//trim(sem_ext)
         call readRealMatrixAsciiDataIO(file_specfem_trace,1,h_p,errmsg,2)
         if (.level.errmsg == 2) then
            write(6,'(a,a,a)') 'Clean station: ',trim(netstaname),' was not included in the station file'
            call updateErrorMessage(errmsg,0,'Missing station error was resolved',myname)
            cycle
         endif
         tanfdp = dble(h_p(1,1)); dts = h_p(2,1)-h_p(1,1); nsamp = size(h_p,1)
         allocate(usemx(size(h_p,1))); usemx = h_p(:,2)
         deallocate(h_p)
      !
         file_specfem_trace = trim(path_specfem_seis)//trim(netstaname)//'.BXY.'//trim(sem_ext)
         call readRealMatrixAsciiDataIO(file_specfem_trace,1,h_p,errmsg,2)
         if (.level.errmsg == 2) return
         allocate(usemy(size(h_p,2))); usemy = h_p(:,2)
         deallocate(h_p)
      !
         file_specfem_trace = trim(path_specfem_seis)//trim(netstaname)//'.BXZ.'//trim(sem_ext)
         call readRealMatrixAsciiDataIO(file_specfem_trace,1,h_p,errmsg,2)
         if (.level.errmsg == 2) return
         allocate(usemz(size(h_p,2))); usemz = h_p(:,2)
         deallocate(h_p)
      !
         if (abs(nsamp-int((tend-tred)/dt)) > 10) then
            call add(errmsg,2,'TEND-TRED and nsamp of seismograms are inconsistent',myname)
            print *,'TEND-TRED: ',int((tend-tred)/dt),' NSAMP = ',nsamp
            goto 1
         end if
      !
         if (abs(dts/dt-1.d0) > 1.e-5) then
            call add(errmsg,2,'DT and dt of seismograms is inconsistent',myname)
            print *,'DT = ',dt,' dts = ',dts
            goto 1
         end if
      !
      !  start of taper is at tt+twinlen-tred-taplen because SPECFEM synthetics begin at t=tred
      !
         if (apply_taper) then
            ntapstart = nint((tt+twinlen-tred)/dt)+1-ntl
         else
            ntapstart = nsamp-ntl+1
         end if
      !
      !  transform SPECFEM XYZ components to UP, N, E components, create time series
      !
         allocate(usemls(nsamp))
         if (equalString(comp,'UP')) then
            ic = 1
         else if (equalString(comp,'N')) then
            ic = 2
         else if (equalString(comp,'E')) then
            ic = 3
         else
            ic = -1
            call add(errmsg,2,'Unknown component',myname)
            goto 1
         endif
      !
      !  compute ZNE components and apply tail taper
      !
         do j = 1,nsamp
            uh(1) = usemx(j); uh(2) = usemy(j); uh(3) = usemz(j)
            call vectorDbleLCfromRCAxesRotation(0.5d0*mc_pid,uh,ul)
            call vectorDbleGCfromLCAxesRotation(0.5*mc_pid-lat_grid_center,lon_grid_center,ul,ug)
            call vectorDbleLSfromLCAxesRotation(theta_sta,phi_sta,ug,us)
            jh = j-ntapstart
            if (jh < 0) then
               tv = 1.d0
            else if (jh >= 0 .and. jh < ntl) then 
               tv = tapval(jh+1)
            else
               tv = 0.d0
            end if
            usemls(j) = -(-1)**ic*us(ic)*tv    ! minus sign for ic = 2
         enddo
         deallocate(usemx,usemy,usemz)
      !
      !  write tail tapered SPECFEM synthetic to gather file
      !  start time of synthetics is therd+tred, ignore non-zero start time in SPECFEM synthetic files
      !
         file_specfem_trace = trim(path_dft_data)//trim(evid)//'_gather_ttap_syn_'//trim(comp)
         call writeAsciiSynseisIO(2,file_specfem_trace,ntapstart+ntl,therd+tred,sngl(dts),netstaname,comp,&
                                  usemls(1:ntapstart+ntl),ios,op=do_open,cl=.false.)
         if (ios > 0) then
            call add(errmsg,2,'can not open '//trim(file_specfem_trace),myname)
            goto 1
         end if
         do_open = .false.
      !
      !  perform DFT
      !
         if (.not. .isdef.dft) call initiateForwardDFT(dft,dble(dt),nsamp,freq,errmsg)
         if (.level.errmsg == 2) goto 1
         call transformForwardDFT(dft,usemls,spectrum,errmsg)
         if (.level.errmsg == 2) goto 1
         deallocate(usemls)
      !
      ! add this spectrum HDF file
      !
         pathname = trim(evid)//"_"//trim(netstaname)
         d2(:,1) = real(freq)
         d2(:,2) = real(spectrum)
         d2(:,3) = aimag(spectrum)
         call arra%assoc2d(d2)
         call writeArrayHDFWrapper(specid,trim(pathname),arra,errmsg)
         if (.level.errmsg == 2) goto 1
         call arra%deassoc()
         deallocate(spectrum)
      !
      !  store real and imaginary part and compute magnitude and phase
      !
         spsta(jsta,:,1) = d2(:,2)
         spsta(jsta,:,2) = d2(:,3)
         spsta(jsta,:,3) = abs(cmplx(spsta(jsta,:,1),spsta(jsta,:,2)))
         spsta(jsta,:,4) = atan2(spsta(jsta,:,2),spsta(jsta,:,1))/mc_deg2rad
      enddo
      close(2)
      call dealloc(dft)
   !
   !  gather values of re, im, mag and phase versus station index
   !
      file_specfem_trace = trim(path_dft_data)//trim(evid)//'_gather_syn_spec_statidx_'//trim(comp)
      do_open = .true.
      do j = 1,nf
         write(cf,'(i2.2)') ifreq(j)
         call writeAsciiSynseisIO(2,file_specfem_trace,nstat,1.d0,1.0,'real',cf,spsta(:,j,1),ios,op=do_open,cl=.false.)
         call writeAsciiSynseisIO(2,file_specfem_trace,nstat,1.d0,1.0,'imag',cf,spsta(:,j,2),ios,op=.false.,cl=.false.)
         call writeAsciiSynseisIO(2,file_specfem_trace,nstat,1.d0,1.0,'abs ',cf,spsta(:,j,3),ios,op=.false.,cl=.false.)
         call writeAsciiSynseisIO(2,file_specfem_trace,nstat,1.d0,1.0,'phas',cf,spsta(:,j,4),ios,op=.false.,cl=.false.)
         do_open = .false.
         if (ios > 0) then
            call add(errmsg,2,'can not open '//trim(file_specfem_trace),myname)
            goto 1
         end if
      enddo
      close(2)
      deallocate(spsta)
      call dealloc(statlist)
   enddo
!------------------------------------------------------------------------
!  clean up
!
   deallocate(d2,tapval)
   call closeFileHDFWrapper(specid,errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(ttime_table)
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(invbasics)
   call dealloc(errmsg)
!------------------------------------------------------------------------
!  treat errors
!
 1 if (.level.errmsg == 2) then
       call print(errmsg)
   endif
!
end program transformSpecfemSynthetics
