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
!> \brief program for transforming time-domain measured data to ASKI-conform frequency-domain measured data files
!
!  Written for data stored in Ascii synseis files
!
!! \author Florian Schumacher, Wolfgang Friederich
!! \date January 2015, December 2021

program transformMeasuredDataAsciiSeis
   use argumentParser
   use axesRotation
   use string
   use constants
   use hdfWrapper
   use inversionBasics
   use asciiSynseisIO
   use seismicEventList
   use seismicEvent
   use seismicNetwork
   use seismicStation
   use discreteFourierTransform
   use travelTimes
   use errorMessage
   use smartUtils
   implicit none
   type (argument_parser) :: ap
   type (inversion_basics) :: invbasics
   type (seismic_event) :: event
   type (seismic_station) :: station
   type (seismic_network) :: station_list
   type (discrete_fourier_transform) :: dft
   type (any_rank_real_array) :: arra
   type (error_message) :: errmsg
   type (input_parameter) :: inpar
   type (input_parameter), pointer :: inpar_inv
   type (seismic_event_list), pointer :: event_list
   type (raytable_travel_times) :: ttime_table
   real :: dt
   double precision :: df,dtnew,taper_length
   double precision :: tanfdp,therd,twinlen,tsbuf,ttmin,tred,tend,tanew
   double precision :: rearth,xlc,ylc,zlc,xg,yg,zg
   double precision :: theta_sta,phi_sta,rd,thetas,phis,rs,xi,vsed,alt,mean_alt
   double precision :: delta,tt,slow,tsalt
   real, dimension(:), allocatable :: urs,unew,d1
   real, dimension(:,:), allocatable :: d2
   real, dimension(:,:,:), allocatable :: spsta
   double precision, dimension(:), allocatable :: freq,tapval
   complex, dimension(:), allocatable :: spectrum
   integer (kind=8) :: specid,toutid
   integer :: j,nf,nsamp,ios,nsnew,nzmax,ntl,nstat,js
   integer, dimension(:), pointer :: ifreq
   logical :: do_open,write_timing
   character(len=char_len_comp) :: comp,h_comp,cf
   character(len=char_len_sta) :: staname
   character(len=max_length_string) :: main_parfile,path_measured_traces,specfile,evid,netstaname,pathname,ttfile,toutfile
   character(len=max_length_string) :: station_file,path_measured_seis,path_dft_data,file_measured_trace,path_gemini
   character(len=max_length_string) :: timing_file,path_injection_seis,main_path_inversion,iter_parfile,iter_path
   character(len=30) :: myname = 'transformMeasuredDataAsciiSeis'
   character(len=80), dimension(4) :: iter_inpar_keys
   data iter_inpar_keys/'ASKI_DFT_taper_length_kd', 'ASKI_DFT_taper_length_gt', 'DT', 'VP_SEDIMENT'/
!----------------------------------------------------------------------------------------------------------------
   call init(ap,myname,'Fourier transform of measured seismograms to ASKI-conform frequency-domain measured data')
   call addPosarg(ap,'main_parfile','sval','Main parameter file of ASKI inversion')
   call addOption(ap,'-comp',.true.,'Component of measured seismograms (UP,N,E)','sval','UP')
   call addOption(ap,'-nzmax',.true.,'number of zeros of sinc where cutoff is made','ival','5')
   call addOption(ap,'-timing',.false.,'write timing of data to HDF file')
   call parse(ap)
   call document(ap)
   main_parfile = ap.sval.'main_parfile'
   comp = ap.sval.'-comp'
   nzmax = ap.ival.'-nzmax'
   write_timing = ap.optset.'-timing'
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
   event_list => getEventListInversionBasics(invbasics)
   main_path_inversion = inpar_inv.sval.'MAIN_PATH_INVERSION'
   path_measured_seis = main_path_inversion+(inpar_inv.sval.'PATH_MEASURED_SEIS')
   path_injection_seis = main_path_inversion+(inpar_inv.sval.'PATH_INJECTION_SEIS')
   path_gemini = main_path_inversion+(inpar_inv.sval.'PATH_GEMINI')
   rearth = inpar_inv.dval.'REARTH'
   ttfile = trim(path_gemini)+(inpar_inv.sval.'ASKI_TTIME_TABLE_FILE_KD')
!
!  setup frequency discretization to be used in Discrete Fourier Transform objects below
!
   df = dval(inpar_inv,'MEASURED_DATA_FREQUENCY_STEP')
   ifreq => getMeasuredDataFrequencyIndicesInversionBasics(invbasics)
   nf = ival(inpar_inv,'MEASURED_DATA_NUMBER_OF_FREQ')
   allocate(freq(nf))
   freq = (/ (ifreq(j)*df, j=1,nf) /)
   path_dft_data = main_path_inversion+sval(inpar_inv,'PATH_MEASURED_DATA')
   specfile = trim(path_dft_data)//'data_'//trim(comp)//'.hdf'
   toutfile = trim(path_dft_data)//'timing_of_data.hdf'
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
   dtnew = inpar.dval.'DT'
   taper_length = inpar.dval.'ASKI_DFT_taper_length_kd'
   vsed = inpar.dval.'VP_SEDIMENT'
!
!  read travel time table
!
   print *,trim(ttfile)
   call readTableTravelTimes(ttime_table,trim(ttfile),errmsg,parallel=.false.)
   if (.level.errmsg == 2) goto 1
!
!  open HDF output file for writing spectral values
!
   call createFileHDFWrapper(specfile,specid,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  open HDF output file for writing time information
!
   if (write_timing) then
      call createFileHDFWrapper(toutfile,toutid,errmsg)
      if (.level.errmsg == 2) goto 1
   end if
!
!  precalculate front taper values
!
   ntl = nint(taper_length/dtnew)
   call frontTaperSmartUtils(ntl,tapval)
!--------------------------------------------------------------------------------------
!  Loop over events
!--------------------------------------------------------------------------------------
   do while (nextEventSeismicEventList(event_list,event))
      call printSeismicEvent(event)
      evid = .evid.event
      therd = .tanfdp.(.otime.event)
      path_measured_traces = trim(path_measured_seis)//trim(evid)//'/per_trace/'
      station_file = trim(path_measured_seis)//trim(evid)//'/ASKI_clean_station_file'
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
   !  read ASKI_clean_station file, create station list, loop over stations
   !
      call createStationListFromStationFile(station_file,1,'ASKI_stations',station_list,errmsg)
      if (.level.errmsg == 2) goto 1
      nstat = .nstat.station_list
   !
   !  mean of station altitudes
   !
      mean_alt = 0.d0
      do while (nextStationSeismicNetwork(station_list,station,c=js))
         mean_alt = mean_alt+.alt.station
      end do
      mean_alt = mean_alt/nstat
      print *,'Mean station altidue is: ',mean_alt
   !
      allocate(spsta(nstat,nf,4))
      do_open = .true.
      do while (nextStationSeismicNetwork(station_list,station,c=js))
         netstaname = trim(.netcode.station)//'.'//trim(.staname.station)
         file_measured_trace = trim(path_measured_traces)//'data_'//trim(netstaname)//'_'//trim(comp)
         print *,trim(file_measured_trace)
         call readAsciiSynseisIO(1,file_measured_trace,nsamp,tanfdp,dt,staname,h_comp,urs,ios)
         if (ios > 0) then
            call add(errmsg,2,'can not write to '//trim(file_measured_trace),myname)
            goto 1
         end if
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
         call getReceiverTravelTimes(ttime_table,rd,delta,tt,errmsg,slow=slow)
         if (.level.errmsg == 2) goto 1
      !
      !  time correction of data due to station altitude
      !  reference altitude is the average of the station altitudes taken into account by source wavelet fitting
      !  if station is higher, P wave arrives later and must be corrected to the ref altitude by a negative time shift
      !  if station is lower, P wave arrives earlier and must be corrected to the ref altitude by a positive time shift
      !
         alt = .alt.station
         tsalt = (mean_alt-alt)*sqrt(1.d0/vsed**2-(slow/rearth)**2)
         write(6,'(a,f15.5,3f15.2)') 'Time shift for altitude correction: ',tsalt,alt,delta/mc_deg2rad,slow
      !
      !  upsample data to SPECFEM's DT before Fourier transforming
      !  start at the SPECFEM sample closest to and before tanfdp
      !  SPECFEM synthetics start at therd+tred
      !  start of data is later than this because of preprocessing in computeSourceWavelet
      !
         tanew = floor((tanfdp-therd-tred)/dtnew)*dtnew+therd+tred
         call sincInterpolateData(nsamp,tanfdp,dble(dt),urs,tanew,dtnew,nzmax,nsnew,unew)
      !
      !  apply a front taper to interpolated data before DFT
      !
         unew(1:ntl) = unew(1:ntl)*tapval
      !
      !  add interpolated and tapered data to gather file
      !
         file_measured_trace = trim(path_dft_data)//trim(evid)//'_gather_sinc_ftap_data_'//trim(comp)
         call writeAsciiSynseisIO(2,file_measured_trace,nsnew,tanew,sngl(dtnew),netstaname,comp,unew,ios,op=do_open,cl=.false.)
         if (ios > 0) then
            call add(errmsg,2,'can not open '//trim(file_measured_trace),myname)
            goto 1
         end if
         do_open = .false.
      !
      !  correct phase of Fourier transform by difference of data start time to therd+tred,
      !  time shift theorem: f(t-s) => F(om)*exp(-i*om*s)
      !  the time therd + tred serves as time reference for Specfem synthetics, kernel seismograms and data
      !  i.e. zero time in the Fourier transforms corrsponds to therd+tred
      !  no taper here because data are already tapered when written by computeSourceWavelet
      !
         call initiateForwardDFT(dft,dtnew,nsnew,freq,errmsg,tshift = tanew-(therd+tred)+tsalt)
         if (.level.errmsg == 2) goto 1
         call transformForwardDFT(dft,unew,spectrum,errmsg)
         if (.level.errmsg == 2) goto 1
         deallocate(urs,unew)
      !
      ! add this spectrum to HDF file
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
         call dealloc(dft)
      !
      ! add timing to HDF file
      !
         if (write_timing) then
            d1 = [tred,twinlen,tt,tt+twinlen-tred]
            call arra%assoc1d(d1)
            call writeArrayHDFWrapper(toutid,trim(pathname),arra,errmsg)
            if (.level.errmsg == 2) goto 1
            call arra%deassoc()
            deallocate(d1)
         end if
      !
      !  store real and imaginary part and compute magnitude and phase
      !
         spsta(js,:,1) = d2(:,2)
         spsta(js,:,2) = d2(:,3)
         spsta(js,:,3) = abs(cmplx(spsta(js,:,1),spsta(js,:,2)))
         spsta(js,:,4) = atan2(spsta(js,:,2),spsta(js,:,1))/mc_deg2rad
      enddo
      close(2)
   !
   !  gather values of re, im, mag and phase versus station index
   !
      file_measured_trace = trim(path_dft_data)//trim(evid)//'_gather_data_spec_statidx_'//trim(comp)
      do_open = .true.
      do j = 1,nf
         write(cf,'(i2.2)') ifreq(j)
         call writeAsciiSynseisIO(2,file_measured_trace,nstat,1.d0,1.0,'real',cf,spsta(:,j,1),ios,op=do_open,cl=.false.)
         call writeAsciiSynseisIO(2,file_measured_trace,nstat,1.d0,1.0,'imag',cf,spsta(:,j,2),ios,op=.false.,cl=.false.)
         call writeAsciiSynseisIO(2,file_measured_trace,nstat,1.d0,1.0,'abs ',cf,spsta(:,j,3),ios,op=.false.,cl=.false.)
         call writeAsciiSynseisIO(2,file_measured_trace,nstat,1.d0,1.0,'phas',cf,spsta(:,j,4),ios,op=.false.,cl=.false.)
         do_open = .false.
         if (ios > 0) then
            call add(errmsg,2,'can not open '//trim(file_measured_trace),myname)
            goto 1
         end if
      enddo
      close(2)
      deallocate(spsta)
      call dealloc(station_list)
   enddo
!
   call closeFileHDFWrapper(specid,errmsg)
   if (.level.errmsg == 2) goto 1
   if (write_timing) then
      call closeFileHDFWrapper(toutid,errmsg)
   end if
   if (.level.errmsg == 2) goto 1
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  clean up
!
   deallocate(tapval)
   call dealloc(invbasics)
   call dealloc(errmsg)
!------------------------------------------------------------------------
!  treat errors
!
 1 if (.level.errmsg == 2) then
        call print(errmsg)
   endif
!
   contains
end program transformMeasuredDataAsciiSeis

