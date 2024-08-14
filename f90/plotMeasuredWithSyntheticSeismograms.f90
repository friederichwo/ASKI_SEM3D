!----------------------------------------------------------------------------
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
!   Plot measured and synthetic seismograms for given event and available stations
!   into one record section. Synthetics may be SPECFEM seismograms.
!
program plotMeasuredWithSyntheticSeismograms
   use seismicNetwork
   use readEventStationFile
   use asciiSynseisIO
   use timeSeries
   use axesRotation
   use argumentParser
   use basePlotGather
   use interactiveBindingsPlotGather
   use pgPlotWindow
   use inputParameter
   use hdfWrapper
   use string
   use errorMessage
   use mathConstants
   implicit none

   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (input_parameter) :: inpar
   type (seismic_network) :: station_list
   type (seismic_station) :: station
   type (time_series) :: meas_ts,sem_ts
   type (date_time) :: otime
   type (any_rank_real_array) :: arra
   type (base_plot_gather), dimension(2) :: gather
   type (pgplot_window) :: pgwin,pgps
   integer :: ios,i,ic,shift,nstat,nsamp_spec,nsamp_meas,slen
   integer (kind=8) :: fid,grid
   real, dimension(:), pointer :: d
   real, dimension(:,:), pointer :: h_p
   real, dimension(:), allocatable :: h_u,usemls
   real :: aspect,width,dt_spec,dt_meas
   double precision, dimension(3) :: ul,ug,us,uh
   double precision, dimension(:), allocatable :: usemx,usemy,usemz
   double precision :: latsta,lonsta,lat_grid_center,lon_grid_center,tanfdp_spec,tanfdp_meas,tred
   logical :: hardcopy,leave,gemflag
   character(len=max_length_string) :: main_parfile,iter_path,station_file,sem_ext,val,dtstr
   character(len=max_length_string) :: path_measured_traces,path_specfem_seis,path_injection_seis,main_path_inversion
   character(len=max_length_string) :: file_measured_trace,file_specfem_trace,file_injection_seis
   character(len=char_len_evid) :: evid
   character(len=char_len_comp) :: comp,h_comp
   character(len=char_len_sta+char_len_netcode) :: netstaname,h_staname
   character (len=1) :: otherkey
   character(len=36) :: myname = 'plotMeasuredWithSyntheticSeismograms'
   character(len=80), dimension(9) :: main_inpar_keys
   data main_inpar_keys/'ITERATION_STEP_PATH','MAIN_PATH_INVERSION','PARFILE_ITERATION_STEP','CURRENT_ITERATION_STEP',&
                        'PATH_MEASURED_SEIS','INVGRID_CENTER_LAT','INVGRID_CENTER_LON','PATH_INJECTION_SEIS',&
                        'PATH_KERNEL_DISPLACEMENTS'/
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
!  read commandline arguments
!
   call init(ap,myname,"Plot measured together with Specfem synthetic seismograms")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addPosarg(ap,'evid','sval',"Event Id of seismograms to be plotted. This option must be set.")
   call addOption(ap,'-comp',.true.,"Component of seismograms, as available in measurements",'sval','UP')
   call addOption(ap,'-sem',.true.,"SPECFEM synthetics file extension to be used (semd,semv,sema)",'sval','semd')
   call addOption(ap,'-gemini',.false.,"compare with convolved Gemini synthetics instead")
   call addOption(ap,'-w',.true.,'plot window width','rval','10.')
   call addOption(ap,'-a',.true.,'plot window aspect ratio','rval','0.65')
   call parse(ap)
   call document(ap)
   main_parfile = ap.sval.'main_parfile'
   evid = ap.sval.'evid'
   comp = ap.sval.'-comp'
   if (.not. equalString(comp,'UP') .and. .not. equalString(comp,'N') .and. .not. equalString(comp,'E')) then
      call add(errmsg,2,'Only UP, N and E components currently implemented',myname)
      goto 1
   end if
   sem_ext = ap.sval.'-sem'
   gemflag = ap.optset.'-gemini'
   width = ap.rval.'-w'
   aspect = ap.rval.'-a'
   call dealloc(ap)
!------------------------------------------------------------------------
!  read input parameters from main parameter file
!
   call createKeywordsInputParameter(inpar,main_inpar_keys)
   call readSubroutineInputParameter(inpar,1,trim(main_parfile),errmsg)
   if (.level.errmsg == 2) goto 1
   main_path_inversion = inpar.sval.'MAIN_PATH_INVERSION'
   write(iter_path,"(2a,i3.3,a)") trim(main_path_inversion),&
         trim(inpar.sval.'ITERATION_STEP_PATH'), inpar.ival.'CURRENT_ITERATION_STEP','/'
   path_measured_traces = trim(inpar.sval.'PATH_MEASURED_SEIS')//trim(evid)//'/per_trace/'
   path_measured_traces = main_path_inversion+path_measured_traces
   path_injection_seis = trim(inpar.sval.'PATH_INJECTION_SEIS')
   path_specfem_seis = trim(iter_path)//trim(inpar.sval.'PATH_KERNEL_DISPLACEMENTS')//trim(evid)//'_OUTPUT_FILES/'
   station_file = trim(inpar.sval.'PATH_MEASURED_SEIS')//trim(evid)//'/ASKI_clean_station_file'
   station_file = main_path_inversion+station_file
   lat_grid_center = inpar.dval.'INVGRID_CENTER_LAT'
   lon_grid_center = inpar.dval.'INVGRID_CENTER_LON'
   call dealloc(inpar)
   print *,'Center of Specfem Box: ',lat_grid_center,lon_grid_center
!
!  read ASKI_clean_station file, create station list
!
   call createStationListFromStationFile(station_file,1,'ASKI_stations',station_list,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  read reduction time and origin time of SPECFEM synthetics from injection seismogram
!
   file_injection_seis = trim(path_injection_seis)//trim(evid)//'_seis.hdf'
   print *,'Get timing info from injection seismogram file for this event from: ',trim(file_injection_seis)
   call openFileRoHDFWrapper(file_injection_seis,fid,errmsg)
   if (.level.errmsg == 2) goto 1
   call readArrayAttributeHDFWrapper(fid,'timing',arra,errmsg)
   if (.level.errmsg == 2) goto 1
   d => arra%get1d()
   print *,'Timing of injection seismograms: ',d
   tred = d(4)
   deallocate(d)
   call openGroupHdfWrapper(fid,'Event',grid,errmsg)
   call readStringAttributeHDFWrapper(grid,'Date_and_time',val,slen,errmsg)
   if (.level.errmsg == 2) goto 1
   dtstr = val(1:slen)
   print *,'Event origin time: ',trim(dtstr)
   call new(otime,dtstr)
   call closeGroupHdfWrapper(grid,errmsg)
   if (.level.errmsg == 2) goto 1
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  iterate over stations in list
!
   nstat = .nstat.station_list
   call gather(1)%allocate(nstat)
   call gather(2)%allocate(nstat)
   do while (nextStationSeismicNetwork(station_list,station))
      netstaname = trim(.netcode.station)//'.'//trim(.staname.station)
      latsta = .latrad.station; lonsta = .lonrad.station
   !
   !  read synthetic SPECFEM traces (all three components)
   !  check before if file exists
   !
      file_specfem_trace = trim(path_specfem_seis)//trim(netstaname)//'.BXX.'//trim(sem_ext)
      inquire(file = trim(file_specfem_trace),exist = leave)
      if (.not. leave) cycle
      call readRealMatrixAsciiDataIO(file_specfem_trace,1,h_p,errmsg,2)
      if (.level.errmsg == 2) goto 1
      tanfdp_spec = dble(h_p(1,1)); dt_spec = h_p(2,1)-h_p(1,1); nsamp_spec = size(h_p,1)
      allocate(usemx(size(h_p,1))); usemx = h_p(:,2)
      deallocate(h_p)
      print *,trim(file_specfem_trace)
      print *,'tanf = ',tanfdp_spec,' dt = ',dt_spec,' nsamp = ',nsamp_spec
      print *,'Station location: ',latsta/mc_deg2radd,lonsta/mc_deg2radd
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
   !  transform SPECFEM XYZ components to UP, N, E components, create time series and add to plot gather
   !  shift start time of synthetics by tred+otime
   !
      allocate(usemls(nsamp_spec))
      if (equalString(comp,'UP')) ic = 1
      if (equalString(comp,'N')) ic = 2
      if (equalString(comp,'E')) ic = 3
      do i = 1,nsamp_spec
         uh(1) = usemx(i); uh(2) = usemy(i); uh(3) = usemz(i)
         call vectorDbleLCfromRCAxesRotation(0.5d0*mc_pid,uh,ul)
         call vectorDbleGCfromLCAxesRotation((90.d0-lat_grid_center)*mc_deg2radd,lon_grid_center*mc_deg2radd,ul,ug)
         call vectorDbleLSfromLCAxesRotation(0.5d0*mc_pid-latsta,lonsta,ug,us)
         usemls(i) = -(-1)**ic*us(ic)    ! minus sign for ic = 2
      enddo
      deallocate(usemx,usemy,usemz)
   !
   !  define start time of Specfem synthetic at source time + reduction time
   !  ignore small negative start time in Specfem synthetics files
   !  this leads to a perfect fit with the Gemini synthetics
   !
      call createDPFromDataTimeSeries(sem_ts,nsamp_spec,tred+(.tanfdp.otime),dble(dt_spec),usemls)
      call gather(2)%addTraceAndTag(sem_ts,trim(.staname.station)//' '//trim(comp),errmsg)
      if (.level.errmsg == 2) goto 1
      deallocate(usemls)
   !
   !  read measured trace, create time series and add to plot gather
   !
      if (gemflag) then
         file_measured_trace = trim(path_measured_traces)//'conv_'//trim(netstaname)//'_'//trim(comp)
      else
         file_measured_trace = trim(path_measured_traces)//'data_'//trim(netstaname)//'_'//trim(comp)
      endif
      call readAsciiSynseisIO(1,file_measured_trace,nsamp_meas,tanfdp_meas,dt_meas,h_staname,h_comp,h_u,ios)
      if (ios > 0) then
         call add(errmsg,2,'can not open '//trim(file_measured_trace),myname)
         goto 1
      end if
      if (.not. equalString(h_staname,netstaname)) then
         call add(errmsg,2,'Name of station in file '//trim(h_staname)//' inconsistent with file name',myname)
         goto 1
      end if
      if (.not. equalString(h_comp,comp)) then
         call add(errmsg,2,'Component of trace in file '//trim(h_staname)//' inconsistent with file name',myname)
         goto 1
      end if
      print *,trim(file_measured_trace)
      print *,'tanf = ',tanfdp_meas,' dt = ',dt_meas,' nsamp = ',nsamp_meas
      call createDPFromDataTimeSeries(meas_ts,nsamp_meas,tanfdp_meas,dble(dt_meas),h_u)
      call gather(1)%addTraceAndTag(meas_ts,trim(.staname.station)//' '//trim(comp),errmsg)
      if (.level.errmsg == 2) goto 1
      deallocate(h_u)
   enddo
!
!  setup of plot gathers
!
   call gather(1)%setOffsetTraceIndices()
   call gather(1)%setup(ci = 1)
!
   call gather(2)%setOffsetTraceIndices()
   call gather(2)%setup(ci = 2)
!
!  open plot window and display gather
!
   call createPgPlotWindow(pgwin,width = width, aspect = aspect)
!
!  display seismogram gathers and enable user interaction
!
   call docuInteractiveBindingsPlotGather()
   shift = 0; leave = .false.
   do while (.not. leave)
      call gather(1)%display()
      call gather(2)%overlay(shift)
      call interactWithPlotGather(gather,leave,shift,hardcopy,otherkey)
      if (leave) exit
   !
   !  produce postscript output
   !
      if (hardcopy) then
         call new(pgps,plotfile = 'plotgather.ps')
         call gather(1)%display()
         call gather(2)%overlay(shift)
         call dealloc(pgps); call select_win(pgwin)
      endif
   enddo
   call closeEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
!
!  treat exceptions
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif
end program plotMeasuredWithSyntheticSeismograms
