! ===============================================================================
!  Basis module for plotting seismogram gathers
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!  Defines a base type "basePlotGather" which specifies common interfaces
!  for procedures required to visualize a seismogram section. Type extension
!  must be used to plot seismograms coming from various sources. These
!  must supply a routine to create a gather, i.e. reading data and giving
!  values to the plot_gather type components.
!-----------------------------------------------------------------------
module basePlotGather
    use string
    use timeSeries
    use fourierSpectrum
    use errorMessage
    implicit none
    type :: base_plot_gather
        integer :: nseis                                              ! true number of traces in gather
        type (time_series), dimension(:), allocatable :: seis         ! array of time series
        type (fourier_spectrum), dimension(:), allocatable :: spec    ! array of Fourier spectra
        real, dimension(:), allocatable :: offset                     ! y-offsets of seismograms
        character (len=12), dimension(:), allocatable :: tag          ! textual description of seismogram
        real, dimension(:), pointer :: smax                           ! current normalization of seismograms
        real, dimension(:), pointer :: somax                          ! original normallization of seismograms
        real :: xmin,xmax,ymin,ymax                                   ! current limits of plot area
        real :: xomin,xomax,yomin,yomax                               ! original limits of plot area
        real :: h,oh,ch                                               ! height for one trace, original height, character height
        real :: tref,pred                                             ! zero time, reduction slowness
        integer :: ls,lw,ci                                           ! line style, line width, color index
        character (len=max_length_string) :: xlab,ylab,title          ! labels of x and y axes, plot title
        logical :: specflag                                           ! plot spectra instead of time series
        character (len=1) :: specmode                                 ! spectra plottig mode: amp or phase
        logical :: plottag                                            ! plot station tags
        logical :: txflag                                             ! plot in time-offset mode
        logical :: pickflag                                           ! plot also picks if available
        logical :: time_after_trace_start                             ! zero of time axis is at trace begin / else midnight
        real, dimension(:,:), allocatable :: tpick                    ! times for pick, earliest and latest pick of each trace
        character (len=max_length_string) :: pickfile                 ! file where picks will be saved to
!
    contains
        procedure :: allocate => allocateBasePlotGather                                    ! allocation of arrays
        procedure :: addTraceAndTag => addTraceAndTagBasePlotGather                        ! aadd traces with tags
        procedure :: computeSpectra => computeSpectraBasePlotGather                        ! calculate spectra from time series
        procedure :: display => displayBasePlotGather                                      ! display plot section
        procedure :: overlay => overlayBasePlotGather                                      ! overlay graph onto existing axes
        procedure :: initXLimits => initXLimitsBasePlotGather                              ! initial limits for x-axis
        procedure :: initYLimits => initYLimitsBasePlotGather                              ! initial limits for y-axis
        procedure :: initNormalization => initNormalizationBasePlotGather                  ! initialize normalization factors for traces
        procedure :: writePicks => writePicksBasePlotGather                                ! write picks to pickfile
        procedure :: dealloc => deallocBasePlotGather                                       ! deallocate
        procedure :: setup => setupBasePlotGather                                          ! basic setup of any plot section object, all type components below tag
        procedure :: plot => plotBasePlotGather                                            ! actually plot time series, called by display and overlay
        procedure :: plotSpectra => plotSpectraBasePlotGather                               ! plot spectra
        procedure :: invokeAmpSpectra => invokeAmpSpectraBasePlotGather                    ! switch to amplitude spectra plotting
        procedure :: invokePhaseSpectra => invokePhaseSpectraBasePlotGather                ! switch to phase spectra plotting
        procedure :: invokeTimeSeries => invokeTimeSeriesBasePlotGather                    ! switch to time series plotting
        procedure :: setTref => setTrefBasePlotGather                                      ! set zero time of seismograms
        procedure :: resetTref => resetTrefBasePlotGather                                  ! set back zero time of seismograms to zero
        procedure :: setReductionSlowness => setReductionSlownessBasePlotGather            ! set a slowness for a velocity reduced plot
        procedure :: activateTag => activateTagBasePlotGather                              ! seismogram tag will be written
        procedure :: deactivateTag => deactivateTagBasePlotGather                          ! seismogram tag will not be written
        procedure :: setXYLimits => setXYLimitsBasePlotGather                              ! modify current axis limits
        procedure :: resetXYLimits => resetXYLimitsBasePlotGather                          ! go back to original axis limits
        procedure :: setOffsetTraceIndices => setOffsetTraceIndicesBasePlotGather          ! initialize or reset offsets to trace indices
        procedure :: setScale => setScaleBasePlotGather                                    ! set amplitude scaling factor
        procedure :: setNormalization => setNormalizationBasePlotGather                    ! set trace specific normalization factors
        procedure :: setGlobalNormalization => setGlobalNormalizationBasePlotGather        ! set normalization factor applied to all traces
        procedure :: resetNormalization => resetNormalizationBasePlotGather                ! reset normalization factors to initial values
        procedure :: setTitle => setTitleBasePlotGather                                    ! set title of plot
        procedure :: getXmax => getXmaxBasePlotGather                                      ! get max value of x-axis
        procedure :: getXmin => getXminBasePlotGather                                      ! get min value of x-axis
        procedure :: getYmax => getYmaxBasePlotGather                                      ! get max value of y-axis
        procedure :: getYmin => getYminBasePlotGather                                      ! get min value of y-axis
        procedure :: getAxesLimits => getAxesLimitsBasePlotGather                          ! get min value of y-axis
        procedure :: getNormalization => getNormalizationBasePlotGather                    ! get array of current normalization factors
        procedure :: getOriginalNormalization => getOriginalNormalizationBasePlotGather    ! get array of original normalization factors
        procedure :: addPick => addPickBasePlotGather                                      ! add a pick to the tpick array
    end type base_plot_gather
!
    contains
!---------------------------------------------------------------------------
!  Pre-allocate arrays of type, add time series later
!
  subroutine allocateBasePlotGather(this,nseis)
     class (base_plot_gather) :: this
     integer :: nseis
     this%nseis = 0
     allocate(this%seis(nseis),this%tag(nseis),this%offset(nseis),this%smax(nseis),this%somax(nseis))
   end subroutine allocateBasePlotGather
!---------------------------------------------------------------------------
!  add a time series to gather
!
  subroutine addTraceAndTagBasePlotGather(this,ts,tag,errmsg)
     class (base_plot_gather) :: this
     type (time_series) :: ts
     character (len=*) :: tag
     type (error_message) :: errmsg
     if (size(this%seis) > this%nseis) then
        this%nseis = this%nseis+1
        this%seis(this%nseis) = ts
        this%tag(this%nseis) = tag
     else
        call add(errmsg,2,'no preallocated space for traces available','addTraceAndTagBasePlotGather')
        return
     endif
  end subroutine addTraceAndTagBasePlotGather
!---------------------------------------------------------------------------
!  Basic common setup for plot_section object
!
  subroutine setupBasePlotGather(this,lw,ci,ls,ch,xlab,ylab,title,txmode,pickfile,tzero)
     class (base_plot_gather) :: this
     integer, optional :: ls,lw,ci
     character (len=*), optional :: xlab,ylab,title,pickfile
     real, optional :: ch
     logical, optional :: txmode,tzero
     integer :: i,j,id
     logical :: exflag
     real :: dmin,dmax
   !
     this%tref = 0.0; this%pred = 0.0; this%specflag = .false.
     this%specmode = 'a'; this%plottag = .true.
   !
   !  deal with pick file
   !
     if (present(tzero)) then; this%time_after_trace_start = tzero; else; this%time_after_trace_start = .false.; endif
     if (present(txmode)) then; this%txflag = txmode; else; this%txflag = .false.; endif
     if (present(pickfile)) then
        this%pickfile = pickfile
     else
        this%pickflag = .false.
        this%pickfile = 'none'
     endif
     if (.not.(this%pickfile.equal.'none')) then
        this%pickflag = .true.
        allocate(this%tpick(this%nseis,3))
        inquire(file = pickfile, exist = exflag)
        if (exflag) then
           open(1,file=pickfile)
           do i = 1,this%nseis
              read(1,*) id,(this%tpick(id,j),j=1,3)
           enddo
           close(1)
        else
           this%tpick = 0.0
        endif
     endif
   !
   !  width of individual trace panels
   !
     dmax = maxval(this%offset); dmin = minval(this%offset)
     if(this%nseis == 1) then; this%h = 1.; else; this%h = (dmax-dmin)/(this%nseis-1); endif
     this%oh = this%h
   !
   !  normalization and axis limits
   !
     call this%initNormalization(); call this%initXLimits(); call this%initYLimits()
   !
   !  line width, line style and line color and character height
   !      
     if (present(lw)) then; this%lw = lw; else; this%lw=2; endif
     if (present(ci)) then; this%ci = ci; else; this%ci=1; endif
     if (present(ls)) then; this%ls = ls; else; this%ls=1; endif
     if (present(ch)) then; this%ch = ch; else; this%ch=1.0; endif
   !
   !  axis labels and title
   !
     if (this%txflag) then
        if (present(xlab)) then; this%xlab = xlab; else; this%xlab = 'Trace number'; endif
        if (present(ylab)) then; this%ylab = ylab; else; this%ylab = 'Time'; endif
     else
        if (present(xlab)) then; this%xlab = xlab; else; this%xlab = 'Time'; endif
        if (present(ylab)) then; this%ylab = ylab; else; this%ylab = 'Trace number'; endif
     endif
     if (present(title)) then; this%title = title; else; this%title=' '; endif
  end subroutine setupBasePlotGather
!--------------------------------------------------------------------
!  deallocate arrays
!
  subroutine deallocBasePlotGather(this)
     class (base_plot_gather) :: this
     integer :: j
     if (allocated(this%tag)) deallocate(this%tag)
     if (allocated(this%offset)) deallocate(this%offset)
     if (associated(this%smax)) deallocate(this%smax)
     if (associated(this%somax)) deallocate(this%somax)
     if (allocated(this%seis)) then
        do j=1,this%nseis; call dealloc(this%seis(j)); enddo
        deallocate(this%seis)
     endif
     if (allocated(this%spec)) then
        do j=1,this%nseis; call dealloc(this%spec(j)); enddo
        deallocate(this%spec)
     endif
     if (allocated(this%tpick)) deallocate(this%tpick)
  end subroutine deallocBasePlotGather
!----------------------------------------------------------------
!  display seismogram section
!
  subroutine displayBasePlotGather(this)
     class (base_plot_gather) :: this
   !
     call pgslw(1); call pgsci(1); call pgscf(2); call pgsch(1.25)
     call pgask(.false.)
     call pgeras
     call pgsvp(0.1,0.9,0.1,0.9)
     call pgswin(this%xmin,this%xmax,this%ymin,this%ymax)
     call pgslw(3); call pgsch(this%ch)
     if (this%specflag) then
        call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
     else
        if (this%txflag) then
           call pgbox('BCNST',0.0,0,'BCNST',0.0,0)
        else
           call pgtbox('ZXOHBCNST',0.0,0,'BCNST',0.0,0)
!           call pgtbox('BCNST',0.0,0,'BCNST',0.0,0)
        endif
     endif
     call pglab(this%xlab,this%ylab,this%title)
     if (this%specflag) then
        call this%plotSpectra(0)
     else
        call this%plot(0)
     endif
  end subroutine displayBasePlotGather
!------------------------------------------------------------------
!  overlay a section on a section plot
!
  subroutine overlayBasePlotGather(this,shift)
     class (base_plot_gather) :: this
     integer :: shift
   !
     if (this%specflag) then
        call this%plotSpectra(shift)
     else
        call this%plot(shift)
     endif
  end subroutine overlayBasePlotGather
!----------------------------------------------------------------------
!  do plotting of the time series (for internal use)
!
  subroutine plotBasePlotGather(this,shift)
     class (base_plot_gather) :: this
     integer :: i,j,shift,nsamp
     real, dimension(:), allocatable :: x
     real, dimension(:), pointer :: pdata
     real :: off,tanf,dt
   !
     nullify(pdata)
     call pgsch(this%ch); call pgslw(this%lw); call pgsci(this%ci)
     do i = 1,this%nseis
        off = 0.2*this%h*shift
        if (this%txflag) then
           if(this%offset(i) > this%xmax .or. this%offset(i) < this%xmin) cycle
        else
           if(this%offset(i) > this%ymax .or. this%offset(i) < this%ymin) cycle
        endif
        nsamp = .nsamp.(this%seis(i)); tanf = .tanf.(this%seis(i)); dt = .dt.(this%seis(i))
        allocate(x(nsamp))
        if (this%time_after_trace_start) then
           do j=1,nsamp
              x(j) = -this%tref-this%pred*this%offset(i)+(j-1)*dt
           enddo
        else
           do j=1,nsamp
              x(j) = tanf-this%tref-this%pred*this%offset(i)+(j-1)*dt
           enddo
        endif
        pdata => .trace.this%seis(i)
        if (this%txflag) then
           call pgline(nsamp,this%offset(i)+off-pdata(1:nsamp)/this%smax(i)*0.5*this%h,x)
           if(this%plottag) call pgptxt(this%offset(i)+off,1.02*this%ymax-.02*this%ymin,90.0,0.0,trim(this%tag(i)))
           if (this%pickflag) then
              call pgsch(1.5)
              do j = 1,3
                 call pgsci(1+j); call pgpt1(this%offset(i)+off,this%tpick(i,j),-4)
              enddo
              call pgsci(this%ci); call pgsch(this%ch)
           endif
        else
           call pgline(nsamp,x,this%offset(i)+off+pdata(1:nsamp)/this%smax(i)*0.5*this%h)
           if(this%plottag) call pgtext(1.02*this%xmax-.02*this%xmin,this%offset(i)+off,trim(this%tag(i)))
           if (this%pickflag) then
              call pgsch(1.5)
              do j = 1,3
                 call pgsci(1+j); call pgpt1(this%tpick(i,j),this%offset(i)+off,-4)
              enddo
              call pgsci(this%ci); call pgsch(this%ch)
           endif
        endif
        deallocate(x)
     enddo
  end subroutine plotBasePlotGather
!---------------------------------------------------------------------------
!> \brief Calculate spectra of time series
!
  subroutine computeSpectraBasePlotGather(this)
     class (base_plot_gather) :: this
     integer :: i
     if (.not. allocated(this%spec)) then
        allocate(this%spec(size(this%seis)))
        do i = 1,this%nseis
           call createFromTimeSeriesFourierSpectrum(this%spec(i),this%seis(i))
        enddo
     endif
  end subroutine computeSpectraBasePlotGather
!----------------------------------------------------------------------
!  just do the plotting of the spectra (for internal use)
!
  subroutine plotSpectraBasePlotGather(this,shift)
     class (base_plot_gather) :: this
     integer :: nspec,i,j,shift,nf
     real, dimension(:), allocatable :: x
     real, dimension(:), pointer :: pdata
     real :: off,fmin,df
   !
     call pgslw(3); call pgsch(this%ch)
     call pgslw(this%lw)
     call pgsci(this%ci)
     nspec = this%nseis
     do i=1,nspec
        off = 0.2*this%h*shift
        if(this%offset(i) > this%ymax .or. this%offset(i) < this%ymin) cycle
        nf = .nf.(this%spec(i)); fmin = .fmin.(this%spec(i)); df = .df.(this%spec(i))
        allocate(x(nf))
        forall (j = 1:nf) x(j) = fmin+(j-1)*df
        if (this%specmode == 'a') then
           pdata => .amp.this%spec(i)
           call pgline(nf,x,this%offset(i)+off+pdata(1:nf)/this%smax(i)*this%h)
        else if (this%specmode == 'p') then
           pdata => .phase.this%spec(i)
           call pgline(nf,x,this%offset(i)+off+pdata(1:nf)/this%smax(i)*0.5*this%h)
        endif
        call pgtext(1.02*this%xmax-.02*this%xmin,this%offset(i)+off,trim(this%tag(i)))
        deallocate(x); deallocate(pdata)
     enddo
  end subroutine plotSpectraBasePlotGather
!----------------------------------------------------------------------
!> \brief Set smax for section
!
  subroutine initNormalizationBasePlotGather(this)
     class (base_plot_gather) :: this
     integer :: j
     real, dimension(:), pointer :: pdata
   !
     nullify(pdata)
     if (.not. this%specflag) then
        do j = 1,this%nseis
           this%smax(j) = absmaxTimeSeries(this%seis(j))
           this%somax(j) = this%smax(j)
        enddo
     else
        do j = 1,this%nseis
           if (this%specmode == 'a') then; pdata => .amp.this%spec(j); this%smax(j) = maxval(pdata); deallocate(pdata); endif
           if (this%specmode == 'p') then; pdata => .phase.this%spec(j); this%smax(j) = maxval(abs(pdata)); deallocate(pdata); endif
           this%somax(j) = this%smax(j)
        enddo
     endif
  end subroutine initNormalizationBasePlotGather
!----------------------------------------------------------------------
!> \brief Set x limits for section
!
  subroutine initXLimitsBasePlotGather(this)
     class (base_plot_gather) :: this
     integer :: j
     real :: dmin,dmax
   !
     this%xmin = 100000000.
     this%xmax =-100000000.
     if ((.not. this%specflag)) then
        if (.not.(this%txflag)) then
           do j = 1,this%nseis
              if (this%time_after_trace_start) then
                 this%xmin = min(this%xmin,-this%tref-this%pred*this%offset(j))
                 this%xmax = max(this%xmax,.tend.this%seis(j)-.tanf.this%seis(j)-this%tref-this%pred*this%offset(j))
              else
                 this%xmin = min(this%xmin,.tanf.this%seis(j)-this%tref-this%pred*this%offset(j))
                 this%xmax = max(this%xmax,.tend.this%seis(j)-this%tref-this%pred*this%offset(j))
              endif
           enddo
        else
           dmax = maxval(this%offset); dmin = minval(this%offset)
           this%xmin = dmin-0.5*this%h; this%xmax = dmax+0.5*this%h
        endif
     else
        do j = 1,this%nseis
           this%xmin = min(this%xmin,.fmin.this%spec(j))
           this%xmax = max(this%xmax,.fmax.this%spec(j))
        enddo
     endif
     this%xomin = this%xmin; this%xomax = this%xmax
  end subroutine initXLimitsBasePlotGather
!----------------------------------------------------------------------
!  set y limits for plot
!
  subroutine initYLimitsBasePlotGather(this)
     class (base_plot_gather) :: this
     real :: dmax,dmin
     integer :: j
   !
     if ((.not. this%specflag)) then
        if (.not.(this%txflag)) then
           dmax = maxval(this%offset); dmin = minval(this%offset)
           this%ymin = dmin-0.5*this%h; this%ymax = dmax+0.5*this%h
        else
           this%ymin = 100000000.
           this%ymax =-100000000.
           do j = 1,this%nseis
              if (this%time_after_trace_start) then
                 this%ymin = min(this%ymin,-this%tref-this%pred*this%offset(j))
                 this%ymax = max(this%ymax,.tend.this%seis(j)-.tanf.this%seis(j)-this%tref-this%pred*this%offset(j))
              else
                 this%ymin = min(this%ymin,.tanf.this%seis(j)-this%tref-this%pred*this%offset(j))
                 this%ymax = max(this%ymax,.tend.this%seis(j)-this%tref-this%pred*this%offset(j))
              endif
           enddo
        endif
     else if (this%specflag .and. this%specmode == 'a') then
        dmax = maxval(this%offset); dmin = minval(this%offset)
        this%ymin = dmin-0.25*this%h; this%ymax = dmax+this%h
     endif
     this%yomax = this%ymax; this%yomin = this%ymin
  end subroutine initYLimitsBasePlotGather
!----------------------------------------------------------------------
!> \brief Switch to amplitude spectra plotting
!
  subroutine invokeAmpSpectraBasePlotGather(this)
     class (base_plot_gather) :: this
     this%specflag = .true.
     this%specmode = 'a'
     call this%initNormalization()
     call this%initXLimits()
     call this%initYLimits()
     this%xlab = 'Frequency'
  end subroutine invokeAmpSpectraBasePlotGather
!----------------------------------------------------------------------
!> \brief Switch to phase spectra plotting
!
  subroutine invokePhaseSpectraBasePlotGather(this)
     class (base_plot_gather) :: this
     this%specflag = .true.
     this%specmode = 'p'
     call this%initNormalization()
     call this%initXLimits()
     call this%initYLimits()
     this%xlab = 'Frequency'
  end subroutine invokePhaseSpectraBasePlotGather
!----------------------------------------------------------------------
!> \brief Switch back section plotting
!
  subroutine invokeTimeSeriesBasePlotGather(this)
     class (base_plot_gather) :: this
     this%specflag = .false.
     this%xlab = 'Time'
     call this%initNormalization()
     call this%initXLimits()
     call this%initYLimits()
  end subroutine invokeTimeSeriesBasePlotGather
!----------------------------------------------------------------------
!  set tref
!
  subroutine setTrefBasePlotGather(this,tref)
     class (base_plot_gather) :: this
     real :: tref
     this%tref = tref
     this%xmax = this%xmax-tref
     this%xmin = this%xmin-tref
  end subroutine setTrefBasePlotGather
!----------------------------------------------------------------------
!  unset tref
!
  subroutine resetTrefBasePlotGather(this)
     class (base_plot_gather) :: this
     this%xmax = this%xmax+this%tref
     this%xmin = this%xmin+this%tref
     this%tref = 0.0
  end subroutine resetTrefBasePlotGather
!----------------------------------------------------------------------
!  set pred
!
  subroutine setReductionSlownessBasePlotGather(this,pred)
     class (base_plot_gather) :: this
     real :: pred
     this%pred = pred
  end subroutine setReductionSlownessBasePlotGather
!----------------------------------------------------------------------
!  set plot tag flag
!
  subroutine activateTagBasePlotGather(this)
     class (base_plot_gather) :: this
     this%plottag = .true.
  end subroutine activateTagBasePlotGather
!----------------------------------------------------------------------
!  unset plot tag flag
!
  subroutine deactivateTagBasePlotGather(this)
     class (base_plot_gather) :: this
     this%plottag = .false.
  end subroutine deactivateTagBasePlotGather
!----------------------------------------------------------------------
!  set XY axes limits
!
  subroutine setXYLimitsBasePlotGather(this,x1,x2,y1,y2)
     class (base_plot_gather) :: this
     real :: x1,x2,y1,y2
     this%xmin = x1; this%xmax = x2; this%ymin = y1; this%ymax = y2
  end subroutine setXYLimitsBasePlotGather
!----------------------------------------------------------------------
!  set original axes limits
!
  subroutine resetXYLimitsBasePlotGather(this)
     class (base_plot_gather) :: this
     this%xmin = this%xomin; this%xmax = this%xomax; this%ymin = this%yomin
     this%ymax = this%yomax
     this%h = this%oh
  end subroutine resetXYLimitsBasePlotGather
!----------------------------------------------------------------------
!  set plot scale
!
  subroutine setScaleBasePlotGather(this,s)
     class (base_plot_gather) :: this
     real :: s
     this%h = this%h*s
  end subroutine setScaleBasePlotGather
!----------------------------------------------------------------------
!  set normalization factor of all traces to values in array seismax
!
  subroutine setNormalizationBasePlotGather(this,seismax)
     class (base_plot_gather) :: this
     real, dimension(:) :: seismax
     this%smax = seismax
  end subroutine setNormalizationBasePlotGather
!----------------------------------------------------------------------
!  set seismogram normalization to one value
!
  subroutine setGlobalNormalizationBasePlotGather(this,s)
     class (base_plot_gather) :: this
     real :: s
     this%smax = s
  end subroutine setGlobalNormalizationBasePlotGather
!----------------------------------------------------------------------
!  reset seismogram maximum to individual maximum
!
  subroutine resetNormalizationBasePlotGather(this)
     class (base_plot_gather) :: this
     this%smax = this%somax
  end subroutine resetNormalizationBasePlotGather
!----------------------------------------------------------------------
!  set offsets to data block index
!  recalculate y-limits
!
  subroutine setOffsetTraceIndicesBasePlotGather(this)
     class (base_plot_gather) :: this
     integer :: i
     this%offset = (/ (i, i=1,this%nseis) /)
     call this%initYLimits()
  end subroutine setOffsetTraceIndicesBasePlotGather
!----------------------------------------------------------------------
!  set title
!
  subroutine setTitleBasePlotGather(this,title)
     class (base_plot_gather) :: this
     character (len=*) :: title
     this%title = title
  end subroutine setTitleBasePlotGather
!----------------------------------------------------------------------
!  get xmin
!
  function getXminBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real :: res
     res = this%xmin
  end function getXminBasePlotGather
!----------------------------------------------------------------------
!  get xmax
!
  function getXmaxBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real :: res
     res = this%xmax
  end function getXmaxBasePlotGather
!----------------------------------------------------------------------
!  get ymin
!
  function getYminBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real :: res
     res = this%ymin
  end function getYminBasePlotGather
!----------------------------------------------------------------------
!  get ymax
!
  function getYmaxBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real :: res
     res = this%ymax
  end function getYmaxBasePlotGather
!----------------------------------------------------------------------
!  get four axes limits
!
  subroutine getAxesLimitsBasePlotGather(this,xmin,xmax,ymin,ymax)
     class (base_plot_gather), intent(in) :: this
     real :: xmin,xmax,ymin,ymax
     xmin = this%xmin; xmax = this%xmax; ymin = this%ymin; ymax = this%ymax
  end subroutine getAxesLimitsBasePlotGather
!----------------------------------------------------------------------
!  get smax
!
  function getNormalizationBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real, dimension(:), pointer :: res
     res => this%smax
  end function getNormalizationBasePlotGather
!----------------------------------------------------------------------
!  get somax
!
  function getOriginalNormalizationBasePlotGather(this) result(res)
     class (base_plot_gather), intent(in) :: this
     real, dimension(:), allocatable :: res
     res = this%somax
  end function getOriginalNormalizationBasePlotGather
!------------------------------------------------------------------------
!  add pick to tpick-array
!
  subroutine addPickBasePlotGather(this,x,y,k)
     class (base_plot_gather) :: this
     real :: x,y
     integer :: k
     if (.not.(this%pickflag)) then
        print *,'Pick mode is not activated!'
        return
     endif
     if (this%txflag) then
        this%tpick(nint(x),k) = y
        print *,'Pick type ',k,' set for trace index ',nint(x),' : ',y
     else
        this%tpick(nint(y),k) = x
        print *,'Pick type ',k,' set for trace index ',nint(y),' : ',x
     endif
  end subroutine addPickBasePlotGather
!----------------------------------------------------------------------------
!  write picks to file
!
  subroutine writePicksBasePlotGather(this)
     class (base_plot_gather), intent(in) :: this
     integer :: i,j
     if (this%pickflag) then
        open(1,file = this%pickfile)
        do i = 1,this%nseis
           write(1,*) i,(this%tpick(i,j),j = 1,3)
        enddo
        close(1)
     endif
  end subroutine writePicksBasePlotGather
!
 end module basePlotGather
