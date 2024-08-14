! ===============================================================================
!  User interaction with gather plot
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
!---------------------------------------------------------------------------------
!  User interaction with a gather plot by mouse or keyboard action.
!  Generic module working with any type of plot gather extending
!  basePlotGather.
!--------------------------------------------------------------------
module interactiveBindingsPlotGather
    use basePlotGather
    use changeAxesLimitsPgplot
    implicit none
!
contains
!-------------------------------------------------------------------
!  Allow interactive modification of gather plot.
!  Should be called within a loop to allow repeated user
!  initeraction with the plot
!  gather:     polymorphic gather plot object (in)
!  leave:      flag signifying exiting interactive mode (out)
!  shift:      plot overlaying traces slightly shifted (out)
!  hardcopy:   flag signifying activation of hardcopy mode (out)
!--------------------------------------------------------------------
     subroutine interactWithPlotGather(gather,leave,shift,hardcopy,otherkey)
     class (base_plot_gather), dimension(:) :: gather
     character (len=1) :: ch,otherkey
     integer :: nsec,i,shift,ios
     real :: x1,x2,y1,y2,xa,ya,xmin,xmax,ymin,ymax,s
     real, dimension(:), pointer :: res
     character (len=5) :: sscale
     logical :: leave,success,hardcopy
   !
     leave = .false.
     hardcopy = .false.
     otherkey = 'x'
     nsec = size(gather)
     call gather(1)%getAxesLimits(xmin,xmax,ymin,ymax)
     call pgcurs(xa,ya,ch)
     select case (ch)
   !
   !  time-offset window, time window, or offset window
   !
     case ('A','D','X')
        call pickAxesLimitsPgplot(ch,xa,ya,x1,x2,y1,y2,success)
        if (success) then
           if (ch == 'A') then
              do i=1,nsec
                 call gather(i)%setXYLimits(x1,x2,y1,y2)
              enddo
           else if (ch == 'D') then
              do i=1,nsec
                 call gather(i)%setXYLimits(x1,x2,ymin,ymax)
              enddo
           else if (ch == 'X') then
              do i=1,nsec
                 call gather(i)%setXYLimits(xmin,xmax,y1,y2)
              enddo
           endif
        endif
   !
   !  zoom in, zoom out, scroll left, right, up, down
   !
     case ('z','Z','l','r','u','d')
        call shiftAxesLimitsPgplot(ch,xmin,xmax,ymin,ymax)
        do i=1,nsec
           call gather(i)%setXYLimits(xmin,xmax,ymin,ymax)
        enddo
   !
   !  go back to original axes limits and amplitude scaling
   !
     case ('o')
        do i=1,nsec
           call gather(i)%resetTref()
           call gather(i)%resetXYLimits()
        enddo
   !
   !  scale trace amplitudes, enter s and then a.b, where 0 < a,b < 9
   !
     case ('s')
        do i=1,3
           call pgcurs(xa,ya,ch)
           sscale(i:i) = ch
        enddo
        read(sscale,'(f3.1)',iostat = ios) s
        if (ios > 0) continue
        do i=1,nsec
           call gather(i)%setScale(s)
        enddo
   !
   !  modify normalization
   !
     case ('m')
        call pgcurs(xa,ya,ch)
        if(ch == '1') then                       ! all gathers normalized as first one
           res => gather(1)%getNormalization()
           do i = 2,nsec
              call gather(i)%setNormalization(res)
           enddo
        else if(ch == 'i') then                 ! all traces normalized individually
           do i = 1,nsec
              call gather(i)%initNormalization()
           enddo
        else if(ch == 'a') then                 ! all traces normalized to max of gather 1
           res => gather(1)%getNormalization()
           do i = 1,nsec
              call gather(i)%setGlobalNormalization(maxval(res))
           enddo
        endif
   !
   !  slightly shift overlaying trace
   !
     case ('S'); shift = abs(shift-1)
   !
   !  plot amplitude spectrum of traces
   !
     case ('F')
        do i = 1,nsec
           call gather(i)%computeSpectra()
           call gather(i)%invokeAmpSpectra()
        enddo
   !
   !  plot phase spectrum of traces
   !
     case ('P')
        do i = 1,nsec
           call gather(i)%computeSpectra()
           call gather(i)%invokePhaseSpectra()
        enddo
   !
   !  back to time series plotting
   !
     case ('f')
        do i = 1,nsec
           call gather(i)%invokeTimeSeries()
        enddo
   !
   !  set reduction slowness: t' = t-p*offset
   !  use 5 digits to specify value
   !
     case ('p')
        do i = 1,5
           call pgcurs(xa,ya,ch)
           sscale(i:i) = ch
        enddo
        read(sscale,*,iostat = ios) s
        if (ios > 0) continue
        print *,'Reduction slowness: ',s
        do i = 1,nsec
           call gather(i)%setReductionSlowness(s)
        enddo
   !
   !  activate name tagging
   !
     case ('n')
        do i = 1,nsec
           call gather(i)%activateTag()
        enddo
   !
   !  deactivate name tagging
   !
     case ('N')
        do i = 1,nsec
           call gather(i)%deactivateTag()
        enddo
   !
   !  organize picking of onsets
   !
     case ('k')
     case (' ')
        call pgband(7,0,x1,y1,x2,y2,ch)
        if (ch == 'w') then
           call gather(1)%writePicks()
        else
           read(ch,'(i1)') i
           call gather(1)%addPick(x2,y2,i)
        endif
   !
   !  produce postscript output
   !
     case ('h'); hardcopy = .true.
   !
   !  print documentation of interactive options to screem
   !
     case ('?'); call docuInteractiveBindingsPlotGather()
   !
   !  exit interactive mode
   !
     case ('x'); leave = .true.
   !
   !  some other key was entered
   !
     case default; otherkey = ch
     end select
     end subroutine interactWithPlotGather
!----------------------------------------------------------------------
!  Documentation of interactive options to be printed on screen
!
  subroutine docuInteractiveBindingsPlotGather()
     print '(50(1h-))'
     print *,'          INTERACTIVE PLOT CONTROL'
     print *,''
     print *,'The appearance of the plot may be modified interactively by pressing keys or using the mouse:'
     print *,'?:                   display this help on the screen'
     print *,'x:                   exit interactive mode and finish program'
     print *,'Left mouse (A):      define rectangular zoom window '
     print *,'Middle mouse (D):    define horizontal zoom window'
     print *,'Right mouse (X):     define vertical zoom window'
     print *,'u:                   scroll upwards'
     print *,'d:                   scroll downwards'
     print *,'l:                   scroll to the left'
     print *,'r:                   scroll to the right'
     print *,'z:                   zoom in'
     print *,'Z:                   zoom out'
     print *,'o:                   go back to initial appearance of plot'
     print *,'sa.b                 scale seismograms by a factor a.b (e.g. 5.2 or 0.2)'
     print *,'m1:                  normalize traces of all displayed sections the same way as in section 1'
     print *,'mi:                  normalize traces individually'
     print *,'ma:                  normalize traces to absolute maximum of all traces of section 1'
     print *,'S:                   slightly shift offset of traces of section 2 relative to those of section 1'
     print *,'h:                   generate a postscript file from the current appearance of the plot'
     print *,'F:                   plot amplitude spectra instead of time series'
     print *,'P:                   plot phase spectra instead of time series'
     print *,'f:                   go back to display of time series'
     print *,'pddddd:              enter a reduction slowness s such that tnew = told - s*distance'
     print *,'                     use five digits to specify the value, use dot for floating point values'
     print *,'                     Example: p3.125 or p0.001'
     print *,'n:                   display name tags of traces'
     print *,'N:                   hide name tags of traces'
     print *,'k:                   initialize picking, a cross hair is diaplyed to select pick'
     print *,'e:                   in pick mode: set earliest pick using the cross hair'
     print *,'l:                   in pick mode: set latest pick using cross hair'
     print *,'p:                   in pick mode: set onset time using cross hair'
     print *,'w:                   in pick mode: write picks to a file with name taken from command line'
     print '(50(1h-))'
  end subroutine docuInteractiveBindingsPlotGather
!
end module interactiveBindingsPlotGather

