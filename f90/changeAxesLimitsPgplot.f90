! ===============================================================================
!  Routines to change plot axes limits interactively
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
!----------------------------------------------------------------------
!  Routines to interactively change plot axes limits of potentially any 
!  figure based on the pgplot library. Routines return new values for
!  the axes limits as chosen interactively by the user.
!----------------------------------------------------------------------
module changeAxesLimitsPgplot
    implicit none
contains
!----------------------------------------------------------------------------------
!  General bindings for interactively selecting new axes limits
!  Assume that the lower left limits have already been picked
!  using either left, middle or right mouse or keys A,D or X.
!  ch:                   character specifying initial action (A, D or X) (in)
!  xa,ya:                coordinates of lower left corner of previous pick (in)
!  xmin,xmax,ymin,ymax:  new extremal values for plot range (inout)
!  success:              flag indicating whether selection was successful
!
  subroutine pickAxesLimitsPgplot(ch,xa,ya,xmin,xmax,ymin,ymax,success)
     real :: xa,ya,xmin,xmax,ymin,ymax
     logical :: success
     real :: x,y
     character (len=1) :: ch,ch1
   !
     success = .false.
     select case (ch)
     case ('A')                                       ! select axes limits by picking corners of a rectangular box
        call pgslw(1); call pgsci(5)
        call pgband(2,0,xa,ya,x,y,ch1)
        if (ch1 == ch) then
           xmin = xa; ymin = ya; xmax = x; ymax = y
           success = .true.
        endif
     case ('D')                                       ! only select new limits for x-axis
        call pgslw(1); call pgsci(5)
        call pgband(4,0,xa,ya,x,y,ch1)
        if (ch1 == ch) then
           xmin = xa; xmax = x
           success = .true.
        endif
     case ('X')                                       ! only select new limits for y-axis
        call pgslw(1); call pgsci(5)
        call pgband(3,0,xa,ya,x,y,ch1)
        if (ch1 == ch) then
           ymin = ya; ymax = y
           success = .true.
        endif
     end select
  end subroutine pickAxesLimitsPgplot
!---------------------------------------------------------------------------------
!  Shift or enlarge or shrink axes limits
!  ch:                  character specifying action
!  xmin,xmax,ymin,ymax: old and new limits (inout)
!
  subroutine shiftAxesLimitsPgplot(ch,xmin,xmax,ymin,ymax)
     character (len=*) :: ch
     real :: xmin,xmax,ymin,ymax
     real :: xc,yc,hx,hy
   !
     select case (ch)
     case ('z')                                          !  zoom in
        xc = 0.5*(xmin+xmax); yc = 0.5*(ymin+ymax)     
        hx = xmax-xmin; hy = ymax-ymin
        xmin = xc-0.75*hx; xmax = xc+0.75*hx
        ymin = yc-0.75*hy; ymax = yc+0.75*hy
     case ('Z')                                          !  zoom out
        xc = 0.5*(xmin+xmax); yc = 0.5*(ymin+ymax)     
        hx = xmax-xmin; hy = ymax-ymin
        xmin = xc-0.25*hx; xmax = xc+0.25*hx
        ymin = yc-0.25*hy; ymax = yc+0.25*hy
     case ('l')                                          !  scroll left
        hx = xmax-xmin
        xmin = xmin - 0.25*hx 
        xmax = xmax - 0.25*hx
     case ('r')                                          !  scroll right
        hx = xmax-xmin
        xmin = xmin + 0.25*hx
        xmax = xmax + 0.25*hx
     case ('u')                                          !  scroll up
        hy = ymax-ymin
        ymin = ymin + 0.25*hy
        ymax = ymax + 0.25*hy
     case ('d')                                          !  scroll down
        hy = ymax-ymin
        ymin = ymin - 0.25*hy
        ymax = ymax - 0.25*hy
     end select
  end subroutine shiftAxesLimitsPgplot
!
end module changeAxesLimitsPgplot
