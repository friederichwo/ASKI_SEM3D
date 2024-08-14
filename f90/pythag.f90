!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------
!	calculate sqrt(x**2+y**2) in a robust way
!--------------------------------------------------------
    double precision function pythag(x,y)
    double precision ::  x,y,ax,ay,r
!
!  absolute values
!
    ax = abs(x)
    ay = abs(y)
!    
!  case where abs(x) >= abs(y)
!
    if (ax .ge. ay) then
        if (ax .eq. 0.d0) then
            pythag = 0.d0
        else
            r = ay/ax
            pythag = ax*sqrt(r*r+1.)
        endif
!
!  case where abs(x) < abs(y)
!
    else
        r = ax/ay
        pythag = ay*sqrt(r*r+1.)
    endif
    return
    end function pythag
