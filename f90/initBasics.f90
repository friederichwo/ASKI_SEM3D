!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
program initBasics
   use inversionBasics
   use iterationStepBasics
   use argumentParser
   use string
   use errorMessage
   implicit none
!
   type (argument_parser) :: ap
   character(len=max_length_string) :: main_parfile
   type (error_message) :: errmsg
   character(len=10) :: myname = 'initBasics'
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------------
   call new(errmsg,myname)
!------------------------------------------------------------------------
!  preliminary processing
!
   call init(ap,myname,"Initiating and testing all basic requirements for ASKI programs (parameter files, event "//&
        "and station list, inversion grid, wavefield points, integration weights, reference model)")
   call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
   end if
   call document(ap)
!
   main_parfile = ap.sval.'main_parfile'
!------------------------------------------------------------------------
!  setup basics
!
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  clean up
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(ap)
 !
!  error treatment
!
 1 if (.level.errmsg == 2) then
      call print(errmsg)
   end if
end program initBasics
