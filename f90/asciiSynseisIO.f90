! =================================================================================
!  Routines for reading and writing synthetic seismograms in ASCII format
! =================================================================================
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
!--------------------------------------------------------------------------------------
!  Routines for reading and writing synthetic seismograms in ASCII format
!-------------------------------------------------------------------------------------
module asciiSynseisIO
    use string
    implicit none
 contains
!------------------------------------------------------------------------
!  Write synthetic seismograms to file in ASCII format
!  lu:          file unit number (in)
!  filename:    name of file (in)
!  nsamp:       number of samples (in)
!  tanfdp:      start time in seconds after midnight with dp accuracy (in)
!  dt:          sampling interval (in)
!  staname:     station name (in)
!  comp:        component (in)
!  urs:         synthetic data (in)
!  ios:         IO status (out)
!  op:          if true, file is opened (in)
!  cl:          if true, file is closed (in)
!
  subroutine writeAsciiSynseisIO(lu,filename,nsamp,tanfdp,dt,staname,comp,urs,ios,op,cl)
     integer :: lu,ios
     character (len=*) :: filename,staname,comp
     integer :: nsamp
     real :: dt
     double precision :: tanfdp
     real, dimension(:) :: urs
     logical, optional :: op,cl
     logical :: do_open,do_close
   !
   !  open and close file on demand, useful for repeated writes
   !
     if (.not. present(op)) then; do_open  = .true.; else; do_open  = op; endif 
     if (.not. present(cl)) then; do_close = .true.; else; do_close = cl; endif
   !
     if (do_open) then
        open(lu, file = filename,status = 'unknown',action = 'write',form = 'formatted',iostat = ios) 
        if (ios /= 0) return
     endif
     write(lu,'(a)') trim(staname)
     write(lu,'(a)') trim(comp)
     write(lu,*) nsamp,dt,tanfdp
     write(lu,*) urs
     if (do_close) close(lu)
  end subroutine writeAsciiSynseisIO
!------------------------------------------------------------------------
!  Read synthetic seismograms from file in ASCII format
!  lu:          file unit number (in)
!  filename:    name of file (in)
!  nsamp:       number of samples (out)
!  tanfdp:      start time in seconds after midnight with dp accuracy (out)
!  dt:          sampling interval (out)
!  staname:     station name (out)
!  comp:        component (out)
!  urs:         synthetic data (out, allocatable)
!  ios:         IO status (out)
!  op:          if true, file is opened (in)
!  cl:          if true, file is closed (in)
!
  subroutine readAsciiSynseisIO(lu,filename,nsamp,tanfdp,dt,staname,comp,urs,ios,op,cl)
     integer :: lu,ios
     character (len=*) :: filename,staname,comp
     integer :: nsamp
     real :: dt
     double precision :: tanfdp
     real, dimension(:), allocatable :: urs
     logical, optional :: op,cl
     logical :: do_open,do_close
   !
   !  open and close file on demand, useful for repeated reads
   !
     if (.not. present(op)) then; do_open  = .true.; else; do_open  = op; endif 
     if (.not. present(cl)) then; do_close = .true.; else; do_close = cl; endif
   !
     if (do_open) then
        open(lu, file = filename,status = 'old',action = 'read',form = 'formatted',iostat = ios) 
        if (ios > 0) return
     endif
     read(lu,'(a)',iostat = ios) staname
     if (ios < 0) return
     read(lu,'(a)') comp
     read(lu,*) nsamp,dt,tanfdp
     allocate(urs(nsamp))
     read(lu,*) urs
     if (do_close) close(lu)
  end subroutine readAsciiSynseisIO
!----------------------------------------------------------------------------------
!  Find out number of traces contained in data file
!  lu:          file unit number (in)
!  filename:    name of file (in)
!  ios:         error indicator (out)
!
  function getTraceCountAsciiSynseisIO(lu,filename,ios) result(nt)
     integer :: lu,ios
     character (len=*) :: filename
     integer :: nt
     character (len=max_length_string) :: staname
   !
     nt = 0
     open(lu, file = filename,status = 'old',action = 'read',form = 'formatted',iostat = ios) 
     if (ios > 0) return
     ios = 0
     do while (ios == 0)
        read(lu,'(a)',iostat = ios) staname
        if (ios < 0) then                    ! if end of file, close and leave
           close(lu); return
        endif
        read(lu,'(1x)')
        read(lu,'(1x)')
        read(lu,'(1x)')
        nt = nt+1                     ! staname successfully read means one trace more
     enddo
  end function getTraceCountAsciiSynseisIO
!
end module asciiSynseisIO
