!----------------------------------------------------------------------------
!   Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option)
!   any later version.
!
!   ASKI version 1.2 is distributed in the hope that it
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!----------------------------------------------------------------
module globalMpiInfo
    implicit none
    integer :: myrank,numtasks
end module globalMpiInfo