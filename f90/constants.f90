!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------------------
!--------------------------------------------------
!  Module with various constants used in ASKI package
!--------------------------------------------------
module constants
!
! length of diverse names
!
   integer, parameter :: char_len_comp = 2
   integer, parameter :: char_len_sta = 8
   integer, parameter :: char_len_netcode = 2
   integer, parameter :: char_len_evid = 17
   integer, parameter :: char_len_par = 6
   integer, parameter :: char_len_pmtrz = 15
   integer, parameter :: char_len_netsta = char_len_sta+char_len_netcode+1
!
!  max length of a string
!
   integer, parameter :: max_length_string = 400
!
!  length of ASKI output ID (used in SPECFEM)
!
   integer, parameter :: char_len_aski_output_id = 17
!
!  valid components for Green tensor
!
   character(len=char_len_comp), dimension(9), parameter :: valid_components = ['CX','CY','CZ','N ','S ','E ','W ','UP','DO']
end module constants
