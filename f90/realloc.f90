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
!--------------------------------------------------------------
!  enlarge or downsize an array in memory
!--------------------------------------------------------------
module realloc
   interface reallocate
      module procedure reallocate_real_ar
      module procedure reallocate_double_precision_ar
      module procedure reallocate_complex_ar
      module procedure reallocate_char_ar
      module procedure reallocate_int_ar
      module procedure reallocate_logical_ar
      module procedure reallocate_real_p
      module procedure reallocate_double_precision_p
      module procedure reallocate_complex_p
      module procedure reallocate_char_p
      module procedure reallocate_int_p
      module procedure reallocate_logical_p
      module procedure reallocate_real_2d
      module procedure reallocate_double_precision_2d
      module procedure reallocate_complex_2d
      module procedure reallocate_int_2d
      module procedure reallocate_logical_2d
   end interface reallocate
contains
!--------------------------------------------------------------------
!                     ALLOCATABLE ARRAYS
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_real_ar(p,n) result(res)              ! reallocate real
   real, dimension(:), allocatable :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_real_ar
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_double_precision_ar(p,n) result(res)              ! reallocate double
   double precision, dimension(:), allocatable :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_double_precision_ar
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_complex_ar(p,n) result(res)              ! reallocate complex
   complex, dimension(:), allocatable :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_complex_ar
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_char_ar(p,n) result(res)              ! reallocate char
   character(len=*), dimension(:), allocatable :: p
   character(len=len(p)), dimension(:), allocatable :: res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_char_ar
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_int_ar(p,n) result(res)              ! reallocate int
   integer, dimension(:), allocatable :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_int_ar
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_logical_ar(p,n) result(res)              ! reallocate logical
   logical, dimension(:), allocatable :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. allocated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_logical_ar
!--------------------------------------------------------------------
   function reallocate_real_p(p, n) result(res)               ! reallocate REAL
   real, pointer, dimension(:) :: p, res
   integer, intent(in) :: n
   integer :: nold, ierr
   allocate(res(1:n), stat=ierr)
   if(ierr /= 0) stop "allocate error"
   if(.not. associated(p)) return
   nold = min(size(p), n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_real_p
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_double_precision_p(p,n) result(res)              ! reallocate double
   double precision, dimension(:), pointer :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. associated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_double_precision_p
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_complex_p(p,n) result(res)              ! reallocate complex
   complex, dimension(:), pointer :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. associated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_complex_p
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_char_p(p,n) result(res)              ! reallocate char
   character(len=*), dimension(:), pointer :: p
   character(len=len(p)), dimension(:), pointer :: res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. associated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_char_p
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_int_p(p,n) result(res)              ! reallocate int
   integer, dimension(:), pointer :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. associated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_int_p
!--------------------------------------------------------------------
!  reallocate array p to new dimension n
!  which may be greater or less than size(p)
!
   function reallocate_logical_p(p,n) result(res)              ! reallocate logical
   logical, dimension(:), pointer :: p, res
   integer, intent(in) :: n
   integer :: nold
   allocate(res(1:n))
   if(.not. associated(p)) return
   nold = min(size(p),n)
   res(1:nold) = p(1:nold)
   deallocate(p)
   end function reallocate_logical_p
!--------------------------------------------------------------------
!     2D
!--------------------------------------------------------------------
   function reallocate_double_precision_2d(p, n, k)    ! reallocate DOUBLE PRECISION 2D
   double precision, pointer, dimension(:,:) :: p, reallocate_double_precision_2d
   integer, intent(in) :: n,k
   integer :: nold, kold, ierr
   allocate(reallocate_double_precision_2d(1:n,1:k), stat=ierr)
   if(ierr /= 0) stop "allocate error"
   if(.not. associated(p)) return
   nold = min(size(p,1), n)
   kold = min(size(p,2), k)
   reallocate_double_precision_2d(1:nold,1:kold) = p(1:nold,1:kold)
   deallocate(p)
   end function reallocate_double_precision_2d
!-----------------------------------------------------------------------
   function reallocate_real_2d(p, n, k)    ! reallocate real 2D
   real, pointer, dimension(:,:) :: p, reallocate_real_2d
   integer, intent(in) :: n,k
   integer :: nold, kold, ierr
   allocate(reallocate_real_2d(1:n,1:k), stat=ierr)
   if(ierr /= 0) stop "allocate error"
   if(.not. associated(p)) return
   nold = min(size(p,1), n)
   kold = min(size(p,2), k)
   reallocate_real_2d(1:nold,1:kold) = p(1:nold,1:kold)
   deallocate(p)
   end function reallocate_real_2d
!-----------------------------------------------------------------------
   function reallocate_int_2d(p, n, k)    ! reallocate integer 2D
   integer, pointer, dimension(:,:) :: p, reallocate_int_2d
   integer, intent(in) :: n,k
   integer :: nold, kold, ierr
   allocate(reallocate_int_2d(1:n,1:k), stat=ierr)
   if(ierr /= 0) stop "allocate error"
   if(.not. associated(p)) return
   nold = min(size(p,1), n)
   kold = min(size(p,2), k)
   reallocate_int_2d(1:nold,1:kold) = p(1:nold,1:kold)
   deallocate(p)
   end function reallocate_int_2d
!-----------------------------------------------------------------------
   function reallocate_logical_2d(p, n, k)    ! reallocate logical 2D
   logical, pointer, dimension(:,:) :: p, reallocate_logical_2d
   integer, intent(in) :: n,k
   integer :: nold, kold, ierr
   allocate(reallocate_logical_2d(1:n,1:k), stat=ierr)
   if(ierr /= 0) stop "allocate error"
   if(.not. associated(p)) return
   nold = min(size(p,1), n)
   kold = min(size(p,2), k)
   reallocate_logical_2d(1:nold,1:kold) = p(1:nold,1:kold)
   deallocate(p)
   end function reallocate_logical_2d
!-----------------------------------------------------------------------
    function reallocate_complex_2d(p, n, k)    ! reallocate complex 2D
    complex, pointer, dimension(:,:) :: p, reallocate_complex_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_complex_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_complex_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_complex_2d
!
end module realloc
