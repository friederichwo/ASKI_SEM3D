! ====================================================================================
!----------------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!  Sorting of a 2D-array along a fixed row or column
!  Carries along an index array
!
module heapSort
    implicit none
    interface heapSort2D
        module procedure heapSort2DReal
        module procedure heapSort2DDble
    end interface heapSort2D
    interface siftdownHeapSort2D
        module procedure siftdownHeapSort2DReal
        module procedure siftdownHeapSort2DDble
    end interface siftdownHeapSort2D
contains
!---------------------------------------------------------------------
!       REAL Version
!
!  Perform heap sorting on a row or column of a 2D array.
!  Also adjust elements of other rows or columns in the same way.
!  Also provide an index array giving the original positions of the
!  sorted elements.
!  a:        2D real array (inout)
!  idx:      index array (inout)
!  dimsort:  dimension along which sorting is done (either 1=down or 2=right)
!  is:       row or column controlling the sorting process
!
    subroutine heapSort2DReal(a,idx,dimsort,is)
    real, dimension(:,:) :: a
    integer, dimension(:) :: idx
    integer :: dimsort,is
    integer :: i,n,itmp
    real, dimension(:), allocatable :: atmp
!
! checks if anything to do
!
    n = size(idx)
    if (n < 2) return
!
!  initialize idx to 1,2,3,4,.....
!
    forall (i = 1:n) idx(i) = i
!
! builds heap
!
    do i = n/2,1,-1
       call siftdownHeapSort2D(i,n,a,idx,dimsort,is)
    enddo
!
! sorts array
!
    allocate(atmp(n))
    if (dimsort == 2) then
       do i = n,2,-1
          atmp = a(:,1); a(:,1) = a(:,i); a(:,i) = atmp
          itmp = idx(1); idx(1) = idx(i); idx(i) = itmp
          call siftdownHeapSort2D(1,i-1,a,idx,dimsort,is)
       enddo
    else
       do i = n,2,-1
          atmp = a(1,:); a(1,:) = a(i,:); a(i,:) = atmp
          itmp = idx(1); idx(1) = idx(i); idx(i) = itmp
          call siftdownHeapSort2D(1,i-1,a,idx,dimsort,is)
       enddo
    endif
    deallocate(atmp)
    end subroutine heapSort2DReal
!--------------------------------------------------------------------
!  Sift down algorithm
!
    subroutine siftdownHeapSort2DReal(start,bottom,a,idx,dimsort,is)
    integer, intent(in) :: start
    integer, intent(in) :: bottom
    real, dimension(:,:) :: a
    integer, dimension(:) :: idx
    integer :: is,dimsort
    integer :: i, j
    real, dimension(:), allocatable :: atmp
    integer :: itmp

    allocate(atmp(size(a,dimsort)))
    i = start
    if (dimsort == 2) then
       atmp = a(:,i)
    else
       atmp = a(i,:)
    endif
    itmp = idx(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
       if (j < bottom) then
          if (dimsort == 2) then
             if (a(is,j) <= a(is,j+1)) j = j+1
          else
             if (a(j,is) <= a(j+1,is)) j = j+1
          endif
      endif
    !
    ! checks if section already smaller than initial value
    !
      if (dimsort == 2) then
         if (a(is,j) < atmp(is)) exit
         a(:,i) = a(:,j)
      else
         if (a(j,is) < atmp(is)) exit
         a(i,:) = a(j,:)
      endif
      idx(i) = idx(j)
      i = j
      j = 2 * i
    enddo
    if (dimsort == 2) then
       a(:,i) = atmp
    else
       a(i,:) = atmp
    endif
    idx(i) = itmp
    deallocate(atmp)
    end subroutine siftdownHeapSort2DReal
!---------------------------------------------------------------------
!       DOUBLE PRECISION version
!
!  Perform heap sorting on a row or column of a 2D array.
!  Also adjust elements of other rows or columns in the same way.
!  Also provide an index array giving the original positions of the
!  sorted elements.
!  a:        2D dble array (inout)
!  idx:      index array (inout)
!  dimsort:  dimension along which sorting is done (either 1=down or 2=right)
!  is:       row or column controlling the sorting process
!
    subroutine heapSort2DDble(a,idx,dimsort,is)
    double precision, dimension(:,:) :: a
    integer, dimension(:) :: idx
    integer :: dimsort,is
    integer :: i,n,itmp
    double precision, dimension(:), allocatable :: atmp
!
! checks if anything to do
!
    n = size(idx)
    if (n < 2) return
!
!  initialize idx to 1,2,3,4,.....
!
    forall (i = 1:n) idx(i) = i
!
! builds heap
!
    do i = n/2,1,-1
       call siftdownHeapSort2D(i,n,a,idx,dimsort,is)
    enddo
!
! sorts array
!
    allocate(atmp(n))
    if (dimsort == 2) then
       do i = n,2,-1
          atmp = a(:,1); a(:,1) = a(:,i); a(:,i) = atmp
          itmp = idx(1); idx(1) = idx(i); idx(i) = itmp
          call siftdownHeapSort2D(1,i-1,a,idx,dimsort,is)
       enddo
    else
       do i = n,2,-1
          atmp = a(1,:); a(1,:) = a(i,:); a(i,:) = atmp
          itmp = idx(1); idx(1) = idx(i); idx(i) = itmp
          call siftdownHeapSort2D(1,i-1,a,idx,dimsort,is)
       enddo
    endif
    deallocate(atmp)
    end subroutine heapSort2DDble
!--------------------------------------------------------------------
!  Sift down algorithm
!
    subroutine siftdownHeapSort2DDble(start,bottom,a,idx,dimsort,is)
    integer, intent(in) :: start
    integer, intent(in) :: bottom
    double precision, dimension(:,:) :: a
    integer, dimension(:) :: idx
    integer :: is,dimsort
    integer :: i, j
    double precision, dimension(:), allocatable :: atmp
    integer :: itmp

    allocate(atmp(size(a,dimsort)))
    i = start
    if (dimsort == 2) then
       atmp = a(:,i)
    else
       atmp = a(i,:)
    endif
    itmp = idx(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
       if (j < bottom) then
          if (dimsort == 2) then
             if (a(is,j) <= a(is,j+1)) j = j+1
          else
             if (a(j,is) <= a(j+1,is)) j = j+1
          endif
      endif
    !
    ! checks if section already smaller than initial value
    !
      if (dimsort == 2) then
         if (a(is,j) < atmp(is)) exit
         a(:,i) = a(:,j)
      else
         if (a(j,is) < atmp(is)) exit
         a(i,:) = a(j,:)
      endif
      idx(i) = idx(j)
      i = j
      j = 2 * i
    enddo
    if (dimsort == 2) then
       a(:,i) = atmp
    else
       a(i,:) = atmp
    endif
    idx(i) = itmp
    deallocate(atmp)
    end subroutine siftdownHeapSort2DDble
end module heapSort
