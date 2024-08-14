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
!   Define (an)elastic material properties
!
!   Organizes (an)elastic properties of earth models by defining
!   the number of properties n used in a property set and the names of the properties,
!   whereby each property is assigned a specific index between 1 and n
!
!   \author Florian Schumacher
!   \date Nov 2015
!----------------------------------------------------------------------------
!   Defines a correlation between the physical properties in the earth and the
!   properties we invert for. For example, we may only invert for vp but assume that
!   a change in vp may also imply a change of vs and rho in the earth. This means,
!   we pull along density and shear wave velocity with P-wave velocity.
!   If we assume such a correlation, it has consequences for kernel computation as
!   well as for the model update.
!
!   Assume that the true perturbations in the earth are given by delrho, delvp and delvs,
!   while our inverted properties are drho,dvp and dvs. We can define a correlation between
!   RELATIVE model changes as follows:
!   delrho/rho = b(1,1)*drho/rho + b(1,2)*dvp/vp + b(1,3)*dvs/vs
!   delvp/vp   = b(2,1)*drho/rho + b(2,2)*dvp/vp + b(2,3)*dvs/vs
!   delvs/vs   = b(3,1)*drho/rho + b(3,2)*dvp/vp + b(3,3)*dvs/vs.
!
!   From this, we obtain a correlation between absolute model changes as follows
!   delrho = a(1,1)*drho + a(1,2)*dvp + a(1,3)*dvs
!   delvp  = a(2,1)*drho + a(2,2)*dvp + a(2,3)*dvs
!   delvs  = a(3,1)*drho + a(3,2)*dvp + a(3,3)*dvs,
!   with a(i,j) = b(i,j)*p_i/p_j, where p_i is any of the three parameters.
!   Thus, if our inversion result is for example dvp we should use the a(i,2)-coefficients
!   to update the parameters.
!
!   For kernel computation, we have a relation between change of displacement delu and
!   the model changes delrho, delvp and delvs:
!   delu = Krho*delrho+Kvp*delvp+Kvs*delvs.
!
!   Inserting our correlation with the properties inverted for, we get:
!   delu =   Krho*(a(1,1)*drho + a(1,2)*dvp + a(1,3)*dvs)
!          + Kvp *(a(2,1)*drho + a(2,2)*dvp + a(2,3)*dvs)
!          + Kvs *(a(3,1)*drho + a(3,2)*dvp + a(3,3)*dvs)
!        =   (a(1,1)*Krho+a(2,1)*Kvp+a(3,1)*Kvs)*drho
!          + (a(1,2)*Krho+a(2,2)*Kvp+a(3,2)*Kvs)*dvp
!          + (a(1,3)*Krho+a(2,3)*Kvp+a(3,3)*Kvs)*dvs
!        =  Krho'*drho + Kvp'*dvp + Kvs'*dvs
!
!   with the new kernels given by
!   Krho' = a(1,1)*Krho+a(2,1)*Kvp+a(3,1)*Kvs
!   Kvp'  = a(1,2)*Krho+a(2,2)*Kvp+a(3,2)*Kvs
!   Kvs'  = a(1,3)*Krho+a(2,3)*Kvp+a(3,3)*Kvs
!
!   We specify the correlation information through a file built as follows:
!   Header Line: name of property set
!   Line 1: property name followed by correlation coefficients  (b(1,j),j = 1,nprop)
!   Line 2: property name followed by correlation coefficients  (b(2,j),j = 1,nprop)
!   Line 3: property name followed by correlation coefficients  (b(3,j),j = 1,nprop)
!   Line X: property name followed by correlation coefficients  (b(X,j),j = 1,nprop)
!   The coefficients are stored in a matrix.
!   The file should be placed in the main inversion folder to ensure consistency over
!   several iterations.
!
!   Example 1: We invert for vp only and want to drag rho and vs along:
!   rho  0.0   0.3   0.0      ->   delrho/rho = 0.3*dvp/vp      and Krho' = 0
!   vp   0.0   1.0   0.0      ->   delvp/vp   = dvp/vp          and Kvp'  = 0.3*rho/vp*Krho + Kvp + 1.8*vs/vp*Kvs
!   vs   0.0   1.8   0.0      ->   delvs/vs   = 1.8*dvp/vp      and Kvs'  = 0
!
!   Example 2: We invert for vp only and want to leave rho and vs unchanged:
!   rho  0.0   0.0   0.0      ->   delrho/rho = 0.0             and Krho' = 0
!   vp   0.0   1.0   0.0      ->   delvp/vp   = dvp/vp          and Kvp'  = Kvp
!   vs   0.0   0.0   0.0      ->   delvs/vs   = 0.0             and Kvs'  = 0
!
!  We can find out which kernels we need to compute from this matrix. If any entry in the
!  i-th row is non-zero, we need the kernel K_i.
!  We can find out which properties we want to evaluate: if any entry in the i-th row is
!  non-zero, we want the i-th property.
! ------------------------------------------------------------------------
module propertySet
   use constants
   use string
   implicit none
   interface dealloc
      module procedure deallocPropertySet
   end interface
   interface operator (.setname.)
      module procedure getNamePropertySet
   end interface
   interface operator (.name.)
      module procedure getNameFromIndexPropertySet
   end interface
   interface operator (.index.)
      module procedure getIndexFromNamePropertySet
   end interface
   interface operator (.nprop.)
      module procedure getNpropPropertySet
   end interface
!
   type property_set
      character(len=14) :: name                                                ! name of property set
      integer :: numprop                                                       ! number of properties
      character(len=char_len_par), dimension(:), allocatable :: properties     ! names of properties
      double precision, dimension(:,:), allocatable :: corrmat                 ! correlation matrix
   end type
!
contains
! -------------------------------------------------------------------------------------------
!  Create an isoVelocitySI property set
!
   subroutine createIsoVelocitySIPropertySet(this)
      type (property_set) :: this
      integer :: i
      this%name = 'isoVelocitySI'
      this%properties = ['rho','vp ','vs ']
      this%numprop = 3
      allocate(this%corrmat(3,3))
      this%corrmat = 0.d0
      forall (i = 1:3) this%corrmat(i,i) = 1.d0
   end subroutine createIsoVelocitySIPropertySet
! -------------------------------------------------------------------------------------------
!  Create an new property set
!
   subroutine createPropertySet(this,setname,props)
      type (property_set) :: this
      character(len=*) :: setname
      character(len=char_len_par), dimension(:) :: props
      integer :: i
      this%name = trim(setname)
      this%properties = props
      this%numprop = size(props)
      allocate(this%corrmat(this%numprop,this%numprop))
      this%corrmat = 0.d0
      forall (i = 1:this%numprop) this%corrmat(i,i) = 1.d0
   end subroutine createPropertySet
! -------------------------------------------------------------------------------------------
!  Deallocate property set
!
   subroutine deallocPropertySet(this)
      type (property_set) :: this
      if (allocated(this%properties)) deallocate(this%properties)
      if (allocated(this%corrmat)) deallocate(this%corrmat)
   end subroutine deallocPropertySet
! -------------------------------------------------------------------
!  Read correlation matrix from main ASKI input parameters
!
   subroutine readCorrmatPropertySet(this,lu,filename,errmsg)
      type (property_set) :: this
      integer :: lu
      character(len=*) :: filename
      type(error_message) :: errmsg
      integer :: ios,j,iprop
      character(len=char_len_par) :: prop
      character(len=max_length_string) :: line
      character(len=21) :: myname = 'readCorrmatPropertySet'
   !
      open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
      if(ios /= 0) then
         call add(errmsg,2,'can not open file '//trim(filename),myname)
         return
      endif
   !
   !  read name of property set and check if it fits
   !
      read(lu,'(a)',iostat=ios) line
      if (ios /= 0) then
         call add(errmsg,2,'problems reading file '//trim(filename),myname)
         return
      endif
      if (.not. equalString(line,this%name)) then
         call add(errmsg,2,'property set '//trim(line)//' in file does not fit',myname)
         return
      endif
   !
   !  read correlation matrix
   !  each line contains property name and row of correlation matrix for this property
   !  correlation matrix is filled according to order of properties define upon construction
   !
      do j = 1,this%numprop
         read(lu,'(a)',iostat=ios) line
         if (ios /= 0) then
            call add(errmsg,2,'problems reading file '//trim(filename),myname)
            return
         endif
         prop = adjustl(line(1:index(line," ")-1))
         if (.not. isValidNamePropertySet(this,prop)) then
            call add(errmsg,2,'invald property name in '//trim(filename),myname)
            return
         end if
         iprop = this.index.prop
         read(line(index(line," "):),*) this%corrmat(iprop,:)
      enddo
      close(lu)
   end subroutine readCorrmatPropertySet
!----------------------------------------------------------------------------------------------
!  Check validity of correlation matrix given the properties to be inverted
!  If a property is not inverted for, its column in the correlation matrix must be all zero
!  If a property is inverted for, its diagonal entry must not be zero
!  prop:      properties to be inverted for
!
   function checkCorrmatPropertySet(this,prop,errmsg) result(res)
      type (property_set) :: this
      character(len=char_len_par), dimension(:) :: prop
      type (error_message) :: errmsg
      logical :: res
      integer :: j,idx
      character(len=22) :: myname = 'checkCorrmatPropertySet'
   !
   !  first check validity of properties to be inverted for
   !
      do j = 1,size(prop)
         if (.not. isValidNamePropertySet(this,prop(j))) then
            res = .false.
            call add(errmsg,2,'Inverted property name not in property set',myname)
            return
         end if
      end do
   !
   !  now check diagonal element in correlation matrix for properties inverted for
   !  diagonal element must not be zero
   !
      do j = 1,size(prop)
         idx = getIndexFromNamePropertySet(this,prop(j))
         if (abs(this%corrmat(idx,idx)) < 1.d-6) then
            res = .false.
            call add(errmsg,2,'Corrmat diagonal entry for property '//trim(prop(j))//' is zero',myname)
            return
         end if
      end do
   !
   !  now check properties not inverted for
   !  their columns must be all zero
   !
      do j = 1,this%numprop
         if (.not. any(prop == this%properties(j))) then
            if (maxval(abs(this%corrmat(:,j))) > 1.d-6) then
               res = .false.
               call add(errmsg,2,'Corrmat column of not inverted-for property '//this%properties(j)//' is non-zero',myname)
               return
            end if
         end if
      end do
      res = .true.
   end function checkCorrmatPropertySet
! ---------------------------------------------------------------------------------------------
!  Get number of properties
!
   function getNpropPropertySet(this) result(n)
      type (property_set), intent(in) :: this
      integer :: n
      n = this%numprop
   end function getNpropPropertySet
! ---------------------------------------------------------------------
!  Get pointer to property names
!
   function getNamesPropertySet(this) result(res)
      type (property_set), target :: this
      character(len=char_len_par), dimension(:), pointer :: res
      res => this%properties
   end function getNamesPropertySet
! ---------------------------------------------------------------------
!  Get pointer to the correlation matrix
!
   function getCorrmatPropertySet(this) result(res)
      type (property_set), target :: this
      double precision, dimension(:,:), pointer :: res
      res => this%corrmat
   end function getCorrmatPropertySet
!------------------------------------------------------------------------
!  Get the name of the property set
!
   function getNamePropertySet(this) result(name)
      type (property_set), intent(in) :: this
      character(len=char_len_pmtrz):: name
      name = this%name
   end function getNamePropertySet
!------------------------------------------------------------------------
!  For given index get the name of the property
!  For performance, no error check here
!
   function getNameFromIndexPropertySet(this,idx) result(name)
      type (property_set), intent(in) :: this
      integer, intent(in) :: idx
      character(len=char_len_par):: name
      name = this%properties(idx)
   end function getNameFromIndexPropertySet
! --------------------------------------------------------------------------
!  Get index from property name
!
   function getIndexFromNamePropertySet(this,name) result(idx)
      type (property_set), intent(in) :: this
      character(len=char_len_par), intent(in) :: name
      integer :: idx
      idx = findloc(this%properties,name,1)
   end function getIndexFromNamePropertySet
!------------------------------------------------------------------------
!  Check if given property name is defined in property set
!
   function isValidNamePropertySet(this,prop) result(l)
      type (property_set) :: this
      character(len=*) :: prop
      logical :: l
      l = any(this%properties == prop)
   end function isValidNamePropertySet
!------------------------------------------------------------------------
!  Iterate over properties of a given property set
!
   logical function nextNamePropertySet(this,prop,iprop)
      type (property_set) :: this
      character(len=char_len_par), optional :: prop
      integer, optional :: iprop
      integer :: call_count = 0
      save :: call_count
   !
      call_count = call_count+1
      if (call_count <= this%numprop) then
         if (present(prop)) prop = getNameFromIndexPropertySet(this,call_count)
         if (present(iprop)) iprop = call_count
         nextNamePropertySet = .true.
      else
         call_count = 0
         if (present(prop)) prop = ''
         if (present(iprop)) iprop = 0
         nextNamePropertySet = .false.
      endif
   end function nextNamePropertySet
! ---------------------------------------------------------------------
!  Sort incoming set of properties following order in definition of property set
!  Performs in-place sorting and returns modified property list
!  Assumes that incoming properties are distinct and unique
!
   subroutine sortNamesPropertySet(this,prop)
      type (property_set) :: this
      character(len=char_len_par), dimension(:) :: prop
      character(len=char_len_par), dimension(:), allocatable :: h_prop
      character(len=char_len_par) :: p
      integer :: np,ip,j
   !
      np = size(prop)
      allocate(h_prop(np))
      ip = 0
      do while(nextNamePropertySet(this,p))
         do j = 1,np
            if (equalString(p,prop(j))) then
               ip = ip+1
               h_prop(ip) = p
               exit
            end if
         end do
      end do
      prop = h_prop
      deallocate(h_prop)
   end subroutine sortNamesPropertySet
! ---------------------------------------------------------------------
!  Based on the correlation matrix, return property names for which
!  kernels need to be calculated
!  If there is non-zero entry in ROW i, we compute kernel i
!
   subroutine getRequiredKernelsByNamePropertySet(this,names)
      type (property_set) :: this
      character(len=char_len_par), dimension(:), allocatable :: names
      integer :: nr,i
      character(len=char_len_par), dimension(:), allocatable :: h_names
   !
      allocate(h_names(this%numprop))
      nr = 0
      do i = 1,this%numprop
         if (count((abs(this%corrmat(i,:)) > 1.e-6),1)  > 0) then
            nr = nr+1
            h_names(nr) = this%properties(i)
         end if
      end do
      allocate(names(nr))
      names = h_names(1:nr)
      deallocate(h_names)
   end subroutine getRequiredKernelsByNamePropertySet
! ---------------------------------------------------------------------
!  Based on the correlation matrix, return property names for which
!  kernels are finally written to file
!  If there is a non-zero entry in COLUMN i, we output kernel i
!
   subroutine getOutputKernelsByNamePropertySet(this,names)
      type (property_set) :: this
      character(len=char_len_par), dimension(:), allocatable :: names
      integer :: nr,i
      character(len=char_len_par), dimension(:), allocatable :: h_names
   !
      allocate(h_names(this%numprop))
      nr = 0
      do i = 1,this%numprop
         if (count((abs(this%corrmat(:,i)) > 1.e-6),1)  > 0) then
            nr = nr+1
            h_names(nr) = this%properties(i)
         end if
      end do
      allocate(names(nr))
      names = h_names(1:nr)
      deallocate(h_names)
   end subroutine getOutputKernelsByNamePropertySet
!
end module propertySet
