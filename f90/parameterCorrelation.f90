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
!   Defines a correlation between the physical properties in the earth and the
!   properties we invert for. For example, we may only invert for vp but assume that
!   a change in vp may also imply a change of vs and rho in the earth. This means,
!   we pull along density and shear wave velocity with P-wave velocity.
!   If we assume such a correlation, it has consequences for kernel computation as
!   well as for the model update.
!   Assume that the true perturbations in the earth are given by delrho, delvp and delvs,
!   while our inverted properties are drho,dvp and dvs. We can define a correlation between
!   RELATIVE model changes as follows:
!   delrho/rho = b(1,1)*drho/rho + b(1,2)*dvp/vp + b(1,3)*dvs/vs
!   delvp/vp   = b(2,1)*drho/rho + b(2,2)*dvp/vp + b(2,3)*dvs/vs
!   delvs/vs   = b(3,1)*drho/rho + b(3,2)*dvp/vp + b(3,3)*dvs/vs.
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
!   Inserting our correlation with the properties inverted for, we get:
!   delu =   Krho*(a(1,1)*drho + a(1,2)*dvp + a(1,3)*dvs)
!          + Kvp *(a(2,1)*drho + a(2,2)*dvp + a(2,3)*dvs)
!          + Kvs *(a(3,1)*drho + a(3,2)*dvp + a(3,3)*dvs)
!        =   (a(1,1)*Krho+a(2,1)*Kvp+a(3,1)*Kvs)*drho
!          + (a(1,2)*Krho+a(2,2)*Kvp+a(3,2)*Kvs)*dvp
!          + (a(1,3)*Krho+a(2,3)*Kvp+a(3,3)*Kvs)*dvs
!        =  Krho'*drho + Kvp'*dvp + Kvs'*dvs
!   with the new kernels given by
!   Krho' = a(1,1)*Krho+a(2,1)*Kvp+a(3,1)*Kvs
!   Kvp'  = a(1,2)*Krho+a(2,2)*Kvp+a(3,2)*Kvs
!   Kvs'  = a(1,3)*Krho+a(2,3)*Kvp+a(3,3)*Kvs
!
!   We specify the correlation information through a file built as follows:
!   Header Line: model parameterization
!   Line 1: property name followed by correlation coefficients  (b(1,j),j = 1,npar)
!   Line 2: property name followed by correlation coefficients  (b(2,j),j = 1,npar)
!   Line 3: property name followed by correlation coefficients  (b(3,j),j = 1,npar)
!   Line X: property name followed by correlation coefficients  (b(X,j),j = 1,npar)
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
!  We can find out which parameters we want to evaluate: if any entry in the i-th row is
!  non-zero, we want the i-th parameter.
!!
!! \author Florian Schumacher, Wolfgang Friederich
!! \date Nov 2015, 2021
!!
module parameterCorrelation
!
  use modelParametrization
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateParameterCorrelation; end interface
  interface correlateParameters
     module procedure getCorrByIndexParameterCorrelation
     module procedure getCorrByNameParameterCorrelation
  end interface correlateParameters
!
  type parameter_correlation
     character(len=char_len_pmtrz) :: model_parametrization
     double precision, dimension(:,:), pointer :: corr_coef => null()            ! correlation matrix B
  end type parameter_correlation
!
contains
!
!------------------------------------------------------------------------
!  Read correlation coefficients from file
!
  subroutine createParameterCorrelation(this,filename,lu,parametrization,errmsg)
    type (parameter_correlation) :: this
    character(len=*) :: filename,parametrization
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=26) :: myname = 'createParameterCorrelation'
    character(len=400) :: errstr,line
    character (char_len_par) :: param
    integer :: nparam,iparam
    integer :: ios,j
!
    call addTrace(errmsg,myname)
    call deallocateParameterCorrelation(this)
!
    open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not open file, iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    read(lu,'(a)',iostat=ios) line         ! parameterization for correlation
    if(ios /= 0) goto 2
    if (trim(adjustl(line)) .ne. trim(adjustl(parametrization))) then
       write(errstr,*) "parameterization of inversion conflicts with that of correlation"
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
!
    if(.not.validModelParametrization(parametrization)) then
       call add(errmsg,2,"incoming model parametrization '"//trim(parametrization)//"' is not valid",myname)
       goto 2
    end if
!
    ! allocate for this model parametrization
    nparam = numberOfParamModelParametrization(parametrization)
    allocate(this%corr_coef(nparam,nparam))
    this%corr_coef = 0.
    this%model_parametrization = parametrization
!
! read parameter correlation file
!
    do j = 1,nparam
       read(lu,"(a400)",iostat=ios) line                     ! parameters with corr_coef
       param = adjustl(line(1:index(line," ")-1))
       if(validParamModelParametrization(parametrization,param)) then
          iparam = indexOfParamModelParametrization(parametrization,param)
          read(line(index(line," "):),*) this%corr_coef(iparam,:)
       endif
    enddo
    return
2   call add(errmsg,2,"object could not be created",myname)
    call deallocateParameterCorrelation(this)
  end subroutine createParameterCorrelation
!------------------------------------------------------------------------
  subroutine deallocateParameterCorrelation(this)
    type (parameter_correlation) :: this
    if (associated(this%corr_coef)) deallocate(this%corr_coef)
    this%model_parametrization = ''
  end subroutine deallocateParameterCorrelation
!------------------------------------------------------------------------
!  check if any non-self correlation is to be done
!  by looking for non-zero off-diagonal elements of correlation matrix
!  If param is present, only check for correlation with this parameter
!
  function anyParameterCorrelation(this,param) result(res)
    type (parameter_correlation) :: this
    logical :: res
    character (len=*), optional :: param
    integer :: i,j,iparam,npar
    !
    res = .false.
    npar = size(this%corr_coef,1)
    if (present(param)) then
       iparam = indexOfParamModelParametrization(this%model_parametrization,param)
       do i = 1,npar
          if (i .ne. iparam .and. abs(this%corr_coef(i,iparam)) .gt. 1.e-6) then
             res = .true.
             exit
          endif
       enddo
    else
        do i = 1,npar
           do j = 1,npar
              if (i .ne. j .and. abs(this%corr_coef(i,j)) .gt. 1.e-6) then
                 res = .true.
              endif
              exit
           enddo
        enddo
    endif
  end function anyParameterCorrelation
!------------------------------------------------------------------------
! return element of correlation matrix b(i,j)
!
  function getCorrByIndexParameterCorrelation(this,i,j,errmsg) result(c_corr)
    type (parameter_correlation) :: this
    type (error_message) :: errmsg
    character(len=34) :: myname = 'getCorrByIndexParameterCorrelation'
    integer :: i,j
    double precision :: c_corr
    c_corr = 0.d0
    if(i < 1 .or. i > size(this%corr_coef,1)) then
       call add(errmsg,2,"invalid first parameter index",myname)
       return
    endif
    if(j < 1 .or. j > size(this%corr_coef,1)) then
       call add(errmsg,2,"invalid second parameter index",myname)
       return
    endif
    c_corr = this%corr_coef(i,j)
  end function getCorrByIndexParameterCorrelation
!------------------------------------------------------------------------
! return element of correlation matrix by giving first and second property name
! corresponding to indices i and j in b(i,j)
!
  function getCorrByNameParameterCorrelation(this,name_param_i,name_param_j,errmsg) result(c_corr)
    type (parameter_correlation) :: this
    type (error_message) :: errmsg
    character(len=*) :: name_param_i,name_param_j
    double precision :: c_corr
    integer :: i,j
    i = indexOfParamModelParametrization(this%model_parametrization,name_param_i)
    j = indexOfParamModelParametrization(this%model_parametrization,name_param_j)
    c_corr = getCorrByIndexParameterCorrelation(this,i,j,errmsg)
  end function getCorrByNameParameterCorrelation
!------------------------------------------------------------------------
!> \brief for a given model parameter, iterate over all parameters which are correlated to it
!! \param  model parametrization
!! \param param this parameter correlation
!! \param 
!! \return logical value which is false if there is no next model parameter
!
  logical function nextCorrelationParameter(this,name_param_main,name_param_corr,c_corr)
    type (parameter_correlation) :: this
    character(len=*) :: name_param_main,name_param_corr
    double precision, optional :: c_corr
    integer :: iparam_corr,iparam
    integer :: count = 0
    integer :: iparam_main = 0
    save :: count,iparam_main
!
    ! COUNT IS NOT A CALL COUNT, but rather is the index iparam_corr of the latest found parameter
    ! which is to be correlated to param_main (or 0 on first call), hence this index can also jump from
    ! one call to another:  
    !   e.g. if for a given main parameter the parameters 1 , 2 and 5 correlate,
    !   then count will have values 1, 2 and 5 on exits of this routine
!
    if(iparam_main==0) then
       ! if this is the first call, memorize the main parameter (index)
       iparam_main = indexOfParamModelParametrization(this%model_parametrization,name_param_main)
    else
       ! otherwise, if the main parameter has changed, return false and set everything back to start
       if(iparam_main /= indexOfParamModelParametrization(this%model_parametrization,name_param_main)) goto 100
    end if
!
    ! now find the next parameter which is to be correlated to the param_main
    ! EXCLUDING the case of self-correlation:
    ! on entry, count is the index iparam_corr found at the last call (or 0 at first call)
    ! so start loop on parameter indices at count+1
    iparam_corr = -1
    do iparam = count+1,numberOfParamModelParametrization(this%model_parametrization)
        if (iparam_main .ne. iparam .and. abs(this%corr_coef(iparam,iparam_main)) .gt. 1.e-6) then
           iparam_corr = iparam
           exit
        endif
    end do
    ! if there was no correlation index found, set everything back to start
    if(iparam_corr == -1) goto 100
!
    ! otherwise store the relevant values
    name_param_corr = getParamFromIndexModelParametrization(this%model_parametrization,iparam_corr)
    if(present(c_corr)) c_corr = this%corr_coef(iparam_corr,iparam_main)
    ! memorize the index of the parameter to be correlated. 
    count = iparam_corr
!
    ! if routine comes here, everything went alright, so return memorizing
    nextCorrelationParameter = .true.
    return
!
    ! if there is a goto 100, the outome is negative, so set everything back to start
100 count = 0
    iparam_main = 0
    name_param_corr = ''
    if(present(c_corr)) c_corr = 0.d0
    nextCorrelationParameter = .false.
  end function nextCorrelationParameter
!
end module parameterCorrelation
