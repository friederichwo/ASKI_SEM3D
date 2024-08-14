!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2022 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   Module handling earth models resulting from ASKI full waveform inversion
!
module kernelInvertedModel
   use inversionGrid
   use mathConstants
   use askiBackgroundModel
   use fmtomoModel
   use propertySet
   use hdfWrapper
   use axesRotation
   use globalHdfInfo

   implicit none
   interface dealloc; module procedure deallocKernelInvertedModel; end interface
   type kernel_inverted_model
      integer :: ncell                                               ! number of grid cells kim is defined on
      integer :: nprop                                               ! number of model properties
      character(len=char_len_par), dimension(:), allocatable :: prop ! property names
      character(len=max_length_string) :: propsetname                ! name property set
      real, dimension(:,:), pointer :: model_values                  ! model values on invgrid (ncell,nprop)
      character(len=max_length_string) :: value_kind                 ! kind of values: absval, absdev, reldev
   end type kernel_inverted_model
!
   contains
!----------------------------------------------------------------------------
!  Create kernel inverted model from solution vector of kernel linear system
!  prop:              ordered list of properties to be written (nprop)
!  values:            values of model (ncell,nprop)
!  vkind:             kind of values: absval, absdev, reldev
!
   subroutine createFromValuesKernelInvertedModel(this,prop,values,propsetname,vkind,errmsg)
      type (kernel_inverted_model) :: this
      character(len=char_len_par), dimension(:) :: prop
      real, dimension(:,:), pointer :: values
      character(len=*) :: propsetname, vkind
      type (error_message) :: errmsg
      character(len=35) :: myname = 'createFromValuesKernelInvertedModel'
   !
      if (size(values,2) .ne. size(prop)) then
         call add(errmsg,2,'Inconsistent dimension of model values and properties',myname)
         return
      end if
   !
      this%ncell = size(values,1)
      this%nprop = size(prop)
      this%model_values => values
      allocate(this%prop(this%nprop))
      this%prop = prop
      this%propsetname = propsetname
      this%value_kind = vkind
      this%model_values => values
   end subroutine createFromValuesKernelInvertedModel
!--------------------------------------------------------------------------
!  Create model perturbations from FMTOMO model and ASKI background model
!  for use with a SPECFEM forward run.
!  Assumes that Fmtomo model is either vp or vs but not both
!  Returns values for properties rho,vp,vs using correlation if specified
!
   subroutine createFromFmtomoAbmKernelInvertedModel(this,abm,fmtomo,invgrid,propset,thetac,phic,errmsg)
      type (kernel_inverted_model) :: this
      type (aski_background_model) :: abm
      type (fmtomo_model) :: fmtomo
      class (inversion_grid) :: invgrid
      type (property_set) :: propset
      double precision :: thetac,phic
      type (error_message) :: errmsg
      integer :: ip,icell,i,j
      double precision :: qmu,qkappa
      double precision, dimension(3) :: abmprop                             ! rho,vp,vs
      double precision :: xlc,ylc,zlc,xg,yg,zg,r,delta,xi,lat,vel,reldev
      double precision, dimension(3) :: cc
      double precision, dimension(:,:), pointer :: corrmat
      character(len=38) :: myname = 'createFromFmtomoAbmKernelInvertedModel'
   !
      this%ncell = invgrid%getNcellLocal()
      this%propsetname = .setname.propset
      this%value_kind = 'absdev'
      if (.not. equalString(this%propsetname,'isoVelocitySI')) then
         call add(errmsg,2,'We assume here that property set is isoVelocitySI',myname)
         return
      end if
   !
   !  prepare property correlation
   !
      call getRequiredKernelsByNamePropertySet(propset,this%prop)
      this%nprop = size(this%prop)
      corrmat => getCorrmatPropertySet(propset)
      ip = 2
      if (getVeltypeFmtomoModel(fmtomo) == 'vs') ip = 3
      allocate(this%model_values(this%ncell,this%nprop))
   !
   !  Go through inversion grid cells
   !
      do while (invgrid%nextCell(cc,icell))
      !
      !  transform cell center to global spherical coordinates
      !
         call coordinatesLCfromRCAxesRotation(mc_pid/2.,cc(1),cc(2),abm%rearth+cc(3),xlc,ylc,zlc)
         call coordinatesGCfromLCAxesRotation(thetac,phic,xlc,ylc,zlc,xg,yg,zg)
         call coordinatesLSfromLCAxesRotation(xg,yg,zg,r,delta,xi)
         lat = 0.5d0*mc_pid-delta
      !
      !  evaluate models there
      !
         call evalASKIBackgroundModel(abm,abm%rearth-r,abmprop(1),abmprop(2),abmprop(3),qmu,qkappa,errmsg)
         if (.level.errmsg == 2) return
         vel = evalFmtomoModel(fmtomo,r,lat,xi,errmsg)
         if (.level.errmsg == 2) return
         if (vel < 0.d0) then                       ! cell is outside Fmtonmo grid, assume background model there
            reldev = 0.d0
         else
            reldev = vel/abmprop(2)-1.d0
            if (ip == 3) reldev = vel/abmprop(3)-1.d0
         end if
      !
      !  apply property correlation
      !
         do i = 1,this%nprop
            j = getIndexFromNamePropertySet(propset,this%prop(i))
            this%model_values(icell,i) = abmprop(j)*corrmat(j,ip)*reldev
         end do
      end do
   end subroutine createFromFmtomoAbmKernelInvertedModel
!--------------------------------------------------------------------------
!  Create kernel inverted model from ASKI background model
!  to obtain 1D reference model values.
!  Returns values for properties rho,vp,vs
!
   subroutine createFromAbmKernelInvertedModel(this,abm,invgrid,thetac,phic,errmsg)
      type (kernel_inverted_model) :: this
      type (aski_background_model) :: abm
      class (inversion_grid) :: invgrid
      double precision :: thetac,phic
      type (error_message) :: errmsg
      integer :: icell
      double precision :: rho,vp,vs,qmu,qkappa
      double precision :: xlc,ylc,zlc,xg,yg,zg,r,delta,xi
      double precision, dimension(3) :: cc
      character(len=38) :: myname = 'createFromAbmKernelInvertedModel'
   !
      this%ncell = invgrid%getNcellLocal()
      this%propsetname = 'isoVelocitySI'
      this%value_kind = 'absval'
      this%nprop = 3
      allocate(this%prop(this%nprop))
      this%prop = ['rho','vp ','vs ']
      allocate(this%model_values(this%ncell,this%nprop))
   !
   !  Go through inversion grid cells
   !
      do while (invgrid%nextCell(cc,icell))
      !
      !  transform cell center to global spherical coordinates
      !
         call coordinatesLCfromRCAxesRotation(mc_pid/2.,cc(1),cc(2),abm%rearth+cc(3),xlc,ylc,zlc)
         call coordinatesGCfromLCAxesRotation(thetac,phic,xlc,ylc,zlc,xg,yg,zg)
         call coordinatesLSfromLCAxesRotation(xg,yg,zg,r,delta,xi)
      !
      !  evaluate models there
      !
         call evalASKIBackgroundModel(abm,abm%rearth-r,rho,vp,vs,qmu,qkappa,errmsg)
         if (.level.errmsg == 2) return
      !
      !  apply property correlation
      !
         this%model_values(icell,1) = rho
         this%model_values(icell,2) = vp
         this%model_values(icell,3) = rho
      end do
   end subroutine createFromAbmKernelInvertedModel
!-------------------------------------------------------------------------
!  Deallocate kernel inverted model
!
   subroutine deallocKernelInvertedModel(this)
      type (kernel_inverted_model) :: this
      if (allocated(this%prop)) deallocate(this%prop)
      if (associated(this%model_values)) deallocate(this%model_values)
   end subroutine deallocKernelInvertedModel
!----------------------------------------------------------------------------
!  Write kernel inverted model to HDF file
!  Cell information is needed for later reading by SPECFEM's generate_databases
!
   subroutine writeHDFKernelInvertedModel(this,filename,invgrid,errmsg)
      type (kernel_inverted_model) :: this
      character(len=*) :: filename
      class (inversion_grid) :: invgrid
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      double precision, dimension(:,:), pointer :: cc
      real, dimension(:), pointer :: rad
      real, dimension(:,:), pointer :: rcc
      integer, dimension(:,:), pointer :: nb
      integer(kind=8) :: fid
      integer :: i
      character(len=max_length_string) :: propstring
   !
      call createFileHDFWrapper(trim(filename),fid,errmsg)
      if (.level.errmsg == 2) return
      call writeStringAttributeHDFWrapper(fid,'property_set_name',trim(this%propsetname),errmsg)
      if (.level.errmsg == 2) return
      propstring = trim(this%prop(1))
      do i = 2,this%nprop
         propstring = trim(propstring)//':'//trim(this%prop(i))
      end do
      call writeStringAttributeHDFWrapper(fid,'properties',trim(propstring),errmsg)
      if (.level.errmsg == 2) return
   !
   !  kind of values (absval, absdev, reldev)
   !
      call writeStringAttributeHDFWrapper(fid,'kind_values',trim(this%value_kind),errmsg)
      if (.level.errmsg == 2) return
   !
   !  absolute model values
   !
      call arra%assoc2d(this%model_values)
      call writeArrayHDFWrapper(fid,'model_values',arra,errmsg,xferprp = hdf_xferprp)
      call arra%deassoc()
      if (.level.errmsg == 2) return
   !
   !  cell centers
   !
      call invgrid%getCellCenters(cc)
      allocate(rcc(3,this%ncell))
      rcc = real(cc)
      call arra%assoc2d(rcc)
      call writeArrayHDFWrapper(fid,'cell_centers',arra,errmsg,xferprp = hdf_xferprp)
      call arra%deassoc()
      deallocate(rcc); nullify(cc)
      if (.level.errmsg == 2) return
   !
   !  cell radii
   !
      call invgrid%getRealCellRadii(rad)
      call arra%assoc1d(rad)
      call writeArrayHDFWrapper(fid,'cell_radii',arra,errmsg,xferprp = hdf_xferprp)
      call arra%deassoc()
      deallocate(rad)
   !
   !  cell face neighbours
   !
      call invgrid%getFaceNeighbours(nb)
      call aria%assoc2d(nb)
      call writeArrayHDFWrapper(fid,'neighbours',aria,errmsg,xferprp = hdf_xferprp)
      call aria%deassoc()
      nullify(nb)

      call closeFileHDFWrapper(fid,errmsg)
      if (.level.errmsg == 2) return
   end subroutine writeHDFKernelInvertedModel
!---------------------------------------------------------------------------
!  Read kernel inverted model
!  Not suitable for reading a kim-file by SPECFEM
!
   subroutine readHDFKernelInvertedModel(this,filename,errmsg)
      type (kernel_inverted_model) :: this
      character(len=*) :: filename
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      integer(kind=8) :: fid
      integer :: i,slen
      character(len=max_length_string), dimension(:), pointer :: res
      character(len=max_length_string) :: propstring,cval
      character(len=26) :: myname = 'readHDFKernelInvertedModel'
   !
      call openFileRoHDFWrapper(filename,fid,errmsg)
      if (.level.errmsg == 2) return
      call readStringAttributeHDFWrapper(fid,'property_set_name',cval,slen,errmsg)
      if (.level.errmsg == 2) return
      if (slen > char_len_pmtrz) then
         call add(errmsg,2,'returned string length for property set name greater than assumed length',myname)
         return
      end if
      this%propsetname = cval(1:slen)
   !
      call readStringAttributeHDFWrapper(fid,'properties',propstring,slen,errmsg)
      if (.level.errmsg == 2) return
      res => getWordsString(propstring(1:slen),':',errmsg)
      if (.level.errmsg == 2) return
      this%nprop = size(res)
      allocate(this%prop(this%nprop))
      do i = 1,this%nprop
         this%prop(i) = trim(res(i))
      end do
      deallocate(res)
   !
      call readStringAttributeHDFWrapper(fid,'kind_values',cval,slen,errmsg)
      if (.level.errmsg == 2) return
      this%value_kind = cval(1:slen)
   !
      call readArrayHDFWrapper(fid,'model_values',arra,errmsg,xferprp = hdf_xferprp)
      if (.level.errmsg == 2) return
      this%model_values => arra%get2d()
      call arra%deassoc()
      this%ncell = size(this%model_values,1)
   !
      call closeFileHDFWrapper(fid,errmsg)
      if (.level.errmsg == 2) return
   end subroutine readHDFKernelInvertedModel
!------------------------------------------------------------------------
!  Iterate over properties of a given kernel inverted model
!
   logical function nextPropKernelInvertedModel(this,prop,iprop)
      type (kernel_inverted_model) :: this
      character(len=char_len_par) :: prop
      integer, optional :: iprop
      integer :: call_count = 0
      save :: call_count
   !
      call_count = call_count+1
      if (call_count <= this%nprop) then
         prop = this%prop(call_count)
         if (present(iprop)) iprop = call_count
         nextPropKernelInvertedModel = .true.
      else
         call_count = 0
         prop = ''
         if (present(iprop)) iprop = 0
         nextPropKernelInvertedModel = .false.
      endif
   end function nextPropKernelInvertedModel
!---------------------------------------------------------------
!  get pointer to absolute model values of given property
!
   function getValuesPropKernelInvertedModel(this,j) result(res)
      type (kernel_inverted_model), intent(in), target :: this
      integer :: j
      real, dimension(:), pointer :: res
      res => this%model_values(:,j)
   end function getValuesPropKernelInvertedModel
!---------------------------------------------------------------
!  return a double precision flattened array of model values
!
   function doubleFlattenedValuesKernelInvertedModel(this) result(res)
      type (kernel_inverted_model), intent(in) :: this
      double precision, dimension(:), pointer :: res
      integer :: j,ip,ic
   !
      allocate(res(this%nprop*this%ncell))
      j = 1
      do ip = 1,this%nprop
         do ic = 1,this%ncell
            res(j) = this%model_values(ic,ip)
            j = j+1
         end do
      enddo
   end function doubleFlattenedValuesKernelInvertedModel

end module kernelInvertedModel
