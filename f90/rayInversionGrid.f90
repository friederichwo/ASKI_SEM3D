!----------------------------------------------------------------------------
!   Copyright 2023 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   Inversion grid based on FM3D propagation grid
!   Elementary cells of the propagation grid are used as inversion grid cells
!   Cell center is set at grid point with plus/minus 1/2 grid spacing around.
!   With this definition the number of cells is equal to the mumber of grid points.
!   There is only one wavefield point in a cell which is the cell center.
!
module rayInversionGrid
!
   use inversionGrid
   use inputParameter
   use mathConstants
   use smartUtils
   use axesRotation
   use errorMessage
   use string
   use globalMpiInfo
   implicit none
!
   type, extends(inversion_grid) :: ray_inversion_grid
      private
      integer :: nr,nlat,nlon                                        ! number of grid points along each dimension
      double precision :: dr,dlat,dlon                               ! grid spacing along each dimension (m,rad,rad)
      double precision :: rmin,latmin,lonmin                         ! minimum value for each space coordinate (m,rad,rad)
      double precision :: rearth                                     ! earth radius in m
      double precision, dimension(:,:), allocatable :: model_values  ! model values on grid (ncell,3)
   contains
      procedure :: create => createRayInversionGrid
      procedure :: dealloc => deallocateRayInversionGrid
      procedure :: getModelValues => getModelValuesOnRayInversionGrid
      procedure :: getRealCellVolumes => getRealCellVolumesRayInversionGrid
      procedure :: getWavefieldPoints => getWavefieldPointsRayInversionGrid
      procedure :: getWpModelValues => getWpModelValuesRayInversionGrid
      procedure :: getIntegrationWeights => getIntegrationWeightsRayInversionGrid
      procedure :: getGlobalDimensions => getGlobalDimensionsRayInversionGrid
      procedure :: getGridCenter => getGridCenterRayInversionGrid
      procedure :: getEarthRadius => getEarthRadiusRayInversionGrid
   end type ray_inversion_grid
!
   contains
!------------------------------------------------------------------------
!  Create ray inversion grid
!  This routine should always be called from iteration step basics
!  path:           path to propgrid file
!  errmsg:         error message
!  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
!
   subroutine createRayInversionGrid(this,iter_path,do_hyperslab,inpar_inv,errmsg)
      class (ray_inversion_grid) :: this
      character(len=*) :: iter_path
      logical :: do_hyperslab
      type (input_parameter) :: inpar_inv
      type (error_message) :: errmsg
      ! local
      integer :: icell,i,j,k,ios,ic,l
      double precision :: thetac,phic,rearth,xg,yg,zg,xlc,ylc,zlc,xr,yr,zr
      double precision :: cc(3),cor(3,8)
      character(len=22) :: myname = 'createRayInversionGrid'
   !----------------------------------------------------------------------------------
   !  read propgrid file (units are km and degrees), better HDF with hyperslab reading
   !
      open(unit=1,file=trim(iter_path//"propgrid.in"),status='unknown',form='formatted',action='read',iostat=ios)
      if (ios /= 0) goto 1
      read(1,*,iostat=ios) this%nr,this%nlat,this%nlon
      if (ios /= 0) goto 1
      read(1,*,iostat=ios) this%dr,this%dlat,this%dlon
      if (ios /= 0) goto 1
      read(1,*,iostat=ios) zmin,this%latmin,this%lonmin              ! surface at z=0
      if (ios /= 0) goto 1
   !
   !  also read model values on grid here
   !
      close(1)
   !
   !  convert to SI and radians
   !
      rearth = inpar_inv.dval.'REARTH'
      this%dr = this%dr*1.d3
      this%rmin = rearth-zmin*1.d3
      this%dlat = this%dlat*mc_deg2radd
      this%dlon = this%dlon*mc_deg2radd
      this%latmin = this%latmin*mc_deg2radd
      this%lonmin = this%lonmin*mc_deg2radd
      thetac = inpar_inv.dval.'INVGRID_CENTER_LAT'
      thetac = 0.5*mc_pid-thetac*mc_deg2radd
      phic = (inpar_inv.dval.'INVGRID_CENTER_LON')*mc_deg2radd
      this%rearth = rearth
   !
   !  cell center is only wavefield point, hence ncell = nwp and ngll = 1
   !
      this%ngll = 1
      this%ncell_all = this%nr*this%nlat*this%nlon
      this%nwp_all = this%ncell_all
   !
   !  distribute cells among processes (used for kernel computation)
   !  equal sharing of invgrid cells, remainder is distributed again starting from rank zero
   !  regular share: nreg = ncell/numtasks
   !  for 0 < myrank < nremain = mod(ncell,numtasks) take one more
   !
      if (do_hyperslab) then
         call shareLinearArraySmartUtils(this%ncell_all,numtasks,myrank,this%ncell_local,this%offcell)
         this%nwp_local = this%ncell_local
      else
         this%ncell_local = this%ncell_all
         this%nwp_local = this%nwp_all
         this%offcell = 0
      end if
   !
   !  order of cells is first along longitude, then along increasing latitude and then along increasing radius
   !  given index of grid point i,j,k the global cell index is ic = k+(j-1)*nlon+(i-1)*nlon*nlat
   !  given index of cell ic, then (with r < 1):
   !  ic/(nlat*nlon) = i-1+r => i = int(ic/(nlat*nlon))+1
   !  nremain = ic-(i-1)*nlat*nlon = k+(j-1)*nlon
   !  nremain/nlon = j-1+r => j = int(nremain/nlon)+1
   !  nremain-(j-1)*nlon = k
   !
      allocate(this%cell_center(3,this%ncell_local))
      allocate(this%cell_corner(3,8,this%ncell_local))
      allocate(this%face_neighbour(7,this%ncell_local))
      do icell = 1,this%ncell_local
         ic = this%offcell+icell
         call gridIndicesFromCellIndexSmartUtils(ic,this%nlon,this%nlat,i,j,k)
         cc(1) = this%rmin+(i-1)*this%dr
         cc(2) = this%latmin+(j-1)*this%dlat
         cc(3) = this%lonmin+(k-1)*this%dlon
      !
      !  convert from global spherical to box centered Cartesian coordinates
      !  (thetac,phic):: position of box center
      !
         call coordinatesLCfromLSAxesRotation(cc(1),0.5d0*mc_pid-cc(2),cc(3),xg,yg,zg)
         call coordinatesLCfromGCAxesRotation(thetac,phic,xg,yg,zg,xlc,ylc,zlc)
         call coordinatesRCfromLCAxesRotation(mc_pid/2.,xlc,ylc,zlc,xr,yr,zr)
         zr = zr-rearth
         this%cell_center(:,icell) = [xr,yr,zr]
      !
      !  compute cell corners
      !
         cor(1,1:4) = cc(1)-0.5d0*this%dr              ! corners at r = rc-0.5*dr
         cor(1,5:8) = cc(1)+0.5d0*this%dr              ! corners at r = rc+0.5*dr
         cor(2,[1,2,5,6]) = cc(2)-0.5d0*this%dlat      ! corners in the south
         cor(2,[3,4,7,8]) = cc(2)+0.5d0*this%dlat      ! corners in the north
         cor(3,[1,4,5,8]) = cc(3)-0.5d0*this%dlon      ! corners in the west
         cor(3,[2,3,6,7]) = cc(3)+0.5d0*this%dlon      ! corners in the east
      !
      !  convert corners from global spherical to box centered Cartesian coordinates
      !  (thetac,phic):: position of box center
      !
         do l = 1,8
            call coordinatesLCfromLSAxesRotation(cor(1,l),0.5d0*mc_pid-cor(2,l),cor(3,l),xg,yg,zg)
            call coordinatesLCfromGCAxesRotation(thetac,phic,xg,yg,zg,xlc,ylc,zlc)
            call coordinatesRCfromLCAxesRotation(0.5d0*mc_pid,xlc,ylc,zlc,xr,yr,zr)
            zr = zr-rearth
            this%cell_corner(:,l,icell) = [xr,yr,zr]
         end do
      !
      !  determine number and indices of face neighbours
      !  first treat cells in corners of grid with 3 neighbours
      !
         if (i == 1 .and. j == 1 .and. k == 1) then                                             ! DSW corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (i == this%nr .and. j == 1 .and. k == 1) then                                  ! USW corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
         else if (i == 1 .and. j == this%nlat .and. k == 1) then                                ! DNW corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (i == this%nr .and. j == this%nlat .and. k == 1) then                          ! UNW corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
         else if (i == 1 .and. j == 1 .and. k == this%nlon) then                                ! DSE corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (i == this%nr .and. j == 1 .and. k == this%nlon) then                          ! USE corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
         else if (i == 1 .and. j == this%nlat .and. k == this%nlon) then                        ! DNE corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (i == this%nr .and. j == this%nlat .and. k == this%nlon) then                  ! UNE corner
            this%face_neighbour(1,icell) = 3
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
      !
      !  now treat cells on faces of grid away from corners with 5 neighbours
      !
         else if (i == 1 .and. j > 1 .and. j < this%nlat .and. k > 1 .and. k < this%nlon) then       ! bottom face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (i == this%nr .and. j > 1 .and. j < this%nlat .and. k > 1 .and. k < this%nlon) then  ! top face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
         else if (j == 1 .and. i > 1 .and. i < this%nr .and. k > 1 .and. k < this%nlon) then          ! south face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (j == this%nlat .and. i > 1 .and. i < this%nr .and. k > 1 .and. k < this%nlon) then  ! north face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (k == 1 .and. i > 1 .and. i < this%nr .and. j > 1 .and. j < this%nlat) then          ! west face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(6,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         else if (k == this%nlon .and. i > 1 .and. i < this%nr .and. j > 1 .and. j < this%nlat) then   ! east face
            this%face_neighbour(1,icell) = 5
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(6,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
      !
      !  treat all the inner cells with 6 neighbours
      !
         else                                                                                   ! inner cells
            this%face_neighbour(1,icell) = 6
            this%face_neighbour(2,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k-1,this%nlon,this%nlat)
            this%face_neighbour(3,icell) = cellIndexFromGridIndicesSmartUtils(i,j,k+1,this%nlon,this%nlat)
            this%face_neighbour(4,icell) = cellIndexFromGridIndicesSmartUtils(i,j-1,k,this%nlon,this%nlat)
            this%face_neighbour(5,icell) = cellIndexFromGridIndicesSmartUtils(i,j+1,k,this%nlon,this%nlat)
            this%face_neighbour(6,icell) = cellIndexFromGridIndicesSmartUtils(i-1,j,k,this%nlon,this%nlat)
            this%face_neighbour(7,icell) = cellIndexFromGridIndicesSmartUtils(i+1,j,k,this%nlon,this%nlat)
         endif
      end do
1     if (ios /= 0) then
         call add(errmsg,2,'problems reading fmtomo model file',myname)
      end if
   end subroutine createRayInversionGrid
!------------------------------------------------------------------------
!  deallocate inversion grid
!
   subroutine deallocateRayInversionGrid(this)
      class (ray_inversion_grid) :: this
      if (allocated(this%model_values)) deallocate(this%model_values)
      if (allocated(this%cell_corner)) deallocate(this%cell_corner)
      if (allocated(this%cell_center)) deallocate(this%cell_center)
      if (allocated(this%face_neighbour)) deallocate(this%face_neighbour)
   end subroutine deallocateRayInversionGrid
!------------------------------------------------------------------------
!  get global dimensions of inversion grid
!
   subroutine getGlobalDimensionsRayInversionGrid(this,n1,n2,n3)
      class (ray_inversion_grid) :: this
      integer :: n1,n2,n3
      n1 = this%nr
      n2 = this%nlat
      n3 = this%nlon
   end subroutine getGlobalDimensionsRayInversionGrid
!------------------------------------------------------------------------
!  get coordinates of center of inversion grid
!
   subroutine getGridCenterRayInversionGrid(this,lat,lon)
      class (ray_inversion_grid) :: this
      double precision :: lat,lon
      lat = 0.5d0*(this%latmin+(this%nlat-1)*this%dlat)
      lon = 0.5d0*(this%lonmin+(this%nlon-1)*this%dlon)
   end subroutine getGridCenterRayInversionGrid
!------------------------------------------------------------------------
!  get earth radius
!
   function getEarthRadiusRayInversionGrid(this) result(r)
      class (ray_inversion_grid) :: this
      double precision :: r
      r = this%rearth
   end function getEarthRadiusRayInversionGrid
!------------------------------------------------------------------------
!  compute volumes of cells as r^2 sin\theta dr d\theta d\phi
!
   subroutine getRealCellVolumesRayInversionGrid(this,res)
      class (ray_inversion_grid) :: this
      real, dimension(:), pointer :: res
      integer :: ioff,joff,koff,iend,jend,kend,ic,i,j,k
      double precision :: r,stheta
!
      allocate(res(this%ncell_local))
      call gridIndicesFromCellIndexSmartUtils(this%offcell,this%nlon,this%nlat,ioff,joff,koff)
      call gridIndicesFromCellIndexSmartUtils(this%offcell+this%ncell_local,this%nlon,this%nlat,iend,jend,kend)
!
      ic = 0
      do i = ioff+1,iend
         r = this%rmin+(i-1)*this%dr
         do j = joff+1,jend
            stheta = dcos(this%latmin+(j-1)*this%dlat)
            do k = koff+1,kend
               ic = ic+1
               res(ic) = r**2*stheta*this%dr*this%dlat*this%dlon
            end do
         end do
      end do
   end subroutine getRealCellVolumesRayInversionGrid
!------------------------------------------------------------------------
!  get pointers to wavefield points
!
   subroutine getWavefieldPointsRayInversionGrid(this,x,y,z)
      class (ray_inversion_grid), target :: this
      double precision, dimension(:), pointer :: x,y,z
      x => this%cell_center(1,:); y => this%cell_center(2,:); z => this%cell_center(3,:)
   end subroutine getWavefieldPointsRayInversionGrid
!------------------------------------------------------------------------
!  get pointers to model values on wavefield points
!
   subroutine getWpModelValuesRayInversionGrid(this,rho,vp,vs)
      class (ray_inversion_grid), target :: this
      double precision, dimension(:), pointer :: rho,vp,vs
      rho => this%model_values(:,1); vp => this%model_values(:,2); vs => this%model_values(:,3)
   end subroutine getWpModelValuesRayInversionGrid
!------------------------------------------------------------------------
!  Return model values on inversion grid
!  If in multi_process hyperslab mode, return property on my share of invgrid
!  Values should then be inserted into global grid using this%offcell
!  In not in hyperslab mode, values on full grid are returned
!
   function getModelValuesOnRayInversionGrid(this) result(res)
      class (ray_inversion_grid), target :: this
      double precision, dimension(:,:), pointer :: res
      res => this%model_values
   end function getModelValuesOnRayInversionGrid
!------------------------------------------------------------------------
!  get pointer to integration weights
!  weight:   on input:  sum of wavenumber vectors for cells (ncell,3)
!  weight:   on output: integration weights for cells (ncell,3)
!
   subroutine getIntegrationWeightsRayInversionGrid(this,weight)
      class (ray_inversion_grid), target :: this
      double precision, dimension(:,:), pointer :: weight
      integer :: i,j,k,ic
      double precision :: r,teta,x,y,z
   !
      do ic = 1,this%ncell_local
         call gridIndicesFromCellIndexSmartUtils(ic,this%nlon,this%nlat,i,j,k)
         r = this%rmin+(i-1)*this%dr
         teta = 0.5d0*mc_pid-(this%latmin+(j-1)*this%dlat)
         x = 0.5d0*weight(1,ic)*this%dr
         y = 0.5d0*weight(2,ic)*r*this%dlat
         z = 0.5d0*weight(3,ic)*r*dsin(teta)*this%dlon
         weight(1,ic) = sincSmartUtils(x)*this%dr
         weight(2,ic) = sincSmartUtils(y)*r*this%dlat
         weight(3,ic) = sincSmartUtils(z)*r*dsin(teta)*this%dlon
      end do
   end subroutine getIntegrationWeightsRayInversionGrid
!
end module rayInversionGrid
