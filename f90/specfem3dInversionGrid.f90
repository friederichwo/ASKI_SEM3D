!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!   inversion grid consisting of SPECFEM3D elements as inversion grid cells
!
module specfem3dInversionGrid
!
   use mathConstants
   use inversionGrid
   use axesRotation
   use specfem3dForASKIFiles
   use heapSort
   use locatePoint
   use errorMessage
   use inputParameter
   use string
   implicit none
!
   type, extends(inversion_grid) ::  specfem3d_inversion_grid
      private
      ! Global grid dimensions
      integer :: nx,ny,nz                                           ! number of elements in x,y,z direction
      ! Geographical coordinates of grid center and earth radius
      double precision :: lat_grid_center,lon_grid_center           ! lat and lon in radians
      double precision :: rearth                                    ! radius im m
      ! Wavefield points
      double precision, dimension(:), allocatable :: x              ! x coordinates of points (proc-specific)
      double precision, dimension(:), allocatable :: y              ! y coordinates of points (proc-specific)
      double precision, dimension(:), allocatable :: z              ! z coordinates of points (proc-specific)
      ! Model on wavefield points
      double precision, dimension(:), allocatable :: rho            ! density at points (proc-specific)
      double precision, dimension(:), allocatable :: vp             ! vp at points (proc-specific)
      double precision, dimension(:), allocatable :: vs             ! vs at points (proc-specific)
      ! GLL points
      integer :: ngllxyz                                            ! number of wavefield points per dimension (all equal)
      double precision, dimension(:), allocatable :: jacobian       ! jacobian values at wp points (proc-specific)
      ! Cells
      double precision, dimension(:,:), allocatable :: integration_weight   ! integration weights (ngll,ncell_local)
      double precision, dimension(:,:), allocatable :: model_values         ! model values on grid (ncell,3)
      ! for finding the cell index from SPECFEM coordinates
      integer, dimension(:), allocatable :: sortidx                 ! index array relating true depth sorted to original cells
      double precision, dimension(:), allocatable :: zlevel         ! true negative depth values of bottom of cell layers
   !
   contains
      procedure :: create => createSpecfem3dInversionGrid
      procedure :: dealloc => deallocateSpecfem3dInversionGrid
      procedure :: getIntegrationWeights => getIntegrationWeightsSpecfem3dInversionGrid
      procedure :: getModelValues => getModelValuesOnSpecfem3dInversionGrid
      procedure :: getRealCellVolumes => getRealCellVolumesSpecfem3dInversionGrid
      procedure :: getWavefieldPoints => getWavefieldPointsSpecfem3dInversionGrid
      procedure :: getWpModelValues => getWpModelValuesSpecfem3dInversionGrid
      procedure :: getGlobalDimensions => getGlobalDimensionsSpecfem3dInversionGrid
      procedure :: getGridCenter => getGridCenterSpecfem3dInversionGrid
      procedure :: getEarthRadius => getEarthRadiusSpecfem3dInversionGrid
      procedure :: convertSphericalCoordinates => convertSphericalCoordinatesSpecfem3dInversionGrid
      procedure :: createLocator => createLocatorSpecfem3dInversionGrid
      procedure :: findCell => findCellSpecfem3dInversionGrid
      procedure :: isBelow => isBelowSpecfem3dInversionGrid
   end type specfem3d_inversion_grid
!
!  GLL weights  for GLL quadrature of orders 5 taken from thesis:
!   Schuberth, B., "The Spectral Element Method for Seismic Wave Propagation - Theory,
!   Implementation and Comparison to Finite Difference Methods", Diploma Thesis, Dept.
!   for Earth and Environmental Sciences, Ludwig-Maximilians-Universität München, (2003)
!
   double precision, dimension(5) :: wgll = (/ 0.100000000, 0.544444444, 0.711111111, 0.544444444, 0.100000000 /)
!
   contains
!------------------------------------------------------------------------
!  Create specfem3d inversion grid
!  This routine should always be called from iteration step basics
!  path:       path to ASKI main file
!  errmsg:         error message
!  do_hyperslab:   read inversion grid wavefield point related arrays in hyperslab mode
!
   subroutine createSpecfem3dInversionGrid(this,iter_path,do_hyperslab,inpar_inv,errmsg)
      class (specfem3d_inversion_grid) :: this
      character(len=*) :: iter_path
      logical :: do_hyperslab
      type (input_parameter) :: inpar_inv
      type (error_message) :: errmsg
      ! local
      integer, dimension(:), allocatable :: id
      integer :: NGLLX,NGLLY,NGLLZ
      integer :: specfem_version,type_invgrid
      integer :: icell,i,j,k,ip,ishift,nf,nproc
      integer, dimension(8) :: p_shift
      integer, dimension(:), allocatable  :: jf
      double precision :: df
      character(len=max_length_string) :: path
      character(len=28) :: myname = 'createSpecfem3dInversionGrid'
   !---------------------------------------------------------------------------------
   !  read ASKI main file
   !
      path = trim(iter_path)//trim(inpar_inv.sval.'PATH_ASKI_MAIN_FILES')
      call readSpecfem3dForASKIMainFileHDF(trim(path)//"ASKI.main",&
           specfem_version,nproc,type_invgrid,this%nwp_all,this%nwp_local,this%offwp,&
           df,nf,jf,NGLLX,NGLLY,NGLLZ,this%ncell_all,this%face_neighbour,&
           do_hyperslab,errmsg,this%x,this%y,this%z,this%rho,this%vp,this%vs,this%jacobian)
      if(.level.errmsg == 2) return
      deallocate(jf)
      this%ngll = NGLLX*NGLLY*NGLLZ
      this%ngllxyz = NGLLX
   !
   !  set global grid dimensions
   !
      id = ivec(inpar_inv,'INVGRID_DIMENSIONS',3)
      this%nx = id(1); this%ny = id(2); this%nz = id(3)
      deallocate(id)
   !
   !  set geographical coordinates of grid center and earth radius
   !
      this%lat_grid_center = (inpar_inv.dval.'INVGRID_CENTER_LAT')*mc_deg2radd
      this%lon_grid_center = (inpar_inv.dval.'INVGRID_CENTER_LON')*mc_deg2radd
      this%rearth = inpar_inv.dval.'REARTH'
   !
   ! check if ASKI main file fits to this special code
   !
      if (specfem_version /= 2) call add(errmsg,2,"use SPECFEM3D cartesian with this code",myname)
      if (type_invgrid /= 4) call add(errmsg,2,"invalid inversion grid type, should be 4",myname)
      if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) call add(errmsg,2,"NGLLX/Y/Z should be 5", myname)
      if (this%ncell_all /= this%nwp_all/this%ngll) then
         call add(errmsg,2,"wavefield points, GLL-points and cells are inconsistent", myname)
      endif
      if(.level.errmsg == 2) return
   !
   !  cells are counted according to the ordering of the wavefield points in the ASKI main file
   !  calculate number of proc-specific cells and associated offset in total cell count
   !
      this%ncell_local = this%nwp_local/this%ngll
      this%offcell = this%offwp/this%ngll
   !
   ! define cell center: for each dimension (i.e. x,y,z) use the center gll point (case ngll = 5, k=63)
   !
      allocate(this%cell_center(3,this%ncell_local))
      k = NGLLX*NGLLY*(NGLLZ-1)/2+NGLLX*(NGLLY-1)/2+(NGLLX+1)/2
      do icell = 1,this%ncell_local
         ishift = (icell-1)*this%ngll    ! compute index in arrays x,y,z, where this element starts (minus 1)
         this%cell_center(1,icell) = this%x(k+ishift)
         this%cell_center(2,icell) = this%y(k+ishift)
         this%cell_center(3,icell) = this%z(k+ishift)
      enddo
   !
   !  compute cell corners
   !
      allocate(this%cell_corner(3,8,this%ncell_local))
      p_shift(1) = 1
      p_shift(2) = NGLLX
      p_shift(3) = NGLLX*NGLLY
      p_shift(4) = (NGLLY-1)*NGLLX + 1
      do j = 1,4
         p_shift(j+4) = (NGLLZ-1)*NGLLX*NGLLY + p_shift(j)      ! above point j
      end do
      do icell = 1,this%ncell_local
         ishift = (icell-1)*this%ngll          ! compute index in arrays x,y,z, where this element starts (minus 1)
         do j = 1,8
            this%cell_corner(1,j,icell) = this%x(ishift+p_shift(j))
            this%cell_corner(2,j,icell) = this%y(ishift+p_shift(j))
            this%cell_corner(3,j,icell) = this%z(ishift+p_shift(j))
         enddo
      enddo
   !
   ! compute integration weights (cell GLL-points in IJK-order)
   !
      allocate(this%integration_weight(this%ngll,this%ncell_local))
      do icell = 1,this%ncell_local
         ishift = (icell-1)*this%ngll
         ip = 0
         do k = 1,NGLLZ
            do j = 1,NGLLY
               do i = 1,NGLLX
                  ip = ip + 1
                  this%integration_weight(ip,icell) = wgll(i)*wgll(j)*wgll(k)*this%jacobian(ishift+ip)
               enddo
            enddo
         enddo
      enddo
   !
   !  compute model values at cell centers of inversion grid
   !
      allocate(this%model_values(this%ncell_local,3))
      k = NGLLX*NGLLY*(NGLLZ-1)/2+NGLLX*(NGLLY-1)/2+(NGLLX+1)/2
      do j = 1,this%ncell_local
         this%model_values(j,1) = this%rho(k)
         this%model_values(j,2) = this%vp(k)
         this%model_values(j,3) = this%vs(k)
         k = k+this%ngll
      enddo
   end subroutine createSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  deallocate simple layered Cartesian inversion grid
!
   subroutine deallocateSpecfem3dInversionGrid(this)
      class (specfem3d_inversion_grid) :: this
      if (allocated(this%x)) deallocate(this%x)
      if (allocated(this%y)) deallocate(this%y)
      if (allocated(this%z)) deallocate(this%z)
      if (allocated(this%rho)) deallocate(this%rho)
      if (allocated(this%vp)) deallocate(this%vp)
      if (allocated(this%vs)) deallocate(this%vs)
      if (allocated(this%jacobian)) deallocate(this%jacobian)
      if (allocated(this%model_values)) deallocate(this%model_values)
      if (allocated(this%cell_corner)) deallocate(this%cell_corner)
      if (allocated(this%cell_center)) deallocate(this%cell_center)
      if (allocated(this%face_neighbour)) deallocate(this%face_neighbour)
      if (allocated(this%integration_weight)) deallocate(this%integration_weight)
      if (allocated(this%zlevel)) deallocate(this%zlevel)
      if (allocated(this%sortidx)) deallocate(this%sortidx)
   end subroutine deallocateSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get global dimensions of inversion grid
!
   subroutine getGlobalDimensionsSpecfem3dInversionGrid(this,n1,n2,n3)
      class (specfem3d_inversion_grid) :: this
      integer :: n1,n2,n3
      n1 = this%nx
      n2 = this%ny
      n3 = this%nz
   end subroutine getGlobalDimensionsSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get coordinates of center of inversion grid
!
   subroutine getGridCenterSpecfem3dInversionGrid(this,lat,lon)
      class (specfem3d_inversion_grid) :: this
      double precision :: lat,lon
      lat = this%lat_grid_center
      lon = this%lon_grid_center
   end subroutine getGridCenterSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get earth radius
!
   function getEarthRadiusSpecfem3dInversionGrid(this) result(r)
      class (specfem3d_inversion_grid) :: this
      double precision :: r
      r = this%rearth
   end function getEarthRadiusSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get pointer to integration weights
!
   subroutine getIntegrationWeightsSpecfem3dInversionGrid(this,weight)
      class (specfem3d_inversion_grid), target :: this
      double precision, dimension(:,:), pointer :: weight
      weight => this%integration_weight
   end subroutine getIntegrationWeightsSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  Return model values on inversion grid
!  If in multi_process hyperslab mode, return property on my share of invgrid
!  Values should then be inserted into global grid using this%offcell
!  In not in hyperslab mode, values on full grid are returned
!
   function getModelValuesOnSpecfem3dInversionGrid(this) result(res)
      class (specfem3d_inversion_grid), target :: this
      double precision, dimension(:,:), pointer :: res
      res => this%model_values
   end function getModelValuesOnSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get volume of inversion grid cells
!
   subroutine getRealCellVolumesSpecfem3dInversionGrid(this,res)
      class (specfem3d_inversion_grid) :: this
      real, dimension(:), pointer :: res
      integer :: ix,iy,iz,ip,icell
   !
      allocate(res(this%ncell_local))
      res = 0.
      ! set start index for accessing this%jacobian to last point index of previous cell
      ip = 0
      do icell = 1,this%ncell_local
         do iz = 1,this%ngllxyz
            do iy = 1,this%ngllxyz
               do ix = 1,this%ngllxyz
                  ip = ip + 1
                  res(icell) = res(icell) + wgll(ix)*wgll(iy)*wgll(iz)*this%jacobian(ip)
               end do
            end do
         end do
      end do
   end subroutine getRealCellVolumesSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get pointers to wavefield points
!
   subroutine getWavefieldPointsSpecfem3dInversionGrid(this,x,y,z)
      class (specfem3d_inversion_grid), target :: this
      double precision, dimension(:), pointer :: x,y,z
      x => this%x; y => this%y; z => this%z
   end subroutine getWavefieldPointsSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  get pointers to model values on wavefield points
!
   subroutine getWpModelValuesSpecfem3dInversionGrid(this,rho,vp,vs)
      class (specfem3d_inversion_grid), target :: this
      double precision, dimension(:), pointer :: rho,vp,vs
      rho => this%rho; vp => this%vp; vs => this%vs
   end subroutine getWpModelValuesSpecfem3dInversionGrid
!-----------------------------------------------------------------------
!  convert spherical coordinates to SPECFEM coordinates
!
   subroutine convertSphericalCoordinatesSpecfem3dInversionGrid(this,r,lat,lon,x,y,z)
      class (specfem3d_inversion_grid) :: this
      double precision :: r,lat,lon,x,y,z
      double precision :: xg,yg,zg,xr,yr,zr
   !
      call coordinatesLCfromLSAxesRotation(r,0.5*mc_pid-lat,lon,xg,yg,zg)
      call coordinatesLCfromGCAxesRotation(0.5*mc_pid-this%lat_grid_center,this%lon_grid_center,xg,yg,zg,xr,yr,zr)
      call coordinatesRCfromLCAxesRotation(0.5*mc_pid,xr,yr,zr,x,y,z)
      z = z-this%rearth
   end subroutine convertSphericalCoordinatesSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  Create specfem3d inversion grid locator
!    invgrid:     inversion grid
!    rearth:      earth radius in meters
!    nx,ny,nz:    number of cells per dimension
!    errmsg:      error message
!
   subroutine createLocatorSpecfem3dInversionGrid(this,errmsg)
      class (specfem3d_inversion_grid) :: this
      type (error_message) :: errmsg
      double precision, dimension(:,:), allocatable :: co_sorted
      integer :: i,ncell,nh
      double precision :: rg,rearth
   !
      rearth = this%rearth
      ncell = this%ncell_all
      nh = this%nx*this%ny
   !
   !  replace z-coordinate of cell corner by true (negative) depth
   !  using rg = sqrt(x**2+y**2+(z+rearth)**2)
   !
      allocate(co_sorted(3,ncell))
      co_sorted = this%cell_corner(:,1,:)
      do i = 1,ncell
         rg = hypot( hypot(co_sorted(1,i),co_sorted(2,i)), co_sorted(3,i)+rearth)
         co_sorted(3,i) = rg-rearth
      enddo
   !
   !  sort co_sorted(3,ncell) according to values of true negative depth, overwrites co_sorted
   !
      allocate(this%sortidx(ncell))
      call heapSort2D(co_sorted,this%sortidx,dimsort = 2,is = 3)
   !
   !  set (negative) depth levels of inversion grid
   !
      allocate(this%zlevel(this%nz))
      do i = 1,this%nz
         this%zlevel(i) = co_sorted(3,1+(i-1)*nh)
         print *,'zlevel: ',i,this%zlevel(i)
      end do
      deallocate(co_sorted)
   end subroutine createLocatorSpecfem3dInversionGrid
!------------------------------------------------------------------------
!  Find cell containing a point given in SPECFEM coordinates
!  x, y, z:     SPECFEM coordinates of point
!  icell:       (out) index of cell
!
   subroutine findCellSpecfem3dInversionGrid(this,x,y,z,icell,errmsg)
      class (specfem3d_inversion_grid), target :: this
      double precision :: x,y,z
      integer :: icell
      type (error_message) :: errmsg
      integer, dimension(:), pointer :: idx
      double precision :: r
      integer :: ic1,ic2,ic,lev,nh
      character(len=37) :: myname = 'findCellSpecfem3dInversionGrid'
   !
      if (.not. allocated(this%sortidx)) then
         call add(errmsg,2,'cell locator was not initiated',myname)
         return
      end if
   !
      r = hypot( hypot(x,y), z+this%rearth)
      if (r > this%rearth*(1.d0+3.d-6)) then
         call add(errmsg,2,'point above earth surface',myname)
         return
      else
         r = min(this%rearth,r)
      end if
   !
      lev = locate(r-this%rearth,this%nz,this%zlevel)           ! zlevel(lev) < r-rearth <= zlevel(lev+1)
      if (lev == 0) then                                        ! r is below bottom level
         call add(errmsg,2,'point below deepest level',myname)
         return
      end if
   !
      nh = this%nx*this%ny
      ic1 = (lev-1)*nh+1                              ! range of sorted cell indices in level lev
      ic2 = lev*nh
   !
   !  cell center within given level closest to current regular grid point
   !  remapped to original cell index
   !
      idx => this%sortidx(ic1:ic2)
      ic = minloc(hypot(this%cell_center(1,idx)-x,this%cell_center(2,idx)-y),1)
      icell = this%sortidx(ic1-1+ic)
   end subroutine findCellSpecfem3dInversionGrid
!---------------------------------------------------------------
!  Is a given true negative depth below the inversion grid
!
   function isBelowSpecfem3dInversionGrid(this,z,errmsg) result(res)
      class (specfem3d_inversion_grid) :: this
      type (error_message) :: errmsg
      double precision :: z
      logical :: res
   !
      if (.not. allocated(this%zlevel)) then
         call add(errmsg,2,'cell locator was not initiated','isBelowSpecfem3dInversionGrid')
         res = .true.
         return
      end if
   !
      if (z < this%zlevel(1)) then
         res = .true.
      else
         res = .false.
      end if
   end function isBelowSpecfem3dInversionGrid
!
end module specfem3dInversionGrid
