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
!   Abstract base class for inversion grids
!     - define real grid type as extension to inversion_grid
!     - real grid type provides definitions of deferred procedures
!     - routines creating real grid type allocate inversion_grid to real grid type
!     - for example:
!         + declare <class (inversion_grid), allocatable :: invgrid> in the creating routine
!         + then for example: <allocate(real_inversion_grid :: invgrid> and
!         + call invgrid%createRealGrid(......)
!     - pointers and dummy arguments of invgrid used elsewhere should be declared
!       as <class (inversion_grid), pointer :: invgrid>
!
module inversionGrid
   use errorMessage
   use inputParameter
   implicit none
   type, abstract :: inversion_grid
      integer :: nwp_all                                                ! total number of wavefield points in inversion grid
      integer :: nwp_local                                              ! number of wavefield points owned by one process
      integer :: offwp                                                  ! offset of proc-specific points in full wp array
      integer :: ncell_all                                              ! total number of cells of inversion grid
      integer :: ncell_local                                            ! number of inversion cells owned by one process
      integer :: offcell                                                ! offset of proc-specific cells in full cell array
      integer :: ngll                                                   ! number of wavefield points belonging to cell
      double precision, dimension(:,:), allocatable :: cell_center      ! (3,ncell_local)-array: xyz-coords of cell centers
      double precision, dimension(:,:,:), allocatable :: cell_corner    ! (3,8,ncell_local)-array: xyz-coords of cell corners
      integer, dimension(:,:), allocatable :: face_neighbour            ! contains cell indices of its neighbours (7,ncell_local)
   contains
      procedure (createInversionGrid), deferred :: create
      procedure (deallocateInversionGrid), deferred :: dealloc
      procedure (getIntegrationWeightsInversionGrid), deferred :: getIntegrationWeights
      procedure (getModelValuesOnInversionGrid), deferred :: getModelValues
      procedure (getRealCellVolumesInversionGrid), deferred :: getRealCellVolumes
      procedure (getWavefieldPointsInversionGrid), deferred :: getWavefieldPoints
      procedure (getWpModelValuesInversionGrid), deferred :: getWpModelValues
      procedure (getGlobalDimensionsInversionGrid), deferred :: getGlobalDimensions
      procedure (getEarthRadiusInversionGrid), deferred :: getEarthRadius
      procedure (getGridCenterInversionGrid), deferred :: getGridCenter
      procedure :: getNwpAll => getNwpAllInversionGrid
      procedure :: getNwpLocal => getNwpLocalInversionGrid
      procedure :: getWpOffset => getWpOffsetInversionGrid
      procedure :: getNcellAll => getNcellAllInversionGrid
      procedure :: getNcellLocal => getNcellLocalInversionGrid
      procedure :: getCellOffset => getCellOffsetInversionGrid
      procedure :: getNgll => getNgllInversionGrid
      procedure :: getCellCenters => getCellCentersInversionGrid
      procedure :: getSelectedCellCenter => getSelectedCellCenterInversionGrid
      procedure :: getCellCentersOfFaceNeigbours => getCellCentersOfFaceNeigboursInversionGrid
      procedure :: getCellCorners => getCellCornersInversionGrid
      procedure :: getFaceNeighbours => getFaceNeighboursInversionGrid
      procedure :: cellHasLateralBoundaryFace => cellHasLateralBoundaryFaceInversionGrid
      procedure :: getSelectedFaceNeighbours => getSelectedFaceNeighboursInversionGrid
      procedure :: getRealCellRadii => getRealCellRadiiInversionGrid
      procedure :: getGeometryVtk => getGeometryVtkInversionGrid
      procedure :: nextCell => nextCellInversionGrid
   end type inversion_grid
!------------------------------------------------------------------------
!  Create inversion grid
!  This routine should always be called from iteration step basics
!  path:       path to file with grid specifications
!  do_hyperslab:   read wavefield point related arrays in hyperslab mode
!  inpar_inv:      input parameters of main inversion file
!  errmsg:         error message
!
   abstract interface
      subroutine createInversionGrid(this,iter_path,do_hyperslab,inpar_inv,errmsg)
         import inversion_grid
         import error_message
         import input_parameter
         class (inversion_grid) :: this
         character(len=*) :: iter_path
         logical :: do_hyperslab
         type (input_parameter) :: inpar_inv
         type (error_message) :: errmsg
      end subroutine createInversionGrid
   end interface
!------------------------------------------------------------------
!  deallocate inversion grid
!
   abstract interface
      subroutine deallocateInversionGrid(this)
         import inversion_grid
         class (inversion_grid) :: this
      end subroutine deallocateInversionGrid
   end interface
!------------------------------------------------------------------------
!  Return model values on inversion grid
!  If in multi_process hyperslab mode, return property on my share of invgrid
!  Values should then be inserted into global grid using this%offcell
!  In in no_hyperslab single process mode, values on full grid are returned
!  res:     pointer to property values
!
   abstract interface
      function getModelValuesOnInversionGrid(this) result(res)
         import inversion_grid
         class (inversion_grid), target :: this
         double precision, dimension(:,:), pointer :: res
      end function getModelValuesOnInversionGrid
   end interface
!------------------------------------------------------------------------
!  get pointer to integration weights
!
   abstract interface
      subroutine getIntegrationWeightsInversionGrid(this,weight)
         import inversion_grid
         class (inversion_grid), target :: this
         double precision, dimension(:,:), pointer :: weight
      end subroutine getIntegrationWeightsInversionGrid
   end interface
!------------------------------------------------------------------------
! Interface for obtaining complex system matrix
! x:      x-value where sysmat is computed
! sysmat: system matrix returned by procedure
! ios:    error indicator
!
   abstract interface
      subroutine getRealCellVolumesInversionGrid(this,res)
         import inversion_grid
         class (inversion_grid) :: this
         real, dimension(:), pointer :: res
      end subroutine getRealCellVolumesInversionGrid
   end interface
!------------------------------------------------------------------------
!  get pointers to wavefield points
!
   abstract interface
      subroutine getWavefieldPointsInversionGrid(this,x,y,z)
         import inversion_grid
         class (inversion_grid), target :: this
         double precision, dimension(:), pointer :: x,y,z
      end subroutine getWavefieldPointsInversionGrid
   end interface
!------------------------------------------------------------------------
!  get pointers to model values
!
   abstract interface
      subroutine getWpModelValuesInversionGrid(this,rho,vp,vs)
         import inversion_grid
         class (inversion_grid), target :: this
         double precision, dimension(:), pointer :: rho,vp,vs
      end subroutine getWpModelValuesInversionGrid
   end interface
!------------------------------------------------------------------------
!  get global dimensions of inversion grd
!
   abstract interface
      subroutine getGlobalDimensionsInversionGrid(this,n1,n2,n3)
         import inversion_grid
         class (inversion_grid) :: this
         integer :: n1,n2,n3
      end subroutine getGlobaldimensionsInversionGrid
   end interface
!------------------------------------------------------------------------
!  get earth radius
!
   abstract interface
      function getEarthRadiusInversionGrid(this) result(r)
         import inversion_grid
         class (inversion_grid) :: this
         double precision :: r
      end function getEarthRadiusInversionGrid
   end interface
!------------------------------------------------------------------------
!  get coordinates of center of inversion grid
!
   abstract interface
      subroutine getGridCenterInversionGrid(this,lat,lon)
         import inversion_grid
         class (inversion_grid) :: this
         double precision :: lat,lon
      end subroutine getGridCenterInversionGrid
   end interface
!
   contains
!------------------------------------------------------------------------
!  get total number of wavefield points in inversion grid
!
   integer function getNwpAllInversionGrid(this) result(res)
      class (inversion_grid), intent(in) :: this
      res = this%nwp_all
   end function getNwpAllInversionGrid
!------------------------------------------------------------------------
!  get number of proc-specific wavefield points
!
   integer function getNwpLocalInversionGrid(this) result(res)
      class (inversion_grid), intent(in) :: this
      res = this%nwp_local
   end function getNwpLocalInversionGrid
!------------------------------------------------------------------------
!  get offsett of proc-specific wavefield points
!
   integer function getWpOffsetInversionGrid(this) result(res)
      class (inversion_grid), intent(in) :: this
      res = this%offwp
   end function getWpOffsetInversionGrid
!------------------------------------------------------------------------
!  get ngll
!
   integer function getNgllInversionGrid(this) result(ngll)
      class (inversion_grid), intent(in) :: this
      ngll = this%ngll
   end function getNgllInversionGrid
!------------------------------------------------------------------------
!  get total number of cells of inversion grid
!
   integer function getNcellAllInversionGrid(this) result(ncell)
      class (inversion_grid), intent(in) :: this
      ncell = this%ncell_all
   end function getNcellAllInversionGrid
!------------------------------------------------------------------------
!  get number of cells owned by this process
!
   integer function getNcellLocalInversionGrid(this) result(res)
      class (inversion_grid), intent(in) :: this
      res = this%ncell_local
   end function getNcellLocalInversionGrid
!------------------------------------------------------------------------
!  get cell offset for this process
!
   integer function getCellOffsetInversionGrid(this) result(res)
      class (inversion_grid), intent(in) :: this
      res = this%offcell
   end function getCellOffsetInversionGrid
!------------------------------------------------------------------------
!  get pointer to cell centers
!
   subroutine getCellCentersInversionGrid(this,cc)
      class (inversion_grid), target :: this
      double precision, dimension(:,:), pointer :: cc
      cc => this%cell_center
   end subroutine getCellCentersInversionGrid
!----------------------------------------------------------------------
!  get coordinates of cell center of given cell
!
   subroutine getSelectedCellCenterInversionGrid(this,icell,xc,yc,zc)
      class (inversion_grid) :: this
      integer :: icell
      double precision :: xc,yc,zc
   !
      xc = this%cell_center(1,icell)
      yc = this%cell_center(2,icell)
      zc = this%cell_center(3,icell)
   end subroutine getSelectedCellCenterInversionGrid
!------------------------------------------------------------------------
!  get coordinates of cell centers of neighbouring cells
!  xc,yc,zc:   dummy array arguments at least dimensioned to 6 on calling side
!
   subroutine getCellCentersOfFaceNeigboursInversionGrid(this,icell,nb,xc,yc,zc)
      class (inversion_grid) :: this
      integer :: icell,nb
      double precision, dimension(:) :: xc,yc,zc
      integer :: i,idx
   !
      nb = this%face_neighbour(1,icell)
      do i = 1,nb
         idx = this%face_neighbour(i+1,icell)
         xc(i) = this%cell_center(1,idx)
         yc(i) = this%cell_center(2,idx)
         zc(i) = this%cell_center(3,idx)
      end do
   end subroutine getCellCentersOfFaceNeigboursInversionGrid
!------------------------------------------------------------------------
!  get pointer to cell corners
!
   subroutine getCellCornersInversionGrid(this,corner)
      class (inversion_grid), target :: this
      double precision, dimension(:,:,:), pointer :: corner
      corner => this%cell_corner
   end subroutine getCellCornersInversionGrid
!------------------------------------------------------------------------
!  get pointer to face neighbours
!
   subroutine getFaceNeighboursInversionGrid(this,nb)
      class (inversion_grid), target :: this
      integer, dimension(:,:), pointer :: nb
      nb => this%face_neighbour
   end subroutine getFaceNeighboursInversionGrid
!------------------------------------------------------------------------
!  get indices of face neighbours of a selected cell
!  nbidx:   dummy array arguments at least dimensioned to 6 on calling side
!
   subroutine getSelectedFaceNeighboursInversionGrid(this,icell,nb,nbidx)
      class (inversion_grid), target :: this
      integer :: icell,nb
      integer, dimension(:) :: nbidx
   !
      nb = this%face_neighbour(1,icell)
      nbidx(1:nb) = this%face_neighbour(2:nb+1,icell)
   end subroutine getSelectedFaceNeighboursInversionGrid
!----------------------------------------------------------------------------
!  check if cell has a lateral boundary face pointing either S, E, N or W
!  do this by finding out if the neighbour in any of these directions is missing
!  iface: 1=S, 2=E, 3=N, 4=W, 5=bot, 6=top
!
   function cellHasLateralBoundaryFaceInversionGrid(this,icell) result(res)
      class (inversion_grid), target :: this
      integer :: icell
      logical :: res
      integer, dimension(6) :: dir_of_face, face
      integer :: nb,ib,idir
      double precision, dimension(3) :: cc,d
      double precision, dimension(3,6) :: cb
   !
      call getCellCentersOfFaceNeigboursInversionGrid(this,icell,nb,cb(1,:),cb(2,:),cb(3,:))
      if (nb == 6) then; res = .false.; return; endif
   !
   !  cell is at a boundary
      res = .true.
      dir_of_face = [-2,1,2,-1,-3,3]                ! face directions with sign (S,E,N,W,B,T)
      face = 0
      call getSelectedCellCenterInversionGrid(this,icell,cc(1),cc(2),cc(3))
   !
   !  find out which faces have neighbours
      do ib = 1,nb
         d = cb(:,ib)-cc(:)                               ! coordinate differences
         idir = maxloc(abs(d),1)                          ! greatest difference and its sign indicate face normal
         if (d(idir) < 0.d0) idir = -idir
         face(ib) = findloc(dir_of_face,idir,1)           ! this face has a neighbour
      end do
   !
   !  check if any of the 1,2,3,4 faces do not have neighbours
      if (any(face .ne. 1) .or. any(face .ne. 2) .or. any(face .ne. 3) .or. any(face .ne.4)) res = .false.
   !
   end function cellHasLateralBoundaryFaceInversionGrid
!------------------------------------------------------------------------
!  get radii of inversion grid cells (single precision)
!  rad radii of cells
!
   subroutine getRealCellRadiiInversionGrid(this,rad)
      class (inversion_grid) :: this
      real, dimension(:), pointer :: rad
      double precision, dimension(3,8) :: p
      integer :: i,j
   !
   !  p(:,i) is the vector pointing from cell center to i'th cell corner
   !
      allocate(rad(this%ncell_local))
      do i = 1,this%ncell_local
         rad(i) =  0.0
         do j = 1,8
            p(:,j) = this%cell_corner(:,j,i) - this%cell_center(:,i)
            rad(i) = max(rad(i),sum(p(:,j)*p(:,j)))
         enddo
         rad(i) = sqrt(rad(i))
      end do
   end subroutine getRealCellRadiiInversionGrid
!------------------------------------------------------------------------
! return geometry information on cells for vtk output
!
   subroutine getGeometryVtkInversionGrid(this,geometry_type,scale,points,cell_connectivity,cell_type)
      class (inversion_grid) :: this
      integer :: geometry_type
      double precision :: scale
      double precision, dimension(:,:), pointer :: points
      integer, dimension(:), pointer :: cell_connectivity,cell_type
      integer :: i,icell
   !
      nullify(points,cell_connectivity,cell_type)
      select case(geometry_type)
      case(0) ! CELLS
         allocate(points(3,8*this%ncell_local))
         allocate(cell_type(this%ncell_local))
         allocate(cell_connectivity((8+1)*this%ncell_local))
         cell_type = 12
         do icell = 1,this%ncell_local
            points(1:3,(icell-1)*8+1:icell*8) = this%cell_corner(1:3,1:8,icell)
            cell_connectivity((icell-1)*9+1) = 8
            cell_connectivity((icell-1)*9+2:icell*9) = (/ (i,i=(icell-1)*8,icell*8-1) /) !(note, that for vtk format, point indices have offset 0!)
         end do
      case(1) ! CELL CENTERS
         allocate(points(3,this%ncell_local))
         points = this%cell_center
      end select
   !
      points = points*scale
   end subroutine getGeometryVtkInversionGrid
!------------------------------------------------------------------------
!  Iterate over cells of inversion grid
!  If in multi_process hyperslab mode, loop through my cells
!  optionally returns cell center and cell index
!
   logical function nextCellInversionGrid(this,cc,icell)
      class (inversion_grid) :: this
      double precision, dimension(:), optional :: cc
      integer, optional :: icell
      integer :: call_count = 0
      save :: call_count
   !
      call_count = call_count+1
      if (call_count <= this%ncell_local) then
         if (present(cc)) cc = this%cell_center(:,call_count)
         if (present(icell)) icell = call_count
         nextCellInversionGrid = .true.
      else
         call_count = 0
         if (present(cc)) cc = 0.d0
         if (present(icell)) icell = 0
         nextCellInversionGrid = .false.
      endif
   end function nextCellInversionGrid

end module inversionGrid
