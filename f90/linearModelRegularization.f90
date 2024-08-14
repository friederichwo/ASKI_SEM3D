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
!   Define and handle linear equations, serving as smoothing or damping constraints,
!   to be added to the kernel equations.
!
!   As the regularization equations will generally be very sparse, avoid storing zeros
!   by only storing relevant coefficients of an equation
!   (along with the indices of the associated variables of that equation)
!
module linearModelRegularization
!
   use mpi
   use inversionGrid
   use errorMessage
   use globalMpiInfo
!
   implicit none
!
   interface dealloc; module procedure deallocateLinearModelRegularization; end interface
   interface operator (.neq.); module procedure getNeqtotalLinearModelRegularization ; end interface
!
!  Type contains indices and values which define regularization equations
!  Only deal with equations of the current process
!  Note: this type does not have a solution vector or related vectors like gradient and p.
!        They are kept by kernelLinearSystem
!
   type linear_model_regularization
      private
      integer :: nprop                                                   ! number of different property names in model space
      integer :: offcell                                                 ! offset of cell range handled by this process
      integer :: ncell                                                   ! number of inversion grid cells of this process
      integer :: neq                                                     ! number of equations handled by this process
   !  smoothing
      double precision, dimension(:), allocatable :: vscal_smooth_mantle ! scaling values for mantle smoothing (one per property)
      double precision, dimension(:), allocatable :: vscal_smooth_crust  ! scaling values for crust smoothing (one per property)
      double precision :: dsmboost                                       ! boosting factor for depth smoothing
      integer, dimension(:,:), allocatable :: smidx                      ! indices of model values with non-zero entry (cell+6 neighbours)
      double precision, dimension(:,:), allocatable :: smval             ! values of non-zero entries
      double precision, dimension(:), allocatable :: smrhs               ! rhs values of equations
      double precision, dimension(:), allocatable :: smq                 ! q values of equations
   !  damping
      double precision, dimension(:), allocatable :: vscal_damp_mantle   ! scaling values for mantle damping (one per property)
      double precision, dimension(:), allocatable :: vscal_damp_crust    ! scaling values for crust damping (one per property)
      double precision :: boundchoke                                     ! damping weight for model boundaries
      integer, dimension(:), allocatable :: dpidx                        ! index of model value with non-zero entry
      double precision, dimension(:), allocatable :: dpval               ! value of non-zero entry
      double precision, dimension(:), allocatable :: dprhs               ! rhs values of equations
      double precision, dimension(:), allocatable :: dpq                 ! q values of equations
   end type linear_model_regularization
!
   contains
!
!------------------------------------------------------------------------------------
!  Initialize regularization object
!  nprop: number of model properties to be inverted for
!  vscal_smooth: smoothing scaling values for these properties
!  vscal_damp:   dampiing scaling values for these properties
!  first_cell:   start of grid cell range handled by this process
!  ncell:        number of grid cells handled by this process
!
   subroutine initiateLinearModelRegularization(this,nprop,vscal_smooth_mantle,vscal_damp_mantle,&
                                                vscal_smooth_crust,vscal_damp_crust,dsmboost,&
                                                boundchoke,offcell,ncell,errmsg)
      type (linear_model_regularization) :: this
      double precision, dimension(:) :: vscal_smooth_mantle,vscal_smooth_crust
      double precision, dimension(:) :: vscal_damp_mantle, vscal_damp_crust
      double precision :: dsmboost,boundchoke
      integer :: nprop,offcell,ncell
      type (error_message) :: errmsg
   !
      this%nprop = nprop
      this%offcell = offcell
      this%ncell = ncell
      this%neq = this%nprop*this%ncell               ! number of my equations for smoothing or damping, respectively
   !
   !  smoothing
   !
      allocate(this%vscal_smooth_mantle(nprop),this%vscal_smooth_crust(nprop))
      this%vscal_smooth_mantle = vscal_smooth_mantle
      this%vscal_smooth_crust = vscal_smooth_crust
      this%dsmboost = dsmboost
      allocate(this%smidx(7,this%neq))               ! cell plus max 6 neighbours
      allocate(this%smval(7,this%neq))
      allocate(this%smrhs(this%neq))
      allocate(this%smq(this%neq))
   !
   !  damping
   !
      allocate(this%vscal_damp_mantle(nprop),this%vscal_damp_crust(nprop))
      this%vscal_damp_mantle = vscal_damp_mantle
      this%vscal_damp_crust = vscal_damp_crust
      this%boundchoke = boundchoke
      allocate(this%dpidx(this%neq))
      allocate(this%dpval(this%neq))
      allocate(this%dprhs(this%neq))
      allocate(this%dpq(this%neq))
   end subroutine initiateLinearModelRegularization
!------------------------------------------------------------------------------------
!  Compute Laplace smoothing (average over neighbours) regularization equations
!
   subroutine smoothingLinearModelRegularization(this,invgrid,crdepth)
      type (linear_model_regularization) :: this
      class (inversion_grid) :: invgrid
      integer :: iprop,icell,nnb,ib,j,ieq,ncell_all
      integer, dimension(:,:), pointer :: face_nb
      double precision :: f,div,xc,yc,zc,depth,rearth,wsmooth,crdepth
   !
      ncell_all = invgrid%getNcellAll()
      call invgrid%getFaceNeighbours(face_nb)
      ieq = 0                                                       ! my local equation counter
      do iprop = 1,this%nprop
         do icell = this%offcell+1,this%offcell+this%ncell          ! my share of global cell range
            ieq = ieq+1
            this%smidx(:,ieq) = icell+(iprop-1)*ncell_all           ! indices of missing neighbours set to center cell
            this%smval(:,ieq) = 0.d0                                ! values of missing neighbours stay zero
            nnb = face_nb(1,icell)                                  ! number of face neighbours

            call invgrid%getSelectedCellCenter(icell,xc,yc,zc)      ! get depth of cell for depth dependent smoothing
            rearth = invgrid%getEarthRadius()
            depth = rearth-hypot(hypot(xc,yc),zc+rearth)

            f = 1.d0
            if (nnb == 6) f = this%dsmboost                         ! enhanced depth smoothing for inner cells
            div = 4.+2.*f
            if (depth > crdepth) then
               wsmooth = this%vscal_smooth_mantle(iprop)
            else
               wsmooth = this%vscal_smooth_crust(iprop)
            endif

            do j = 1,nnb
               ib = face_nb(j+1,icell)                                      ! global cell index of neighbour j
               this%smidx(j,ieq) = ib+(iprop-1)*ncell_all                   ! model value index with non-zero entry
               if (j > 4) then                                              ! top and bottom face
                  this%smval(j,ieq) = f/div*wsmooth                         ! value of non-zero entry
               else                                                         ! other faces
                  this%smval(j,ieq) = 1.d0/div*wsmooth
               end if
            end do
            this%smidx(nnb+1,ieq) = icell+(iprop-1)*ncell_all                ! center cell
            this%smval(nnb+1,ieq) = -1.d0*wsmooth
            this%smrhs(ieq) = 0.d0
         end do
      enddo
   end subroutine smoothingLinearModelRegularization
!------------------------------------------------------------------------------------
!  Compute damping regularization equations
!
   subroutine dampingLinearModelRegularization(this,invgrid,crdepth,wx,wy)
      type (linear_model_regularization) :: this
      class (inversion_grid) :: invgrid
      double precision :: crdepth,wx,wy
      integer :: ncell_all
      integer :: iprop,icell,ieq
      integer, dimension(:,:), pointer :: face_nb
      double precision :: xc,yc,zc,rearth,depth,wdamp,beta,phi,xpm,ypm,r
   !
      call invgrid%getFaceNeighbours(face_nb)
      ieq = 0
      do iprop = 1,this%nprop
         do icell = this%offcell+1,this%offcell+this%ncell              ! my share of global cell range
            ieq = ieq+1
            this%dpidx(ieq) = icell+(iprop-1)*ncell_all

            call invgrid%getSelectedCellCenter(icell,xc,yc,zc)          ! get depth of cell for depth dependent smoothing
            rearth = invgrid%getEarthRadius()
            r = hypot(hypot(xc,yc),zc+rearth)
            depth = rearth-r
            beta = asin(yc/r)                                           ! get pseudo mesh coordinates of point
            phi = asin(xc/(r*cos(beta)))
            xpm = rearth*phi
            ypm = rearth*beta
            if (depth > crdepth) then
               wdamp = this%vscal_damp_mantle(iprop)
            else
               wdamp = this%vscal_damp_crust(iprop)
            endif

            if (invgrid%cellHasLateralBoundaryFace(icell)) then
!            if (abs(xpm) > 0.25*wx .or. abs(ypm) > 0.25*wy) then
                this%dpval(ieq) = +this%boundchoke*wdamp
            else
                this%dpval(ieq) = 1.d0*wdamp
            end if
            this%dprhs(ieq) = 0.d0
         end do
      enddo
   end subroutine dampingLinearModelRegularization
!------------------------------------------------------------------------
!  compute in place rhs - A*sol (using my values of rhs only)
!
   subroutine rhsMinusMatrixDotSolLinearModelRegularization(this,sol)
      type (linear_model_regularization) :: this
      double precision, dimension(:) :: sol
      double precision :: sum
      integer :: i,j
   !
      do i = 1,this%neq
         sum = 0.d0
         do j = 1,7
            sum = sum+this%smval(j,i)*sol(this%smidx(j,i))
         end do
         this%smrhs(i) = this%smrhs(i)-sum
         this%dprhs(i) = this%dprhs(i)-this%dpval(i)*sol(this%dpidx(i))
      end do
   end subroutine rhsMinusMatrixDotSolLinearModelRegularization
!------------------------------------------------------------------------
!  compute norm of my share of rhs and then do allreduce
!
   subroutine computeNormSquaredRhsLinearModelRegularization(this,norm_sq,errmsg)
      type (linear_model_regularization) :: this
      double precision :: norm_sq_sm,norm_sq_dp,norm_sq_loc,norm_sq
      type (error_message) :: errmsg
      integer :: ios
      double precision :: DDOT
      external DDOT
   !
      norm_sq_loc = DDOT(this%neq,this%smrhs,1,this%smrhs,1)
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq_sm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of smrhs','computeNormSquaredRhsLinearModelRegularization')
         return
      end if
      norm_sq_loc = DDOT(this%neq,this%dprhs,1,this%dprhs,1)
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq_dp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of dprhs','computeNormSquaredRhsLinearModelRegularization')
         return
      end if
      norm_sq = norm_sq_sm+norm_sq_dp
   end subroutine computeNormSquaredRhsLinearModelRegularization
!------------------------------------------------------------------------
!  add regularization part of g = A^T rhs to that of kernel linear system
!
   subroutine addGradIsMatrixTransRhsLinearModelRegularization(this,g,errmsg)
      type (linear_model_regularization) :: this
      double precision, dimension(:) :: g
      double precision, dimension(:), allocatable :: g_loc,g_glob
      type (error_message) :: errmsg
      integer :: ios,nmval,i,j,k
   !
      nmval = size(g)
      allocate(g_loc(nmval),g_glob(nmval))
      g_loc = 0.d0
   !
   !  Go through my columns of matrix A^T and multiply each of them by the associated element of rhs
   !  and sum up over columns
   !  Do this for smoothing with 7 entries per column and damping with one entry per column
   !  In case of missing neighbours in smoothing, a column entry may occur more than once
   !  but since smval=0 in this case it does not change the sum
   !
      do i = 1,this%neq
         do j = 1,7                                             ! treat smoothing
            k = this%smidx(j,i)
            g_loc(k) = g_loc(k)+this%smval(j,i)*this%smrhs(i)
         end do
         k = this%dpidx(i)                                      ! treat damping
         g_loc(k) = g_loc(k)+this%dpval(i)*this%dprhs(i)
      end do
   !
      call MPI_ALLREDUCE(g_loc,g_glob,nmval,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather results of A^T*rhs','addGradIsMatrixTransRhsLinearModelRegularization')
         return
      end if
      deallocate(g_loc)
   !
   !  add g_glob to gradient from kernel linear system
   !
      g = g+g_glob
      deallocate(g_glob)
   end subroutine addGradIsMatrixTransRhsLinearModelRegularization
!------------------------------------------------------------------------
!  compute my share of q = A*p
!
   subroutine qIsMatrixDotPLinearModelRegularization(this,p)
      type (linear_model_regularization) :: this
      double precision, dimension(:) :: p
      double precision :: sum
      integer :: i,j
   !
      do i = 1,this%neq
         sum = 0.d0
         do j = 1,7
            sum = sum+this%smval(j,i)*p(this%smidx(j,i))
         end do
         this%smq(i) = sum
         this%dpq(i) = this%dpval(i)*p(this%dpidx(i))
      end do
   end subroutine qIsMatrixDotPLinearModelRegularization
!------------------------------------------------------------------------
!  compute norm of my share of q and then do allreduce
!
   subroutine computeNormSquaredQLinearModelRegularization(this,norm_sq,errmsg)
      type (linear_model_regularization) :: this
      double precision :: norm_sq_sm,norm_sq_dp,norm_sq_loc,norm_sq
      type (error_message) :: errmsg
      integer :: ios
      double precision :: DDOT
      external DDOT
   !
      norm_sq_loc = DDOT(this%neq,this%smq,1,this%smq,1)
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq_sm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of smq','computeNormSquaredQLinearModelRegularization')
         return
      end if
      norm_sq_loc = DDOT(this%neq,this%dpq,1,this%dpq,1)
      call MPI_ALLREDUCE(norm_sq_loc,norm_sq_dp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ios)
      if (ios .ne. MPI_SUCCESS) then
         call add(errmsg,2,'could not gather norm of dpq','computeNormSquaredQLinearModelRegularization')
         return
      end if
      norm_sq = norm_sq_sm+norm_sq_dp
   end subroutine computeNormSquaredQLinearModelRegularization
!------------------------------------------------------------------------
!  update rhs = rhs - alfa * q
!
   subroutine updateRhsLinearModelRegularization(this,alfa)
      type (linear_model_regularization) :: this
      double precision :: alfa
   !
      call DAXPY(this%neq,-alfa,this%smq,1,this%smrhs,1)
      call DAXPY(this%neq,-alfa,this%dpq,1,this%dprhs,1)
   end subroutine updateRhsLinearModelRegularization
!------------------------------------------------------------------------
!  Return total number of regularization equations
!
  function getNeqTotalLinearModelRegularization(this) result(neq)
    type (linear_model_regularization), intent(in) :: this
    integer :: neq
    neq = this%neq
  end function getNeqTotalLinearModelRegularization
!------------------------------------------------------------------------
!  Deallocate object
!
   subroutine deallocateLinearModelRegularization(this)
      type (linear_model_regularization) :: this
      if(allocated(this%smidx)) deallocate(this%smidx)
      if(allocated(this%smval)) deallocate(this%smval)
      if(allocated(this%smrhs)) deallocate(this%smrhs)
      if(allocated(this%smrhs)) deallocate(this%smq)
      if(allocated(this%dpidx)) deallocate(this%dpidx)
      if(allocated(this%dpval)) deallocate(this%dpval)
      if(allocated(this%dprhs)) deallocate(this%dprhs)
      if(allocated(this%dprhs)) deallocate(this%dpq)
      if (allocated(this%vscal_smooth_mantle)) deallocate(this%vscal_smooth_mantle)
      if (allocated(this%vscal_damp_mantle)) deallocate(this%vscal_damp_mantle)
      if (allocated(this%vscal_smooth_crust)) deallocate(this%vscal_smooth_crust)
      if (allocated(this%vscal_damp_crust)) deallocate(this%vscal_damp_crust)
   end subroutine deallocateLinearModelRegularization
end module linearModelRegularization

