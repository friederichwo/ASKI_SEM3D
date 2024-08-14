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
!
!  New main program for writing VTK files on 1 CPU only
!
program writeVtkFiles
   use inversionBasics
   use iterationStepBasics
   use vtkFile
   use specfem3dInversionGrid
   use string
   use argumentParser
   use fileUnitHandler
   use errorMessage
   use mpiSupport
   use globalMpiInfo
   implicit none
!
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics
   type (vtk_info) :: wp_vtk,invgrid_vtk
   type (mpi_support) :: mpisup
   class (inversion_grid), allocatable :: invgrid
   type (input_parameter), pointer :: inpar_inv,inpar_iter
   real, dimension(:), pointer :: rdata
   double precision, dimension(:,:), pointer :: weight
   double precision, dimension(:), allocatable :: sumweight
   double precision, dimension(:), pointer :: rho,vp,vs
   double precision :: vtk_scale
   integer :: j,ia,ie,ngll,ncell
   logical :: no_wp
   character(len=10) :: myname = 'writeVtkFiles'
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: vtk_format,vtkpath,vtk_geometry_type
!-----------------------------------------------------------------------------
!  initialise MPI
!
   call new(mpisup)
   myrank = .myrank.mpisup
   numtasks = .numtasks.mpisup
!------------------------------------------------------------------------
!  preliminary processing
!
   call init(ap,myname,"Write diverse inversion grid informaton to VTK files")
   call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
   call addOption(ap,'-nowp',.false.,'If set, skip VTk files on wavefield points')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      call abort(mpisup)
   end if
   call document(ap)
   main_parfile = ap.sval.'main_parfile'
   no_wp = ap.optset.'-nowp'
   call dealloc(ap)
!------------------------------------------------------------------------
!  setup basics
!
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   call initiateIterationStepBasics(iterbasics,invbasics,1,errmsg)
   if (.level.errmsg == 2) goto 1
!
   inpar_inv => getInputParameterInversionBasics(invbasics)
   inpar_iter =>  getInputParameterIterationStepBasics(iterbasics)
   vtk_format = trim(inpar_inv.sval.'DEFAULT_VTK_FILE_FORMAT')
   vtkpath = trim(.iterpath.iterbasics)//trim(inpar_inv.sval.'PATH_VTK_FILES')//'stats'
   vtk_scale = inpar_iter.dval.'VTK_COORDS_SCALING_FACTOR'
   vtk_geometry_type = inpar_iter.sval.'VTK_GEOMETRY_TYPE'
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(.iterpath.invbasics,.false.,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
   ncell = invgrid%getNcellLocal()
   call invgrid%getWpModelValues(rho,vp,vs)
!
!  write inversion grid to vtk file
!
   print *,'write inversion grid to vtk file'
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeHeaderVtkFile(invgrid_vtk,vtkpath,1,errmsg)
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
!
!  write sum of integration weights per invgrid cell (should correspond to volume) to invgridVtkFile
!
   print *,'write sum of integration weights (volume) to vtk file'
   call invgrid%getIntegrationWeights(weight)
   allocate(rdata(ncell),sumweight(ncell))
   do j = 1,ncell
      sumweight(j) = sum(weight(:,j))
      rdata(j) = sumweight(j)
   enddo
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(invgrid_vtk,trim(vtkpath)//"_sum_intw",1,rdata,errmsg,data_name='sum_of_integration_weights')
   if (.level.errmsg == 2) return
   deallocate(rdata)
   call dealloc(invgrid_vtk)
!
! write cell volume per invgrid cell to invgridVtkFile
!
   print *,'write cell volume to vtk file'
   call invgrid%getRealCellVolumes(rdata)
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(invgrid_vtk,trim(vtkpath)//"_cell_vol",1,rdata,errmsg,data_name='cell_volume')
   if (.level.errmsg == 2) return
   deallocate(rdata)
   call dealloc(invgrid_vtk)
!
!  rho on grid
!
   print *,'write rho on grid points to vtk file'
   ngll = invgrid%getNgll()
   allocate(rdata(ncell))
   ia = 1; ie = ngll
   do j = 1,ncell
      rdata(j) = sum(weight(:,j)*rho(ia:ie))/sumweight(j)
      ia = ia+ngll; ie = ie+ngll
   enddo
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(invgrid_vtk,trim(vtkpath)//'_rho_on_grid',1,rdata,errmsg,data_name='rho_on_grid')
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
!
!  vp on grid
!
   print *,'write vp on grid points to vtk file'
   ia = 1; ie = ngll
   do j = 1,ncell
      rdata(j) = sum(weight(:,j)*vp(ia:ie))/sumweight(j)
      ia = ia+ngll; ie = ie+ngll
   enddo
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(invgrid_vtk,trim(vtkpath)//'_vp_on_grid',1,rdata,errmsg,data_name='vp_on_grid')
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
!
!  vs on grid
!
   print *,'write vs on grid points to vtk file'
   ia = 1; ie = ngll
   do j = 1,ncell
      rdata(j) = sum(weight(:,j)*vs(ia:ie))/sumweight(j)
      ia = ia+ngll; ie = ie+ngll
   enddo
   call initiateVtkFile(invgrid_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,1,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(invgrid_vtk,trim(vtkpath)//'_vs_on_grid',1,rdata,errmsg,data_name='vs_on_grid')
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
   deallocate(rdata)
!
   if (no_wp) goto 2
!
!  write wavefield points to vtk file
!
   print *,'write wavefield points to vtk file'
   call initiateVtkFile(wp_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,2,errmsg)
   if (.level.errmsg == 2) return
   call writeHeaderVtkFile(wp_vtk,vtkpath,1,errmsg)
   if (.level.errmsg == 2) return
   call dealloc(wp_vtk)
!
! write rho, vp, vs on wavefield points
!
   print *,'write density on wp points to vtk file'
   allocate(rdata(size(rho)))
   rdata = rho
   call initiateVtkFile(wp_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,2,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(wp_vtk,trim(vtkpath)//'_rho_on_wp',1,rdata,errmsg,data_name='rho_on_wavefield_points')
   if (.level.errmsg == 2) return
   call dealloc(wp_vtk)
   !
   print *,'write vp on wp points to vtk file'
   rdata = vp
   call initiateVtkFile(wp_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,2,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(wp_vtk,trim(vtkpath)//'_vp_on_wp',1,rdata,errmsg,data_name='vp_on_wavefield_points')
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
   !
   print *,'write vs on wp points to vtk file'
   rdata = vs
   call initiateVtkFile(wp_vtk,invgrid,vtk_format,vtk_geometry_type,vtk_scale,2,errmsg)
   if (.level.errmsg == 2) return
   call writeRealDataVtkFile(wp_vtk,trim(vtkpath)//'_vs_on_wp',1,rdata,errmsg,data_name='vs_on_wavefield_points')
   if (.level.errmsg == 2) return
   call dealloc(invgrid_vtk)
   deallocate(rdata)

 1 if (.level.errmsg == 2) then
      call print(errmsg)
      call abort(mpisup)
   end if
 2 call dealloc(mpisup)
   call dealloc(iterbasics)
   call dealloc(invbasics)
   call invgrid%dealloc(); deallocate(invgrid)
end program writeVtkFiles
