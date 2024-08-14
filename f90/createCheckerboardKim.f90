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
!
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!   Create a 3D checkerboard kernel inverted model (defined on inversion grid)
!
program createCheckerboardKim
   use mathConstants
   use inversionBasics
   use specfem3dInversionGrid
   use kernelInvertedModel
   use askiBackgroundModel
   use regularSphericalGrid
   use inputParameter
   use argumentParser
   use hdfWrapper
   use errorMessage
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (input_parameter), pointer :: inpar_inv
   class (inversion_grid), allocatable :: invgrid
   type (kernel_inverted_model) :: kim,abkim
   type (aski_background_model) :: abm
   type (regular_spherical_grid) :: rsg
   real, dimension(:), allocatable :: d
   real, dimension(:,:), pointer :: mval
   double precision :: dxpm,dypm,dzpm,dlon,dlat,dr,rmin,latmin,lonmin,wxpm,wypm,wzpm
   double precision :: r,lat,lon,rp,rearth,xc,yc,zc,plon,plat,lat_box_center,lon_box_center,thetac,reldev,m1d
   integer, dimension(3) :: griddim,blockdims
   integer :: it,iprop,iprop_abm,nr,nlat,nlon,icell,ncell,k,j,i,br,blat,blon,ib,jb,kb,kbshift
   logical :: ko,jo,io,ke,je,ie,k2o,j2o,i2o,k2e,j2e,i2e
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: main_path,iter_path,outpath
   character(len=max_length_string) :: abmfile,specfem_input
   character(len=char_len_par), dimension(1) :: prop
   character(len=4) :: cbldims
   character(len=21) :: myname = 'createCheckerboardKim'
!
   call new(errmsg,myname)
!
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 1
! ------------------------------------------------------------------------
!  command line processing
!
   call init(ap,myname,"Create a checkerboard model on the inversion grid")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-prop',.true.,'model property','sval','vp')
   call addOption(ap,'-reldev',.true.,'relative deviation','dval','0.05')
   call addOption(ap,'-it',.true.,'take inversion grid from specified iteration (def = current)','ival','0')
   call addOption(ap,'-blockdims',.true.,'block dimensions (br,blat,blon)','ivec','2 2 2')
   call addOption(ap,'-vblshift',.true.,'shift blocks vertically by given number of cell layers (def = 0)','ival','0')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
   call document(ap)
!
   main_parfile = ap.sval.'main_parfile'
   prop(1) = ap.sval.'-prop'
   reldev = ap.dval.'-reldev'
   it = ap.ival.'-it'
   blockdims = ap.ivec.'-blockdims'
   br = blockdims(1); blat = blockdims(2); blon = blockdims(3)
   write(cbldims,'(a1,3i1)') '_',br,blat,blon
   kbshift = ap.ival.'-vblshift'
   print *,'Pattern shifted down by ',kbshift,' layers' 
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
   wxpm = inpar_inv.dval.'ASKI_wx'
   wypm = inpar_inv.dval.'ASKI_wy'
   wzpm = inpar_inv.dval.'ASKI_wz'
   d = rvec(inpar_inv,'PSEUDO_MESH_SPACING',3)
   dxpm = d(1); dypm = d(2); dzpm = d(3); deallocate(d)
   griddim = ivec(inpar_inv,'INVGRID_DIMENSIONS',3)
   specfem_input = inpar_inv.sval.'PATH_SPECFEM_INPUT'
   abmfile = trim(specfem_input)//trim(inpar_inv.sval.'FILE_ASKI_BACKGROUND_MODEL')
!
!  take model perturbations and inversion grid either from current iteration or from specified one
!
   if (it > 0) then
      write(iter_path,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it,'/'
   else
      iter_path = .iterpath.invbasics
      it = inpar_inv.ival.'CURRENT_ITERATION_STEP'
   end if
!
!  read inversion grid
!
   invgrid_type = trim(inpar_inv.sval.'INVERSION_GRID')
   if (equalString(invgrid_type,'specfem3d')) then
      allocate(specfem3d_inversion_grid :: invgrid)
      call invgrid%create(iter_path,.false.,inpar_inv,errmsg)
      if (.level.errmsg == 2) return
   else
      call add(errmsg,2,'Inversion grid type not implemented',myname)
      return
   end if
!
!  rearth and
!  coordinates of pole the grid coordinates refer to in range (-pi,pi)
!
   rearth = invgrid%getEarthRadius()
   call invgrid%getGridCenter(lat_box_center,lon_box_center)
   plon = mc_pid+lon_box_center
   if (plon > mc_pid) plon = plon-mc_two_pid
   plat = 0.5*mc_pid-lat_box_center
   thetac = plat

   !  path for kim output
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
!
!  get ASKI background model
!
   call readASKIBackgroundModel(abm,abmfile,1,rearth,errmsg)
   if (.level.errmsg == 2) goto 1
   call createFromAbmKernelInvertedModel(abkim,abm,invgrid,thetac,lon_box_center,errmsg)
   if (.level.errmsg == 2) goto 1
!
   iprop = 1
   iprop_abm = findloc(abkim%prop,prop(1),1)
   print *,'Model property ',trim(prop(1)),' has index ',iprop_abm,' in ASKI background model'
!
!  generate equator centered SPECFEM pseudo mesh as regular spherical grid
!
   dlon = dxpm/rearth
   dlat = dypm/rearth
   dr = dzpm
   nlon = griddim(1)
   nlat = griddim(2)
   nr = griddim(3)
   write(6,'(a,3i8)') 'Dimensions of pseudo mesh: ',nlon,nlat,nr
   lonmin = (-0.5*wxpm+0.5*dxpm)/rearth
   latmin = (-0.5*wypm+0.5*dypm)/rearth
   rmin = rearth-(wzpm-0.5*dzpm)
   call rsg%create(1,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'reldev','cb'//cbldims)
!
!  step through elements ony by one
!  here we transform from the SPECFEM Cartesian element coordinates
!  to their spherical coordinates in the equator centered pseudo mesh
!
   ncell = invgrid%getNcellAll()
   allocate(mval(ncell,1))
!
!  set all zeros first and then fill the positives and negatives
!
   mval = 0.0
   do while (invgrid%nextCell(icell = icell))
      m1d = abkim%model_values(icell,iprop_abm)
      call invgrid%getSelectedCellCenter(icell,xc,yc,zc)
      r = hypot(hypot(zc+rearth,xc),yc)
      lat = dasin(yc/r)
      rp = r*dcos(lat)
      lon = dasin(xc/rp)
      call rsg%getIndicesClosest(r,lat,lon,k,j,i)
   !
   !  block indices of cell and if it is odd or even
   !
      kb = (k-1+kbshift)/br+1
      jb = (j-1)/blat+1
      ib = (i-1)/blon+1
      ke = (mod(kb,2) == 0); ko = .not. ke
      je = (mod(jb,2) == 0); jo = .not. je
      ie = (mod(ib,2) == 0); io = .not. ie
      if (ko) then; k2e = (mod((kb+1)/2,2) == 0); k2o = .not. k2e; endif
      if (ke) then; k2e = (mod(kb/2,2) == 0); k2o = .not. k2e; endif
      if (jo) then; j2e = (mod((jb+1)/2,2) == 0); j2o = .not. j2e; endif
      if (je) then; j2e = (mod(jb/2,2) == 0); j2o = .not. j2e; endif
      if (io) then; i2e = (mod((ib+1)/2,2) == 0); i2o = .not. i2e; endif
      if (ie) then; i2e = (mod(ib/2,2) == 0); i2o = .not. i2e; endif
   !
   !  set either positive or negative perturbations
   !
      if ((ko .and. jo .and. io)) then
         if ((j2o .and. i2o) .or. (j2e .and. i2e)) then
            if (k2o) then; mval(icell,iprop) = +reldev*m1d; else; mval(icell,iprop) = -reldev*m1d; endif
         else if ((j2o .and. i2e) .or. (j2e .and. i2o)) then
            if (k2o) then; mval(icell,iprop) = -reldev*m1d; else; mval(icell,iprop) = +reldev*m1d; endif
         end if
      endif
   enddo
!
   call createFromValuesKernelInvertedModel(kim,prop,mval,'isoVelocitySI','absdev',errmsg)
   if (.level.errmsg == 2) goto 1
   call writeHDFKernelInvertedModel(kim,trim(outpath)//'checkerboard_'//trim(prop(1))//cbldims//'_absdev.kim',invgrid,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   call dealloc(kim)
   call dealloc(abm)
   call dealloc(abkim)
   call invgrid%dealloc()
   deallocate(invgrid)
   call rsg%dealloc()
   call dealloc(invbasics)
   call dealloc(errmsg)
   call dealloc(ap)
!------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif
end program createCheckerboardKim
