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
!   Evaluate kernel inverted model defined on partitioned SPECFEM mesh
!   at element centers of pseudo mesh defined at the equator
!   which forms a regular spherical grid
!   with spacings given by the mesh spacings and number of grid points given
!   by the mesh dimensions.
!
program throwKimOnSphericalPseudoMesh
   use mathConstants
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use kernelInvertedModel
   use askiBackgroundModel
   use regularSphericalGrid
   use inputParameter
   use argumentParser
   use hdfWrapper
   use errorMessage
   use smartUtils
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (input_parameter), pointer :: inpar_inv
   class (inversion_grid), allocatable :: invgrid
   type (kernel_inverted_model) :: kim,krm,abkim,kimcs
   type (aski_background_model) :: abm
   type (regular_spherical_grid) :: rsg
   real, dimension(:), allocatable :: d
   double precision, dimension(:), allocatable :: mravg,reldevmin,reldevmax
   double precision :: dxpm,dypm,dzpm,dlon,dlat,dr,rmin,latmin,lonmin,wxpm,wypm,wzpm
   double precision :: r,lat,lon,rp,rearth,xc,yc,zc,plon,plat,lat_box_center,lon_box_center,thetac
   double precision :: mval,refval,m1d,colsum
   integer, dimension(3) :: griddim
   integer :: it,itref,iprop,iprop_ref,iprop_abm,nr,nlat,nlon,icell,ierr,k,j,i,ic
   logical :: excs
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: main_path,iter_path,outpath,iter_path_ref,outpath_ref,rsgpath
   character(len=max_length_string) :: modelfile,refmodelfile,colsumfile,vgridfile,abmfile,specfem_input
   character(len=char_len_par) :: prop
   character(len=15) :: myname = 'throwKimOnSphericalPseudoMesh'
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
   call init(ap,myname,"Throw kernel inverted model on regular spherical SPECFEM pseudo mesh")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-prop',.true.,'model property to be plotted','sval','vp')
   call addOption(ap,'-it',.true.,'take model perturbationsfrom specified iteration (def = current)','ival','0')
   call addOption(ap,'-itref',.true.,'subtract model perturbations from specified iteration','ival','0')
   call addOption(ap,'-modelfile',.true.,'name of modelfile under OUTPUT_FILES','sval','total_absdev_model.kim')
   call addOption(ap,'-rsgpath',.true.,'subpath of OUTPUT_FILES where rsg output goes','sval','kimrsg/')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   prop = ap.sval.'-prop'
   it = ap.ival.'-it'
   itref = ap.ival.'-itref'
   modelfile = ap.sval.'-modelfile'
   rsgpath = ap.sval.'-rsgpath'
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   iter_path = .iterpath.invbasics
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
!
!  take reference model perturbations either from specified iteration or none
!  default is none
!
   if (itref == 0) then
      iter_path_ref = 'not_used/'
   else
      write(iter_path_ref,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),itref,'/'
   end if
!
!  paths and files for model input and vgrid output
!
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   modelfile = outpath + modelfile
   colsumfile = outpath + 'column_sums.kim'
   outpath_ref = iter_path_ref+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   refmodelfile = outpath_ref + 'total_absdev_model.kim'
   inquire(file = colsumfile, exist = excs)
!
!  create a subfolder for regular grid output and later the extracted slices
!
   rsgpath = trim(outpath)//trim(rsgpath)
   call system('mkdir -p '//trim(rsgpath),ierr)
   write(vgridfile,'(a,a4,i2.2,a5,i2.2,a)') trim(rsgpath)//'pseudo_mesh_'//trim(prop),'_it_',it,'_itr_',itref,'.hdf'
   print *,'Model taken from :',trim(modelfile)
   print *,'Reference model taken from: ',trim(refmodelfile)
   if (excs) print *,'Column sums taken from: ',trim(colsumfile)
   print *,'Values on pseudo mesh written to: ',trim(vgridfile)
!
!  get ASKI background model
!
   call readASKIBackgroundModel(abm,abmfile,1,rearth,errmsg)
   if (.level.errmsg == 2) goto 1
   call createFromAbmKernelInvertedModel(abkim,abm,invgrid,thetac,lon_box_center,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  read model and reference model perturbations and column sums
!
   call readHDFKernelInvertedModel(kim,modelfile,errmsg)
   if (.level.errmsg == 2) goto 1
   if (excs) then
       call readHDFKernelInvertedModel(kimcs,colsumfile,errmsg)
       if (.level.errmsg == 2) goto 1
   end if
   if (itref > 0) then
      call readHDFKernelInvertedModel(krm,refmodelfile,errmsg)
      if (.level.errmsg == 2) goto 1
   end if
   if (.not. any(kim%prop == prop)) then
      call add(errmsg,2,'desired property was not inverted for',myname)
      goto 1
   end if
   iprop = findloc(kim%prop,prop,1)
   print *,'Model property ',trim(prop),' has index ',iprop,' in model'
   print *,'Number of inversion grid cells: ',kim%ncell
   if (itref > 0) then
      iprop_ref = findloc(krm%prop,prop,1)
      print *,'Model property ',trim(prop),' has index ',iprop_ref,' in reference model'
   end if
   iprop_abm = findloc(abkim%prop,prop,1)
   print *,'Model property ',trim(prop),' has index ',iprop_abm,' in ASKI background model'
!
!  generate spherical pseudo mesh from cell centers
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
   write(6,'(a,2f12.3,f15.0)') 'Adjusted min-limits of pseudo mesh: ',lonmin/mc_deg2rad,latmin/mc_deg2rad,rmin
   write(6,'(a,2f12.3,f15.0)') 'Adjusted max-limits of pseudo mesh: ',&
                               (lonmin+(nlon-1)*dlon)/mc_deg2rad,(latmin+(nlat-1)*dlat)/mc_deg2rad,rmin+(nr-1)*dr
!
   if (excs) then
      call rsg%create(4,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'absval:absdev:reldev:colsum','kim')
   else
      call rsg%create(3,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'absval:absdev:reldev','kim')
   end if
!
!  step through elements ony by one
!  here we transform from the SPECFEM Cartesian element coordinates 
!  to their spherical coordinates in the equator centered pseudo mesh
!
   allocate(mravg(nr),reldevmin(nr),reldevmax(nr))
   mravg = 0.d0; reldevmin = 1000.d0; reldevmax = -100.d0
   do icell = 1,kim%ncell
      mval = kim%model_values(icell,iprop)             ! total perturbations after it
      if (excs) colsum = kimcs%model_values(icell,1)
      if (itref > 0) then
         refval = krm%model_values(icell,iprop_ref)    ! total perturbations after itref
      else
         refval = 0.d0                                 ! or zero
      endif
      m1d = abkim%model_values(icell,iprop_abm)
      call invgrid%getSelectedCellCenter(icell,xc,yc,zc)
      r = hypot(hypot(zc+rearth,xc),yc)
      lat = dasin(yc/r)
      rp = r*dcos(lat)
      lon = dasin(xc/rp)
      call rsg%getIndicesClosest(r,lat,lon,k,j,i)
   !
   !  compute absval, absdev and reldev of model at regular grid point
   !  ic is flattened index of grid point
   !  also add to average of model deviations at given radius
   !
      ic = i+nlon*(j-1+(k-1)*nlat)
      rsg%field(ic,1) = mval+m1d
      rsg%field(ic,2) = mval-refval
      rsg%field(ic,3) = (mval-refval)/m1d
      if (excs) rsg%field(ic,4) = colsum
      mravg(k) = mravg(k)+mval-refval
      reldevmin(k) = min(rsg%field(ic,3), reldevmin(k))
      reldevmax(k) = max(rsg%field(ic,3), reldevmax(k))
   !
      if (icell == 1) write(6,'(a8,a12,2a12,3a6,a12,2a12)') 'icell','r','lat','lon','k','j','i','rg','latg','long'
      if (icell == 1 .or. icell == 21000 .or. icell == 51000 .or. icell == 87000) then
         write(6,'(i8,f12.1,2f12.6,3i6,f12.1,2f12.6)') icell,r,lat,lon,k,j,i,rmin+(k-1)*dr,latmin+(j-1)*dlat,lonmin+(i-1)*dlon
      end if
   enddo
!
!  write mean relative deviations to file and subtract it from field values
!
   write(6,'(a12,3a16)') 'Depth','Mean deviation','Min reldev','Max reldev'
   do k = 1,nr
      write(6,'(f12.1,3f16.4)') (rearth-rsg%rmin-(k-1)*rsg%dr)*1.d-3,mravg(k)/(nlat*nlon),reldevmin(k),reldevmax(k)
   end do
!
!  write regular grid to HDF file
!
   call rsg%writeHDF(trim(vgridfile),errmsg)
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   call dealloc(kim)
   if (excs) call dealloc(kimcs)
   if (itref > 0) call dealloc(krm)
   call dealloc(abm)
   call dealloc(abkim)
   call invgrid%dealloc()
   deallocate(invgrid)
   deallocate(mravg)
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

end program throwKimOnSphericalPseudoMesh
