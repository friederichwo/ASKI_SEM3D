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
program throwKimOnVgrid
   use mathConstants
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use kernelInvertedModel
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
   type (kernel_inverted_model) :: kim,krm
   type (regular_spherical_grid) :: rsg
   double precision, dimension(7) :: xc,yc,zc,w,mv,refmv
   double precision :: rg,latg,long,dlon,dlat,dr,hwlon,hwlat,rmin,latmin,lonmin
   double precision :: lat_grid_center,lon_grid_center,rearth,x,y,z,plat,plon,mval,refval
   integer, dimension(6) :: nbidx
   integer :: it,itref,iprop,iprop_ref,nr,nlat,nlon,icell,nb,ip,ierr
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: main_path,iter_path,outpath,iter_path_ref,outpath_ref,rsgpath
   character(len=max_length_string) :: modelfile,refmodelfile,vgridfile
   character(len=char_len_par) :: prop
   character(len=15) :: myname = 'throwKimOnVgrid'
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
   call init(ap,myname,"Throw kernel inverted model on regular spherical velocity grid")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-hwlat',.true.,'latitude half width of grid around inversion grid center (deg)','dval','6.0')
   call addOption(ap,'-hwlon',.true.,'longitude half width of grid around inversion grid center (deg)','dval','11.0')
   call addOption(ap,'-dlat',.true.,'latitude spacing of velocity grid (deg)','dval','0.1')
   call addOption(ap,'-rmin',.true.,'bottom radius of velocity grid (km)','dval','5771')
   call addOption(ap,'-dr',.true.,'depth spacing (km)','dval','10')
   call addOption(ap,'-prop',.true.,'model property to be plotted','sval','vp')
   call addOption(ap,'-it',.true.,'take model from specified iteration','ival','0')
   call addOption(ap,'-itref',.true.,'take reference model from specified iteration','ival','1')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   hwlat = (ap.dval.'-hwlat')*mc_deg2radd
   hwlon = (ap.dval.'-hwlon')*mc_deg2radd
   dlat = (ap.dval.'-dlat')*mc_deg2radd
   rmin = (ap.dval.'-rmin')*1.d3
   dr = (ap.dval.'-dr')*1.d3
   prop = ap.sval.'-prop'
   it = ap.ival.'-it'
   itref = ap.ival.'-itref'
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   iter_path = .iterpath.invbasics
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
!
!  take model and inversion grid either from current iteration or from specified one
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
   rearth = invgrid%getEarthRadius()
   call invgrid%getGridCenter(lat_grid_center,lon_grid_center)
!
!  create locator and check rmin
!
   select type (invgrid)
      class is (specfem3d_inversion_grid)
         call invgrid%createLocator(errmsg)
         if (.level.errmsg == 2) return
         if (invgrid%isBelow(rmin-rearth,errmsg)) then
            call add(errmsg,2,'minimum radius of spherical grid is below deepest inversion grid level',myname)
            print *,'rmin = ',rmin
            goto 1
         endif
   end select
!
!  take reference model either from same iteration as model or from specified one
!  default is the reference model of the first iteration (i.e. 1D reference model)
!
   if (itref == 0) then
      iter_path_ref = iter_path
      itref = it
   else
      write(iter_path_ref,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),itref,'/'
   end if
!
!  paths and files for model input and vgrid output
!
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   outpath_ref = iter_path_ref+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   modelfile = outpath + 'inverted_abs_model.kim'
   refmodelfile = outpath_ref + 'reference_abs_model.kim'
!
!  create a subfolder for regular grid output and later the extracted slices
!
   rsgpath = trim(outpath)//'kimrsg/'
   call system('mkdir -p '//trim(rsgpath),ierr)
   write(vgridfile,'(a,a4,i2.2,a5,i2.2,a)') trim(rsgpath)//'rsg_'//trim(prop),'_it_',it,'_itr_',itref,'.hdf'
   print *,'Model taken from :',trim(modelfile)
   print *,'Reference model taken from: ',trim(refmodelfile)
   print *,'Vgrid file written to: ',trim(vgridfile)
!
!  read model and reference model values
!
   call readHDFKernelInvertedModel(kim,modelfile,errmsg)
   if (.level.errmsg == 2) goto 1
   call readHDFKernelInvertedModel(krm,refmodelfile,errmsg)
   if (.level.errmsg == 2) goto 1
   if (.not. any(kim%prop == prop)) then
      call add(errmsg,2,'desired property was not inverted for',myname)
      goto 1
   end if
   iprop = findloc(kim%prop,prop,1)
   iprop_ref = findloc(krm%prop,prop,1)
   print *,'Model property ',trim(prop),' has index ',iprop,' in model'
   print *,'Model property ',trim(prop),' has index ',iprop_ref,' in reference model'
!
!  generate vgrid and convert to SPECFEM box Cartesian coordinates
!  find closest inversion grid cell and its neighbours
!  calculate model values using Shepard interpolation
!
   dlon = dlat/cos(lat_grid_center)                        ! choose dlon to get same arclength along latitude circle
   nr = ceiling((rearth-rmin)/dr)+1
   nlat = 2.*ceiling(hwlat/dlat)+1
   nlon = 2.*ceiling(hwlon/dlon)+1
   write(6,'(a,3i8)') 'Dimensions of vgrid: ',nlon,nlat,nr
   latmin = lat_grid_center-0.5*(nlat-1)*dlat
   lonmin = lon_grid_center-0.5*(nlon-1)*dlon
   rmin = rearth-(nr-1)*dr
   plon = 0.d0
   plat = 0.5*mc_pid
   write(6,'(a,2f12.3,f15.0)') 'Adjusted min-limits of vgrid: ',lonmin/mc_deg2rad,latmin/mc_deg2rad,rmin
   write(6,'(a,2f12.3,f15.0)') 'Adjusted max-limits of vgrid: ',&
                               (lonmin+(nlon-1)*dlon)/mc_deg2rad,(latmin+(nlat-1)*dlat)/mc_deg2rad,rmin+(nr-1)*dr
!
   call rsg%create(3,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'absval:absdev:reldev','vgrid')
   do while (rsg%nextPoint(ip,rg,latg,long))
   !
   !  conversion to SPECFEM box centered Cartesian coordinates
   !  and finding of cell index
   !
      select type (invgrid)
         class is (specfem3d_inversion_grid)
            call invgrid%convertSphericalCoordinates(rg,latg,long,x,y,z)
            call invgrid%findCell(x,y,z,icell,errmsg)
            if (.level.errmsg == 2) goto 1
      end select
   !
   !  central cell and neighbours for interpolation
   !  pack neighbour points into one array with closest cell
   !
      call invgrid%getSelectedCellCenter(icell,xc(1),yc(1),zc(1))
      call invgrid%getSelectedFaceNeighbours(icell,nb,nbidx)
      call invgrid%getCellCentersOfFaceNeigbours(icell,nb,xc(2:),yc(2:),zc(2:))
      call shepardInterpolationSmartUtils(x,y,z,xc(1:nb+1),yc(1:nb+1),zc(1:nb+1),w(1:nb+1))
   !
   !  get model and reference model values
   !
      mv(1) = kim%model_values(icell,iprop)
      mv(2:nb+1) = kim%model_values(nbidx(1:nb),iprop)
      refmv(1) = krm%model_values(icell,iprop_ref)
      refmv(2:nb+1) = krm%model_values(nbidx(1:nb),iprop_ref)
      mval = sum(w(1:nb+1)*mv(1:nb+1))
      refval = sum(w(1:nb+1)*refmv(1:nb+1))
      if (mod(ip,100000) == 0) write(6,'(i12,2f12.1,f12.1,2f15.3)') ip,mval,refval,rg,latg/mc_deg2rad,long/mc_deg2rad
   !
   !  compute absval, absdev and reldev of model at regular grid point
   !
      rsg%field(ip,1) = mval
      rsg%field(ip,2) = mval-refval
      rsg%field(ip,3) = (mval-refval)/refval
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
   call dealloc(krm)
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
!
end program throwKimOnVgrid
