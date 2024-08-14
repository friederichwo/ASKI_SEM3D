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
program throwKernelOnInvgrid
   use mathConstants
   use inversionBasics
   use iterationStepBasics
   use propertySet
   use specfem3dInversionGrid
   use spectralWaveformKernel
   use regularSphericalGrid
   use inputParameter
   use argumentParser
   use hdfWrapper
   use errorMessage
   use smartUtils
   use globalHdfInfo
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (spectral_waveform_kernel) :: kernel
   type (property_set), pointer :: propset
   type (input_parameter), pointer :: inpar_inv
   class (inversion_grid), allocatable :: invgrid
   type (regular_spherical_grid) :: rsg
   real, dimension(:), pointer :: pre,pim
   double precision, dimension(7) :: xc,yc,zc,w,mv
   double precision :: rg,latg,long,dlon,dlat,dr,hwlon,hwlat,rmin,latmin,lonmin
   double precision :: lat_grid_center,lon_grid_center,rearth,x,y,z,plat,plon
   integer, dimension(6) :: nbidx
   integer :: it,ifreq,nr,nlat,nlon,icell,nb,ip,ierr
   integer(kind=8) :: fid
   character(len=max_length_string) :: main_parfile,invgrid_type
   character(len=max_length_string) :: main_path,iter_path,outpath,outsub
   character(len=max_length_string) :: kernelbase,kernelfile,vgridfile
   character(len=char_len_evid) :: evid
   character(len=char_len_netsta) :: nsname
   character(len=char_len_par) :: prop
   character(len=char_len_comp) :: comp
   character(len=15) :: myname = 'throwKernelOnVgrid'
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
   call init(ap,myname,"Throw kernel  on regular spherical velocity grid")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addPosarg(ap,'evid','sval','event id of the path')
   call addPosarg(ap,'netstaname','sval','net.station name of the path')
   call addOption(ap,'-hwlat',.true.,'latitude half width of grid around inversion grid center (deg)','dval','6.0')
   call addOption(ap,'-hwlon',.true.,'longitude half width of grid around inversion grid center (deg)','dval','11.0')
   call addOption(ap,'-dlat',.true.,'latitude spacing of velocity grid (deg)','dval','0.1')
   call addOption(ap,'-rmin',.true.,'bottom radius of velocity grid (km)','dval','5771')
   call addOption(ap,'-dr',.true.,'depth spacing (km)','dval','10')
   call addOption(ap,'-comp',.true.,'receiver component of the kernel','sval','UP')
   call addOption(ap,'-prop',.true.,'kernel property to be plotted','sval','vp')
   call addOption(ap,'-ifreq',.true.,'frequency index of kernel','ival','10')
   call addOption(ap,'-it',.true.,'take kernel from specified iteration','ival','0')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   evid = ap.sval.'evid'
   nsname = ap.sval.'netstaname'
   hwlat = (ap.dval.'-hwlat')*mc_deg2radd
   hwlon = (ap.dval.'-hwlon')*mc_deg2radd
   dlat = (ap.dval.'-dlat')*mc_deg2radd
   rmin = (ap.dval.'-rmin')*1.d3
   dr = (ap.dval.'-dr')*1.d3
   prop = ap.sval.'-prop'
   comp = ap.sval.'-comp'
   ifreq = ap.ival.'-ifreq'
   it = ap.ival.'-it'
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   iter_path = .iterpath.invbasics
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
   propset => getPropertySetInversionBasics(invbasics)
!
!  take model either from current iteration or from specified one
!
   if (it > 0) then
      write(iter_path,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),it,'/'
   else
      iter_path = .iterpath.invbasics
      it = inpar_inv.ival.'CURRENT_ITERATION_STEP'
   end if
!
!  paths and files for kernel input and vgrid output
!
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   kernelbase = trim(iter_path)//trim(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')//'spectral_kernel_'//trim(.setname.propset)
   write(kernelfile,'(a,i6.6,a)') trim(kernelbase)//'_jf',ifreq,'.hdf'
   write(outsub,'(a,i3.3,a)') 'kernelrsg_'//evid(1:15)//'_'//trim(nsname)//'_'//trim(comp)//'_f',ifreq,'/'
   call system('mkdir -p '//trim(outpath)//trim(outsub),ierr)
   write(vgridfile,'(a,i2.2,a)') trim(outpath)//trim(outsub)//'rsg_'//trim(prop)//'_it_',it,'.hdf'
   print *,'Kernel taken from :',trim(kernelfile)
   print *,'Vgrid file written to: ',trim(vgridfile)
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
!  read waveform sensitivity kernel
!
   call openFileRoHDFWrapper(kernelfile,fid,errmsg)
   if (.level.errmsg == 2) goto 1
   call readMetaSpectralWaveformKernel(kernel,fid,ifreq,errmsg)
   if (.level.errmsg == 2) goto 1
   call readSpectralWaveformKernel(kernel,fid,evid,nsname,errmsg)
   if (.level.errmsg == 2) goto 1
   call getValuesByComponentAndPropertySpectralWaveformKernel(kernel,comp,prop,pre,pim)
   if (.not. associated(pre)) then
      call add(errmsg,2,"no sensitivity kernel for "//trim(comp)//"_"//trim(prop)//" available",myname)
      goto 1
   end if
   call closeFileHDFWrapper(fid,errmsg)
!
!  generate vgrid and convert to SPECFEM box Cartesian coordinates
!  find closest inversion grid cell and its neighbours
!  calculate kernel values using Shepard interpolation
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
!
   call rsg%create(3,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'kernel_re:kernel_im:kernel_abs','kernel')
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
      if (mod(ip,100000) == 0) then
         write(6,'(2i8,f10.0,2f7.2,3f10.0,2e15.3)') ip,icell,rg,latg/mc_deg2radd,long/mc_deg2radd,x,y,z,hypot(pre(icell),pim(icell))
      endif
   !
   !  central cell and neighbours for interpolation
   !  pack neighbour points into one array with closest cell
   !
      call invgrid%getSelectedCellCenter(icell,xc(1),yc(1),zc(1))
      call invgrid%getSelectedFaceNeighbours(icell,nb,nbidx)
      call invgrid%getCellCentersOfFaceNeigbours(icell,nb,xc(2:),yc(2:),zc(2:))
      call shepardInterpolationSmartUtils(x,y,z,xc(1:nb+1),yc(1:nb+1),zc(1:nb+1),w(1:nb+1))
   !
   !  real part
   !
      mv(1) = pre(icell)
      mv(2:nb+1) = pre(nbidx(1:nb))
      rsg%field(ip,1) = sum(w(1:nb+1)*mv(1:nb+1))
   !
   !  imaginary part
   !
      mv(1) = pim(icell)
      mv(2:nb+1) = pim(nbidx(1:nb))
      rsg%field(ip,2) = sum(w(1:nb+1)*mv(1:nb+1))
   !
   !  absolute value
   !
      mv(1) = hypot(pre(icell),pim(icell))
      mv(2:nb+1) = hypot(pre(nbidx(1:nb)),pim(nbidx(1:nb)))
      rsg%field(ip,3) = sum(w(1:nb+1)*mv(1:nb+1))
   end do
!
!  write regular grid to HDF file
!
   call rsg%writeHDF(trim(vgridfile),errmsg)
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   call dealloc(kernel)
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
end program throwKernelOnInvgrid
