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
!   Evaluate kernel values defined on partitioned SPECFEM mesh
!   at element centers of spherical pseudo mesh defined at the equator
!   which forms a regular spherical grid
!   with spacings given by the mesh spacings and number of grid points given
!   by the mesh dimensions.
!
program throwKernelOnSphericalPseudoMesh
   use mathConstants
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use spectralWaveformKernel
   use regularSphericalGrid
   use propertySet
   use inputParameter
   use argumentParser
   use hdfWrapper
   use errorMessage
   use string
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (input_parameter), pointer :: inpar_inv
   class (inversion_grid), allocatable :: invgrid
   type (spectral_waveform_kernel) :: skernel
   type (regular_spherical_grid) :: rsg
   type (property_set), pointer :: propset
   real, dimension(:), allocatable :: d
   real, dimension(:), pointer :: pre,pim
   double precision :: dxpm,dypm,dzpm,dlon,dlat,dr,rmin,latmin,lonmin,wxpm,wypm,wzpm
   double precision :: r,lat,lon,rp,rearth,xc,yc,zc,plon,plat,lat_box_center,lon_box_center
   integer (kind=8) :: fid
   integer, dimension(3) :: griddim
   integer :: it,nr,nlat,nlon,icell,ierr,k,j,i,ic,ncell,ifreq
   character(len=max_length_string) :: main_parfile,invgrid_type,property_set_name
   character(len=max_length_string) :: main_path,iter_path,outpath,rsgpath
   character(len=max_length_string) :: kernelpath,kernelfile,rsgfile
   character(len=char_len_comp) :: comp
   character(len=char_len_par) :: prop
   character(len=char_len_evid) :: eventid
   character(len=char_len_netsta) :: nsname
   character(len=15) :: myname = 'throwKernelOnSphericalPseudoMesh'
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
   call init(ap,myname,"Throw kernel on regular spherical SPECFEM pseudo mesh")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addPosarg(ap,'eventid','sval','ID of event')
   call addPosarg(ap,'nsname','sval','Net.Station name')
   call addOption(ap,'-ifreq',.true.,'frequency index','ival','10')
   call addOption(ap,'-comp',.true.,'kernel component to be plotted','sval','UP')
   call addOption(ap,'-prop',.true.,'kernel property to be plotted','sval','vp')
   call addOption(ap,'-it',.true.,'take kernel from specified iteration','ival','0')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   eventid = ap.sval.'eventid'
   nsname = ap.sval.'nsname'
   ifreq = ap.ival.'-ifreq'
   comp = ap.sval.'-comp'
   prop = ap.sval.'-prop'
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
   property_set_name = .setname.propset
   wxpm = inpar_inv.dval.'ASKI_wx'
   wypm = inpar_inv.dval.'ASKI_wy'
   wzpm = inpar_inv.dval.'ASKI_wz'
   d = rvec(inpar_inv,'PSEUDO_MESH_SPACING',3)
   dxpm = d(1); dypm = d(2); dzpm = d(3); deallocate(d)
   griddim = ivec(inpar_inv,'INVGRID_DIMENSIONS',3)
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
!  rearth and
!  coordinates of pole the grid coordinates refer to in range (-pi,pi)
!
   rearth = invgrid%getEarthRadius()
   call invgrid%getGridCenter(lat_box_center,lon_box_center)
   plon = mc_pid+lon_box_center
   if (plon > mc_pid) plon = plon-mc_two_pid
   plat = 0.5*mc_pid-lat_box_center
!
!  paths and files for kernel input and rs-grid output
!
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   kernelpath = iter_path+(inpar_inv.sval.'PATH_SENSITIVITY_KERNELS')
   write(kernelfile,'(a,i6.6,a)') trim(kernelpath)//'spectral_kernel_'//trim(property_set_name)//'_'//'jf',ifreq,'.hdf'
!
!  create a subfolder for regular grid output and later the extracted slices
!
   rsgpath = trim(outpath)//'kernelrsg_'//eventid(1:15)//'_'//trim(nsname)//'/'
   call system('mkdir -p '//trim(rsgpath),ierr)
   write(rsgfile,'(a,a2,i3.3,a3,i3.3,a)') trim(rsgpath)//'pseudo_mesh_'//trim(prop),'_f',ifreq,'_it',it,'.hdf'
   print *,'Kernel taken from :',trim(kernelfile)
   print *,'Values on pseudo mesh written to: ',trim(rsgfile)
!
!  read kernel
!
   call openFileRoHDFWrapper(kernelfile,fid,errmsg)
   if (.level.errmsg == 2) goto 1
   call readMetaSpectralWaveformKernel(skernel,fid,ifreq,errmsg)
   if (.level.errmsg == 2) goto 1
   call readSpectralWaveformKernel(skernel,fid,eventid,nsname,errmsg)
   if (.level.errmsg == 2) goto 1
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) return
   if (.not. any(.prop.skernel == prop)) then
      call add(errmsg,2,'desired kernel property was not computed',myname)
      goto 1
   end if
   if (.not. any(.comp.skernel == comp)) then
      call add(errmsg,2,'desired kernel component was not computed',myname)
      goto 1
   end if
   ncell = invgrid%getNcellAll()
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
   call rsg%create(3,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'real:imag:abs','kernel')
!
!  step through elements ony by one
!  here we transform from the SPECFEM Cartesian element coordinates
!  to their spherical coordinates in the equator centered pseudo mesh
!
   do icell = 1,ncell
      call getValuesByComponentAndPropertySpectralWaveformKernel(skernel,comp,prop,pre,pim)
      call invgrid%getSelectedCellCenter(icell,xc,yc,zc)
      r = hypot(hypot(zc+rearth,xc),yc)
      lat = dasin(yc/r)
      rp = r*dcos(lat)
      lon = dasin(xc/rp)
      call rsg%getIndicesClosest(r,lat,lon,k,j,i)
   !
   !  compute absval, absdev and reldev of model at regular grid point
   !  ic is flattened index of grid point
   !
      ic = i+nlon*(j-1+(k-1)*nlat)
      rsg%field(ic,1) = pre(icell)
      rsg%field(ic,2) = pim(icell)
      rsg%field(ic,3) = hypot(pre(icell), pim(icell))
   !
      if (icell == 1) write(6,'(a8,a12,2a12,3a6,a12,2a12)') 'icell','r','lat','lon','k','j','i','rg','latg','long'
      if (icell == 1 .or. icell == 21000 .or. icell == 51000 .or. icell == 87000) then
         write(6,'(i8,f12.1,2f12.6,3i6,f12.1,2f12.6)') icell,r,lat,lon,k,j,i,rmin+(k-1)*dr,latmin+(j-1)*dlat,lonmin+(i-1)*dlon
      end if
   enddo
!
!  write regular grid to HDF file
!
   call rsg%writeHDF(trim(rsgfile),errmsg)
   if (.level.errmsg == 2) return
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   call dealloc(skernel)
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

end program throwKernelOnSphericalPseudoMesh
