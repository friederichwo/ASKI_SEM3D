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
program throwWavefieldOnSphericalPseudoMesh
   use mathConstants
   use inversionBasics
   use iterationStepBasics
   use specfem3dInversionGrid
   use specfem3dKernelWavefield
   use rayKernelWavefield
   use regularSphericalGrid
   use inputParameter
   use argumentParser
   use hdfWrapper
   use mpiSupport
   use errorMessage
   use string
   use globalHdfInfo
!
   implicit none
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (inversion_basics) :: invbasics
   type (input_parameter), pointer :: inpar_inv
   type (mpi_support) :: mpisup
   class (inversion_grid), allocatable :: invgrid
   class (kernel_wavefield), allocatable :: kwev,kwsta
   type (regular_spherical_grid) :: rsg
   real, dimension(:), allocatable :: d
   double complex, dimension(:,:), pointer :: u
   double complex, dimension(:,:,:), pointer :: g
   double precision :: dxpm,dypm,dzpm,dlon,dlat,dr,rmin,latmin,lonmin,wxpm,wypm,wzpm
   double precision :: df,r,lat,lon,rp,rc,latc,lonc,rearth,xc,yc,zc,plon,plat,lat_box_center,lon_box_center
   integer, dimension(3) :: griddim
   integer :: it,nr,nlat,nlon,icell,ierr,kk,jj,ii,k,j,i,ic,ncell,ifreq,ip1,ip2,ip,ishift,ngll,ngll3,nwp,comp
   logical :: hres,gt
   character(len=max_length_string) :: main_parfile,invgrid_type,forward_method
   character(len=max_length_string) :: main_path,iter_path,outpath,rsgpath
   character(len=max_length_string) :: filebase,rsgfile
   character(len=char_len_evid) :: eventid
   character(len=char_len_sta+char_len_netcode+1) :: nsname
   character(len=15) :: myname = 'throwWavefieldOnSphericalPseudoMesh'
!
   call new(errmsg,myname)
!
!  initialise MPI
!
   call new(mpisup)
   myrank = .myrank.mpisup
   numtasks = .numtasks.mpisup
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
   call init(ap,myname,"Throw wavefield on regular spherical SPECFEM pseudo mesh")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addPosarg(ap,'eventid','sval','ID of event (kd) or network station name (kgt)')
   call addOption(ap,'-ifreq',.true.,'frequency index','ival','10')
   call addOption(ap,'-comp',.true.,'wavefield component to be plotted','ival','3')
   call addOption(ap,'-it',.true.,'take wavefield from specified iteration','ival','0')
   call addOption(ap,'-hres',.false.,'use 3 of 5 GLL points')
   call addOption(ap,'-gt',.false.,'extract Green tensor wavefields (force UP')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   main_parfile = ap.sval.'main_parfile'
   ifreq = ap.ival.'-ifreq'
   comp = ap.ival.'-comp'
   it = ap.ival.'-it'
   hres = ap.optset.'-hres'
   gt = ap.optset.'-gt'
   if (gt) then
      nsname = ap.sval.'eventid'
      eventid = 'none'
   else
      nsname = 'none'
      eventid = ap.sval.'eventid'
   end if
! ------------------------------------------------------------------------
!  get info from main_parfile
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   iter_path = .iterpath.invbasics
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
   forward_method = inpar_inv.sval.'FORWARD_METHOD'
   wxpm = inpar_inv.dval.'ASKI_wx'
   wypm = inpar_inv.dval.'ASKI_wy'
   wzpm = inpar_inv.dval.'ASKI_wz'
   d = rvec(inpar_inv,'PSEUDO_MESH_SPACING',3)
   if (hres) then
      dxpm = 0.5*d(1); dypm = 0.5*d(2); dzpm = 0.5*d(3); deallocate(d)
   else
      dxpm = d(1); dypm = d(2); dzpm = d(3); deallocate(d)
   end if
   griddim = ivec(inpar_inv,'INVGRID_DIMENSIONS',3)
   if (hres) griddim = 2*griddim+1
!
!  take wavefield and inversion grid either from current iteration or from specified one
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
!  paths and files for wavefield input and rs-grid output
!
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   if (gt) then
      filebase = iter_path+(inpar_inv.sval.'PATH_KERNEL_GREEN_TENSORS')
   else
      filebase = iter_path+(inpar_inv.sval.'PATH_KERNEL_DISPLACEMENTS')
   end if
!
!  create a subfolder for regular grid output and later the extracted slices
!
   if (gt) then
      rsgpath = trim(outpath)//'kgtrsg_'//trim(nsname)//'/'
   else
      rsgpath = trim(outpath)//'kdrsg_'//eventid(1:15)//'/'
   end if
   call system('mkdir -p '//trim(rsgpath),ierr)
   if (hres) then
      write(rsgfile,'(a,a2,i3.3,a3,i3.3,a)') trim(rsgpath)//'hres_pseudo_mesh','_f',ifreq,'_it',it,'.hdf'
   else
      write(rsgfile,'(a,a2,i3.3,a3,i3.3,a)') trim(rsgpath)//'pseudo_mesh','_f',ifreq,'_it',it,'.hdf'
   end if
   print *,'Wavefields taken from :',trim(filebase)
   print *,'Values on pseudo mesh written to: ',trim(rsgfile)
!
!  set the type of kernel displacement and kernel Green tensor
!  according to forward method
!
   if (equalString(forward_method,'specfem3d')) then
      allocate(specfem3d_kernel_wavefield :: kwev)
      allocate(specfem3d_kernel_wavefield :: kwsta)
   else if (equalString(forward_method,'fm3d')) then
      allocate(ray_kernel_wavefield :: kwev)
      allocate(ray_kernel_wavefield :: kwsta)
   else
      call add(errmsg,2,'Forward method '//trim(forward_method)//' not implemented',myname)
      goto 1
   end if
!
!  read wavefield
!
   df = .df.invbasics
   nwp = invgrid%getNwpAll()
   if (gt) then
      call kwsta%readGreenTensor(nwp,ifreq,df,filebase,nsname,['UP'],.false.,errmsg)
   else
      call kwev%readDisplacement(nwp,ifreq,df,filebase,eventid,.false.,errmsg)
   end if
   if (.level.errmsg == 2) goto 1
   if (gt) then
      g => kwsta%getGreenTensor()
      u => g(:,:,1)
   else
      u => kwev%getDisplacement()
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
   if (hres) then
      lonmin = -0.5*wxpm/rearth
      latmin = -0.5*wypm/rearth
      rmin = rearth-wzpm
   else
      lonmin = (-0.5*wxpm+0.5*dxpm)/rearth
      latmin = (-0.5*wypm+0.5*dypm)/rearth
      rmin = rearth-(wzpm-0.5*dzpm)
   end if
   write(6,'(a,2f12.3,f15.0)') 'Adjusted min-limits of pseudo mesh: ',lonmin/mc_deg2rad,latmin/mc_deg2rad,rmin
   write(6,'(a,2f12.3,f15.0)') 'Adjusted max-limits of pseudo mesh: ',&
                               (lonmin+(nlon-1)*dlon)/mc_deg2rad,(latmin+(nlat-1)*dlat)/mc_deg2rad,rmin+(nr-1)*dr
!
   call rsg%create(3,plat,plon,nr,nlat,nlon,dr,dlat,dlon,rmin,latmin,lonmin,rearth,'real:imag:abs','wavefield')
!
!  step through elements ony by one
!  here we transform from the SPECFEM Cartesian element/GLL coordinates
!  to their spherical coordinates in the equator centered pseudo mesh
!
   ngll3 = invgrid%getNgll()
   ngll = 5
   do icell = 1,ncell
      ishift = (icell-1)*ngll3
      call invgrid%getSelectedCellCenter(icell,xc,yc,zc)
      rc = hypot(hypot(zc+rearth,xc),yc)
      latc = dasin(yc/rc)
      rp = rc*dcos(latc)
      lonc = dasin(xc/rp)
   !
   !  go through every second GLL-point in each dimension
   !  touching the 27 points at the corners, centers of the edges and center of the element
   !  point count from bottom left front running through x, y, and z in that order.
   !  ip = ngll**2*(k-1)+(j-1)*ngll+i
   !
      if (hres) then
         do kk = 1,ngll,2
            r = rc-dr+(kk-1)*dr*0.5
            ip1 = ishift+ngll**2*(kk-1)
            do jj = 1,ngll,2
               lat = latc-dlat+(jj-1)*dlat*0.5
               ip2 = ip1+(jj-1)*ngll
               do ii = 1,ngll,2
                  lon = lonc-dlon+(ii-1)*dlon*0.5
                  ip = ip2+ii
                  call rsg%getIndicesClosest(r,lat,lon,k,j,i)
                  ic = i+nlon*(j-1+(k-1)*nlat)
                  rsg%field(ic,1) = real(u(ip,comp),4)
                  rsg%field(ic,2) = real(aimag(u(ip,comp)),4)
                  rsg%field(ic,3) = real(abs(u(ip,comp)),4)
               end do
            end do
         end do
      else
         call rsg%getIndicesClosest(rc,latc,lonc,k,j,i)
         ip = ishift+2*ngll**2+2*ngll+3
         ic = i+nlon*(j-1+(k-1)*nlat)
         rsg%field(ic,1) = real(u(ip,comp),4)
         rsg%field(ic,2) = real(aimag(u(ip,comp)),4)
         rsg%field(ic,3) = real(abs(u(ip,comp)),4)
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
   call kwev%dealloc()
   call kwsta%dealloc()
   deallocate(kwev,kwsta)
   call invgrid%dealloc()
   deallocate(invgrid)
   call rsg%dealloc()
   call dealloc(invbasics)
   call dealloc(errmsg)
   call dealloc(ap)
   call dealloc(mpisup)
!------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   endif

end program throwWavefieldOnSphericalPseudoMesh
