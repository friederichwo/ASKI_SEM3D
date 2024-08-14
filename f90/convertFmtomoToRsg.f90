!----------------------------------------------------------------------------
!   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!   Convert data in fmtomo format into a HDF regular spherical grid file
!
program convertFmtomoToRsg
   use regularSphericalGrid
   use fmtomoModel
   use argumentParser
   use string
   use errorMessage
   use globalHdfInfo

   implicit none

   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (fmtomo_model) :: fmtomo,fmtomoref
   type (regular_spherical_grid) :: rsg
   integer :: nrnew,nabsurf
   double precision :: rearth,rmax
   double precision, dimension(:), allocatable :: vel,velref
   character(len=max_length_string) :: vgridfile,vgridreffile,rsgfile,field_name
   character(len=18) :: myname = 'convertFmtomoToRsg'
!-------------------------------------------------------------------------------
!  open HDF environment
!
   call new(errmsg,myname)
   call openEnvironmentHDFWrapper(errmsg)
   if (.level.errmsg == 2) goto 1
   call setXferprpIndependentHDFWrapper(hdf_xferprp,errmsg)
   if (.level.errmsg == 2) goto 1
!------------------------------------------------------------------------
!  command line processing
!
   call init(ap,myname,"Convert a FMTOMO 3D model to a RSG HDF file")
   call addPosarg(ap,"vgrid","sval","name of fmtomo 3D model")
   call addPosarg(ap,"vgridref","sval","name of fmtomo 1D reference model")
   call addPosarg(ap,"rsgfile","sval","Name of RSG otput file")
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then; call print(.errmsg.ap); call usage(ap); goto 1; end if
   call document(ap)
!
   vgridfile = ap.sval.'vgrid'
   vgridreffile = ap.sval.'vgridref'
   rsgfile = ap.sval.'rsgfile'
   call dealloc(ap)
!
!  read FMOTMO model
!
   print *,'Take fmtomo model from file: ',trim(vgridfile)
   print *,'Take fmtomo reference model from file: ',trim(vgridreffile)
   print *,'Write rsg to file: ',trim(rsgfile)
   call readFmtomoModel(fmtomo,1,vgridfile,errmsg)
   if (.level.errmsg == 2) goto 1
   call readFmtomoModel(fmtomoref,1,vgridreffile,errmsg)
   if (.level.errmsg == 2) goto 1
!
!  create regular spherical grid
!  remove points baove earth's surface
!
   field_name = 'reldev'
   rearth = 6371000.d0
   rmax = fmtomo%rmin+(fmtomo%nr-1)*fmtomo%dr
   nabsurf = nint((rmax-rearth)/fmtomo%dr)
   nrnew = fmtomo%nr-nabsurf
   call rsg%create(1,0.5*mc_pid,0.d0,nrnew,fmtomo%nlat,fmtomo%nlon,fmtomo%dr,fmtomo%dlat,fmtomo%dlon,&
                   fmtomo%rmin,fmtomo%latmin,fmtomo%lonmin,6371000.d0,field_name,'fmtomo')
!
!  fill field array of RSG
!
   vel = reshape(fmtomo%vel(:,:,1:nrnew),[rsg%npts])
   velref = reshape(fmtomoref%vel(:,:,1:nrnew),[rsg%npts])
   rsg%field(:,1) = (vel-velref)/velref
!
!  write to HDF
!
   call rsg%writeHDF(rsgfile,errmsg)
   if (.level.errmsg == 2) goto 1
!
   deallocate(vel,velref)
!
!  error handling
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   end if
!
end program convertFmtomoToRsg
