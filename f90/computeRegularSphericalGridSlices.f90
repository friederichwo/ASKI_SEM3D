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
!   Compute slices from a regular spherical grid
!   Slices are defined in a parameter file
!-------------------------------------------------------------------------------
program computeRegularSphericalGridSlices
!
   use regularSphericalGrid
   use axesRotation
   use hdfWrapper
   use inputParameter
   use argumentParser
   use errorMessage
   use string
   use globalHdfInfo
!
   implicit none
   type (input_parameter) :: slpar
   type (argument_parser) :: ap
   type (error_message) :: errmsg
   type (regular_spherical_grid) :: rsg,hslice,vslice,gc
   integer :: nhsl,nvsl,i
   integer, dimension(:), pointer :: slnlat,slnlon,slnr,slndel
   double precision, dimension(:), pointer :: lat1,lon1,lat2,lon2,slrad
   character(len=500), dimension(:), pointer :: slname
   character(len=max_length_string) :: slicefile,rsgfile,filename
   character(len=33) :: myname = 'computeRegularSphericalGridSlices'
   character (len=80), dimension(15) :: slpar_keys
   data slpar_keys/'HORIZONTAL_SLICES', 'HS_NUM', 'HS_NAME', 'HS_RAD', 'HS_NLAT', 'HS_NLON', &
                   'VERTICAL_SLICES',  'VS_NUM', 'VS_NAME', 'VS_NR', 'VS_NDEL', &
                   'VS_LAT1', 'VS_LAT2', 'VS_LON1', 'VS_LON2'/
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
!
!  command line processing
!
   call init(ap,myname,"Compute some slices through a regular spherical grid")
   call addPosarg(ap,'rsgfile','sval','name of RSG file')
   call addOption(ap,'-slicedefs',.true.,'Text file defining slices','sval','kim_slicedef.txt')
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap); call usage(ap)
      call add(errmsg,2,'Command line could not be read successfully',myname)
      goto 1
   endif
!
   rsgfile = ap.sval.'rsgfile'
   slicefile = ap.sval.'-slicedefs'
   print *,'Read 3D regular spherical grid from ',trim(rsgfile)
!
!  read regular grid file
!
   call readHDFRegularSphericalGrid(rsg,trim(rsgfile),errmsg)
   if (.level.errmsg == 2) goto 1
!
!  read slice parameter file
!
   call createKeywordsInputParameter(slpar,slpar_keys)
   call readSubroutineInputParameter(slpar,1,trim(slicefile),errmsg)
   if (.level.errmsg == 2) return
   !
   !  horizontal slices first
   !
   nhsl = slpar.ival.'HS_NUM'
   if ((slpar.lval.'HORIZONTAL_SLICES') .and. nhsl > 0) then
      slname => svecp(slpar,'HS_NAME',nhsl)
      slrad => dvecp(slpar,'HS_RAD',nhsl)
      slnlat => ivecp(slpar,'HS_NLAT',nhsl)
      slnlon => ivecp(slpar,'HS_NLON',nhsl)
      do i = 1,nhsl
         call rsg%extractHorizontalSlice(slrad(i),slnlat(i),slnlon(i),hslice,slname(i))
         filename = getAheadLastSeparatorString(rsgfile,'.',errmsg)
         call appendString(filename,'_hslice_'//trim(slname(i))//'.hdf',errmsg)
         if (.level.errmsg == 2) goto 1
         print *,'write slice to: ',trim(filename)
         call writeHDFRegularSphericalGrid(hslice,filename,errmsg)
         if (.level.errmsg == 2) goto 1
         call hslice%dealloc()
      end do
      deallocate(slname,slrad,slnlat,slnlon)
   endif
   !
   !  vertical slices second
   !
   nvsl = slpar.ival.'VS_NUM'
   if ((slpar.lval.'VERTICAL_SLICES') .and. nvsl > 0) then
      slname => svecp(slpar,'VS_NAME',nvsl)
      slnr => ivecp(slpar,'VS_NR',nvsl)
      slndel => ivecp(slpar,'VS_NDEL',nvsl)
      lat1 => dvecp(slpar,'VS_LAT1',nvsl)
      lat2 => dvecp(slpar,'VS_LAT2',nvsl)
      lon1 => dvecp(slpar,'VS_LON1',nvsl)
      lon2 => dvecp(slpar,'VS_LON2',nvsl)
      do i = 1,nvsl
         call rsg%extractVerticalSlice(lat1(i)*mc_deg2radd,lon1(i)*mc_deg2radd,lat2(i)*mc_deg2radd,lon2(i)*mc_deg2radd,&
                                       slnr(i),slndel(i),slname(i),vslice,gc)
         filename = getAheadLastSeparatorString(rsgfile,'.',errmsg)
         call appendString(filename,'_vslice_'//trim(slname(i))//'.hdf',errmsg)
         if (.level.errmsg == 2) goto 1
         print *,'write slice to: ',trim(filename)
         call writeHDFRegularSphericalGrid(vslice,filename,errmsg)
         if (.level.errmsg == 2) goto 1

         filename = getAheadLastSeparatorString(rsgfile,'.',errmsg)
         call appendString(filename,'_great_circle_'//trim(slname(i))//'.hdf',errmsg)
         if (.level.errmsg == 2) goto 1
         print *,'write slice to: ',trim(filename)
         call writeHDFRegularSphericalGrid(gc,filename,errmsg)
         if (.level.errmsg == 2) goto 1
         call vslice%dealloc()
         call gc%dealloc()
      end do
      deallocate(slname,slnr,slndel,lat1,lat2,lon1,lon2)
   endif
!
!  clean up
!
   call closeEnvironmentHDFWrapper(errmsg)
   call rsg%dealloc()
   call dealloc(errmsg)
   call dealloc(ap)
!
!  error treatment
!
1  if (.level.errmsg == 2) then
      call print(errmsg)
   end if
!
end program computeRegularSphericalGridSlices
