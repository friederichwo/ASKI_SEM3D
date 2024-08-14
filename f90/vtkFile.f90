!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!   Module to write data living on innversion grid or wavefield points to vkt output
!
!   The inversion grid geometry (as type USTRUCTURED_GRID) or wavefield points
!   alongside with scalar data living on the inversion grid is written
!   to binary or ascii vtk files. Complex data is handled as 2 component
!   scalar float data. As an option, multiple files containing data
!   w.r.t. some index (frequency, time) may be written having the same
!   file base name followed by an index, in order to be considered by
!   Paraview as a sequence of data.
!
!  \author Florian Schumacher
!  \date Nov 2015
!
module vtkFile
!
   use inversionGrid
   use errorMessage
!
   implicit none
   interface dealloc; module procedure deallocateVtkFile; end interface
!
!  general file, geometry and cell information of vtk file
   type vtk_info
      private
      integer :: mode                                      ! either 1 = invgrid or 2 = wavefield points
      integer :: geometry_type                             ! either 0 = CELLS or 1 = CELL CENTERS
      logical is_ascii                                     ! ascii (true) or binary (false) format of vtk file
      integer :: ndata                                     ! number of data handled in this vtk file
      ! vtk points
      integer :: npoints                                   ! number of points in vtk file
      double precision, dimension(:,:), pointer :: points  ! POINTS geometry for UNSTRUCTURED_GRID
      ! vtk cells
      integer :: ncells                                    ! number of inversion grid cells handled in this vtk file
      integer, dimension(:), pointer :: cell_connectivity  ! array of point indices defining the vtk cells
      integer, dimension(:), pointer :: cell_type          ! vtk cell type for each cell (usually all equal)
   end type vtk_info
!
   contains
!------------------------------------------------------------------------
!  Initiate geometry structure of vtk file
!  Define format (ascii /binary), geometry and scale
!  and get the cell geometry from module inversion grid.
!  invgrid:            inversion grid
!  vtk_format:         'ASCII' or 'BINARY' indicating the vtk file format
!  vtk_geometry_type:  either CELLS or CELL_CENTERS
!  vtk_scale:          scaling of VTK coordinates
!  mode:               1 = inverson grid, 2 = wavefield points
!  errmsg              error message
!
   subroutine initiateVtkFile(this,invgrid,vtk_format,vtk_geometry_type,vtk_scale,mode,errmsg)
      type (vtk_info) :: this
      class (inversion_grid) :: invgrid
      character(len=*) :: vtk_format
      character(len=*) :: vtk_geometry_type
      double precision :: vtk_scale
      integer :: mode
      type (error_message) :: errmsg
      double precision, dimension(:), pointer :: x,y,z
      character(len=11) :: myname = 'initVtkFile'
   !
      if(trim(vtk_format) /= 'ASCII' .and. trim(vtk_format) /= 'BINARY') then
         call add(errmsg,2," invalid vtk format ",myname)
         return
      endif
   !
      select case(trim(vtk_geometry_type))
      case('CELLS')
         this%geometry_type = 0
      case('CELL_CENTERS')
         this%geometry_type = 1
      case default
         this%geometry_type = -1
      end select
      if (this%geometry_type < 0) then
         call add(errmsg,2,'invalid geometry type',myname)
         return
      endif
   !
      this%is_ascii = (trim(vtk_format) == 'ASCII')
      this%mode = mode
      if (mode == 1) then
         call invgrid%getGeometryVtk(this%geometry_type,vtk_scale,this%points,&
              this%cell_connectivity,this%cell_type)
         select case (this%geometry_type)
         case(0) ! VOLUMETRIC CELLS
            this%npoints = size(this%points,2)
            this%ncells = size(this%cell_type)
            this%ndata = this%ncells
         case(1) ! CELL CENTER POINTS
            this%npoints = size(this%points,2)
            this%ncells = 0
            this%ndata = this%npoints
         case default
            call add(errmsg,2,'Invalid geomtry type',myname)
            return
         end select
      else if(mode == 2) then
         call invgrid%getWavefieldPoints(x,y,z)
         this%npoints = size(x)
         this%ndata = this%npoints
         this%ncells = 0
         allocate(this%points(3,size(x)))
         this%points(1,:) = x
         this%points(2,:) = y
         this%points(3,:) = z
         nullify(x,y,z)
         this%cell_connectivity => null()
         this%cell_type => null()
      else
         call add(errmsg,2,'invalid mode, either 1 = invgrid or 2 = wavefiels points',myname)
      end if
   end subroutine initiateVtkFile
!------------------------------------------------------------------------
!  deallocate object
!
   subroutine deallocateVtkFile(this)
      type (vtk_info) :: this
      if(associated(this%points)) deallocate(this%points)
      if(associated(this%cell_connectivity)) deallocate(this%cell_connectivity)
      if(associated(this%cell_type)) deallocate(this%cell_type)
   end subroutine deallocateVtkFile
!------------------------------------------------------------------------
!  open vtk file to write
!  lu file unit
!  file_index optional index of file (will be appended to filename bas
!  errmsg error message
!
   subroutine openVtkFile(this,basename,lu,filename_extension,errmsg)
      type (vtk_info) :: this
      character (len=*) :: basename
      integer :: lu
      character (len=*) :: filename_extension
      type (error_message) :: errmsg
      character(len=11) :: myname = 'openVtkFile'
      character (len=400) :: vtk_file
      integer :: ios
   !
   ! create filename from basename and (possibly empty) filename extension plus '.vtk' extension
      vtk_file = trim(basename)//trim(filename_extension)//'.vtk'
   !
   ! open file. According to value of open_status, an existing file is overwritten
      if(this%is_ascii) then
         open(unit=lu,file=trim(vtk_file),form='FORMATTED',status='unknown',action='WRITE',iostat=ios)
         if(ios/=0) call add(errmsg,2,"could not open ascii file '"//trim(vtk_file)//"'to write",myname)
      else ! this%is_ascii
         open(unit=lu,file=trim(vtk_file),form='UNFORMATTED',access='STREAM',status='unknown',action='WRITE',&
              convert='BIG_ENDIAN',iostat=ios)
         if(ios/=0) call add(errmsg,2,"could not open binary file '"//trim(vtk_file)//"'to write",myname)
      endif ! this%is_ascii
   end subroutine openVtkFile
!------------------------------------------------------------------------
!  write vtk header, points and cell structure to open vtk file
!  lu file unit of file (MUST BE OPENED ALREADY!)
!  errmsg error message
!
   subroutine writeHeaderGeometryVtkFile(this,basename,lu,errmsg)
      type (vtk_info) :: this
      character (len=*) :: basename
      integer :: lu
      type(error_message) :: errmsg
      integer :: ios
      character(len=25) :: vtk_dataset_type
      character (len=1), parameter :: eol_char = char(10)
      logical :: err
      character(len=26) :: myname = 'writeHeaderGeometryVtkFile'
   !
      if (this%mode == 1) then
         select case(this%geometry_type)
         case(0) ! VOLUMETRIC CELLS
            vtk_dataset_type = 'DATASET UNSTRUCTURED_GRID'
         case(1) ! CELL CENTER POINTS
            vtk_dataset_type = 'DATASET POLYDATA'
         end select
      else
         vtk_dataset_type = 'DATASET POLYDATA'
      end if
   !
   ! remember with err if there was an error somewhere
      err = .false.
      if(this%is_ascii) then
         ! vkt Header
         write(unit=lu,fmt='(a)',iostat=ios) '# vtk DataFile Version 3.1'  ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(a)',iostat=ios) trim(basename)                ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(a)',iostat=ios) 'ASCII'                       ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(a)',iostat=ios) trim(vtk_dataset_type)        ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing vtk Header',myname)
            return
         endif
      else ! this%is_ascii
         ! vtk Header
         write(unit=lu,iostat=ios) '# vtk DataFile Version 3.1'//eol_char  ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) trim(basename)//eol_char                ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) 'BINARY'//eol_char                      ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) trim(vtk_dataset_type)//eol_char        ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing vtk Header',myname)
            return
         endif
      endif ! this%is_ascii
   !
      if (this%mode == 2) then
         call writePointsGeometryVtkFile(this,lu,errmsg,myname)
         call writeVerticesGeometryVtkFile(this,lu,errmsg,myname)
      else
         select case(this%geometry_type)
         case(0) ! VOLUMETRIC CELLS
            call writePointsGeometryVtkFile(this,lu,errmsg,myname)
            call writeCellsGeometryVtkFile(this,lu,errmsg,myname)
         case(1) ! CELL CENTER POINTS
            call writePointsGeometryVtkFile(this,lu,errmsg,myname)
            call writeVerticesGeometryVtkFile(this,lu,errmsg,myname)
         end select
      end if
   end subroutine writeHeaderGeometryVtkFile
!------------------------------------------------------------------------
   subroutine writePointsGeometryVtkFile(this,lu,errmsg,myname)
      type (vtk_info) :: this
      integer :: lu
      type(error_message) :: errmsg
      character(len=*) :: myname
      integer :: ios
      character (len=500) :: string
      logical :: err
      character (len=1), parameter :: eol_char = char(10)
   !
      err = .false.
      if(this%is_ascii) then
         ! POINTS
         write(unit=lu,fmt='(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'  ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(3e14.6e2)',iostat=ios) this%points                     ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing POINTS',myname)
            return
         endif
      else ! this%is_ascii
         ! POINTS
         write(string,'(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'
         write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) this%points                             ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing POINTS',myname)
            return
         endif
      end if ! this%is_ascii
   end subroutine writePointsGeometryVtkFile
!------------------------------------------------------------------------
   subroutine writeVerticesGeometryVtkFile(this,lu,errmsg,myname)
      type (vtk_info) :: this
      integer :: lu
      type(error_message) :: errmsg
      character(len=*) :: myname
      integer :: ios,i
      character (len=500) :: string
      logical :: err
      character (len=1), parameter :: eol_char = char(10)
   !
      err = .false.
      if(this%is_ascii) then
         ! VERTICES
         write(unit=lu,fmt='(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints     ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(2i12)',iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)            ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing VERTICES',myname)
            return
         endif
      else ! this%is_ascii
         ! VERTICES
         write(string,'(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints
         write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)          ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing VERTICES',myname)
            return
         endif
      end if ! this%is_ascii
   end subroutine writeVerticesGeometryVtkFile
!------------------------------------------------------------------------
   subroutine writeCellsGeometryVtkFile(this,lu,errmsg,myname)
      type (vtk_info) :: this
      integer :: lu
      type(error_message) :: errmsg
      character(len=*) :: myname
      ! logical
      integer :: ios
      character (len=500) :: string
      logical :: err
      character (len=1), parameter :: eol_char = char(10)
!
      err = .false.
      if(this%is_ascii) then
         ! CELL CONNECTIVITY
         write(unit=lu,fmt='(a,2i12)',iostat=ios) 'CELLS ',this%ncells,size(this%cell_connectivity)  ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(i12)',iostat=ios) this%cell_connectivity                                ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing CELLS',myname)
            return
         endif
         ! CELL TYPES
         write(unit=lu,fmt='(a,i12)',iostat=ios) 'CELL_TYPES ',this%ncells  ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(i12)',iostat=ios) this%cell_type               ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing CELL_TYPES',myname)
            return
         endif
      else ! this%is_ascii
         ! CELL CONNECTIVITY
         write(string,'(a,2i12)',iostat=ios)'CELLS ',this%ncells,size(this%cell_connectivity)
         write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) this%cell_connectivity                  ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing CELLS',myname)
            return
         endif
         ! CELL TYPES
         write(string,'(a,i12)',iostat=ios)'CELL_TYPES ',this%ncells
         write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) this%cell_type                          ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing CELL_TYPES',myname)
            return
         endif
      end if ! this%is_ascii
   end subroutine writeCellsGeometryVtkFile
!------------------------------------------------------------------------
!  open file, write header and geometry and one component float scalar valued cell data to vtk file
!  The number of incoming real data values must match the number this%ndata
!  of this vtk_info object, because here scalar CELL_DATA (i.e. one value per inversion grid
!  cell center point) is added to the vtk file.
!  First a file is opened calling openVtkFile and header and point and cell
!  geometry information is written to that file calling writeHeaderGeometryVtkFile.
!  Then the incoming data values are added to the vtk file as scalar valued float cell data.
!  lu:               file unit
!  data:             data values to be added to vtk file
!  data_name:        optional name of the data (by default 'data')
!  file_index:       optional index of file (will be appended to filename base)
!  fname_extension:  optional character string as file name extension IN FRONT OF file index
!  errmsg:           error message
!
   subroutine writeRealDataVtkFile(this,basename,lu,data,errmsg,data_name,file_index,fname_extension)
      ! incoming
      type (vtk_info) :: this
      character (len=*) :: basename
      integer :: lu
      real, dimension(:) :: data
      character (len=*), optional :: data_name
      integer, optional :: file_index
      character (len=*), optional :: fname_extension
      type (error_message) :: errmsg
     ! local
      character(len=10) :: vtk_data_type
      character (len=500) :: string
      character (len=500) :: filename_extension
      character (len=1), parameter :: eol_char = char(10)
      integer :: ios,ndata_expect
      logical :: err
      character(len=20) :: myname = 'writeRealDataVtkFile'
   !
      if (this%mode == 2) then
         vtk_data_type = 'POINT_DATA'
      else
         select case(this%geometry_type)
         case(0) ! VOLUMETRIC CELLS
            ! set variables for the vtk file line defining the type of data
            vtk_data_type = 'CELL_DATA'
         case(1) ! CELL CENTER POINTS
            ! set variables for the vtk file line defining the type of data
            vtk_data_type = 'POINT_DATA'
         end select
      end if
      ndata_expect = this%ndata
   !
   ! check size of incoming data
      if(size(data) /= ndata_expect) then
         call add(errmsg,2,'incoming data values do not match cells or wp points',myname)
         return
      endif
   !
      filename_extension = ''
      if(present(fname_extension)) then
         filename_extension = trim(filename_extension)//trim(fname_extension)
      endif
      if(present(file_index)) then
         write(string,"('_',i6.6)") file_index
         filename_extension = trim(filename_extension)//trim(string)
      endif
      ! open vtk file to write
      call openVtkFile(this,basename,lu,trim(filename_extension),errmsg)
      if(.level.errmsg == 2) return
   !
   !  write header and geometry information to file
      call writeHeaderGeometryVtkFile(this,basename,lu,errmsg)
      if(.level.errmsg == 2) return
   !
   !  write data to file
      err = .false.
      if(this%is_ascii) then
         write(unit=lu,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata    ; err = err.or.(ios/=0)
         if(present(data_name)) then
            string = 'SCALARS '//trim(data_name)//' float 1'
         else
            string = 'SCALARS data float 1'
         endif
         write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(e14.6e2)', iostat=ios) data                          ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing '//trim(vtk_data_type),myname)
            close(lu)
            return
         endif
      else
         write(string,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata
         write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
         if(present(data_name)) then
            string = 'SCALARS '//trim(data_name)//' float 1'
         else
            string = 'SCALARS data float 1'
         endif
         write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) data                                        ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing '//trim(vtk_data_type),myname)
            close(lu)
            return
         endif
      endif
      !
      ! close file
      close(lu)
   end subroutine writeRealDataVtkFile
!------------------------------------------------------------------------
!   open file, write header and geometry and two component float scalar valued cell data to vtk file
!   The number of incoming complex data values must match the number this%ndata
!   of this vtk_info object, because here two component scalar CELL_DATA (i.e. one complex value per inversion grid
!   cell center point) is added to the vtk file.
!   First a file is opened calling openVtkFile and header and point and cell
!   geometry information is written to that file calling writeHeaderGeometryVtkFile.
!   Then the incoming data values are added to the vtk file as two component scalar valued float cell data.
!   lu:               file unit
!   data:             data values to be added to vtk file
!   data_name:        optional name of the data (by default 'data')
!   file_index:       optional index of file (will be appended to filename base)
!   fname_extension:  optional character string as file name extension IN FRONT OF file index: this%basename//fname_extension//file_indx
!   errmsg:           error message
!
   subroutine writeComplexDataVtkFile(this,basename,lu,datare,dataim,errmsg,data_name,file_index,fname_extension)
      type (vtk_info) :: this
      character (len=*) :: basename
      integer :: lu
      real, dimension(:) :: datare,dataim
      character (len=*), optional :: data_name
      integer, optional :: file_index
      character (len=*), optional :: fname_extension
      type (error_message) :: errmsg
      ! local
      character(len=23) :: myname = 'writeComplexDataVtkFile'
      character(len=10) :: vtk_data_type
      character (len=500) :: string
      character (len=500) :: filename_extension
      character (len=1), parameter :: eol_char = char(10)
      integer :: ios,ndata_expect,i
      logical :: err
   !
      if (this%mode == 2) then
         vtk_data_type = 'POINT_DATA'
      else
         select case(this%geometry_type)
         case(0) ! VOLUMETRIC CELLS
            ! set variables for the vtk file line defining the type of data
            vtk_data_type = 'CELL_DATA'
         case(1) ! CELL CENTER POINTS
            ! set variables for the vtk file line defining the type of data
            vtk_data_type = 'POINT_DATA'
         end select
      end if

   !  by default, expect as many data as cells are handled in this vtk file
      ndata_expect = this%ndata
   !
   !  check size of incoming data
      if(size(datare) /= ndata_expect) then
         call add(errmsg,2,'incoming data values do not match invgrid points',myname)
         return
      endif
!
      filename_extension = ''
      if(present(fname_extension)) then
         filename_extension = trim(filename_extension)//trim(fname_extension)
      end if
      if(present(file_index)) then
         write(string,"('_',i6.6)") file_index
         filename_extension = trim(filename_extension)//trim(string)
      endif
   ! open vtk file to write
      call openVtkFile(this,basename,lu,trim(filename_extension),errmsg)
      if(.level.errmsg == 2) return
   !
   ! write header and geometry information to file
      call writeHeaderGeometryVtkFile(this,basename,lu,errmsg)
      if(.level.errmsg == 2) return
   !
   ! write data to file
      err = .false.
      if(this%is_ascii) then
         write(unit=lu,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata         ; err = err.or.(ios/=0)
         if(present(data_name)) then
            string = 'SCALARS '//trim(data_name)//' float 2'
         else
            string = 'SCALARS data float 2'
         endif
         write(unit=lu,fmt='(a)',iostat=ios) trim(string)                            ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'                  ; err = err.or.(ios/=0)
         write(unit=lu,fmt='(2e14.6e2)', iostat=ios) (datare(i),dataim(i),i=1,this%ndata)         ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing data',myname)
            close(lu)
            return
         endif
      else
         write(string,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata
         write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
         if(present(data_name)) then
            string = 'SCALARS '//trim(data_name)//' float 2'
         else
            string = 'SCALARS data float 2'
         endif
         write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) (datare(i),dataim(i),i=1,this%ndata)        ; err = err.or.(ios/=0)
         write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
         if(err) then ! if any of the above ios were /= 0
            call add(errmsg,2,'there was an error writing data',myname)
            close(lu)
            return
         endif
      endif
      !
      ! close file
      close(lu)
   end subroutine writeComplexDataVtkFile
!------------------------------------------------------------------------
!  open file, write only header and geometry without data to vtk file
!  lu file unit
!  errmsg error message
!
   subroutine writeHeaderVtkFile(this,basename,lu,errmsg)
      type (vtk_info) :: this
      character (len=*) :: basename
      integer :: lu
      type (error_message) :: errmsg
      character (len=1), parameter :: eol_char = char(10)
   !
   ! open vtk file to write
   !
      if (this%mode == 1) then
         call openVtkFile(this,basename,lu,'_invgrid',errmsg)
      else
         call openVtkFile(this,basename,lu,'_wp',errmsg)
      end if
      if(.level.errmsg == 2) return
    !
    ! write header and geometry information to file
      call writeHeaderGeometryVtkFile(this,basename,lu,errmsg)
      if(.level.errmsg == 2) return
   !
   ! close file
      close(lu)
   end subroutine writeHeaderVtkFile
!
end module vtkFile
