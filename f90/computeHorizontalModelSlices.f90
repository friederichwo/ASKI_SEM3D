program computeHorizontalModelSlices
   use inversionBasics
   use kernelInvertedModel
   use inputParameter
   use argumentParser
   use anyRankRealArray
   use anyRankIntegerArray
   use string
   use errorMessage
   use mpiSupport
   use smartUtils
   use globalMpiInfo
   use globalHdfInfo
   implicit none

   type (argument_parser) :: ap
   type (inversion_basics) :: invbasics
   type (kernel_inverted_model) :: kim,krm
   type (any_rank_real_array) :: arra
   type (any_rank_integer_array) :: aria
   type (mpi_support) :: mpisup
   type (error_message) :: errmsg
   type (input_parameter), pointer :: inpar_inv
   double precision, dimension(7) :: mv,refmv,w,xc,yc,zc
   double precision, dimension(:), pointer :: depths_km,dxy
   double precision :: rearth,wxpm,wypm,wzpm,dxpm,dypm,dzpm,xs,ys,zs,tol,dxs,dys
   double precision :: xtl,xtr,ytd,ytu
   real, dimension(:,:,:), allocatable :: mval,absdev,reldev,xval,yval
   real, dimension(:,:), pointer :: cell_centers
   real, dimension(:), allocatable :: d
   integer(kind=8) :: fid
   integer :: itref,ncell,nsl,imax,jmax,nb,nc,j,k,jsl,il,ir,jd,ju,i,icell,iprop,iprop_ref
   integer, dimension(:), pointer :: nbidx
   integer, dimension(:,:), pointer :: neighbours
   integer, dimension(:), allocatable :: id
   character(len=char_len_par) :: prop
   character(len=max_length_string) :: main_parfile,ext,main_path,iter_path,iter_path_ref,outpath,outpath_ref
   character(len=max_length_string) :: modelfile,refmodelfile,slicefile
   character(len=28) :: myname = 'computeHorizontalModelSlices'
! -------------------------------------------------------------------------------
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
   call init(ap,myname,"Compute horizontal (depth) slices through earth model")
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-depths',.true.,'array of depths in km of slices','dvec','150 180 210 240')
   call addOption(ap,'-dxy',.true.,'lateral spacings of slice samples in pseudo mesh in km','dvec','7.5 7.5')
   call addOption(ap,'-prop',.true.,'model property to be plotted','sval','vp')
   call addOption(ap,'-ext',.true.,'extension for output file','sval','v01')
   call addOption(ap,'-itref',.true.,'take reference model from specified iteration','ival','0')
   call parse(ap)
   if (myrank == 0) then
      call document(ap)
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap); call usage(ap)
         call add(errmsg,2,'Command line oculd not be read successfully',myname)
         goto 1
      endif
   endif
! -------------------------------------------------------------
!  get values of positional arguments
!
   main_parfile = ap.sval.'main_parfile'
   depths_km =>  getDvecArgumentParser(ap,'-depths')
   dxy => getDvecArgumentParser(ap,'-dxy')
   dxs = 1.d3*dxy(1); dys = 1.d3*dxy(2)
   prop = ap.sval.'-prop'
   ext = ap.sval.'-ext'
   itref = ap.ival.'-itref'
   deallocate(dxy)
! ------------------------------------------------------------------------
!  setup inversion basics
!
   call init(invbasics,main_parfile,1,errmsg)
   if (.level.errmsg == 2) goto 1
   iter_path = .iterpath.invbasics
   inpar_inv => getInputParameterInversionBasics(invbasics)
   main_path = inpar_inv.sval.'MAIN_PATH_INVERSION'
   outpath = iter_path+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   if (itref == 0) then
      iter_path_ref = iter_path
   else
      write(iter_path_ref,"(2a,i3.3,a)") trim(main_path),&
            trim(inpar_inv.sval.'ITERATION_STEP_PATH'),itref,'/'
   end if
   outpath_ref = iter_path_ref+(inpar_inv.sval.'PATH_OUTPUT_FILES')
   modelfile = outpath + 'inverted_abs_model.kim'
   refmodelfile = outpath_ref + 'reference_abs_model.kim'
   slicefile = trim(outpath)//'depth_slices_'//trim(prop)//'_'//trim(ext)//'.hdf'
   print *,'Model taken from ',trim(modelfile)
   print *,'Reference model taken from ',trim(refmodelfile)
   print *,'Interpolated model written to ',trim(slicefile)
!
!  pseudo-mesh properties
!
   rearth = inpar_inv.rval.'REARTH'
   wxpm = inpar_inv.dval.'ASKI_wx'
   wypm = inpar_inv.dval.'ASKI_wy'
   wzpm = inpar_inv.dval.'ASKI_wz'
   d = dvec(inpar_inv,'PSEUDO_MESH_SPACING',3)
   dxpm = d(1); dypm = d(2); dzpm = d(3); deallocate(d)
   tol = dxpm*1e-3
!
!  read model values
!
   call readHDFKernelInvertedModel(kim,modelfile,errmsg,cc=cell_centers,nb=neighbours)
   if (.level.errmsg == 2) goto 1
   call readHDFKernelInvertedModel(krm,refmodelfile,errmsg)
   if (.level.errmsg == 2) goto 1
   if (.not. any(kim%prop == prop)) then
      call add(errmsg,2,'desired property was not inverted for',myname)
      goto 1
   end if
!
!  index of property in kim and krm
!
   if (.not. any(kim%prop == prop)) then
      call add(errmsg,2,'Desired property was not inverted for',myname)
      goto 1
   end if
   if (.not. any(krm%prop == prop)) then
      call add(errmsg,2,'Desired property not in reference model',myname)
      goto 1
   end if
   iprop = findloc(kim%prop,prop,1)
   iprop_ref = findloc(krm%prop,prop,1)
   print *,'Model property ',trim(prop),' has index ',iprop,' in model'
   print *,'Model property ',trim(prop),' has index ',iprop_ref,' in reference model'
!
!  step through elements ony by one
!
   ncell = kim%ncell
   nsl = size(depths_km)
   imax = floor(wxpm/dxs)+1
   jmax = floor(wypm/dys)+1
   print *,'Dimensions of slice: ',imax,jmax
!
   allocate(mval(nsl,jmax,imax))
   allocate(absdev(nsl,jmax,imax))
   allocate(reldev(nsl,jmax,imax))
   allocate(xval(nsl,jmax,imax))
   allocate(yval(nsl,jmax,imax))
!
   do icell = 1,ncell
      nb = neighbours(1,icell)
      nbidx => neighbours(2:nb+1,icell)
      nc = nb+1
   !
      mv(1) = kim%model_values(icell,iprop)
      mv(2:nb+1) = kim%model_values(nbidx,iprop)
      refmv(1) = krm%model_values(icell,iprop_ref)                    ! reference model values
      refmv(2:nb+1) = krm%model_values(nbidx,iprop_ref)
   !
      do j = 1,nc
         if (j == 1) then
            k = icell
         else
            k = nbidx(j-1)
         end if
         xc(j) = cell_centers(1,k)
         yc(j) = cell_centers(2,k)
         zc(j) = cell_centers(3,k)
         call mapSphericalChunkToPseudoMesh(xc(j),yc(j),zc(j),rearth)
      enddo
   !
   !  identify sampling points in element and calculate interpolated values there
   !
      do jsl = 1,nsl
         zs = -depths_km(jsl)*1.d3
         if (zs-zc(1) > -0.5*dzpm .and. zs-zc(1) .le. 0.5*dzpm) then       ! element in depth slice
            xtl = xc(1)-0.5*dxpm+0.5*wxpm                                  ! cell margin coordinates shifted by 0.5*wx or 0.5*wy
            xtr = xtl+dxpm
            ytd = yc(1)-0.5*dypm+0.5*wypm
            ytu = ytd+dypm
 
            il = ceiling(xtl/dxs-tol/dxs)                           ! il*dxs >=  xtl-tol
            ir = floor(xtr/dxs+tol/dxs)                             ! ir*dxs <=  xtr+tol
            jd = ceiling(ytd/dys-tol/dys)                           ! jd*dys >=  ytd-tol
            ju = floor(ytu/dys+tol/dys)                             ! ju*dys <=  ytu+tol

            do i = il,ir                                                       ! Double loop over sampling points
               do j = jd,ju
                  xs = i*dxs-0.5*wxpm                                          ! xy-coordinates
                  ys = j*dys-0.5*wypm
                  call shepard_interpolation(xs,ys,zs,xc(1:nc),yc(1:nc),zc(1:nc),w(1:nc))
                  mval(jsl,j+1,i+1) = sum(w(1:nc)*mv(1:nc))
                  absdev(jsl,j+1,i+1) = sum(w(1:nc)*(mv(1:nc)-refmv(1:nc)))
                  reldev(jsl,j+1,i+1) = sum(w(1:nc)*(mv(1:nc)/refmv(1:nc)-1.0))
                  xval(jsl,j+1,i+1) = xs
                  yval(jsl,j+1,i+1) = ys
               end do
            end do                          ! end sampling points
         end if                             ! end in slice
      end do                                ! end slice
   end do                                   ! end icell
!
!  write HDF output file
!
   call createFileHDFWrapper(trim(slicefile),fid,errmsg)
   if (.level.errmsg == 2) return
!
   id = [imax+1,jmax+1,nsl]
   call aria%assoc1d(id)
   call writeArrayAttributeHDFWrapper(fid,"dimensions",aria,errmsg)
   call aria%deassoc(); deallocate(id)
   if (.level.errmsg == 2) return
!
   d = depths_km*1.e3
   call arra%assoc1d(d)
   call writeArrayAttributeHDFWrapper(fid,"depths",arra,errmsg)
   call aria%deassoc(); deallocate(d)
   if (.level.errmsg == 2) return
!
   call arra%assoc3d(mval)
   call writeArrayHDFWrapper(fid,'model',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) return
!
   call arra%assoc3d(absdev)
   call writeArrayHDFWrapper(fid,'absdev',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) return
!
   call arra%assoc3d(reldev)
   call writeArrayHDFWrapper(fid,'reldev',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) return
!
   call arra%assoc3d(xval)
   call writeArrayHDFWrapper(fid,'x',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) return
!
   call arra%assoc3d(yval)
   call writeArrayHDFWrapper(fid,'y',arra,errmsg,xferprp = hdf_xferprp)
   call arra%deassoc()
   if (.level.errmsg == 2) return
!
   call closeFileHDFWrapper(fid,errmsg)
   if (.level.errmsg == 2) return
!
!  clean up
!
   call dealloc(kim)
   call dealloc(krm)
   deallocate(mval,absdev,reldev,xval,yval,depths_km)
!------------------------------------------------------------------------------
!  error handling
!
1  if (.level.errmsg == 2) then
      if (myrank == 0) call print(errmsg)
      call abort(mpisup)
   endif
 !
   contains
!------------------------------------------------------------------------
!  map spherical chunk to cartesian box (MD, WF)
!  assumes that origin is at the equator at lon=0 and lat=0 and at surface
!  x,y,z: cartesian coordinates of point in spherial chunk on input
!  x,y,z: cartesian coordinates of point in cartesian box on output
!
   subroutine mapSphericalChunkToPseudoMesh(x,y,z,re)
   implicit none
   double precision :: x,y,z,re
   double precision :: r,rp,lon,lat,xs,ys,zs
 !
   xs = x; ys = y; zs = z
   r = dsqrt((zs+re)**2+xs**2+ys**2)
   lat = dasin(ys/r)
   rp = r*dcos(lat)
   lon = dasin(xs/rp)
   x = re*lon
   y = re*lat
   z = r-re
   end subroutine mapSphericalChunkToPseudoMesh
!---------------------------------------------------------------------------------
! x,y,z:   coordinates of GLL point (in)
! cc:      coordinates of interpolation anchors (3,nc) (in)
! w:       interpolation weigths (nc) (out)
!---------------------------------------------------------------------------------
   subroutine shepard_interpolation(x,y,z,xc,yc,zc,w)
      implicit none
      double precision :: x,y,z
      double precision, dimension(:) :: xc,yc,zc
      double precision, dimension(:) :: w
      ! local variables
      double precision :: r,h,sum_s
      double precision, dimension(:), allocatable :: t,s,d
      integer :: nc,i,j,imin
      double precision :: dmin,dmax
   !
   ! calculate the distances of (x,y,z) to interpolation anchors
   !
      nc = size(xc)
      w = 0.d0
      allocate(d(nc))
      d = sqrt((xc-x)**2+(yc-y)**2+(zc-z)**2)
      dmin = minval(d)
      imin = minloc(d,1)
      dmax = maxval(d)

   ! if the minimum distance to a cell is close to zero, just take the value there.
      if(dmin/dmax < 1.d-6) then
         w(imin) = 1.d0
         deallocate(d)
         return
      end if

   ! choose radius r, beyond which the influence of a cell center will be absolutely zero
   ! choose here dmax
      r = 1.5*dmax

   ! compute values s_i
      allocate(s(nc))
      do i = 1,nc
         if(d(i) < r/3.d0) then
            s(i) = 1.d0/d(i)                        ! near zero was checked above
         else
            h = d(i)/r - 1.d0
            s(i) = 27.*h*h/(4.d0*r)
         end if
      enddo

   ! define interpolation weights from s and
   ! direction factors t (interpolation function f_3 in the paper)
      sum_s = sum(s)
      allocate(t(nc))
      do i = 1,nc
         t(i) = 0.
         do j = 1,nc
            if (j == i) cycle
            h = (xc(i)-x)*(xc(j)-x) + (yc(i)-y)*(yc(j)-y) + (zc(i)-z)*(zc(j)-z)
            t(i) = t(i) + s(j)*(1.-h/(d(i)*d(j)))
         end do ! j
         t(i) = t(i)/sum_s
      end do ! i
      w = s*s*(1.d0+t)
      w = w/sum(w)
      deallocate(s,t,d)
   end subroutine shepard_interpolation
end program computeHorizontalModelSlices
