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
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!   Module which computes integrated spectral waveform sensitivity kernels
!
!   This module computes waveform sensitivity kernels on the
!   wavefield points at the given frequencies for property sets consistent with
!   module propertySet and, on demand, integrates them onto the inversion grid
!
!   Authors: Florian Schumacher, Wolfgang Friederich
!   2015-2023
!
module spectralWaveformKernel
!
   use propertySet
   use inversionGrid
   use specfem3dKernelWavefield
   use rayKernelWavefield
   use hdfWrapper
   use realloc
   use mathConstants
   use errorMessage
   use globalHdfInfo
   use globalMpiInfo
!
   implicit none
!
   interface dealloc; module procedure deallocateSpectralWaveformKernel; end interface
   interface operator (.prop.); module procedure getPropertiesSpectralWaveformKernel; end interface
   interface operator (.nprop.); module procedure getNpropOutSpectralWaveformKernel; end interface
   interface operator (.nk.); module procedure getNkSpectralWaveformKernel; end interface
   interface operator (.df.); module procedure getDfSpectralWaveformKernel; end interface
   interface operator (.ncomp.); module procedure getNcompSpectralWaveformKernel; end interface
   interface operator (.comp.); module procedure getComponentsSpectralWaveformKernel; end interface
   interface operator (.oninvgrid.); module procedure getOnInvgridSpectralWaveformKernel; end interface
!
   type spectral_waveform_kernel
      private
      integer :: ncell                                                  ! number of PE-specific inversion cells
      integer :: ncell_all                                              ! total number of inversion grid cells
      integer :: offcell                                                ! offset of cells of this process
      integer :: nwp                                                    ! number of PE-specific wavefield points
      integer :: nwp_all                                                ! total number of wavefield points
      integer :: offwp                                                  ! offset of wavefield points of this process
      integer :: ncomp                                                  ! number of components in this kernel (i.e. force components)
      integer :: nprop_out                                              ! number of model properties for which kernels are written
      character(len=char_len_par), dimension(:), allocatable :: prop_out  ! names of model properties for which kernels are written
      integer :: nprop                                                  ! number of model properties for which kernels are required
      character(len=char_len_par), dimension(:), allocatable :: prop    ! names of model properties for which kernels are required
      character(len=char_len_comp), dimension(:), allocatable :: comp   ! names of force components in kernel
      logical :: on_invgrid = .true.                                    ! pre-integrated values on inversion grid
      double precision :: df                                            ! df
      double precision :: rsctap                                        ! residual scattering time kernel taper length
      real, dimension(:,:,:), pointer :: kernel_re => null()            ! real part of kernel values (nk,nprop,ncomp)
      real, dimension(:,:,:), pointer :: kernel_im => null()            ! imag part of kernel values (nk,nprop,ncomp)
   end type spectral_waveform_kernel
!
contains
!------------------------------------------------------------------------
!  initiate computation of spectral waveform kernel
!
   subroutine initSpectralWaveformKernel(this,ncell,ncell_all,offcell,nwp,nwp_all,offwp,prop_req,prop_out,df,rsctap,on_wp)
      type (spectral_waveform_kernel) :: this
      integer :: ncell,ncell_all,offcell,nwp,nwp_all,offwp
      character(len=char_len_par), dimension(:) ::  prop_req,prop_out
      double precision :: df,rsctap
      logical :: on_wp
   !
   !  set object members
   !
      this%ncell = ncell
      this%ncell_all = ncell_all
      this%offcell = offcell
      this%nwp = nwp
      this%nwp_all = nwp_all
      this%offwp = offwp
   !
      this%nprop = size(prop_req)                                             ! required kernels
      allocate(this%prop(this%nprop))
      this%prop = prop_req
      this%nprop_out = size(this%prop_out)
      allocate(this%prop_out(this%nprop_out))
      this%prop_out = prop_out
   !
      this%on_invgrid = .not. on_wp
      this%df = df
      this%rsctap = rsctap
   end subroutine initSpectralWaveformKernel
!------------------------------------------------------------------------
!  deallocate object
!
   subroutine deallocateMetaSpectralWaveformKernel(this)
      type (spectral_waveform_kernel) :: this
      if(allocated(this%prop)) deallocate(this%prop)
      if(allocated(this%prop_out)) deallocate(this%prop_out)
   end subroutine deallocateMetaSpectralWaveformKernel
!------------------------------------------------------------------------
!  deallocate kernel arrays
!
   subroutine deallocateKernelSpectralWaveformKernel(this)
      type (spectral_waveform_kernel) :: this
      if (allocated(this%comp)) deallocate(this%comp)
      if(associated(this%kernel_re)) deallocate(this%kernel_re)
      if(associated(this%kernel_im)) deallocate(this%kernel_im)
   end subroutine deallocateKernelSpectralWaveformKernel
!------------------------------------------------------------------------
!  deallocate everything
!
   subroutine deallocateSpectralWaveformKernel(this)
      type (spectral_waveform_kernel) :: this
      if(allocated(this%prop)) deallocate(this%prop)
      if(allocated(this%prop_out)) deallocate(this%prop_out)
      if (allocated(this%comp)) deallocate(this%comp)
      if(associated(this%kernel_re)) deallocate(this%kernel_re)
      if(associated(this%kernel_im)) deallocate(this%kernel_im)
   end subroutine deallocateSpectralWaveformKernel
!------------------------------------------------------------------------
!   Compute spectral waveform kernel at current frequency
!   and optionally pre-integrate to inversion grid
!   kwf:        kernel wavefields
!   invgrid:    inversion grid object
!   uf_mdata:   unit factor of the measured data
!   propset:    propserty set
!   comp:       force components for Green tensor
!   jf:         frequency index of current frequency
!   petsta:     phase end time of data and synthetics at station (relative to start of injection seismograms)
!   errmsg:     error message
!
   subroutine computeSpectralWaveformKernel(this,kwf,invgrid,uf_mdata,propset,comp,jf,petsta,errmsg)
      type (spectral_waveform_kernel) :: this
      class (kernel_wavefield) :: kwf
      class (inversion_grid) :: invgrid
      double precision :: uf_mdata,petsta
      type (property_set) :: propset
      character(len=char_len_comp), dimension(:)  ::  comp
      integer :: jf
      type (error_message) :: errmsg
      integer :: nwp,ncomp
      double precision, dimension(:), allocatable :: uf_equation
      double complex :: omega
      character(len=400) :: errstr
      character (len=38) :: myname = 'computeSpectralWaveformKernel'

      nwp = kwf%getNwp()
      if (nwp /= this%nwp) then
         print *,'Rank: ',myrank,' kwf%nwp: ',nwp,' invgrid%nwp: ',this%nwp
         call add(errmsg,2,'# of wavefield points in kernel wavefield inconsistent with inversion grid',myname)
         return
      end if
   !
      this%ncomp = size(comp)
      allocate(this%comp(this%ncomp))
      this%comp = comp
      ncomp = kwf%getNcomp()
      if (ncomp /= this%ncomp) then
         call add(errmsg,2,'Force components in kernel wavefield inconsistent with dmspace',myname)
         return
      end if
   !
   !  allocate space for output kernels
   !
      if (.not. this%on_invgrid) then
         allocate(this%kernel_re(this%nwp,this%nprop_out,this%ncomp))
         allocate(this%kernel_im(this%nwp,this%nprop_out,this%ncomp))
      else
         allocate(this%kernel_re(this%ncell,this%nprop_out,this%ncomp))
         allocate(this%kernel_im(this%ncell,this%nprop_out,this%ncomp))
      end if
   !
   ! Compute kernel factors which accounts for the rest of the kernel equation, outside the kernel expressions,
   ! i.e. pre-integration, values of inverted model (i.e. need a vector of factors) and measured data
   ! These factors account for (1) unit factor of measured data (dividing by this factor)
   !
      allocate(uf_equation(this%nprop))
      uf_equation = 1.d0 / uf_mdata
   !
   ! now that everything seems alright, set current frequency index in this object (can be used by other routines below)
   ! furthermore, compute the complex angular frequency accounting for any possible imaginary part that a method assumes
   !
      omega = 2.d0*mc_pid*jf*this%df
   !
      select case(trim(.setname.propset))
      case('isoVelocitySI')
         select type (kwf)
         class is (specfem3d_kernel_wavefield)
            call computeIsoVelocitySpecfem3dSpectralWaveformKernel(this,invgrid,propset,kwf,uf_equation,omega,petsta,errmsg)
         class is (ray_kernel_wavefield)
            call computeIsoVelocityRaySpectralWaveformKernel(this,invgrid,propset,kwf,uf_equation,omega,petsta,errmsg)
         class default
            call add(errmsg,2,'kernel_wavefield-type not implemented',myname)
            return
         end select
      case default
         write(errstr,*) "there are no routines implemented to compute spectral waveform kernels of property set '",&
              trim(.setname.propset),"'"
         call add(errmsg,2,errstr,myname)
         return
      end select
!
      if(allocated(uf_equation)) deallocate(uf_equation)
   end subroutine computeSpectralWaveformKernel
!------------------------------------------------------------------------
!  compute kernel values for poperty set of type 'isoVelocity'
!  this:          spectral waveform kernel
!  invgrid:       inversion grid
!  kd:            kernelDisplacement object
!  kgt:           kernelGreenTensor object
!  uf_equation:   array of unit factors, one factor for each prameter
!  omega:         angular frequency which corresponds to current frequency (can include imaginary part used by method)
!  errmsg:        error message
!
   subroutine computeIsoVelocitySpecfem3dSpectralWaveformKernel(this,invgrid,propset,kwf,uf_equation,omega,petsta,errmsg)
      type (spectral_waveform_kernel) :: this
      class (inversion_grid) :: invgrid
      type (property_set) :: propset
      class (specfem3d_kernel_wavefield) :: kwf
      double precision, dimension(:) :: uf_equation
      double complex :: omega
      double precision :: petsta
      type (error_message):: errmsg
      double complex, dimension(:,:,:), allocatable :: kernel_on_wp
      double complex, dimension(:,:), pointer :: ustr,u
      double complex, dimension(:,:,:), pointer :: gstr,g
      double precision, dimension(:), pointer :: petu,petg
   !
   ! get wavefields and strains and phase end times
   !
      u => kwf%getDisplacement()
      ustr => kwf%getStrain()
      g => kwf%getGreenTensor()
      gstr => kwf%getGreenStrain()
      petu => kwf%getDisplacementPhaseEndTime()
      petg => kwf%getGreenTensorPhaseEndTime()
   !
      call computeIsoVelocityOnWpSpectralWaveformKernel(this,u,ustr,petu,g,gstr,petg,invgrid,propset,uf_equation,&
         omega,petsta,kernel_on_wp,errmsg)
      if (.not. this%on_invgrid) then
         this%kernel_re = real(kernel_on_wp,4)
         this%kernel_im = real(aimag(kernel_on_wp),4)
      else
         call integrateSpecfem3dSpectralWaveformKernel(this,kernel_on_wp,invgrid)
      endif
   !
      deallocate(kernel_on_wp)
      if(associated(u)) nullify(u)
      if(associated(ustr)) nullify(ustr)
      if(associated(g)) nullify(g)
      if(associated(gstr)) nullify(gstr)
   !
   end subroutine computeIsoVelocitySpecfem3dSpectralWaveformKernel
!------------------------------------------------------------------------
!   given kernel values on wavefield points, do integration onto inversion grid
!   this spectral waveform kernel
!   kernel_on_wp kernel values computed on wavefield points, as returned by private routines of this module
!   invgrid inverasion grid
!   errmsg error message
!   myname name of routine which called this routine
!
   subroutine integrateSpecfem3dSpectralWaveformKernel(this,kernel_on_wp,invgrid)
      type (spectral_waveform_kernel) :: this
      double complex, dimension(:,:,:) :: kernel_on_wp
      class (inversion_grid) :: invgrid
      ! local
      integer :: icell,iprop,icomp,ip,ngll3
      double precision, dimension(:,:), pointer :: weight
      double complex :: temp
   !
      call invgrid%getIntegrationWeights(weight)
      ngll3 = invgrid%getNgll()
      do icomp = 1,this%ncomp
         do iprop = 1,this%nprop_out
            ip = 0
            do icell = 1,this%ncell
               temp = sum(weight(:,icell)*kernel_on_wp(ip+1:ip+ngll3,iprop,icomp))
               this%kernel_re(icell,iprop,icomp) = real(temp,4)
               this%kernel_im(icell,iprop,icomp) = real(aimag(temp),4)
               ip = ip+ngll3
            end do
         end do
      end do
   end subroutine integrateSpecfem3dSpectralWaveformKernel
!-------------------------------------------------------------------------------------------------------------------
!  compute kernel values for poperty set of type 'isoVelocity'
!  this:          spectral waveform kernel
!  invgrid:       inversion grid
!  kwf:           kernel wavefield object
!  uf_equation:   array of unit factors, one factor for each prameter
!  omega:         angular frequency which corresponds to current frequency (can include imaginary part used by method)
!  errmsg:        error message
!
   subroutine computeIsoVelocityRaySpectralWaveformKernel(this,invgrid,propset,kwf,uf_equation,omega,petsta,errmsg)
      type (spectral_waveform_kernel) :: this
      class (inversion_grid) :: invgrid
      type (property_set) :: propset
      class (ray_kernel_wavefield) :: kwf
      double precision, dimension(:) :: uf_equation
      double complex :: omega
      double precision :: petsta
      type (error_message):: errmsg
      integer :: js,jr,ip,ic
      double complex, dimension(:,:), pointer :: u
      double complex, dimension(:,:), pointer :: ustr
      double complex, dimension(:,:,:), pointer :: g
      double complex, dimension(:,:,:), pointer :: gstr
      double precision, dimension(:), pointer :: petu,petg
      double complex, dimension(:,:,:), allocatable :: kernel_on_wp
      double precision, dimension(:,:), pointer :: wksum
   !   character (len=43) :: myname = 'computeIsoVelocityRaySpectralWaveformKernel'
   !
      petu => kwf%getDisplacementPhaseEndTime()
      petg => kwf%getGreenTensorPhaseEndTime()
      do js = 1,kwf%getNps()
         u => getPhaseDisplacementRayKernelWavefield(kwf,js)
         ustr => getPhaseStrainRayKernelWavefield(kwf,js)
         do jr = 1,kwf%getNpr()
            g => getPhaseGreenTensorRayKernelWavefield(kwf,jr)
            gstr => getPhaseGreenStrainRayKernelWavefield(kwf,jr)
            call computeIsoVelocityOnWpSpectralWaveformKernel(this,u,ustr,petu,g,gstr,petg,invgrid,propset,uf_equation,&
                 omega,petsta,kernel_on_wp,errmsg)
            if (.level.errmsg == 2) return
            call kwf%getWavenumberSum(js,jr,real(omega),wksum)
            call invgrid%getIntegrationWeights(wksum)
            forall (ip = 1:this%nprop_out, ic = 1:this%ncomp)
               this%kernel_re(:,ip,ic) = this%kernel_re(:,ip,ic)&
                     +real(kernel_on_wp(:,ip,ic),4)*wksum(1,:)*wksum(2,:)*wksum(3,:)
               this%kernel_im(:,ip,ic) = this%kernel_im(:,ip,ic)&
                     +real(aimag(kernel_on_wp(:,ip,ic)),4)*wksum(1,:)*wksum(2,:)*wksum(3,:)
            end forall
            deallocate(kernel_on_wp)
         end do
      end do
   end subroutine computeIsoVelocityRaySpectralWaveformKernel
!-------------------------------------------------------------------------------------------------------------------
!  compute kernel values for property set of type 'isoVelocity' on wavefield points
!  this:          spectral waveform kernel
!  nwp:           number of wavefield points
!  ncomp:         number of force components for Green tensor
!  u,ustr:        source displacement and strain
!  g,gstr:        Green displacement and strain
!  invgrid:       inversion grid
!  uf_equation:   array of unit factors, one factor for each parameter
!  omega:         angular frequency which corresponds to current frequency
!  errmsg:        error message
!
   subroutine computeIsoVelocityOnWpSpectralWaveformKernel(this,u,ustr,petu,g,gstr,petg,invgrid,propset,uf_equation,&
               omega,petsta,kernel_on_wp,errmsg)
      type (spectral_waveform_kernel) :: this
      double complex, dimension(:,:) :: ustr,u
      double complex, dimension(:,:,:) :: gstr,g
      double precision, dimension(:), pointer :: petu,petg
      class (inversion_grid) :: invgrid
      type (property_set) :: propset
      double precision, dimension(:) :: uf_equation
      double precision :: petsta
      double complex :: omega
      double complex, dimension(:,:,:), allocatable :: kernel_on_wp
      type (error_message):: errmsg
      integer :: iprop,icomp,idx
      logical :: any_rho,any_vp,any_vs
      logical, dimension(:), allocatable :: rsc_in_taper
      double precision, dimension(:), allocatable :: rsc
      double precision, dimension(:), pointer :: rho,vp,vs,h_p
      double precision, dimension(:,:), pointer :: corrmat
      double complex, dimension(:,:), allocatable :: krho,klambda,kmu
      double complex, dimension(:,:), allocatable :: krhop,kvp,kvs
      character (len=44) :: myname = 'computeIsoVelocityOnWpSpectralWaveformKernel'
   !
      any_rho = any(this%prop == 'rho')
      any_vp = any(this%prop == 'vp')
      any_vs = any(this%prop == 'vs')
   !
   !  get values of properties on wavefield points
   !
      call invgrid%getWpModelValues(rho,vp,vs)
   !
   !  compute kernels
   !
      if(any_rho) then
         allocate(krho(this%nwp,this%ncomp))
         do icomp = 1,this%ncomp
            krho(:,icomp) = (omega*omega) * &
                 ( u(:,1)*g(:,1,icomp) + u(:,2)*g(:,2,icomp) + u(:,3)*g(:,3,icomp))
         end do
      else
         allocate(krho(1,1))           ! dummy allocation
      end if
   !
   !  lambda kernel is required for ALL rho-, vp- and vs-isoVelocity kernels
   !
      allocate(klambda(this%nwp,this%ncomp))
      do icomp = 1,this%ncomp
         klambda(:,icomp) = -1.d0 * (ustr(:,1)+ustr(:,2)+ustr(:,3)) * &
              (gstr(:,1,icomp)+gstr(:,2,icomp)+gstr(:,3,icomp))
      end do
   !
   !  mu kernel is required only for rho- and vs-isoVelocity kernels
   !
      if (any_rho .or. any_vs) then
         allocate(kmu(this%nwp,this%ncomp))
         do icomp = 1,this%ncomp
            kmu(:,icomp) = &
               -2.d0*(ustr(:,1)*gstr(:,1,icomp) + ustr(:,2)*gstr(:,2,icomp) + ustr(:,3)*gstr(:,3,icomp)) &
               -4.d0*(ustr(:,4)*gstr(:,4,icomp) + ustr(:,5)*gstr(:,5,icomp) + ustr(:,6)*gstr(:,6,icomp))
         end do
      else
         allocate(kmu(1,1))           ! dummy allocation
      end if
   !
   !  compute isoVelocity kernels via linearized relations between properties ro,lam,mu and ro,vp,vs
   !  do this for the required kernel properties
   !
      do iprop = 1,this%nprop
         if (equalString(this%prop(iprop),'rho')) then
            allocate(krhop(this%nwp,this%ncomp))
            do icomp = 1,this%ncomp
               krhop(:,icomp) = uf_equation(iprop) * ( krho(:,icomp) + &
                    (vp**2-2.d0*vs**2) * klambda(:,icomp) + vs**2 * kmu(:,icomp) )
            end do
         else if (equalString(this%prop(iprop),'vp')) then
            allocate(kvp(this%nwp,this%ncomp))
            do icomp = 1,this%ncomp
               kvp(:,icomp) = uf_equation(iprop) * 2.d0*rho*vp*klambda(:,icomp)
            end do
         else if (equalString(this%prop(iprop),'vs')) then
            allocate(kvs(this%nwp,this%ncomp))
            do icomp = 1,this%ncomp
               kvs(:,icomp) = uf_equation(iprop) * (2.d0*rho*vs*(kmu(:,icomp)-2.d0*klambda(:,icomp)))
            end do
         else
            call add(errmsg,2,'Invalid kernel property requested',myname)
            return
         endif
      end do
   !
   !  dummy allocation to avoid compiler warnings
   !
      if (.not. allocated(krhop)) allocate(krhop(1,1))
      if (.not. allocated(kvp)) allocate(kvp(1,1))
      if (.not. allocated(kvs)) allocate(kvs(1,1))
   !
   !  calculate the output kernels really needed for doing inversion on invgrid
   !  as determined by the property correlation matrix
   !  Krho' = a(1,1)*Krho+a(2,1)*Kvp+a(3,1)*Kvs
   !  Kvp'  = a(1,2)*Krho+a(2,2)*Kvp+a(3,2)*Kvs
   !  Kvs'  = a(1,3)*Krho+a(2,3)*Kvp+a(3,3)*Kvs
   !  with a(i,j) = b(i,j)*p_i/p_j, where p_i is any of the three parameters.
   !
   !  residual scattering time: rsc = petsta-(petu+petg)
   !
      allocate(kernel_on_wp(this%nwp,this%nprop_out,this%ncomp))
      allocate(rsc_in_taper(this%nwp),rsc(this%nwp))
      rsc = petsta-petu-petg
      rsc_in_taper = (rsc .gt. -this%rsctap .and. rsc < 0.d0)
      corrmat => getCorrmatPropertySet(propset)
      do icomp = 1,this%ncomp
         do iprop = 1,this%nprop_out
            idx = propset.index.(this%prop_out(iprop))
            if (idx == 0) then
               call add(errmsg,2,'Requested model property not found in property set',myname)
               return
            end if
            if (equalString(this%prop_out(iprop),'rho')) then
               h_p => rho
            else if (equalString(this%prop_out(iprop),'vp')) then
               h_p => vp
            else if (equalString(this%prop_out(iprop),'vs')) then
               h_p => vs
            else
               h_p => null()
               call add(errmsg,2,'Invalid model property requested',myname)
               return
            end if
         !
            kernel_on_wp(:,iprop,icomp) = 0.d0
            if (abs(corrmat(1,idx)) > 1.e-6) then      ! contribution by rho-kernel to idx-kernel
               if (idx == 1) then
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(1,idx)*krhop(:,icomp)
               else
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(1,idx)*rho/h_p*krhop(:,icomp)
               end if
            end if
            if (abs(corrmat(2,idx)) > 1.e-6) then      ! contribution by vp-kernel to idx-kernel
               if (idx == 2) then
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(2,idx)*kvp(:,icomp)
               else
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(2,idx)*vp/h_p*kvp(:,icomp)
               end if
            end if
            if (abs(corrmat(3,idx)) > 1.e-6) then      ! contribution by vs-kernel to idx-kernel
               if (idx == 3) then
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(3,idx)*kvs(:,icomp)
               else
                  kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp) + corrmat(3,idx)*vs/h_p*kvs(:,icomp)
               endif
            end if
         !
         !  restrict kernel to positive residual scattering time
         !  petu+petg < petsta                  --> take kernel as is  --> rsc > 0
         !  petsta < petu+petg < petsta+taper   --> taper kernel       --> -taper < rsc < 0
         !  petu+petg > petsta+taper            --> kernel zero        --> rsc < -taper
         !
            where(rsc < -this%rsctap) kernel_on_wp(:,iprop,icomp) = 0.d0
            where(rsc_in_taper) kernel_on_wp(:,iprop,icomp) = kernel_on_wp(:,iprop,icomp)*(1.d0+rsc/this%rsctap)
         end do
      end do
   !
      deallocate(krho,klambda,kmu)
      deallocate(krhop,kvp,kvs)
      deallocate(rsc,rsc_in_taper)
      if(associated(rho)) nullify(rho)
      if(associated(vp)) nullify(vp)
      if(associated(vs)) nullify(vs)

   end subroutine computeIsoVelocityOnWpSpectralWaveformKernel
!------------------------------------------------------------------------
!  write meta information of waveform kernel to file
!
   subroutine writeMetaSpectralWaveformKernel(this,propset,fid,jf,npath,errmsg)
      type (spectral_waveform_kernel), target :: this
      type (property_set) :: propset
      integer (kind=8) :: fid
      integer :: jf,npath
      type (error_message) :: errmsg
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer, dimension(:), allocatable :: id
      real, dimension(:), allocatable :: d
      character(len=max_length_string) :: ps_written,ps_computed
      integer :: i
   !
      call writeStringAttributeHDFWrapper(fid,'property_set_name',trim(.setname.propset),errmsg)
      if (.level.errmsg == 2) return
   !
      ps_computed = trim(this%prop(1))
      do i = 2,this%nprop
         ps_computed = trim(ps_computed)//':'//trim(this%prop(i))
      end do
      call writeStringAttributeHDFWrapper(fid,'properties_computed',trim(ps_computed),errmsg)
      if (.level.errmsg == 2) return
   !
      ps_written = trim(this%prop_out(1))
      do i = 2,this%nprop_out
         ps_written = trim(ps_written)//':'//trim(this%prop_out(i))
      end do
      call writeStringAttributeHDFWrapper(fid,'properties_written',trim(ps_written),errmsg)
      if (.level.errmsg == 2) return
   !
      if (.not. this%on_invgrid) then
         id = [this%nwp_all,jf,npath,0]
      else
         id = [this%ncell_all,jf,npath,1]
      end if
      call aria%assoc1d(id)
      call writeArrayAttributeHDFWrapper(fid,"nkall_ifreq_npath_ongrid",aria,errmsg)
      call aria%deassoc(); deallocate(id)
      if (.level.errmsg == 2) return
   !
      d = [real(this%df)]; call arra%assoc1d(d)
      call writeArrayAttributeHDFWrapper(fid,"df",arra,errmsg)
      call arra%deassoc(); deallocate(d)
      if (.level.errmsg == 2) return
   end subroutine writeMetaSpectralWaveformKernel
!------------------------------------------------------------------------
!  read meta information of waveform kernels from file
!  assumes that everything is read by one process
!  hyperslab reading is not implemented
!
   subroutine readMetaSpectralWaveformKernel(this,fid,jf,errmsg,npath,propsetname)
      type (spectral_waveform_kernel) :: this
      integer (kind=8) :: fid
      integer :: jf
      type (error_message) :: errmsg
      integer, optional :: npath
      character (len=max_length_string), optional :: propsetname
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      integer, dimension(:), pointer :: id
      integer :: i,slen
      real, dimension(:), pointer :: d
      character(len=max_length_string) :: propstring,cval
      character(len=max_length_string), dimension(:), pointer :: res
      character (len=30) :: myname = 'readMetaSpectralWaveformKernel'
   !
      if (present(propsetname)) then
         call readStringAttributeHDFWrapper(fid,'property_set_name',cval,slen,errmsg)
         if (.level.errmsg == 2) return
         if (slen > char_len_pmtrz) then
            call add(errmsg,2,'returned string length for property set name greater than assumed length',myname)
            return
         end if
         propsetname = cval(1:slen)
      end if
   !
      call readStringAttributeHDFWrapper(fid,'properties_computed',propstring,slen,errmsg)
      if (.level.errmsg == 2) return
      res => getWordsString(propstring(1:slen),':',errmsg)
      if (.level.errmsg == 2) return
      this%nprop = size(res)
      allocate(this%prop(this%nprop))
      do i = 1,this%nprop
         this%prop(i) = trim(res(i))
      end do
      deallocate(res)
   !
      call readStringAttributeHDFWrapper(fid,'properties_written',propstring,slen,errmsg)
      if (.level.errmsg == 2) return
      res => getWordsString(propstring(1:slen),':',errmsg)
      if (.level.errmsg == 2) return
      this%nprop_out = size(res)
      allocate(this%prop_out(this%nprop_out))
      do i = 1,this%nprop_out
         this%prop_out(i) = trim(res(i))
      end do
   !
      call readArrayAttributeHDFWrapper(fid,"nkall_ifreq_npath_ongrid",aria,errmsg)
      if (.level.errmsg == 2) return
      id => aria%get1d()
      this%on_invgrid = (id(4) > 0)
      if (.not. this%on_invgrid) then                 ! cell information not required
         this%nwp_all = id(1)
         this%nwp = this%nwp_all
         this%offwp = 0
         this%ncell_all = -1
         this%ncell = this%ncell_all
         this%offcell = -1
      else                                            ! wp information not required
         this%ncell_all = id(1)
         this%ncell = this%ncell_all
         this%offcell = 0
         this%nwp_all = -1
         this%nwp = this%nwp_all
         this%offwp = -1
      end if
      if (id(2) /= jf) then
         call add(errmsg,2,'desired frequency index and frequency index of waveform kernel file are inconsistent',myname)
         return
      end if
      if (present(npath)) npath = id(3)
      deallocate(id)
   !
      call readArrayAttributeHDFWrapper(fid,"df",arra,errmsg)
      if (.level.errmsg == 2) return
      d => arra%get1d()
      this%df = d(1)
      deallocate(d)
   end subroutine readMetaSpectralWaveformKernel
!------------------------------------------------------------------------
!  write values of waveform kernel to file
!
   subroutine writeSpectralWaveformKernel(this,fid,evid,netstaname,errmsg)
      type (spectral_waveform_kernel) :: this
      integer (kind=8) :: fid
      character (len=*) :: evid,netstaname
      type (error_message) :: errmsg
      ! local
      integer (kind=8), dimension(3) :: dims3d,offset3d,count3d
      integer (kind=8) :: dsetre,dsetim,dsp
      integer :: icomp,ierr,nprop
      type (any_rank_real_array) :: arra
      character(len=char_len_sta+char_len_evid+1) :: path
      character(len=max_length_string) :: compstring
      character(len=27) :: myname = 'writeSpectralWaveformKernel'
   !
      nprop = size(this%kernel_re,2)
      if (.not. this%on_invgrid) then
         dims3d = [this%nwp_all,nprop,this%ncomp]
         offset3d = [this%offwp,0,0]
         count3d = [this%nwp,nprop,this%ncomp]
      else
         dims3d = [this%ncell_all,nprop,this%ncomp]
         offset3d = [this%offcell,0,0]
         count3d = [this%ncell,nprop,this%ncomp]
      end if
   !
   !  create data sets for waveform kernels with proper dimensions
   !  and write to HDF
   !
      path = trim(evid)//'_'//trim(netstaname)
      call h5screate_simple_f(3,dims3d,dsp,ierr)
      if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
      call h5dcreate_f(fid,'kernel_real_'//trim(path),H5T_NATIVE_REAL,dsp,dsetre,ierr)
      if (ierr < 0) then; print *,'h5dcreate_f kernel_real'; goto 1; endif
      call arra%assoc3d(this%kernel_re)
      call writeArrayHDFWrapper(fid,'kernel_real_'//trim(path),arra,errmsg,xferprp = hdf_xferprp,&
                                ds = dsetre,offset = offset3d,count = count3d)
      call arra%deassoc()
      if (.level.errmsg == 2) return
   !
      call h5dcreate_f(fid,'kernel_imag_'//trim(path),H5T_NATIVE_REAL,dsp,dsetim,ierr)
      if (ierr < 0) then; print *,'h5dcreate_f kernel_imag'; goto 1; endif
      call arra%assoc3d(this%kernel_im)
      call writeArrayHDFWrapper(fid,'kernel_imag_'//trim(path),arra,errmsg,xferprp = hdf_xferprp,&
                                ds = dsetim,offset = offset3d,count = count3d)
      if (.level.errmsg == 2) return
      call arra%deassoc()
   !
      compstring = trim(this%comp(1))
      do icomp = 2,this%ncomp
         compstring = trim(compstring)//':'//trim(this%comp(icomp))
      end do
      call writeStringAttributeHDFWrapper(dsetre,'components',compstring,errmsg)
      if (.level.errmsg == 2) return
   !
      call h5dclose_f(dsetre,ierr)
      call h5dclose_f(dsetim,ierr)
1     if (ierr < 0) then
         call add(errmsg,2,'Error in HDF call',myname)
         return
      endif
   end subroutine writeSpectralWaveformKernel
! ---------------------------------------------------------------------------------
!  read values of waveform kernel from file, written for reading by one process only
!
   subroutine readSpectralWaveformKernel(this,fid,evid,netstaname,errmsg)
      type (spectral_waveform_kernel) :: this
      integer (kind=8) :: fid
      character (len=*) :: evid,netstaname
      type (error_message) :: errmsg
      integer (kind=8) :: dsetre
      integer :: i,ierr,slen
      type (any_rank_real_array) :: arra
      character(len=max_length_string) :: compstring,path
      character(len=max_length_string), dimension(:), pointer :: res
      character(len=26) :: myname = 'readSpectralWaveformKernel'
   !
   !  read data sets for waveform kernels with proper dimensions
   !
      path = trim(evid)//'_'//trim(netstaname)
      call readArrayHDFWrapper(fid,'kernel_real_'//trim(path),arra,errmsg,xferprp = hdf_xferprp,ds = dsetre)
      if (.level.errmsg == 2) return
      this%kernel_re => arra%get3d()
      call arra%deassoc()
      call readArrayHDFWrapper(fid,'kernel_imag_'//trim(path),arra,errmsg,xferprp = hdf_xferprp)
      if (.level.errmsg == 2) return
      this%kernel_im => arra%get3d()
      call arra%deassoc()
   !
      call readStringAttributeHDFWrapper(dsetre,'components',compstring,slen,errmsg)
      if (.level.errmsg == 2) return
      res => getWordsString(compstring(1:slen),':',errmsg)
      if (.level.errmsg == 2) return
      this%ncomp = size(res)
      allocate(this%comp(this%ncomp))
      do i = 1,this%ncomp
         this%comp(i) = trim(res(i))
      end do
      deallocate(res)
      call h5dclose_f(dsetre,ierr)
      if (ierr < 0) goto 1
   !
1     if (ierr < 0) then
         call add(errmsg,2,'Error in HDF call',myname)
      endif
   end subroutine readSpectralWaveformKernel
!------------------------------------------------------------------------
!  return pointer to array of property names
!
   function getPropertiesSpectralWaveformKernel(this) result(res)
      type (spectral_waveform_kernel), intent(in), target :: this
      character(len=char_len_par), dimension(:), pointer :: res
      res => this%prop_out
   end function getPropertiesSpectralWaveformKernel
!------------------------------------------------------------------------
!  return number of properties written to file
!
   function getNpropOutSpectralWaveformKernel(this) result(res)
      type (spectral_waveform_kernel), intent(in) :: this
      integer :: res
      res = this%nprop_out
    end function getNpropOutSpectralWaveformKernel
!------------------------------------------------------------------------
!  return number of proc-specific inversion grid cells of kernel
!
   function getNkSpectralWaveformKernel(this) result(res)
      type (spectral_waveform_kernel), intent(in) :: this
      integer :: res
      if (this%on_invgrid) then
         res = this%ncell
      else
         res = this%nwp
      end if
    end function getNkSpectralWaveformKernel
!------------------------------------------------------------------------
!  return df of kernel
!
   function getDfSpectralWaveformKernel(this) result(df)
      type (spectral_waveform_kernel), intent(in) :: this
      real :: df
      df = this%df
   end function getDfSpectralWaveformKernel
!------------------------------------------------------------------------
!  return number of receiver components of kernel
!
   function getNcompSpectralWaveformKernel(this) result(ncomp)
      type (spectral_waveform_kernel), intent(in) :: this
      integer :: ncomp
      ncomp = this%ncomp
   end function getNcompSpectralWaveformKernel
!------------------------------------------------------------------------
!  return pointer to array of receiver components of kernel
!
   function getComponentsSpectralWaveformKernel(this) result(comp)
      type (spectral_waveform_kernel), intent(in), target :: this
      character(len=char_len_comp), dimension(:), pointer :: comp
      comp => this%comp
   end function getComponentsSpectralWaveformKernel
!------------------------------------------------------------------------
!  return this%on_invgrid (indicating that the values of this are on inversion grid, pre-integrated)
!
   function getOnInvgridSpectralWaveformKernel(this) result(l)
      type (spectral_waveform_kernel), intent(in) :: this
      logical :: l
      l = this%on_invgrid
   end function getOnInvgridSpectralWaveformKernel
!------------------------------------------------------------------------
!  return pointer to real part of all kernel values (nk*nprop*ncomp)
!
   subroutine getAllValuesSpectralWaveformKernel(this,pre,pim)
      type (spectral_waveform_kernel) :: this
      real, dimension(:,:,:), pointer :: pre,pim
      pre => this%kernel_re
      pim => this%kernel_im
   end subroutine getAllValuesSpectralWaveformKernel
!------------------------------------------------------------------------
!  get kernel values for given component
!
   subroutine getValuesByCompSpectralWaveformKernel(this,comp,pre,pim)
      type (spectral_waveform_kernel) :: this
      character (len=*) :: comp
      real, dimension(:,:), pointer :: pre,pim
      integer :: icomp
      pre => null(); pim => null()
      do icomp = 1,this%ncomp
         if (this%comp(icomp) == comp) then
            pre => this%kernel_re(:,:,icomp)
            pim => this%kernel_im(:,:,icomp)
            return
         end if
      end do ! icomp
   end subroutine getValuesByCompSpectralWaveformKernel
!------------------------------------------------------------------------
!  return kernel values for given property
!
   subroutine getValuesByPropertySpectralWaveformKernel(this,prop,pre,pim)
      type (spectral_waveform_kernel) :: this
      character(len=*) :: prop
      real, dimension(:,:), pointer :: pre,pim
      integer :: iprop
      pre => null(); pim => null()
      do iprop = 1,this%nprop_out
         if (this%prop_out(iprop) == prop) then
            pre => this%kernel_re(:,iprop,:)
            pim => this%kernel_im(:,iprop,:)
            return
         end if
      end do
   end subroutine getValuesByPropertySpectralWaveformKernel
!------------------------------------------------------------------------
!  return kernel values for given property and component
!
   subroutine getValuesByComponentAndPropertySpectralWaveformKernel(this,comp,prop,pre,pim)
      type (spectral_waveform_kernel) :: this
      character(len=*) :: prop,comp
      real, dimension(:), pointer :: pre,pim
      integer :: i,iprop,icomp
      iprop = 0; icomp = 0
      do i = 1,this%nprop_out
         if (trim(this%prop_out(i)) == trim(prop)) iprop = i
      enddo
      do i = 1,this%ncomp
         if (trim(this%comp(i)) == trim(comp)) icomp = i
      enddo
      if (icomp == 0 .or. iprop == 0) then
         pre => null(); pim => null()
      else
         pre => this%kernel_re(:,iprop,icomp)
         pim => this%kernel_im(:,iprop,icomp)
      end if
   end subroutine getValuesByComponentAndPropertySpectralWaveformKernel
!
end module spectralWaveformKernel
