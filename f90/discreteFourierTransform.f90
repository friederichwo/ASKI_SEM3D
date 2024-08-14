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
!   module for discrete Fourier transform
!
!  In particular developed for transform from time to frequency domain
!  in case of a sparse frequency spectrum (few frequency samples, in general
!  not equidistant). A "direct" discrete Fourier transform is applied, explicitely
!  computing coefficients efactors(omega=2*pi*f,t) = exp ( - i * omega * t ) with
!  which a Rieman-sum type of integration is conducted, (hence, very basic after all)
!
!  \author Florian Schumacher, Wolfgang Friederich
!  \date August 2013
!
module discreteFourierTransform
!
   use errorMessage
   use mathConstants
!
   implicit none
!
   interface transformForwardDFT
      module procedure dbleTransformTracesForwardDFT
      module procedure dbleTransformOneTraceForwardDFT
      module procedure singleTransformTracesForwardDFT
      module procedure singleTransformOneTraceForwardDFT
   end interface transformForwardDFT
   interface dealloc
      module procedure deallocateDFT
   end interface dealloc
   interface operator (.isdef.)
      module procedure isDefinedDFT
   end interface
!
   type discrete_fourier_transform
      private
      integer :: nf                                  ! number of frequencies
      logical :: is_defined = .false.                ! signals state of dft object (initialized with .false. upon declaration)
      double complex, dimension(:), allocatable :: f ! frequencies for which (rows of) efactors are computed
      double precision :: dt                         ! time step of time series
      integer :: nt                                  ! total number of time samples (length of t, columns of efactors)
      double precision, dimension(:), allocatable :: t  ! times for which (columns of) efactors are computed
      double complex, dimension(:,:), allocatable :: efactors  ! complex exponential coefficients for forward Fourier transform
   end type discrete_fourier_transform
!
contains
!
!-----------------------------------------------------------------------------------
!  initiate complex exponential coefficients for direct forward Fourier transform
!  this DFT object
!  dt time step of the time series which will be transformed
!  f vector of real frequency values, at which the spectrum of the Fourier transform should be evaluated
!  errmsg error message
!  taper_frac value between 0. and 1., defining a tail portion of the timeseries which are tapered by a cos hanning taper
!
   subroutine initiateForwardDFT(this,dt,nsamp,f,errmsg,taper_frac,tshift)
      type (discrete_fourier_transform) :: this
      double precision :: dt
      integer :: nsamp
      double precision, dimension(:) :: f
      type (error_message) :: errmsg
      double precision, optional :: taper_frac,tshift
      logical :: apply_taper
      double precision :: tlen,wtaper,h_tshift
      integer :: ntaper,it,jf
      character(len=18) :: myname = 'initiateForwardDFT'
   !
      if (present(taper_frac)) then
         if(taper_frac <= 0. .or. taper_frac > 1.) then
            call add(errmsg,2,'taper_frac not between 0 and 1',myname)
            return
         end if
         apply_taper = .true.
      else
         apply_taper = .false.
      endif
   !
      if (present(tshift)) then
         h_tshift = tshift
      else
         h_tshift = 0.d0
      end if
   !
      this%nt = nsamp
      this%dt = dt
      allocate(this%t(this%nt))
      this%t = (/ (dble(it-1)*this%dt, it = 1,nsamp) /)
   !
      this%nf = size(f)
      allocate(this%f(this%nf))
      this%f = f
      allocate(this%efactors(this%nf,this%nt))
      do it = 1,this%nt
         do jf = 1,this%nf
            this%efactors(jf,it) = cdexp(-mc_cid*mc_two_pid*this%f(jf)*(this%t(it)+h_tshift))*dt
         end do
      end do
   !
      if (apply_taper) then
         tlen = this%dt*dble(this%nt-1)
         wtaper = tlen*dble(taper_frac)
         ntaper = wtaper/this%dt
         if (ntaper > 0) then
            do it = this%nt-ntaper+1,this%nt
               this%efactors(:,it) = this%efactors(:,it)*0.5d0*(1.d0-dcos(mc_pid*this%dt*dble(this%nt-it)/wtaper))
            enddo
         endif
      endif
      this%is_defined = .true.
   !
   end subroutine initiateForwardDFT
! ------------------------------------------------------------------------
!  deallocate object of type discrete_fourier_transform
!  this DFT object
!
   subroutine deallocateDFT(this)
      type (discrete_fourier_transform) :: this
      this%is_defined = .false.
      if (allocated(this%t)) deallocate(this%t)
      if (allocated(this%f)) deallocate(this%f)
      if (allocated(this%efactors)) deallocate(this%efactors)
   end subroutine deallocateDFT
!------------------------------------------------------------------------
!  Transform more than one time series to frequency domain
!  Using the complex exponential coefficients defined in routine initateForwardDFT
!  each incoming trace is transformed by a matrix vector multiplication efactors*trace
!  which essentially computes the sum of the time samples weighted by respective exponential
!  coefficients, which also include (simple constant) integeration weights, such that the
!  actual Fourier integral is computed.
!  this DFT object
!  traces incoming (nt,ntrace)-array of traces
!  spectra (nf,ntrace)-array containing allocated space for the spectra of the transformed traces
!  errmsg error message
!
   subroutine dbleTransformTracesForwardDFT(this,traces,spectra,errmsg)
      type (discrete_fourier_transform) :: this
      double precision, dimension(:,:) :: traces
      double complex, dimension(:,:), allocatable :: spectra
      type (error_message) :: errmsg
      character(len=29) :: myname = 'dbleTransformTracesForwardDFT'
      integer :: ntrace,nt
!
      nt = size(traces,1)
      ntrace = size(traces,2)
      if (nt /= this%nt) then
         call add(errmsg,2,"Number of samples in traces inconsistent with DFT initiation",myname)
         return
      endif
      allocate(spectra(this%nf,ntrace))
      spectra = matmul(this%efactors,traces)
   end subroutine dbleTransformTracesForwardDFT
!------------------------------------------------------------------------
!  Transform more than one time series to frequency domain
!  Version with single precision trace and spectrum
!
   subroutine singleTransformTracesForwardDFT(this,traces,spectra,errmsg)
      type (discrete_fourier_transform) :: this
      real, dimension(:,:) :: traces
      complex, dimension(:,:), allocatable :: spectra
      type (error_message) :: errmsg
      character(len=31) :: myname = 'singleTransformTracesForwardDFT'
      integer :: ntrace,nt
!
      nt = size(traces,1)
      ntrace = size(traces,2)
      if (nt /= this%nt) then
         call add(errmsg,2,"Number of samples in traces inconsistent with DFT initiation",myname)
         return
      endif
      allocate(spectra(this%nf,ntrace))
      spectra = matmul(this%efactors,traces)
   end subroutine singleTransformTracesForwardDFT
!------------------------------------------------------------------------
!  Transform exactly one time series to frequency domain
!  Using the complex exponential coefficients defined in routine initateForwardDFT
!  the incoming trace is transformed by a matrix vector multiplication efactors*trace
!  which essentially computes the sum of the time samples weighted by respective exponential
!  coefficients, which also include (simple constant) integeration weights, such that the
!  actual Fourier integral is computed.
!  this DFT object
!  trace incoming (nt)-array containing time series
!  spectrum (nf)-array containing allocated specfe for the spectrum of the transformed trace
!  errmsg error message
!
   subroutine dbleTransformOneTraceForwardDFT(this,trace,spectrum,errmsg)
      type (discrete_fourier_transform) :: this
      double precision, dimension(:) :: trace
      double complex, dimension(:), allocatable :: spectrum
      type (error_message) :: errmsg
      character(len=31) :: myname = 'dbleTransformOneTraceForwardDFT'
      integer :: nt
!
      nt = size(trace)
      if (nt /= this%nt) then
         call add(errmsg,2,"Number of samples in trace inconsistent with DFT initiation",myname)
         return
      endif
      allocate(spectrum(this%nf))
      spectrum = matmul(this%efactors,trace)
   end subroutine dbleTransformOneTraceForwardDFT
!------------------------------------------------------------------------
!  Transform exactly one time series to frequency domain
!  Version with single precision trace and spectrum
!
   subroutine singleTransformOneTraceForwardDFT(this,trace,spectrum,errmsg)
      type (discrete_fourier_transform) :: this
      real, dimension(:) :: trace
      complex, dimension(:), allocatable :: spectrum
      type (error_message) :: errmsg
      character(len=33) :: myname = 'singleTransformOneTraceForwardDFT'
      integer :: nt
!
      nt = size(trace)
      if (nt /= this%nt) then
         call add(errmsg,2,"Number of samples in trace inconsistent with DFT initiation",myname)
         return
      endif
      allocate(spectrum(this%nf))
      spectrum = matmul(this%efactors,trace)
   end subroutine singleTransformOneTraceForwardDFT
!--------------------------------------------------------------------------
!  Function that returns if DFT object is defined
!
   function isDefinedDFT(this) result(res)
      type(discrete_fourier_transform), intent(in) :: this
      logical :: res
      res = this%is_defined
   end function isDefinedDFT
!
end module discreteFourierTransform
