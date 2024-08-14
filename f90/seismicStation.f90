! ===============================================================================
!  Module describing a seismic station
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   ASKI version 1.2 is distributed in the hope that it
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!----------------------------------------------------------------------------
!--------------------------------------------------------------------
!> \brief  Module for seismic station information
!--------------------------------------------------------------------
 module seismicStation
    use mathConstants
    use constants
    use dateTime
    use errorMessage
    implicit none
    interface new; module procedure createSeismicStation; end interface
    interface operator (.csys.); module procedure getCoordinateSystemSeismicStation; end interface
    interface operator (.staname.); module procedure getStationNameSeismicStation; end interface
    interface operator (.nsname.); module procedure getNetstaNameSeismicStation; end interface
    interface operator (.netcode.); module procedure getNetcodeSeismicStation; end interface
    interface operator (.loc.); module procedure getLocationSeismicStation; end interface
    interface operator (.inst.); module procedure getInstnameSeismicStation; end interface
    interface operator (.lat.); module procedure getLatitudeSeismicStation; end interface
    interface operator (.lon.); module procedure getLongitudeSeismicStation; end interface
    interface operator (.latrad.); module procedure getRadianLatitudeSeismicStation; end interface
    interface operator (.lonrad.); module procedure getRadianLongitudeSeismicStation; end interface
    interface operator (.alt.); module procedure getAltitudeSeismicStation; end interface
    interface operator (.opstart.); module procedure getOperationStartSeismicStation; end interface
    interface operator (.opend.); module procedure getOperationEndSeismicStation; end interface
    type seismic_station
        private
        character (len=1) :: csys            ! coordinate system of coordinates
        character (len=char_len_netcode) :: netcode         ! network code
        character (len=char_len_sta) :: staname   ! station name
        character (len=2) :: location        ! location code
        character (len=5) :: inst_name       ! short name of instrument (e.g. CMG-3)
        double precision :: lat, lon         ! latitude and longitude (degrees) or x and y (meter) of station
        double precision :: alt              ! altitude of station (meter)
        type (date_time) :: tin              ! date and time of begin of operation
        type (date_time) :: tout             ! date and time of end of operation
    end type
!
 contains
!--------------------------------------------------------------------------------------
!> \brief Create station info object
!
    subroutine createSeismicStation(this,csys,staname,lat,lon,alt,inst_name,tin,tout,netcode,location)
    type (seismic_station) :: this
    character (len=*) :: staname,csys
    double precision :: lat,lon
    character (len = *), optional :: inst_name,netcode,location
    double precision, optional :: alt
    type (date_time), optional :: tin,tout
!
    this%csys = csys; this%lat = lat; this%lon = lon; this%staname = staname
!
    if (present(inst_name)) then; this%inst_name = inst_name; else; this%inst_name = 'undef'; endif
    if (present(netcode)) then; this%netcode = netcode; else; this%netcode = 'XX'; endif
    if (present(location)) then; this%location = location; else; this%location = ''; endif
    if (present(alt)) then; this%alt = alt; else; this%alt = 0.0; endif
    if (present(tin)) then; this%tin = tin; else; call new(this%tin,1900,1,1,0,0,0,0); endif
    if (present(tout)) then; this%tout = tout; else; call new(this%tout,2999,12,31,23,59,59,999999999); endif
    end subroutine createSeismicStation
!---------------------------------------------------------------------------
!> \brief Extend a pointer array of station info objects
!! If object contains pointers to array, the copied ones should be unlinked
!! and the no more used ones deallocated (example miniSeed.f90)
!
    function extendArraySeismicStation(array,n) result(newarray)
    type (seismic_station), dimension(:), pointer :: array
    type (seismic_station), dimension(:), pointer :: newarray
    integer :: n,nold
!
    allocate(newarray(n))
    if (.not. associated(array)) return
    nold = min(size(array),n)
    newarray(1:nold) = array(1:nold)
    deallocate(array)
    end function extendArraySeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station altitude
!
    function getAltitudeSeismicStation(this) result(alt)
    type (seismic_station), intent(in) :: this
    double precision :: alt
    alt = this%alt
    end function getAltitudeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get coordinate system
!
    function getCoordinateSystemSeismicStation(this) result(res)
    type (seismic_station), intent(in) :: this
    character (len=1) :: res
    res = this%csys
    end function getCoordinateSystemSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get instrument name
!
    function getInstnameSeismicStation(this) result(res)
    type (seismic_station), intent(in) :: this
    character (len=5) :: res
    res = this%inst_name
    end function getInstnameSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station latitude
!
    function getLatitudeSeismicStation(this) result(lat)
    type (seismic_station), intent(in) :: this
    double precision ::lat
    lat = this%lat
    end function getLatitudeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station latitude in rad
!
    function getRadianLatitudeSeismicStation(this) result(lat)
    type (seismic_station), intent(in) :: this
    double precision ::lat
    lat = this%lat*mc_deg2radd
    end function getRadianLatitudeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station longitude
!
    function getLongitudeSeismicStation(this) result(lon)
    type (seismic_station), intent(in) :: this
    double precision ::lon
    lon = this%lon
    end function getLongitudeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station longitude in rad
!
    function getRadianLongitudeSeismicStation(this) result(lon)
    type (seismic_station), intent(in) :: this
    double precision ::lon
    lon = this%lon*mc_deg2radd
    end function getRadianLongitudeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get location code
!
    function getLocationSeismicStation(this) result(location)
    type (seismic_station), intent(in) :: this
    character (len=2) :: location
    location = this%location
    end function getLocationSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get netcode
!
    function getNetcodeSeismicStation(this) result(netcode)
    type (seismic_station), intent(in) :: this
    character (len=6) :: netcode
    netcode = this%netcode
    end function getNetcodeSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get operation start time
!
    function getOperationStartSeismicStation(this) result(tin)
    type (seismic_station), intent(in) :: this
    type (date_time) :: tin
    tin = this%tin
    end function getOperationStartSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get operation end time
!
    function getOperationEndSeismicStation(this) result(tout)
    type (seismic_station), intent(in) :: this
    type (date_time) :: tout
    tout = this%tout
    end function getOperationEndSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get station name with underscores replaced by blanks
!
    function getStationNameSeismicStation(this) result(staname)
    type (seismic_station), intent(in) :: this
    character (len=char_len_sta) :: staname
    staname = this%staname
    if(char_len_sta<4) return
    if (this%staname(4:4) == '_') staname(4:4) = ' '    ! replace underscores by blanks
    if(char_len_sta<5) return
    if (this%staname(5:5) == '_') staname(5:5) = ' ' 
    end function getStationNameSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Get net.station name with underscores replaced by blanks
!
    function getNetstaNameSeismicStation(this) result(nsname)
    type (seismic_station), intent(in) :: this
    character (len=char_len_sta+char_len_netcode+1) :: nsname
    character (len=char_len_sta) :: staname
    staname = this%staname
    if (char_len_sta > 3) then
       if (this%staname(4:4) == '_') staname(4:4) = ' '    ! replace underscores by blanks
    end if
    if(char_len_sta > 4) then
       if (this%staname(5:5) == '_') staname(5:5) = ' '
    end if
    nsname = trim(this%netcode)//'.'//trim(staname)
    end function getNetstaNameSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Is station a seismometer ?
!
    logical function isSeismometerSeismicStation(this)
    type (seismic_station), intent(in) :: this
    isSeismometerSeismicStation = .false.
    select case (this%inst_name)
    case ('SEISM','CMG-3','L4-3D','STS-2','LE-3D','CMG3Q','CMG3C','CMG40','TRI40','STS-1')
        isSeismometerSeismicStation = .true.
    end select 
    end function isSeismometerSeismicStation
!-------------------------------------------------------------------------------------
!> \brief Is station a hydrophone ?
!
    logical function isHydrophoneSeismicStation(this)
    type (seismic_station), intent(in) :: this
    isHydrophoneSeismicStation = .false.
    select case (this%inst_name)
    case ('xxhyd')
        isHydrophoneSeismicStation = .true.
    end select 
    end function isHydrophoneSeismicStation
!---------------------------------------------------------------------------
!> \brief Print seismic station
!
    subroutine printSeismicStation(this)
    type (seismic_station) :: this
    write(*,*) this%staname,this%netcode,this%lat,this%lon,this%alt
    end subroutine printSeismicStation
!
 end module
