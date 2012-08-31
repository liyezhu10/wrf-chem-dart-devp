! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Created by Rafael Rosolem (09/30/2011) for COSMOS based on file by Ave Arellano

! BEGIN DART PREPROCESS KIND LIST
! COSMOS_NEUTRON_INTENSITY,    KIND_NEUTRON_INTENSITY
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_COSMOS_mod, only : read_neutron_intensity, &
!                                 write_neutron_intensity, &
!                          get_expected_neutron_intensity, &
!                           set_obs_def_neutron_intensity, &
!                           interactive_neutron_intensity
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! TJH Questions for Nancy: AFAICT none of this should be public ... 

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(COSMOS_NEUTRON_INTENSITY)
!         call get_expected_neutron_intensity(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(COSMOS_NEUTRON_INTENSITY)
!         call read_neutron_intensity(obs_def%key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(COSMOS_NEUTRON_INTENSITY)
!         call write_neutron_intensity(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(COSMOS_NEUTRON_INTENSITY)
!         call interactive_neutron_intensity(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS SET_OBS_DEF_NEUTRON_INTENSITY
!      case(COSMOS_NEUTRON_INTENSITY)
!         call set_obs_def_neutron_intensity(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_NEUTRON_INTENSITY


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_COSMOS_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             logfileunit, get_unit, open_file, close_file
use     location_mod, only : location_type, set_location, get_location, &
                             vert_is_height,   VERTISHEIGHT,            &
                             vert_is_level,    VERTISLEVEL
use        model_mod, only : model_interpolate
use     obs_kind_mod, only : KIND_GEOPOTENTIAL_HEIGHT, KIND_SOIL_MOISTURE

implicit none
private

public :: set_obs_def_neutron_intensity,  &
          get_expected_neutron_intensity, &
          interactive_neutron_intensity,  &
          read_neutron_intensity,         &
          write_neutron_intensity

integer, parameter :: nlayers        = 4 ! Number of soil layers used in Noah model
integer            :: num_neutron_intensity = 0 ! observation counter

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision $", &
   revdate  = "$Date$"

character(len=256) :: string1
logical, save      :: module_initialized = .false.

!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module
!
! a DART tradition
call register_module(source, revision, revdate)

module_initialized = .true.

end subroutine initialize_module


 subroutine read_neutron_intensity(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_neutron_intensity(key, ifile, fform)
!
integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
integer           :: keyin
integer           :: counts1
character(len=32) :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
!      read(ifile) keyin
   CASE DEFAULT
!      read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key     = counts1
call set_obs_def_neutron_intensity(key)

end subroutine read_neutron_intensity



 subroutine write_neutron_intensity(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_neutron_intensity(key, ifile, fform)

integer,          intent(in)           :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)                      :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
!         write(ifile) key
   CASE DEFAULT
!         write(ifile, *) key
END SELECT

end subroutine write_neutron_intensity



 subroutine get_expected_neutron_intensity(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_neutron_intensity(state, location, key, val, istatus)
! uses a weighting function calculated by COSMIC (COsmic-ray Soil
! Moisture Interaction Code)

real(r8),            intent(in)  :: state(:) ! state vector
type(location_type), intent(in)  :: location ! location of obs
integer,             intent(in)  :: key      ! obs key
real(r8),            intent(out) :: val      ! value of obs
integer,             intent(out) :: istatus  ! status of the calculation

!========================================================================================
! COsmic-ray Soil Moisture Interaction Code (COSMIC) - Version 1.5
!
! W. James Shuttleworth and Rafael Rosolem - January/2012
! Fortran code developed by Rafael Rosolem
!========================================================================================
! Updates:
! 01/20/2012 - Version 1.0: * Original version based on SPAM
! 01/28/2012 - Version 1.1: * Some parameters are re-defined for better physical realism
! 02/17/2012 - Version 1.2: * After contribution from all angles are taken, need to
!                             multiply fastflux by 2/PI
!                           * Angle increments can be specified here (ideg)
!                             resolution
! 02/29/2012 - Version 1.3: * Reduced number of parameters based on relationship of ns
!                             and nw (now given as alpha = ns/nw)
! 04/03/2012 - Version 1.4: * Soil thickness (i.e., input.dat file) needs to be specified
!                             at finer resolution (i.e., 0.1 cm)
! 04/04/2012 - Version 1.5  * Now the contributions to soil and water densities/mass are
!                             taken to be at the center of a given soil layer
!=================================================================================
! COSMIC: Variables list
!=================================================================================

!rr: Use 1 mm layer increments (down to 3 meters) to compute the weighting function
!rr: (as originally done for COSMIC), and then compute the cumulative weights for
!rr: individual soil layers

integer, parameter :: nlyr=3000 ! Total number of soil layers - each 1mm for 3m

real(r8) :: bd     = 0.0_r8 ! Dry soil bulk density (g/m3)
real(r8) :: vwclat = 0.0_r8 ! Volumetric "lattice" water content (m3/m3)
real(r8) :: N      = 0.0_r8 ! High energy neutron flux (-)
real(r8) :: alpha  = 0.0_r8 ! Ratio of Fast Neutron Creation Factor (Soil to Water), alpha (-)
real(r8) :: L1     = 0.0_r8 ! High Energy Soil Attenuation Length (g/cm2)
real(r8) :: L2     = 0.0_r8 ! High Energy Water Attenuation Length (g/cm2)
real(r8) :: L3     = 0.0_r8 ! Fast Neutron Soil Attenuation Length (g/cm2)
real(r8) :: L4     = 0.0_r8 ! Fast Neutron Water Attenuation Length (g/cm2)
real(r8) :: zdeg
real(r8) :: zrad
real(r8) :: ideg
real(r8) :: costheta
real(r8) :: dtheta
real(r8) :: totflux     ! Total flux of above-ground fast neutrons

real(r8), dimension(:), allocatable :: dz          ! Soil layers (cm)
real(r8), dimension(:), allocatable :: zthick      ! Soil layer thickness (cm)
real(r8), dimension(:), allocatable :: vwc         ! Volumetric Water Content (m3/m3)
real(r8), dimension(:), allocatable :: isoimass    ! Integrated dry soil mass above layer (g)
real(r8), dimension(:), allocatable :: iwatmass    ! Integrated water mass above layer (g)
real(r8), dimension(:), allocatable :: hiflux      ! High energy neutron flux
real(r8), dimension(:), allocatable :: fastpot     ! Fast neutron source strength of layer
real(r8), dimension(:), allocatable :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
real(r8), dimension(:), allocatable :: idegrad     ! Integrated neutron degradation factor (-)
real(r8), dimension(:), allocatable :: fastflux    ! Contribution to above-ground neutron flux

!rr: Not needed for DART
!rr: real(r8), dimension(:), allocatable :: normfast ! Normalized contribution to neutron flux (-) [weighting factors]

real(r8), parameter   :: h2odens = 1000.0_r8 ! Density of water (g/cm3)

integer,  parameter   :: maxlayers = 1000 ! more than the maximum # of model soil layers
real(r8), allocatable :: layerz(:)        ! original soil layer depths
real(r8), allocatable :: soil_moisture(:) ! original soil layer moistures

integer  :: angle, angledz, maxangle  ! loop indices for an integration interval
integer  :: i, zi, nlevels, iunit
real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat, obs_val
type(location_type) :: loc

!=================================================================================

if ( .not. module_initialized ) call initialize_module

val = 0.0_r8 ! set return value early

!=================================================================================
! Determine the number of soil layers and their depths
! using only the standard DART interfaces to model
! set the locations for each of the model levels

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

nlevels = 0
COUNTLEVELS : do i = 1,maxlayers 
   loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)
   call model_interpolate(state,loc,KIND_GEOPOTENTIAL_HEIGHT,obs_val,istatus)
   if (istatus /= 0) exit COUNTLEVELS
   nlevels = nlevels + 1
enddo COUNTLEVELS

if (nlevels == maxlayers) then
   write(string1,*) 'FAILED to determine number of soil layers in model.'
   call error_handler(E_ERR,'get_expected_neutron_intensity',string1,source,revision,revdate)
else
   write(*,*)'get_expected_neutron_intensity: we have ',nlevels,' soil layers'
endif

!=================================================================================
! Now actually find the depths at each level
! While we're at it, might as well get the soil moisture there, too.

allocate(layerz(nlevels),soil_moisture(nlevels))

FINDLEVELS : do i = 1,nlevels 
   loc = set_location(loc_lon, loc_lat, real(i,r8), VERTISLEVEL)
   call model_interpolate(state,loc,KIND_GEOPOTENTIAL_HEIGHT,layerz(i),istatus)
   ! not checking this error code because it worked just a few lines earlier

   loc = set_location(loc_lon, loc_lat, layerz(i), VERTISHEIGHT)
   call model_interpolate(state,loc,KIND_SOIL_MOISTURE,obs_val,istatus)

   if (istatus /= 0) then
      write(string1,*) 'FAILED to determine soil moisture for layer',i
      call error_handler(E_ERR,'get_expected_neutron_intensity',string1,source,revision,revdate)
   endif

   soil_moisture(i) = obs_val 
   
enddo FINDLEVELS

!rr: DART soil layers are negative and given in meters while COSMIC needs them to be
!rr: positive and in centimeters
if     ( all(layerz < 0.0_r8) ) then
   layerz(:) = -1.0_r8 * 100.0_r8 * layerz(:)
elseif ( any(layerz < 0.0_r8) )then
   write(string1,*) 'unusual values of soil layers in model.'
   call error_handler(E_ERR,'get_expected_neutron_intensity',string1,source,revision,revdate)
endif

!=================================================================================
! COSMIC: Site specific-parameters
!=================================================================================

!rr: #####
!rr: I will need to change find a way to read those parameters externally
!rr: #####

! SANTA RITA SITE
bd     = 1.4620_r8
vwclat = 0.0366_r8
N      = 399.05099359_r8
alpha  = 0.25985853017_r8
L1     = 161.98621864_r8
L2     = 129.14558985_r8
L3     = 114.04156391_r8
L4     = 3.8086191933_r8

!=================================================================================
! COSMIC: Allocate arrays and initialize variables
!=================================================================================

allocate(dz(nlyr), zthick(nlyr), vwc(nlyr), &
         hiflux(nlyr), fastpot(nlyr), h2oeffdens(nlyr), &
        idegrad(nlyr), fastflux(nlyr), &
        isoimass(nlyr), iwatmass(nlyr))

totflux = 0.0_r8

do i = 1,nlyr

   dz(i)         = real(i,r8)/10.0_r8 ! 0.1 cm intervals ... for now.
   zthick(i)     = 0.0_r8
   h2oeffdens(i) = 0.0_r8
   vwc(i)        = 0.0_r8
   isoimass(i)   = 0.0_r8
   iwatmass(i)   = 0.0_r8
   hiflux(i)     = 0.0_r8
   fastpot(i)    = 0.0_r8
   idegrad(i)    = 0.0_r8
   fastflux(i)   = 0.0_r8

    !rr: Not needed for DART
    !rr: normfast(i)    = 0.0_r8

enddo

!=================================================================================
! Get soil moisture from individual model layers and assign them to
! 1 mm intervals (down to 3 meters)

do i = 1,nlyr

   SOIL : do zi = 1,nlevels
      if (dz(i) >= layerz(nlevels)) then
         vwc(i) = soil_moisture(nlevels)
         exit SOIL
      elseif (dz(i) <= layerz(zi)) then
         vwc(i) = soil_moisture(zi) ! soil moisture (m3/m3)
         exit SOIL
      endif
   enddo SOIL

enddo

deallocate(layerz ,soil_moisture)

!=================================================================================
! COSMIC: Neutron flux calculation
!=================================================================================

! At some point, you might want to tinker around with non-uniform 
! soil layer thicknesses.
zthick(1) = dz(1) - 0.0_r8 ! Surface layer
do i = 2,nlyr
   zthick(i) = dz(i) - dz(i-1) ! Remaining layers
enddo

if ( 2 == 1 ) then ! TJH DEBUG OUTPUT BLOCK 
   iunit = open_file('cosmos_layers.txt',form='formatted',action='write')
   do i = 1,nlyr
      write(iunit,*)dz(i),vwc(i),zthick(i)
   enddo
   call close_file(iunit)
endif

! Angle distribution parameters (HARDWIRED)
!rr: Using 0.5 deg angle intervals appears to be sufficient
!rr: (smaller angles increase the computing time for COSMIC)
ideg     = 0.5_r8                   ! ideg ultimately controls the number of trips through
angledz  = nint(ideg*10.0_r8)       ! the ANGLE loop. Make sure the 10.0_r8 is enough
maxangle = 900 - angledz            ! to create integers with no remainder
dtheta   = ideg*(PI/180.0_r8)

if ( real(angledz,r8) /= ideg*10.0_r8 ) then
   write(string1,*) 'ideg*10.0 must result in an integer - it results in ',ideg*10.0_r8
   call error_handler(E_ERR,'get_expected_neutron_intensity',string1,source,revision,revdate)
endif

do i = 1,nlyr

   ! High energy neutron downward flux
   ! The integration is now performed at the node of each layer (i.e., center of the layer)

   h2oeffdens(i) = ((vwc(i)+vwclat)*h2odens)/1000.0_r8

   if(i > 1) then
      ! Assuming an area of 1 cm2
      isoimass(i) = isoimass(i-1) + bd*(0.5_r8*zthick(i-1))*1.0_r8 + &
                                    bd*(0.5_r8*zthick(i  ))*1.0_r8
      ! Assuming an area of 1 cm2
      iwatmass(i) = iwatmass(i-1) + h2oeffdens(i-1)*(0.5_r8*zthick(i-1))*1.0_r8 + &
                                    h2oeffdens(i  )*(0.5_r8*zthick(i  ))*1.0_r8
   else
      isoimass(i) =            bd*(0.5_r8*zthick(i))*1.0_r8 ! Assuming an area of 1 cm2
      iwatmass(i) = h2oeffdens(i)*(0.5_r8*zthick(i))*1.0_r8 ! Assuming an area of 1 cm2
   endif

   hiflux( i) = N*exp(-(isoimass(i)/L1 + iwatmass(i)/L2) )
   fastpot(i) = zthick(i)*hiflux(i)*(alpha*bd + h2oeffdens(i))

   ! This second loop needs to be done for the distribution of angles for fast neutron release
   ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
   ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.  

   do angle=0,maxangle,angledz
      zdeg     = real(angle,r8)/10.0_r8   ! 0.0  0.5  1.0  1.5 ...
      zrad     = (zdeg*PI)/180.0_r8
      costheta = cos(zrad)

      ! Angle-dependent low energy (fast) neutron upward flux
      fastflux(i) = fastflux(i) + fastpot(i)*exp(-(isoimass(i)/L3 + iwatmass(i)/L4)/costheta)*dtheta
   enddo

   ! After contribution from all directions are taken into account,
   ! need to multiply fastflux by 2/PI

   fastflux(i) = (2.0_r8/PI)*fastflux(i)

   ! Low energy (fast) neutron upward flux
   totflux = totflux + fastflux(i)

enddo

!=================================================================================

! ... and finally calculated what the neutron intensity (which is basically totflux)
val = totflux

deallocate(dz, zthick, vwc, hiflux, fastpot, &
           h2oeffdens, idegrad, fastflux, isoimass, iwatmass)

! assume all is well for now
istatus = 0

return
end subroutine get_expected_neutron_intensity



 subroutine interactive_neutron_intensity(key)
!----------------------------------------------------------------------
!subroutine interactive_neutron_intensity(key)
!
integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! Increment the index
num_neutron_intensity = num_neutron_intensity + 1
key = num_neutron_intensity

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_COSMOS observation'

end subroutine interactive_neutron_intensity



 subroutine set_obs_def_neutron_intensity(key)
!----------------------------------------------------------------------
! Allows passing of obs_def special information

integer, intent(in) :: key

if ( .not. module_initialized ) call initialize_module

end subroutine set_obs_def_neutron_intensity


end module obs_def_COSMOS_mod

! END DART PREPROCESS MODULE CODE

