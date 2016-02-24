! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Fortran has a limit of 32 characters for variable names. Hence,
! each column can be at most 32 characters wide.
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! BEGIN DART PREPROCESS KIND LIST
! SAT_TEMPERATURE,                 KIND_TEMPERATURE,                COMMON_CODE
! SAT_TEMPERATURE_ELECTRON,        KIND_TEMPERATURE_ELECTRON,       COMMON_CODE
! SAT_TEMPERATURE_ION,             KIND_TEMPERATURE_ION,            COMMON_CODE
! SAT_DENSITY_NEUTRAL_O3P,         KIND_DENSITY_NEUTRAL_O3P,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_O2,          KIND_DENSITY_NEUTRAL_O2,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2,          KIND_DENSITY_NEUTRAL_N2,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N4S,         KIND_DENSITY_NEUTRAL_N4S,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_NO,          KIND_DENSITY_NEUTRAL_NO,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2D,         KIND_DENSITY_NEUTRAL_N2D,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2P,         KIND_DENSITY_NEUTRAL_N2P,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_H,           KIND_DENSITY_NEUTRAL_H,          COMMON_CODE
! SAT_DENSITY_NEUTRAL_HE,          KIND_DENSITY_NEUTRAL_HE,         COMMON_CODE
! SAT_DENSITY_NEUTRAL_CO2,         KIND_DENSITY_NEUTRAL_CO2,        COMMON_CODE
! SAT_DENSITY_NEUTRAL_O1D,         KIND_DENSITY_NEUTRAL_O1D,        COMMON_CODE
! SAT_DENSITY_ION_O4SP,            KIND_DENSITY_ION_O4SP,           COMMON_CODE
! SAT_DENSITY_ION_O2P,             KIND_DENSITY_ION_O2P,            COMMON_CODE
! SAT_DENSITY_ION_N2P,             KIND_DENSITY_ION_N2P,            COMMON_CODE
! SAT_DENSITY_ION_NP,              KIND_DENSITY_ION_NP,             COMMON_CODE
! SAT_DENSITY_ION_NOP,             KIND_DENSITY_ION_NOP,            COMMON_CODE
! SAT_DENSITY_ION_O2DP,            KIND_DENSITY_ION_O2DP,           COMMON_CODE
! SAT_DENSITY_ION_O2PP,            KIND_DENSITY_ION_O2PP,           COMMON_CODE
! SAT_DENSITY_ION_HP,              KIND_DENSITY_ION_HP,             COMMON_CODE
! SAT_DENSITY_ION_HEP,             KIND_DENSITY_ION_HEP,            COMMON_CODE
! SAT_DENSITY_ION_E,               KIND_DENSITY_ION_E,              COMMON_CODE
! SAT_VELOCITY_U,                  KIND_VELOCITY_U,                 COMMON_CODE
! SAT_VELOCITY_V,                  KIND_VELOCITY_V,                 COMMON_CODE
! SAT_VELOCITY_W,                  KIND_VELOCITY_W,                 COMMON_CODE
! SAT_VELOCITY_U_ION,              KIND_VELOCITY_U_ION,             COMMON_CODE
! SAT_VELOCITY_V_ION,              KIND_VELOCITY_V_ION,             COMMON_CODE
! SAT_VELOCITY_W_ION,              KIND_VELOCITY_W_ION,             COMMON_CODE
! SAT_VELOCITY_VERTICAL_O3P,       KIND_VELOCITY_VERTICAL_O3P,      COMMON_CODE
! SAT_VELOCITY_VERTICAL_O2,        KIND_VELOCITY_VERTICAL_O2,       COMMON_CODE
! SAT_VELOCITY_VERTICAL_N2,        KIND_VELOCITY_VERTICAL_N2,       COMMON_CODE
! SAT_VELOCITY_VERTICAL_N4S,       KIND_VELOCITY_VERTICAL_N4S,      COMMON_CODE
! SAT_VELOCITY_VERTICAL_NO,        KIND_VELOCITY_VERTICAL_NO,       COMMON_CODE
! SAT_F107,                        KIND_1D_PARAMETER,               COMMON_CODE
! SAT_RHO,                         KIND_DENSITY
! GPS_PROFILE,                     KIND_ELECTRON_DENSITY,           COMMON_CODE
! GND_GPS_VTEC,		           KIND_GND_GPS_VTEC
! CHAMP_DENSITY,                   KIND_DENSITY
! MIDAS_TEC,                       KIND_VERTICAL_TEC
! SSUSI_O_N2_RATIO,                KIND_O_N2_COLUMN_DENSITY_RATIO
! GPS_VTEC_EXTRAP,                 KIND_VERTICAL_TEC,               COMMON_CODE
! SUPER_DARN_ELECTRIC_POTENTIAL,   KIND_ELECTRIC_POTENTIAL,         COMMON_CODE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_upper_atm_mod, only : get_expected_upper_atm_density
!  use obs_def_upper_atm_mod, only : get_expected_gnd_gps_vtec
!  use obs_def_upper_atm_mod, only : get_expected_vtec
!  use obs_def_upper_atm_mod, only : get_expected_O_N2_ratio
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(SAT_RHO) 
!      call get_expected_upper_atm_density(state_handle, ens_size, location, expected_obs, istatus)
! case(CHAMP_DENSITY) 
!      call get_expected_upper_atm_density(state_handle, ens_size, location, expected_obs, istatus)
! case(MIDAS_TEC) 
!      call get_expected_vtec(state_handle, ens_size, location, expected_obs, istatus)
! case(GND_GPS_VTEC)
!      call get_expected_gnd_gps_vtec(state_handle, ens_size, location, expected_obs, istatus)
! case(SSUSI_O_N2_RATIO)
!      call get_expected_O_N2_ratio(state_handle, ens_size, location, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(SAT_RHO) 
!      continue
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!      continue
! case(GND_GPS_VTEC)
!      continue
! case(SSUSI_O_N2_RATIO)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_upper_atm_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_TEMPERATURE, &
                             KIND_PRESSURE, &
                             KIND_DENSITY, &
                             KIND_DENSITY_ION_E, &
                             KIND_GND_GPS_VTEC, &
                             KIND_GEOPOTENTIAL_HEIGHT, &
                             KIND_GEOMETRIC_HEIGHT, &
                             KIND_O_N2_COLUMN_DENSITY_RATIO
use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private
public :: get_expected_upper_atm_density, &
          get_expected_gnd_gps_vtec, &
          get_expected_vtec, &
          get_expected_O_N2_ratio

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

real(r8), PARAMETER :: N2_molar_mass = 28.0_r8
real(r8), PARAMETER :: O_molar_mass  = 16.0_r8
real(r8), PARAMETER :: O2_molar_mass = 32.0_r8
real(r8), PARAMETER :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
integer,  PARAMETER :: MAXLEVELS = 100 ! more than max levels expected in the model 

character(len=512) :: string1, string2, string3

contains

!-----------------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!-----------------------------------------------------------------------------
!> @todo Test RMA.
! Given DART state vector and a location, 
! it computes thermospheric neutral density [Kg/m3] 
! The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_upper_atm_density(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val(ens_size)
integer,            intent(out) :: istatus(ens_size)

real(r8) :: mmro1(ens_size), mmro2(ens_size) ! mass mixing ratio 
real(r8) :: mass_reciprocal(ens_size), pressure(ens_size), temperature(ens_size)
integer  :: this_istatus(ens_size)
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0

! Some models (i.e. GITM) have density as part of the state.
! If it is available, just return it. If density is not state,
! then we need to create it from its constituents.

call interpolate(state_handle, ens_size, location, KIND_DENSITY, obs_val, istatus)
if(any(istatus == 0)) return ! density is part of the state


! This part was implemented for TIEGCM. Check the units for use with
! other models.
istatus(:) = 0
call interpolate(state_handle, ens_size, location, KIND_ATOMIC_OXYGEN_MIXING_RATIO, mmro1, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, KIND_MOLEC_OXYGEN_MIXING_RATIO, mmro2, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, KIND_PRESSURE, pressure, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return

call interpolate(state_handle, ens_size, location, KIND_TEMPERATURE, temperature, this_istatus)
call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
if (return_now) return


! density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
! where M is the mean molar mass 
! 1/M = sum(wi/Mi) where wi are mass mixing fractions and Mi are individual molar masses

where (istatus == 0) 
   mass_reciprocal = mmro1/O_molar_mass + mmro2/O2_molar_mass + &
                    (1.0_r8-mmro1-mmro2)/N2_molar_mass

   obs_val = pressure / mass_reciprocal / temperature / universal_gas_constant 
endwhere

end subroutine get_expected_upper_atm_density

!-----------------------------------------------------------------------------

! Given DART state vector and a location, 
! it computes ground GPS vertical total electron content
! The istatus variable should be returned as 0 unless there is a problem
!> @todo Is the logic correct in this code on the Trunk
!>  Should you return from the subroutine instead of exiting 
!> the loop at exit LEVELS
subroutine get_expected_gnd_gps_vtec(state_handle, ens_size, location, obs_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val(ens_size)
integer,            intent(out) :: istatus(ens_size)

! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be in the
! integrated column from an instrument looking straight down at the tangent point.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

integer  :: nAlts, iAlt, this_istatus(ens_size)
real(r8), dimension(ens_size, MAXLEVELS) :: ALT, IDensityS_ie  ! num_ens by num levels
real(r8) :: loc_vals(3)
real(r8) :: tec(ens_size)
type(location_type) :: probe
logical  :: return_now

if ( .not. module_initialized ) call initialize_module

istatus = 0     ! must be 0 to use track_status()
obs_val = MISSING_R8

loc_vals = get_location(location)

nAlts = 0
LEVELS: do iAlt=1, size(ALT)+1
   ! loop over levels.  if we get to one more than the allocated array size,
   ! this model must have more levels than we expected.  increase array sizes,
   ! recompile, and try again.
   if (iAlt > size(ALT)) then
      write(string1,'(''more than '',i4,'' levels in the model.'')') MAXLEVELS
      string2='increase MAXLEVELS in obs_def_upper_atm_mod.f90, rerun preprocess and recompile.'
      string3='increase ALT, IDensityS_ie array sizes in code and recompile'
      call error_handler(E_ERR, 'get_expected_gnd_gps_vtec', string1, &
           source, revision, revdate, text2=string2, text3=string3)
   endif

   ! At each altitude interpolate the 2D IDensityS_ie to the lon-lat where data 
   ! point is located. After this loop we will have a column centered at the data 
   ! point's lon-lat and at all model altitudes.
   probe = set_location(loc_vals(1), loc_vals(2), real(iAlt, r8), VERTISLEVEL) !probe is where we have data 

   call interpolate(state_handle, ens_size, probe, KIND_DENSITY_ION_E, IDensityS_ie(:, iAlt), this_istatus) 
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit LEVELS

   call interpolate(state_handle, ens_size, probe, KIND_GEOPOTENTIAL_HEIGHT, ALT(:, iAlt), this_istatus) 
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit LEVELS
   
   nAlts = nAlts+1
enddo LEVELS

if (nAlts == 0) return

tec=0.0_r8 !start with zero for the summation

do iAlt = 1, nAlts-1 !approximate the integral over the altitude as a sum of trapezoids
   !area of a trapezoid: A = (h2-h1) * (f2+f1)/2
   where (istatus == 0) &
      tec = tec + ( ALT(:, iAlt+1)-ALT(:, iAlt) )  * ( IDensityS_ie(:, iAlt+1)+IDensityS_ie(:, iAlt) ) /2.0_r8
enddo
where (istatus == 0) &
   obs_val = tec * 10.0**(-16) !units of TEC are "10^16" #electron/m^2 instead of just "1" #electron/m^2

! return code set by track_status

end subroutine get_expected_gnd_gps_vtec

!-----------------------------------------------------------------------------

! Given DART state vector and a location, 
! it computes thermospheric neutral density [Kg/m3] 
! The istatus variable should be returned as 0 unless there is a problem

subroutine get_expected_vtec(state_handle, ens_size, location, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),           intent(out) :: expected_obs(ens_size)
integer,            intent(out) :: istatus(ens_size)


if ( .not. module_initialized ) call initialize_module

call error_handler(E_ERR, 'get_expected_vtec', 'routine needs to be written', &
           source, revision, revdate)

expected_obs = missing_r8
istatus = 1

end subroutine get_expected_vtec

!-----------------------------------------------------------------------------

! First, find the number of levels in the model.
! Then, loop down through the levels to create a top-down vertical profile.
!       As we do that, we accumulate the amount of N2 and O, stopping when
!       the N2 reaches 10^21 M^-2. This will probably mean only using part
!       of the 'last' layer.

subroutine get_expected_O_N2_ratio(state_handle, ens_size, location, obs_val, istatus)
 
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val(ens_size)
integer,            intent(out) :: istatus(ens_size)

real(r8) :: loc_array(3)
real(r8) :: loc_lon, loc_lat
type(location_type) :: loc
logical  :: return_now

real(r8), parameter :: Max_N2_column_density = 1.0E21_r8

real(r8) :: N2_total(ens_size)
real(r8) :: O_total(ens_size)

real(r8) :: O_mmr(ens_size, MAXLEVELS)
real(r8) :: O2_mmr(ens_size, MAXLEVELS)
real(r8) :: pressure(ens_size, MAXLEVELS)
real(r8) :: temperature(ens_size, MAXLEVELS)
real(r8) :: heights(ens_size, MAXLEVELS)
real(r8) :: thickness(ens_size, MAXLEVELS)
real(r8) :: O_integrated(ens_size)
real(r8) :: N2_integrated(ens_size)

real(r8), allocatable :: N2_mmr(:, :)
real(r8), allocatable :: mbar(:, :)
real(r8), allocatable :: N2_number_density(:, :)
real(r8), allocatable :: total_number_density(:, :)
real(r8), allocatable :: O_number_density(:, :)

real(r8), PARAMETER :: k_constant = 1.381e-23_r8 ! m^2 * kg / s^2 / K
integer :: ilayer, nlevels, nilevels
integer :: this_istatus(ens_size)
real(r8) :: layerfraction(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus = 0
obs_val = MISSING_R8

call error_handler(E_ERR, 'get_expected_O_N2_ratio', 'routine not tested', &
           source, revision, revdate, &
           text2='routine in obs_def/obs_def_upper_atm_mod.f90', &
           text3='test and inform the DART development team. Thanks -- Tim.')

if ( .not. module_initialized ) call initialize_module

loc_array = get_location(location) ! loc is in DEGREES
loc_lon   = loc_array(1)
loc_lat   = loc_array(2)

! some variables are defined on interface layers

nilevels = 0
heights = 0.0_r8

FILLINTERFACES : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_handle, ens_size, loc, KIND_GEOMETRIC_HEIGHT, heights(:, ilayer), istatus)
   if (any(istatus /= 0)) exit FILLINTERFACES

   nilevels = nilevels + 1

enddo FILLINTERFACES


if (nilevels == 0) return

istatus(:) = 0
thickness = 0.0_r8
thickness(:, 1:nilevels-1) = heights(:, 2:nilevels) - heights(:, 1:nilevels-1)

! Some variables are defined on midpoints of the layers

nlevels = 0

FILLMIDPOINTS : do ilayer = 1,MAXLEVELS

   loc = set_location(loc_lon, loc_lat, real(ilayer,r8), VERTISLEVEL)

   call interpolate(state_handle, ens_size, loc, KIND_PRESSURE, pressure(:, ilayer), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (any(istatus /= 0)) exit FILLMIDPOINTS

   call interpolate(state_handle, ens_size, loc, KIND_TEMPERATURE, temperature(:, ilayer), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, loc, KIND_ATOMIC_OXYGEN_MIXING_RATIO, O_mmr(:, ilayer), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (return_now) return

   call interpolate(state_handle, ens_size, loc, KIND_MOLEC_OXYGEN_MIXING_RATIO, O2_mmr(:, ilayer), this_istatus)
   call track_status(ens_size, this_istatus, obs_val, istatus, return_now)
   if (return_now) return

   nlevels = nlevels + 1

enddo FILLMIDPOINTS

if (nlevels == 0) return

! Check to make sure we have more interfaces than layers.
!> @TODO should this be an error instead of a message?

if (nilevels /= (nlevels+1)) then
   write(string1,*)'Require there to be 1 more interfaces than midpoints.'
   write(string2,*)'Found ',nilevels,' interface layers.'
   write(string3,*)'Found ',nlevels,' midpoint layers.'
   call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
              source, revision, revdate, text2=string2,text3=string3)
   where (istatus == 0) obs_val = missing_r8
   where (istatus == 0) istatus = 11
   return
endif

! calculate what we can using array notation

allocate(N2_mmr(ens_size, nlevels), mbar(ens_size, nlevels), total_number_density(ens_size, nlevels), &
         O_number_density(ens_size, nlevels), N2_number_density(ens_size, nlevels))

N2_mmr = 1.0_r8 -  O_mmr(:, 1:nlevels) - O2_mmr(:, 1:nlevels)
  mbar = 1.0_r8/( O2_mmr(:, 1:nlevels)/O2_molar_mass + &
                   O_mmr(:, 1:nlevels)/ O_molar_mass + &
                  N2_mmr(:, 1:nlevels)/N2_molar_mass )

! O_mmr and N2_mmr defined at midpoints, heights defined at interfaces, so the
! calculated thicknesses apply directly to the O and N2 densities.

total_number_density = pressure(:, 1:nlevels) / (k_constant * temperature(:, 1:nlevels))

 O_number_density =  O_mmr(:, 1:nlevels) * mbar /  O_molar_mass * total_number_density 
N2_number_density = N2_mmr(:, 1:nlevels) * mbar / N2_molar_mass * total_number_density

if ( 1 == 2 ) then ! DEBUG BLOCK NOT IN USE
   write(*,*)
   do ilayer = nlevels,1,-1
      write(*,*)'DEBUG level, thickness, ens member 1: ',ilayer, thickness(1, ilayer), &
            O_number_density(1, ilayer), N2_number_density(1, ilayer), &
                 temperature(1, ilayer), total_number_density(1, ilayer), &
                 O2_mmr(1, ilayer), O_mmr(1, ilayer), N2_mmr(1, ilayer), mbar(1, ilayer)
   enddo
   write(*,*)
endif

N2_total = 0.0_r8
 O_total = 0.0_r8

TOPDOWN : do ilayer = nlevels,1,-1

   if (ilayer == 1) then
      write(string1,*)'Integrated all the way down to the surface.'
      write(string2,*)'Still do not have ',Max_N2_column_density,' nitrogen molecules per m^2'
      call error_handler(E_MSG,'get_expected_O_N2_ratio', string1, &
                 source, revision, revdate, text2=string2)
      where (istatus == 0) istatus = 2
      return
   endif

   ! integrate over layer thickness
   O_integrated  =  O_number_density(:, ilayer) * thickness(:, ilayer)
   N2_integrated = N2_number_density(:, ilayer) * thickness(:, ilayer)

   where ((N2_total+N2_integrated) >= Max_N2_column_density) 
      ! only store part of the final layer so as not to overshoot 10^21 m^-2
      ! Let y2 == N2_total, y = Max_N2_column_density, y1 = N2_total + N2_integrated
      ! the layer fraction is (y - y2)/(y1-y2) 
      ! (Max_N2_column_density - N2_total)/(N2_total + N2_integrated - N2_total)
      layerfraction = (Max_N2_column_density - N2_total) / N2_integrated
      N2_total = N2_total + N2_integrated*layerfraction
       O_total =  O_total +  O_integrated*layerfraction
   elsewhere
      N2_total = N2_total + N2_integrated
       O_total =  O_total +  O_integrated
   endwhere

   if (any((N2_total+N2_integrated) >= Max_N2_column_density)) exit TOPDOWN

enddo TOPDOWN

where (istatus == 0) obs_val = O_total / N2_total

deallocate(N2_mmr, mbar, total_number_density, O_number_density, N2_number_density)

end subroutine get_expected_O_N2_ratio


end module obs_def_upper_atm_mod
! END DART PREPROCESS MODULE CODE      

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
