! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! An example of a simple forward operator that involves more than
! just interpolating directly from a state vector in a model.
!
! This section defines a specific type in the left column and
! can be any string you want to use for an observation.  The
! right column must be a generic kind that already exists in
! the obs_kind/DEFAULT_obs_kind_mod.F90 file.

! BEGIN DART PREPROCESS KIND LIST
! AIRNOW_PM25,         KIND_PM25
! END DART PREPROCESS KIND LIST

! This section will be added to the main obs_def_mod.f90 that
! is going to be generated, to allow it to call the code we
! are defining here.

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_AIRNOW_PM25_mod, only : get_expected_airnow_pm25
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! This section will be dropped into a large case statement in the
! main obs_def_mod.f90 code to control what happens with each
! observation type that is processed.

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(AIRNOW_PM25)
!        call get_expected_airnow_pm25(state, location, obs_val, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! The next few sections do nothing because there is no additional
! data to read, write, or prompt for.  But there still needs to be a
! case statement in the large select, so they must be here.

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(AIRNOW_PM25)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(AIRNOW_PM25)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(AIRNOW_PM25)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! This is the code that implements the forward operator.
! Define a module, and make public anything that will be called
! from the main obs_def_mod.f90 file.  Here it is just the
! get_expected routine.  There isn't any initialization needed
! but the stub is there; it could read a namelist if there are
! any run-time options to be set.

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_AIRNOW_PM25_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type
use  assim_model_mod, only : interpolate
use     obs_kind_mod

implicit none
private

public :: get_expected_airnow_pm25

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_AIRNOW_PM25_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.
logical       :: use_log_pm25
namelist /obs_def_AIRNOW_PM25_nml/ use_log_pm25

contains

! ---------------------------------------------------

subroutine initialize_module
! Handle any module initialization tasks
integer ::     iunit, rc

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.
use_log_pm25=.false.
! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_AIRNOW_PM25_nml", iunit)
read(iunit, nml = obs_def_AIRNOW_PM25_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_AIRNOW_PM25_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_AIRNOW_PM25_nml)
if (do_nml_term()) write(     *     , nml=obs_def_AIRNOW_PM25_nml)

end subroutine initialize_module

! ---------------------------------------------------

subroutine get_expected_airnow_pm25(state_vector, location, pm25, istatus)  
 real(r8),            intent(in)  :: state_vector(:)
 type(location_type), intent(in)  :: location
 integer,             intent(out) :: istatus
 real(r8),            intent(out) :: pm25 ! (ug/m^3 dry air)

! Forward operator for PM25.  The argument list to this routine
! must match the call in the GET_EXPECTED_OBS_FROM_DEF section above.

real(r8) :: qvapor   ! water vapor mixing ratio
real(r8) :: p25      ! p25 
real(r8) :: sulf     ! Sulphate 
real(r8) :: BC1      ! Hydrophobic Black Carbon 
real(r8) :: BC2      ! Hydrophilic Black Carbon 
real(r8) :: OC1      ! Hydrophobic Organic Carbon
real(r8) :: OC2      ! Hydrophilic Organic Carbon 
real(r8) :: DUST1    ! Dust 1
real(r8) :: DUST2    ! Dust 2
real(r8) :: DUST3    ! Dust 3
real(r8) :: DUST4    ! Dust 4
real(r8) :: DUST5    ! Dust 5
real(r8) :: SS1      ! Sea Salt 1
real(r8) :: SS2      ! Sea Salt 2
real(r8) :: SS3      ! Sea Salt 3
real(r8) :: SS4      ! Sea Salt 4
real(r8) :: Theta    ! perturbation potential temperature
real(r8) :: T        ! temperature
real(r8) :: Tv       ! virtual temperature
real(r8) :: P        ! perturbation pressure
real(r8) :: rho_d    ! density dry air
real(r8) :: rho_m    ! density moist air
real(r8) :: Rd       ! gas constant dry air
real(r8) :: Cp       ! heat capacity dry air
real(r8) :: Pa_to_torr  ! convert Pa to torr
!
 namelist /obs_def_AIRNOW_PM25_nml/ use_log_pm25

if ( .not. module_initialized ) call initialize_module

! water vapor mixing ratio (kg/kg) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_VAPOR_MIXING_RATIO, qvapor, istatus)
if (istatus /= 0) then
   qvapor = missing_r8
   pm25 = missing_r8
   return
endif
!
! sulfate (ppmv) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_SO4, sulf, istatus)
if (istatus /= 0) then
   sulf = missing_r8
   pm25 = missing_r8
   return
endif
!
! dust (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_DST01, DUST1, istatus)
if (istatus /= 0) then
   DUST1 = missing_r8
   pm25 = missing_r8
   return
endif
!
! dust (ug/kg - dray air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_DST02, DUST2, istatus)
if (istatus /= 0) then
   DUST2 = missing_r8
   pm25 = missing_r8
   return
endif
!
! dust (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_DST03, DUST3, istatus)
if (istatus /= 0) then
   DUST3 = missing_r8
   pm25 = missing_r8
   return
endif
!
! dust (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_DST04, DUST4, istatus)
if (istatus /= 0) then
   DUST4 = missing_r8
   pm25 = missing_r8
   return
endif
!
! dust (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_DST05, DUST5, istatus)
if (istatus /= 0) then
   DUST5 = missing_r8
   pm25 = missing_r8
   return
endif
!
! hydrophilic black carbon (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_BC1, BC1, istatus)
if (istatus /= 0) then
   BC1 = missing_r8
   pm25 = missing_r8
   return
endif
!
! hydrophobic black carbon (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_BC2, BC2, istatus)
if (istatus /= 0) then
   BC2 = missing_r8
   pm25 = missing_r8
   return
endif
!
! hydrophilic organic carbon (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_OC1, OC1, istatus)
if (istatus /= 0) then
   OC1 = missing_r8
   pm25 = missing_r8
   return
endif
!
! hydrophobic organic carbon (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_OC2, OC2, istatus)
if (istatus /= 0) then
   OC2 = missing_r8
   pm25 = missing_r8
   return
endif
!
! sea salt (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_SSLT01, SS1, istatus)
if (istatus /= 0) then
   SS1 = missing_r8
   pm25 = missing_r8
   return
endif
!
! sea salt (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_SSLT02, SS2, istatus)
if (istatus /= 0) then
   SS2 = missing_r8
   pm25 = missing_r8
   return
endif
!
! sea salt (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_SSLT03, SS3, istatus)
if (istatus /= 0) then
   SS3 = missing_r8
   pm25 = missing_r8
   return
endif
!
! sea salt (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_SSLT04, SS4, istatus)
if (istatus /= 0) then
   SS4 = missing_r8
   pm25 = missing_r8
   return
endif
!
! perturbation theta (K) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_POTENTIAL_TEMPERATURE, Theta, istatus)
if (istatus /= 0) then
   Theta = missing_r8
   pm25 = missing_r8
   return
endif
!
! perturbation pressure (Pa) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_PRESSURE, P, istatus)
if (istatus /= 0) then
   P = missing_r8
   pm25 = missing_r8
   return
endif
!
! p25 (ug/kg-dry air) at this location - this calls the model_mod code.
call interpolate(state_vector, location, KIND_P25, p25, istatus)
if (istatus /= 0) then
   p25 = missing_r8
   pm25 = missing_r8
   return
endif
!
! The actual forward operator computation.  This is the value that
! will be returned.  istatus (the return code) of 0 is good,
! return any value > 0 for error.  (values < 0 reserved for
! system use.)
Rd = 287.05
Cp = 1006.0
T = Theta / (100000./P)**(Rd/Cp)
Tv = Theta / (100000./P)**(Rd/Cp) * qvapor/(1.+qvapor)
rho_d = P/(Rd*T)
rho_m = P/(Rd*Tv)
Pa_to_torr = 133.322
!
! Expected pm25 (ug/kg)
if (use_log_pm25) then
   pm25 = (exp(p25) + 1.375*exp(sulf) + exp(BC1) + exp(BC2) + 1.8*(exp(OC1)+exp(OC2)) + &
        exp(DUST1) + .286*exp(DUST2) + exp(SS1) + .942*exp(SS2))
else
   pm25 = (p25 + 1.375*sulf + BC1 + BC2 + 1.8*(OC1+OC2) + &
        DUST1 + .286*DUST2 + SS1 + .942*SS2)
endif
!
! Expected pm25 (ug/m^3) (LC) Should be this
pm25 = pm25*rho_d
!
! Expected pm25 (ug/m^3) (STP)
!pm25 = pm25 * T/298.15 *760./P*pa_to_torr

istatus = 0
    
end subroutine get_expected_airnow_pm25

! ---------------------------------------------------

end module obs_def_AIRNOW_PM25_mod
! END DART PREPROCESS MODULE CODE

