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
! MODIS_AOD_RETRIEVAL,         KIND_AOD
! END DART PREPROCESS KIND LIST

! This section will be added to the main obs_def_mod.f90 that
! is going to be generated, to allow it to call the code we
! are defining here.

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_MODIS_AOD_mod, only : get_expected_modis_aod
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! This section will be dropped into a large case statement in the
! main obs_def_mod.f90 code to control what happens with each
! observation type that is processed.

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(MODIS_AOD_RETRIEVAL)
!        call get_expected_modis_aod(state, location, obs_val, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! The next few sections do nothing because there is no additional
! data to read, write, or prompt for.  But there still needs to be a
! case statement in the large select, so they must be here.

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(MODIS_AOD_RETRIEVAL)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(MODIS_AOD_RETRIEVAL)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(MODIS_AOD_RETRIEVAL)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! This is the code that implements the forward operator.
! Define a module, and make public anything that will be called
! from the main obs_def_mod.f90 file.  Here it is just the
! get_expected routine.  There isn't any initialization needed
! but the stub is there; it could read a namelist if there are
! any run-time options to be set.

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_MODIS_AOD_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE, VERTISUNDEF
use  assim_model_mod, only : interpolate
use     obs_kind_mod
implicit none
private

public :: get_expected_modis_aod

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'obs_def_MODIS_AOD_mod.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

logical, save :: module_initialized = .false.
logical       :: use_log_aod
namelist /obs_def_MODIS_AOD_nml/ use_log_aod

contains

! ---------------------------------------------------

subroutine initialize_module
! Handle any module initialization tasks
integer ::     iunit, rc

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.
use_log_aod=.false.
! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_MODIS_AOD_nml", iunit)
read(iunit, nml = obs_def_MODIS_AOD_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_MODIS_AOD_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_MODIS_AOD_nml)
if (do_nml_term()) write(     *     , nml=obs_def_MODIS_AOD_nml)

end subroutine initialize_module

! ---------------------------------------------------

subroutine get_expected_modis_aod(state_vector, location, modis_aod, istatus)  
 real(r8),            intent(in)  :: state_vector(:)
 type(location_type), intent(in)  :: location
 real(r8),            intent(out) :: modis_aod ! (dimensionless)
 integer,             intent(out) :: istatus

 integer, parameter :: wrf_nlev=33
 integer, parameter :: aod_nlev=wrf_nlev-1

! Forward operator for MODIS AOD.  The argument list to this routine
! must match the call in the GET_EXPECTED_OBS_FROM_DEF section above.

 type(location_type)  :: nloc,nploc
 integer  :: ilev
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
 real(r8) :: PH_dn    ! perturbation geopotential
 real(r8) :: PH_up    ! perturbation geopotential
 real(r8) :: rho_d    ! density dry air
 real(r8) :: rho_m    ! density moist air
 real(r8) :: Rd       ! gas constant dry air
 real(r8) :: Cp       ! heat capacity dry air
 real(r8) :: Pa_to_torr  ! convert Pa to torr
 real(r8) :: grav     ! gravity
 real(r8) :: fac      ! units conversion factor
 real(r8) :: mloc(3)
!
 namelist /obs_def_MODIS_AOD_nml/ use_log_aod

 if ( .not. module_initialized ) call initialize_module

 fac=1.e-6
 grav = 9.8      ! m/s^2
 Rd = 287.05     ! J/kg
 Cp = 1006.0     ! J/kg/K
 Pa_to_torr = 133.322
 modis_aod = 0.  ! dimensionless

 mloc = get_location(location)
 if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
 elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
 endif

 do ilev=1,aod_nlev
    mloc(3)=real(ilev)
    nloc = set_location(mloc(1), mloc(2), mloc(3), VERTISLEVEL)
    mloc(3)=real(ilev+1)
    nploc = set_location(mloc(1), mloc(2), mloc(3), VERTISLEVEL)
!
! water vapor mixing ratio (kg/kg) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_VAPOR_MIXING_RATIO, qvapor, istatus)
    if (istatus /= 0) then
       qvapor = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, qvapor ',qvapor
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_DST01, DUST1, istatus)
    if (istatus /= 0) then
       DUST1 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, DUST1 ',DUST1
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_DST02, DUST2, istatus)
    if (istatus /= 0) then
       DUST2 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, DUST2 ',DUST2
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_DST03, DUST3, istatus)
    if (istatus /= 0) then
       DUST3 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, DUST3 ',DUST3
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_DST04, DUST4, istatus)
    if (istatus /= 0) then
       DUST4 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, DUST4 ',DUST4
!
! dust (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_DST05, DUST5, istatus)
    if (istatus /= 0) then
       DUST5 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, DUST5 ',DUST5
!
! hydrophilic black carbon (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_BC1, BC1, istatus)
    if (istatus /= 0) then
       BC1 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, BC1 ',BC1
!
! hydrophobic black carbon (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_BC2, BC2, istatus)
    if (istatus /= 0) then
       BC2 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, BC2 ',BC2
!
! hydrophilic organic carbon (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_OC1, OC1, istatus)
    if (istatus /= 0) then
       OC1 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, OC1 ',OC1
!
! hydrophobic organic carbon (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_OC2, OC2, istatus)
    if (istatus /= 0) then
       OC2 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, OC2 ',OC2
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_SSLT01, SS1, istatus)
    if (istatus /= 0) then
       SS1 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, SS1 ',SS1
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_SSLT02, SS2, istatus)
    if (istatus /= 0) then
       SS2 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, SS2 ',SS2
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_SSLT03, SS3, istatus)
    if (istatus /= 0) then
       SS3 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, SS3 ',SS3
!
! sea salt (ug/kg - dry air) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_SSLT04, SS4, istatus)
    if (istatus /= 0) then
       SS4 = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, SS4 ',SS4
!
! theta (K) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_POTENTIAL_TEMPERATURE, Theta, istatus)
    if (istatus /= 0) then
       Theta = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, Theta ',Theta
!
! pressure (Pa) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_PRESSURE, P, istatus)
    if (istatus /= 0) then
       P = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, P ',P
!
! sulfate (kg/kg) at this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_SO4, sulf, istatus)
    if (istatus /= 0) then
       sulf = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, sulf ',sulf
!
! geopotential height (m) below this location - this calls the model_mod code.
    call interpolate(state_vector, nloc, KIND_GEOPOTENTIAL_HEIGHT, PH_dn, istatus)
    if (istatus /= 0) then
       PH_dn = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, PH_dn ',PH_dn
!
! geopotential height (m) above this location - this calls the model_mod code.
    call interpolate(state_vector, nploc, KIND_GEOPOTENTIAL_HEIGHT, PH_up, istatus)
    if (istatus /= 0) then
       PH_up = missing_r8
       modis_aod = missing_r8
       return
    endif
!    print *, 'ilev, PH_up ',PH_up
!
! p25 at this location - this calls the model_mod code.
!    call interpolate(state_vector, location, KIND_P25, p25, istatus)
!    if (istatus /= 0) then
!       p25 = missing_r8
!       modis_aod = missing_r8
!       return
!    endif
!    print *, 'ilev, p25 ',p25
!
! The actual forward operator computation.  This is the value that
! will be returned.  istatus (the return code) of 0 is good,
! return any value > 0 for error.  (values < 0 reserved for
! system use.)
!
    T = Theta / (100000./(P**(Rd/Cp)))
    Tv = Theta / (100000./(P**(Rd/Cp))) &
       * qvapor/(1.+qvapor)
    rho_d = P/(Rd*T)
    rho_m = P/(Rd*Tv)
!
! Expected Modis AOD
    if (use_log_aod) then
       modis_aod = modis_aod + (3.673*exp(sulf) + 2.473*(exp(OC1)+exp(OC2)) + 8.99*(exp(BC1)+exp(BC2)) + &
                 2.548*exp(SS1) + .889*exp(SS2) + .227*exp(SS3) + .096*exp(SS4) + 1.596*exp(DUST1) + &
                 .507*exp(DUST2) + .275*exp(DUST3) + .140*exp(DUST4) + .078*exp(DUST5)) * &
                 rho_d * (PH_up-PH_dn) * fac
    else
       modis_aod = modis_aod + (3.673*sulf + 2.473*(OC1+OC2) + 8.99*(BC1+BC2) + &
                 2.548*SS1 + .889*SS2 + .227*SS3 + .096*SS4 + 1.596*DUST1 + &
                 .507*DUST2 + .275*DUST3 + .140*DUST4 + .078*DUST5) * &
                 rho_d * (PH_up-PH_dn) * fac
    endif 
 enddo

 istatus = 0
    
end subroutine get_expected_modis_aod

! ---------------------------------------------------

end module obs_def_MODIS_AOD_mod
! END DART PREPROCESS MODULE CODE

