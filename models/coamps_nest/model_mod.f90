! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

module model_mod

!------------------------------
! MODULE:       model_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!               Manhattan (updated jun 2017)
!
! Module for the interface between DART and the U. S. Navy's COAMPS
! mesoscale model.  COAMPS was developed by the Naval Research Laboratory,
! Monterey, California.  COAMPS is a registered trademark of the Naval
! Research Laboratory.
!------------------------------ 

use coamps_domain_mod,   only : coamps_domain,            &
                                initialize_domain,        &
                                dump_domain_info,         &
                                get_domain_nest,          &
                                get_nest_count,           &
                                get_domain_num_levels,    &
                                get_domain_wsigma,        &
                                get_domain_msigma,        &
                                get_domain_dsigmaw,       &
                                latlon_to_nest_point,     &
                                nest_point_to_latlon


use coamps_statevec_mod, only : state_vector,             &
                                state_iterator,           &
                                get_iterator,             &
                                has_next,                 &
                                get_next,                 &
                                initialize_state_vector,  &
                                find_state_variable,      &
                                get_num_fields,           &
                                get_var_by_index,         &
                                dump_state_vector,        &
                                get_total_size,           &
                                construct_domain_info


use coamps_statevar_mod, only : state_variable,           &
                                get_var_substate,         &
                                get_nest_number,          &
                                get_vert_type,            &
                                get_vert_loc,             &
                                get_state_begin,          &
                                get_state_end,            &
                                is_var_at_index,          &
                                dump_state_variable,      &
                                get_var_kind

use coamps_nest_mod,     only : coamps_nest,              &
                                get_terrain,              &
                                get_terrain_height_at_points, &
                                get_nest_i_width,         &
                                get_nest_j_width,         &
                                get_nest_delta_x,         &
                                get_nest_delta_y,         &
                                get_nest,                 &
                                make_nest_point,          &
                                nest_point,               &
                                get_nest_latlon,          &
                                dump_nest_info,          &
                                nest_index_1d_to_3d

use coamps_intrinsic_mod, only : vor,              &
                                 z2zint

use coamps_interp_mod,   only : interpolate,              &
                                set_interp_diag

use coamps_util_mod,     only : check_alloc_status,       &
                                set_debug_level,          &
                                timestamp_message,        &
                                dump_data_file,           &
                                check_dealloc_status,     &
                                HDF5_FILE_NAME

use coamps_netcdf_mod,   only : nc_write_prognostic_atts, &
                                nc_write_prognostic_data

use coamps_translate_mod, only : initialize_translator,        &
                                 finalize_translator,          &
                                 record_hdf_varnames,          &
       generate_coamps_varnames => generate_coamps_filenames,  &
                                 get_dtg

!#!    use coamps_pert_mod,     only : perturb_state

use location_mod,        only : get_close_type,                &
                                get_dist,                      &
                                get_location,                  &
                                location_type,                 &
                       loc_get_close_state => get_close_state, &
                       loc_get_close_obs   => get_close_obs,   &
                                set_location,                  &
                                vertical_localization_on,      &
                                set_vertical_localization_coord, &
                                query_location,                &
                                VERTISLEVEL,                   &
                                VERTISPRESSURE,                &
                                VERTISHEIGHT,                  &
                                VERTISSURFACE,                 &
                                VERTISUNDEF,                   &
                                VERTISSCALEHEIGHT

use obs_kind_mod

use ensemble_manager_mod, only : ensemble_type

use time_manager_mod,    only : set_time, set_time_missing,    &
                                set_date, get_date, time_type, &
                                set_calendar_type,             &
                                print_date, print_time,        &
                                operator(+), operator(-)

use types_mod,           only : MISSING_R8, MISSING_I, DEG2RAD, &
                                r8, i8

use utilities_mod,       only : check_namelist_read,           &
                                do_output,                     &
                                E_ERR,                         &
                                E_MSG,                         &
                                error_handler,                 &
                                find_namelist_in_file,         &
                                get_unit,                      &
                                register_module,               &
                                string_to_real,                &
                                string_to_logical,             &
                                to_upper,                      &
                                do_nml_file,                   &
                                do_nml_term,                   &
                                nmlfileunit

use default_model_mod,  only :  init_conditions,               &
                                init_time,                     &
                                adv_1step,                     &
                                nc_write_model_vars

use state_structure_mod, only : add_domain,                    &
                                get_domain_size,               &
                                state_structure_info

use netcdf_utilities_mod, only : nc_add_global_attribute,      &
                                 nc_add_global_creation_time,  &
                                 nc_check

use netcdf
use typesizes

implicit none
private

!-------------------------------------------------------
! BEGIN PUBLIC INTERFACE - updated for Manhattan release
!-------------------------------------------------------

! Initialization/finalization
public :: static_init_model 
public :: end_model      

! NetCDF diagnostics
public :: nc_write_model_atts 

! Ensemble generation
public :: pert_model_copies 

! Forward operator
public :: model_interpolate

! Localization
public :: get_close_obs
public :: get_close_state

! Vertical conversion if multiple choices for vert coordinate
public :: convert_vertical_obs
public :: convert_vertical_state

! Information about model setup
public :: get_model_size 
public :: get_state_meta_data 
public :: shortest_time_between_assimilations
public :: get_coamps_domain   ! not required by DART

! Null interfaces - code in other modules
public :: init_conditions
public :: init_time      
public :: adv_1step 
public :: nc_write_model_vars

! Time management - must conform to COAMPS file format for time
public :: read_model_time
public :: write_model_time

!------------------------------
! END PUBLIC INTERFACE
!------------------------------

!------------------------------
! BEGIN EXTERNAL INTERFACE
!------------------------------
  !  [none]
!------------------------------
! END EXTERNAL INTERFACE
!------------------------------

!------------------------------
! BEGIN TYPES AND CONSTANTS
!------------------------------
!  [none]
!------------------------------
! END TYPES AND CONSTANTS
!------------------------------

!------------------------------
! BEGIN MODULE VARIABLES
!------------------------------

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
      
character(len=512) :: string1, string2
logical, save :: module_initialized = .false.

! DART state vector contents are computed from coamps_statevec_mod and state.vars
integer, parameter :: MAX_STATE_VARIABLES = 40

character(len=NF90_MAX_NAME) :: var_names(MAX_STATE_VARIABLES)   = ' '
logical  ::                   update_list(MAX_STATE_VARIABLES)   = .FALSE.
integer  ::                     kind_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::                    clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

! Main model_mod namelist - not too much here as we read most of
! the data we need in from the COAMPS files themselves
character(len=10) :: cdtg                 = '1999083100' ! Date-time group
integer           :: y_bound_skip         = 3            ! How many x and y boundary
integer           :: x_bound_skip         = 3            ! points to skip when
                                                         ! perturbing the model
                                                         ! state
logical           :: need_mean            = .true.       ! Do we need the ensemble
                                                         ! mean for for forward
                                                         ! operator computation?
logical           :: output_interpolation = .false. 
integer           :: debug                = 0            ! increase for debug messages
integer           :: assimilation_period_days = 0
integer           :: assimilation_period_seconds = 216000

namelist /model_nml/ cdtg, y_bound_skip, x_bound_skip, need_mean, &
                     output_interpolation, debug, &
                     assimilation_period_days, assimilation_period_seconds

! Locations of state variables
integer, dimension(:), allocatable :: all_vars

! Grid information structure
type(coamps_domain) :: domain
type(state_vector)  :: state_definition   ! ENTIRE COAMP_NEST STATE
type(state_vector)  :: state_layout_3D    ! Just the variables for DART

! Ensemble mean
real(kind=r8), dimension(:), allocatable :: ensemble_mean

! default to localizing in pressure.  override with namelist
integer :: vertical_localization_type = VERTISPRESSURE

integer :: nfields
integer :: domid

!------------------------------
! END MODULE VARIABLES
!------------------------------

contains

!------------------------------
! BEGIN PUBLIC ROUTINES
!------------------------------

!-------------------------------------------------------------------------------
!> One-time initialization of the model.  For COAMPS, this:
!>  1. Reads in the model_mod namelist
!>  2. Initializes the pressure levels for the state vector
!>  3. Generate the location data for each member of the state
!>  4. Queries the template file 'dart_vector.nc' to glean variable sizes
!>  PARAMETERS
!>   [none]
!>
!> latitu_sfc_000000_000000_1a0201x0204_2013011000_00000000_fcstfld
!> longit_sfc_000000_000000_1a0201x0204_2013011000_00000000_fcstfld

subroutine static_init_model()

character(len=*), parameter :: STATE_VEC_DEF_FILE = 'state.vars'
character(len=*), parameter :: routine = 'static_init_model'

integer :: ivar, nvars
integer :: inest, numnests
integer(i8) :: model_size

if (module_initialized) return ! only need to do this once

call register_module(source, revision, revdate)

module_initialized = .true.

call set_calendar_type('Gregorian')
call read_model_namelist()
call set_debug_level(debug)

! the domain information is reported upon initialization
call initialize_domain(HDF5_FILE_NAME, cdtg, domain)

call set_interp_diag(output_interpolation)

! 'state_definition' contains state vector necessary for entire coamps_nest
call initialize_state_vector(state_definition, STATE_VEC_DEF_FILE, domain)

! 'state_layout_3D' contains state vector necessary for DART
call initialize_state_vector(state_layout_3D, STATE_VEC_DEF_FILE, domain, .true.)
call allocate_metadata_arrays()
call populate_metadata_arrays()
call initialize_translator()
call generate_coamps_varnames()
call record_hdf_varnames(state_layout_3D)

if (debug > 0 .and. do_output()) call dump_state_vector(state_layout_3D)

numnests = get_nest_count(domain)  ! DART 'domain' is a COAMPS 'nest'

model_size = 0_i8
NESTLOOP : do inest = 1,numnests

   call construct_domain_info(state_layout_3D, inest, var_names, kind_list, clamp_vals, update_list, nvars)

   domid = add_domain('dart_vector.nc', nvars, var_names, kind_list, clamp_vals, update_list )

   ! print information in the state structure

   if (debug > 0 .and. do_output()) call state_structure_info(domid)

   model_size = model_size + get_domain_size(domid)

enddo NESTLOOP

if (debug > 0 .and. do_output()) then
write(string1, *)'static_init_model: model_size = ', model_size
call error_handler(E_MSG, routine, string1)
endif
call finalize_translator()

end subroutine static_init_model


!-------------------------------------------------------------------------------
!> Clean up the workspace once the program is ending
!>  PARAMETERS
!>   [none]

subroutine end_model()
    real(kind=r8) x

    x = MISSING_R8

    call deallocate_metadata_arrays()
end subroutine end_model


!-------------------------------------------------------------------------------
!> Write model-specific global attributes to a NetCDF file.
!>  PARAMETERS
!>   IN  ncFileID          numeric ID of an *open* NetCDF file
!>   OUT ierr              0 if writing was successful

subroutine nc_write_model_atts( ncFileID, domain_id ) 
integer, intent(in) :: ncFileID
integer, intent(in) :: domain_id

call nc_check(nf90_redef(ncFileID),'nc_write_model_atts', 'redef')

call nc_write_prognostic_atts(ncFileID, state_layout_3D, define_vars=.false.)

end subroutine nc_write_model_atts


!-------------------------------------------------------------------------------
!> Perturb the model state, field by field.  This can be done 3
!> different ways:
!>  1. No perturbation
!>  2. Uniform perturbation - each element of the field has the
!>                            same additive perturbation
!>  3. Individual perturbation - each element of the field has a 
!>                               different additive perturbation
!> The perturbation magnitude and option are supplied out of the
!> dynamic restart vector definition - this allows us to supply a
!> variance appropriate for each type of variable at each level.
!>  PARAMETERS
!>   IN  state             DART state vector
!>   OUT pert_state        state vector after perturbations
!>   OUT interf_provided   true if this routine did the perturbation

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_copies


!#!    subroutine pert_model_state(state, pert_state, interf_provided)
!#!        real(kind=r8), dimension(:), intent(in)  :: state
!#!        real(kind=r8), dimension(:), intent(out) :: pert_state
!#!        logical,                     intent(out) :: interf_provided
!#!
!#!        call error_handler(E_MSG, 'pert_model_state', 'Perturbing '// &
!#!                           'model state', source, revision, revdate) 
!#!        
!#!        pert_state = perturb_state(state, state_definition, x_bound_skip, y_bound_skip)
!#!        interf_provided = .true.
!#!    end subroutine pert_model_state
    
    ! model_interpolate
    ! -----------------
    ! Given the DART state vector, a location, and a raw variable type
    ! interpolates the state to that location
    ! This implementation currently only supports pressure coordinates.
    !  PARAMETERS
    !   IN  x                 full state vector
    !   IN  location          DART location structure for interpolation
    !                         target location
    !   IN  obs_kind          raw variable kind
    !   OUT obs_val           interpolated value
    !   OUT interp_status     status of interpolation (0 is success)
    !
    ! ORIGINAL:
    !subroutine model_interpolate(x, location, obs_kind, obs_val, interp_status)
    !    real(r8), dimension(:), intent(in)  :: x
    !    type(location_type),    intent(in)  :: location
    !    integer,                intent(in)  :: obs_kind
    !    real(r8),               intent(out) :: obs_val
    !    integer,                intent(out) :: interp_status
    !
    ! NEW WITH MANHATTAN:
    !  PARAMETERS
    !   IN  state_handle      full state vector handle
    !   IN  ens_size          ensemble size
    !   IN  location          DART location structure for interpolation
    !                         target location
    !   IN  obs_kind          raw variable quantity
    !   OUT obs_val           interpolated value (now ens_size array)
    !   OUT interp_status     all 0s if interpolation was successful
    !                         900 unspecified failure
    !                         901 where the location is not in domain or on an unsupported level type
    !                         902 where there are not enough vertical levels
    !                         903 where the location is too high or too low (extrapolation)
    !                         904 where unable to interpolate to a single level
    !                         915 where altimeter is unrealistic
    !                         999 vortex obs ... untested at this point, skipping

    subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, interp_status)
        type(ensemble_type),    intent(in)  :: state_handle
        integer,                intent(in)  :: ens_size
        type(location_type),    intent(in)  :: location
        integer,                intent(in)  :: obs_kind
        real(r8),               intent(out) :: expected_obs(ens_size)
        integer,                intent(out) :: interp_status(ens_size)

!       logical :: interp_worked(ens_size)
        logical :: in_domain

        type(coamps_nest)      :: nest
        type(location_type)    :: cur_loc 
        type(nest_point)       :: nest_pt
        type(state_variable)   :: u_var, v_var

        integer                :: loc_which
        integer                :: i, j, nx, ny, nz
        integer,  dimension(2) :: ij
        real(r8), dimension(3) :: loc_array
        real(r8)               :: delx, dely, ztop

        real(r8), dimension(1) :: terrain

        !real(r8), allocatable  :: heights(:, :)
        real(r8), allocatable  :: heights(:)
        real(r8), allocatable  :: fm(:, :)
        real(r8), allocatable  :: vort_z(:, :, :)
        real(r8), allocatable  :: vort_z_zlev(:, :)
        real(r8), allocatable  :: vort_z_smooth(:, :)
        logical,  allocatable  :: mask(:,:)

        real(r8), parameter    :: vort_srch_radius = 5.0_r8
        real(r8), parameter    :: vrtx_scale = 200000.0_r8
        real(r8), parameter    :: vort_max_lev = 100.0_r8

        integer,  parameter    :: USE_MISSING_VALUE = 1
        integer,  parameter    :: DART_LOC_LON  = 1
        integer,  parameter    :: DART_LOC_LAT  = 2
        integer,  parameter    :: DART_LOC_VERT = 3

        character(len=*), parameter :: routine = 'model_interpolate'
        integer                     :: alloc_status, dealloc_status

        expected_obs(:) = MISSING_R8
        interp_status(:) = 900

        select case (obs_kind)    
        case (QTY_VORTEX_LAT, QTY_VORTEX_LON)

          ! this currently HAS NOT BEEN CONVERTED
          interp_status(:) = 999
          return
         
!#!           ! obs_loc
!#!           loc_array = get_location(location)
!#!           loc_which = nint(query_location(location))
!#! 
!#!           ztop = get_domain_wsigma(domain, 1)
!#! 
!#!           call latlon_to_nest_point(domain, loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON),&
!#!                                     nest_pt, in_domain)
!#!           if(.not. in_domain) then
!#!             obs_val = MISSING_R8
!#!             interp_status = 1
!#!             return
!#!           end if
!#! 
!#!           nest = get_nest(nest_pt)
!#! 
!#!           nx = get_nest_i_width(nest)
!#!           ny = get_nest_j_width(nest)
!#!           nz = get_domain_num_levels(domain)
!#! 
!#!           delx = get_nest_delta_x(nest)
!#!           dely = get_nest_delta_y(nest)
!#! 
!#!           ! allocate arrays
!#!           allocate(vort_z(nx, ny, nz), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: vort_z')
!#! 
!#!           allocate(vort_z_zlev(nx, ny), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: vort_z_zlev')
!#! 
!#!           allocate(vort_z_smooth(nx, ny), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: vort_z_smooth')
!#! 
!#!           !allocate(heights(nx*ny, nz), stat=alloc_status )
!#!           allocate(heights(nz), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: heights')
!#! 
!#! 
!#!           ! (1) Get u and v location in the state vector 
!#!           u_var = find_state_variable(state_layout_3D, nest, &
!#!                    QTY_U_WIND_COMPONENT, .false., 'M', nz)
!#! 
!#!           v_var = find_state_variable(state_layout_3D, nest, &
!#!                    QTY_V_WIND_COMPONENT, .false., 'M', nz)
!#! 
!#!           ! (2) Calculate vorticity
!#!           call timestamp_message('VRTX: Before calculate vorticity')
!#!           allocate(fm(nx, ny), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: fm')
!#! 
!#!           fm(:,:) = 1.0_r8
!#! 
!#!           ! get_var_substate() uses old x() array
!#!           !call vor(reshape(get_var_substate(u_var, x),(/nx, ny, nz/)),                         &
!#!           !         reshape(get_var_substate(v_var, x),(/nx, ny, nz/)),                         &
!#!           !         nx, ny, nz, delx, dely, fm,                              &
!#!           !         get_domain_msigma(domain), get_domain_wsigma(domain),    &
!#!           !         get_domain_dsigmaw(domain), get_terrain(nest), 1, vort_z )
!#! 
!#!           deallocate(fm, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: fm')
!#! 
!#!           call timestamp_message('VRTX: After calculate vorticity')
!#! 
!#!           ! (3) Compute sigma level heights.  Since we want to interpolate to a 
!#!           ! a height above ground level we will not add the surface height.
!#!           !call timestamp_message('VRTX: Before calculate heights')
!#!           !heights = spread(get_domain_msigma(domain), 1, nx*ny) * &
!#!           !          ( 1 - spread(reshape(get_terrain(nest)/ztop, (/nx*ny/)), 2, nz) )
!#!           !call timestamp_message('VRTX: After calculate heights')
!#! 
!#!           ! (4) Interpolate vorticity to the deisred level
!#!           call timestamp_message('VRTX: Before interpolate vorticity')
!#!           do i=1,nx ; do j=1,ny
!#! 
!#!             call get_terrain_height_at_points(nest, (/i/), (/j/), terrain) 
!#! 
!#!             heights(:) = get_domain_msigma(domain) * (1 - terrain(1)/ztop)
!#! 
!#!             call z2zint(vort_z(i, j, :), vort_z_zlev(i, j), heights(:),   &
!#!                         (/vort_max_lev/), nz, 1, 1, USE_MISSING_VALUE, MISSING_R8 )
!#! 
!#!           end do ; end do
!#! 
!#! !          call z2zint(reshape(vort_z, (/nx*ny, nz/)), reshape(vort_z_zlev, (/nx*ny, nz/)), &
!#! !                      reshape(heights,(/nx*ny, nz/)), (/vort_max_lev/), nz, 1, nx*ny,      &
!#! !                      USE_MISSING_VALUE, MISSING_R8)
!#!           call timestamp_message('VRTX: After interpolate vorticity')
!#! 
!#!           ! (5) Smooth vorticity
!#!           call timestamp_message('VRTX: Before smooth vorticity')
!#!           call kernal_smoother(vort_z_zlev, vort_z_smooth, &
!#!                                ceiling(vrtx_scale/delx), nx, ny)
!#! 
!#!           call timestamp_message('VRTX: After smooth vorticity')
!#!          
!#!           ! (6) Search for vorticity maximum within a given radius of obs.
!#!           call timestamp_message('VRTX: Before maximum vorticity search')
!#! 
!#!           allocate(mask(nx,ny), stat=alloc_status )
!#!           call check_alloc_status(alloc_status, routine, source, revision, &
!#!                                   revdate, 'model_interpolate: mask')
!#! 
!#!           do i=1,nx ; do j=1,ny
!#!             call nest_point_to_latlon(domain, make_nest_point(nest, i, j),  &
!#!                                       loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON))
!#! 
!#!             cur_loc = set_location(loc_array(DART_LOC_LON), loc_array(DART_LOC_LAT), &
!#!                                    loc_array(DART_LOC_VERT), loc_which)
!#! 
!#!             mask(i,j) = get_dist(location, cur_loc, no_vert=.true.) &
!#!                                  <= (vort_srch_radius*DEG2RAD)
!#!           end do ; end do
!#!           ij = maxloc(vort_z_smooth, mask)
!#! 
!#!           call timestamp_message('VRTX: After maximum vorticity search')
!#! 
!#!           ! (7) Set desired return value
!#!           call nest_point_to_latlon(domain, make_nest_point(nest, ij(1), ij(2)),  &
!#!                                     loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON))
!#!           if(obs_kind == QTY_VORTEX_LAT) then
!#!             obs_val = loc_array(DART_LOC_LAT)
!#!           else
!#!             obs_val = loc_array(DART_LOC_LON)
!#!           endif
!#!           interp_worked(:) = .true.
!#! 
!#!           ! deallocate arrays
!#!           deallocate(vort_z_smooth, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: vort_z_smooth')
!#! 
!#!           deallocate(heights, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: heights')
!#! 
!#!           deallocate(vort_z, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: vort_z')
!#! 
!#!           deallocate(vort_z_zlev, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: vort_z_zlev')
!#! 
!#!           deallocate(mask, stat=dealloc_status)
!#!           call check_dealloc_status(dealloc_status, routine, source, revision, &
!#!                                     revdate, 'model_interpolate: mask')

        case default 

          ! ORIGINAL:
          !call interpolate(x, domain, state_definition, location, &
          !                 obs_kind, obs_val, interp_worked)
          ! NEW:
          do i = 1, ens_size
             call interpolate(state_handle, ens_size, i, domain, state_definition, &
                       location, obs_kind, expected_obs(i), interp_status(i))
          enddo

        end select

        return
    end subroutine model_interpolate

    ! get_close_state
    ! -------------
    ! Gets the number of close state locations.  Wrapper for location
    ! module's get_close_state subroutine.

    subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                               num_close, close_ind, dist, ens_handle)
    
    ! The specific type of the base observation, plus the generic kinds list
    ! for either the state or obs lists are available if a more sophisticated
    ! distance computation is needed.
    
    type(get_close_type),          intent(in)  :: gc
    type(location_type),           intent(inout)  :: base_loc, locs(:)
    integer,                       intent(in)  :: base_type, loc_qtys(:)
    integer(i8),                   intent(in)  :: loc_indx(:)
    integer,                       intent(out) :: num_close, close_ind(:)
    real(r8),            optional, intent(out) :: dist(:)
    type(ensemble_type), optional, intent(in)  :: ens_handle
    
    call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, &
                             loc_indx, num_close, close_ind, dist, ens_handle)
    
    end subroutine get_close_state


    ! get_close_obs
    ! -------------
    ! Gets the number of close observations.  Wrapper for location
    ! module's get_close_obs subroutine.
    !  PARAMETERS
    !   IN  gc                get_close_type accelerator
    !   IN  base_obs_loc      location of the base observation
    !   IN  obs_loc           location of all the observations
    !   IN  base_obs_kind     raw type of the base observation
    !   IN  obs_kind          raw type of all the observations
    !   OUT num_close         how many observations are close to the
    !                         base observation
    !   OUT close_ind         which of the observations are close to
    !                         the base observation
    !   OUT dist              OPTIONAL distance from the observations
    !                         to the base observation
    !   IN  ens_handle        handle to ensemble mean
    subroutine get_close_obs(gc, base_obs_loc, base_obs_type, obs_locs, loc_qtys, loc_types, &
                             num_close, close_ind, dist, ens_mean_handle)
    
    type(get_close_type),          intent(in)  :: gc
    type(location_type),           intent(inout) :: base_obs_loc, obs_locs(:)
    integer,                       intent(in)  :: base_obs_type, loc_qtys(:), loc_types(:)
    integer,                       intent(out) :: num_close
    integer,                       intent(out) :: close_ind(:)
    real(r8),            optional, intent(out) :: dist(:)
    type(ensemble_type), optional, intent(in)  :: ens_mean_handle

    integer                                 :: t_ind, istatus1(1), istatus2(1), k
    integer                                 :: base_which, local_obs_which
    real(r8), dimension(3)                  :: base_array, local_obs_array
    real(r8)                                :: obs_val(1)
    type(location_type)                     :: local_obs_loc

    ! Initialize variables to missing status
    num_close = 0 ; close_ind = -99 ; dist = 1.0e9
    istatus1  = 0 ; istatus2  = 0

    base_which = nint(query_location(base_obs_loc)) 
    base_array = get_location(base_obs_loc)

    ! Consider horizontal localization first.  Also consider undefined vertical coordinate.
    if(.not. vertical_localization_on() .or. base_which == VERTISUNDEF) then
        call loc_get_close_obs(gc, base_obs_loc, base_obs_type, obs_locs, loc_qtys, &
                               loc_types, num_close, close_ind, dist, ens_mean_handle)
    else
      ! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

      if (base_which /= VERTISLEVEL) then
        if(base_which == VERTISSURFACE) then
          base_array(3)=0.0_r8 ; base_which=VERTISLEVEL
        else
          call model_interpolate(ens_mean_handle, 1, base_obs_loc, QTY_VERTLEVEL, obs_val, istatus1)
          base_array(3)=obs_val(1) ; base_which=VERTISLEVEL
        endif
        base_obs_loc = set_location(base_array(1),  base_array(2), base_array(3), base_which)
      elseif (base_array(3) == MISSING_R8) then
        istatus1 = 1
      endif

      if (istatus1(1) == 0) then
      ! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
      ! This way, we are decreasing the number of distance computations that will follow.
      ! This is a horizontal-distance operation and we don't need to have the relevant vertical
      ! coordinate information yet (for obs_loc).
        call loc_get_close_obs(gc, base_obs_loc, base_obs_type, obs_locs, loc_qtys, &
                               loc_types, num_close, close_ind)

        ! Loop over potentially close subset of obs priors or state variables
        do k = 1, num_close
          t_ind           = close_ind(k)
          local_obs_loc   = obs_locs(t_ind)
          local_obs_which = nint(query_location(local_obs_loc))
          local_obs_array = get_location(local_obs_loc)

          ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
          ! This should only be necessary for obs priors, as state location information already
          ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
          if (local_obs_which /= VERTISLEVEL) then
            if(local_obs_which == VERTISSURFACE) then
              local_obs_array(3)=obs_val(1) ; local_obs_which=VERTISLEVEL
            elseif(local_obs_which == VERTISUNDEF) then
              local_obs_which=VERTISUNDEF
            else
              call model_interpolate(ens_mean_handle, 1, obs_locs(t_ind), QTY_VERTLEVEL, obs_val, istatus2)
              local_obs_array(3)=obs_val(1) ; local_obs_which=VERTISLEVEL
            end if

            ! Store the "new" location into the original full local array
            local_obs_loc = set_location(local_obs_array(1),local_obs_array(2), &
                                         local_obs_array(3),local_obs_which)
            obs_locs(t_ind) = local_obs_loc
          endif

          ! Compute distance - set distance to a very large value if vert coordinate is missing
          ! or vert_interpolate returned error (istatus2=1)
          if ((local_obs_array(3) == missing_r8) .or. (istatus2(1) == 1)) then
            dist(k) = 1.0e9
          else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_type, loc_qtys(t_ind))
          endif

        end do
      end if
    end if

    end subroutine get_close_obs

    ! get_model_size
    ! --------------
    ! Returns the size of the DART state vector
    !  PARAMETERS
    !   OUT get_model_size    length of the DART state vector
    function get_model_size()
        integer(i8) :: get_model_size

        get_model_size = int(get_total_size(state_definition), i8)
    end function get_model_size

    ! get_state_meta_data
    ! -------------------
    ! Get the location and variable type for a particular index in
    ! the state vector
    !  PARAMETERS
    !   IN  index_in          position in state vector to query
    !   OUT location          DART location_type for that index
    !   OUT var_type          OPTIONAL numeric variable type

    !>@todo FIXME at present, the grid is being read in from the hdf5 file,
    !> but a lot of the interpolation routines etc depend on grid parameters
    !> from a the datahd structure ... To be consistent, should revert to 
    !> earlier form where the grid is computed from the datahd info. 
    !> The interplation takes care of the stagger but the locations are 
    !> ONLY for the mass points, so the localization distances are off by
    !> some amount.

    subroutine get_state_meta_data(index_in, location, var_type)
        integer(i8),                   intent(in)  :: index_in
        type(location_type), optional, intent(out) :: location
        integer,             optional, intent(out) :: var_type

        type(state_variable)  :: cur_var
        type(coamps_nest)     :: cur_nest
        integer, dimension(3) :: ijk
        integer               :: index_loc
        real(kind=r8)         :: lat, lon

        cur_var = get_var_by_index(state_layout_3D, all_vars(index_in))

        if (present(location)) then 
            cur_nest = get_domain_nest(domain, get_nest_number(cur_var))

            index_loc = index_in - get_state_begin(cur_var) + 1
            ijk = nest_index_1d_to_3d(cur_nest, index_loc)

            ! note there is nothing here that could account for stagger.
            call get_nest_latlon(cur_nest, ijk(1), ijk(2), lat, lon)        

            location = set_location(lon, lat, &
                                    get_vert_loc(cur_var, domain, ijk(3)), &
                                    get_vert_type(cur_var))
        end if

        if (present(var_type)) var_type = get_var_kind(cur_var)

    end subroutine get_state_meta_data


!--------------------------------------------------------------------
! Seems like the strategy is to localize using sigma levels. If that is
! mandatory, then convert_vertical_obs must convert pressure,height, etc to
! sigma levels.

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: status(:)

status(:) = 0

call error_handler(E_ERR,'convert_vertical_obs','routine not written')

end subroutine convert_vertical_obs


!--------------------------------------------------------------------
! Seems like the strategy is to localize using sigma levels. If that is
! mandatory, then convert_vertical_state does nothing since the coamps state is
! already on sigma levels.

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

istatus = 0

! nothing to do

end subroutine convert_vertical_state


!--------------------------------------------------------------------
!> Returns the smallest increment in time that the model is capable
!> of advancing the state.
!>  PARAMETERS
!>   OUT shortest_time_between_assimilations  model timestep as a DART time_type

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations !> model forecast length

shortest_time_between_assimilations = &
     set_time(assimilation_period_seconds,assimilation_period_days)

end function shortest_time_between_assimilations



    ! get_coamps_domain
    ! ---------
    ! Returns the module defined coamps_domain
    function get_coamps_domain()
      type(coamps_domain) :: get_coamps_domain
      get_coamps_domain = domain
    end function get_coamps_domain


!--------------------------------------------------------------------
!>
!> read the time from the cdtg character string in the model_nml namelist.
!> and add in the 'latest' forecast 'tau' (the forecast length is 'tau' - in hours)
!> This is a required interface that has a mandatory argument.
!> In this case, the mandatory argument is not needed.

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

! local variables
type(time_type) :: base_time
integer, allocatable :: seconds(:)
integer :: year, month, day, hour
integer :: forecast_hours

if ( .not. module_initialized ) call static_init_model()

! cdtg has the format '2013011000'

read(cdtg,'(i4,i2,i2,i2)')year,month,day,hour
base_time = set_date(year, month, day, hour)

!>@todo Get the forecast tau from someplace
forecast_hours = 0
read_model_time = base_time + set_time(forecast_hours*3600, days=0)

if (debug > 0 .and. do_output()) then
   call print_date(base_time, 'read_model_time:simulation starting date')
   call print_time(base_time, 'read_model_time:simulation starting time')
   call print_date(read_model_time, 'read_model_time:current model date')
   call print_time(read_model_time, 'read_model_time:current model time')
endif

end function read_model_time


!-----------------------------------------------------------------------
!>
!> write model time to netCDF file when creating files from scratch

subroutine write_model_time(ncid, dart_time)
integer,         intent(in) :: ncid
type(time_type), intent(in) :: dart_time

integer :: iyear, imonth, iday, ihour, iminute, isecond

call nc_check(nf90_redef(ncid), "write_model_time", "redef")

call nc_add_global_attribute(ncid, 'date-time-group', cdtg, &
                       context='model_mod:write_model_time:dtg')

call get_date(dart_time, iyear, imonth, iday, ihour, iminute, isecond)

write(string1,'(I4,"-",I2.2,"-",I2.2,1x,2(I2.2,":"),I2.2)') &
                         iyear, imonth, iday, ihour, iminute, isecond

call nc_add_global_attribute(ncid, 'model_valid_time', trim(string1), &
                       context='model_mod:write_model_time')

call nc_check(nf90_enddef(ncid),"write_model_time", "Enddef" )

end subroutine write_model_time

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! read_model_namelist
    ! -------------------
    ! Read in parameters from the model_nml namelist
    subroutine read_model_namelist

        character(len=*), parameter :: NAMELIST_FILE  = 'input.nml'
        character(len=*), parameter :: MODEL_NAMELIST = 'model_nml'

        integer :: io_status
        integer :: nml_unit

        call find_namelist_in_file(NAMELIST_FILE, MODEL_NAMELIST, nml_unit)
        read (nml_unit,nml=model_nml,iostat=io_status)
        call check_namelist_read(nml_unit, io_status, MODEL_NAMELIST)

        ! Record the namelist values used for the run
        if (do_nml_file()) write(nmlfileunit, nml=model_nml)
        if (do_nml_term()) write(     *     , nml=model_nml)

    end subroutine read_model_namelist

    ! allocate_metadata_arrays
    ! ------------------------
    ! Initialize storage for state vector locations, types, and possibly the
    ! ensemble mean
    subroutine allocate_metadata_arrays

        character(len=*), parameter :: routine = 'allocate_metadata_arrays'
        integer                     :: alloc_status

        allocate( all_vars(get_model_size()), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'all_vars')

        if (need_mean) then
            allocate( ensemble_mean(get_model_size()), stat=alloc_status )
            call check_alloc_status(alloc_status, routine, source, revision, &
            revdate, 'ensemble_mean')
        end if

    end subroutine allocate_metadata_arrays

    ! deallocate_metadata_arrays
    ! --------------------------
    ! Finalize storage for state vector lcoations, types and possibly the
    ! ensemble mean
    subroutine deallocate_metadata_arrays

        character(len=*), parameter :: routine = 'deallocate_metadata_arrays'
        integer                     :: dealloc_status

        deallocate(all_vars, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                  revdate, 'all_vars')

        if (need_mean) then
            deallocate(ensemble_mean, stat=dealloc_status)
            call check_dealloc_status(dealloc_status, routine, source,       &
                                      revision, revdate, 'all_locs')
        end if
    end subroutine deallocate_metadata_arrays

    ! populate_metadata_arrays
    ! ------------------------
    ! Populates the various arrays of metadata for the model state
    subroutine populate_metadata_arrays()

        type(state_variable), pointer :: cur_var
        type(state_iterator)          :: iterator
        integer :: state_ii

        state_ii = 0
        iterator  = get_iterator(state_layout_3D)
        do while (has_next(iterator))
            state_ii = state_ii + 1
            cur_var   =>  get_next(iterator)
            all_vars(get_state_begin(cur_var):get_state_end(cur_var)) = state_ii
        end do

    end subroutine populate_metadata_arrays

    ! index_to_var
    ! ------------------------
    ! returns the state_variable at a given index point
    function index_to_var(index_in)
      integer,            intent(in) :: index_in
      type(state_variable)           :: index_to_var

      type(state_iterator) :: var_iterator
      integer              :: ii

      var_iterator  = get_iterator(state_layout_3D)
      do while (has_next(var_iterator))
         index_to_var = get_next(var_iterator)
         if(is_var_at_index(index_to_var, index_in)) exit 
      end do

    end function index_to_var

    subroutine kernal_smoother(fin, fout, half_width, nx, ny)
      real(kind=r8), intent(in)  :: fin(:, :)
      real(kind=r8), intent(out) :: fout(:, :)
      integer,       intent(in)  :: half_width
      integer,       intent(in)  :: nx
      integer,       intent(in)  :: ny

      real(kind=r8), allocatable :: kernal(:,:)
      real(kind=r8) :: LL, kernal_sum
      integer :: i, j, ii, jj, n_width

        character(len=*), parameter :: routine = 'kernal_smoother'
        integer                     :: alloc_status, dealloc_status

      n_width = 2*half_width + 1

      allocate(kernal(n_width, n_width), stat=alloc_status )
      call check_alloc_status(alloc_status, routine, source, revision, revdate, 'kernal')

      ! Setup kernal
      kernal_sum = 0.0_r8 
      do i=1,n_width ; do j=1,n_width
        LL = sqrt(real((half_width + 1 - j)**2 + (half_width + 1 - i)**2, kind=r8))  
        if(LL <= real(half_width,kind=r8)) then
          kernal(i, j) = 1.0_r8 - LL/real(half_width,kind=r8)
        else
          kernal(i, j) = 0.0_r8
        end if
        kernal_sum = kernal(i, j) + kernal_sum
      end do ; end do

      kernal(:,:) = kernal(:,:) / kernal_sum

      fout(:, :) = 0.0_r8
      do i=1+half_width,nx-half_width ; do j=1+half_width,ny-half_width
        do ii=-half_width,half_width ; do jj=-half_width,half_width
          fout(i, j) = fout(i, j) + fin(i+ii, j+jj)*kernal(ii+half_width+1, jj+half_width+1)
        end do ; end do
      end do ; end do

      deallocate(kernal, stat=dealloc_status )
      call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'kernal')

    end subroutine kernal_smoother

    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module model_mod

