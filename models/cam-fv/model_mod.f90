! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
! $Id$
!----------------------------------------------------------------
!>
!> this is the interface between the cam-fv atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.
!>
!----------------------------------------------------------------

module model_mod

!>@todo fixme fill in the actual names we use after we've gotten
!>further into writing the coded

use             types_mod, only : MISSING_R8, MISSING_I, i8, r8, vtablenamelength, &
                                  gravity, DEG2RAD
use      time_manager_mod, only : set_time, time_type, set_date, &
                                  set_calendar_type, get_date
use          location_mod, only : location_type, set_vertical, set_location, &
                                  get_location,get_close_obs, get_close_state, & 
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, &
                                  VERTISSCALEHEIGHT, query_location
use         utilities_mod, only : find_namelist_in_file, check_namelist_read, &
                                  string_to_logical, string_to_real,& 
                                  logfileunit, do_nml_file, do_nml_term, &
                                  register_module, error_handler, &
                                  file_exist, to_upper, E_ERR, E_MSG
use          obs_kind_mod, only : QTY_SURFACE_ELEVATION, QTY_PRESSURE, &
                                  QTY_GEOMETRIC_HEIGHT, QTY_VERTLEVEL, &
                                  QTY_SURFACE_PRESSURE, QTY_PRECIPITABLE_WATER, &
                                  QTY_TEMPERATURE, QTY_SPECIFIC_HUMIDITY, &
                                  QTY_GEOPOTENTIAL_HEIGHT,  &
                                  get_index_for_quantity, get_num_quantities
use     mpi_utilities_mod,  only : my_task_id
use        random_seq_mod,  only : random_seq_type, init_random_seq, random_gaussian
use  ensemble_manager_mod,  only : ensemble_type, get_my_num_vars, get_my_vars
use distributed_state_mod,  only : get_state
use   state_structure_mod,  only : add_domain, get_dart_vector_index, get_domain_size, &
                                   get_dim_name, get_kind_index, get_num_dims, &
                                   get_num_variables, get_varid_from_kind, &
                                   get_model_variable_indices, state_structure_info
use  netcdf_utilities_mod,  only : nc_get_variable, nc_get_variable_size, &
                                   nc_add_attribute_to_variable, &
                                   nc_define_integer_variable, &
                                   nc_define_real_variable, &
                                   nc_add_global_creation_time, &
                                   nc_add_global_attribute, &
                                   nc_define_dimension, nc_put_variable, &
                                   nc_sync, nc_enddef, nc_redef, nc_open_readonly, &
                                   nc_close, nc_variable_exists
use        chem_tables_mod, only : init_chem_tables, finalize_chem_tables, chem_convert_factor
!#!use       location_io_mod
use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, finalize_quad_interp, &
                                   quad_lon_lat_locate, quad_lon_lat_evaluate, &
                                   GRID_QUAD_IRREG_SPACED_REGULAR,  &
                                   QUAD_LOCATED_CELL_CENTERS
use     default_model_mod,  only : adv_1step, init_time, init_conditions, &
                                   nc_write_model_vars

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the dart code.

! routines in this list have code in this module
public :: static_init_model,                   &
          get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   & 
          shortest_time_between_assimilations, &
          nc_write_model_atts,                 &
          write_model_time,                    & 
          read_model_time,                     &
          end_model

! code for these routines are in other modules
public :: nc_write_model_vars,           &
          pert_model_copies,             & 
          adv_1step,                     &
          init_time,                     &
          init_conditions,               & 
          convert_vertical_obs,          & ! todo
          convert_vertical_state,        & ! doing
          get_close_obs,                 & ! todo
          get_close_state                  ! todo

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! maximum number of fields you can list to be perturbed
! to generate an ensemble if starting from a single state.
integer, parameter :: MAX_PERT = 100

! model_nml namelist variables and default values
character(len=256) :: cam_template_filename           = 'caminput.nc'
character(len=256) :: cam_phis_filename               = 'camphis.nc'
character(len=32)  :: vertical_localization_coord     = 'PRESSURE'
integer            :: assimilation_period_days        = 0
integer            :: assimilation_period_seconds     = 21600
logical            :: use_log_vertical_scale          = .false.
integer            :: no_assim_above_this_model_level = 4
logical            :: use_damping_ramp_at_model_top   = .false.  
integer            :: debug_level                     = 0
logical            :: suppress_grid_info_in_output    = .false.
logical            :: custom_routine_to_generate_ensemble = .false.
character(len=32)  :: fields_to_perturb(MAX_PERT)     = "QTY_TEMPERATURE"
real(r8)           :: perturbation_amplitude(MAX_PERT)= 0.00001_r8
logical            :: using_chemistry                 = .false.

! state_variables defines the contents of the state vector.
! each line of this input should have the form:
!
!    netcdf_variable_name, dart_quantity, clamp_min, clamp_max, update_variable
!
! all items must be strings (even if numerical values).
! for no clamping, use the string 'NA'
! to have the assimilation change the variable use 'UPDATE', else 'NO_UPDATE'

integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: num_state_table_columns = 5
character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLES * &
                                                   num_state_table_columns ) = ' '

namelist /model_nml/  &
   cam_template_filename,           &
   cam_phis_filename,               &
   vertical_localization_coord,     &
   state_variables,                 &
   assimilation_period_days,        &
   assimilation_period_seconds,     &
   no_assim_above_this_model_level, &
   use_damping_ramp_at_model_top,   &
   suppress_grid_info_in_output,    &
   custom_routine_to_generate_ensemble, &
   fields_to_perturb,               &
   perturbation_amplitude,          &
   debug_level

! global variables
character(len=512) :: string1, string2
logical, save      :: module_initialized = .false.

! domain id for the cam model.  this allows us access to all of the state structure
! info and is require for getting state variables.
integer :: domain_id

! this one must match the threed_sphere codes for VERTISxx
! default to pressure (2)
integer :: vert_local_coord = VERTISPRESSURE

!> Everything needed to describe a variable. Basically all the metadata from
!> a netCDF file is stored here as well as all the information about where
!> the variable is stored in the DART state vector.
!>

type cam_1d_array
   integer  :: nsize
   real(r8), allocatable :: vals(:)
   !>@todo FIXME do we need a string name here anymore?
end type

type cam_grid
   type(cam_1d_array) :: lon
   type(cam_1d_array) :: lat
   type(cam_1d_array) :: slon
   type(cam_1d_array) :: slat
   type(cam_1d_array) :: lev
   type(cam_1d_array) :: ilev
   type(cam_1d_array) :: gw
   type(cam_1d_array) :: hyai
   type(cam_1d_array) :: hybi
   type(cam_1d_array) :: hyam
   type(cam_1d_array) :: hybm
   type(cam_1d_array) :: P0
end type

type(cam_grid) :: grid_data


integer, parameter :: STAGGER_NONE = -1
integer, parameter :: STAGGER_U    =  1
integer, parameter :: STAGGER_V    =  2
integer, parameter :: STAGGER_W    =  3 
integer, parameter :: STAGGER_UV   =  4

type cam_stagger
   integer, allocatable :: qty_stagger(:)
end type

type(cam_stagger) :: grid_stagger

! Surface potential; used for calculation of geometric heights.
real(r8), allocatable :: phis(:, :)

! Horizontal interpolation code.  Need a handle for nonstaggered, U and V.
type(quad_interp_handle) :: interp_nonstaggered, &
                            interp_u_staggered, &
                            interp_v_staggered


contains


!-----------------------------------------------------------------------
! All the required interfaces are first.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> Called to do one time initialization of the model.
!> In this case, it reads in the grid information, the namelist
!> containing the variables of interest, where to get them, their size,
!> their associated DART KIND, etc.
!>
!> In addition to harvesting the model metadata (grid,
!> desired model advance step, etc.), it also fills a structure
!> containing information about what variables are where in the DART
!> framework.

subroutine static_init_model()

integer :: iunit, io
integer :: nfields

if ( module_initialized ) return

! The Plan:
!
! * read in the grid sizes from grid file
! * allocate space, and read in actual grid values
! * figure out model timestep
! * Compute the model size.

! Record version info
call register_module(source, revision, revdate)

module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(logfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('GREGORIAN')

call read_grid_info(cam_template_filename, grid_data)

!>@todo do we need to map_qtys here?
!>@todo do we need to set the model top related stuff here?

! set_cam_variable_info() fills var_names, kind_list, clamp_vals, update_list
! from the &model_mod_nml state_variables

call set_cam_variable_info(state_variables, nfields)

! convert from string to integer
call set_vert_localization(vertical_localization_coord, vert_local_coord)

! if you are going to have chemistry variables in the model state, set
! this namelist variable so we can initialize the proper tables
if (using_chemistry) call init_chem_tables()

end subroutine static_init_model


!-----------------------------------------------------------------------
!> Returns the size of the DART state vector (i.e. model) as an integer.
!>

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(domain_id)

end function get_model_size



!-----------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.
!>
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type the DART KIND at that index
!>

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type

! Local variables

integer  :: iloc, vloc, jloc
integer  :: myvarid, myqty, nd

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=myvarid, kind_index=myqty)

nd = get_num_dims(domain_id, myvarid)

location = get_location_from_index(iloc, jloc, vloc, myqty, nd)

! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
!> given the (i,j,k) indices into a field in the state vector,
!> and the quantity, and the dimensionality of the field (2d, 3d),
!> compute the location of that item.  

function get_location_from_index(i, j, k, q, nd)
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
integer, intent(in) :: q
integer, intent(in) :: nd
type(location_type) :: get_location_from_index

real(r8) :: slon_val
real(r8) :: use_vert_val
integer  :: use_vert_type

! full 3d fields are returned with lon/lat/level.
! 2d fields are either surface fields, or if they
! are column integrated values then they are 'undefined'
! in the vertical.

if (nd == 3) then
   use_vert_type = VERTISLEVEL
   use_vert_val  = real(k,r8)
else
   if (q == QTY_SURFACE_ELEVATION .or. q == QTY_SURFACE_PRESSURE) then
      use_vert_type = VERTISSURFACE
      use_vert_val  = MISSING_R8  ! could also be surface elevation:
      !use_vert_val  = phis(lon_index, lat_index) / gravity
   else
      use_vert_type = VERTISUNDEF
      use_vert_val  = MISSING_R8
   endif
endif

! the horizontal location depends on whether this quantity is on the
! mass point grid or staggered in either lat or lon.  

select case (grid_stagger%qty_stagger(q))
  case (STAGGER_U)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%slat%vals(j), &
                                          use_vert_val, use_vert_type)

  case (STAGGER_V)
   ! the first staggered longitude is negative.  dart requires lons
   ! be between 0 and 360.
   slon_val = grid_data%slon%vals(i)
   if (slon_val < 0) slon_val = slon_val + 360.0_r8
   get_location_from_index = set_location(slon_val, &
                                          grid_data%lat%vals(j), &
                                          use_vert_val, use_vert_type)
   
  !>@todo not sure what to do yet. ? +-1/2 ?
  case (STAGGER_W)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          use_vert_val - 0.5_r8, use_vert_type)
  ! no stagger - cell centers
  case default
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          use_vert_val, use_vert_type)

end select

end function get_location_from_index

!-----------------------------------------------------------------------
!> this routine should be called to get a value from an unstaggered grid
!> that corresponds to a staggered grid.  e.g. you need the surface pressurre
!> under a V wind point.

subroutine get_staggered_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, &
                                         lev_index, stagger_qty, vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: lon_index
integer,             intent(in) :: lat_index
integer,             intent(in) :: lev_index
integer,             intent(in) :: stagger_qty
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status

integer :: next_lat, prev_lon, stagger
real(r8) :: vals_bot(ens_size), vals_top(ens_size)

stagger = grid_stagger%qty_stagger(stagger_qty)

!> latitudes:  staggered value N is between N and (N + 1) on the unstaggered grid
!> longitudes: staggered value N is between N and (N - 1) on the unstaggered grid

select case (stagger)
  case (STAGGER_U)
   call stagger_quad(lon_index, lat_index, prev_lon, next_lat)

   call get_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, vals_bot, my_status)
   call get_values_from_qty(ens_handle, ens_size, qty, lon_index, next_lat,  lev_index, vals_top, my_status)

   vals = (vals_bot + vals_top) * 0.5_r8

  case (STAGGER_V)
   call stagger_quad(lon_index, lat_index, prev_lon, next_lat)

   call get_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, vals_bot, my_status)
   call get_values_from_qty(ens_handle, ens_size, qty, prev_lon,  lat_index, lev_index, vals_top, my_status)

   vals = (vals_bot + vals_top) * 0.5_r8

  ! no stagger - cell centers, or W stagger
  case default
   call get_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, vals, my_status)

end select

end subroutine get_staggered_values_from_qty


!-----------------------------------------------------------------------
!>
!> 

subroutine get_values_from_qty(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, vals, my_status)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: qty
integer,             intent(in) :: lon_index
integer,             intent(in) :: lat_index
integer,             intent(in) :: lev_index
real(r8),            intent(out) :: vals(ens_size)
integer,             intent(out) :: my_status

integer :: varid
integer(i8) :: state_indx

varid      = get_varid_from_kind(domain_id, qty)
if (varid < 0) then
   vals(:) = MISSING_R8
   my_status = 12
   return
endif

state_indx = get_dart_vector_index(lon_index, lat_index, lev_index, domain_id, varid)
vals(:)    = get_state(state_indx, ens_handle)

end subroutine get_values_from_qty


!-----------------------------------------------------------------------
!> this routine takes care of getting the actual state values.  get_state()
!> communicates with other MPI tasks and can be expensive.
!>
!> all ensemble members have the same horizontal location, but different 
!> ensemble members could have different vertical locations and
!> so be between different vertical layers.  this code tries to do the fewest
!> calls to get_state by only calling it for levels that are actually needed
!> and setting all members with those same levels in a single pass.
!> 

subroutine get_values_from_varid(ens_handle, ens_size, lon_index, lat_index, lev_index, &
                                 varid, vals, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,  intent(in)  :: ens_size
integer,  intent(in)  :: lon_index
integer,  intent(in)  :: lat_index
integer,  intent(in)  :: lev_index(ens_size)
integer,  intent(in)  :: varid
real(r8), intent(out) :: vals(ens_size)
integer,  intent(out) :: my_status(ens_size)

integer(i8) :: state_indx
integer  :: i, j
real(r8) :: temp_vals(ens_size) 
logical  :: member_done(ens_size)


! as we get the values for each ensemble member, we set the 'done' flag
! and a good return code. 
my_status(:) = 16
member_done(:) = .false.

! start with lev_index(1).  get the vals into a temp var.  
! run through 2-N. any other member that has the same level 
! set the outgoing values.  keep a separate flag for which 
! member(s) have been done.  skip to the next undone member 
! and get the state for that level.  repeat until all
! levels done.

do i=1, ens_size

   if (member_done(i)) cycle

   state_indx = get_dart_vector_index(lon_index, lat_index, lev_index(i), domain_id, varid)
   if (state_indx < 0) then
      !>@todo FIXME this shouldn't happen, right?  should it call the error handler instead?
      my_status(:) = 15
      return
   endif

   !print*, 'state_indx ', state_indx
   temp_vals(:) = get_state(state_indx, ens_handle)    ! all the ensemble members for level (i)
   ! if (lev_index(1) /= lev_index(2) .and. all(temp_vals > 0.00000001)) then
   ! endif

   ! start at i, because my ensemble member is clearly at this level.
   ! then continue on to see if any other members are also at this level.
   do j=i, ens_size
      if (member_done(j)) cycle
         print*, 'lev_index(1) ', i, lev_index(i)
         print*, 'lev_index(2) ', j, lev_index(j)
         print*, 'tem_values(:)', temp_vals(j)
      
      if (lev_index(j) == lev_index(i)) then
         vals(j) = temp_vals(j)
         member_done(j) = .true.
         my_status(j) = 0
         print*, 'tem_values(:)', temp_vals(j)
      endif
      
   enddo
enddo

end subroutine get_values_from_varid

!-----------------------------------------------------------------------
!> this is just for 3d fields

subroutine get_values_from_nonstate_fields(ens_handle, ens_size, lon_index, lat_index, &
                                           lev_index, obs_quantity, vals, my_status)
type(ensemble_type),  intent(in)  :: ens_handle
integer,              intent(in)  :: ens_size
integer,              intent(in)  :: lon_index
integer,              intent(in)  :: lat_index
integer,              intent(in)  :: lev_index(ens_size)
integer,              intent(in)  :: obs_quantity
real(r8),             intent(out) :: vals(ens_size)
integer,              intent(out) :: my_status(ens_size)

integer  :: imember
real(r8) :: vals_array(grid_data%lev%nsize,ens_size)

character(len=*), parameter :: routine = 'get_values_from_nonstate_fields:'

vals(:) = MISSING_R8

select case (obs_quantity)
   case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT) 
      if (obs_quantity == QTY_PRESSURE) then
         call cam_pressure_levels(ens_handle, ens_size, &
                                  lon_index, lat_index, grid_data%lev%nsize, &
                                  obs_quantity, vals_array, my_status)
      else
         call cam_height_levels(ens_handle, ens_size, &
                                lon_index, lat_index, grid_data%lev%nsize, &
                                obs_quantity, vals_array, my_status) 
      endif

      if (any(my_status /= 0)) return

      do imember=1,ens_size
         vals(imember) = vals_array(lev_index(imember), imember)
      enddo

   case (QTY_VERTLEVEL)
      vals(:)      = lev_index(:)
      my_status(:) = 0

   case default
      write(string1,*)'contact dart support. unexpected error for quantity ', obs_quantity
      call error_handler(E_MSG,routine,string1,source,revision,revdate)

end select

end subroutine get_values_from_nonstate_fields

!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> to the given location.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_qty the DART KIND of interest
!> @param interp_vals the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>
!> istatus = 2    asked to interpolate an unknown/unsupported quantity
!> istatus = 3    cannot locate horizontal quad
!> istatus = 4    cannot locate enclosing vertical levels
!> istatus = 5    cannot retrieve state vector values
!> istatus = 6    cannot do vertical interpolation for bottom layer
!> istatus = 7    cannot do vertical interpolation for top layer
!> istatus = 8    cannot interpolate in the quad to get the values
!> istatus = 9    cannot get vertical levels for an obs on model levels
!> istatus = 10   cannot get vertical levels for an obs on pressure levels
!> istatus = 11   cannot get vertical levels for an obs on height levels
!> istatus = 12   cannot get values from obs quantity
!> istatus = 13   can not interpolate values of this quantity
!> istatus = 14   can not interpolate values above a given level
!> istatus = 15   can not get indices from given state vector index
!> istatus = X
!> istatus = 99   unknown error - shouldn't happen
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_vals, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: interp_vals(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'model_interpolate:'

integer  :: varid, icorner, numdims, which_vert, status1
integer  :: lon_bot, lat_bot, lon_top, lat_top
integer  :: four_lons(4), four_lats(4)
integer  :: level_one_array(ens_size), status_array(ens_size)
integer  :: four_bot_levs(4, ens_size), four_top_levs(4, ens_size)
real(r8) :: lon_fract, lat_fract
real(r8) :: lon_lat_vert(3)
real(r8) :: botvals(ens_size), topvals(ens_size)
real(r8) :: four_vert_fracts(4, ens_size), quad_vals(4, ens_size)
type(quad_interp_handle) :: interp_handle

if ( .not. module_initialized ) call static_init_model


! Successful istatus is 0
interp_vals(:) = MISSING_R8
istatus(:)     = 99

! do we know how to interpolate this quantity?
call ok_to_interpolate(obs_qty, varid, status1)

if (status1 /= 0) then  
   ! the interpolation code should not print or error out if
   ! it is asked to interpolate a quantity it cannot do.  it is just
   ! a failed forward operator.  (some forward operators try one quantity
   ! and if that fails, then they compute what they need using different quantities.)
   ! if you set the debug level up high enough, we can print out info about
   ! failed forward operators.
   if(debug_level > 12) then
      write(string1,*)'did not find observation quantity ', obs_qty, ' in the state vector'
      call error_handler(E_MSG,routine,string1,source,revision,revdate)
   endif
   istatus(:) = 2   ! this quantity not in the state vector
   return
endif

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location)) 

! get the grid handle for the right staggered grid
interp_handle = get_interp_handle(obs_qty)

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2) , &
                         four_lons, four_lats, lon_fract, lat_fract, status1)

if (status1 /= 0) then
   istatus(:) = 3  ! cannot locate enclosing horizontal quad
   return
endif

! need to consider the case for 2d vs 3d variables
numdims = get_dims_from_qty(obs_qty, varid)

!>@todo FIXME need to be refactored
! if variable is in the state
!    do 2d and 3d

! if variable is not in the state
!    do 2d and 3d

! now here potentially we have different results for different
! ensemble members.  the things that can vary are dimensioned by ens_size.

if (numdims > 2 ) then

   ! build 4 columns to find vertical level numbers
   do icorner=1, 4
      call find_vertical_levels(state_handle, ens_size, &
                                four_lons(icorner), four_lats(icorner), lon_lat_vert(3), &
                                which_vert, obs_qty, &
                                four_bot_levs(icorner, :), four_top_levs(icorner, :), &
                                four_vert_fracts(icorner, :), status_array)
      

      !>@todo FIXME should we let the process continue if at least one
      !>member has failed?  pro: save work  con: don't get forward operator
      !>values for members that could compute them
      !>(this is true for all the subsequent returns from this routine)

      !if (any(four_bot_levs(:, 1) /= four_bot_levs(:, 2))) then
      !   print*, 'four_bot_levs(:, 1)', four_bot_levs(:, 1)
      !   print*, 'four_bot_levs(:, 2)', four_bot_levs(:, 2)
      !endif
      !
      !if (any(four_bot_levs(:, 1) /= four_bot_levs(:, 2))) then
      !   print*, 'four_top_levs(:, 1)', four_top_levs(:, 1)
      !   print*, 'four_top_levs(:, 2)', four_top_levs(:, 2)
      !endif
      if (any(status_array /= 0)) then
         istatus(:) = 4   ! cannot locate enclosing vertical levels  !>@todo FIXME use where statements?
         return
      endif
      
      ! if we are avoiding assimilating obs above a given level, test here and return
      ! if any of the bottom corners are above the limit (meaning the obs is at least
      ! in the layer above the given cutoff.)
      !
      ! level 1 is top, so test that the level numbers are *smaller* than the limit.
      ! (meaning the obs is above the given limit in at least one ensemble member)
      !print*, 'levs', four_bot_levs(icorner,:), no_assim_above_this_model_level

      if (no_assim_above_this_model_level > 0) then
         if (any(four_bot_levs(icorner,:) <= no_assim_above_this_model_level)) then
            print*, 'istatus is 14'
            istatus(:) = 14
            return
         endif
      endif
   enddo
   
   ! we have all the indices and fractions we could ever want.
   ! now get the data values at the bottom levels, the top levels, 
   ! and do vertical interpolation to get the 4 values in the columns.
   ! the final horizontal interpolation will happen later.
      
   if (varid > 0) then

      call get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                four_bot_levs, four_top_levs, four_vert_fracts, &
                                varid, quad_vals, status_array)

   else ! get 3d special variables in another ways ( like QTY_PRESSURE )

      call get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                   four_bot_levs, four_top_levs, four_vert_fracts, &
                                   obs_qty, quad_vals, status_array)

   endif

      ! if we are avoiding assimilating obs above a given level, test here and return
      ! if any of the bottom corners are above the limit (meaning the obs is at least
      ! in the layer above the given cutoff.)
      !
      ! level 1 is top, so test that the level numbers are *smaller* than the limit.
      ! (meaning the obs is above the given limit in at least one ensemble member)
   do icorner=1, 4
      !print*, 'levs', four_bot_levs(icorner,:), no_assim_above_this_model_level

      if (no_assim_above_this_model_level > 0) then
         if (any(four_bot_levs(icorner,:) <= no_assim_above_this_model_level)) then
            print*, 'istatus is 15'
            istatus(:) = 15
            return
         endif
      endif
   enddo

   if (any(status_array /= 0)) then
      istatus(:) = status_array(:)
      return
   endif

else ! 2 dimensional variables

   if (varid > 0) then
      level_one_array(:) = 1
      do icorner=1, 4
         call get_values_from_varid(state_handle,  ens_size, & 
                                    four_lons(icorner), four_lats(icorner), &
                                    level_one_array, varid, quad_vals(icorner,:), status_array)
      enddo
   else ! special 2d case
      do icorner=1, 4
         call get_quad_corners(ens_size, four_lons(icorner), four_lats(icorner), &
                               obs_qty, obs_qty, quad_vals(icorner,:), status1)
      enddo
   endif

endif

! do the horizontal interpolation for each ensemble member
call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, ens_size, &
                           quad_vals, interp_vals, status_array)
if (any(status_array /= 0)) then
   istatus(:) = 8   ! cannot evaluate in the quad
   return
endif

! all interp values should be set by now.  set istatus
istatus(:) = 0

end subroutine model_interpolate

!-----------------------------------------------------------------------
!>

subroutine get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                 four_bot_levs, four_top_levs, four_vert_fracts, &
                                 varid, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: four_bot_levs(4, ens_size), four_top_levs(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: varid
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: botvals(ens_size), topvals(ens_size)

character(len=*), parameter :: routine = 'get_four_state_values:'

do icorner=1, 4
   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_bot_levs(icorner, :), varid, botvals, &
                              my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 5   ! cannot retrieve bottom values
      return
   endif

   call get_values_from_varid(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_top_levs(icorner, :), varid, topvals, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 6   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, botvals, topvals, four_vert_fracts(icorner, :), &
                    quad_vals(icorner, :), my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 7   ! cannot do vertical interpolation
      return
   endif
enddo


end subroutine get_four_state_values

!-----------------------------------------------------------------------
!>

subroutine get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                 four_bot_levs, four_top_levs, four_vert_fracts, &
                                 obs_qty, quad_vals, my_status)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
integer,             intent(in) :: four_lons(4), four_lats(4)
integer,             intent(in) :: four_bot_levs(4, ens_size), four_top_levs(4, ens_size)
real(r8),            intent(in) :: four_vert_fracts(4, ens_size)
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: quad_vals(4, ens_size) !< array of interpolated values
integer,            intent(out) :: my_status(ens_size)

integer  :: icorner
real(r8) :: botvals(ens_size), topvals(ens_size)

character(len=*), parameter :: routine = 'get_four_nonstate_values:'

do icorner=1, 4
   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_bot_levs(icorner, :), obs_qty, botvals, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 5   ! cannot retrieve bottom values
      return
   endif

   call get_values_from_nonstate_fields(state_handle,  ens_size, &
                              four_lons(icorner), four_lats(icorner), &
                              four_top_levs(icorner, :), obs_qty, topvals, my_status)
   if (any(my_status /= 0)) then
      my_status(:) = 6   ! cannot retrieve top values
      return
   endif

   call vert_interp(ens_size, botvals, topvals, four_vert_fracts(icorner, :), &
                    quad_vals(icorner, :), my_status)

   if (any(my_status /= 0)) then
      my_status(:) = 7   ! cannot do vertical interpolation
      return
   endif
enddo

end subroutine get_four_nonstate_values

!-----------------------------------------------------------------------
!> figure out whether this is a 2d or 3d field based on the quantity.
!> if this field is in the state vector, use the state routines.
!> if it's not, there are cases for known other quantities we can
!> interpolate and return.  add any new non-state fields here.

function get_dims_from_qty(obs_quantity, var_id)
integer, intent(in) :: obs_quantity
integer, intent(in) :: var_id
integer :: get_dims_from_qty

character(len=*), parameter :: routine = 'get_dims_from_qty:'

if (var_id > 0) then
   get_dims_from_qty = get_num_dims(domain_id,var_id)
else
   select case (obs_quantity)
      case (QTY_SURFACE_ELEVATION)
         get_dims_from_qty = 2
      case (QTY_PRESSURE)
         get_dims_from_qty = 3
      case default 
         write(string1, *) 'we can not interpolate qty', obs_quantity, &
                           ' if the dimension is not known'
         call error_handler(E_ERR,routine, string1,source,revision,revdate)
    end select
endif

end function get_dims_from_qty

!-----------------------------------------------------------------------
!> return 0 (ok) if we know how to interpolate this quantity.
!> if it is a field in the state, return the variable id from
!> the state structure.  if not in the state, varid will return -1

subroutine ok_to_interpolate(obs_qty, varid, my_status)
integer, intent(in)  :: obs_qty
integer, intent(out) :: varid
integer, intent(out) :: my_status

! See if the state contains the obs quantity 
varid = get_varid_from_kind(domain_id, obs_qty)

! in the state vector
if (varid > 0) then
   my_status = 0
   return
endif
   

! add any quantities that can be interpolated to this list if they
! are not in the state vector.
select case (obs_qty)
   case (QTY_SURFACE_ELEVATION, &
         QTY_PRESSURE, &
         QTY_SPECIFIC_HUMIDITY, &
         QTY_GEOPOTENTIAL_HEIGHT,  &
         QTY_PRECIPITABLE_WATER, &
         QTY_GEOMETRIC_HEIGHT, &
         QTY_VERTLEVEL)
      my_status = 0
   case default
      my_status = 13
end select


end subroutine ok_to_interpolate


!-----------------------------------------------------------------------
!>
!>  This is for 2d special observations quantities not in the state

subroutine get_quad_corners(ens_size, lon_index, lat_index, obs_quantity, stagger_qty, vals, my_status)
integer,  intent(in) :: ens_size
integer,  intent(in) :: lon_index
integer,  intent(in) :: lat_index
integer,  intent(in) :: obs_quantity
integer,  intent(in) :: stagger_qty
real(r8), intent(out) :: vals(ens_size) 
integer,  intent(out) :: my_status

character(len=*), parameter :: routine = 'get_quad_corners'

integer :: stagger, prev_lon, next_lat
real(r8) :: vals_bot(ens_size), vals_top(ens_size)

stagger = grid_stagger%qty_stagger(stagger_qty)

!>@todo FIXME we need to look at the stagger of the obs_quantity here
!>(don't we?)

select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)

     select case (stagger)
       case (STAGGER_U)
          call stagger_quad(lon_index, lat_index, prev_lon, next_lat)
          vals_bot(:) = phis(lon_index, lat_index) 
          vals_top(:) = phis(lon_index, next_lat) 
     
        vals = (vals_bot + vals_top) * 0.5_r8 
     
       case (STAGGER_V)
          call stagger_quad(lon_index, lat_index, prev_lon, next_lat)
          vals_bot(:) = phis(lon_index, lat_index) 
          vals_top(:) = phis(prev_lon,  lat_index) 
     
        vals = (vals_bot + vals_top) * 0.5_r8
     
       ! no stagger - cell centers, or W stagger
       case default
  
        vals = phis(lon_index, lat_index)
  
     end select
    
     !>@todo FIXME:
     ! should this be using gravity at the given latitude? 
     vals = vals / gravity

   case default 
      write(string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)

end select

my_status = 0

end subroutine get_quad_corners


!-----------------------------------------------------------------------
!>  interpolate in the vertical between 2 arrays of items.

subroutine vert_interp(nitems, botvals, topvals, vert_fracts, out_vals, my_status)
integer,  intent(in)  :: nitems
real(r8), intent(in)  :: botvals(nitems)
real(r8), intent(in)  :: topvals(nitems)
real(r8), intent(in)  :: vert_fracts(nitems)
real(r8), intent(out) :: out_vals(nitems)
integer,  intent(out) :: my_status(nitems)

!>@todo FIXME does this need a status return if it cannot fail?

! vert_fracts: 1 is 100% of the bottom level and 
!              0 is 100% of the top level

out_vals(:) = (botvals(:)* vert_fracts(:)) + (topvals(:) * (1.0_r8-vert_fracts(:)))
my_status(:) = 0

end subroutine vert_interp

!-----------------------------------------------------------------------
!> given lon/lat indices, add one to lat and subtract one from lon
!> check for wraparound in lon, and north pole at lat.
!> intent is that you give the indices into the staggered grid
!> and the return values are the indices in the original unstaggered
!> grid that you need to compute the midpoints for the staggers.
!>@todo FIXME this needs a picture or ascii art

subroutine stagger_quad(lon_index, lat_index, prev_lon, next_lat)
integer, intent(in)  :: lon_index
integer, intent(in)  :: lat_index
integer, intent(out) :: prev_lon
integer, intent(out) :: next_lat

   next_lat = lat_index+1
   if (next_lat > grid_data%lat%nsize) next_lat = grid_data%lat%nsize

   prev_lon = lon_index-1
   if (prev_lon < 1) prev_lon = grid_data%lon%nsize

end subroutine stagger_quad


!-----------------------------------------------------------------------
!> given a lon/lat index number, a quantity and a vertical value and type,
!> return which two levels these are between and the fraction across.
!> 

subroutine find_vertical_levels(ens_handle, ens_size, lon_index, lat_index, vert_val, &
                                which_vert, obs_qty, bot_levs, top_levs, vert_fracts, my_status)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index 
integer,             intent(in)  :: lat_index
real(r8),            intent(in)  :: vert_val
integer,             intent(in)  :: which_vert
integer,             intent(in)  :: obs_qty
integer,             intent(out) :: bot_levs(ens_size)
integer,             intent(out) :: top_levs(ens_size)
real(r8),            intent(out) :: vert_fracts(ens_size)
integer,             intent(out) :: my_status(ens_size)

character(len=*), parameter :: routine = 'find_vertical_levels:'

integer  :: bot1, top1, imember, nlevels, level_one, status1, k
real(r8) :: fract1
real(r8) :: surf_pressure (  ens_size )
real(r8) :: pressure_array( grid_data%lev%nsize, ens_size )
real(r8) :: height_array  ( grid_data%lev%nsize, ens_size )

! assume the worst
bot_levs(:)    = MISSING_I
top_levs(:)    = MISSING_I
vert_fracts(:) = MISSING_R8
my_status(:)   = 98

!>@todo FIXME we need nlevels everywhere - should we have a global var for it?

! number of vertical levels (midlayer points)
nlevels   = grid_data%lev%nsize
level_one = 1

select case (which_vert)

   case(VERTISPRESSURE)
      ! construct a pressure column here and find the model levels
      ! that enclose this value
      call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                         lon_index, lat_index, level_one, obs_qty, &
                                         surf_pressure, status1)

      call build_cam_pressure_column(ens_size, surf_pressure, nlevels, pressure_array)

      do imember=1, ens_size
         call pressure_to_level(nlevels, pressure_array(:, imember), vert_val, &
                                bot_levs(imember), top_levs(imember), &
                                vert_fracts(imember), my_status(imember))

      enddo

   case(VERTISHEIGHT)
      ! construct a height column here and find the model levels
      ! that enclose this value
      !@>todo put in arguments and write height_to_level
      call cam_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, obs_qty, &
                             height_array, my_status)
      if (any(my_status /= 0)) return   !>@todo FIXME let successful members continue?
      do imember=1, ens_size
         !>@todo FIXME somewhere cull out unique levels and only get_state() for those. (see wrf)
         call height_to_level(nlevels, height_array(:, imember), vert_val, &
                              bot_levs(imember), top_levs(imember), vert_fracts(imember), my_status(imember))
         do k = 1,ens_size
            print*, 'bot_levs(i), top_levs(i), vert_fracts()', bot_levs(imember), top_levs(imember), vert_fracts(imember)
         enddo
      enddo
      

      !write(string1, *) 'we have not written the code yet for vertical type: ', which_vert
      !call error_handler(E_ERR,routine,string1,source,revision,revdate)

   case(VERTISLEVEL)
      ! this routine returns false if the level number is out of range.
      if (range_set(vert_val, nlevels, bot1, top1, fract1)) then
         my_status(:) = 8
         return
      endif

      ! because we're given a model level as input, all the ensemble
      ! members have the same outgoing values.
      bot_levs(:) = bot1
      top_levs(:) = top1
      vert_fracts(:) = fract1

   ! 2d fields
   case(VERTISSURFACE)
   case(VERTISUNDEF)  
      !>@todo FIXME  OR  all levels = 1?  for a 2d field and get_state_index()?
      bot_levs(:) = nlevels
      top_levs(:) = nlevels - 1
      vert_fracts(:) = 1.0_r8

   case default
      !>@todo FIXME: do nothing or error out here?

      write(string1, *) 'unsupported vertical type: ', which_vert
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
      
end select

end subroutine find_vertical_levels

!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

subroutine cam_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, height_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: qty
real(r8),            intent(out) :: height_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer  :: k, level_one, imember, status1
real(r8) :: surface_elevation(1)
real(r8) :: temperature(ens_size), specific_humidity(ens_size), surface_pressure(ens_size)
real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom

!>@todo this should come from a model specific constant module.
!> the forward operators and model_mod should use it.
real(r8), parameter :: rd = 287.05_r8 ! dry air gas constant
real(r8), parameter :: rv = 461.51_r8 ! wet air gas constant
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

! this is for surface obs
level_one = 1

!@>todo make into a subroutine get_val or something similar
! get the surface pressure from the ens_handle
call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                   lon_index, lat_index, level_one, qty, surface_pressure, status1)


! get the surface elevation from the phis, including stagger if needed
call get_quad_corners(1, lon_index, lat_index, QTY_SURFACE_ELEVATION, qty, surface_elevation, status1)

print*, 'surface_elevation ', surface_elevation

do k = 1, nlevels
   ! temperature
   call get_staggered_values_from_qty(ens_handle, ens_size, QTY_TEMPERATURE, &
                                     lon_index, lat_index, k, qty, temperature, status1)

   ! specific humidity
   call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, &
                                     lon_index, lat_index, k, qty, specific_humidity, status1)
   
   !>tv == virtual temperature.
   tv(k,:) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))

enddo

! compute the height columns for each ensemble member
do imember = 1, ens_size
   call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), tv(:, imember), &
                      height_array(:, imember))  ! can pass in variable_r here
enddo

! convert entire array to geometric height (from potential height)
call gph2gmh(height_array, grid_data%lat%vals(lat_index))

if (debug_level > 100) then
   do k = 1,nlevels
     do imember = 1, ens_size
       print *, "member, level, height: ", imember, k, height_array(k, imember)
     enddo
   enddo
endif

my_status(:) = 0

end subroutine cam_height_levels

!-----------------------------------------------------------------------
!> Compute the pressures at the layer midpoints

subroutine build_cam_pressure_column(ens_size, surface_pressure, n_levels, pressure_array)

integer,            intent(in)  :: ens_size
real(r8),           intent(in)  :: surface_pressure(:)   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: pressure_array(:,:)

integer :: j, k

! Set midpoint pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do j=1, ens_size
   do k=1,n_levels
      pressure_array(k, j) = grid_data%hyam%vals(k) * grid_data%P0%vals(1) + &
                             grid_data%hybm%vals(k) * surface_pressure(j)
   enddo
enddo

end subroutine build_cam_pressure_column


!-----------------------------------------------------------------------
!> Compute column of pressures at the layer midpoints for the given 
!> surface pressure. 
!>
!>@todo FIXME unlike some other things - you could pass in an upper and lower
!>level number and compute the pressure only at the levels between those.
!>this isn't a column that has to be built from the surface up or top down.
!>each individual pressure can be computed independently given the surface pressure.

subroutine cam_p_col_midpts(surface_pressure, n_levels, pressure_array)

real(r8),           intent(in)  :: surface_pressure   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: pressure_array(n_levels)

integer :: k

! Set midpoint pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do k=1, n_levels
   pressure_array(k) = grid_data%hyam%vals(k) * grid_data%P0%vals(1) + &
                       grid_data%hybm%vals(k) * surface_pressure
enddo

end subroutine cam_p_col_midpts

!-----------------------------------------------------------------------
!> Compute columns of pressures at the layer interfaces for the given number 
!> of surface pressures. 
!>
!>@todo FIXME see comment above in cam_p_col_midpts()

subroutine cam_p_col_intfcs( surface_pressure, n_levels, pressure_array)

real(r8),           intent(in)  :: surface_pressure   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: pressure_array(n_levels)

integer :: k

! Set interface pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do k=1, n_levels
   pressure_array(k) = grid_data%hyai%vals(k) * grid_data%P0%vals(1) + &
                       grid_data%hybi%vals(k) * surface_pressure
enddo

end subroutine cam_p_col_intfcs

!-----------------------------------------------------------------------
!> Compute columns of heights at the layer midpoints for the given number 
!> of surface pressures and surface elevations.  (should this one compute
!> those itself?  if so give it the lon/lat indx and staggers.)

!> since it's no more work, should we compute the interface heights as well here?
!> return both, or make the arrays optional and return one or both depending
!> on what's asked for?

subroutine cam_h_col_midpts(num_cols, surface_pressure, n_levels, height_array)

integer,            intent(in)  :: num_cols
real(r8),           intent(in)  :: surface_pressure(num_cols)   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: height_array(n_levels, num_cols)  ! in meters

integer :: j, k
real(r8) :: surface_elev
real(r8) :: p_midpts(n_levels,   num_cols)  ! in pascals, layer midpoints
real(r8) :: p_intfcs(n_levels+1, num_cols)  ! in pascals, layer interfaces

! Set midpoint heights.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

! unlike pressure, we need to start at the surface and work our way up.
! we can't compute a height in the middle of the column.

do j=1, num_cols
   call cam_p_col_midpts(surface_pressure(j), n_levels, p_midpts(:,j))
   call cam_p_col_intfcs(surface_pressure(j), n_levels, p_intfcs(:,j))
   ! get surface elevation at this location...
   surface_elev = MISSING_R8
   height_array(1, 1) = surface_elev 
   do k=1, n_levels
      ! add on lower half of this layer and record the height of the midpoint,
      ! then add the upper half of this layer.  each layer has a different 
      ! virtual temperature (at the midpt) so we can't compute midpoint to midpoint.
      height_array(k, j) = -1
      !height_array(k, j) = grid_data%hyam%vals(k) * grid_data%P0%vals(1) + &
      !                     grid_data%hybm%vals(k) * surface_pressure(j)
   enddo
enddo

end subroutine cam_h_col_midpts

!-----------------------------------------------------------------------
!> Compute columns of heights at the layer interfaces for the given number 
!> of surface pressures. 

subroutine cam_h_col_intfcs(num_cols, surface_pressure, n_levels, height_array)

integer,            intent(in)  :: num_cols
real(r8),           intent(in)  :: surface_pressure(num_cols)   ! in pascals
integer,            intent(in)  :: n_levels
real(r8),           intent(out) :: height_array(n_levels, num_cols)  ! in meters

integer :: j, k

! Set interface heights.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do j=1, num_cols
   do k=1, n_levels
      height_array(k, j) = -1
      !height_array(k, j) = grid_data%hyai%vals(k) * grid_data%P0%vals(1) + &
      !                     grid_data%hybi%vals(k) * surface_pressure(j)
   enddo
enddo

end subroutine cam_h_col_intfcs

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> 1 is top, N is bottom, bot is the lower level, top is the upper level
!> so top value will be smaller than bot.  fract = 0 is the full top,
!> fract = 1 is the full bot.  return non-zero if value outside the valid range.

subroutine pressure_to_level(nlevels, pressures, p_val, &
                             bot_lev, top_lev, fract, my_status)

integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: pressures(:)
real(r8), intent(in)  :: p_val
integer,  intent(out) :: bot_lev
integer,  intent(out) :: top_lev
real(r8), intent(out) :: fract  
integer,  intent(out) :: my_status

character(len=*), parameter :: routine = 'pressure_to_level:'

integer :: this_lev

bot_lev = MISSING_I
top_lev = MISSING_I
fract   = MISSING_R8

if (p_val < pressures(1) .or. p_val > pressures(nlevels)) then
   my_status = 10
   return
endif

! search from the top down.  smaller pressures are higher
! than larger ones.  in theory you can never fall out of
! this loop without setting the top/bot because we've tested
! already for p_val out of range.
levloop: do this_lev = 2, nlevels
   if (p_val >= pressures(this_lev)) cycle levloop

   top_lev = this_lev - 1
   bot_lev = this_lev
   ! this is where we do the vertical fraction in either linear or log scale
   if (use_log_vertical_scale) then
      fract = (log(p_val) - log(pressures(top_lev))) / (log(pressures(bot_lev)) - log(pressures(top_lev)))
   else
      fract = (p_val - pressures(top_lev)) / (pressures(bot_lev) - pressures(top_lev))
   endif
   my_status = 0
   return
enddo levloop

! you shouldn't get here
if (bot_lev == MISSING_I) then
  write(string1,*) 'should not happen - contact DART support'
  write(string2,*) 'pressure value ', p_val, ' was not found in pressure column'
  call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine pressure_to_level

!-----------------------------------------------------------------------
!>
!> return the level indices and fraction across the level.
!> 1 is top, N is bottom, bot is the lower level, top is the upper level
!> so top value will be larger than bot.  fract = 0 is the full top,
!> fract = 1 is the full bot.  return non-zero if value outside the valid range.
!>

subroutine height_to_level(nlevels, heights, h_val, &
                           bot_lev, top_lev, fract, my_status)

integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: heights(:)
real(r8), intent(in)  :: h_val
integer,  intent(out) :: bot_lev
integer,  intent(out) :: top_lev
real(r8), intent(out) :: fract  
integer,  intent(out) :: my_status

character(len=*), parameter :: routine = 'height_to_level:'

integer :: this_lev

bot_lev = MISSING_I
top_lev = MISSING_I
fract   = MISSING_R8

if (h_val > heights(1) .or. h_val < heights(nlevels)) then
   my_status = 11
   return
endif

! search from the top down.  larger heights are higher
! than smaller ones.  in theory you can never fall out of
! this loop without setting the top/bot because we've tested
! already for h_val out of range.
levloop: do this_lev = 2, nlevels
   if (h_val < heights(this_lev)) cycle levloop
   top_lev = this_lev - 1
   bot_lev = this_lev
   fract = (h_val - heights(top_lev)) / (heights(bot_lev) - heights(top_lev))
   my_status = 0
   return
enddo levloop

! you shouldn't get here
if (bot_lev == MISSING_I) then
  write(string1,*) 'should not happen - contact DART support'
  write(string2,*) 'height value ', h_val, ' was not found in height column'
  call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine height_to_level

!-----------------------------------------------------------------------
!> in cam level 1 is at the model top, level N is the lowest level
!> our convention in this code is:  between levels a fraction of 0
!> is 100% the top level, and fraction of 1 is 100% the bottom level.
!> the top level is always closer to the model top and so has a *smaller*
!> level number than the bottom level. 

function range_set(vert_value, valid_range, bot, top, fract)
real(r8), intent(in)  :: vert_value
integer,  intent(in)  :: valid_range
integer,  intent(out) :: bot
integer,  intent(out) :: top
real(r8), intent(out) :: fract
logical               :: range_set

integer :: integer_level
real(r8) :: fract_level

! be a pessimist, then you're never disappointed
range_set = .false.
bot = MISSING_I
top = MISSING_I
fract = MISSING_R8

integer_level = floor(vert_value)         ! potential top
fract_level = vert_value - integer_level 

! cam levels start at the top so level 1 is
! the highest level and increase on the way down.

!>@todo FIXME might want to add debugging print

!>might want to allow extrapolation - which means
!>allowing out of range values here and handling
!>them correctly in the calling and vert_interp() code.

! out of range checks
if (vert_value < 1.0_r8 .or. vert_value > valid_range) return

if (vert_value /= valid_range) then
   top = integer_level
   bot = top + 1
   fract = fract_level
else   
   ! equal to the bottom level
   top = integer_level - 1
   bot = integer_level
   fract = 1.0_r8
endif

range_set = .true.

end function range_set


!-----------------------------------------------------------------------
!> based on the stagger that corresponds to the given quantity,
!> return the handle to the interpolation grid


function get_interp_handle(obs_quantity)
integer, intent(in)      :: obs_quantity
type(quad_interp_handle) :: get_interp_handle

character(len=*), parameter :: routine = 'get_interp_handle:'

select case (grid_stagger%qty_stagger(obs_quantity))
   case ( STAGGER_U ) 
      get_interp_handle = interp_u_staggered
   case ( STAGGER_V ) 
      get_interp_handle = interp_v_staggered
   case ( STAGGER_NONE )
      get_interp_handle = interp_nonstaggered
   case ( STAGGER_W ) 
      write(string1,*) 'w stagger -- not supported yet'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   case ( STAGGER_UV ) 
      write(string1,*) 'uv stagger -- not supported yet'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   case default
      write(string1,*) 'unknown stagger -- this should never happen'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end select
                      
end function get_interp_handle

!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>
!>

function shortest_time_between_assimilations()

character(len=*), parameter :: routine = 'shortest_time_between_assimilations:'

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = set_time(assimilation_period_seconds, &
                                               assimilation_period_days)

write(string1,*)'assimilation period is ',assimilation_period_days,   ' days ', &
                                          assimilation_period_seconds,' seconds'
call error_handler(E_MSG,routine,string1,source,revision,revdate)

end function shortest_time_between_assimilations




!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! deallocate arrays from grid and anything else

call free_cam_grid(grid_data)

call finalize_quad_interp(interp_nonstaggered)
call finalize_quad_interp(interp_u_staggered)
call finalize_quad_interp(interp_v_staggered)

if (using_chemistry) call finalize_chem_tables()

end subroutine end_model


!-----------------------------------------------------------------------
!>
!> Writes the model-specific attributes to a DART 'diagnostic' netCDF file.
!> This includes coordinate variables and some metadata, but NOT the
!> actual DART state.
!>
!> @param ncid    the netCDF handle of the DART diagnostic file opened by
!>                assim_model_mod:init_diag_output

subroutine nc_write_model_atts(ncid, dom_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: dom_id    ! not used since there is only one domain

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

character(len=*), parameter :: routine = 'nc_write_model_atts'

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------
call nc_redef(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "CAM")

! this option is for users who want the smallest output
! or diagnostic files - only the state vector data will
! be written.   otherwise, if you want to plot this data
! the rest of this routine writes out enough grid info
! to make the output file look like the input.
if (suppress_grid_info_in_output) then
   call nc_enddef(ncid)
   return
endif

!----------------------------------------------------------------------------
! Output the grid variables.
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_define_dimension(ncid, 'lon',  grid_data%lon%nsize,  routine)
call nc_define_dimension(ncid, 'lat',  grid_data%lat%nsize,  routine)
call nc_define_dimension(ncid, 'slon', grid_data%slon%nsize, routine)
call nc_define_dimension(ncid, 'slat', grid_data%slat%nsize, routine)
call nc_define_dimension(ncid, 'lev',  grid_data%lev%nsize,  routine)
call nc_define_dimension(ncid, 'ilev', grid_data%ilev%nsize, routine)
call nc_define_dimension(ncid, 'gw',   grid_data%gw%nsize,   routine)
call nc_define_dimension(ncid, 'hyam', grid_data%hyam%nsize, routine)
call nc_define_dimension(ncid, 'hybm', grid_data%hybm%nsize, routine)
call nc_define_dimension(ncid, 'hyai', grid_data%hyai%nsize, routine)
call nc_define_dimension(ncid, 'hybi', grid_data%hybi%nsize, routine)

!----------------------------------------------------------------------------
! Create the Coordinate Variables and the Attributes
! The contents will be written in a later block of code.
!----------------------------------------------------------------------------

! U,V Grid Longitudes
call nc_define_real_variable(     ncid, 'lon', (/ 'lon' /),                 routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'long_name', 'longitude',    routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'units',     'degrees_east', routine)


call nc_define_real_variable(     ncid, 'slon', (/ 'slon' /),                       routine)
call nc_add_attribute_to_variable(ncid, 'slon', 'long_name', 'staggered longitude', routine)
call nc_add_attribute_to_variable(ncid, 'slon', 'units',     'degrees_east',        routine)

! U,V Grid Latitudes
call nc_define_real_variable(     ncid, 'lat', (/ 'lat' /),                  routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'long_name', 'latitude',      routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'units',     'degrees_north', routine)


call nc_define_real_variable(     ncid, 'slat', (/ 'slon' /),                      routine)
call nc_add_attribute_to_variable(ncid, 'slat', 'long_name', 'staggered latitude', routine)
call nc_add_attribute_to_variable(ncid, 'slat', 'units',     'degrees_north',      routine)

! Vertical Grid Latitudes
call nc_define_real_variable(     ncid, 'lev', (/ 'lev' /),                                                     routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'long_name', 'hybrid level at midpoints (1000*(A+B))',           routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'units',          'level',                                       routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'positive',       'down',                                        routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'standard_name',  'atmosphere_hybrid_sigma_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'lev', 'formula_terms',  'a: hyam b: hybm p0: P0 ps: PS',                routine)


call nc_define_real_variable(     ncid, 'ilev', (/ 'ilev' /),                                                    routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'long_name',      'hybrid level at interfaces (1000*(A+B))',     routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'units',          'level',                                       routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'positive',       'down',                                        routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'standard_name',  'atmosphere_hybrid_sigma_pressure_coordinate', routine)
call nc_add_attribute_to_variable(ncid, 'ilev', 'formula_terms',  'a: hyai b: hybi p0: P0 ps: PS',               routine)

! Hybrid Coefficients
call nc_define_real_variable(     ncid, 'hyam', (/ 'lev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hyam', 'long_name', 'hybrid A coefficient at layer midpoints', routine)

call nc_define_real_variable(     ncid, 'hybm', (/ 'lev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hybm', 'long_name', 'hybrid B coefficient at layer midpoints', routine)


call nc_define_real_variable(     ncid, 'hyai', (/ 'ilev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hyai', 'long_name', 'hybrid A coefficient at layer interfaces', routine)


call nc_define_real_variable(     ncid, 'hybi', (/ 'ilev' /),                                            routine)
call nc_add_attribute_to_variable(ncid, 'hybi', 'long_name', 'hybrid B coefficient at layer interfaces', routine)

! Gaussian Weights
call nc_define_real_variable(     ncid, 'gw', (/ 'lat' /),                  routine)
call nc_add_attribute_to_variable(ncid, 'gw', 'long_name', 'gauss weights', routine)

call nc_define_real_variable(ncid, 'P0', 0, routine)
call nc_add_attribute_to_variable(ncid, 'P0', 'long_name', 'reference pressure', routine)
call nc_add_attribute_to_variable(ncid, 'P0', 'units',     'Pa',                 routine)

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_enddef(ncid)

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_put_variable(ncid, 'lon',  grid_data%lon%vals,  routine)
call nc_put_variable(ncid, 'lat',  grid_data%lat%vals,  routine)
call nc_put_variable(ncid, 'slon', grid_data%slon%vals, routine)
call nc_put_variable(ncid, 'slat', grid_data%slat%vals, routine)
call nc_put_variable(ncid, 'lev',  grid_data%lev%vals,  routine)
call nc_put_variable(ncid, 'ilev', grid_data%ilev%vals, routine)
call nc_put_variable(ncid, 'gw',   grid_data%gw%vals,   routine)
call nc_put_variable(ncid, 'hyam', grid_data%hyam%vals, routine)
call nc_put_variable(ncid, 'hybm', grid_data%hybm%vals, routine)
call nc_put_variable(ncid, 'hyai', grid_data%hyai%vals, routine)
call nc_put_variable(ncid, 'hybi', grid_data%hybi%vals, routine)
call nc_put_variable(ncid, 'P0',   grid_data%P0%vals,   routine)

call nc_sync(ncid)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
!> writes CAM's model date and time of day into file.  CAM uses
!> integer date values and interger time of day measured in seconds
!>
!> @param ncid         name of the file
!> @param model_time   the current time of the model state
!>

subroutine write_model_time(ncid, model_time)
integer,         intent(in) :: ncid
type(time_type), intent(in) :: model_time

integer :: iyear, imonth, iday, ihour, iminute, isecond
integer :: cam_date(1), cam_tod(1)

character(len=*), parameter :: routine = 'write_model_time'

if ( .not. module_initialized ) call static_init_model

call get_date(model_time, iyear, imonth, iday, ihour, iminute, isecond)

cam_date = iyear*10000 + imonth*100 + iday
cam_tod  = ihour*3600  + iminute*60 + isecond

! if the file doesn't already have a "date" variable, so we make one
if (.not. nc_variable_exists(ncid, "date")) then
   call error_handler(E_MSG, routine,'"date" variable not found in file ', &
                      source, revision, revdate, text2='creating one')

   call nc_redef(ncid)
   call nc_define_integer_variable(ncid, 'date', (/ 'time' /), routine)
   call nc_enddef(ncid)
   call nc_put_variable(ncid, 'date', cam_date, routine)
endif

! if the file doesn't already have a "datesec" variable, so we make one
if (.not. nc_variable_exists(ncid, "datesec")) then

   call error_handler(E_MSG, routine,'"datesec" variable not found in file ', &
                      source, revision, revdate, text2='creating one')

   call nc_redef(ncid)
   call nc_define_integer_variable(ncid, 'datesec', (/ 'time' /), routine)
   call nc_enddef(ncid)
   call nc_put_variable(ncid, 'datesec', cam_tod,  routine)
endif

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> Read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid
integer :: timesize
integer :: cam_date, cam_tod
integer :: iyear, imonth, iday, ihour, imin, isec, rem

character(len=*), parameter :: routine = 'read_model_time'

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) trim(filename), ' does not exist.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif

ncid = nc_open_readonly(filename, routine)

! CAM initial files have two variables of length 
! 'time' (the unlimited dimension): date, datesec

call nc_get_variable(ncid, 'date',    cam_date)
call nc_get_variable(ncid, 'datesec', cam_tod)

! 'date' is YYYYMMDD 
! 'cam_tod' is seconds of current day
iyear  = cam_date / 10000
rem    = cam_date - iyear*10000
imonth = rem / 100
iday   = rem - imonth*100

ihour  = cam_tod / 3600
rem    = cam_tod - ihour*3600
imin   = rem / 60
isec   = rem - imin*60

! some cam files are from before the start of the gregorian calendar.
! since these are 'arbitrary' years, just change the offset.
if (iyear < 1601) then
   write(string1,*)' '
   write(string2,*)'WARNING - ',trim(filename),' changing year from ', &
                   iyear,'to',iyear+1601

   call error_handler(E_MSG, routine, string1, source, revision, &
                      revdate, text2=string2,text3='to make it a valid Gregorian date.')

   write(string1,*)' '
   call error_handler(E_MSG, routine, string1, source, revision)
   iyear = iyear + 1601
endif

read_model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

call nc_close(ncid, routine)

end function read_model_time

!--------------------------------------------------------------------
!> if the namelist is set to not use this custom routine, the default
!> dart routine will add 'pert_amp' of noise to every field in the state
!> to generate an ensemble from a single member.  if it is set to true
!> this routine will be called.  the pert_amp will be ignored, and the
!> given list of quantities will be perturbed by the given amplitude
!> (which can be different for each field) to generate an ensemble.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle 
integer,   intent(in) :: ens_size
real(r8),  intent(in) :: pert_amp   ! ignored in this version
logical,  intent(out) :: interf_provided

type(random_seq_type) :: seq

integer :: iloc, jloc, vloc, var_id, myqty
integer :: max_qtys, j

integer(i8) :: i, items
integer(i8), allocatable :: my_vars(:)

logical,  allocatable :: do_these_qtys(:)
real(r8), allocatable :: perturb_by(:)

character(len=*), parameter :: routine = 'pert_model_copies:'

! set by namelist to select using the default routine or the code here
if (custom_routine_to_generate_ensemble) then
   interf_provided = .true.
else
   interf_provided = .false.
   return
endif

! make sure each task is using a different random sequence
call init_random_seq(seq, my_task_id())

max_qtys = get_num_quantities()
allocate(do_these_qtys(max_qtys), perturb_by(max_qtys))

do_these_qtys(:) = .false.
perturb_by(:)    = 0.0_r8

do i=1, MAX_PERT
   if (fields_to_perturb(i) == '') exit
 
   myqty = get_index_for_quantity(fields_to_perturb(i))
   if (myqty < 0) then
      string1 = 'unrecognized quantity name in "fields_to_perturb" list: ' // &
                trim(fields_to_perturb(i))
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   do_these_qtys(i) = .true.
   perturb_by(i)    = perturbation_amplitude(i)
enddo

! get the global index numbers of the part of the state that 
! we have in this task.
items = get_my_num_vars(state_ens_handle)
allocate(my_vars(items))
call get_my_vars(state_ens_handle, my_vars)

do i=1, items
   call get_model_variable_indices(my_vars(i), iloc, jloc, vloc, kind_index=myqty)

   ! if myqty is in the namelist, perturb it.  otherwise cycle
   if (.not. do_these_qtys(myqty)) cycle
  
   do j=1, ens_size
      state_ens_handle%copies(j, i) = random_gaussian(seq, state_ens_handle%copies(j, i), perturb_by(i))
   enddo

enddo

deallocate(my_vars)
deallocate(do_these_qtys, perturb_by)

end subroutine pert_model_copies


!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!>
!>@param variable_array  the list of variables and kinds from model_mod_nml
!>@param nfields         the number of variable/KIND pairs specified

subroutine set_cam_variable_info( variable_array, nfields )

character(len=*), intent(in)  :: variable_array(:)
integer,          intent(out) :: nfields

character(len=*), parameter :: routine = 'set_cam_variable_info:'

integer :: i

!>@todo FIXME what should these be?  hardcode to 128 just for a hack.
character(len=128) :: varname    ! column 1, NetCDF variable name
character(len=128) :: dartstr    ! column 2, DART KIND
character(len=128) :: minvalstr  ! column 3, Clamp min val
character(len=128) :: maxvalstr  ! column 4, Clamp max val
character(len=128) :: updatestr  ! column 5, Update output or not

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  :: update_list(MAX_STATE_VARIABLES)   = .FALSE.
integer  ::   kind_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::  clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8


nfields = 0
ParseVariables : do i = 1, MAX_STATE_VARIABLES

   varname   = variable_array(num_state_table_columns*i-4)
   dartstr   = variable_array(num_state_table_columns*i-3)
   minvalstr = variable_array(num_state_table_columns*i-2)
   maxvalstr = variable_array(num_state_table_columns*i-1)
   updatestr = variable_array(num_state_table_columns*i  )

   if ( varname == ' ' .and. dartstr == ' ' ) exit ParseVariables ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml:model "state_variables" not fully specified'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(3A)') 'there is no obs_kind <', trim(dartstr), '> in obs_kind_mod.f90'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   call to_upper(minvalstr)
   call to_upper(maxvalstr)
   call to_upper(updatestr)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstr)
   clamp_vals(i,2) = string_to_real(maxvalstr)
   update_list( i) = string_to_logical(updatestr, 'UPDATE')

   nfields = nfields + 1

enddo ParseVariables

if (nfields == MAX_STATE_VARIABLES) then
   write(string1,'(2A)') 'WARNING: There is a possibility you need to increase ', &
                         'MAX_STATE_VARIABLES in the global variables in model_mod.f90'

   write(string2,'(A,i4,A)') 'WARNING: you have specified at least ', nfields, &
                             ' perhaps more'

   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

domain_id = add_domain(cam_template_filename, nfields, var_names, kind_list, &
                       clamp_vals, update_list )

!>@todo JPH we may need to call map_qtys. where do we call this?
call fill_cam_stagger_info(grid_stagger)

if (debug_level > 2) call state_structure_info(domain_id)

end subroutine set_cam_variable_info


!-----------------------------------------------------------------------
!>
!> Fill the qty_stagger array to tell what type of stagger each variable 
!> has. This will be useful for interpolating observations.
!> This currently doesn't support both slon/slat stagger - but cam-fv 
!> doesn't have any fields like that.
!>

subroutine fill_cam_stagger_info(stagger)
type(cam_stagger), intent(inout) :: stagger

integer :: ivar, jdim, qty_index

allocate(stagger%qty_stagger(0:get_num_quantities()))

stagger%qty_stagger = STAGGER_NONE

do ivar = 1, get_num_variables(domain_id)
   do jdim = 1, get_num_dims(domain_id, ivar)

      if (get_dim_name(domain_id, ivar, jdim) == 'slat') then
         qty_index = get_kind_index(domain_id, ivar) 
         stagger%qty_stagger(qty_index) = STAGGER_U
      endif

      if (get_dim_name(domain_id, ivar, jdim) == 'slon') then
         qty_index = get_kind_index(domain_id, ivar)
         stagger%qty_stagger(qty_index) = STAGGER_V
      endif

      if (get_dim_name(domain_id, ivar, jdim) == 'ilev') then
         qty_index = get_kind_index(domain_id, ivar)
         stagger%qty_stagger(qty_index) = STAGGER_W
      endif

   enddo
enddo

end subroutine fill_cam_stagger_info


!-----------------------------------------------------------------------
!> Read in the grid information from the given CAM restart file.
!> Note that none of the data will be used from this file; just the
!> grid size and locations.  Also read in the elevation information
!> from the "PHIS' file.

subroutine read_grid_info(grid_file, grid)
character(len=*), intent(in)  :: grid_file
type(cam_grid),   intent(out) :: grid

! Get the grid info plus additional non-state arrays
call get_cam_grid(grid_file, grid)

! This non-state variable is used to compute surface elevation.
call read_cam_phis_array(cam_phis_filename)

! Set up the interpolation structures for later 
call setup_interpolation(grid)

end subroutine read_grid_info


!-----------------------------------------------------------------------
!> Read the data from the various cam grid arrays 
!>
!>@todo FIXME not all of these are used.  can we either
!> not read them in, or make them optional?  this does affect
!> what we can write out in the diagnostic file.  if we have
!> to have them in the diag files then we have to read them all
!> even if we never use them.  both ilev and gw currently fall
!> into this category.
!> 

subroutine get_cam_grid(grid_file, grid)
character(len=*), intent(in)  :: grid_file
type(cam_grid), intent(out) :: grid

character(len=*), parameter :: routine = 'get_cam_grid:'

integer :: ncid

! put this in a subroutine that deals with the grid
ncid = nc_open_readonly(grid_file, routine)

call fill_cam_1d_array(ncid, 'lon',  grid%lon)
call fill_cam_1d_array(ncid, 'lat',  grid%lat)
call fill_cam_1d_array(ncid, 'lev',  grid%lev)
call fill_cam_1d_array(ncid, 'ilev', grid%ilev) ! for staggered vertical grid
call fill_cam_1d_array(ncid, 'slon', grid%slon)
call fill_cam_1d_array(ncid, 'slat', grid%slat)
call fill_cam_1d_array(ncid, 'gw',   grid%gw)   ! gauss weights
call fill_cam_1d_array(ncid, 'hyai', grid%hyai)
call fill_cam_1d_array(ncid, 'hybi', grid%hybi)
call fill_cam_1d_array(ncid, 'hyam', grid%hyam)
call fill_cam_1d_array(ncid, 'hybm', grid%hybm)

! P0 is a scalar with no dimensionality
call fill_cam_0d_array(ncid, 'P0',   grid%P0) 

call nc_close(ncid, routine)

end subroutine get_cam_grid


!-----------------------------------------------------------------------
!>
!> allocate space for a scalar variable and read values into the grid_array
!>   


subroutine fill_cam_1d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

integer :: i, per_line
!>@todo need to check that this exists
call nc_get_variable_size(ncid, varname, grid_array%nsize)
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals)

!>@todo FIXME this should be an array_dump() routine
!> in a utilities routine somewhere. 
if (debug_level > 10) then
   !call array_dump2(varname, grid_array%vals(:,:), nper_linei, nsize)
   !subroutine
   per_line = 5
   print*, 'variable name ', trim(varname)
   do i=1, grid_array%nsize, per_line
      print*,  grid_array%vals(i:min(grid_array%nsize,i+per_line-1))
   enddo
endif

end subroutine fill_cam_1d_array


!-----------------------------------------------------------------------
!>
!> allocate space for a scalar variable and read values into the grid_array
!>   


subroutine fill_cam_0d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

grid_array%nsize = 1
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals)

if (debug_level > 10) then
   print*, 'variable name ', trim(varname), grid_array%vals
endif

end subroutine fill_cam_0d_array

!-----------------------------------------------------------------------
!>
!> free space in the various grid arrays
!> 

subroutine free_cam_grid(grid)

type(cam_grid), intent(inout) :: grid

call free_cam_1d_array(grid%lon)
call free_cam_1d_array(grid%lat)
call free_cam_1d_array(grid%lev)
call free_cam_1d_array(grid%ilev) 
call free_cam_1d_array(grid%slon)
call free_cam_1d_array(grid%slat)
call free_cam_1d_array(grid%gw)
call free_cam_1d_array(grid%hyai)
call free_cam_1d_array(grid%hybi)
call free_cam_1d_array(grid%hyam)
call free_cam_1d_array(grid%hybm)

call free_cam_1d_array(grid%P0) 

deallocate(phis)

end subroutine free_cam_grid


!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine free_cam_1d_array(grid_array)
type(cam_1d_array), intent(inout) :: grid_array

deallocate(grid_array%vals)
grid_array%nsize = -1

end subroutine free_cam_1d_array

!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine set_vert_localization(typename, vcoord)
character(len=*), intent(in)  :: typename
integer,          intent(out) :: vcoord

character(len=*), parameter :: routine = 'set_vert_localization'

character(len=32) :: ucasename

ucasename = typename
call to_upper(ucasename)

select case (ucasename)
  case ("PRESSURE")
     vcoord = VERTISPRESSURE
  case ("HEIGHT")
     vcoord = VERTISHEIGHT
  case ("SCALEHEIGHT", "SCALE_HEIGHT")
     vcoord = VERTISSCALEHEIGHT
  case ("LEVEL", "MODEL_LEVEL")
     vcoord = VERTISLEVEL
  case default
     write(string1,*)'unrecognized vertical localization coordinate type: '//trim(typename)
     write(string2,*)'valid values are: PRESSURE, HEIGHT, SCALEHEIGHT, LEVEL'
     call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
end select
 
end subroutine set_vert_localization

!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine setup_interpolation(grid)
type(cam_grid), intent(in) :: grid

!>@todo FIXME the cam fv grid is really evenly spaced in lat and lon,
!>even though they provide full lon() and lat() arrays.  the deltas
!>between each pair would be faster

! mass points at cell centers
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%lon%nsize, grid%lat%nsize, QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_nonstaggered)
call set_quad_coords(interp_nonstaggered, grid%lon%vals, grid%lat%vals)

! U stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%lon%nsize, grid%slat%nsize, QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_staggered)
call set_quad_coords(interp_u_staggered, grid%lon%vals, grid%slat%vals)

! V stagger
call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid%slon%nsize, grid%lat%nsize, QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_staggered)
call set_quad_coords(interp_v_staggered, grid%slon%vals, grid%lat%vals)

end subroutine setup_interpolation

!-----------------------------------------------------------------------
!>
!> 

subroutine read_cam_phis_array(phis_filename)
character(len=*),   intent(in)    :: phis_filename

character(len=*), parameter :: routine = 'read_cam_phis_array'

integer :: ncid, nsize(2)

ncid = nc_open_readonly(phis_filename, routine)

call nc_get_variable_size(ncid, 'PHIS', nsize(:))
allocate( phis(nsize(1), nsize(2)) )

call nc_get_variable(ncid, 'PHIS', phis)

call nc_close(ncid, routine)

end subroutine read_cam_phis_array

!-----------------------------------------------------------------------

subroutine build_heights(n_levels,p_surf,h_surf,vertual_temp,height_midpts,height_interf,variable_r)

integer,  intent(in)  :: n_levels                            ! Number of vertical levels
real(r8), intent(in)  :: p_surf                              ! Surface pressure (pascals)
real(r8), intent(in)  :: h_surf                              ! Surface height (m)
real(r8), intent(in)  :: vertual_temp( n_levels)             ! Vertual Temperature
real(r8), intent(out) :: height_midpts(n_levels)             ! Geopotential height at midpoints, top to bottom
real(r8), intent(out), optional :: height_interf(n_levels+1) ! Geopotential height at interfaces, top to bottom
real(r8), intent(in),  optional :: variable_r(   n_levels)   ! Dry air gas constant, if varies, top to bottom

! Local variables
!>@todo have a model constants module?
real(r8), parameter :: const_r = 287.04_r8    ! Different than model_heights ! dry air gas constant.
real(r8), parameter :: g0 = 9.80616_r8        ! Different than model_heights:gph2gmh:G !
real(r8) :: dist_from_interf(n_levels)         ! pressure profile
real(r8) :: pterm(n_levels)
real(r8) :: log_points_mid(n_levels+1)

integer  :: i,k,l
real(r8) :: temp, hybrid_base_pressure

! an array now: real(r8), parameter :: rbyg=r/g0
real(r8) :: const_r_g0(n_levels)      ! rbyg=r/g0
real(r8) :: midpts(n_levels)          ! midpoints  in column
real(r8) :: interf(n_levels)          ! interfaces in column
real(r8) :: log_p_midpoints(n_levels) ! log of pressure at layer midpoints
real(r8) :: log_p_interf(n_levels+1)  ! log of pressure at layer interfaces
real(r8) :: geo_hgt_intf(n_levels+1)  ! Geopotential height at interfaces, top to bottom

!state of dry air dry air Air pressure
!    p=rho*R*T 
! where rho   is the mass density 
!       alpha is specific volume
!       R      gas constant and 
!       T     temperature
! gravity force must exactly be exactly balanced by  by the vertical component of the pressure gradiant force
!    dp/dz = -rho*g    
! hydrostatic balance
!    p(z) = /int_z^/infinity rho*g
! finite difference forms of the Dyn operator, including semi implicit time integration are :

! normally we only need the midpoint heights, but if its useful to have
! the heights on the interface levels, here it is.  


! Invert vertical loop
! Compute top only if top interface pressure is nonzero.
! Implemented by setting loop limit klim
!
! hyba, hybb are bottom to top, starting at ground.
! log_points_mid(i,k) is the mid-point pressure of layer k.
! SHEA MODIFICATION

if (present(variable_r)) then
   do i=1, n_levels
      const_r_g0(i) = variable_r(i) / g0
   enddo
else
      const_r_g0(:) = const_r / g0
endif

call cam_p_col_midpts(p_surf, n_levels,   midpts)
call cam_p_col_intfcs(p_surf, n_levels+1, interf)

if (debug_level > 200) then
   do i=1, n_levels
     print *, 'midpts, interf: ', i, log_p_midpoints(i), log_p_interf(i)
   enddo
   print *, 'interf: ', n_levels+1, log_p_interf(n_levels+1)
endif

log_p_midpoints = log(midpts)
log_p_interf    = log(interf)
hybrid_base_pressure = grid_data%P0%vals(1)

do k = n_levels+1, 1, -1
   i = n_levels-k+2
   ! hybrid_As coord coeffs for P0 reference pressure term in plevs_cam
   ! hybrid_Bs coord coeffs for surf pressure term in plevs_cam (in same format as hybrid_As)
   temp = hybrid_base_pressure*grid_data%hyam%vals(k) + p_surf*grid_data%hybm%vals(k)
   if (temp > 0.0_r8) THEN 
      log_points_mid(k) = log(temp)
   else 
      log_points_mid(k) = 0.0_r8
   endif
enddo

! Initialize height_midpts to sum of ground height and thickness of
! top half-layer  (i.e. (phi)sfc in equation 1.14)
!
!        (height_midpts(1)=top  ->  height_midpts(n_levels)=bottom
! 
!        Eq 3.a.109.2  where l=K,k<K  h(k,l) = 1/2 * ln [  p(k+1) / p(k) ]

do k = 2,n_levels - 1
   pterm(k) = const_r_g0(k)*vertual_temp(k)*0.5_r8*(log_points_mid(k+1)-log_points_mid(k-1))
enddo

! Initialize height_midpts to sum of ground height and thickness of top half layer
! this is NOT adding the thickness of the 'top' half layer.
!    it's adding the thickness of the half layer at level K,
do k = 1,n_levels - 1
   height_midpts(k) = h_surf + const_r_g0(k)*vertual_temp(k)*0.5_r8*(log_points_mid(k+1)-log_points_mid(k))

enddo

! Add thickness of the remaining full layers
!  (i.e., integrate from ground to highest layer interface)
! 
! Eqs 1.14 & 3.a.109.3 where l>K, k<K
!                          h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]

height_midpts(n_levels) = h_surf + const_r_g0(n_levels)*vertual_temp(n_levels)* &
                          (log_p_interf(n_levels)-log_points_mid(n_levels))

! THIS is adding the half layer at the BOTTOM.
!
! Add thickness of the remaining full layers
!  (i.e., integrate from ground to highest layer interface)
! 
! Eqs 1.14 & 3.a.109.3 where l>K, k<K
!                          h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]
do k = 1,n_levels - 1
    height_midpts(k) = height_midpts(k) + const_r_g0(k)*vertual_temp(k)*log_p_interf(k)-0.5_r8* &
                       (log_points_mid(n_levels-1)+log_points_mid(n_levels))
enddo


! Add thickness of the remaining full layers
! (i.e., integrate from ground to highest layer interface)
!
!       Eqs 1.14 & 3.a.109.3 where l>K, k<K
!                                h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]
!
! 3.a.109.3

do k = 1,n_levels - 2
   do l = k+1, n_levels-1
      height_midpts(k) = height_midpts(k) + pterm(l)
   enddo
enddo

end subroutine build_heights

!-----------------------------------------------------------------------
!>  Convert a 2d array of geopotential altitudes to mean sea level altitudes.

subroutine gph2gmh(h, lat)

real(r8), intent(inout) :: h(:,:)    ! geopotential altitude in m
real(r8), intent(in)    :: lat       ! latitude in degrees.

real(r8), parameter ::  be = 6356751.6_r8        ! min earth radius, m
real(r8), parameter ::  ae = 6378136.3_r8        ! max earth radius, m
real(r8), parameter ::  pi = 3.14159265358979_r8
real(r8), parameter ::  G = 9.80665_r8 ! WMO reference g value, m/s**2, at 45.542N(S)

real(r8) :: g0
real(r8) :: r0
real(r8) :: latr

integer :: i, j

latr = lat * DEG2RAD  ! convert to radians
call compute_surface_gravity(latr, g0)

! compute local earth's radius using ellipse equation

r0 = sqrt( ae**2 * cos(latr)**2 + be**2 * sin(latr)**2)

! Compute altitude above sea level
do j=1, size(h, 2)
   do i=1, size(h, 1)
      h(i,j) = (r0 * h(i,j)) / (((g0*r0)/G) - h(i,j))
   enddo
enddo

end subroutine gph2gmh

!-----------------------------------------------------------------------
!> This subroutine computes the Earth's gravity at any latitude.
!> The model assumes the Earth is an oblate spheriod rotating at 
!> the Earth's spin rate.  The model was taken from 
!> "Geophysical Geodesy, Kurt Lambeck, 1988".
!>
!>  input:    xlat, latitude in radians
!>  output:   galt, gravity at the given lat, m/sec**2
!>
!> taken from code from author Bill Schreiner, 5/95
!>
!>

subroutine compute_surface_gravity(xlat,galt)

real(r8), intent(in)  :: xlat
real(r8), intent(out) :: galt

real(r8),parameter :: xmu = 398600.4415_r8         ! km^3/s^2
real(r8),parameter :: ae  = 6378.1363_r8           ! km
real(r8),parameter :: f   = 1.0_r8/298.2564_r8
real(r8),parameter :: xm  = 0.003468_r8            !
real(r8),parameter :: f2  = 5.3481622134089e-03_r8 ! f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
real(r8),parameter :: f4  = 2.3448248012911e-05_r8 ! f4 = -f**2* 0.50 + 5.0* 0.50*f*xm

! gravity at the equator, km/s2
real(r8), parameter :: ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0_r8/14.0_r8*xm*f)


! compute gravity at any latitude, km/s2
galt = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)

! convert to meters/s2
galt = galt*1000.0_r8

end subroutine compute_surface_gravity

!-----------------------------------------------------------------------
!> This subroutine computes converts vertical state
!>
!>  in:    ens_handle  - mean ensemble handle
!>  in:    num         - number of locations
!>  inout: locs(:)     - locations
!>  in:    loc_qtys(:) - location quantities
!>  in:    loc_indx(:) - location index
!>  in:    which_vert  - vertical location to convert
!>  out:   istatus     - return status 0 is a successful conversion
!>

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

character(len=*), parameter :: routine = 'convert_vertical_state'

integer :: current_vert_type, ens_size, i

ens_size = 1

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if ( current_vert_type == which_vert ) cycle

   select case (which_vert)
      case (VERTISPRESSURE)
         call convert_to_pressure(    ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISHEIGHT)
         call convert_to_height(      ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISLEVEL)
         call convert_to_level(                   ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case (VERTISSCALEHEIGHT)
         call convert_to_scaleheight( ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i) )
      case default
         write(string1,*)'unable to convert vertical state "', which_vert, '"'
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end select
enddo

istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

subroutine  convert_to_pressure(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

real(r8) :: pressure

!>@todo FIXME move up to convert_vertical_state ?
call state_level_to_pressure(ens_handle, ens_size, location_indx, qty, pressure)

call set_vertical(location, pressure, VERTISPRESSURE)

end subroutine  convert_to_pressure

!--------------------------------------------------------------------

subroutine state_level_to_pressure(ens_handle, ens_size, location_indx, qty, pressure)
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer(i8),         intent(in)  :: location_indx
integer,             intent(in)  :: qty
real(r8),            intent(out) :: pressure

integer  :: iloc, jloc, vloc
integer  :: my_status(ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)

call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                         qty, pressure_array, my_status)

pressure = pressure_array(vloc)

end subroutine state_level_to_pressure


!--------------------------------------------------------------------

subroutine convert_to_height(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

!>@todo FIXME move up to convert_vertical_state ?
call state_level_to_height(ens_handle, ens_size, location, location_indx, qty)

end subroutine  convert_to_height

!--------------------------------------------------------------------

subroutine state_level_to_height(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty


integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)

! build a height column and a pressure column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_height_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       qty, height_array, my_status) 

!>@todo FIXME this can only be used if ensemble size is 1
call set_vertical(location, height_array(vloc,1), VERTISHEIGHT)

end subroutine state_level_to_height

!--------------------------------------------------------------------

subroutine  convert_to_scaleheight(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

call state_level_to_scaleheight(ens_handle, ens_size, location, location_indx, qty)

end subroutine  convert_to_scaleheight

!--------------------------------------------------------------------

subroutine convert_to_level(ens_size, location, location_indx, qty)
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

call state_level_to_level(ens_size, location, location_indx, qty)

end subroutine convert_to_level

!--------------------------------------------------------------------

subroutine state_level_to_level(ens_size, location, location_indx, qty)
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc
integer  :: my_status(ens_size)

!>@todo FIXME qty is currently unused.  if we need it, its here.
!>if we really don't need it, we can remove it.  all the other
!>corresponding routines like this use it.

call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call set_vertical(location, real(vloc,r8), VERTISLEVEL)

end subroutine state_level_to_level



!--------------------------------------------------------------------

subroutine state_level_to_scaleheight(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty

integer  :: iloc, jloc, vloc, level_one, status1, my_status(ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)
real(r8) :: surface_pressure(1), scaleheight

level_one = 1

! build a height column and a scaleheight column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

! get the surface pressure from the ens_handle
call get_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                         iloc, jloc, level_one, surface_pressure, status1)
   
call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                         qty, pressure_array, my_status)
   
scaleheight = log(pressure_array(vloc)/surface_pressure(1))

call set_vertical(location, scaleheight, VERTISSCALEHEIGHT)

end subroutine state_level_to_scaleheight

!--------------------------------------------------------------------
!> this only works for ens_size = 1

subroutine obs_height_to_pressure(ens_handle, ens_size, location, location_indx, qty)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx
integer,             intent(in)    :: qty


integer  :: iloc, jloc, vloc, my_status(ens_size), i, j
real(r8) :: height_array(  grid_data%lev%nsize, ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize, ens_size)

!>@todo FIXME this isn't done but is not called yet.
!>it needs to search for the actual height in the first column
!>and then use that index into the second.

! build a height column and a pressure column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_height_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       qty, height_array, my_status) 

print*, iloc, jloc
do j = 1, ens_size
   print*, j
   print*, height_array(:,j)
enddo

call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       qty, pressure_array, my_status)

call set_vertical(location, pressure_array(vloc,1), VERTISPRESSURE)

end subroutine obs_height_to_pressure


!-----------------------------------------------------------------------
!> Compute the pressure values at midpoint levels
!>
!> this version does all ensemble members at once.

subroutine cam_pressure_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, pressure_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
integer,             intent(in)  :: qty
real(r8),            intent(out) :: pressure_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: level_one, status1, imember
real(r8)    :: surface_pressure(ens_size)

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                         lon_index, lat_index, level_one, qty, surface_pressure, status1)

call build_cam_pressure_column(ens_size, surface_pressure, grid_data%lev%nsize, &
                               pressure_array)
my_status(:) = 0

end subroutine cam_pressure_levels

!--------------------------------------------------------------------

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, my_status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer,             intent(in)    :: loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: my_status(:)

character(len=*), parameter :: routine = 'convert_vertical_obs'

integer :: current_vert_type, i, ens_size


ens_size = 1

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if ( current_vert_type == which_vert ) cycle

   select case (which_vert)
      case (VERTISPRESSURE)
         call convert_obs_to_pressure(ens_handle, ens_size, locs(i))
      case (VERTISHEIGHT)
         call convert_obs_to_height(ens_handle, ens_size, locs(i))
      case (VERTISLEVEL)
         call convert_obs_to_level(ens_handle, ens_size, locs(i))
      case (VERTISSCALEHEIGHT)
         call convert_obs_to_scaleheight(ens_handle, ens_size, locs(i))
      case default
         write(string1,*)'unable to convert vertical obs "', which_vert, '"'
         call error_handler(E_MSG,routine,string1,source,revision,revdate)
   end select
enddo

my_status(:) = 0

end subroutine convert_vertical_obs

!--------------------------------------------------------------------

subroutine convert_obs_to_pressure(ens_handle, ens_size, location)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

call obs_vertical_to_pressure(ens_handle, ens_size, location)

end subroutine  convert_obs_to_pressure

!--------------------------------------------------------------------

subroutine obs_vertical_to_pressure(ens_handle, ens_size, location)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

integer  :: my_status(ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)

call model_interpolate(ens_handle, ens_size, location, QTY_PRESSURE, pressure_array(:), my_status(:))

call set_vertical(location, pressure_array(1), VERTISPRESSURE)

end subroutine obs_vertical_to_pressure

!--------------------------------------------------------------------

subroutine convert_obs_to_height(ens_handle, ens_size, location)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

call obs_vertical_to_height(ens_handle, ens_size, location)

end subroutine  convert_obs_to_height

!--------------------------------------------------------------------

subroutine obs_vertical_to_height(ens_handle, ens_size, location)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

integer  :: my_status(    ens_size )
real(r8) :: height_array( ens_size )

call model_interpolate(ens_handle, ens_size, location, &
                       QTY_GEOMETRIC_HEIGHT, height_array(:), my_status(:))

call set_vertical(location, height_array(1), VERTISHEIGHT)

end subroutine obs_vertical_to_height

!--------------------------------------------------------------------

subroutine convert_obs_to_level(ens_handle, ens_size, location)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

call obs_vertical_to_level(ens_handle, ens_size, location)

end subroutine  convert_obs_to_level

!--------------------------------------------------------------------

subroutine obs_vertical_to_level(ens_handle, ens_size, location)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

integer  :: my_status(ens_size)
real(r8) :: level_array(ens_size)

call model_interpolate(ens_handle, ens_size, location, QTY_VERTLEVEL, level_array(:), my_status(:))

call set_vertical(location, level_array(1), VERTISLEVEL)

end subroutine obs_vertical_to_level

!--------------------------------------------------------------------

subroutine convert_obs_to_scaleheight(ens_handle, ens_size, location)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

call obs_vertical_to_scaleheight(ens_handle, ens_size, location)

end subroutine  convert_obs_to_scaleheight

!--------------------------------------------------------------------

subroutine obs_vertical_to_scaleheight(ens_handle, ens_size, location)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location

integer  :: my_status(ens_size)
real(r8) :: pressure_array(ens_size), surface_pressure_array(ens_size)
real(r8) :: scaleheight

call model_interpolate(ens_handle, ens_size, location, QTY_PRESSURE,         pressure_array(:),         my_status(:))
call model_interpolate(ens_handle, ens_size, location, QTY_SURFACE_PRESSURE, surface_pressure_array(:), my_status(:))

scaleheight = log(pressure_array(1)/surface_pressure_array(1))

call set_vertical(location, scaleheight, VERTISSCALEHEIGHT)

end subroutine obs_vertical_to_scaleheight

!--------------------------------------------------------------------

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
