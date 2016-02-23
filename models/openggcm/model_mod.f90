! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between the openggcm ocean model and DART.

! Modules that are absolutely required for use are listed
use        types_mod,    only : r4, r8, i4, i8, SECPERDAY, MISSING_R8, rad2deg, PI
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date,                           &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)
use     location_mod, only : location_type, get_dist, get_close_maxdist_init,  &
                             get_close_obs_init, set_location,                 &
                             VERTISUNDEF, VERTISHEIGHT, get_location,          &
                             vert_is_height, vert_is_level, vert_is_surface,   &
                             loc_get_close_obs => get_close_obs, get_close_type
use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text
use     obs_kind_mod, only : KIND_ELECTRON_DENSITY, KIND_ELECTRIC_POTENTIAL,   &
                             get_raw_obs_kind_index, get_raw_obs_kind_name,    &
                             paramname_length 
use mpi_utilities_mod, only: my_task_id, task_count
use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian
use ensemble_manager_mod,  only : ensemble_type, map_pe_to_task, get_copy_owner_index, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_model_variable_indices, &
                                  get_num_variables, get_index_start, &
                                  get_num_dims, get_domain_size

use dart_time_io_mod,      only : write_model_time

use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,                &
          adv_1step,                     &
          get_state_meta_data,           &
          model_interpolate,             &
          get_model_time_step,           &
          static_init_model,             &
          end_model,                     &
          init_time,                     &
          init_conditions,               &
          nc_write_model_atts,           &
          nc_write_model_vars,           &
          pert_model_copies,             &
          get_close_maxdist_init,        &
          get_close_obs_init,            &
          get_close_obs,                 &
          query_vert_localization_coord, &
          vert_convert,                  &
          construct_file_name_in,        &
          read_model_time,               &
          write_model_time



! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! message strings
character(len=512) :: string1
character(len=512) :: string2
character(len=512) :: msgstring

logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10 
integer, parameter :: num_state_table_columns = 3
! NOTE: may need to increase character length if netcdf variables are
! larger than paramname_length = 32.
character(len=paramname_length) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
logical :: update_var_list( max_state_variables )

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_KIND_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

! things which can/should be in the model_nml
character(len=NF90_MAX_NAME) :: openggcm_template
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
character(len=paramname_length) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 0   ! turn up for more and more debug messages

namelist /model_nml/  &
   openggcm_template,           &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   model_state_variables,       &
   debug

!------------------------------------------------------------------
! fields in the state vector
!------------------------------------------------------------------

! Number of fields in the state vector
integer :: nfields

! Grid parameters - the values will be read from the
! openggcm netcdf restart file

! the size of the grid
integer :: nlon, nlat, nvert

real(r8), allocatable :: grid_longitude(:), grid_latitude(:)
real(r8), allocatable :: levels(:)

type(time_type) :: model_time, model_timestep

integer(i8) :: model_size    ! the state vector length


!------------------------------------------------


! global domain id to be used by routines in state_structure_mod
integer :: domain_id

contains

!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine static_init_model()

! Called to do one time initialization of the model. In this case,
! it reads in the grid information.

integer :: iunit, io, i, ncid
integer :: ss, dd

! The Plan:
!
!   read in the grid sizes from the horiz grid file and the vert grid file
!   horiz is netcdf, vert is ascii
!  
!   allocate space, and read in actual grid values
!
!   figure out model timestep.  FIXME: from where?
!
!   Compute the model size.
!
!   set the index numbers where the field types change
!

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)


! Set the time step ... causes openggcm namelists to be read.
! Ensures model_timestep is multiple of 'ocean_dynamics_timestep'

model_timestep = set_time(0,0)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(msgstring,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)


! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

!> @TODO FIXME

ncid = get_grid_template_fileid(openggcm_template)

call get_grid_sizes(ncid, nlon, nlat, nvert)

! Allocate space for grid variables. 
allocate(grid_longitude(nlon), grid_latitude(nlat))
allocate(levels(nvert))

call read_horiz_grid(ncid)
call read_vert_levels(ncid)

! verify that the model_state_variables namelist was filled in correctly.  
! returns variable_table which has variable names, kinds and update strings.
call verify_state_variables(model_state_variables, nfields, variable_table, state_kinds_list, update_var_list)

!> @todo  - need input filename, hardcode for now to openggcm.nc
domain_id = add_domain(openggcm_template, nfields, &
                       var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                       update_list = update_var_list(1:nfields))

model_size = get_domain_size(domain_id)
if (do_output()) write(*,*) 'model_size = ', model_size


end subroutine static_init_model

!------------------------------------------------------------

subroutine init_conditions(x)
 real(r8), intent(out) :: x(:)

! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

character(len=128) :: msgstring2, msgstring3

msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = 'use openggcm_to_dart to generate an initial state'
call error_handler(E_ERR,'init_conditions', &
                  'ERROR!!  openggcm model has no built-in default state', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

subroutine adv_1step(x, time)
 real(r8),        intent(inout) :: x(:)
 type(time_type), intent(in)    :: time

! If the model could be called as a subroutine, does a single
! timestep advance.  openggcm cannot be called this way, so fatal error
! if this routine is called.

call error_handler(E_ERR,'adv_1step', &
                  'openggcm model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

end subroutine adv_1step

!------------------------------------------------------------------

function get_model_size()
 integer(i8) :: get_model_size

! Returns the size of the model as an integer. Required for all
! applications.


if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------

subroutine init_time(time)
 type(time_type), intent(out) :: time

! Companion interface to init_conditions. Returns a time that is
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.

character(len=128) :: msgstring2, msgstring3

msgstring2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
msgstring3 = 'use openggcm_to_dart to generate an initial state which contains a timestamp'
call error_handler(E_ERR,'init_time', &
                  'ERROR!!  openggcm model has no built-in default time', &
                  source, revision, revdate, &
                  text2=msgstring2, text3=msgstring3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)

 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 integer,            intent(out) :: istatus(ens_size)
 real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values

! Model interpolate will interpolate any state variable
! the given location given a state vector. The type of the variable being
! interpolated is obs_type since normally this is used to find the expected
! value of an observation at some location. The interpolated value is 
! returned in interp_val and istatus is 0 for success.

! Local storage
real(r8)       :: loc_array(3), llon, llat, lheight
integer(i8)    :: base_offset
integer        :: ind
integer        :: hgt_bot, hgt_top
real(r8)       :: hgt_fract
integer        :: hstatus

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs(:) = MISSING_R8     ! the DART bad value flag
istatus(:) = 99                ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if (debug > 1) print *, 'requesting interpolation of ', obs_type, ' at ', llon, llat, lheight

if( vert_is_height(location) ) then
   ! Nothing to do 
elseif ( vert_is_surface(location) ) then
   ! Nothing to do 
elseif (vert_is_level(location)) then
   ! convert the levels index to an actual depth 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(levels)) ) then 
      istatus = 11
      return
   else
      lheight = levels(ind)
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 17
   return
endif

! The following kinds are either in the state vector (so you
! can simply interpolate to find the value) or they are a simple
! transformation of something in the state vector.

SELECT CASE (obs_type)
   CASE (KIND_ELECTRON_DENSITY)
      base_offset = get_index_start(domain_id, get_varid_from_kind(obs_type))

   CASE DEFAULT
      ! Not a legal type for interpolation, return istatus error
      istatus = 15
      return

END SELECT

if( vert_is_surface(location) ) then
   call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, 1, expected_obs, istatus)

   if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus
   return
endif

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, nvert, levels, hgt_bot, hgt_top, hgt_fract, hstatus)
if(hstatus /= 0) then
   istatus = 12
   return
endif

! do a 2d interpolation for the value at the bottom level, then again for
! the top level, then do a linear interpolation in the vertical to get the
! final value.  this sets both interp_val and istatus.
call do_interp(state_handle, ens_size, base_offset, hgt_bot, hgt_top, hgt_fract, &
               llon, llat, obs_type, expected_obs, istatus)
if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus

end subroutine model_interpolate

!------------------------------------------------------------------

!------------------------------------------------------------------
!> Is height ens_size? Should quad status be ens_size?
subroutine lon_lat_interpolate(state_handle, ens_size, offset, lon, lat, var_type, height, expected_obs, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer(i8),         intent(in)  :: offset ! Not sure if this is the best way to do this
real(r8),            intent(in)  :: lon, lat
integer,             intent(in)  :: var_type, height
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! Subroutine to interpolate to a lon lat location given the state vector
! for that level, x. This works just on one horizontal slice.
! NOTE: Using array sections to pass in the x array may be inefficient on some
! compiler/platform setups. Might want to pass in the entire array with a base
! offset value instead of the section if this is an issue.
! This routine works for either the dipole or a regular lat-lon grid.
! Successful interpolation returns istatus=0.

! Local storage
integer  :: lat_bot, lat_top, lon_bot, lon_top, num_inds, start_ind
real(r8) :: x_corners(4), y_corners(4)
real(r8) :: p(4,ens_size), xbot(ens_size), xtop(ens_size)
real(r8) :: lon_fract, lat_fract

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

!> @TODO FIXME this is where we need to find the lat/lon corners
! Getting corners for accurate interpolation
!#! call get_quad_corners(grid_longitude, lon_bot, lat_bot, x_corners)
!#! call get_quad_corners(grid_latitude, lon_bot, lat_bot, y_corners)

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > nlat) then
   istatus = 2
   return
endif

! Watch for wraparound in longitude
lon_top = lon_bot + 1
if(lon_top > nlon) lon_top = 1

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
p(1, :) = get_val(lon_bot, lat_bot, nlon, state_handle, offset, ens_size, var_type, height)
p(2, :) = get_val(lon_top, lat_bot, nlon, state_handle, offset, ens_size, var_type, height)

p(3, :) = get_val(lon_top, lat_top, nlon, state_handle, offset, ens_size, var_type, height)
p(4, :) = get_val(lon_bot, lat_top, nlon, state_handle, offset, ens_size, var_type, height)

! Rectangular biliear interpolation
xbot = p(1, :) + lon_fract * (p(2, :) - p(1, :))
xtop = p(4, :) + lon_fract * (p(3, :) - p(4, :))
! Now interpolate in latitude
expected_obs = xbot + lat_fract * (xtop - xbot)

end subroutine lon_lat_interpolate

!------------------------------------------------------------

function get_val(lon_index, lat_index, nlon, state_handle, offset, ens_size, var_type, height)
 integer,             intent(in)  :: lon_index, lat_index, nlon, var_type, height
 type(ensemble_type), intent(in)  :: state_handle
 integer(i8),         intent(in)  :: offset
 integer,             intent(in)  :: ens_size

 real(r8)    :: get_val(ens_size)
 integer(i8) :: state_index

! Returns the value from a single level array given the lat and lon indices

if ( .not. module_initialized ) call static_init_model

!> @TODO FIXME is this still right?
! state index must be 8byte integer
state_index = int(lat_index - 1,i8)*int(nlon,i8) + int(lon_index,i8) + int(offset-1,i8)

! The x above is only a horizontal slice, not the whole state.   HK WHY -1?
get_val = get_state(state_index, state_handle)

end function get_val

!------------------------------------------------------------

subroutine lon_bounds(lon, nlons, lon_array, bot, top, fract)
 real(r8),    intent(in) :: lon
 integer,     intent(in) :: nlons
 real(r8),    intent(in) :: lon_array(:, :)
 integer,    intent(out) :: bot, top
 real(r8),   intent(out) :: fract

! Given a longitude lon, the array of longitudes for grid boundaries, and the
! number of longitudes in the grid, returns the indices of the longitude
! below and above the location longitude and the fraction of the distance
! between. It is assumed that the longitude wraps around for a global grid. 
! Since longitude grids are going to be regularly spaced, this could be made more efficient.
! Algorithm fails for a silly grid that has only two longitudes separated by 180 degrees.

! Local storage
integer  :: i
real(r8) :: dist_bot, dist_top

if ( .not. module_initialized ) call static_init_model

! This is inefficient, someone could clean it up since longitudes are regularly spaced
! But note that they don't have to start at 0
do i = 2, nlons
   dist_bot = lon_dist(lon, lon_array(i - 1, 1))
   dist_top = lon_dist(lon, lon_array(i, 1))
   if(dist_bot <= 0 .and. dist_top > 0) then
      bot = i - 1
      top = i
      fract = abs(dist_bot) / (abs(dist_bot) + dist_top)
      return
   endif
enddo

! Falling off the end means it's in between; wraparound
bot = nlons
top = 1
dist_bot = lon_dist(lon, lon_array(bot, 1))
dist_top = lon_dist(lon, lon_array(top, 1)) 
fract = abs(dist_bot) / (abs(dist_bot) + dist_top)

end subroutine lon_bounds

!-------------------------------------------------------------

subroutine lat_bounds(lat, nlats, lat_array, bot, top, fract, istatus)
 real(r8),   intent(in) :: lat
 integer,    intent(in) :: nlats
 real(r8),   intent(in) :: lat_array(:, :)
 integer,   intent(out) :: bot, top
 real(r8),  intent(out) :: fract
 integer,   intent(out) :: istatus

! Given a latitude lat, the array of latitudes for grid boundaries, and the
! number of latitudes in the grid, returns the indices of the latitude
! below and above the location latitude and the fraction of the distance
! between. istatus is returned as 0 unless the location latitude is 
! south of the southernmost grid point (1 returned) or north of the 
! northernmost (2 returned). If one really had lots of polar obs would 
! want to worry about interpolating around poles.

! Local storage
integer :: i

if ( .not. module_initialized ) call static_init_model

! Success should return 0, failure a positive number.
istatus = 0

! Check for too far south or north
if(lat < lat_array(1, 1)) then
   istatus = 1
   return
else if(lat > lat_array(1, nlats)) then
   istatus = 2
   return
endif

! In the middle, search through
do i = 2, nlats
   if(lat <= lat_array(1, i)) then
      bot = i - 1
      top = i
      fract = (lat - lat_array(1, bot)) / (lat_array(1, top) - lat_array(1, bot))
      return
   endif
enddo

! Shouldn't get here. Might want to fail really hard through error handler
istatus = 40

end subroutine lat_bounds

!------------------------------------------------------------------

function lon_dist(lon1, lon2)
 real(r8), intent(in) :: lon1, lon2
 real(r8)             :: lon_dist

! Returns the smallest signed distance between lon1 and lon2 on the sphere
! If lon1 is less than 180 degrees east of lon2 the distance is negative
! If lon1 is less than 180 degrees west of lon2 the distance is positive

if ( .not. module_initialized ) call static_init_model

lon_dist = lon2 - lon1
if(lon_dist >= -180.0_r8 .and. lon_dist <= 180.0_r8) then
   return
else if(lon_dist < -180.0_r8) then
   lon_dist = lon_dist + 360.0_r8
else
   lon_dist = lon_dist - 360.0_r8
endif

end function lon_dist

!------------------------------------------------------------

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)
 real(r8),    intent(in) :: lheight
 integer,     intent(in) :: nheights
 real(r8),    intent(in) :: hgt_array(nheights)
 integer,    intent(out) :: bot, top
 real(r8),   intent(out) :: fract
 integer,    intent(out) :: istatus

! Local variables
integer   :: i

if ( .not. module_initialized ) call static_init_model

! Succesful istatus is 0
! Make any failure here return istatus in the 20s
istatus = 0

! The level array contains the depths of the center of the vertical grid boxes
! In this case (unlike how we handle the MIT depths), positive is really down.
! FIXME: in the MIT model, we're given box widths and we compute the centers,
! and we computed them with larger negative numbers being deeper.  Here,
! larger positive numbers are deeper.

! It is assumed that the top box is shallow and any observations shallower
! than the depth of this box's center are just given the value of the
! top box.
if(lheight <= hgt_array(1)) then
   top = 1
   bot = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   fract = 1.0_r8 
if (debug > 7) print *, 'above first level in height'
if (debug > 7) print *, 'hgt_array, top, bot, fract=', hgt_array(1), top, bot, fract
   return
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is shallower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i -1
      bot = i
      fract = (hgt_array(bot) - lheight) / (hgt_array(bot) - hgt_array(top))
if (debug > 7) print *, 'i, hgt_array, top, bot, fract=', i, hgt_array(i), top, bot, fract
      return
   endif
enddo

! Falling off the end means the location is lower than the deepest height
istatus = 20

end subroutine height_bounds

!------------------------------------------------------------------

function get_model_time_step()
 type(time_type) :: get_model_time_step

! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step

!------------------------------------------------------------------

subroutine get_state_meta_data(state_handle, index_in, location, var_type)
 type(ensemble_type), intent(in)  :: state_handle
 integer(i8),         intent(in)  :: index_in
 type(location_type), intent(out) :: location
 integer,             intent(out), optional :: var_type

! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

real(r8) :: lat, lon, height
integer  :: lon_index, lat_index, height_index, local_var, var_id

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, height_index, var_id=var_id)
! get from state structure
call get_state_kind(var_id, local_var)

lon = grid_longitude(lon_index)
lat = grid_latitude(lat_index)

if (local_var == KIND_ELECTRIC_POTENTIAL) then
   height   = 0.0_r8
   location = set_location(lon, lat, height, VERTISUNDEF)
else
   height   = levels(height_index)
   location = set_location(lon, lat, height, VERTISHEIGHT)
endif

if (debug > 5) print *, 'lon, lat, height = ', lon, lat, height

if (present(var_type)) then
   var_type = local_var
endif

end subroutine get_state_meta_data

!--------------------------------------------------------------------

function get_varid_from_kind(dart_kind)

integer, intent(in) :: dart_kind
integer             :: get_varid_from_kind

! given a kind, return what variable number it is

integer :: i

do i = 1, get_num_variables(domain_id)
   if (dart_kind == state_kinds_list(i)) then
      get_varid_from_kind = i
      return
   endif
end do

write(string1, *) 'Kind ', dart_kind, ' not found in state vector'
write(string2, *) 'AKA ', get_raw_obs_kind_name(dart_kind), ' not found in state vector'
call error_handler(E_MSG,'get_varid_from_kind', string1, &
                   source, revision, revdate, text2=string2)

get_varid_from_kind = -1

end function get_varid_from_kind


!------------------------------------------------------------------

subroutine get_state_kind(var_ind, var_type)
 integer, intent(in)  :: var_ind
 integer, intent(out) :: var_type

! Given an integer index into the state vector structure, returns the kind,
! and both the starting offset for this kind, as well as the offset into
! the block of this kind.

if ( .not. module_initialized ) call static_init_model

var_type = state_kinds_list(var_ind)

end subroutine get_state_kind

!------------------------------------------------------------------

subroutine end_model()

! Shutdown and clean-up.

! assume if one is allocated, they all were.  if no one ever
! called the init routine, don't try to dealloc something that
! was never alloc'd.
if (allocated(grid_latitude)) deallocate(grid_latitude, grid_longitude)
if (allocated(levels))        deallocate(levels)

end subroutine end_model

!------------------------------------------------------------------

function get_grid_template_fileid(filename)

character(len=*), intent(in) :: filename
integer get_grid_template_fileid

call nc_check( nf90_open(filename, NF90_NOWRITE, get_grid_template_fileid), &
                  'get_grid_template_fileid', 'open '//trim(filename))

end function get_grid_template_fileid

!------------------------------------------------------------------

subroutine get_grid_sizes(ncFileID, nlon, nlat, nvert)

integer, intent(in)  :: ncFileID
integer, intent(out) :: nlon
integer, intent(out) :: nlat
integer, intent(out) :: nvert

! netcdf variables
integer :: DimID

nvert = 1

call nc_check(NF90_inq_dimid(ncid=ncFileID, name='nphi', dimid=DimID), &
                           'get_grid_sizes','inq_dimid nphi')

call nc_check(NF90_inquire_dimension(ncFileID, DimID, len=nlon), &
              'get_grid_sizes', 'inqure_dimension nphi')

call nc_check(NF90_inq_dimid(ncid=ncFileID, name='nthe', dimid=DimID), &
                           'get_grid_sizes','inq_dimid nthe')

call nc_check(NF90_inquire_dimension(ncFileID, DimID, len=nlat), &
              'get_grid_sizes', 'inqure_dimension nthe')

end subroutine get_grid_sizes

!------------------------------------------------------------------

subroutine read_horiz_grid(ncFileID)

integer, intent(in)  :: ncFileID

! netcdf variables
integer :: VarID

call nc_check(NF90_inq_varid(ncFileID, 'nphi', VarID), &
              'read_horiz_grid', 'nphi inq_varid')

call nc_check(NF90_get_var(ncFileID, VarID, grid_longitude), &
              'read_horiz_grid', 'nphi get_var')

call nc_check(NF90_inq_varid(ncFileID, 'nthe', VarID), &
              'read_horiz_grid', 'nthe inq_varid')

call nc_check(NF90_get_var(ncFileID, VarID, grid_latitude), &
              'read_horiz_grid', 'nthe get_var')

grid_latitude(:) = 90.0_r8 - grid_latitude(:)

end subroutine read_horiz_grid

!------------------------------------------------------------------

subroutine read_vert_levels(ncFileID)

integer, intent(in)  :: ncFileID

! netcdf variables
integer :: VarID

!>@todo FIXME
!#! call nc_check(NF90_inq_varid(ncFileID, 'xxxx', VarID), &
!#!               'read_vert_levels', 'xxxx inq_varid')
!#! 
!#! call nc_check(NF90_get_var(ncFileID, VarID, levels), &
!#!               'read_vert_levels', 'xxxx get_var')
!#! 

levels(1) = 1

end subroutine read_vert_levels

!------------------------------------------------------------------

function nc_write_model_atts( ncFileID, model_mod_writes_state_variables ) result (ierr)
 integer, intent(in)  :: ncFileID      ! netCDF file identifier
 logical, intent(out) :: model_mod_writes_state_variables
 integer              :: ierr          ! return value of function

! TJH -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector.
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode 
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NvertDimID
integer :: lonVarID, latVarID, levelVarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: nmlVarID
logical :: has_openggcm_namelist

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

character(len=128)  :: filename

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

! have dart write out the state variables using the state structure
model_mod_writes_state_variables = .true. 

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'creation put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'source put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'revision put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'revdate put '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'openggcm' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

!----------------------------------------------------------------------------
! We need to output grid information
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', &
       len = nlon, dimid = NlonDimID),'nc_write_model_atts', 'lon def_dim '//trim(filename))
call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', &
       len = nlat, dimid = NlatDimID),'nc_write_model_atts', 'lon def_dim '//trim(filename))
call nc_check(nf90_def_dim(ncid=ncFileID, name='lev', &
       len = nvert, dimid = NvertDimID),'nc_write_model_atts', 'lev def_dim '//trim(filename))

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! Grid Longitudes
call nc_check(nf90_def_var(ncFileID,name='grid_longitude', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=lonVarID),&
              'nc_write_model_atts', 'grid_longitude def_var '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitudes of U,V grid'), &
              'nc_write_model_atts', 'grid_longitude long_name '//trim(filename))

! Grid Latitudes
call nc_check(nf90_def_var(ncFileID,name='grid_latitude', xtype=nf90_real, &
              dimids=(/ NlonDimID, NlatDimID /), varid=latVarID),&
              'nc_write_model_atts', 'grid_latitude def_var '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitudes of U,V grid'), &
              'nc_write_model_atts', 'grid_latitude long_name '//trim(filename))

! Levels
call nc_check(nf90_def_var(ncFileID,name='levels', xtype=nf90_real, &
              dimids=NvertDimID, varid= levelVarID), &
              'nc_write_model_atts', 'levels def_var '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levelVarID, 'long_name', 'depth at grid edges'), &
              'nc_write_model_atts', 'levels long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levelVarID, 'cartesian_axis', 'Z'),   &
              'nc_write_model_atts', 'levels cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levelVarID, 'units', 'meters'),  &
              'nc_write_model_atts', 'levels units '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levelVarID, 'positive', 'down'),  &
              'nc_write_model_atts', 'levels units '//trim(filename))
call nc_check(nf90_put_att(ncFileID, levelVarID, 'comment', &
               'more positive is closer to the center of the earth'),  &
              'nc_write_model_atts', 'levels comment '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_check(nf90_put_var(ncFileID, lonVarID, grid_longitude ), &
             'nc_write_model_atts', 'grid_longitude put_var '//trim(filename))
call nc_check(nf90_put_var(ncFileID, latVarID, grid_latitude ), &
             'nc_write_model_atts', 'grid_latitude put_var '//trim(filename))
call nc_check(nf90_put_var(ncFileID, levelVarID, levels ), &
             'nc_write_model_atts', 'levels put_var '//trim(filename))

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_openggcm_namelist) then
   call file_to_text('openggcm_in', textblock)
   call nc_check(nf90_put_var(ncFileID, nmlVarID, textblock ), &
                 'nc_write_model_atts', 'put_var nmlVarID')
   deallocate(textblock)
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts

!------------------------------------------------------------------

function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
 integer,                intent(in) :: ncFileID      ! netCDF file identifier
 real(r8), dimension(:), intent(in) :: statevec
 integer,                intent(in) :: copyindex
 integer,                intent(in) :: timeindex
 integer                            :: ierr          ! return value of function

! DEFINE ANY EXTRA VARIABLES OTHER THAN THOSE IN THE STATE

if ( .not. module_initialized ) call static_init_model

end function nc_write_model_vars

!------------------------------------------------------------------

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

 type(ensemble_type), intent(inout) :: state_ens_handle
 integer,             intent(in)    :: ens_size
 real(r8),            intent(in)    :: pert_amp
 logical,             intent(out)   :: interf_provided

! Perturbs state copies for generating initial ensembles.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

integer     :: var_type
integer     :: j,i 
integer(i8) :: dart_index

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq
logical :: lanai_bitwise

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.
lanai_bitwise = .true.

if (lanai_bitwise) then
   call pert_model_copies_bitwise_lanai(state_ens_handle, ens_size)
else

   ! Initialize random number sequence
   call init_random_seq(random_seq, my_task_id())

   ! only perturb the actual ocean cells; leave the land and
   ! ocean floor values alone.
   do i=1,state_ens_handle%my_num_vars
      dart_index = state_ens_handle%my_vars(i)
      do j=1, ens_size
         state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
            state_ens_handle%copies(j,i), &
            model_perturbation_amplitude)
      enddo
   enddo

endif


end subroutine pert_model_copies

!------------------------------------------------------------------
!> Perturb the state such that the perturbation is bitwise with
!> a perturbed Lanai.
!> Note:
!> * This is not bitwise with itself (like Lanai) if task_count < ens_size
!> * This is very slow because you have a loop around the length of the
!> state inside a loop around the ensemble.
!> If a task has more than one copy then the random number
!> sequence continues from the end of one copy to the start of the other.
subroutine pert_model_copies_bitwise_lanai(ens_handle, ens_size)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size

type(random_seq_type) :: r(ens_size)
integer     :: i ! loop variable
integer(i8) :: j ! loop variable
real(r8)    :: random_number
integer     :: var_type
integer     :: sequence_to_use
integer     :: owner, owners_index
integer     :: copy_owner

do i = 1, min(ens_size, task_count())
   call init_random_seq(r(i), map_pe_to_task(ens_handle,i-1)) ! my_task_id in Lanai
   sequence_to_use = i 
enddo


do i = 1, ens_size

   call get_copy_owner_index(i, copy_owner, owners_index) ! owners index not used. Distribution type 1
   sequence_to_use = copy_owner + 1
   do j = 1, ens_handle%num_vars

     random_number = random_gaussian(r(sequence_to_use), 0.0_r8, model_perturbation_amplitude)
     call get_var_owner_index(j, owner, owners_index)
     if (ens_handle%my_pe==owner) then
        ens_handle%copies(i, owners_index) = ens_handle%copies(i, owners_index) + random_number
     endif

   enddo

enddo

end subroutine pert_model_copies_bitwise_lanai

!------------------------------------------------------------------

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs, obs_kind, num_close, close_ind, dist, state_handle)
 type(ensemble_type),               intent(in) :: state_handle
 type(get_close_type),              intent(in) :: gc
 type(location_type),               intent(in) :: base_obs_loc
 integer,                           intent(in) :: base_obs_kind
 type(location_type), dimension(:), intent(in) :: obs
 integer,             dimension(:), intent(in) :: obs_kind
 integer,                           intent(out):: num_close
 integer,             dimension(:), intent(out):: close_ind
 real(r8),            dimension(:), intent(out):: dist !does this need to be optional? It is not in WRF

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

integer :: t_ind, k

! Initialize variables to missing status

num_close = 0
close_ind = -99
!if (present(dist)) dist = 1.0e9   !something big and positive (far away)
dist = 1.0e9   !something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet 
! (for obs).

call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs, obs_kind, &
                       num_close, close_ind)

end subroutine get_close_obs

!------------------------------------------------------------------

subroutine do_interp(state_handle, ens_size, base_offset, hgt_bot, hgt_top, hgt_fract, &
                     llon, llat, obs_type, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_Size
integer(i8),         intent(in) :: base_offset
integer,             intent(in) :: hgt_bot, hgt_top
real(r8),            intent(in) :: hgt_fract, llon, llat
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size)
integer,           intent(out) :: istatus(ens_size)
 
! do a 2d horizontal interpolation for the value at the bottom level, 
! then again for the top level, then do a linear interpolation in the 
! vertical to get the final value.

integer(i8) :: offset
real(r8)    :: bot_val(ens_size), top_val(ens_size)
integer     :: e
integer     :: temp_status(ens_size)


! Find the base location for the bottom height and interpolate horizontally 
!  on this level.  Do bottom first in case it is below the ocean floor; can
!  avoid the second horizontal interpolation.
offset = base_offset + (hgt_bot - 1) * nlon * nlat
if (debug > 6) &
   print *, 'bot, field, abs offset: ', hgt_bot, base_offset, offset

call lon_lat_interpolate(state_handle, ens_size, offset, llon, llat, obs_type, hgt_bot, bot_val, temp_status)
! Failed istatus from interpolate means give up
istatus = temp_status
if(all(istatus /= 0)) return
if (debug > 6) &
   print *, 'bot_val = ', bot_val

! Find the base location for the top height and interpolate horizontally 
!  on this level.
offset = base_offset + (hgt_top - 1) * nlon * nlat
if (debug > 6) &
   print *, 'top, field, abs offset: ', hgt_top, base_offset, offset

call lon_lat_interpolate(state_handle, ens_size, offset, llon, llat, obs_type, hgt_top, top_val, temp_status)
do e = 1, ens_size
   if(temp_status(e) /= 0) istatus(e) = temp_status(e)
enddo
! Failed istatus from interpolate means give up
if(all(istatus /= 0)) return
if (debug > 6) &
   print *, 'top_val = ', top_val

! Then weight them by the vertical fraction and return
expected_obs = bot_val + hgt_fract * (top_val - bot_val)
if (debug > 2) print *, 'do_interp: interp val = ',expected_obs


end subroutine do_interp

!--------------------------------------------------------------------
!> construct restart file name for reading
function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_in

! stub is found in input.nml io_filename_nml
! restart files typically are of the form
! openggcm.r0001.nc
! write(construct_file_name_in, '(A, i4.4, A)') trim(stub), copy, ".nc"
write(construct_file_name_in, '(A, A)') trim(stub), ".nc"


end function construct_file_name_in

!--------------------------------------------------------------------
!> read the time from template file
function read_model_time(filename)

character(len=1024) :: filename
type(time_type) :: read_model_time

! netcdf variables
integer :: ncFileID, VarID

! fractional days
real(r8) :: days

! interger conversion of fractional days
integer :: seconds, whole_days

type(time_type) :: base_time

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',msgstring,source,revision,revdate)
endif

call nc_check( nf90_open(filename, NF90_NOWRITE, ncFileID), &
                  'read_model_time', 'open '//filename )

call nc_check(NF90_inq_varid(ncFileID, 'time', VarID), &
              'read_model_time', 'time inq_varid')

call nc_check(NF90_get_var(ncFileID, VarID, days), &
              'read_model_time', 'time get_var')

!>@todo FIXME : Should be grabbing this information from the attributes
!>              harcoded for now.
base_time = set_date(1900,1,1,0,0)

whole_days = floor(days)
seconds    = nint( (days - whole_days)*86400 )

read_model_time = base_time + set_time(seconds, whole_days)

end function read_model_time

!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod
function query_vert_localization_coord()

integer :: query_vert_localization_coord

query_vert_localization_coord = 1 ! any old value

end function query_vert_localization_coord

!--------------------------------------------------------------------
!> This is used in the filter_assim. The vertical conversion is done using the 
!> mean state.
!> Calling this is a waste of time
subroutine vert_convert(state_handle, location, obs_kind, istatus)

type(ensemble_type), intent(in)  :: state_handle
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
integer,             intent(out) :: istatus

istatus = 0

end subroutine vert_convert

!------------------------------------------------------------------
!> Verify that the namelist was filled in correctly, and check
!> that there are valid entries for the dart_kind. 
!> Returns a table with columns:  
!>
!>    netcdf_variable_name ; dart_kind_string ; update_string
!>
subroutine verify_state_variables( state_variables, ngood, table, kind_list, update_var )

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: kind_list(:)   ! kind number
logical, optional, intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0

if ( state_variables(1) == ' ' ) then ! no model_state_variables namelist provided
   string1 = 'model_nml:model_state_variables not specified using default variables'
   call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
endif

MyLoop : do i = 1, nrows

   varname = trim(state_variables(3*i -2))
   dartstr = trim(state_variables(3*i -1))
   update  = trim(state_variables(3*i   ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' ) then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   kind_list(i) = get_raw_obs_kind_index(dartstr)
   if( kind_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif
   
   ! Make sure the update variable has a valid name

   if ( present(update_var) )then
      SELECT CASE (update)
         CASE ('UPDATE')
            update_var(i) = .true.
         CASE ('NO_COPY_BACK')
            update_var(i) = .false.
         CASE DEFAULT
            write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
            write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
            call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate, text2=string2)
      END SELECT
   endif

   ! Record the contents of the DART state vector

   if (do_output()) then
      write(string1,'(A,I2,6A)') 'variable ',i,' is ',trim(varname), ', ', trim(dartstr), ', ', trim(update)
      call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine verify_state_variables

!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
