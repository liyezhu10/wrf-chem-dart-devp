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
                             VERTISHEIGHT, get_location, vert_is_height,       &
                             vert_is_level, vert_is_surface,                   &
                             loc_get_close_obs => get_close_obs, get_close_type
use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text
use     obs_kind_mod, only : KIND_TEMPERATURE, KIND_SALINITY, KIND_DRY_LAND,   &
                             KIND_U_CURRENT_COMPONENT,KIND_V_CURRENT_COMPONENT,&
                             KIND_SEA_SURFACE_HEIGHT, KIND_SEA_SURFACE_PRESSURE,&
                             KIND_POTENTIAL_TEMPERATURE, get_raw_obs_kind_index,&
                             get_raw_obs_kind_name, paramname_length 
use mpi_utilities_mod, only: my_task_id, task_count
use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian
use      dart_openggcm_mod, only: set_model_time_step,                              &
                             get_horiz_grid_dims, get_vert_grid_dim,           &
                             read_horiz_grid, read_topography, read_vert_grid, &
                             get_openggcm_restart_filename

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


! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.
public :: get_gridsize, restart_file_to_sv, sv_to_restart_file, &
          get_openggcm_restart_filename, test_interpolation

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
logical  :: output_state_vector = .true.
integer  :: assimilation_period_days = 1
integer  :: assimilation_period_seconds = 0
real(r8) :: model_perturbation_amplitude = 0.2
character(len=paramname_length) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 0   ! turn up for more and more debug messages

namelist /model_nml/  &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   model_state_variables,       &
   debug

!------------------------------------------------------------------
!
! The default DART state vector (control vector) will consist of:  S, T, U, V, PSURF
! (Salinity, Temperature, U velocity, V velocity, Sea Surface Height).
! S, T are 3D arrays, located at cell centers.  U,V are at grid cell corners.
! PSURF is a 2D field (X,Y only).  The Z direction is downward. 
! 
! Additional variables can be read into the state vector using the 
! model_state_variables namelist by specifying the netcdf variable name
! dart kind string and an update string.  Currently the update string
! is not being used.
!------------------------------------------------------------------

! Number of fields in the state vector
integer :: nfields

! Grid parameters - the values will be read from the
! openggcm netcdf restart file

! the size of the grid
integer :: nlon, nlat, nvert

real(r8) :: grid_longitude(:), grid_latitude(:)
real(r8) :: levels(:)

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

integer :: iunit, io
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

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(msgstring,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',msgstring,source,revision,revdate)


! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

call get_horiz_grid_dims(nlon, nlat)
call get_vert_grid_dim(nvert)

! Allocate space for grid variables. 
allocate(grid_longitude(nlon), grid_latitude(nlat))
allocate(levels(nvert))

! Fill them in. get the size from the array size.
call read_horiz_grid(grid_longitude, grid_latitude)
call read_vert_grid(levels)

if (debug > 2) call write_grid_netcdf() ! DEBUG only

! verify that the model_state_variables namelist was filled in correctly.  
! returns variable_table which has variable names, kinds and update strings.
call verify_state_variables(model_state_variables, nfields, variable_table, state_kinds_list, update_var_list)

!> @todo  - need input filename, hardcode for now to openggcm.nc
domain_id = add_domain('openggcm.nc', nfields, &
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
logical        :: convert_to_ssh
integer        :: e

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
   ! convert the level index to an actual depth 
   ind = nint(loc_array(3))
   if ( (ind < 1) .or. (ind > size(zc)) ) then 
      istatus = 11
      return
   else
      lheight = level(ind)
   endif
else   ! if pressure or undefined, we don't know what to do
   istatus = 17
   return
endif

! kind (in-situ) temperature is a combination of potential temp,
! salinity, and pressure based on depth.  call a routine that
! interpolates all three, does the conversion, and returns the
! sensible/in-situ temperature.
if(obs_type == KIND_TEMPERATURE) then
   ! we know how to interpolate this from potential temp,
   ! salinity, and pressure based on depth.
   call compute_temperature(state_handle, ens_size, llon, llat, lheight, expected_obs, istatus)
   if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus
   return
endif


! The following kinds are either in the state vector (so you
! can simply interpolate to find the value) or they are a simple
! transformation of something in the state vector.

convert_to_ssh = .FALSE.

SELECT CASE (obs_type)
   CASE (KIND_SALINITY,              &
         KIND_POTENTIAL_TEMPERATURE, &
         KIND_U_CURRENT_COMPONENT,   &
         KIND_V_CURRENT_COMPONENT,   &
         KIND_SEA_SURFACE_PRESSURE)
      base_offset = get_index_start(domain_id, get_varid_from_kind(obs_type))

   CASE (KIND_SEA_SURFACE_HEIGHT)
      base_offset = get_index_start(domain_id, get_varid_from_kind(KIND_SEA_SURFACE_PRESSURE))
      convert_to_ssh = .TRUE. ! simple linear transform of PSURF

   CASE DEFAULT
      ! Not a legal type for interpolation, return istatus error
      istatus = 15
      return

END SELECT

! For Sea Surface Height or Pressure don't need the vertical coordinate
! SSP needs to be converted to a SSH if height is required.
if( vert_is_surface(location) ) then
   ! HK CHECK surface observations
   call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, obs_type, 1, expected_obs, istatus)
   do e = 1, ens_size
      if (convert_to_ssh .and. (istatus(e) == 0)) then !HK why check istatus?
         expected_obs(e) = expected_obs(e) / 980.6_r8   ! openggcm uses CGS units
      endif
   enddo

   if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus
   return
endif

! Get the bounding vertical levels and the fraction between bottom and top
call height_bounds(lheight, Nz, ZC, hgt_bot, hgt_top, hgt_fract, hstatus)
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
integer  :: x_ind, y_ind
real(r8) :: x_corners(4), y_corners(4)
real(r8) :: p(4,ens_size), xbot(ens_size), xtop(ens_size)
real(r8) :: lon_fract, lat_fract
logical  :: masked
integer  :: quad_status
integer  :: e

if ( .not. module_initialized ) call static_init_model

! Succesful return has istatus of 0
istatus = 0

! Get the lower left corner for either grid type
if(dipole_grid) then
   ! Figure out which of the regular grid boxes this is in
   call get_reg_box_indices(lon, lat, x_ind, y_ind)

   ! Is this on the U or T grid?
   if(is_on_ugrid(var_type)) then
      ! On U grid
      num_inds =  u_dipole_num  (x_ind, y_ind)
      start_ind = u_dipole_start(x_ind, y_ind)

      ! If there are no quads overlapping, can't do interpolation
      if(num_inds == 0) then
         istatus = 1
         return
      endif

      ! Search the list of quads to see if (lon, lat) is in one
      call get_dipole_quad(lon, lat, ulon, ulat, num_inds, start_ind, &
         u_dipole_lon_list, u_dipole_lat_list, lon_bot, lat_bot, quad_status)
      ! Fail on bad istatus return
      if(quad_status /= 0) then
         istatus = quad_status
         return
      endif

      ! Getting corners for accurate interpolation
      call get_quad_corners(ulon, lon_bot, lat_bot, x_corners)
      call get_quad_corners(ulat, lon_bot, lat_bot, y_corners)

      ! Fail if point is in one of the U boxes that go through the
      ! pole (this could be fixed up if necessary)
      if(lat_bot == u_pole_y .and. (lon_bot == pole_x -1 .or. &
         lon_bot == pole_x)) then
         istatus = 4
         return
      endif

   else
      ! On T grid
      num_inds =  t_dipole_num  (x_ind, y_ind)
      start_ind = t_dipole_start(x_ind, y_ind)

      ! If there are no quads overlapping, can't do interpolation
      if(num_inds == 0) then
         istatus = 1
         return
      endif

      call get_dipole_quad(lon, lat, tlon, tlat, num_inds, start_ind, &
         t_dipole_lon_list, t_dipole_lat_list, lon_bot, lat_bot, quad_status)

      ! Fail on bad istatus return
      if(quad_status /= 0) then
         istatus = quad_status
         return
      endif

      ! Fail if point is in T box that covers pole
      if(lon_bot == pole_x .and. lat_bot == t_pole_y) then
         istatus = 5
         return
      endif

      ! Getting corners for accurate interpolation
      call get_quad_corners(tlon, lon_bot, lat_bot, x_corners)
      call get_quad_corners(tlat, lon_bot, lat_bot, y_corners)

   endif

else
   ! This is an irregular grid
   ! U and V are on velocity grid
   if (is_on_ugrid(var_type)) then
      ! Get the corner indices and the fraction of the distance between
      call get_irreg_box(lon, lat, ulon, ulat, &
         lon_bot, lat_bot, lon_fract, lat_fract, quad_status)
   else
      ! Eta, T and S are on the T grid
      ! Get the corner indices
      call get_irreg_box(lon, lat, tlon, tlat, &
         lon_bot, lat_bot, lon_fract, lat_fract, quad_status)
   endif

   ! Return passing through error status
   if(quad_status /= 0) then
      istatus = quad_status
      return
   endif

endif

! Find the indices to get the values for interpolating
lat_top = lat_bot + 1
if(lat_top > ny) then
   istatus = 2
   return
endif

! Watch for wraparound in longitude
lon_top = lon_bot + 1
if(lon_top > nx) lon_top = 1

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
p(1, :) = get_val(lon_bot, lat_bot, nx, state_handle, offset, ens_size, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif


p(2, :) = get_val(lon_top, lat_bot, nx, state_handle, offset, ens_size, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif



p(3, :) = get_val(lon_top, lat_top, nx, state_handle, offset, ens_size, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

p(4, :) = get_val(lon_bot, lat_top, nx, state_handle, offset, ens_size, var_type, height, masked)
if(masked) then
   istatus = 3
   return
endif

! Full bilinear interpolation for quads
if(dipole_grid) then
   do e = 1, ens_size
      call quad_bilinear_interp(lon, lat, x_corners, y_corners, p(:,e), ens_size, expected_obs(e))
   enddo
else
   ! Rectangular biliear interpolation
   xbot = p(1, :) + lon_fract * (p(2, :) - p(1, :))
   xtop = p(4, :) + lon_fract * (p(3, :) - p(4, :))
   ! Now interpolate in latitude
   expected_obs = xbot + lat_fract * (xtop - xbot)
endif

end subroutine lon_lat_interpolate

!------------------------------------------------------------

function get_val(lon_index, lat_index, nlon, state_handle, offset, ens_size, var_type, height, masked)
 integer,             intent(in)  :: lon_index, lat_index, nlon, var_type, height
 type(ensemble_type), intent(in)  :: state_handle
 integer(i8),         intent(in)  :: offset
 integer,             intent(in)  :: ens_size
 logical,             intent(out) :: masked

 real(r8)    :: get_val(ens_size)
 integer(i8) :: state_index

! Returns the value from a single level array given the lat and lon indices
! 'masked' returns true if this is NOT a valid grid location (e.g. land, or
! below the ocean floor in shallower areas).

if ( .not. module_initialized ) call static_init_model

! check the land/ocean bottom map and return if not valid water cell.
if(is_dry_land(var_type, lon_index, lat_index, height)) then
   masked = .true.
   get_val = MISSING_R8
   return
endif

! state index must be 8byte integer
state_index = int(lat_index - 1,i8)*int(nlon,i8) + int(lon_index,i8) + int(offset-1,i8)

! Layout has lons varying most rapidly
!get_val = x((lat_index - 1) * nlon + lon_index)
! The x above is only a horizontal slice, not the whole state.   HK WHY -1?
get_val = get_state(state_index, state_handle)

! this is a valid ocean water cell, not land or below ocean floor
masked = .false.

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

! The zc array contains the depths of the center of the vertical grid boxes
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

real(r8) :: lat, lon, depth
integer :: lon_index, lat_index, depth_index, local_var, var_id

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, depth_index, var_id=var_id)
call get_state_kind(var_id, local_var)

if (local_var == KIND_SEA_SURFACE_HEIGHT) then
   depth = 0.0_r8
else
   depth = ZC(depth_index)
endif

if (debug > 5) print *, 'lon, lat, depth = ', lon, lat, depth

location = set_location(lon, lat, depth, VERTISHEIGHT)

if (present(var_type)) then
   var_type = local_var
   if(is_dry_land(var_type, lon_index, lat_index, depth_index)) then
      var_type = KIND_DRY_LAND
   endif
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

subroutine get_state_kind_inc_dry(index_in, var_type)
 integer(i8), intent(in)  :: index_in
 integer,     intent(out) :: var_type

! Given an integer index into the state vector structure, returns the
! type, taking into account the ocean bottom and dry land.

integer :: lon_index, lat_index, depth_index, var_id

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, depth_index, var_id=var_id)
call get_state_kind(var_id, var_type)

! if on land or below ocean floor, replace type with dry land.
if(is_dry_land(var_type, lon_index, lat_index, depth_index)) then
   var_type = KIND_DRY_LAND
endif

end subroutine get_state_kind_inc_dry

!------------------------------------------------------------------

subroutine end_model()

! Shutdown and clean-up.

! assume if one is allocated, they all were.  if no one ever
! called the init routine, don't try to dealloc something that
! was never alloc'd.
if (allocated(ULAT)) deallocate(ULAT, ULON, TLAT, TLON, KMT, KMU, HT, HU)
if (allocated(ZC))   deallocate(ZC, ZG, pressure)

end subroutine end_model

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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!----------------------------------------------------------------------
! variables if we just blast out one long state vector
!----------------------------------------------------------------------

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables if we parse the state vector into prognostic variables.
!----------------------------------------------------------------------

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NzDimID
integer :: ulonVarID, ulatVarID, tlonVarID, tlatVarID, ZGVarID, ZCVarID
integer :: KMTVarID, KMUVarID

! for the prognostic variables
integer :: SVarID, TVarID, UVarID, VVarID, PSURFVarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
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

integer     :: i
character(len=128)  :: filename

integer  :: model_size_i4 ! this is for checking model_size

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly
model_mod_writes_state_variables = .true. 

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
                           'nc_write_model_atts','inq_dimid NMLlinelen')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(msgstring,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', msgstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
! JH -- nf90_def_dim is expecting a lenght that is i4.  Here we type cast size and
!    check if the values are the same.  In the case where model_size is larger
!    than the largest i4 integer we error out.
!-------------------------------------------------------------------------------

model_size_i4 = int(model_size,i4) 
if (model_size_i4 /= model_size) then
   write(msgstring,*)'model_size =  ', model_size, ' is too big to write ', &
             ' diagnostic files.'
   call error_handler(E_ERR,'nc_write_model_atts', msgstring, source, revision, revdate)
endif

call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size_i4, &
        dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

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
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('openggcm_in', nlines, linelen)
if (nlines > 0) then
  has_openggcm_namelist = .true.
else
  has_openggcm_namelist = .false.
endif

if (debug > 0)    print *, 'openggcm namelist: nlines, linelen = ', nlines, linelen
  
if (has_openggcm_namelist) then 
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncFileID,name='openggcm_in', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var openggcm_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of openggcm_in namelist'), 'nc_write_model_atts', 'put_att openggcm_in')

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                 dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                 'statevariable def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                 'nc_write_model_atts','statevariable long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                 'nc_write_model_atts', 'statevariable units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1_i8,model_size /)),&
                 'nc_write_model_atts', 'statevariable valid_range '//trim(filename))

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                 dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                 'nc_write_model_atts','state def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                 'nc_write_model_atts', 'state long_name '//trim(filename))

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                 'nc_write_model_atts', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to output the prognostic variables.
   !----------------------------------------------------------------------------
   ! Define the new dimensions IDs
   !----------------------------------------------------------------------------
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='i', &
          len = Nx, dimid = NlonDimID),'nc_write_model_atts', 'i def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='j', &
          len = Ny, dimid = NlatDimID),'nc_write_model_atts', 'j def_dim '//trim(filename))
   call nc_check(nf90_def_dim(ncid=ncFileID, name='k', &
          len = Nz, dimid =   NzDimID),'nc_write_model_atts', 'k def_dim '//trim(filename))
   
   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------


   ! U,V Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='ULON', xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=ulonVarID),&
                 'nc_write_model_atts', 'ULON def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitudes of U,V grid'), &
                 'nc_write_model_atts', 'ULON long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'cartesian_axis', 'X'),  &
                 'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'ULON units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'ULON valid_range '//trim(filename))

   ! U,V Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='ULAT', xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=ulatVarID),&
                 'nc_write_model_atts', 'ULAT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitudes of U,V grid'), &
                 'nc_write_model_atts', 'ULAT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'ULAT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  ulatVarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'ULAT valid_range '//trim(filename))

   ! S,T,PSURF Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='TLON', xtype=nf90_real, &
                 dimids=(/ NlonDimID, NlatDimID /), varid=tlonVarID),&
                 'nc_write_model_atts', 'TLON def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, 'long_name', 'longitudes of S,T,... grid'), &
                 'nc_write_model_atts', 'TLON long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, 'cartesian_axis', 'X'),   &
                 'nc_write_model_atts', 'TLON cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, 'units', 'degrees_east'),  &
                 'nc_write_model_atts', 'TLON units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlonVarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'TLON valid_range '//trim(filename))


   ! S,T,PSURF Grid (center) Latitudes
   call nc_check(nf90_def_var(ncFileID,name='TLAT', xtype=nf90_real, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=tlatVarID), &
                 'nc_write_model_atts', 'TLAT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, 'long_name', 'latitudes of S,T, ... grid'), &
                 'nc_write_model_atts', 'TLAT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, 'cartesian_axis', 'Y'),   &
                 'nc_write_model_atts', 'TLAT cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'TLAT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, tlatVarID, 'valid_range', (/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'TLAT valid_range '//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name='ZG', xtype=nf90_real, &
                 dimids=NzDimID, varid= ZGVarID), &
                 'nc_write_model_atts', 'ZG def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, 'long_name', 'depth at grid edges'), &
                 'nc_write_model_atts', 'ZG long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'ZG cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ZG units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'ZG units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZGVarID, 'comment', &
                  'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'ZG comment '//trim(filename))

   ! Depths
   call nc_check(nf90_def_var(ncFileID,name='ZC',xtype=nf90_real,dimids=NzDimID,varid=ZCVarID), &
                 'nc_write_model_atts', 'ZC def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'long_name', 'depth at grid centroids'), &
                 'nc_write_model_atts', 'ZC long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'ZC cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'ZC units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, ZCVarID, 'comment', &
                  'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'ZC comment '//trim(filename))

   ! Depth mask
   call nc_check(nf90_def_var(ncFileID,name='KMT',xtype=nf90_int, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=KMTVarID), &
                 'nc_write_model_atts', 'KMT def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, 'long_name', 'lowest valid depth index at grid centroids'), &
                 'nc_write_model_atts', 'KMT long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, 'units', 'levels'),  &
                 'nc_write_model_atts', 'KMT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'KMT units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMTVarID, 'comment', &
                  'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'KMT comment '//trim(filename))

   ! Depth mask
   call nc_check(nf90_def_var(ncFileID,name='KMU',xtype=nf90_int, &
                 dimids= (/ NlonDimID, NlatDimID /), varid=KMUVarID), &
                 'nc_write_model_atts', 'KMU def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, 'long_name', 'lowest valid depth index at grid corners'), &
                 'nc_write_model_atts', 'KMU long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, 'units', 'levels'),  &
                 'nc_write_model_atts', 'KMU units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, 'positive', 'down'),  &
                 'nc_write_model_atts', 'KMU units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID, KMUVarID, 'comment', &
                  'more positive is closer to the center of the earth'),  &
                 'nc_write_model_atts', 'KMU comment '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------
   
   !>@todo JH : If we store the variable attributes in a structure we can simplly
   ! loop over all of the variables and output prognostic variables and attributes
   !> For now we are only writting the default variables if they exist.
   if ( get_varid_from_kind(KIND_SALINITY) > 0 ) then
      call nc_check(nf90_def_var(ncid=ncFileID, name='SALT', xtype=nf90_real, &
            dimids = (/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=SVarID),&
            'nc_write_model_atts', 'S def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, SVarID, 'long_name', 'salinity'), &
            'nc_write_model_atts', 'S long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, SVarID, 'units', 'kg/kg'), &
            'nc_write_model_atts', 'S units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, SVarID, 'missing_value', NF90_FILL_REAL), &
            'nc_write_model_atts', 'S missing '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, SVarID, '_FillValue', NF90_FILL_REAL), &
            'nc_write_model_atts', 'S fill '//trim(filename))
   endif

   if ( get_varid_from_kind(KIND_POTENTIAL_TEMPERATURE) > 0 ) then
      call nc_check(nf90_def_var(ncid=ncFileID, name='TEMP', xtype=nf90_real, &
            dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=TVarID),&
            'nc_write_model_atts', 'T def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, TVarID, 'long_name', 'Potential Temperature'), &
            'nc_write_model_atts', 'T long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, TVarID, 'units', 'deg C'), &
            'nc_write_model_atts', 'T units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, TVarID, 'units_long_name', 'degrees celsius'), &
            'nc_write_model_atts', 'T units_long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, TVarID, 'missing_value', NF90_FILL_REAL), &
            'nc_write_model_atts', 'T missing '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, TVarID, '_FillValue', NF90_FILL_REAL), &
            'nc_write_model_atts', 'T fill '//trim(filename))
   endif


   if ( get_varid_from_kind(KIND_U_CURRENT_COMPONENT) > 0 ) then
      call nc_check(nf90_def_var(ncid=ncFileID, name='UVEL', xtype=nf90_real, &
            dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=UVarID),&
            'nc_write_model_atts', 'U def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, UVarID, 'long_name', 'U velocity'), &
            'nc_write_model_atts', 'U long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, UVarID, 'units', 'cm/s'), &
            'nc_write_model_atts', 'U units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, UVarID, 'units_long_name', 'centimeters per second'), &
            'nc_write_model_atts', 'U units_long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, UVarID, 'missing_value', NF90_FILL_REAL), &
            'nc_write_model_atts', 'U missing '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, UVarID, '_FillValue', NF90_FILL_REAL), &
            'nc_write_model_atts', 'U fill '//trim(filename))
   endif


   if ( get_varid_from_kind(KIND_V_CURRENT_COMPONENT) > 0 ) then
      call nc_check(nf90_def_var(ncid=ncFileID, name='VVEL', xtype=nf90_real, &
            dimids=(/NlonDimID,NlatDimID,NzDimID,MemberDimID,unlimitedDimID/),varid=VVarID),&
            'nc_write_model_atts', 'V def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, VVarID, 'long_name', 'V Velocity'), &
            'nc_write_model_atts', 'V long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, VVarID, 'units', 'cm/s'), &
            'nc_write_model_atts', 'V units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, VVarID, 'units_long_name', 'centimeters per second'), &
            'nc_write_model_atts', 'V units_long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, VVarID, 'missing_value', NF90_FILL_REAL), &
            'nc_write_model_atts', 'V missing '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, VVarID, '_FillValue', NF90_FILL_REAL), &
            'nc_write_model_atts', 'V fill '//trim(filename))
   endif

   if ( get_varid_from_kind(KIND_SEA_SURFACE_PRESSURE) > 0 ) then
      call nc_check(nf90_def_var(ncid=ncFileID, name='PSURF', xtype=nf90_real, &
            dimids=(/NlonDimID,NlatDimID,MemberDimID,unlimitedDimID/),varid=PSURFVarID), &
            'nc_write_model_atts', 'PSURF def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'long_name', 'surface pressure'), &
            'nc_write_model_atts', 'PSURF long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'units', 'dyne/cm2'), &
            'nc_write_model_atts', 'PSURF units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, PSURFVarID, 'missing_value', NF90_FILL_REAL), &
            'nc_write_model_atts', 'PSURF missing '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, PSURFVarID, '_FillValue', NF90_FILL_REAL), &
            'nc_write_model_atts', 'PSURF fill '//trim(filename))
   endif

   ! Finished with dimension/variable definitions, must end 'define' mode to fill.

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncFileID, ulonVarID, ULON ), &
                'nc_write_model_atts', 'ULON put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ulatVarID, ULAT ), &
                'nc_write_model_atts', 'ULAT put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, tlonVarID, TLON ), &
                'nc_write_model_atts', 'TLON put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, tlatVarID, TLAT ), &
                'nc_write_model_atts', 'TLAT put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZGVarID, ZG ), &
                'nc_write_model_atts', 'ZG put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, ZCVarID, ZC ), &
                'nc_write_model_atts', 'ZC put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, KMTVarID, KMT ), &
                'nc_write_model_atts', 'KMT put_var '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, KMUVarID, KMU ), &
                'nc_write_model_atts', 'KMU put_var '//trim(filename))

endif

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

! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: VarID

real(r8), dimension(Nx,Ny,Nz) :: data_3d
real(r8), dimension(Nx,Ny)    :: data_2d
character(len=128)  :: filename

integer :: S_index ,T_index ,U_index, V_index, PSURF_index 

if ( .not. module_initialized ) call static_init_model

ierr = -1 ! assume things go poorly

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncFileID', ncFileID

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
              'nc_write_model_vars', 'inquire '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,statevec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   ! Replace missing values (0.0) with netcdf missing value.
   ! Staggered grid causes some logistical problems.
   !----------------------------------------------------------------------------

   !>@todo JH: If we store the variable attributes in a structure we can simplly
   !> loop over all of the variables and output prognostic variables and attributes.
   !> For now we are only writting the default variables if they exist.
   S_index = get_varid_from_kind(KIND_SALINITY)
   if ( S_index > 0 ) then
      !>@todo JH: do not need to use vector_to_prog_var to reshape variables for
      !> netcdf file.  you can simply use the count=(dim1, dim2, dim3) in the 
      !> netcdf optional arguments.
      call vector_to_prog_var(statevec, S_index, data_3d)
      where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
      call nc_check(NF90_inq_varid(ncFileID, 'SALT', VarID), &
                   'nc_write_model_vars', 'S inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                   'nc_write_model_vars', 'S put_var '//trim(filename))
   endif

   T_index = get_varid_from_kind(KIND_POTENTIAL_TEMPERATURE)
   if ( T_index > 0 ) then
      call vector_to_prog_var(statevec, T_index, data_3d)
      where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
      call nc_check(NF90_inq_varid(ncFileID, 'TEMP', VarID), &
                   'nc_write_model_vars', 'T inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                   'nc_write_model_vars', 'T put_var '//trim(filename))
   endif

   U_index = get_varid_from_kind(KIND_U_CURRENT_COMPONENT)
   if ( U_index > 0 ) then
      call vector_to_prog_var(statevec, U_index, data_3d)
      where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
      call nc_check(NF90_inq_varid(ncFileID, 'UVEL', VarID), &
                   'nc_write_model_vars', 'U inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                   'nc_write_model_vars', 'U put_var '//trim(filename))
   endif

   V_index = get_varid_from_kind(KIND_V_CURRENT_COMPONENT)
   if ( V_index > 0 ) then
      call vector_to_prog_var(statevec, V_index, data_3d)
      where (data_3d == 0.0_r8) data_3d = NF90_FILL_REAL
      call nc_check(NF90_inq_varid(ncFileID, 'VVEL', VarID), &
                   'nc_write_model_vars', 'V inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID,VarID,data_3d,start=(/1,1,1,copyindex,timeindex/)),&
                   'nc_write_model_vars', 'V put_var '//trim(filename))
   endif

   PSURF_index = get_varid_from_kind(KIND_SEA_SURFACE_PRESSURE)
   if ( PSURF_index > 0 ) then
      call vector_to_prog_var(statevec, PSURF_index, data_2d)
      where (data_2d == 0.0_r8) data_2d = NF90_FILL_REAL
      call nc_check(NF90_inq_varid(ncFileID, 'PSURF', VarID), &
                   'nc_write_model_vars', 'PSURF inq_varid '//trim(filename))
      call nc_check(nf90_put_var(ncFileID,VarID,data_2d,start=(/1,1,copyindex,timeindex/)),&
                   'nc_write_model_vars', 'PSURF put_var '//trim(filename))
   endif
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

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
real(r8)    :: random_number

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
      call get_state_kind_inc_dry(dart_index, var_type)
      do j=1, ens_size
         if (var_type /= KIND_DRY_LAND) then
            state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
               state_ens_handle%copies(j,i), &
               model_perturbation_amplitude)
   
         endif
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

      call get_state_kind_inc_dry(j, var_type)

      if(var_type /= KIND_DRY_LAND) then
         random_number = random_gaussian(r(sequence_to_use), 0.0_r8, model_perturbation_amplitude)
         call get_var_owner_index(j, owner, owners_index)
         if (ens_handle%my_pe==owner) then
            ens_handle%copies(i, owners_index) = ens_handle%copies(i, owners_index) + random_number
         endif
      endif

   enddo

enddo

end subroutine pert_model_copies_bitwise_lanai

!------------------------------------------------------------------

subroutine restart_file_to_sv(filename, state_vector, model_time)
 character(len=*), intent(in)    :: filename 
 real(r8),         intent(inout) :: state_vector(:)
 type(time_type),  intent(out)   :: model_time

! Reads the current time and state variables from a openggcm restart
! file and packs them into a dart state vector.

! temp space to hold data while we are reading it
real(r8) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)
integer  :: i, j, k, ivar, indx

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname 
integer :: VarID, numdims, dimlen
integer :: ncid, iyear, imonth, iday, ihour, iminute, isecond

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ... 
! Read the time data. 
! Note from Nancy Norton as pertains time:
! "The time recorded in the openggcm2 restart files is the current time,
! which corresponds to the time of the XXXX_CUR variables.
!
! current time is determined from iyear, imonth, iday, and *seconds_this_day*
!
! The ihour, iminute, and isecond variables are used for internal
! model counting purposes, but because isecond is rounded to the nearest
! integer, it is possible that using ihour,iminute,isecond information
! on the restart file to determine the exact curtime would give you a 
! slightly wrong answer."
!
! DART only knows about integer number of seconds, so using the rounded one
! is what we would have to do anyway ... and we already have a set_date routine
! that takes ihour, iminute, isecond information.

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'restart_file_to_sv', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  'restart_file_to_sv', 'get_att iyear')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  'restart_file_to_sv', 'get_att imonth')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  'restart_file_to_sv', 'get_att iday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  'restart_file_to_sv', 'get_att ihour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  'restart_file_to_sv', 'get_att iminute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  'restart_file_to_sv', 'get_att isecond')

! FIXME: we don't allow a real year of 0 - add one for now, but
! THIS MUST BE FIXED IN ANOTHER WAY!
if (iyear == 0) then
  call error_handler(E_MSG, 'restart_file_to_sv', &
                     'WARNING!!!   year 0 not supported; setting to year 1')
  iyear = 1
endif

model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

if (do_output()) &
    call print_time(model_time,'time for restart file '//trim(filename))
if (do_output()) &
    call print_date(model_time,'date for restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the 3d arrays into a single 1d list of numbers.
! These must be a fixed number and in a fixed order.

indx = 1

! Fill the state vector with the variables that are provided.
! The openggcm restart files have two time steps for each variable,
! the variables are named SALT_CUR and SALT_OLD ... for example.
! We are only interested in the CURrent time step.

do ivar=1, nfields

   varname = trim(variable_table(ivar, VAR_NAME_INDEX))
   numdims = get_num_dims(domain_id, ivar)
   string2 = trim(filename)//' '//trim(varname)

   SELECT CASE (numdims)
      CASE (3) ! 3D variable

         do i = 1,numdims
            write(msgstring,'(''inquire dimension'',i2,A)') i,trim(string2)
            call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                  'restart_file_to_sv', msgstring)

            if (dimlen /= size(data_3d_array,i)) then
               write(msgstring,*) trim(string2),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
               call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
            endif
         enddo   

         ! Actually get the variable and stuff it into the array

         call nc_check(nf90_get_var(ncid, VarID, data_3d_array), 'restart_file_to_sv', &
                      'get_var '//trim(varname))

         !>@todo JH: use reshape instead of manualy unfolding the vector
         do k = 1, Nz   ! size(data_3d_array,3)
         do j = 1, Ny   ! size(data_3d_array,2)
         do i = 1, Nx   ! size(data_3d_array,1)
            state_vector(indx) = data_3d_array(i, j, k)
            indx = indx + 1
         enddo
         enddo
         enddo

      CASE (2) ! 2D variable

         do i = 1,numdims
            write(msgstring,'(''inquire dimension'',i2,A)') i,trim(string2)
            call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                  'restart_file_to_sv', msgstring)

            if (dimlen /= size(data_2d_array,i)) then
               write(msgstring,*) trim(string2),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
               call error_handler(E_ERR,'restart_file_to_sv',msgstring,source,revision,revdate)
            endif
         enddo   

         ! Actually get the variable and stuff it into the array

         call nc_check(nf90_get_var(ncid, VarID, data_2d_array), 'restart_file_to_sv', &
                      'get_var '//trim(varname))

         !>@todo JH: use reshape instead of manualy unfolding the vector
         do j = 1, Ny   ! size(data_3d_array,2)
         do i = 1, Nx   ! size(data_3d_array,1)
            state_vector(indx) = data_2d_array(i, j)
            indx = indx + 1
         enddo
         enddo

      CASE DEFAULT

         write(msgstring,*) trim(string2),'numdims ',numdims,' not supported'
         call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)

   END SELECT
enddo

end subroutine restart_file_to_sv

!------------------------------------------------------------------

subroutine sv_to_restart_file(state_vector, filename, statedate)
 real(r8),         intent(in) :: state_vector(:)
 character(len=*), intent(in) :: filename 
 type(time_type),  intent(in) :: statedate

! Writes the current time and state variables from a dart state
! vector (1d fortran array) into a openggcm netcdf restart file.

integer :: iyear, imonth, iday, ihour, iminute, isecond
type(time_type) :: openggcm_time

! temp space to hold data while we are writing it
real(r8) :: data_2d_array(Nx,Ny), data_3d_array(Nx,Ny,Nz)

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname 

integer :: i, ivar, ncid, VarID, numdims, dimlen

!----------------------------------------------------------------------
! Get the show underway
!----------------------------------------------------------------------

if ( .not. module_initialized ) call static_init_model

! Check that the input file exists. 
! make sure the time tag in the restart file matches 
! the current time of the DART state ...

if ( .not. file_exist(filename)) then
   write(msgstring,*)trim(filename),' does not exist. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

call nc_check( nf90_open(trim(filename), NF90_WRITE, ncid), &
                  'sv_to_restart_file', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  'sv_to_restart_file', 'get_att iyear')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  'sv_to_restart_file', 'get_att imonth')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  'sv_to_restart_file', 'get_att iday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  'sv_to_restart_file', 'get_att ihour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  'sv_to_restart_file', 'get_att iminute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  'sv_to_restart_file', 'get_att isecond')

openggcm_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

if ( openggcm_time /= statedate ) then
   call print_time(statedate,'DART current time',logfileunit) 
   call print_time( openggcm_time,'openggcm  current time',logfileunit) 
   call print_time(statedate,'DART current time') 
   call print_time( openggcm_time,'openggcm  current time') 
   write(msgstring,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate) 
endif

if (do_output()) &
    call print_time(openggcm_time,'time of restart file '//trim(filename))
if (do_output()) &
    call print_date(openggcm_time,'date of restart file '//trim(filename))

! fill up the netcdf file with values
do ivar=1, nfields

   varname = trim(variable_table(ivar, VAR_NAME_INDEX))
   numdims = get_num_dims(domain_id, ivar)
   string2 = trim(filename)//' '//trim(varname)

   SELECT CASE (numdims)
      CASE (3) ! 3D variable

         do i = 1,numdims
            write(msgstring,'(''inquire dimension'',i2,A)') i,trim(string2)
            call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                  'sv_to_restart_file', msgstring)

            if (dimlen /= size(data_3d_array,i)) then
               write(msgstring,*) trim(string2),'dim/dimlen',i,dimlen,'not',size(data_3d_array,i)
               call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
            endif
         enddo

         call vector_to_prog_var(state_vector, ivar, data_3d_array)

         ! Actually stuff it into the netcdf file
         call nc_check(nf90_put_var(ncid, VarID, data_3d_array), &
                  'sv_to_restart_file', 'put_var '//trim(string2))

      CASE (2) ! 2D variable

         do i = 1,numdims
            write(msgstring,'(''inquire dimension'',i2,A)') i,trim(string2)
            call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                  'sv_to_restart_file', msgstring)
      
            if (dimlen /= size(data_2d_array,i)) then
               write(msgstring,*) trim(string2),'dim/dimlen',i,dimlen,'not',size(data_2d_array,i)
               call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)
            endif
         enddo
      
         call vector_to_prog_var(state_vector, ivar, data_2d_array)
      
         call nc_check(nf90_put_var(ncid, VarID, data_2d_array), &
                  'sv_to_restart_file', 'put_var '//trim(string2))
      
      CASE DEFAULT

         write(msgstring,*) trim(string2),'numdims ',numdims,' not supported'
         call error_handler(E_ERR,'sv_to_restart_file',msgstring,source,revision,revdate)

   END SELECT

enddo

call nc_check(nf90_close(ncid), 'sv_to_restart_file', 'close '//trim(filename))

end subroutine sv_to_restart_file

!------------------------------------------------------------------

subroutine get_gridsize(num_x, num_y, num_z)
 integer, intent(out) :: num_x, num_y, num_z

! public utility routine.

if ( .not. module_initialized ) call static_init_model

 num_x = Nx
 num_y = Ny
 num_z = Nz

end subroutine get_gridsize

!------------------------------------------------------------------

subroutine write_grid_netcdf()

! Write the grid to a netcdf file for checking.

integer :: ncid, NlonDimID, NlatDimID, NzDimID
integer :: nlon, nlat, nz
integer :: ulatVarID, ulonVarID, TLATvarid, TLONvarid
integer :: ZGvarid, ZCvarid, KMTvarid, KMUvarid

integer :: dimids(2);

if ( .not. module_initialized ) call static_init_model

nlon = size(ULAT,1)
nlat = size(ULAT,2)
nz   = size(ZG)

call nc_check(nf90_create('dart_grid.nc', NF90_CLOBBER, ncid),'write_grid_netcdf')

! define dimensions

call nc_check(nf90_def_dim(ncid, 'i', nlon, NlonDimID),'write_grid_netcdf')
call nc_check(nf90_def_dim(ncid, 'j', nlat, NlatDimID),'write_grid_netcdf')
call nc_check(nf90_def_dim(ncid, 'k',   nz,   NzDimID),'write_grid_netcdf')

dimids(1) = NlonDimID 
dimids(2) = NlatDimID 

! define variables

! FIXME: we should add attributes to say what units the grids are in (degrees).
call nc_check(nf90_def_var(ncid,  'KMT', nf90_int,     dimids,  KMTvarid),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid,  'KMU', nf90_int,     dimids,  KMUvarid),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid, 'ULON', nf90_double,  dimids, ulonVarID),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid, 'ULAT', nf90_double,  dimids, ulatVarID),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid, 'TLON', nf90_double,  dimids, TLONvarid),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid, 'TLAT', nf90_double,  dimids, TLATvarid),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid,   'ZG', nf90_double, NzDimID,   ZGvarid),'write_grid_netcdf')
call nc_check(nf90_def_var(ncid,   'ZC', nf90_double, NzDimID,   ZCvarid),'write_grid_netcdf')

call nc_check(nf90_put_att(ncid,ulonVarID,'long_name','U,V grid lons'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ulatVarID,'long_name','U,V grid lats'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,tlonVarID,'long_name','S,T grid lons'), &
                                                     'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,tlatVarID,'long_name','S,T grid lats'), &
                                                    'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ZCVarID,'long_name','vertical grid centers'), &
                                                    'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ZCVarID,'units','meters'), &
                                                    'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ZGVarID,'long_name','vertical grid bottoms'), &
                                                    'write_grid_netcdf')
call nc_check(nf90_put_att(ncid,ZGVarID,'units','meters'), &
                                                    'write_grid_netcdf')

call nc_check(nf90_enddef(ncid),'write_grid_netcdf')

! fill variables

call nc_check(nf90_put_var(ncid,  KMTvarid,  KMT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid,  KMUvarid,  KMU),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, ulatVarID, ULAT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, ulonVarID, ULON),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, TLATvarid, TLAT),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid, TLONvarid, TLON),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid,   ZGvarid,   ZG),'write_grid_netcdf')
call nc_check(nf90_put_var(ncid,   ZCvarid,   ZC),'write_grid_netcdf')

call nc_check(nf90_close(ncid),'write_grid_netcdf')

end subroutine write_grid_netcdf

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

! Loop over potentially close subset of obs priors or state variables
!if (present(dist)) then
do k = 1, num_close

   t_ind = close_ind(k)

   ! if dry land, leave original 1e9 value.  otherwise, compute real dist.
   if (obs_kind(t_ind) /= KIND_DRY_LAND) then
      dist(k) = get_dist(base_obs_loc,       obs(t_ind), &
                         base_obs_kind, obs_kind(t_ind))
   endif

enddo
!endif

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
offset = base_offset + (hgt_bot - 1) * nx * ny
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
offset = base_offset + (hgt_top - 1) * nx * ny
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

!------------------------------------------------------------------

subroutine insitu_temp(potemp, s, lpres, insitu_t)
 real(r8), intent(in)  :: potemp, s, lpres
 real(r8), intent(out) :: insitu_t

! CODE FROM openggcm MODEL -
! nsc 1 nov 2012:  i have taken the original subroutine with call:
!  subroutine dpotmp(press,temp,s,rp,potemp)
! and removed the original 'press' argument (setting it to 0.0 below)
! and renamed temp -> potemp, and potemp -> insitu_t
! i also reordered the args to be a bit more logical.  now you specify:
! potential temp, salinity, local pressure in decibars, and you get
! back in-situ temperature (called sensible temperature in the atmosphere;
! what a thermometer would measure).  the original (F77 fixed format) code
! had a computed goto which is deprecated/obsolete.  i replaced it with 
! a set of 'if() then else if()' lines.  i did try to not alter the original
! code so much it wasn't recognizable anymore.
!
!  aliciak note: rp = 0 and press = local pressure as function of depth
!  will return potemp given temp.
!  the trick here that if you make rp = local pressure and press = 0.0, 
!  and put potemp in the "temp" variable , it will return insitu temp in the 
!  potemp variable.

! an example figure of the relationship of potential temp and in-situ temp
! at depth:  http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_05.htm
! see the 'potential temperature' section (note graph starts at -1000m)

!     title:
!     *****

!       insitu_temp  -- calculate sensible (in-situ) temperature from 
!                       local pressure, salinity, and potential temperature

!     purpose:
!     *******

!       to calculate sensible temperature, taken from a converter that
!       went from sensible/insitu temperature to potential temperature
!
!       ref: N.P. Fofonoff
!            Deep Sea Research
!            in press Nov 1976

!     arguments:
!     **********

!       potemp     -> potential temperature in celsius degrees
!       s          -> salinity pss 78
!       lpres      -> local pressure in decibars
!       insitu_t   <- in-situ (sensible) temperature (deg c)


!     local variables:
!     *********

integer  :: i,j,n
real(r8) :: dp,p,q,r1,r2,r3,r4,r5,s1,t,x

!     code:
!     ****

      s1 = s - 35.0_r8
      p  = 0.0_r8
      t  = potemp

      dp = lpres - p
      n  = int (abs(dp)/1000.0_r8) + 1
      dp = dp/n

      do i=1,n
         do j=1,4

            r1 = ((-2.1687e-16_r8 * t + 1.8676e-14_r8) * t - 4.6206e-13_r8) * p
            r2 = (2.7759e-12_r8*t - 1.1351e-10_r8) * s1
            r3 = ((-5.4481e-14_r8 * t + 8.733e-12_r8) * t - 6.7795e-10_r8) * t
            r4 = (r1 + (r2 + r3 + 1.8741e-8_r8)) * p + (-4.2393e-8_r8 * t+1.8932e-6_r8) * s1
            r5 = r4 + ((6.6228e-10_r8 * t-6.836e-8_r8) * t + 8.5258e-6_r8) * t + 3.5803e-5_r8

            x  = dp*r5

            if (j == 1) then
               t = t + 0.5_r8 * x
               q = x
               p = p + 0.5_r8 * dp
          
            else if (j == 2) then
               t = t + 0.29298322_r8 * (x-q)
               q = 0.58578644_r8 * x + 0.121320344_r8 * q
   
            else if (j == 3) then
               t = t + 1.707106781_r8 * (x-q)
               q = 3.414213562_r8*x - 4.121320344_r8*q
               p = p + 0.5_r8*dp

            else ! j must == 4
               t = t + (x - 2.0_r8 * q) / 6.0_r8

            endif
   
         enddo ! j loop
      enddo ! i loop

      insitu_t = t

if (debug > 2) print *, 'potential temp, salinity, local pressure -> sensible temp'
if (debug > 2) print *, potemp, s, lpres, insitu_t

!       potemp     -> potential temperature in celsius degrees
!       s          -> salinity pss 78
!       lpres      -> local pressure in decibars
!       insitu_t   <- in-situ (sensible) temperature (deg c)

end subroutine insitu_temp

!------------------------------------------------------------------

subroutine dpth2pres(nd, depth, pressure)
 integer,  intent(in)  :: nd
 real(r8), intent(in)  :: depth(nd)
 real(r8), intent(out) :: pressure(nd)

!  description:
!  this function computes pressure in bars from depth in meters
!  using a mean density derived from depth-dependent global 
!  average temperatures and salinities from levitus 1994, and 
!  integrating using hydrostatic balance.
! 
!  references:
! 
!  levitus, s., r. burgett, and t.p. boyer, world ocean atlas 
!  volume 3: salinity, noaa atlas nesdis 3, us dept. of commerce, 1994.
! 
!  levitus, s. and t.p. boyer, world ocean atlas 1994, volume 4:
!  temperature, noaa atlas nesdis 4, us dept. of commerce, 1994.
! 
!  dukowicz, j. k., 2000: reduction of pressure and pressure
!  gradient errors in ocean simulations, j. phys. oceanogr., submitted.

!  input parameters:
!  nd     - size of arrays
!  depth  - depth in meters. no units check is made

!  output parameters:
!  pressure - pressure in bars 

!  local variables & parameters:
integer :: n
real(r8), parameter :: c1 = 1.0_r8

! -----------------------------------------------------------------------
!  convert depth in meters to pressure in bars
! -----------------------------------------------------------------------

      do n=1,nd
         pressure(n) = 0.059808_r8*(exp(-0.025_r8*depth(n)) - c1)  &
                     + 0.100766_r8*depth(n) + 2.28405e-7_r8*depth(n)**2
      end do

if (debug > 2 .and. do_output()) then
   print *, 'depth->pressure conversion table.  cols are: N, depth(m), pressure(bars)'
   do n=1,nd
      print *, n, depth(n), pressure(n)
   enddo
endif

end subroutine dpth2pres

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
write(construct_file_name_in, '(A, i4.4, A)') trim(stub), copy, ".nc"


end function construct_file_name_in

!--------------------------------------------------------------------
!> read the time from the input file
!> Stolen from openggcm model_mod.f90 restart_to_sv
function read_model_time(filename)

character(len=1024) :: filename
type(time_type) :: read_model_time


integer :: ncid !< netcdf file id
integer :: iyear, imonth, iday, ihour, iminute, isecond

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(msgstring,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',msgstring,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iyear'  , iyear), &
                  'read_model_time', 'get_att iyear')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'imonth' , imonth), &
                  'read_model_time', 'get_att imonth')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iday'   , iday), &
                  'read_model_time', 'get_att iday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'ihour'  , ihour), &
                  'read_model_time', 'get_att ihour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'iminute', iminute), &
                  'read_model_time', 'get_att iminute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'isecond', isecond), &
                  'read_model_time', 'get_att isecond')

! FIXME: we don't allow a real year of 0 - add one for now, but
! THIS MUST BE FIXED IN ANOTHER WAY!
if (iyear == 0) then
  call error_handler(E_MSG, 'read_model_time', &
                     'WARNING!!!   year 0 not supported; setting to year 1')
  iyear = 1
endif

read_model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)


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
   call use_default_state_variables( state_variables )
   string1 = 'model_nml:model_state_variables not specified using default variables'
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
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

! check to see if temp and salinity are both in the state otherwise you will not
! be able to interpolate in XXX subroutine
if ( any(kind_list == KIND_SALINITY) ) then
   ! check to see that temperature is also in the variable list
   if ( .not. any(kind_list == KIND_POTENTIAL_TEMPERATURE) ) then
      write(string1,'(A)') 'in order to compute temperature you need to have both '
      write(string2,'(A)') 'KIND_SALINITY and KIND_POTENTIAL_TEMPERATURE in the model state'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate, text2=string2)
   endif
endif
 
end subroutine verify_state_variables

!------------------------------------------------------------------
!> Default state_variables from the original openggcm model_mod.  Must
!> keep in the same order to be consistent with previous versions.
subroutine use_default_state_variables( state_variables )

character(len=*),  intent(inout) :: state_variables(:)

! strings must all be the same length for the gnu compiler
state_variables( 1:5*num_state_table_columns ) = &
   (/ 'SALT_CUR                  ', 'KIND_SALINITY             ', 'UPDATE                    ', &
      'TEMP_CUR                  ', 'KIND_POTENTIAL_TEMPERATURE', 'UPDATE                    ', &
      'UVEL_CUR                  ', 'KIND_U_CURRENT_COMPONENT  ', 'UPDATE                    ', &
      'VVEL_CUR                  ', 'KIND_V_CURRENT_COMPONENT  ', 'UPDATE                    ', &
      'PSURF_CUR                 ', 'KIND_SEA_SURFACE_PRESSURE ', 'UPDATE                    ' /)

end subroutine use_default_state_variables

!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
