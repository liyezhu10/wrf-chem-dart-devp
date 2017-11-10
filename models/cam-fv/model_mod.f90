! DART software - Copyright UCAR. This open source software is provided
! by ucar, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/dares/dart/dart_download
!
! $id: new_model_mod.f90 12058 2017-11-07 16:42:42z nancy@ucar.edu $
!----------------------------------------------------------------
!>
!> this is the interface between the cam-fv atmosphere model and dart.
!> the required public interfaces and arguments cannot be changed.
!>
!----------------------------------------------------------------

module model_mod

!>@todo fixme fill in the actual names we use after we've gotten
!>further into writing the coded

use             types_mod
use      time_manager_mod
use          location_mod, only : location_type, set_vertical, set_location, &
                                  get_location,get_close_obs, get_close_state, & 
                                  query_location, &
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, &
                                  VERTISSCALEHEIGHT
use         utilities_mod
use          obs_kind_mod
use     mpi_utilities_mod
use        random_seq_mod
use  ensemble_manager_mod
use distributed_state_mod
use   state_structure_mod
use  netcdf_utilities_mod,  only : nc_get_variable, nc_get_variable_size, &
                                   nc_add_attribute_to_variable, &
                                   nc_define_integer_variable, &
                                   nc_define_real_variable, &
                                   nc_add_global_creation_time, &
                                   nc_add_global_attribute, &
                                   nc_define_dimension, nc_put_variable, &
                                   nc_sync, nc_enddef, nc_redef, nc_open_readonly, &
                                   nc_close, nc_variable_exists
use       location_io_mod
use        quad_utils_mod
use     default_model_mod,  only : adv_1step, init_time, init_conditions, &
                                   nc_write_model_vars, pert_model_copies

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the dart code.

! routines in this list have code in this module
public :: static_init_model,                   &
          get_model_size,                      &
          get_state_meta_data,                 &
          model_interpolate,                   & ! big todo
          shortest_time_between_assimilations, &
          nc_write_model_atts,                 &
          write_model_time,                    & ! todo
          read_model_time,                     &
          end_model

! code for these routines are in other modules
public :: nc_write_model_vars,           &
          pert_model_copies,             & ! todo
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

! model_nml namelist variables and default values
character(len=256) :: cam_template_filename           = 'caminput.nc'
character(len=256) :: cam_phis_filename               = 'camphis.nc'
character(len=32)  :: vertical_localization_coord     = 'PRESSURE'
integer            :: assimilation_period_days        = 0
integer            :: assimilation_period_seconds     = 21600
integer            :: no_assim_above_this_model_level = 5
logical            :: use_damping_ramp_at_model_top   = .false.  
integer            :: debug_level                     = 0
logical            :: suppress_grid_info_in_output    = .false.

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
   debug_level


! global variables
character(len=512) :: string1, string2, string3
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
! All the REQUIRED interfaces come first - by convention.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
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
integer :: ncid
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

end subroutine static_init_model


!-----------------------------------------------------------------------
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.
!> Required for all applications.
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
integer  :: myvarid, myqty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=myvarid)

myqty = get_kind_index(domain_id, myvarid)

location = get_location_from_index(iloc, jloc, vloc, myvarid)

! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data

!-----------------------------------------------------------------------

function get_location_from_index(i, j, k, q)
integer, intent(in) :: i
integer, intent(in) :: j
integer, intent(in) :: k
integer, intent(in) :: q
type(location_type) :: get_location_from_index


select case (grid_stagger%qty_stagger(q))
  case (STAGGER_U)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%slat%vals(j), &
                                          real(k,r8), VERTISLEVEL)

  case (STAGGER_V)
   get_location_from_index = set_location(grid_data%slon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          real(k,r8), VERTISLEVEL)
   
  !>@todo not sure what to do yet. ? +-1/2 ?
  case (STAGGER_W)
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          real(k,r8)-0.5_r8, VERTISLEVEL)
  ! no stagger - cell centers
  case default
   get_location_from_index = set_location(grid_data%lon%vals(i), &
                                          grid_data%lat%vals(j), &
                                          real(k,r8), VERTISLEVEL)

end select

end function get_location_from_index

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
state_indx = get_dart_vector_index(lon_index, lat_index, lev_index, domain_id, varid)
vals(:)    = get_state(state_indx, ens_handle)

if (varid < 0) my_status = 12

end subroutine get_values_from_qty


!-----------------------------------------------------------------------
!>
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

!>@todo FIXME add error checking?  is state_indx < 0 how it indicates error?

!>@todo FIXME find unique level indices here (see new wrf code)

state_indx = get_dart_vector_index(lon_index, lat_index, lev_index(1), domain_id, varid)
vals       = get_state(state_indx, ens_handle)

my_status(:) = 0

end subroutine get_values_from_varid

!-----------------------------------------------------------------------
!> this is just for 3d fields

subroutine get_values_from_nonstate_fields(ens_handle, ens_size, lon_index, lat_index, lev_index, obs_quantity, vals, my_status)
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

vals(:) = MISSING_R8

select case (obs_quantity)
   case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT) 
      if (obs_quantity == QTY_PRESSURE) then
         call cam_pressure_levels(ens_handle, ens_size, &
                                  lon_index, lat_index, grid_data%lev%nsize, &
                                  vals_array, my_status)
      else
         call cam_height_levels(ens_handle, ens_size, &
                                lon_index, lat_index, grid_data%lev%nsize, &
                                vals_array, my_status) 
      endif

      if (any(my_status /= 0)) return

      do imember=1,ens_size
         vals(imember) = vals_array(lev_index(imember), imember)
      enddo

   case (QTY_VERTLEVEL)
      vals(:)      = lev_index(:)
      my_status(:) = 0

   case default
      print*, 'should not go here'

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

integer  :: varid
integer  :: lon_bot, lat_bot, lon_top, lat_top
real(r8) :: lon_fract, lat_fract
real(r8) :: lon_lat_vert(3), botvals(ens_size), topvals(ens_size)
integer  :: level_one_array(ens_size)
integer  :: which_vert, status1, status2, status_array(ens_size)
type(quad_interp_handle) :: interp_handle
integer  :: ijk(3), icorner, imember, numdims
integer  :: four_lons(4), four_lats(4)
integer  :: two_bots(2), two_tops(2)
real(r8) :: two_horiz_fracts(2)
integer  :: four_bot_levs(4, ens_size), four_top_levs(4, ens_size)
real(r8) :: four_vert_fracts(4, ens_size)
real(r8) :: quad_vals(4, ens_size)

if ( .not. module_initialized ) call static_init_model

! the overall strategy:
! 1. figure out if the quantity to interpolate is
! in the state.  if so, compute and return the value.
! if not, is it something that is a simple function of
! items in the state that we should compute here instead
! of computing it in a separate forward operator?  ok, we'll
! do it.. else error.
! (if there *are* functions of state vars or others that we
! need to compute here, put the rest of the interp code into
! a separate subroutine for reuse.)
! 2. compute the 4 horizontal i,j indices for the quad corners
! that enclose the obs location.
! 3. for each of the 4 quad corners compute the 2 vertical levels 
! that enclose the obs in the vertical. also return the fraction
! in the vertical - this needs a linear/log option.  (how set?)
! 4. now we have 8 i,j,k index numbers and 3 fractions.
! 5. compute the data values at each of the 4 horizontal i,j quad
! corners, interpolating in the vertical using the k fraction
! 6. compute the final horizontal obs value based on the i,j fractions


! Successful istatus is 0
interp_vals(:) = MISSING_R8
istatus(:)     = 99

call ok_to_interpolate(obs_qty, varid, status1)

! If not, for now return an error.  if there are other quantities
! that we can return that aren't part of the state then add code here.
! (make the rest of this routine into a separate subroutine that can
! be reused?)
if (status1 /= 0) then  
   !>@todo FIXME there may be things we need to compute that
   !> has multiple variables involved
   ! e.g. phis (elevation?)

   ! generally the interpolation code should not print or error out if
   ! it is asked to interpolate a quantity it cannot do.  this is then just
   ! a failed forward operator.  (some forward operators try one quantity
   ! and if that fails, then they compute what they need using different quantities.)
   if(debug_level > 12) then
      write(string1,*)'did not find observation quantity ', obs_qty, ' in the state vector'
      call error_handler(E_MSG,'model_interpolate:',string1,source,revision,revdate)
   endif
   istatus(:) = 2   ! this quantity not in the state vector
   return
endif

! unpack the location type into lon, lat, vert
! also may need the vert type?
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location))  ! default is to return the vertical type

! get the grid handle for the right staggered grid
interp_handle = get_interp_handle(obs_qty)

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2) , &
                         lon_bot, lat_bot, lon_top, lat_top, lon_fract, lat_fract, &
                         status1)

if (status1 /= 0) then
   istatus(:) = 3  ! cannot locate enclosing horizontal quad
   return
endif

!>@todo FIXME move this into quad_lon_lat_locate() interface so these are the
!>items directly returned.
! the order of the returned vertical info is counterclockwise around 
! the quad starting at lon/lat bot:
!  (lon_bot, lat_bot), (lon_top, lat_bot), (lon_top, lat_top), (lon_bot, lat_top)
! stuff this info into arrays of length 4 so we can loop over them easier.

call index_setup(lon_bot, lat_bot, lon_top, lat_top, lon_fract, lat_fract, &
                 four_lons, four_lats, two_horiz_fracts)

! need to consider the case for 2d vs 3d variables
numdims = get_dims_from_qty(obs_qty, varid)

!>@todo FIXME need to be refactored
! if variable is in the state
!    do 2d and 3d

! if variable is not in the state
!    do 2d and 3d

if (numdims > 2 ) then
   ! and now here potentially we have different results for different
   ! ensemble members.  the things that can vary are dimensioned by ens_size.
   do icorner=1, 4
      ! build a column to find vertical level numbers
      !>@doto FIXME this needs an option for linear or log scale.  affects fraction only.
      call find_vertical_levels(state_handle, ens_size, &
                                four_lons(icorner), four_lats(icorner), lon_lat_vert(3), &
                                which_vert, obs_qty, &
                                four_bot_levs(icorner, :), four_top_levs(icorner, :), &
                                four_vert_fracts(icorner, :), status_array)
      
      !>@todo FIXME should we let the process continue if at least one
      !>member has failed?  pro: save work  con: don't get forward operator
      !>values for members that could compute them
      !>(this is true for all the subsequent returns from this routine)
      if (any(status_array /= 0)) then
         istatus(:) = 4   ! cannot locate enclosing vertical levels  !>@todo FIXME use where statements?
         return
      endif
   enddo
   
   if (varid > 0) then
      ! we have all the indices and fractions we could ever want.
      ! now get the data values at the bottom levels, the top levels, 
      ! and do vertical interpolation to get the final result.
      
      do icorner=1, 4
         call get_values_from_varid(state_handle,  ens_size, &
                                    four_lons(icorner), four_lats(icorner), &
                                    four_bot_levs(icorner, :), varid, botvals, status_array)
         if (any(status_array /= 0)) then
            istatus(:) = 5   ! cannot retrieve bottom values
            return
         endif
         
         call get_values_from_varid(state_handle,  ens_size, &
                                    four_lons(icorner), four_lats(icorner), &
                                    four_top_levs(icorner, :), varid, topvals, status_array)
         if (any(status_array /= 0)) then
            istatus(:) = 6   ! cannot retrieve top values
            return
         endif

         call vert_interp(ens_size, botvals, topvals, four_vert_fracts(icorner, :), &
                          quad_vals(icorner, :), status_array)
  
         if (any(status_array /= 0)) then
            istatus(:) = 7   ! cannot do vertical interpolation
            return
         endif
      enddo
   else ! get 3d special variables in another ways ( like QTY_PRESSURE )
      ! fill topvals and botvals somehow
      do icorner=1, 4
         call get_values_from_nonstate_fields(state_handle,  ens_size, &
                                    four_lons(icorner), four_lats(icorner), &
                                    four_bot_levs(icorner, :), varid, botvals, status_array)
         if (any(status_array /= 0)) then
            istatus(:) = 5   ! cannot retrieve bottom values
            return
         endif

         call get_values_from_nonstate_fields(state_handle,  ens_size, &
                                    four_lons(icorner), four_lats(icorner), &
                                    four_top_levs(icorner, :), varid, topvals, status_array)
         if (any(status_array /= 0)) then
            istatus(:) = 6   ! cannot retrieve top values
            return
         endif

         call vert_interp(ens_size, botvals, topvals, four_vert_fracts(icorner, :), &
                          quad_vals(icorner, :), status_array)
      
         if (any(status_array /= 0)) then
            istatus(:) = 7   ! cannot do vertical interpolation
            return
         endif
      enddo
   endif
else ! 2 dimensional variable
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
                               obs_qty, quad_vals(icorner,:), status1)
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
!>  

function get_dims_from_qty(obs_quantity, var_id)
integer, intent(in) :: obs_quantity
integer, intent(in) :: var_id
integer :: get_dims_from_qty

if (var_id > 0) then
   get_dims_from_qty = get_num_dims(domain_id,var_id)
else
   select case (obs_quantity)
      case (QTY_SURFACE_ELEVATION)
         get_dims_from_qty = 2
      case default 
         write(string1, *) 'we can not interpolate qty', obs_quantity, &
                           ' if the dimension is not known'
         call error_handler(E_ERR,'get_dims_from_qty', &
                            string1,source,revision,revdate)
    end select
endif

end function get_dims_from_qty

!-----------------------------------------------------------------------
!>
!>  

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
   case (QTY_SURFACE_ELEVATION, QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT, QTY_VERTLEVEL)
      my_status = 0
   case default
      my_status = 13
end select


end subroutine ok_to_interpolate


!-----------------------------------------------------------------------
!>
!>  This is for 2d special observations quantities not in the state

subroutine get_quad_corners(ens_size, lon_index, lat_index, obs_quantity, vals, my_status)
integer,  intent(in) :: ens_size
integer,  intent(in) :: lon_index
integer,  intent(in) :: lat_index
integer,  intent(in) :: obs_quantity
real(r8), intent(out) :: vals(ens_size) 
integer,  intent(out) :: my_status

character(len=*), parameter :: routine = 'get_quad_corners'

select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)
      vals(:) = phis(lon_index, lat_index) / gravity
   case default 
      write(string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
end select

my_status = 0

end subroutine get_quad_corners


!-----------------------------------------------------------------------
!>
!>  interpolating first in the horizontal then the vertical

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
!>
!> 

! populate arrays with the 4 combinations of bot/top vs lon/lat
! to make it easier below to loop over the corners.
! no rocket science here.

subroutine index_setup(lon_bot, lat_bot, lon_top, lat_top, lon_fract, lat_fract, &
                       four_lons, four_lats, two_horiz_fracts)
integer,  intent(in)  :: lon_bot, lat_bot, lon_top, lat_top
real(r8), intent(in)  :: lon_fract, lat_fract
integer,  intent(out) :: four_lons(4), four_lats(4)
real(r8), intent(out) :: two_horiz_fracts(2)

! order is counterclockwise around the quad:
!  (lon_bot, lat_bot), (lon_top, lat_bot), (lon_top, lat_top), (lon_bot, lat_top)

four_lons(1) = lon_bot
four_lons(2) = lon_top
four_lons(3) = lon_bot
four_lons(4) = lon_top

four_lats(1) = lat_bot
four_lats(2) = lat_bot
four_lats(3) = lat_top
four_lats(4) = lat_top

two_horiz_fracts(1) = lon_fract
two_horiz_fracts(2) = lat_fract

end subroutine index_setup

!-----------------------------------------------------------------------
!>
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

integer :: bot1, top1, imember, nlevels, varid
integer :: level_one, status1
integer(i8) :: state_indx
real(r8) :: fract1
real(r8) :: surf_pressure(ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)

! assume the worst
bot_levs(:) = MISSING_I
top_levs(:) = MISSING_I
vert_fracts(:) = MISSING_R8
my_status(:) = 98

!>@todo FIXME we need nlevels everywhere - should we have a global var for it?

! number of vertical levels (midlayer points)
nlevels = grid_data%lev%nsize
level_one = 1

select case (which_vert)

   case(VERTISPRESSURE)
      ! construct a pressure column here and find the model levels
      ! that enclose this value
      call get_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                               lon_index, lat_index, level_one, &
                               surf_pressure, status1)

      !>@todo FIXME: should we figure out now or later? how many unique levels we have?
      !> for now - do the unique culling later so we don't have to carry that count around.
      do imember=1, ens_size
         call build_cam_pressure_column(surf_pressure(imember), nlevels, &
                                        grid_data%hyam, grid_data%hybm, &
                                        grid_data%P0, pressure_array)

         call pressure_to_level(nlevels, pressure_array, vert_val, &
                                bot_levs(imember), top_levs(imember), &
                                vert_fracts(imember), my_status(imember))

      enddo

   case(VERTISHEIGHT)
      ! construct a height column here and find the model levels
      ! that enclose this value
      !@>todo put in arguments and write height_to_level
      call cam_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, &
                             height_array, my_status)
      if (any(my_status /= 0)) return   !>@todo FIXME let successful members continue?
      do imember=1, ens_size
         !>@todo FIXME somewhere cull out unique levels and only get_state() for those. (see wrf)
         call height_to_level(nlevels, height_array(:, imember), vert_val, &
                              bot_levs(imember), top_levs(imember), vert_fracts(imember), my_status(imember))
      enddo
      

      write(string1, *) 'we have not written the code yet for vertical type: ', which_vert
      call error_handler(E_ERR,'find_vertical_levels', &
                         string1,source,revision,revdate)

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
      call error_handler(E_ERR,'find_vertical_levels', &
                         string1,source,revision,revdate)
      
end select

end subroutine find_vertical_levels

!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

subroutine cam_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, height_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: height_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: k, varid, level_one, imember, status1
real(r8)    :: temperature(ens_size), specific_humidity(ens_size), surface_pressure(ens_size)
real(r8)    :: phi(ens_size, nlevels)
real(r8)    :: surface_elevation
real(r8)    :: tv(ens_size, nlevels+1)  !>@todo FIXME  ????           ! Virtual temperature, top to bottom
integer(i8) :: state_indx

!>@todo this should come from a model specific constant module.
!> the forward operators and model_mod should use it.
real(r8), parameter :: rd = 287.05_r8 ! dry air gas constant
real(r8), parameter :: rv = 461.51_r8 ! wet air gas constant
real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

!>@todo FIXME: these need to be replaced by hyam, hymb, hyai, hybi -- VERY VERY CAREFULLY
real(r8) ::hybrid_As(nlevels+1,2), hybrid_Bs(nlevels+1,2)

! this is for surface obs
level_one = 1

!@>todo make into a subroutine get_val or something similar
! get the surface pressure from the ens_handle
call get_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                         lon_index, lat_index, level_one, surface_pressure, status1)


! get the surface elevation from the phis
surface_elevation = phis(lon_index, lat_index)

do k = 1, nlevels
   ! temperature
   call get_values_from_qty(ens_handle, ens_size, QTY_TEMPERATURE, &
                                     lon_index, lat_index, k, temperature, status1)

   ! specific humidity
   call get_values_from_qty(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, &
                                     lon_index, lat_index, k, specific_humidity, status1)
   
   !>@todo rename tv to something that mens something to users
   tv(:,k) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))
enddo

! need to convert to geopotential height
do imember = 1, ens_size
   !>@todo refacfor to just put out geometric height
   call dcz2(nlevels, surface_pressure(imember), surface_elevation, tv(imember,:), &
             grid_data%P0%vals(1), hybrid_As, hybrid_Bs, phi(imember,:))
   do k = 1,nlevels
      height_array(k, imember) = gph2gmh(phi(imember,k), grid_data%lat%vals(lat_index))
   enddo
enddo

my_status(:) = 0

end subroutine cam_height_levels

!-----------------------------------------------------------------------
!> Compute the pressures at the layer midpoints

subroutine build_cam_pressure_column(surface_pressure, n_levels, hyam, hybm, P0, pressure_array)

real(r8),           intent(in)  :: surface_pressure   ! in pascals
integer,            intent(in)  :: n_levels
type(cam_1d_array), intent(in)  :: hyam
type(cam_1d_array), intent(in)  :: hybm
type(cam_1d_array), intent(in)  :: P0
real(r8),           intent(out) :: pressure_array(:)

integer :: k

! Set midpoint pressures.  This array mirrors the order of the
! cam model levels: 1 is the model top, N is the bottom.

do k=1,n_levels
   pressure_array(k) = hyam%vals(k)*P0%vals(1) + hybm%vals(k)*surface_pressure
enddo

end subroutine build_cam_pressure_column


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

integer :: this_lev, varid

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
   fract = (p_val - pressures(top_lev)) / (pressures(bot_lev) - pressures(top_lev))
   my_status = 0
   return
enddo levloop

! you shouldn't get here
if (bot_lev == MISSING_I) then
  write(string1,*) 'should not happen - contact dart support'
  write(string2,*) 'pressure value ', p_val, ' was not found in pressure column'
  call error_handler(E_ERR,'pressure_to_level', &
                         string1,source,revision,revdate,text2=string2)
endif

end subroutine pressure_to_level

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> 1 is top, N is bottom, bot is the lower level, top is the upper level
!> so top value will be larger than bot.  fract = 0 is the full top,
!> fract = 1 is the full bot.  return non-zero if value outside the valid range.

subroutine height_to_level(nlevels, heights, h_val, &
                           bot_lev, top_lev, fract, my_status)

integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: heights(:)
real(r8), intent(in)  :: h_val
integer,  intent(out) :: bot_lev
integer,  intent(out) :: top_lev
real(r8), intent(out) :: fract  
integer,  intent(out) :: my_status

integer :: this_lev, varid

bot_lev = MISSING_I
top_lev = MISSING_I
fract   = MISSING_R8

if (h_val < heights(1) .or. h_val > heights(nlevels)) then
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
  write(string1,*) 'should not happen - contact dart support'
  write(string2,*) 'height value ', h_val, ' was not found in height column'
  call error_handler(E_ERR,'height_to_level', &
                         string1,source,revision,revdate,text2=string2)
endif

end subroutine height_to_level

!-----------------------------------------------------------------------
!> in cam level 1 is at the model top, level N is the lowest level
!> our convention in this code is:  between levels a fraction of 0
!> is 100% the top level, and fraction of 1 is 100% the bottom level.
!> the top level is always closer to the model top and so has a *smaller*
!> level number than the bottom level. stay alert for this!

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
!>  Next, get the pressures on the levels for this ps
!>  Assuming we'll only need pressures on model mid-point levels, not 
!>  interface levels.  This pressure column will be for the correct grid 
!> for obs_qty, since p_surf was taken
!>      from the grid-correct ps[_stagr] grid

subroutine find_pressure_levels()
real(r8), allocatable :: p_col(:,:)

end subroutine find_pressure_levels


!-----------------------------------------------------------------------
!>
!> 


function get_interp_handle(obs_quantity)
integer, intent(in)     :: obs_quantity
type(quad_interp_handle) :: get_interp_handle

select case (grid_stagger%qty_stagger(obs_quantity))
   case ( STAGGER_U ) 
      get_interp_handle = interp_u_staggered
   case ( STAGGER_V ) 
      get_interp_handle = interp_v_staggered
   case ( STAGGER_NONE )
      get_interp_handle = interp_nonstaggered
   case ( STAGGER_W ) 
      write(string1,*) 'w stagger -- not supported yet'
      call error_handler(E_ERR,'get_interp_handle', &
                         string1,source,revision,revdate)
   case ( STAGGER_UV ) 
      write(string1,*) 'uv stagger -- not supported yet'
      call error_handler(E_ERR,'get_interp_handle', &
                         string1,source,revision,revdate)
   case default
      write(string1,*) 'unknown stagger -- this should never happen'
      call error_handler(E_ERR,'get_interp_handle', &
                         string1,source,revision,revdate)
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

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = set_time(assimilation_period_seconds, &
                                               assimilation_period_days)

write(string1,*)'assimilation period is ',assimilation_period_days,   ' days ', &
                                          assimilation_period_seconds,' seconds'

call error_handler(E_MSG,'shortest_time_between_assimilations:', &
                   string1,source,revision,revdate)

end function shortest_time_between_assimilations




!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

! deallocate arrays from grid and anything else

call free_cam_grid(grid_data)

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
integer, intent(in) :: dom_id

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NzDimID
integer :: ulonVarID, ulatVarID, tlonVarID, tlatVarID, ZGVarID, ZCVarID
integer :: KMTVarID, KMUVarID

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

! Gaussian Wave
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

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
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
integer :: ios, VarID

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
! This code require that the time length be size 1

call nc_get_variable_size(ncid, 'time', timesize)

!>@todo do we really need to ceck this if it is never going to happen.
!#! if (timesize /= 1) then
!#!    write(string1,*) trim(filename),' has',timesize,'times. Require exactly 1.'
!#!    call error_handler(E_ERR, 'read_model_time', string1, source, revision, revdate)
!#! endif

call nc_get_variable(ncid, 'date',    cam_date)
call nc_get_variable(ncid, 'datesec', cam_tod)

! The 'date' is YYYYMMDD ... cam_tod is 'current seconds of current day'
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
      call error_handler(E_ERR,'set_cam_variable_info:',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(3A)') 'there is no obs_kind <', trim(dartstr), '> in obs_kind_mod.f90'
      call error_handler(E_ERR,'set_cam_variable_info:',string1,source,revision,revdate)
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

   call error_handler(E_MSG,'set_cam_variable_info:',string1, &
                      source,revision,revdate,text2=string2)
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
!>

subroutine fill_cam_stagger_info(stagger)
type(cam_stagger), intent(inout) :: stagger

integer :: ivar, jdim, qty_index

allocate(stagger%qty_stagger(0:get_num_quantities()))

stagger%qty_stagger = STAGGER_NONE

do ivar = 1, get_num_variables(domain_id)
   do jdim = 1, get_num_dims(domain_id, ivar)

      if (get_dim_name(domain_id, ivar, jdim) == 'slat') then
         qty_index = get_kind_index(domain_id, ivar) ! get_kind_index actualy returns the qty
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
!>
!> Given a generic CAM restart file, that contains grid information and populate
!> the grid information


subroutine read_grid_info(grid_file, grid)
character(len=*), intent(in)  :: grid_file
type(cam_grid),   intent(out) :: grid

character(len=*), parameter :: routine = 'read_grid_info'

integer :: ncid

! put this in a subroutine that deals with the grid
ncid = nc_open_readonly(grid_file, routine)

! Get the grid info plus additional non-state arrays
call get_cam_grid(ncid, grid)
call read_cam_phis_array(cam_phis_filename)

! Set up the interpolation structure for later 
call setup_interpolation(grid)

call nc_close(ncid, routine)

end subroutine read_grid_info


!-----------------------------------------------------------------------
!>
!> fill in the various cam grid arrays 
!> 

subroutine get_cam_grid(ncid, grid)
integer,        intent(in)  :: ncid
type(cam_grid), intent(out) :: grid

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
     call error_handler(E_ERR,'set_vert_localization:',string1, &
                        source,revision,revdate,text2=string2)
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

subroutine dcz2(kmax,p_surf,h_surf,tv,hprb,hybrid_As,hybrid_Bs,z2)

! Compute geopotential height for a CESM hybrid coordinate column.
! All arrays except hybrid_As, hybrid_Bs are oriented top to bottom.
! hybrid_[AB]s first subscript:
!  = 1 for layer interfaces
!  = 2 for layer midpoints
! hybrid_As coord coeffs for P0 reference pressure term in plevs_cam
! hybrid_Bs coord coeffs for surf pressure term in plevs_cam (in same format as hybrid_As)

integer,  intent(in)  :: kmax                ! Number of vertical levels
real(r8), intent(in)  :: p_surf              ! Surface pressure           (pascals)
real(r8), intent(in)  :: h_surf               ! Surface height (m)
real(r8), intent(in)  :: tv(kmax)            ! Virtual temperature, top to bottom
real(r8), intent(in)  :: hprb                ! Hybrid base pressure       (pascals)
real(r8), intent(in)  :: hybrid_As(kmax+1,2)
real(r8), intent(in)  :: hybrid_Bs(kmax+1,2)
real(r8), intent(out) :: z2(kmax)            ! Geopotential height, top to bottom

! Local variables

!>@todo have a model constants module?
real(r8), parameter :: r = 287.04_r8    ! Different than model_heights ! dry air gas constant.
real(r8), parameter :: g0 = 9.80616_r8  ! Different than model_heights:gph2gmh:G !
real(r8), parameter :: rbyg=r/g0
real(r8) :: pterm(kmax)         ! pressure profile
real(r8) :: pmln(kmax+1)        ! logs of midpoint pressures

integer  :: i,k,l
real(r8) :: ARG

! Compute intermediate quantities using scratch space

! DEBUG: z2 was unassigned in previous code.
z2(:) = MISSING_R8

! Invert vertical loop
! Compute top only if top interface pressure is nonzero.
!
! newFIXME; p_col could be used here, instead of (re)calculating it in ARG
do K = kmax+1, 1, -1
   i = kmax-K+2
   ARG = hprb*hybrid_As(i,2) + p_surf *hybrid_Bs(i,2)
   if (ARG > 0.0_r8) THEN
       pmln(K) = LOG(ARG)
   else
       pmln(K) = 0.0_r8
   endif
enddo

do K = 2,kmax - 1
   pterm(k) = rbyg*tv(k)*0.5_r8* (pmln(k+1)-pmln(k-1))
enddo

! Initialize z2 to sum of ground height and thickness of top half layer
! this is NOT adding the thickness of the 'top' half layer.
!    it's adding the thickness of the half layer at level K,
do K = 1,kmax - 1
   z2(k) = h_surf + rbyg*tv(k)*0.5_r8* (pmln(K+1)-pmln(K))
enddo
z2(kmax) = h_surf + rbyg*tv(kmax)* (log(p_surf*hybrid_Bs(1,1))-pmln(kmax))

! THIS is adding the half layer at the BOTTOM.
do k = 1,kmax - 1
    z2(k) = z2(k) + rbyg*tv(kmax)* (log(p_surf*hybrid_Bs(1,1))-0.5_r8* &
                                       (pmln(kmax-1)+pmln(kmax)))
enddo

! Add thickness of the remaining full layers
! (i.e., integrate from ground to highest layer interface)

do K = 1,kmax - 2
    do L = K+1, kmax-1
       z2(K) = z2(K) + pterm(L)
    enddo
enddo

end subroutine dcz2

!-----------------------------------------------------------------------

function gph2gmh(h, lat)

!  Convert a list of geopotential altitudes to mean sea level altitude.

real(r8), intent(in) :: h         ! geopotential altitude (in m)
real(r8), intent(in) :: lat       ! latitude  of profile in degrees.
real(r8)             :: gph2gmh   ! MSL altitude, in km.

real(r8), parameter ::  be = 6356751.6_r8        ! min earth radius, m
real(r8), parameter ::  ae = 6378136.3_r8        ! max earth radius, m
real(r8), parameter ::  pi = 3.14159265358979_r8

! FIXME; another definition of gravitational acceleration.  
! See g0 and gravity_constant elsewhere.
real(r8), parameter ::  G = 9.80665_r8 ! WMO reference g value, m/s**2, at 45.542N(S)

real(r8) :: g0
real(r8) :: r0
real(r8) :: latr

latr = lat * (pi/180.0_r8)           ! in radians
call compute_gravity(latr, 0.0_r8, g0)

! compute local earth's radius using ellipse equation

r0 = sqrt( ae**2 * cos(latr)**2 + be**2 * sin(latr)**2)

! Compute altitude above sea level
gph2gmh = (r0 * h) / (((g0*r0)/G) - h)

end function gph2gmh

!-----------------------------------------------------------------------

!> This subroutine computes the Earth's gravity at any altitude
!> and latitude.  The model assumes the Earth is an oblate
!> spheriod rotating at a the Earth's spin rate.  The model
!> was taken from "Geophysical Geodesy, Kurt Lambeck, 1988".
!>
!>  input:    xlat, latitude in radians
!>            alt,  altitude above the reference ellipsiod, km
!>  output:   galt, gravity at the given lat and alt, m/sec**2
!>
!> Compute acceleration due to the Earth's gravity at any latitude/altitude
!> author     Bill Schreiner   5/95
!>
!> changed from using kilometers to meters since that is our native unit.
!>

subroutine compute_gravity(xlat,alt,galt)

real(r8), intent(in)  :: xlat
real(r8), intent(in)  :: alt
real(r8), intent(out) :: galt

real(r8),parameter :: xmu = 398600.4415_r8         ! km^3/s^2
real(r8),parameter :: ae  = 6378.1363_r8           ! km
real(r8),parameter :: f   = 1.0_r8/298.2564_r8
real(r8),parameter :: w   = 7.292115e-05_r8        ! rad/s
real(r8),parameter :: xm  = 0.003468_r8            !
real(r8),parameter :: f2  = 5.3481622134089e-03_r8 ! f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
real(r8),parameter :: f4  = 2.3448248012911e-05_r8 ! f4 = -f**2* 0.50 + 5.0* 0.50*f*xm

real(r8) :: ge
real(r8) :: g


! compute gravity at the equator, km/s2
ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0_r8/14.0_r8*xm*f)

! compute gravity at any latitude, km/s2
g = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)

! compute gravity at any latitude and at any height, km/s2
galt = g - 2.0_r8*ge*alt/ae*(1.0_r8 + f + xm + (-3.0_r8*f + 5.0_r8* 0.50_r8*xm)*  &
                          (sin(xlat))**2) + 3.0_r8*ge*alt**2/ae**2
! convert to meters/s2
galt = galt*1000.0_r8

end subroutine compute_gravity

!-----------------------------------------------------------------------

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:)
integer(i8),         intent(in)    :: loc_indx(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: current_vert_type, ens_size, i

ens_size = 1

do i=1,num
   current_vert_type = nint(query_location(locs(i)))

   if ( current_vert_type == which_vert ) cycle

   select case (which_vert)
      case (VERTISPRESSURE)
         call convert_to_pressure(ens_handle, ens_size, locs(i), loc_indx(i))
      case (VERTISHEIGHT)
         call convert_to_height(ens_handle, ens_size, locs(i), loc_indx(i))
      case (VERTISLEVEL)
         call convert_to_level(ens_handle, ens_size, locs(i), loc_indx(i))
      case (VERTISSCALEHEIGHT)
         call convert_to_scaleheight(ens_handle, ens_size, locs(i), loc_indx(i))
      case default
         print*, ' can not convert vert'
   end select
enddo

istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

subroutine  convert_to_pressure(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

real(r8) :: pressure

!>@todo FIXME move up to convert_vertical_state ?
call state_level_to_pressure(ens_handle, ens_size, location_indx, pressure)

call set_vertical(location, pressure, VERTISPRESSURE)

end subroutine  convert_to_pressure

!--------------------------------------------------------------------

subroutine state_level_to_pressure(ens_handle, ens_size, location_indx, pressure)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
integer(i8),         intent(in)    :: location_indx
real(r8),            intent(out)   :: pressure

integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)

! build a height column and a pressure column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                         pressure_array, my_status)

pressure = pressure_array(vloc)

end subroutine state_level_to_pressure


!--------------------------------------------------------------------

subroutine convert_to_height(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

!>@todo FIXME move up to convert_vertical_state ?
call state_level_to_height(ens_handle, ens_size, location, location_indx)

end subroutine  convert_to_height

!--------------------------------------------------------------------

subroutine state_level_to_height(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx


integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)

! build a height column and a pressure column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_height_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       height_array, my_status) 

!>@todo FIXME this can only be used if ensemble size is 1
call set_vertical(location, height_array(vloc,1), VERTISHEIGHT)

end subroutine state_level_to_height

!--------------------------------------------------------------------

subroutine  convert_to_scaleheight(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

call state_level_to_scaleheight(ens_handle, ens_size, location, location_indx)

end subroutine  convert_to_scaleheight

!--------------------------------------------------------------------

subroutine  convert_to_level(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

call state_level_to_level(ens_handle, ens_size, location, location_indx)

end subroutine  convert_to_level

!--------------------------------------------------------------------

subroutine state_level_to_level(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)
real(r8) :: level(grid_data%lev%nsize)

! build a height column and a level column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call set_vertical(location, real(vloc,r8), VERTISLEVEL)

end subroutine state_level_to_level



!--------------------------------------------------------------------

subroutine state_level_to_scaleheight(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx

integer  :: iloc, jloc, vloc, level_one, status1, my_status(ens_size)
real(r8) :: height_array(grid_data%lev%nsize, ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize)
real(r8) :: surface_pressure(1), scaleheight

level_one = 1

! build a height column and a scaleheight column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

! get the surface pressure from the ens_handle
call get_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                         iloc, jloc, level_one, surface_pressure, status1)

call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                         pressure_array, my_status)

scaleheight = log(pressure_array(vloc)/surface_pressure(1))

call set_vertical(location, scaleheight, VERTISSCALEHEIGHT)

end subroutine state_level_to_scaleheight

!--------------------------------------------------------------------
!> this only works for ens_size = 1

subroutine obs_height_to_pressure(ens_handle, ens_size, location, location_indx)
type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: ens_size
type(location_type), intent(inout) :: location
integer(i8),         intent(in)    :: location_indx


integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(  grid_data%lev%nsize, ens_size)
real(r8) :: pressure_array(grid_data%lev%nsize, ens_size)

! build a height column and a pressure column and find the levels?
call get_model_variable_indices(location_indx, iloc, jloc, vloc)

call cam_height_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       height_array, my_status) 

call cam_pressure_levels(ens_handle, ens_size, iloc, jloc, grid_data%lev%nsize, &
                       pressure_array, my_status)

call set_vertical(location, pressure_array(vloc,1), VERTISPRESSURE)

end subroutine obs_height_to_pressure


!-----------------------------------------------------------------------
!> Compute the pressure values at midpoint levels
!>
!> this version does all ensemble members at once.

subroutine cam_pressure_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, pressure_array, my_status) 
type(ensemble_type), intent(in)  :: ens_handle
integer,             intent(in)  :: ens_size
integer,             intent(in)  :: lon_index
integer,             intent(in)  :: lat_index
integer,             intent(in)  :: nlevels
real(r8),            intent(out) :: pressure_array(nlevels, ens_size)
integer,             intent(out) :: my_status(ens_size)

integer     :: level_one, status1, imember
real(r8)    :: surface_pressure(ens_size)

! this is for surface obs
level_one = 1

! get the surface pressure from the ens_handle
call get_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                         lon_index, lat_index, level_one, surface_pressure, status1)

do imember=1, ens_size
   call build_cam_pressure_column(surface_pressure(imember), grid_data%lev%nsize, &
                                  grid_data%hyam, grid_data%hybm, grid_data%P0, &
                                  pressure_array(:,imember))
enddo
my_status(:) = 0

end subroutine cam_pressure_levels

!--------------------------------------------------------------------

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, my_status)

type(ensemble_type), intent(in)    :: ens_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(:)
integer,             intent(in)    :: loc_qtys(:), loc_types(:)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: my_status(:)

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
         print*, ' can not convert vert'
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

integer  :: iloc, jloc, vloc, my_status(ens_size)
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

integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: height_array(ens_size)

call model_interpolate(ens_handle, ens_size, location, QTY_GEOMETRIC_HEIGHT, height_array(:), my_status(:))

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

integer  :: iloc, jloc, vloc, my_status(ens_size)
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

integer  :: iloc, jloc, vloc, my_status(ens_size)
real(r8) :: pressure_array(ens_size)
real(r8) :: surface_pressure_array(ens_size)
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
