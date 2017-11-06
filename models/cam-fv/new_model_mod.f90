! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
!----------------------------------------------------------------
!>
!> This is the interface between the CAM-FV atmosphere model and DART.
!> The required public interfaces and arguments CANNOT be changed.
!>
!----------------------------------------------------------------

module model_mod

!>@todo FIXME fill in the actual names we use after we've gotten
!>further into writing the coded

use             types_mod
use      time_manager_mod
use          location_mod
use         utilities_mod
use          obs_kind_mod
use     mpi_utilities_mod
use        random_seq_mod
use  ensemble_manager_mod
use distributed_state_mod
use   state_structure_mod
use  netcdf_utilities_mod,  only : nc_check, nc_get_variable, nc_get_variable_size, &
                                   nc_add_attribute_to_variable, nc_define_real_variable, &
                                   nc_add_global_creation_time, nc_add_global_attribute, &
                                   nc_define_dimension, nc_put_variable, nc_sync, nc_enddef, &
                                   nc_redef
use       location_io_mod
use        quad_utils_mod
use     default_model_mod,  only : adv_1step, init_time, init_conditions, &
                                   nc_write_model_vars, pert_model_copies
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: static_init_model,             &
          get_model_size,                &
          get_state_meta_data,           &
          model_interpolate,             &
          shortest_time_between_assimilations, &
          convert_vertical_obs,          &
          convert_vertical_state,        &
          nc_write_model_atts,           &
          write_model_time,              &
          read_model_time,               &
          end_model

! code for these routines are in other modules
public :: nc_write_model_vars,           &
          pert_model_copies,             &
          adv_1step,                     &
          init_time,                     &
          init_conditions,               &
          get_close_obs,                 &
          get_close_state

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! global variables
character(len=512) :: string1, string2, string3
logical, save      :: module_initialized = .false.

! model_nml namelist variables
integer  :: assimilation_period_days     = 0
integer  :: assimilation_period_seconds  = 21600
integer  :: vert_localization_coord      = VERTISPRESSURE
integer  :: debug = 0   ! turn up for more and more debug messages
logical  :: minimal_output = .false.
character(len=256) :: cam_template_filename = 'caminput.nc' ! template restart file
!#! character(len=256) :: cam_grid_filename     = 'caminput.nc' ! JPH probably delete?
character(len=256) :: cam_phis_filename     = 'camphis.nc'

namelist /model_nml/  &
!#!    cam_grid_filename,           &
   assimilation_period_days,    &
   assimilation_period_seconds, &
   cam_template_filename,       &
   cam_phis_filename,           &
   vert_localization_coord,     &
   minimal_output,              &
   debug,                       &
   state_variables

! the number of state variables allowed in the state.  this may need to change
! depending on the number of variables that you are using.
integer, parameter :: MAX_STATE_VARIABLES = 100
! column variables should be the same as the cam_template_file and
! have the same shape as the restarts being used. For no clamping use 'NA'
!
!    variable_name, variable_type, clamp_min, clamp_max, update_variable
integer, parameter :: num_state_table_columns = 5
!>@todo JPH what should vtablenamelength be?
character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLES * &
                                                   num_state_table_columns ) = ' '

! domain id for the cam model.  this allows us access to all of the state structure
! info and is require for getting state variables.
integer :: domain_id

!> Everything needed to describe a variable. Basically all the metadata from
!> a netCDF file is stored here as well as all the information about where
!> the variable is stored in the DART state vector.
!>

type cam_1d_array
   integer  :: nsize
   real(r8), allocatable :: vals(:)
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
call read_cam_phis_array(cam_phis_filename)

!>@todo do we need to map_qtys here?
!>@todo do we need to set the model top related stuff here?

! set_cam_variable_info() fills var_names, kind_list, clamp_vals, update_list
! from the &model_mod_nml state_variables
call set_cam_variable_info(state_variables, nfields)

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
!> Model interpolate will interpolate any DART state variable
!> (i.e. S, T, U, V, Eta) to the given location given a state vector.
!> The type of the variable being interpolated is obs_qty since
!> normally this is used to find the expected value of an observation
!> at some location. The interpolated value is returned in interp_vals
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_qty the DART KIND of interest
!> @param interp_val the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_val, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
integer,            intent(out) :: istatus(ens_size)
real(r8),           intent(out) :: interp_val(ens_size) !< array of interpolated values

integer  :: varid
real(r8) :: lon_lat_vert(3)

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_val = MISSING_R8
istatus = 99

lon_lat_vert = get_location(location)

varid = get_varid_from_kind(domain_id, obs_qty)

if (varid < 0) then
   if(debug > 12) then
      write(string1,*)'did not find obs_qty ', obs_qty, ' in the state'
      call error_handler(E_MSG,'model_interpolate:',string1,source,revision,revdate)
   endif
   return
endif

write(string1,*)'model_interpolate should not be called.'
write(string2,*)'we are getting forward observations directly from CAM'
call error_handler(E_MSG,'model_interpolate:',string1,source,revision,revdate, text2=string2)

end subroutine model_interpolate


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
! deallocate(state_loc)

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
integer     :: i

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
if (minimal_output) then
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
! Reference Pressure


end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
!> writes the time of the current state and (optionally) the time
!> to be conveyed to ROMS to dictate the length of the forecast.
!> This file is then used by scripts to modify the ROMS run.
!> The format in the time information is totally at your discretion.
!>
!> @param ncfile_out name of the file
!> @param model_time the current time of the model state
!> @param adv_to_time the time in the future of the next assimilation.
!>

subroutine write_model_time(ncid, model_time, adv_to_time)
integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

integer :: io, varid, seconds, days
type(time_type) :: origin_time, deltatime
real(digits12)  :: run_duration

if ( .not. module_initialized ) call static_init_model

!>@todo need to put code to write cam model time

if (present(adv_to_time)) then
   write(string1,*)'CAM/DART not configured to advance CAM.'
   write(string2,*)'called with optional advance_to_time '
   call error_handler(E_ERR, 'write_model_time', string1, &
              source, revision, revdate, text2=string2)
endif

end subroutine write_model_time

!--------------------------------------------------------------------
!>
!> read the time from the input file
!>
!> @param filename name of file that contains the time
!>

function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type)              :: read_model_time

integer :: ncid
integer :: timesize
integer :: datefull, datesec
integer :: iyear, imonth, iday, ihour, imin, isec, rem

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) trim(filename), ' does not exist.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))

! CAM initial files have two variables of length 
! 'time' (the unlimited dimension): date, datesec
! This code require that the time length be size 1

call nc_get_variable_size(ncid, 'time', timesize)

!>@todo do we really need to ceck this if it is never going to happen.
!#! if (timesize /= 1) then
!#!    write(string1,*) trim(filename),' has',timesize,'times. Require exactly 1.'
!#!    call error_handler(E_ERR, 'read_model_time', string1, source, revision, revdate)
!#! endif


call nc_get_variable(ncid, 'date',    datefull)
call nc_get_variable(ncid, 'datesec', datesec)

! The 'date' is YYYYMMDD ... datesec is 'current seconds of current day'
iyear  = datefull / 10000
rem    = datefull - iyear*10000
imonth = rem / 100
iday   = rem - imonth*100

ihour  = datesec / 3600
rem    = datesec - ihour*3600
imin   = rem / 60
isec   = rem - imin*60

! some cam files are from before the start of the gregorian calendar.
! since these are 'arbitrary' years, just change the offset.
if (iyear < 1601) then
   write(string1,*)' '
   write(string2,*)'WARNING - ',trim(filename),' changing year from ', &
                   iyear,'to',iyear+1601

   call error_handler(E_MSG, 'read_model_time', string1, source, revision, &
                      revdate, text2=string2,text3='to make it a valid Gregorian date.')

   write(string1,*)' '
   call error_handler(E_MSG, 'read_model_time', string1, source, revision)
   iyear = iyear + 1601
endif

read_model_time = set_date(iyear,imonth,iday,ihour,imin,isec)

call nc_check( nf90_close(ncid), 'read_model_time', 'close '//trim(filename))

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

character(len=NF90_MAX_NAME) :: varname    ! column 1, NetCDF variable name
character(len=NF90_MAX_NAME) :: dartstr    ! column 2, DART KIND
character(len=NF90_MAX_NAME) :: minvalstr  ! column 3, Clamp min val
character(len=NF90_MAX_NAME) :: maxvalstr  ! column 4, Clamp max val
character(len=NF90_MAX_NAME) :: updatestr  ! column 5, Update output or not

character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  :: update_list(MAX_STATE_VARIABLES)   = .FALSE.
integer  ::   kind_list(MAX_STATE_VARIABLES)   = MISSING_I
real(r8) ::  clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8


nfields = 0
ParseVariables : do i = 1, MAX_STATE_VARIABLES

   varname   = trim(variable_array(num_state_table_columns*i-4))
   dartstr   = trim(variable_array(num_state_table_columns*i-3))
   minvalstr = trim(variable_array(num_state_table_columns*i-2))
   maxvalstr = trim(variable_array(num_state_table_columns*i-1))
   updatestr = trim(variable_array(num_state_table_columns*i  ))

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

if (debug > 2) call state_structure_info(domain_id)

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

integer :: ncid

! put this in a subroutine that deals with the grid
call nc_check( nf90_open(grid_file, NF90_NOWRITE, ncid), &
               'read_grid_info', 'open '//trim(grid_file))

! Get the grid info
call get_cam_grid(ncid, grid)

! Set up the interpolation structure for later 
call setup_interpolation(grid)

call nc_check( nf90_close(ncid), 'read_grid_info', 'close '//trim(grid_file))

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


subroutine fill_cam_0d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

grid_array%nsize = 1
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals)

if (debug > 10) then
   print*, 'variable name ', trim(varname), grid_array%vals
endif

end subroutine fill_cam_0d_array


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
if (debug > 10) then
   per_line = 5
   print*, 'variable name ', trim(varname)
   do i=1, grid_array%nsize, per_line
      print*,  grid_array%vals(i:min(grid_array%nsize,i+per_line-1))
   enddo
endif

end subroutine fill_cam_1d_array

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


subroutine read_cam_phis_array(phis_filename)
character(len=*),   intent(in)    :: phis_filename

integer :: ncid, nsize(2)

call nc_check( nf90_open(phis_filename, NF90_NOWRITE, ncid), &
               'read_cam_phis_array', 'open '//trim(phis_filename))

call nc_get_variable_size(ncid, 'PHIS', nsize(:))
allocate( phis(nsize(1), nsize(2)) )

call nc_get_variable(ncid, 'PHIS', phis)

call nc_check( nf90_close(ncid), 'read_cam_phis_array', 'close '//trim(phis_filename))

end subroutine read_cam_phis_array

!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
