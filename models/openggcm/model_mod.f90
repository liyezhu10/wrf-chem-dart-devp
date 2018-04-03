! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This is the interface between the OpenGGCM space weather model and DART.

module model_mod

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, i4, i8, SECPERDAY, MISSING_R8, rad2deg, PI, &
                             earth_radius, vtablenamelength, digits12

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time, &
                             print_time, print_date, set_calendar_type,         &
                             operator(*),  operator(+), operator(-),            &
                             operator(>),  operator(<), operator(/),            &
                             operator(/=), operator(<=), GREGORIAN

use     location_mod, only : location_type, set_location, get_location,         & 
                             VERTISUNDEF, VERTISHEIGHT, write_location,         &
                             get_close_type,  get_dist, is_vertical,            &
                             convert_vertical_obs, convert_vertical_state

use    utilities_mod, only : register_module, error_handler,                    &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,       &
                             do_output, to_upper,                               &
                             find_namelist_in_file, check_namelist_read,        &
                             file_exist, find_textfile_dims, file_to_text,      &
                             do_nml_file, do_nml_term, nmlfileunit, open_file,  &
                             close_file

use  netcdf_utilities_mod, only : nc_check, nc_get_variable, &
                                  nc_add_global_attribute, &
                                  nc_add_global_creation_time, &
                                  nc_end_define_mode, &
                                  nc_synchronize_file

use          obs_kind_mod, only : QTY_ELECTRON_DENSITY, QTY_ELECTRIC_POTENTIAL, &
                                  get_index_for_quantity, get_name_for_quantity

use     mpi_utilities_mod, only : my_task_id, task_count

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use  ensemble_manager_mod, only : ensemble_type, map_pe_to_task, get_copy_owner_index, &
                                  get_var_owner_index

use distributed_state_mod, only : get_state

use   state_structure_mod, only : add_domain, get_model_variable_indices,        &
                                  get_num_variables, get_index_start,            &
                                  get_num_dims, get_domain_size, get_kind_index, &
                                  get_varid_from_kind, get_dart_vector_index,    &
                                  get_dim_name, get_variable_name,               &
                                  state_structure_info

use obs_def_utilities_mod, only : track_status

use     default_model_mod, only : get_close_obs, get_close_state, nc_write_model_vars

use              cotr_mod, only : transform, cotr_set, cotr, xyzdeg, degxyz

use typesizes
use netcdf 

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

!> required routines with code in this module
public :: get_model_size,                &
          get_state_meta_data,           &
          model_interpolate,             &
          shortest_time_between_assimilations, &
          static_init_model,             &
          end_model,                     &
          init_time,                     &
          init_conditions,               &
          adv_1step,                     &
          pert_model_copies,             &
          nc_write_model_atts,           &
          read_model_time,               &
          write_model_time

!> required routines where code is in other modules
public :: get_close_obs,                 &
          get_close_state,               &
          nc_write_model_vars,           &
          convert_vertical_obs,          &
          convert_vertical_state

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! Transform type openggcm, this come directly from the 
! openggcm source code that had been modified and put
! as a subroutine in cotr_mod.nml
type(transform) :: openggcm_transform

! Message strings
character(len=512) :: string1
character(len=512) :: string2
character(len=512) :: string3

integer, parameter :: VERT_LEVEL_1 = 1

! Grid types
integer, parameter :: GEOGRAPHIC_GRID = 1
integer, parameter :: MAGNETIC_GRID   = 2

! Grid convertion direction
integer, parameter :: SM_TO_GEO = 1
integer, parameter :: GEO_TO_SM = 2

! Logical to keep track of if we have initialized static_init_model
logical, save :: module_initialized = .false.

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10 
integer, parameter :: num_state_table_columns = 4
! NOTE: may need to increase character length if netcdf variables are
! larger than vtablenamelength = 64.
character(len=vtablenamelength) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
logical ::  update_var_list( max_state_variables )
integer ::   grid_info_list( max_state_variables )
integer ::   dim_order_list( max_state_variables, 3 )

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX   = 1
integer, parameter :: VAR_KIND_INDEX   = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

! identifiers for LAT, LON and HEIGHT
! The 'natural' order of the variables is height varies fastest, then lat, then lon.
! in ncdump  ... (lon,lat,height)
! in fortran ... (height,lat,lon)
integer, parameter :: VAR_HGT_INDEX   = 1
integer, parameter :: VAR_LAT_INDEX   = 2
integer, parameter :: VAR_LON_INDEX   = 3

! things which can/should be in the model_nml
character(len=NF90_MAX_NAME) :: openggcm_template = 'DATA.ionos2.nc'
integer  :: assimilation_period_days     = 1
integer  :: assimilation_period_seconds  = 0
real(r8) :: model_perturbation_amplitude = 0.2
character(len=vtablenamelength) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 0   ! turn up for more and more debug messages

namelist /model_nml/  &
   openggcm_template,           &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   model_state_variables,       &
   debug

! Mon Jul 18 17:00:01 MDT 2016 Jimmy Raeder : 
! I had at the CEDAR meeting some longish discussions with Naomi.
! I learned that the O+ grid is actually different than I thought.
! It is based on dipole flux tubes, such that the footprints of the 
! flux tubes are a regular lat-lon grid (call these dimensions ph and th 
! for longitude and latitude) with indices iph and ith.  Then, in the 
! 3rd dimension (call it z, with index iz) the points are irregularly 
! spaced along flux tubes.  Thus, the grid is topologically a cube 
! (1..nph, 1..nth, 1..nz), but each location depends on all 3 indices, 
! i.e., we have ph(iph,ith,iz), th(iph,ith,iz), and z(iph,ith,iz), 
! and each holds a value, say v(iph,with,iz), O+ for example.  
! Although the grid is static, I suggest that we compute and write 
! out all for DART.

!>@todo since the grid is static, it should be written ONCE by openggcm
!> at startup ... to a separate file from the data.

! Generic grid type to hold CTIM and Magnetic grid information
type grid_type
   ! the size of the grid
   integer :: nlon, nlat, nheight

   ! grid information
   real(r8), allocatable :: longitude(:,:,:)
   real(r8), allocatable ::  latitude(:,:,:)
   real(r8), allocatable ::    height(:,:,:)
   logical :: uses_colatitude

   ! optional conversion grid - 2d arrays: (lon, lat)
   real(r8), allocatable :: conv_2d_lon(:,:), conv_2d_lat(:,:)
end type grid_type


!------------------------------------------------------------------
! Global Variables 
!------------------------------------------------------------------

! Number of fields in the state vector
integer :: nfields

type(grid_type), target :: geo_grid   ! Geometric Grid (oplus)
type(grid_type), target :: mag_grid   ! Magnetic Grid

! Global Time Variables
type(time_type) :: model_time, model_timestep

! The state vector length
integer(i8) :: model_size

! Domain id to be used by routines in state_structure_mod
integer :: domain_id

contains

!------------------------------------------------------------------
!------------------------------------------------------------------

!> Called to do one time initialization of the model. In this case,
!> it reads in the grid information.

subroutine static_init_model()

integer :: iunit, io, ncid
integer :: ss, dd

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
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! set calendar type
call set_calendar_type(GREGORIAN)

model_timestep = set_time(assimilation_period_days,assimilation_period_seconds)

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

! get data dimensions, then allocate space, then open the files
! and actually fill in the arrays.

ncid = get_grid_template_fileid(openggcm_template)

call get_grid_sizes(ncid, geo_grid, 'geo_lon', 'geo_lat','geo_height')
! call get_grid_sizes(ncid, mag_grid, 'ig_lon', 'ig_lat')

! allocate space for geographic and magnetic grids

call allocate_grid_space(geo_grid, conv=.false.)
! call allocate_grid_space(mag_grid, conv=.true.)

! read in geographic and magnetic grids, and for the mag grid read
! in the 2d conversion arrays to go to geographic coords

call read_grid(ncid, geo_grid, 'geo_lon', 'geo_lat', 'geo_height', is_conv=.false., is_co_latitude=.false.)
! call read_grid(ncid, mag_grid, 'ig_lon',  'ig_lat', 'ig_height', is_conv=.false., is_co_latitude=.true.) 

! verify that the model_state_variables namelist was filled in correctly.  
! returns variable_table which has variable names, kinds and update strings, 
! and grid information.

call verify_state_variables(model_state_variables, nfields, variable_table, &
                            state_kinds_list, update_var_list, grid_info_list)

! fill up the state structure with information from the model template file

domain_id = add_domain(openggcm_template, nfields, &
                       var_names   = variable_table  (1:nfields , VAR_NAME_INDEX), &
                       kind_list   = state_kinds_list(1:nfields), &
                       update_list = update_var_list (1:nfields))

if (debug > 0) call state_structure_info(domain_id)

! order the dimensions according to lat, lon and height

call make_dim_order_table(nfields)

! only one domain in the model at the moment.

model_size = get_domain_size(domain_id)
write(string1,*)'model_size = ', model_size
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

! set the transform geo -> magnetic grid
call initialize_openggcm_transform(openggcm_template)

end subroutine static_init_model

!------------------------------------------------------------------

!> Returns the size of the model state vector as an I8 integer.
!> The size is computed in static_init_model() an stored in a
!> module global variable.

function get_model_size()

integer(i8) :: get_model_size  !< state vector length

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------

!> Model interpolate will interpolate any state variable
!> the given location given a state vector. The 'generic kind' of the variable being
!> interpolated is obs_kind since normally this is used to find the expected
!> value of an observation at some location. The interpolated value is 
!> returned in interp_val and istatus is 0 for success.

subroutine model_interpolate(state_handle, ens_size, location, obs_kind, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle !< ensemble handle for data to interpolate in
integer,             intent(in) :: ens_size !< number of ensembles, sets size of expected_obs and istatus arrays
type(location_type), intent(in) :: location !< dart location to interpolate to
integer,             intent(in) :: obs_kind !< dart kind to interpolate
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size) !< array of returned statuses

! Local storage
real(r8)    :: loc_array(3), llon, llat, lheight
integer     :: ind
integer     :: hgt_bot, hgt_top
real(r8)    :: hgt_fract
integer     :: hstatus
type(grid_type), pointer :: mygrid

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

expected_obs(:) = MISSING_R8     ! the DART bad value flag
istatus(:)      = 99             ! unknown error

! Get the individual locations values
loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
lheight = loc_array(3)

if ( is_vertical(location, "UNDEFINED") ) then
   ! this is what we expect and it is ok
elseif ( is_vertical(location, "HEIGHT") ) then
   ! this is what we expect and it is ok
   ! once we write the code to search in the vertical
elseif ( is_vertical(location, "LEVEL") ) then
   write(string1,*)'requesting interp of an obs on level, not supported yet'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)

   !>@todo FIXME something like this
   ! convert the heights index to an actual height 
   ind = nint(loc_array(3))
   if ( ind < 1 .or. ind > geo_grid%nheight ) then 
      istatus = 11
      return
   else
      lheight = ind
   endif

else   ! if pressure or surface we don't know what to do
   write(string1,*)'requesting interp of an obs on pressure or surface, not supported yet'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)

   istatus = 17
   return
endif

! Determine which grid the incoming observation is on. If the observation
! is on the magnetic grid transform it to the geographic grid.

SELECT CASE (get_grid_type(obs_kind))
   CASE (MAGNETIC_GRID)
      call transform_mag_geo(llon, llat, lheight, GEO_TO_SM)
      mygrid => mag_grid

   CASE (GEOGRAPHIC_GRID)
      mygrid => geo_grid

   CASE DEFAULT
      call error_handler(E_ERR, 'model_interpolate', 'unknown grid type, should not happen', &
            source, revision, revdate)
END SELECT


if( is_vertical(location,"UNDEFINED") ) then
   call lon_lat_interpolate(state_handle, ens_size, mygrid, obs_kind, llon, llat, VERT_LEVEL_1, &
                            expected_obs, istatus)

   return
endif

!>@todo  the heights vary at each location ... 
call height_bounds(lheight, mygrid%nheight, mygrid%height, hgt_bot, hgt_top, hgt_fract, hstatus)

if(hstatus /= 0) then
   istatus = 12
   return
endif

! do a 2d interpolation for the value at the bottom level, then again for
! the top level, then do a linear interpolation in the vertical to get the
! final value.  this sets both interp_val and istatus.
call do_interp(state_handle, ens_size, mygrid, hgt_bot, hgt_top, hgt_fract, &
               llon, llat, obs_kind, expected_obs, istatus)

end subroutine model_interpolate

!------------------------------------------------------------------

!> Subroutine to interpolate to a lon lat location given the state handle.
!> Successful interpolation returns istatus=0.

subroutine lon_lat_interpolate(state_handle, ens_size, grid_handle, var_kind, &
                               lon, lat, height_index, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle           !< state ensemble handle
integer,             intent(in)  :: ens_size               !< ensemble size
type(grid_type),     intent(in)  :: grid_handle            !< geo or mag grid
integer,             intent(in)  :: var_kind               !< dart variable kind
real(r8),            intent(in)  :: lon                    !< longitude to interpolate
real(r8),            intent(in)  :: lat                    !< latitude to interpolate
integer,             intent(in)  :: height_index           !< height index to interpolate
real(r8),            intent(out) :: expected_obs(ens_size) !< returned interpolations
integer,             intent(out) :: istatus(ens_size)      !< returned statuses

! Local storage, 
integer  :: lat_bot, lat_top, lon_bot, lon_top
real(r8) :: p(4,ens_size), xbot(ens_size), xtop(ens_size)
real(r8) :: lon_fract, lat_fract

! Succesful return has istatus of 0
istatus = 0

! find the lower and upper indices which enclose the given value
! in this model, the data at lon 0 is replicated at lon 360, so no special
! wrap case is needed.
call lon_bounds(lon, grid_handle, lon_bot, lon_top, lon_fract)

if (grid_handle%uses_colatitude) then
   call colat_bounds(lat, grid_handle, lat_bot, lat_top, lat_fract, istatus(1))
else
   call lat_bounds(lat, grid_handle, lat_bot, lat_top, lat_fract, istatus(1))
endif

if (istatus(1) /= 0) then
   istatus(:) = 18 
   return
endif

! Get the values at the four corners of the box or quad
! Corners go around counterclockwise from lower left
p(1, :) = get_val(lon_bot, lat_bot, height_index, var_kind, state_handle, ens_size)
p(2, :) = get_val(lon_top, lat_bot, height_index, var_kind, state_handle, ens_size)
p(3, :) = get_val(lon_top, lat_top, height_index, var_kind, state_handle, ens_size)
p(4, :) = get_val(lon_bot, lat_top, height_index, var_kind, state_handle, ens_size)

! Rectangular bilinear interpolation
xbot = p(1, :) + lon_fract * (p(2, :) - p(1, :))
xtop = p(4, :) + lon_fract * (p(3, :) - p(4, :))

! Now interpolate in latitude
expected_obs = xbot + lat_fract * (xtop - xbot)

end subroutine lon_lat_interpolate

!------------------------------------------------------------

!> Returns the value given the lon, lat, and height indices

function get_val(lon_index, lat_index, height_index, var_kind, state_handle, ens_size)

integer,             intent(in)  :: lon_index     !< longitude index
integer,             intent(in)  :: lat_index     !< latitude index
integer,             intent(in)  :: height_index  !< height index
integer,             intent(in)  :: var_kind      !< dart variable kind
type(ensemble_type), intent(in)  :: state_handle  !< ensemble handle for state vector
integer,             intent(in)  :: ens_size      !< size of the ensemble

! Local variables
real(r8)    :: get_val(ens_size)
integer(i8) :: state_index
integer     :: var_id

integer, dimension(3) :: dim_index

if (var_kind < 0 ) then
   write(string1,*) 'dart kind < 0 which should not happen'
   write(string2,*) 'the dart kind provided is : ', var_kind
   call error_handler(E_ERR, 'get_val', string1, &
                      source, revision, revdate, text2=string2)
endif   

var_id = get_varid_from_kind(domain_id, var_kind)

! This should take of any any permutation of storage order.
! The use of dim_order_list means any variable can be stored
! in any order.
dim_index(dim_order_list(var_id, VAR_HGT_INDEX)) = height_index
dim_index(dim_order_list(var_id, VAR_LAT_INDEX)) = lat_index
dim_index(dim_order_list(var_id, VAR_LON_INDEX)) = lon_index

! the order of the arguments to get_dart_vector_index() is always
! in the same order as the native storage order.
state_index = get_dart_vector_index(dim_index(1), dim_index(2), dim_index(3), &
                                    domain_id, var_id)

get_val = get_state(state_index, state_handle)

end function get_val

!------------------------------------------------------------

!> Given a longitude lon and a grid handle which contains both the 1D array 
!> of longitudes and the grid longitude size, returns the indices of the grid
!> below and above the location longitude and the fraction of the distance
!> between.  This code assumes that the first and last rows are replicated
!> and identical (e.g. 0 and 360 both have entries in the array)

subroutine lon_bounds(lon, grid_handle, bot, top, fract)

real(r8),        intent(in)  :: lon          !< input longitude
type(grid_type), intent(in)  :: grid_handle  !< handle to either a geo or mag grid
integer,         intent(out) :: bot          !< index of bottom layer
integer,         intent(out) :: top          !< index of top layer
real(r8),        intent(out) :: fract        !< fraction between layers

! Local storage
integer  :: i

call error_handler(E_ERR, 'lon_bounds', 'routine needs to be rewritten for fully 3D grid', &
                   source, revision, revdate )

!todo do i = 2, grid_handle%nlon
!todo    if (lon <= grid_handle%longitude(i)) then
!todo       bot = i-1
!todo       top = i
!todo       fract = (lon - grid_handle%longitude(bot)) / &
!todo               (grid_handle%longitude(top) - grid_handle%longitude(bot))
!todo       return
!todo    endif
!todo enddo
!todo 
!todo write(string1, *) 'looking for lon ', lon
!todo call error_handler(E_ERR, 'lon_bounds', 'reached end of loop without finding lon', &
!todo                    source, revision, revdate, text2=string1)

end subroutine lon_bounds

!-------------------------------------------------------------

!> Given a latitude lat and the grid_handle which contains both the 
!> 1D array of latitudes and the grid latitude count, returns the
!> indices of the grid below and above the location latitude and 
!> the fraction of the distance between. istatus is returned as 0 
!> unless the location latitude is south of the southernmost grid 
!> point (1 returned) or north of the northernmost (2 returned),
!> which may not be possible anymore and possibly could be removed.

subroutine lat_bounds(lat, grid_handle, bot, top, fract, istatus)

real(r8),        intent(in)  :: lat          !< input latitude
type(grid_type), intent(in)  :: grid_handle  !< geo or mag grid
integer,         intent(out) :: bot          !< index of bottom layer
integer,         intent(out) :: top          !< index of top layer
real(r8),        intent(out) :: fract        !< fraction between layers
integer,         intent(out) :: istatus      !< return status

! Local storage
integer :: i

! Success should return 0, failure a positive number.
istatus = 0

call error_handler(E_ERR, 'lat_bounds', 'routine needs to be rewritten for fully 3D grid', &
                   source, revision, revdate )

!todo ! Check for too far south or north
!todo if(lat < grid_handle%latitude(1)) then
!todo    istatus = 1
!todo    return
!todo else if(lat > grid_handle%latitude(grid_handle%nlat)) then
!todo    istatus = 2
!todo    return
!todo endif
!todo 
!todo ! In the middle, search through
!todo do i = 2, grid_handle%nlat
!todo    if(lat <= grid_handle%latitude(i)) then
!todo       bot = i - 1
!todo       top = i
!todo       fract = (lat - grid_handle%latitude(bot)) / &
!todo               (grid_handle%latitude(top) - grid_handle%latitude(bot))
!todo       return
!todo    endif
!todo enddo
!todo 
!todo write(string1, *) 'looking for lat ', lat
!todo call error_handler(E_ERR, 'lat_bounds', 'reached end of loop without finding lat', &
!todo                    source, revision, revdate, text2=string1)

end subroutine lat_bounds

!-------------------------------------------------------------

!> Given a latitude lat, the grid handle which contains the 1d array of 
!> colatitudes for the grid and the grid count, return the indices of
!> the grid below and above the location colatitude and the fraction 
!> of the distance between. colatitudes start at 0 and go to 180, but
!> to be consistent with our locations mod we have already transformed
!> them into 90 to -90.  this routine has to be different because the
!> order of the points is north pole to south, while latitudes are ordered
!> south pole to north.  we have to search in a different order and the
!> test itself is reversed from the lat_bounds() routine.
!> istatus is returned as 0 unless the location latitude is 
!> south of the southernmost grid point (1 returned) or north of the 
!> northernmost (2 returned). given our locations module i believe this
!> test is no longer needed since the grid includes the poles.

subroutine colat_bounds(lat, grid_handle, bot, top, fract, istatus)

real(r8),        intent(in)  :: lat          !< input latitude
type(grid_type), intent(in)  :: grid_handle  !< geo or mag grid
integer,         intent(out) :: bot          !< index of bottom layer
integer,         intent(out) :: top          !< index of top layer
real(r8),        intent(out) :: fract        !< fraction between layers
integer,         intent(out) :: istatus      !< return status

! Local storage
integer :: i

call error_handler(E_ERR, 'colat_bounds', 'routine needs to be rewritten for fully 3D grid', &
                   source, revision, revdate )

! Success should return 0, failure a positive number.
istatus = 0

!todo ! Check for too far south or north
!todo if(lat > grid_handle%latitude(1)) then
!todo    istatus = 1
!todo    return
!todo else if(lat < grid_handle%latitude(grid_handle%nlat)) then
!todo    istatus = 2
!todo    return
!todo endif
!todo 
!todo ! In the middle, search through
!todo do i = 2, grid_handle%nlat
!todo    if(lat >= grid_handle%latitude(i)) then
!todo       bot = i - 1
!todo       top = i
!todo       fract = (lat - grid_handle%latitude(bot)) / &
!todo               (grid_handle%latitude(top) - grid_handle%latitude(bot))
!todo       return
!todo    endif
!todo enddo
!todo 
!todo write(string1, *) 'looking for colat ', lat
!todo call error_handler(E_ERR, 'colat_bounds', 'reached end of loop without finding colat', &
!todo                    source, revision, revdate, text2=string1)

end subroutine colat_bounds

!------------------------------------------------------------

!> find the index top and bottom index for a variable given an lheight and an
!> array of heights.

subroutine height_bounds(lheight, nheights, hgt_array, bot, top, fract, istatus)

real(r8),   intent(in)  :: lheight             !< height location
integer,    intent(in)  :: nheights            !< number of total heights
real(r8),   intent(in)  :: hgt_array(nheights) !< array of heights
integer,    intent(out) :: bot                 !< bottom bounding height
integer,    intent(out) :: top                 !< top bounding height
real(r8),   intent(out) :: fract               !< fraction inbetween
integer,    intent(out) :: istatus             !< return status


! Local variables
integer   :: i

! Succesful istatus is 0
! Make any failure here return istatus in the 20s
istatus = 0

if(lheight <= hgt_array(1)) then
   bot = 1
   top = 2
   ! NOTE: the fract definition is the relative distance from bottom to top
   fract = 1.0_r8 
   return
endif

! Search through the boxes
do i = 2, nheights
   ! If the location is lower than this entry, it must be in this box
   if(lheight < hgt_array(i)) then
      top = i
      bot = i -1
      fract = (lheight - hgt_array(bot)) / (hgt_array(top) - hgt_array(bot))
      return
   endif
enddo

! Falling off the end means the location is higher than the model top.
bot   = -1
top   = -1
fract = -1.0_r8

istatus = 20

end subroutine height_bounds

!------------------------------------------------------------------

!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. In fact this sets the assimilation window size.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations !< returned timestep

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location and (optionally) the generic kind.  For this model
!> we must always return geographic coordinates.  For fields on the
!> magnetic grid this requires a coordinate transformation.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in     !< dart state index of interest
type(location_type), intent(out) :: location     !< location of interest
integer, OPTIONAL,   intent(out) :: var_type     !< optional dart kind return

! Local variables
real(r8) :: lat, lon, height
integer  :: lon_index, lat_index, height_index, local_var, var_id
integer  :: state_loc(3)

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, state_loc(1), state_loc(2), state_loc(3), var_id=var_id)

local_var = get_kind_index(domain_id, var_id)

height_index = state_loc(dim_order_list(var_id, VAR_HGT_INDEX))
lat_index    = state_loc(dim_order_list(var_id, VAR_LAT_INDEX))
lon_index    = state_loc(dim_order_list(var_id, VAR_LON_INDEX))

! we are getting a mapping array between magnetic -> geogrid

if ( get_grid_type(local_var) == MAGNETIC_GRID ) then
   height   = 0.0_r8
   lat = mag_grid%conv_2d_lat(lon_index, lat_index)
   lon = mag_grid%conv_2d_lon(lon_index, lat_index)
else
   height = geo_grid%height(   height_index, lat_index, lon_index)
   lat    = geo_grid%latitude( height_index, lat_index, lon_index)
   lon    = geo_grid%longitude(height_index, lat_index, lon_index)
endif

!>@todo  Here we are ASSUMING that electric potential is a 2D
! variable.  If this is NOT the case FIX HERE!
if (local_var == QTY_ELECTRIC_POTENTIAL) then

   call error_handler(E_ERR,'get_state_meta_data','not tested with QTY_ELECTRIC_POTENTIAL', &
         source, revision, revdate, text2='check entire stream - end-to-end')

   location = set_location(lon, lat, height, VERTISUNDEF)
else
   location = set_location(lon, lat, height, VERTISHEIGHT)
endif

if (present(var_type)) then
   var_type = local_var
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------

!> Shutdown and clean-up.

subroutine end_model()

call deallocate_grid_space(geo_grid)
call deallocate_grid_space(mag_grid)

end subroutine end_model

!------------------------------------------------------------------

!> open and return the netcdf file id of template file

function get_grid_template_fileid(filename)

character(len=*), intent(in) :: filename
integer :: get_grid_template_fileid

call nc_check( NF90_open(filename, NF90_NOWRITE, get_grid_template_fileid), &
                  'get_grid_template_fileid', 'open '//trim(filename))

end function get_grid_template_fileid

!------------------------------------------------------------------

!> get grid sizes given netcdf id, lon, lat and height name. results
!> are store in the provided grid handle

subroutine get_grid_sizes(ncid, grid_handle, lon_name, lat_name, height_name)

integer,                    intent(in)    :: ncid        !< netcdf file id
type(grid_type),            intent(inout) :: grid_handle !< geo or mag grid
character(len=*),           intent(in)    :: lon_name    !< longitude name
character(len=*),           intent(in)    :: lat_name    !< latitude name
character(len=*), OPTIONAL, intent(in)    :: height_name !< height name

grid_handle%nlon = get_dim(ncid,lon_name, 'get_grid_sizes')
grid_handle%nlat = get_dim(ncid,lat_name, 'get_grid_sizes')

if (present(height_name)) then
   grid_handle%nheight = get_dim(ncid, height_name, 'get_grid_sizes')
else
   grid_handle%nheight = 1
endif

end subroutine get_grid_sizes

!------------------------------------------------------------------

!> is_conv:  if true, read the data into the conversion grid.
!> otherwise read into the normal lat/lon arrays.
!> 
!> is_co_latitude:  if true, subtract 90 from the lat values
!> co_latitudes start at 0 at the north pole and go to 180 at the south.
!> "normal" latitudes for us are -90 at the south pole up to 90 at the north.

subroutine read_grid(ncid, grid_handle, lon_name, lat_name, z_name, is_conv, is_co_latitude)

integer,          intent(in)    :: ncid           !< netcdf file id for grid
type(grid_type),  intent(inout) :: grid_handle    !< grid handle to be read into
character(len=*), intent(in)    :: lon_name       !< longitude variable name
character(len=*), intent(in)    :: lat_name       !< latitude variable name
character(len=*), intent(in)    :: z_name         !< height variable name
logical,          intent(in)    :: is_conv        !< fill conversion grid
logical,          intent(in)    :: is_co_latitude !< is grid in co-latitude

integer :: hgtdim, latdim, londim

if (is_conv) then

   call nc_get_variable(ncid, lon_name, grid_handle%conv_2d_lon, 'read_conv_horiz_grid')
   call nc_get_variable(ncid, lat_name, grid_handle%conv_2d_lat, 'read_conv_horiz_grid')

   if (is_co_latitude) then
      grid_handle%conv_2d_lat(:,:) = 90.0_r8 - grid_handle%conv_2d_lat(:,:)
      grid_handle%uses_colatitude = .true.
   else
      grid_handle%uses_colatitude = .false.
   endif

   where(grid_handle%conv_2d_lon < 0.0_r8) grid_handle%conv_2d_lon = grid_handle%conv_2d_lon + 360.0_r8

else

   call nc_get_variable(ncid, lon_name, grid_handle%longitude, 'read_grid')
   call nc_get_variable(ncid, lat_name, grid_handle%latitude,  'read_grid')
   call nc_get_variable(ncid,   z_name, grid_handle%height,    'read_grid')

   if (is_co_latitude) then
      grid_handle%latitude = 90.0_r8 - grid_handle%latitude
      grid_handle%uses_colatitude = .true.
   else
      grid_handle%uses_colatitude = .false.
   endif

   where(grid_handle%longitude < 0.0_r8) grid_handle%longitude = grid_handle%longitude + 360.0_r8

endif

end subroutine read_grid

!------------------------------------------------------------------

!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables for the geometric and
!> magnetic grids.

subroutine nc_write_model_atts( ncid, domain_id )

integer, intent(in)  :: ncid
integer, intent(in)  :: domain_id

!----------------------------------------------------------------------
! local variables 
!----------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! for the dimensions and coordinate variables
integer :: NlonDimID, NlatDimID, NhgtDimID
integer :: geoLonVarID, geoLatVarID, geoHeightVarID
integer :: magLonVarID, magLatVarID, magHeightVarID
integer ::  coLonVarID,  coLatVarID

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=128) :: filename

if ( .not. module_initialized ) call static_init_model

!--------------------------------------------------------------------
! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.
!--------------------------------------------------------------------

write(filename,*) 'ncid', ncid

!-------------------------------------------------------------------------------
! make sure ncid refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                   'nc_write_model_atts', 'inquire '//trim(filename))
call nc_check(nf90_Redef(ncid),'nc_write_model_atts',   'redef '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "openggcm")

!----------------------------------------------------------------------------
! We need to output grid information
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

NlonDimID = set_dim(ncid, 'geo_lon',    geo_grid%nlon,    filename)
NlatDimID = set_dim(ncid, 'geo_lat',    geo_grid%nlat,    filename)
NhgtDimID = set_dim(ncid, 'geo_height', geo_grid%nheight, filename)

!----------------------------------------------------------------------------
! Write out Geographic Grid attributes
!----------------------------------------------------------------------------

call nc_check(NF90_def_var(ncid,name='geo_lon', xtype=NF90_real, &
              dimids=(/NhgtDimID, NlatDimID, NlonDimID /), varid=geoLonVarID),&
              'nc_write_model_atts', 'geo_lon def_var '//trim(filename))
call nc_check(NF90_put_att(ncid,  geoLonVarID, 'long_name', 'geographic longitudes'), &
              'nc_write_model_atts', 'geo_lon long_name '//trim(filename))

! Grid Latitudes
call nc_check(NF90_def_var(ncid,name='geo_lat', xtype=NF90_real, &
              dimids=(/NhgtDimID, NlatDimID, NlonDimID /), varid=geoLatVarID),&
              'nc_write_model_atts', 'geo_lat def_var '//trim(filename))
call nc_check(NF90_put_att(ncid,  geoLatVarID, 'long_name', 'geographic latitudes'), &
              'nc_write_model_atts', 'geo_lat long_name '//trim(filename))

! Heights
call nc_check(NF90_def_var(ncid,name='geo_height', xtype=NF90_real, &
              dimids=(/ NhgtDimID, NlatDimID, NlonDimID /), varid= geoHeightVarID), &
              'nc_write_model_atts', 'geo_height def_var '//trim(filename))
call nc_check(NF90_put_att(ncid, geoHeightVarID, 'long_name', 'heights of geographic grid'), &
              'nc_write_model_atts', 'geo_height long_name '//trim(filename))

!TJH !----------------------------------------------------------------------------
!TJH ! Write out Magnetic Grid attributes
!TJH !----------------------------------------------------------------------------
!TJH 
!TJH NlonDimID = set_dim(ncid, 'mag_lon', mag_grid%nlon,    filename)
!TJH NlatDimID = set_dim(ncid, 'mag_lat', mag_grid%nlat,    filename)
!TJH NhgtDimID = set_dim(ncid, 'mag_hgt', mag_grid%nheight, filename)
!TJH 
!TJH ! Grid Longitudes
!TJH call nc_check(NF90_def_var(ncid,name='mag_lon', xtype=NF90_real, &
!TJH               dimids=NlonDimID, varid=magLonVarID),&
!TJH               'nc_write_model_atts', 'mag_lon def_var '//trim(filename))
!TJH call nc_check(NF90_put_att(ncid,  magLonVarID, 'long_name', 'longitudes mag grid'), &
!TJH               'nc_write_model_atts', 'mag_lon long_name '//trim(filename))
!TJH 
!TJH ! Grid Latitudes
!TJH call nc_check(NF90_def_var(ncid,name='mag_lat', xtype=NF90_real, &
!TJH               dimids=NlatDimID, varid=magLatVarID),&
!TJH               'nc_write_model_atts', 'mag_lat def_var '//trim(filename))
!TJH call nc_check(NF90_put_att(ncid,  magLatVarID, 'long_name', 'latitudes of mag grid'), &
!TJH               'nc_write_model_atts', 'mag_lat long_name '//trim(filename))
!TJH 
!TJH ! Heights
!TJH call nc_check(NF90_def_var(ncid,name='mag_height', xtype=NF90_real, &
!TJH               dimids=NhgtDimID, varid= magHeightVarID), &
!TJH               'nc_write_model_atts', 'mag_height def_var '//trim(filename))
!TJH call nc_check(NF90_put_att(ncid, magHeightVarID, 'long_name', 'height for mag grid'), &
!TJH               'nc_write_model_atts', 'mag_height long_name '//trim(filename))
!TJH 
!TJH !----------------------------------------------------------------------------
!TJH ! Write out Co-Latitude Grid attributes
!TJH !----------------------------------------------------------------------------
!TJH 
!TJH NlonDimID = set_dim(ncid, 'co_lon', mag_grid%nlon,    filename)
!TJH NlatDimID = set_dim(ncid, 'co_lat', mag_grid%nlat,    filename)
!TJH 
!TJH ! Grid Longitudes
!TJH call nc_check(NF90_def_var(ncid,name='co_lon', xtype=NF90_real, &
!TJH               dimids=(/ NlonDimID, NlatDimID /), varid=coLonVarID),&
!TJH               'nc_write_model_atts', 'co_lon def_var '//trim(filename))
!TJH call nc_check(NF90_put_att(ncid,  coLonVarID, 'long_name', 'longitudes co grid'), &
!TJH               'nc_write_model_atts', 'co_lon long_name '//trim(filename))
!TJH 
!TJH ! Grid Latitudes
!TJH call nc_check(NF90_def_var(ncid,name='co_lat', xtype=NF90_real, &
!TJH               dimids=(/ NlonDimID, NlatDimID /), varid=coLatVarID),&
!TJH               'nc_write_model_atts', 'co_lat def_var '//trim(filename))
!TJH call nc_check(NF90_put_att(ncid,  coLatVarID, 'long_name', 'latitudes of co grid'), &
!TJH               'nc_write_model_atts', 'co_lat long_name '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_end_define_mode(ncid, 'nc_write_model_atts enddef "'//trim(filename)//'"')

!----------------------------------------------------------------------------
! Fill Geographic Grid values
!----------------------------------------------------------------------------

call nc_check(NF90_put_var(ncid, geoLonVarID, geo_grid%longitude ), &
             'nc_write_model_atts', 'geo_lon put_var '//trim(filename))
call nc_check(NF90_put_var(ncid, geoLatVarID, geo_grid%latitude ), &
             'nc_write_model_atts', 'geo_lat put_var '//trim(filename))
call nc_check(NF90_put_var(ncid, geoHeightVarID, geo_grid%height ), &
             'nc_write_model_atts', 'geo_height put_var '//trim(filename))

!TJH !----------------------------------------------------------------------------
!TJH ! Fill Magnetic Grid values
!TJH !----------------------------------------------------------------------------
!TJH 
!TJH call nc_check(NF90_put_var(ncid, magLonVarID, mag_grid%longitude ), &
!TJH              'nc_write_model_atts', 'mag_lon put_var '//trim(filename))
!TJH call nc_check(NF90_put_var(ncid, magLatVarID, mag_grid%latitude ), &
!TJH              'nc_write_model_atts', 'mag_lat put_var '//trim(filename))
!TJH call nc_check(NF90_put_var(ncid, magHeightVarID, mag_grid%height ), &
!TJH              'nc_write_model_atts', 'mag_height put_var '//trim(filename))
!TJH 
!TJH !----------------------------------------------------------------------------
!TJH ! Fill Co-Latitude Grid values
!TJH !----------------------------------------------------------------------------
!TJH 
!TJH call nc_check(NF90_put_var(ncid, coLonVarID, mag_grid%conv_2d_lon(:,:)), &
!TJH              'nc_write_model_atts', 'co_lon put_var '//trim(filename))
!TJH call nc_check(NF90_put_var(ncid, coLatVarID, mag_grid%conv_2d_lat(:,:) ), &
!TJH              'nc_write_model_atts', 'co_lat put_var '//trim(filename))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_synchronize_file(ncid, 'nc_write_model_atts')

end subroutine nc_write_model_atts


!------------------------------------------------------------------
!> Perturbs state copies for generating initial ensembles.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle !< state ensemble handle
integer,             intent(in)    :: ens_size !< ensemble size
real(r8),            intent(in)    :: pert_amp !< perterbation amplitude
logical,             intent(out)   :: interf_provided !< have you provided an interface?

! Local Variables
integer     :: i, j
integer(i8) :: dart_index

! Storage for a random sequence for perturbing a single initial state
type(random_seq_type) :: random_seq

if ( .not. module_initialized ) call static_init_model

write(string1,*)'routine not tested'
call error_handler(E_ERR, 'pert_model_copies', string1, source, revision, revdate)

interf_provided = .true.

! Initialize random number sequence
call init_random_seq(random_seq, my_task_id())

! perturb the state using random gaussian noise
do i=1,state_ens_handle%my_num_vars
   dart_index = state_ens_handle%my_vars(i)
   do j=1, ens_size
      ! Since potential and electron density have such radically different values
      ! we weight the standard deviation with the actual state value so that noise 
      ! created is closer to the actual values in the state. NOTE: If the value is
      ! state value is zero the standard deviation will be zero and your values will
      ! not be perturbed.
      state_ens_handle%copies(j,i) = random_gaussian(random_seq, &
         mean               = state_ens_handle%copies(j,i),      &
         standard_deviation = state_ens_handle%copies(j,i)*model_perturbation_amplitude)
   enddo
enddo


end subroutine pert_model_copies

!------------------------------------------------------------------

!> do a 2d horizontal interpolation for the value at the bottom level, 
!> then again for the top level, then do a linear interpolation in the 
!> vertical to get the final value.

subroutine do_interp(state_handle, ens_size, grid_handle, hgt_bot, hgt_top, hgt_fract, &
                     llon, llat, obs_kind, expected_obs, istatus)

type(ensemble_type), intent(in)  :: state_handle           !< state ensemble handle
integer,             intent(in)  :: ens_size               !< ensemble size
type(grid_type),     intent(in)  :: grid_handle            !< geo or mag grid
integer,             intent(in)  :: hgt_bot                !< index to bottom bound
integer,             intent(in)  :: hgt_top                !< index to top bound
real(r8),            intent(in)  :: hgt_fract              !< fraction inbetween top and bottom
real(r8),            intent(in)  :: llon                   !< longitude to interpolate
real(r8),            intent(in)  :: llat                   !< latitude to interpolate
integer,             intent(in)  :: obs_kind               !< dart kind
real(r8),            intent(out) :: expected_obs(ens_size) !< interpolated value
integer,             intent(out) :: istatus(ens_size)      !< status of interpolation

! Local Variables
real(r8)    :: bot_val(ens_size), top_val(ens_size)
integer     :: temp_status(ens_size)
logical     :: return_now

istatus(:) = 0 ! need to start with istatus = 0 for track status to work properly

call lon_lat_interpolate(state_handle, ens_size, grid_handle, obs_kind, &
                         llon, llat, hgt_bot, bot_val, temp_status)
call track_status(ens_size, temp_status, bot_val, istatus, return_now)
if (return_now) return

call lon_lat_interpolate(state_handle, ens_size, grid_handle, obs_kind, &
                         llon, llat, hgt_top, top_val, temp_status)
call track_status(ens_size, temp_status, top_val, istatus, return_now)
if (return_now) return

! Then weight them by the vertical fraction and return
where (istatus == 0) 
   expected_obs = bot_val + hgt_fract * (top_val - bot_val)
elsewhere
   expected_obs = MISSING_R8
endwhere

end subroutine do_interp

!--------------------------------------------------------------------
!> read the current model time from openggcm output file
!> the time is encoded as a 64bit real with units of 
!> 'seconds since 1966-01-01'  (the start of the satellite era)

function read_model_time(filename)

character(len=*), intent(in) :: filename !< file to get time
type(time_type) :: read_model_time !< returned time from file

character(len=*), parameter :: routine = 'read_model_time'

integer         :: ncid, VarID, rc
real(digits12)  :: length_of_run
type(time_type) :: model_time_base, model_time_offset
integer         :: iyear, imonth, iday, isecond
integer(i8)     :: ibigday
character(len=NF90_MAX_NAME) :: unitstring

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR, routine, string1, source, revision, revdate)
endif

rc = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(rc, routine, 'open '//filename )

rc = nf90_inq_varid(ncid, 'time', VarID)
call nc_check(rc, routine, 'time inq_varid')

! read the 64 bit real and form the number of days and number of
! seconds - both (32bit) integers. the 64bit real could overflow
! a single 32bit integer, so must use 'long' integers

rc = nf90_get_var(ncid, VarID, length_of_run)
call nc_check(rc, routine, 'time get_var')

ibigday = int(length_of_run,i8)/86400_i8
iday    = int(ibigday)
isecond = int(length_of_run - real(iday,digits12)*86400.0_digits12)

model_time_offset = set_time(isecond, iday)

! read the units ... and check for valid assumption
rc = nf90_get_att(ncid, VarID, 'units', unitstring)
call nc_check(rc, routine,'get_att time:units')

if (unitstring(1:14) == 'seconds since ') then
   read(unitstring,'(14x,i4,2(1x,i2))',iostat=rc)iyear,imonth,iday
   if (rc /= 0) then
      write(string1,*)'Unable to parse year/month/day from units. Error status was ',rc
      write(string2,*)'expected something like "seconds since YYYY-MM-DD"'
      write(string3,*)'was                     "'//trim(unitstring)//'"'
      call error_handler(E_ERR, routine, string1, &
          source, revision, revdate, text2=string2, text3=string3)
   endif
else
   write(string1,*)'Unable to read time units.'
   write(string2,*)'expected "seconds since YYYY-MM-DD"'
   write(string3,*)'was      "'//trim(unitstring)//'"'
   call error_handler(E_ERR, routine, string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

model_time_base = set_date(iyear, imonth, iday)

read_model_time = model_time_base + model_time_offset

call print_date(read_model_time,'read_model_time:netcdf model date')
call print_time(read_model_time,'read_model_time:DART   model time')

end function read_model_time

!--------------------------------------------------------------------
!> write the current model time to an openggcm output file
!> the time is encoded as a 64bit real with units of 
!> 'seconds since 1966-01-01'  (the start of the satellite era)

subroutine write_model_time(ncid, dart_time)

integer, intent(in) :: ncid    !< file to receive time
type(time_type) :: dart_time   !< model time

character(len=*), parameter :: routine = 'write_model_time'

integer         :: DimID, VarID, rc
type(time_type) :: model_time_base, model_time_offset
real(digits12)  :: rbigday, rsecond
integer         ::    iday, isecond
integer         :: ntimes

if ( .not. module_initialized ) call static_init_model

rc = nf90_Redef(ncid)
call nc_check(rc, routine, "redef")

rc = nf90_inq_dimid(ncid, "time", DimID)
if (rc == NF90_NOERR) then
   rc = nf90_inquire_dimension(ncid, DimID, len=ntimes)
   call nc_check(rc, routine, "inquire_dimension time")
else
   rc = nf90_def_dim(ncid, "time", NF90_UNLIMITED, DimID)
   call nc_check(rc, routine, "def_dim time")
   ntimes = 1
endif

rc = nf90_inq_varid(ncid, "time", VarID)
if (rc /= NF90_NOERR) then

   rc = nf90_def_var(ncid,name="time",xtype=nf90_double,dimids=(/DimID/),varid=VarID)
   call nc_check(rc, routine, "def_var time")

   rc = nf90_put_att(ncid, VarID, "long_name", "valid time of the model state")
   call nc_check(rc, routine, "long_name time")

   rc = nf90_put_att(ncid, VarID, "units", "seconds since 1966-01-01 00:00:00")
   call nc_check(rc, routine, "units long_name")

endif

call nc_check( nf90_Enddef(ncid),routine, "Enddef" )

! read the 64 bit real and form the number of days and number of
! seconds - both (32bit) integers. the 64bit real could overflow
! a single 32bit integer, so must use 'long' integers

model_time_base   = set_date(1966, 1, 1)
model_time_offset = dart_time - model_time_base

call get_time(model_time_offset, isecond, iday)

rbigday = real(   iday,digits12) * 86400.0_digits12
rsecond = real(isecond,digits12) + rbigday

rc = nf90_put_var(ncid, VarID, rsecond, start=(/ ntimes /))
call nc_check(rc, routine, "put_var model_time")

end subroutine write_model_time


!--------------------------------------------------------------------
!> pass the vertical localization coordinate to assim_tools_mod

function query_vert_localization_coord()

integer :: query_vert_localization_coord !< return which height we want to 
                                         !< localize in

query_vert_localization_coord = VERTISHEIGHT

end function query_vert_localization_coord

!------------------------------------------------------------

!> Always an error to call this routine.
!> At present, this is only called if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> This is not possible for a large geophysical model.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:) !< initalize state from scratch

string2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
string3 = 'use openggcm_to_dart to generate an initial state'
call error_handler(E_ERR,'init_conditions', &
                  'ERROR!!  openggcm model has no built-in default state', &
                  source, revision, revdate, &
                  text2=string2, text3=string3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
x = 0.0_r8

end subroutine init_conditions

!------------------------------------------------------------------

!> If the model could be called as a subroutine, does a single
!> timestep advance.  openggcm cannot be called this way, so fatal error
!> if this routine is called.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:) !< state vector
type(time_type), intent(in)    :: time !< time stamp for state_vector

call error_handler(E_ERR,'adv_1step', &
                  'openggcm model cannot be called as a subroutine; async cannot = 0', &
                  source, revision, revdate)

end subroutine adv_1step

!------------------------------------------------------------------

!> Companion interface to init_conditions. Returns a time that is
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_time(time)

type(time_type), intent(out) :: time !< time restart file

string2 = "cannot run perfect_model_obs with 'start_from_restart = .false.' "
string3 = 'use openggcm_to_dart to generate an initial state which contains a timestamp'
call error_handler(E_ERR,'init_time', &
                  'ERROR!!  openggcm model has no built-in default time', &
                  source, revision, revdate, &
                  text2=string2, text3=string3)

! this code never reached - just here to avoid compiler warnings
! about an intent(out) variable not being set to a value.
time = set_time(0,0)

end subroutine init_time

!------------------------------------------------------------------

!> Verify that the namelist was filled in correctly, and check
!> that there are valid entries for the dart_kind. 
!> Returns a table with columns:  
!>
!> netcdf_variable_name ; dart_kind_string ; update_string ; grid_id
!>

subroutine verify_state_variables( state_variables, ngood, table, kind_list, update_var, grid_id )

character(len=*), intent(inout) :: state_variables(:) !< list of state variables and attributes from nml
integer,          intent(out)   :: ngood              !< number of good namelist values or nfields
character(len=*), intent(out)   :: table(:,:)         !< 2d table with information from namelist
integer,          intent(out)   :: kind_list(:)       !< dart kind
logical,          intent(out)   :: update_var(:)      !< list of logical update information
integer,          intent(out)   :: grid_id(:)         !< list of dart kind numbers

! Local Variables
integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update, gridname

nrows = size(table,1)

ngood = 0

if ( state_variables(1) == ' ' ) then ! no model_state_variables namelist provided
   string1 = 'model_nml:model_state_variables not specified using default variables'
   call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
endif

MyLoop : do i = 1, nrows

   varname  = trim(state_variables(4*i-3))
   dartstr  = trim(state_variables(4*i-2))
   update   = trim(state_variables(4*i-1))
   gridname = trim(state_variables(4*i  ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)
   table(i,4) = trim(gridname)

   if ( table(i,1) == ' ' .and. &
        table(i,2) == ' ' .and. &
        table(i,3) == ' ' .and. &
        table(i,4) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. &
        table(i,2) == ' ' .or. &
        table(i,3) == ' ' .or. &
        table(i,4) == ' ') then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   kind_list(i) = get_index_for_quantity(dartstr)
   if( kind_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif
   
   ! Make sure the update variable has a valid name

   SELECT CASE (update)
      CASE ('UPDATE')
         update_var(i) = .true.
      CASE ('NO_COPY_BACK')
         update_var(i) = .false.
      CASE DEFAULT
         write(string1, '(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update), ', ', trim(gridname)
         call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate, text2=string2)
   END SELECT

   ! Make sure the update variable has a valid name

   SELECT CASE (gridname)
      CASE ('GEOGRAPHIC_GRID')
         grid_id(i) = GEOGRAPHIC_GRID
      CASE ('MAGNETIC_GRID')
         grid_id(i) = MAGNETIC_GRID
      CASE DEFAULT
         write(string1, '(A)')  'only GEOGRAPHIC_GRID or MAGNETIC_GRID supported in model_state_variable namelist'
         write(string2,'(8A)') 'you provided : ',&
                               trim(varname), ', ', trim(dartstr), ', ', trim(update), ', ', trim(gridname)
         call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate, text2=string2)
   END SELECT

   ! Record the contents of the DART state vector

   if (do_output()) then
      write(string1,'(A,I2,8A)') 'variable ',i,' is ',trim(varname), ', ', &
                                  trim(dartstr), ', ', trim(update), ', ', trim(gridname)
      call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate)
   endif

   ngood = ngood + 1
enddo MyLoop

end subroutine verify_state_variables

!----------------------------------------------------------------------

!> determine if grid is geographic or magnetic given a dart kind
!> this is based on what the user tells us in the namelist.

function get_grid_type(dart_kind)

integer, intent(in) :: dart_kind !< dart kind
integer :: get_grid_type         !< returned grid

! Local Variables
integer :: i

do i = 1,nfields
   if ( dart_kind == state_kinds_list(i) ) then
      get_grid_type = grid_info_list(i)
      return
   endif
enddo

write(string1,*)' Can not find grid type for : ', get_name_for_quantity(dart_kind)
call error_handler(E_ERR,'get_grid_type',string1,source,revision,revdate)

end function get_grid_type

!----------------------------------------------------------------------

!> Allocate space for grid variables. 
!> cannot be called until the grid sizes are set.

subroutine allocate_grid_space(grid_handle, conv)

type(grid_type), intent(inout) :: grid_handle !< geo or mag grid handle
logical,         intent(in)    :: conv !< if true, the grid has conversion arrays

!>@todo use dim_order_list ... maybe ... at least check for consistency with variable shape
! so that when we use the indices for the variable, we are sure that the indices will work
! on the coordinate arrays

allocate(grid_handle%longitude(grid_handle%nheight, grid_handle%nlat, grid_handle%nlon))
allocate(grid_handle%latitude( grid_handle%nheight, grid_handle%nlat, grid_handle%nlon))
allocate(grid_handle%height(   grid_handle%nheight, grid_handle%nlat, grid_handle%nlon))

!>@todo check if these are nlon,nlat or nlat,nlon ...
if (conv) then
   allocate(grid_handle%conv_2d_lon(grid_handle%nlon, grid_handle%nlat))
   allocate(grid_handle%conv_2d_lat(grid_handle%nlon, grid_handle%nlat))
endif

end subroutine allocate_grid_space

!----------------------------------------------------------------------

!> Deallocate space for grid variables. 

subroutine deallocate_grid_space(grid_handle)

type(grid_type), intent(inout) :: grid_handle !< geo or mag grid handle

if (allocated(grid_handle%longitude))  deallocate(grid_handle%longitude)
if (allocated(grid_handle%latitude))   deallocate(grid_handle%latitude)
if (allocated(grid_handle%height))     deallocate(grid_handle%height)

if (allocated(grid_handle%conv_2d_lon)) deallocate(grid_handle%conv_2d_lon)
if (allocated(grid_handle%conv_2d_lat)) deallocate(grid_handle%conv_2d_lat)

end subroutine deallocate_grid_space

!------------------------------------------------------------------

!>@todo FIXME: the following routines should be in a netcdf utils 
!> module somewhere

!------------------------------------------------------------------

!< returns the length of the dimension from a given a dimension name
!< and netcdf file id

function get_dim(ncid, dim_name, context)

integer,          intent(in)    :: ncid !< netcdf file id
character(len=*), intent(in)    :: dim_name !< dimension name of interest
character(len=*), intent(in)    :: context !< routine calling from
integer :: get_dim !< returns the length of the dimension

! netcdf variables
integer :: DimID, rc

rc = NF90_inq_dimid(ncid=ncid, name=dim_name, dimid=DimID)
call nc_check(rc, trim(context)//' inquiring for dimension '//trim(dim_name))

rc = NF90_inquire_dimension(ncid, DimID, len=get_dim)
call nc_check(rc, trim(context)//' getting length of dimension '//trim(dim_name))

end function get_dim

!------------------------------------------------------------------

!< sets a netcdf file dimension given a dimension name and length.
!< returns the value of the netcdf dimension id.

function set_dim(ncid, dim_name, dim_val, context)

integer,          intent(in)    :: ncid !< netcdf id
character(len=*), intent(in)    :: dim_name !< dimension name
integer,          intent(in)    :: dim_val !< dimension length
character(len=*), intent(in)    :: context !< routine calling from
integer :: set_dim !< netcdf id to the defined dimension


! netcdf variables
integer :: rc

rc = NF90_def_dim(ncid=ncid, name=dim_name, len=dim_val, dimid=set_dim)
call nc_check(rc, trim(context)//' setting dimension '//trim(dim_name))

end function set_dim

!------------------------------------------------------------------

!> initialize a grid transformation type.  the transform type
!> and conversion routines are in the cort_mod.  the grid transformations
!> change with time, so the initialization routine must know the
!> current model time.

subroutine initialize_openggcm_transform(filename)

character(len=*), intent(inout) :: filename !< file with openggcm time information


! Local Variables
type(time_type) :: dart_time
integer :: yr,mo,dy,hr,mn,se

dart_time = read_model_time(filename)

call get_date(dart_time, yr, mo, dy, hr, mn, se)

call cotr_set(yr,mo,dy,hr,mn,real(se,r4),openggcm_transform)

end subroutine initialize_openggcm_transform

!----------------------------------------------------------------------

!> transform grid from geo->magnetic or vice-versa. 
!> currently we only need it to translate from magnetic to geo.
!> direction options GEO_TO_SM or SM_TO_GEO

subroutine transform_mag_geo(llon, llat, lheight, direction)

real(r8), intent(inout) :: llon !< dart longitude [0,360]
real(r8), intent(inout) :: llat !< dart latitude [-90,90]
real(r8), intent(inout) :: lheight !< height
integer,  intent(in)    :: direction


! Local Variables
real(r4) :: xin, xout, yin, yout, zin, zout
character(len=4) :: dirstrin, dirstrout

if (direction == SM_TO_GEO) then
   dirstrin  = 'sm '
   dirstrout = 'geo'
else if (direction == GEO_TO_SM) then
   dirstrin  = 'geo'
   dirstrout = 'sm '
else
  call error_handler(E_ERR, 'transform_mag_geo', &
          'unexpected direction input, should not happen', &
           source, revision, revdate)
endif

! cort_mod is expecting height from center of earth
if (lheight == MISSING_R8) then
   lheight = 0.0_r8 + earth_radius
else
   lheight = lheight + earth_radius
endif

! degxyz requires longitude to be between -180 and 180
if (llon > 180.0_r8) llon = llon - 360.0_r8
! degxyz requires latitude to be between 0 and 180
llat = 90.0_r8 - llat

! transform spherical coordinates to cartesian
call degxyz(lheight, llon, llat, xin, yin, zin)

! transform from geographic to magnetic grid or back
call cotr(openggcm_transform, dirstrin, dirstrout, &
          xin, yin, zin, xout, yout, zout)

! transform cartesian coordinates to spherical 
call xyzdeg(xout, yout, zout, lheight, llon, llat)

! transform back to dart longitude coordinates [0,360]
if (llon < 0.0_r8) llon = llon + 360.0_r8
! transform back to dart latitude coordinates [-90,90]
llat = 90.0_r8 - llat

end subroutine transform_mag_geo

!----------------------------------------------------------------------

!> recording the storage order of the dimensions for each variable.
!>
!> hgt_index is   dim_order_list(VAR_ID, VAR_HGT_INDEX)
!> lat_index is   dim_order_list(VAR_ID, VAR_LAT_INDEX)
!> lon_index is   dim_order_list(VAR_ID, VAR_LON_INDEX)
!>
!> variables without a height dimension, VAR_HGT_INDEX is set to one.

subroutine make_dim_order_table(ngood)
integer, intent(in) :: ngood !< number of good fields

! Local Variables
integer :: ivar, jdim
character(len=NF90_MAX_NAME) :: dimname

! initialize list
dim_order_list(:,:) = 1

do ivar = 1,ngood
   do jdim = 1,get_num_dims(domain_id, ivar)
      dimname = get_dim_name(domain_id, ivar, jdim)
      SELECT CASE (trim(dimname))
         CASE ('cg_height','ig_height','geo_height')
            dim_order_list(ivar, VAR_HGT_INDEX) = jdim
         CASE ('cg_lat','ig_lat','geo_lat')
            dim_order_list(ivar, VAR_LAT_INDEX) = jdim
         CASE ('cg_lon','ig_lon','geo_lon')
            dim_order_list(ivar, VAR_LON_INDEX) = jdim
         CASE DEFAULT
            write(string1,*) 'cannot find dimension ', trim(dimname),&
                               ' for variable', get_variable_name(domain_id, ivar)
            call error_handler(E_ERR,'make_dim_order_table',string1,source,revision,revdate)
      END SELECT
   enddo
enddo

!>@todo what about falling off the list? 

end subroutine make_dim_order_table

!----------------------------------------------------------------------
!------------------------------------------------------------------
! End of model_mod
!------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
