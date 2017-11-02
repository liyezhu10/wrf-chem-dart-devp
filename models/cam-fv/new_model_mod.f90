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
use  netcdf_utilities_mod
use       location_io_mod
use     default_model_mod

use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

! routines in this list have code in this module
public :: get_model_size,                &
          get_state_meta_data,           &
          static_init_model,             &
          model_interpolate,             &
          shortest_time_between_assimilations, &
          convert_vertical_obs,          &
          convert_vertical_state,        &
          end_model,                     &
          nc_write_model_atts,           &
          write_model_time,              &
          read_model_time

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

character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

! things which can/should be in the model_nml
integer  :: assimilation_period_days     = 0
integer  :: assimilation_period_seconds  = 21600
integer  :: vert_localization_coord      = VERTISPRESSURE
integer  :: debug = 0   ! turn up for more and more debug messages
logical  :: minimal_output = .false.
character(len=256) :: cam_grid_filename = 'caminput.nc'
character(len=256) :: cam_phis_filename = 'camphis.nc'

namelist /model_nml/  &
   assimilation_period_days,    &
   assimilation_period_seconds, &
   cam_grid_filename,           &
   cam_phis_filename,           &
   vert_localization_coord,     &
   minimal_output,              &
   debug,                       &
   variables

! DART contents are specified in the input.nml:&model_nml namelist.
!>@todo  NF90_MAX_NAME is 256 ... this makes the namelist output unreadable
integer, parameter :: MAX_STATE_VARIABLES = 8
integer, parameter :: num_state_table_columns = 5
character(len=vtablenamelength) :: variables(MAX_STATE_VARIABLES * num_state_table_columns ) = ' '
character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
logical  ::                   update_list(MAX_STATE_VARIABLES) = .FALSE.
integer  ::                     kind_list(MAX_STATE_VARIABLES) = MISSING_I
real(r8) ::                    clamp_vals(MAX_STATE_VARIABLES,2) = MISSING_R8

integer :: nfields   ! This is the number of variables in the DART state vector.

integer :: domain_id ! global variable for state_structure_mod routines

!> Everything needed to describe a variable. Basically all the metadata from
!> a netCDF file is stored here as well as all the information about where
!> the variable is stored in the DART state vector.
!>

type cam_grid
   integer :: nlon, nlat
   integer :: nslon, nslat, 
   integer :: nlevels   ! is a vertical stagger possible?
   real(r8), allocatable :: lon_vals(:)
   real(r8), allocatable :: lat_vals(:)
   ! something for vertical = hyma, hymb?
end type

type(time_type) :: model_timestep

integer(i8) :: model_size    ! the state vector length

type(cam_grid) :: grid

contains


!-----------------------------------------------------------------------
! All the REQUIRED interfaces come first - by convention.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Returns the size of the DART state vector (i.e. model) as an integer.
!> Required for all applications.
!>

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

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

select case (myqty)
  case U_WIND_COMPONENT:
   location = set_location(grid%slon(iloc), grid%lat(jloc), grid%level(vloc), VERTISLEVEL)

  case V_WIND_COMPONENT:
   location = set_location(grid%lon(iloc), grid%slat(jloc), grid%level(vloc), VERTISLEVEL)

  case default
   location = set_location(grid%lon(iloc), grid%lat(jloc), grid%level(vloc), VERTISLEVEL)

end select

! return state quantity for this index if requested
if (present(var_type)) var_type = myqty

end subroutine get_state_meta_data


!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> (i.e. S, T, U, V, Eta) to the given location given a state vector.
!> The type of the variable being interpolated is obs_type since
!> normally this is used to find the expected value of an observation
!> at some location. The interpolated value is returned in interp_vals
!> and istatus is 0 for success. NOTE: This is a workhorse routine and is
!> the basis for all the forward observation operator code.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_type the DART KIND of interest
!> @param interp_val the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_type, interp_val, istatus)

 type(ensemble_type), intent(in) :: state_handle
 integer,             intent(in) :: ens_size
 type(location_type), intent(in) :: location
 integer,             intent(in) :: obs_type
 integer,            intent(out) :: istatus(ens_size)
 real(r8),           intent(out) :: interp_val(ens_size) !< array of interpolated values

if ( .not. module_initialized ) call static_init_model

! Successful istatus is 0
interp_val = MISSING_R8
istatus = 99

write(string1,*)'model_interpolate should not be called.'
write(string2,*)'we are getting forward observations directly from ROMS'
call error_handler(E_MSG,'model_interpolate:',string1,source,revision,revdate, text2=string2)

end subroutine model_interpolate


!-----------------------------------------------------------------------
!>
!> Returns the the time step of the model; the smallest increment in
!> time that the model is capable of advancing the ROMS state.
!>

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = model_timestep

end function shortest_time_between_assimilations


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
integer :: ss, dd
integer :: ncid

character(len=32) :: calendar

type(time_type) :: model_time

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


! put this in a subroutine that deals with time
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd)
write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model:',string1,source,revision,revdate)

call set_calendar_type( trim(calendar) )


! put this in a subroutine that deals with the grid
call nc_check( nf90_open(trim(roms_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(roms_filename))

! Get the grid info
call get_grid_dimensions()
call get_grid()

call get_time_information(roms_filename, ncid, 'ocean_time', 'ocean_time', &
                          calendar=calendar, last_time=model_time)

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(roms_filename))


! is there a common subroutine outside of the model mod we can call here?

! parse_variable_input() fills var_names, kind_list, clamp_vals, update_list
call parse_variable_input(variables, nfields)

domain_id = add_domain(roms_filename, nfields, &
                    var_names, kind_list, clamp_vals, update_list )

if (debug > 2) call state_structure_info(domain_id)



model_size = get_domain_size(domain_id)

end subroutine static_init_model


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
!> @param ncid the netCDF handle of the DART diagnostic file opened by
!>                 assim_model_mod:init_diag_output
!> @param model_writes_state have the state structure write out all of the
!>                 state variables

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

! for the dimensions and coordinate variables
integer :: nxirhoDimID, nxiuDimID, nxivDimID
integer :: netarhoDimID, netauDimID, netavDimID
integer :: nsrhoDimID, nswDimID
integer :: VarID

! local variables

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

! we only have a netcdf handle here so we do not know the filename
! or the fortran unit number.  but construct a string with at least
! the netcdf handle, so in case of error we can trace back to see
! which netcdf file is involved.

write(filename,*) 'ncid', ncid

! Write Global Attributes

!>@todo FIXME  make writing the grid info optional.
!> based on a namelist setting.  if not writing grid,
!> this routine has nothing to do.

if (minimal_output) return

! add grid info
call nc_redef(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source)
call nc_add_global_attribute(ncid, "model_revision", revision)
call nc_add_global_attribute(ncid, "model_revdate", revdate)

call nc_add_global_attribute(ncid, "model", "ROMS")

! We need to output the grid information
! Define the new dimensions IDs

call nc_check(nf90_def_dim(ncid, name='xi_rho',  len = Nxi_rho, &
     dimid = nxirhoDimID),'nc_write_model_atts', 'xi_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_rho', len = Neta_rho,&
     dimid = netarhoDimID),'nc_write_model_atts', 'eta_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='s_rho',   len = Ns_rho,&
     dimid = nsrhoDimID),'nc_write_model_atts', 's_rho def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='s_w',   len = Ns_w,&
     dimid = nswDimID),'nc_write_model_atts', 's_w def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='xi_u',    len = Nxi_u,&
     dimid = nxiuDimID),'nc_write_model_atts', 'xi_u def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='xi_v',    len = Nxi_v,&
     dimid = nxivDimID),'nc_write_model_atts', 'xi_v def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_u',   len = Neta_u,&
     dimid = netauDimID),'nc_write_model_atts', 'eta_u def_dim '//trim(filename))

call nc_check(nf90_def_dim(ncid, name='eta_v',   len = Neta_v,&
     dimid = netavDimID),'nc_write_model_atts', 'eta_v def_dim '//trim(filename))

! Create the Coordinate Variables and give them Attributes
! The values will be added in a later block of code.

call nc_check(nf90_def_var(ncid,name='lon_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho longitudes'), &
              'nc_write_model_atts', 'lon_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho latitudes'), &
              'nc_write_model_atts', 'lat_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lon_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u longitudes'), &
              'nc_write_model_atts', 'lon_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u latitudes'), &
              'nc_write_model_atts', 'lat_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lon_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
              'nc_write_model_atts', 'lon_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v longitudes'), &
              'nc_write_model_atts', 'lon_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'lon_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='lat_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
              'nc_write_model_atts', 'lat_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v latitudes'), &
              'nc_write_model_atts', 'lat_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
              'nc_write_model_atts', 'lat_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_rho', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_rho def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_rho long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_rho units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_u', xtype=nf90_double, &
              dimids=(/ nxiuDimID, netauDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_u def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_u long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_u units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_v', xtype=nf90_double, &
              dimids=(/ nxivDimID, netavDimID, nsrhoDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_v def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_v long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_v units '//trim(filename))

call nc_check(nf90_def_var(ncid,name='z_w', xtype=nf90_double, &
              dimids=(/ nxirhoDimID, netarhoDimID, nswDimID /), varid=VarID),&
              'nc_write_model_atts', 'z_w def_var '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
              'nc_write_model_atts', 'z_w long_name '//trim(filename))
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
              'nc_write_model_atts', 'z_w units '//trim(filename))

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_enddef(ncid)

! Fill the coordinate variable values

! the RHO grid

call nc_check(NF90_inq_varid(ncid, 'lon_rho', VarID), &
              'nc_write_model_atts', 'lon_rho inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, TLON ), &
             'nc_write_model_atts', 'lon_rho put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_rho', VarID), &
              'nc_write_model_atts', 'lat_rho inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, TLAT ), &
             'nc_write_model_atts', 'lat_rho put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'z_rho', VarID), &
              'nc_write_model_atts', 'z_rho inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, TDEP ), &
             'nc_write_model_atts', 'z_rho put_var '//trim(filename))

! the U grid

call nc_check(NF90_inq_varid(ncid, 'lon_u', VarID), &
              'nc_write_model_atts', 'lon_u inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, ULON ), &
             'nc_write_model_atts', 'lon_u put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_u', VarID), &
              'nc_write_model_atts', 'lat_u inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, ULAT ), &
             'nc_write_model_atts', 'lat_u put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'z_u', VarID), &
              'nc_write_model_atts', 'z_u inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, UDEP ), &
             'nc_write_model_atts', 'z_u put_var '//trim(filename))

! the V grid

call nc_check(NF90_inq_varid(ncid, 'lon_v', VarID), &
              'nc_write_model_atts', 'lon_v inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, VLON ), &
             'nc_write_model_atts', 'lon_v put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'lat_v', VarID), &
              'nc_write_model_atts', 'lat_v inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, VLAT ), &
             'nc_write_model_atts', 'lat_v put_var '//trim(filename))

call nc_check(NF90_inq_varid(ncid, 'z_v', VarID), &
              'nc_write_model_atts', 'z_v inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, VDEP ), &
             'nc_write_model_atts', 'z_v put_var '//trim(filename))

! the W grid

call nc_check(NF90_inq_varid(ncid, 'z_w', VarID), &
              'nc_write_model_atts', 'z_w inq_varid '//trim(filename))
call nc_check(nf90_put_var(ncid, VarID, WDEP ), &
             'nc_write_model_atts', 'z_w put_var '//trim(filename))

! Flush the buffer and leave netCDF file open
call nc_sync(ncid)


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

if (present(adv_to_time)) then
   string3 = time_to_string(adv_to_time)
   write(string1,*)'ROMS/DART not configured to advance ROMS.'
   write(string2,*)'called with optional advance_to_time of'
   call error_handler(E_ERR, 'write_model_time', string1, &
              source, revision, revdate, text2=string2,text3=string3)
endif

! If the ocean_time variable exists, we are updating a ROMS file,
! if not ... must be updating a DART diagnostic file.

io = nf90_inq_varid(ncid,'ocean_time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'ocean_time', 'ocean_time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
endif

io = nf90_inq_varid(ncid,'time',varid)
if (io == NF90_NOERR) then
   call get_time_information('unknown', ncid, 'time', 'time', &
                myvarid=varid, origin_time=origin_time)
   deltatime = model_time - origin_time
   call get_time(deltatime, seconds, days)
   run_duration = real(days,digits12)*86400.0_digits12 + real(seconds,digits12)
   call nc_check(nf90_put_var(ncid, varid, run_duration), 'write_model_time', 'put_var')
   return
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

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'read_model_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))

call get_time_information(filename, ncid, 'ocean_time', 'ocean_time', last_time=read_model_time)

call nc_check( nf90_close(ncid), 'read_model_time', 'close '//trim(filename))

end function read_model_time



!-----------------------------------------------------------------------
! The remaining (private) interfaces come last.
! None of the private interfaces need to call static_init_model()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>
!> Set the desired minimum model advance time. This is generally NOT the
!> dynamical timestep of the model, but rather the shortest forecast length
!> you are willing to make. This impacts how frequently the observations
!> may be assimilated.
!>

function set_model_time_step()

type(time_type) :: set_model_time_step

! assimilation_period_seconds, assimilation_period_days are from the namelist

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!-----------------------------------------------------------------------
!>
!> Read the grid dimensions from the ROMS grid netcdf file.
!> By reading the dimensions first, we can use them in variable
!> declarations later - which is faster than using allocatable arrays.
!>

subroutine get_grid_dimensions()

integer :: ncid

! Read the (static) grid dimensions from the ROMS grid file.

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
              'get_grid_dimensions', 'open '//trim(roms_filename))

Nxi_rho   = get_dimension_length(ncid, 'xi_rho',   roms_filename)
Nxi_u     = get_dimension_length(ncid, 'xi_u',     roms_filename)
Nxi_v     = get_dimension_length(ncid, 'xi_v',     roms_filename)
Neta_rho  = get_dimension_length(ncid, 'eta_rho',  roms_filename)
Neta_u    = get_dimension_length(ncid, 'eta_u',    roms_filename)
Neta_v    = get_dimension_length(ncid, 'eta_v',    roms_filename)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_filename))

! Read the vertical dimensions from the dedicated file.

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
               'get_grid_dimensions', 'open '//trim(roms_filename))

Ns_rho    = get_dimension_length(ncid, 's_rho',    roms_filename)
Ns_w      = get_dimension_length(ncid, 's_w'  ,    roms_filename)

call nc_check(nf90_close(ncid), &
              'get_grid_dimensions','close '//trim(roms_filename))

Nx =  Nxi_rho  ! Setting the nominal value of the 'global' variables
Ny = Neta_rho  ! Setting the nominal value of the 'global' variables
Nz =   Ns_rho  ! Setting the nominal value of the 'global' variables

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
!>
!> Read the actual grid values from the ROMS netcdf file.
!>
!>@todo FIXME:  the original implementation opened 3 different files
!> to get the grid info - the namelist was:
!>    roms_ini_filename            = '../data/wc13_ini.nc'
!>    grid_definition_filename     = '../data/wc13_grd.nc'
!>    depths_definition_filename   = '../data/wc13_depths.nc'
!>
!> these have been consolidated by hernan for the santa cruz version
!> into a single file.  check with the other rutgers folks to see if
!> they still need to open 3 different files.  if so, we might need
!> to restore the 3 namelist items and we can use the same file for
!> all 3 types of grid info in the first case, and 3 different files
!> for the second case.
!>

subroutine get_grid()

integer  :: ncid, VarID

real(r8), parameter :: all_land = 0.001_r8

if (.not. allocated(ULAT)) allocate(ULAT(Nxi_u, Neta_u))
if (.not. allocated(ULON)) allocate(ULON(Nxi_u, Neta_u))
if (.not. allocated(UDEP)) allocate(UDEP(Nxi_u, Neta_u, Nz))

if (.not. allocated(VLAT)) allocate(VLAT(Nxi_v, Neta_v))
if (.not. allocated(VLON)) allocate(VLON(Nxi_v, Neta_v))
if (.not. allocated(VDEP)) allocate(VDEP(Nxi_v, Neta_v, Nz))

if (.not. allocated(TLAT)) allocate(TLAT(Nxi_rho, Neta_rho))
if (.not. allocated(TLON)) allocate(TLON(Nxi_rho, Neta_rho))
if (.not. allocated(TDEP)) allocate(TDEP(Nxi_rho, Neta_rho, Nz))
if (.not. allocated(WDEP)) allocate(WDEP(Nxi_rho, Neta_rho, Ns_w))

! Read the vertical information from the (separate) roms_filename

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_u', VarID), &
      'get_grid', 'inq_varid z_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, UDEP), &
      'get_grid', 'get_var z_u '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_w', VarID), &
      'get_grid', 'inq_varid z_w '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, WDEP), &
      'get_grid', 'get_var z_w '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_v', VarID), &
      'get_grid', 'inq_varid z_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VDEP), &
      'get_grid', 'get_var z_v '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'z_rho', VarID), &
      'get_grid', 'inq_varid z_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TDEP), &
      'get_grid', 'get_var z_rho '//trim(roms_filename))

call nc_check(nf90_close(ncid), &
             'get_var','close '//trim(roms_filename))

! Read the rest of the grid information from the traditional grid file

call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
      'get_grid', 'open '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_rho', VarID), &
   'get_grid', 'inq_varid lon_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TLON), &
      'get_grid', 'get_var lon_rho '//trim(roms_filename))

where (TLON < 0.0_r8) TLON = TLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_rho', VarID), &
      'get_grid', 'inq_varid lat_rho '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, TLAT), &
      'get_grid', 'get_var lat_rho '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_u', VarID), &
      'get_grid', 'inq_varid lon_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, ULON), &
      'get_grid', 'get_var lon_u '//trim(roms_filename))

where (ULON < 0.0_r8) ULON = ULON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_u', VarID), &
      'get_grid', 'inq_varid lat_u '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, ULAT), &
      'get_grid', 'get_var lat_u '//trim(roms_filename))

call nc_check(nf90_inq_varid(ncid, 'lon_v', VarID), &
      'get_grid', 'inq_varid lon_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VLON), &
      'get_grid', 'get_var lon_v '//trim(roms_filename))

where (VLON < 0.0_r8) VLON = VLON + 360.0_r8

call nc_check(nf90_inq_varid(ncid, 'lat_v', VarID), &
      'get_grid', 'inq_varid lat_v '//trim(roms_filename))
call nc_check(nf90_get_var( ncid, VarID, VLAT), &
      'get_grid', 'get_var lat_v '//trim(roms_filename))

! Be aware that all the depths are negative values.
! The surface of the ocean is 0.0, the deepest is a big negative value.

if (do_output() .and. debug > 0) then
    write(string1,*)'    min/max ULON ',minval(ULON), maxval(ULON)
    write(string2,*)    'min/max ULAT ',minval(ULAT), maxval(ULAT)
    write(string3,*)    'min/max UDEP ',minval(UDEP), maxval(UDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)

    write(string1,*)'    min/max VLON ',minval(VLON), maxval(VLON)
    write(string2,*)    'min/max VLAT ',minval(VLAT), maxval(VLAT)
    write(string3,*)    'min/max VDEP ',minval(VDEP), maxval(VDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)

    write(string1,*)'    min/max TLON ',minval(TLON), maxval(TLON)
    write(string2,*)    'min/max TLAT ',minval(TLAT), maxval(TLAT)
    write(string3,*)    'min/max TDEP ',minval(TDEP), maxval(TDEP)
    call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)
endif

end subroutine get_grid


!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!>
!>@param state_variables the list of variables and kinds from model_mod_nml
!>@param ngood the number of variable/KIND pairs specified

subroutine parse_variable_input( state_variables, ngood )

character(len=*), intent(in)  :: state_variables(:)
integer,          intent(out) :: ngood

integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 5   change to updateable

ngood = 0
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname      = trim(state_variables(num_state_table_columns*i-4))
   dartstr      = trim(state_variables(num_state_table_columns*i-3))
   minvalstring = trim(state_variables(num_state_table_columns*i-2))
   maxvalstring = trim(state_variables(num_state_table_columns*i-1))
   state_or_aux = trim(state_variables(num_state_table_columns*i  ))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or. dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
      call error_handler(E_ERR,'parse_variable_input:',string1,source,revision,revdate)
   endif

   ! Make sure DART kind is valid

   if( get_index_for_quantity(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_input:',string1,source,revision,revdate)
   endif

   call to_upper(minvalstring)
   call to_upper(maxvalstring)
   call to_upper(state_or_aux)

   var_names(   i) = varname
   kind_list(   i) = get_index_for_quantity(dartstr)
   clamp_vals(i,1) = string_to_real(minvalstring)
   clamp_vals(i,2) = string_to_real(maxvalstring)
   update_list( i) = string_to_logical(state_or_aux, 'UPDATE')

   ngood = ngood + 1

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: There is a possibility you need to increase ''MAX_STATE_VARIABLES'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_input:',string1,source,revision,revdate,text2=string2)
endif

end subroutine parse_variable_input


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
