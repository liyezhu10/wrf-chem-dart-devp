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
use  netcdf_utilities_mod,  only : nc_check, nc_get_variable, nc_get_variable_size
use       location_io_mod
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

call read_grid_info(cam_template_filename, cam_phis_filename, grid_data)

! is there a common subroutine outside of the model mod we can call here?

! set_cam_variable_info() fills var_names, kind_list, clamp_vals, update_list
! from the &model_mod_nml variables
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

select case (grid_stagger%qty_stagger(myqty))
  case (STAGGER_U)
   location = set_location(grid_data%lon%vals(iloc), &
                           grid_data%slat%vals(jloc), &
                           real(vloc,r8), VERTISLEVEL)

  case (STAGGER_V)
   location = set_location(grid_data%slon%vals(iloc), &
                           grid_data%lat%vals(jloc), &
                           real(vloc,r8), VERTISLEVEL)
   
  !>@todo not sure what to do yet. ? +-1/2 ?
  case (STAGGER_W)
   location = set_location(grid_data%lon%vals(iloc), &
                           grid_data%lat%vals(jloc), &
                           real(vloc,r8), VERTISLEVEL)
  case default
   location = set_location(grid_data%lon%vals(iloc), &
                           grid_data%lat%vals(jloc), &
                           real(vloc,r8), VERTISLEVEL)

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

!#! ! Write Global Attributes
!#! 
!#! !>@todo FIXME  make writing the grid info optional.
!#! !> based on a namelist setting.  if not writing grid,
!#! !> this routine has nothing to do.
!#! 
!#! if (minimal_output) return
!#! 
!#! ! add grid info
!#! call nc_redef(ncid)
!#! 
!#! call nc_add_global_creation_time(ncid)
!#! 
!#! call nc_add_global_attribute(ncid, "model_source", source)
!#! call nc_add_global_attribute(ncid, "model_revision", revision)
!#! call nc_add_global_attribute(ncid, "model_revdate", revdate)
!#! 
!#! call nc_add_global_attribute(ncid, "model", "CAM")
!#! 
!#! ! We need to output the grid information
!#! ! Define the new dimensions IDs
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='xi_rho',  len = Nxi_rho, &
!#!      dimid = nxirhoDimID),'nc_write_model_atts', 'xi_rho def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='eta_rho', len = Neta_rho,&
!#!      dimid = netarhoDimID),'nc_write_model_atts', 'eta_rho def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='s_rho',   len = Ns_rho,&
!#!      dimid = nsrhoDimID),'nc_write_model_atts', 's_rho def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='s_w',   len = Ns_w,&
!#!      dimid = nswDimID),'nc_write_model_atts', 's_w def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='xi_u',    len = Nxi_u,&
!#!      dimid = nxiuDimID),'nc_write_model_atts', 'xi_u def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='xi_v',    len = Nxi_v,&
!#!      dimid = nxivDimID),'nc_write_model_atts', 'xi_v def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='eta_u',   len = Neta_u,&
!#!      dimid = netauDimID),'nc_write_model_atts', 'eta_u def_dim '//trim(filename))
!#! 
!#! call nc_check(nf90_def_dim(ncid, name='eta_v',   len = Neta_v,&
!#!      dimid = netavDimID),'nc_write_model_atts', 'eta_v def_dim '//trim(filename))
!#! 
!#! ! Create the Coordinate Variables and give them Attributes
!#! ! The values will be added in a later block of code.
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lon_rho', xtype=nf90_double, &
!#!               dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lon_rho def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho longitudes'), &
!#!               'nc_write_model_atts', 'lon_rho long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
!#!               'nc_write_model_atts', 'lon_rho units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lat_rho', xtype=nf90_double, &
!#!               dimids=(/ nxirhoDimID, netarhoDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lat_rho def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'rho latitudes'), &
!#!               'nc_write_model_atts', 'lat_rho long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
!#!               'nc_write_model_atts', 'lat_rho units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lon_u', xtype=nf90_double, &
!#!               dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lon_u def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u longitudes'), &
!#!               'nc_write_model_atts', 'lon_u long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
!#!               'nc_write_model_atts', 'lon_u units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lat_u', xtype=nf90_double, &
!#!               dimids=(/ nxiuDimID, netauDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lat_u def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'u latitudes'), &
!#!               'nc_write_model_atts', 'lat_u long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
!#!               'nc_write_model_atts', 'lat_u units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lon_v', xtype=nf90_double, &
!#!               dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lon_v def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v longitudes'), &
!#!               'nc_write_model_atts', 'lon_v long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_east'), &
!#!               'nc_write_model_atts', 'lon_v units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='lat_v', xtype=nf90_double, &
!#!               dimids=(/ nxivDimID, netavDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'lat_v def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'v latitudes'), &
!#!               'nc_write_model_atts', 'lat_v long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'degrees_north'), &
!#!               'nc_write_model_atts', 'lat_v units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='z_rho', xtype=nf90_double, &
!#!               dimids=(/ nxirhoDimID, netarhoDimID, nsrhoDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'z_rho def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
!#!               'nc_write_model_atts', 'z_rho long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
!#!               'nc_write_model_atts', 'z_rho units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='z_u', xtype=nf90_double, &
!#!               dimids=(/ nxiuDimID, netauDimID, nsrhoDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'z_u def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
!#!               'nc_write_model_atts', 'z_u long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
!#!               'nc_write_model_atts', 'z_u units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='z_v', xtype=nf90_double, &
!#!               dimids=(/ nxivDimID, netavDimID, nsrhoDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'z_v def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
!#!               'nc_write_model_atts', 'z_v long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
!#!               'nc_write_model_atts', 'z_v units '//trim(filename))
!#! 
!#! call nc_check(nf90_def_var(ncid,name='z_w', xtype=nf90_double, &
!#!               dimids=(/ nxirhoDimID, netarhoDimID, nswDimID /), varid=VarID),&
!#!               'nc_write_model_atts', 'z_w def_var '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'z at rho'), &
!#!               'nc_write_model_atts', 'z_w long_name '//trim(filename))
!#! call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'), &
!#!               'nc_write_model_atts', 'z_w units '//trim(filename))
!#! 
!#! ! Finished with dimension/variable definitions, must end 'define' mode to fill.
!#! 
!#! call nc_enddef(ncid)
!#! 
!#! ! Fill the coordinate variable values
!#! 
!#! ! the RHO grid
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lon_rho', VarID), &
!#!               'nc_write_model_atts', 'lon_rho inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, TLON ), &
!#!              'nc_write_model_atts', 'lon_rho put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lat_rho', VarID), &
!#!               'nc_write_model_atts', 'lat_rho inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, TLAT ), &
!#!              'nc_write_model_atts', 'lat_rho put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'z_rho', VarID), &
!#!               'nc_write_model_atts', 'z_rho inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, TDEP ), &
!#!              'nc_write_model_atts', 'z_rho put_var '//trim(filename))
!#! 
!#! ! the U grid
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lon_u', VarID), &
!#!               'nc_write_model_atts', 'lon_u inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, ULON ), &
!#!              'nc_write_model_atts', 'lon_u put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lat_u', VarID), &
!#!               'nc_write_model_atts', 'lat_u inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, ULAT ), &
!#!              'nc_write_model_atts', 'lat_u put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'z_u', VarID), &
!#!               'nc_write_model_atts', 'z_u inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, UDEP ), &
!#!              'nc_write_model_atts', 'z_u put_var '//trim(filename))
!#! 
!#! ! the V grid
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lon_v', VarID), &
!#!               'nc_write_model_atts', 'lon_v inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, VLON ), &
!#!              'nc_write_model_atts', 'lon_v put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'lat_v', VarID), &
!#!               'nc_write_model_atts', 'lat_v inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, VLAT ), &
!#!              'nc_write_model_atts', 'lat_v put_var '//trim(filename))
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'z_v', VarID), &
!#!               'nc_write_model_atts', 'z_v inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, VDEP ), &
!#!              'nc_write_model_atts', 'z_v put_var '//trim(filename))
!#! 
!#! ! the W grid
!#! 
!#! call nc_check(NF90_inq_varid(ncid, 'z_w', VarID), &
!#!               'nc_write_model_atts', 'z_w inq_varid '//trim(filename))
!#! call nc_check(nf90_put_var(ncid, VarID, WDEP ), &
!#!              'nc_write_model_atts', 'z_w put_var '//trim(filename))
!#! 
!#! ! Flush the buffer and leave netCDF file open
!#! call nc_sync(ncid)


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
!#! 
!#! integer  :: ncid, VarID
!#! 
!#! real(r8), parameter :: all_land = 0.001_r8
!#! 
!#! if (.not. allocated(ULAT)) allocate(ULAT(Nxi_u, Neta_u))
!#! if (.not. allocated(ULON)) allocate(ULON(Nxi_u, Neta_u))
!#! if (.not. allocated(UDEP)) allocate(UDEP(Nxi_u, Neta_u, Nz))
!#! 
!#! if (.not. allocated(VLAT)) allocate(VLAT(Nxi_v, Neta_v))
!#! if (.not. allocated(VLON)) allocate(VLON(Nxi_v, Neta_v))
!#! if (.not. allocated(VDEP)) allocate(VDEP(Nxi_v, Neta_v, Nz))
!#! 
!#! if (.not. allocated(TLAT)) allocate(TLAT(Nxi_rho, Neta_rho))
!#! if (.not. allocated(TLON)) allocate(TLON(Nxi_rho, Neta_rho))
!#! if (.not. allocated(TDEP)) allocate(TDEP(Nxi_rho, Neta_rho, Nz))
!#! if (.not. allocated(WDEP)) allocate(WDEP(Nxi_rho, Neta_rho, Ns_w))
!#! 
!#! ! Read the vertical information from the (separate) roms_filename
!#! 
!#! call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
!#!       'get_grid', 'open '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'z_u', VarID), &
!#!       'get_grid', 'inq_varid z_u '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, UDEP), &
!#!       'get_grid', 'get_var z_u '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'z_w', VarID), &
!#!       'get_grid', 'inq_varid z_w '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, WDEP), &
!#!       'get_grid', 'get_var z_w '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'z_v', VarID), &
!#!       'get_grid', 'inq_varid z_v '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, VDEP), &
!#!       'get_grid', 'get_var z_v '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'z_rho', VarID), &
!#!       'get_grid', 'inq_varid z_rho '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, TDEP), &
!#!       'get_grid', 'get_var z_rho '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_close(ncid), &
!#!              'get_var','close '//trim(roms_filename))
!#! 
!#! ! Read the rest of the grid information from the traditional grid file
!#! 
!#! call nc_check(nf90_open(trim(roms_filename), nf90_nowrite, ncid), &
!#!       'get_grid', 'open '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lon_rho', VarID), &
!#!    'get_grid', 'inq_varid lon_rho '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, TLON), &
!#!       'get_grid', 'get_var lon_rho '//trim(roms_filename))
!#! 
!#! where (TLON < 0.0_r8) TLON = TLON + 360.0_r8
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lat_rho', VarID), &
!#!       'get_grid', 'inq_varid lat_rho '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, TLAT), &
!#!       'get_grid', 'get_var lat_rho '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lon_u', VarID), &
!#!       'get_grid', 'inq_varid lon_u '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, ULON), &
!#!       'get_grid', 'get_var lon_u '//trim(roms_filename))
!#! 
!#! where (ULON < 0.0_r8) ULON = ULON + 360.0_r8
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lat_u', VarID), &
!#!       'get_grid', 'inq_varid lat_u '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, ULAT), &
!#!       'get_grid', 'get_var lat_u '//trim(roms_filename))
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lon_v', VarID), &
!#!       'get_grid', 'inq_varid lon_v '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, VLON), &
!#!       'get_grid', 'get_var lon_v '//trim(roms_filename))
!#! 
!#! where (VLON < 0.0_r8) VLON = VLON + 360.0_r8
!#! 
!#! call nc_check(nf90_inq_varid(ncid, 'lat_v', VarID), &
!#!       'get_grid', 'inq_varid lat_v '//trim(roms_filename))
!#! call nc_check(nf90_get_var( ncid, VarID, VLAT), &
!#!       'get_grid', 'get_var lat_v '//trim(roms_filename))
!#! 
!#! ! Be aware that all the depths are negative values.
!#! ! The surface of the ocean is 0.0, the deepest is a big negative value.
!#! 
!#! if (do_output() .and. debug > 0) then
!#!     write(string1,*)'    min/max ULON ',minval(ULON), maxval(ULON)
!#!     write(string2,*)    'min/max ULAT ',minval(ULAT), maxval(ULAT)
!#!     write(string3,*)    'min/max UDEP ',minval(UDEP), maxval(UDEP)
!#!     call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)
!#! 
!#!     write(string1,*)'    min/max VLON ',minval(VLON), maxval(VLON)
!#!     write(string2,*)    'min/max VLAT ',minval(VLAT), maxval(VLAT)
!#!     write(string3,*)    'min/max VDEP ',minval(VDEP), maxval(VDEP)
!#!     call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)
!#! 
!#!     write(string1,*)'    min/max TLON ',minval(TLON), maxval(TLON)
!#!     write(string2,*)    'min/max TLAT ',minval(TLAT), maxval(TLAT)
!#!     write(string3,*)    'min/max TDEP ',minval(TDEP), maxval(TDEP)
!#!     call error_handler(E_MSG,'get_grid',string1, text2=string2, text3=string3)
!#! endif
!#! 
end subroutine get_grid


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
MyLoop : do i = 1, MAX_STATE_VARIABLES

   varname   = trim(variable_array(num_state_table_columns*i-4))
   dartstr   = trim(variable_array(num_state_table_columns*i-3))
   minvalstr = trim(variable_array(num_state_table_columns*i-2))
   maxvalstr = trim(variable_array(num_state_table_columns*i-1))
   updatestr = trim(variable_array(num_state_table_columns*i  ))

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml:model "variables" not fully specified'
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

enddo MyLoop

if (nfields == MAX_STATE_VARIABLES) then
   write(string1,'(2A)') 'WARNING: There is a possibility you need to increase ', &
                         'MAX_STATE_VARIABLES in the global variables in model_mod.f90'

   write(string2,'(A,i4,A)') 'WARNING: you have specified at least ', nfields, &
                             ' perhaps more'

   call error_handler(E_MSG,'set_cam_variable_info:',string1, &
                      source,revision,revdate,text2=string2)
endif

!>TODO JPH: do we need another namelist
! JPH cam_template_filename comes from the namelist and should look like
! a generic restart file.
domain_id = add_domain(cam_template_filename, nfields, var_names, kind_list, &
                       clamp_vals, update_list )

call fill_cam_stagger_info(grid_stagger)

if (debug > 2) call state_structure_info(domain_id)

end subroutine set_cam_variable_info

!-----------------------------------------------------------------------
!>
!> Fill table to tell what type of stagger the variable has
!>


subroutine fill_cam_stagger_info(stagger)
type(cam_stagger), intent(inout) :: stagger

integer :: ivar, jdim, qty_index

allocate(stagger%qty_stagger(0:get_num_quantities()))

stagger%qty_stagger = STAGGER_NONE

do ivar = 1, get_num_variables(domain_id)
   do jdim = 1, get_num_dims(domain_id, ivar)

      if (get_dim_name(domain_id, ivar, jdim) == 'slat') then
         qty_index = get_kind_index(domain_id, ivar) ! qty is kind
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
!> Read type(cam_grid) information
!>
!>
!>
!>
!>
!>


subroutine read_grid_info(grid_file, phis_file, grid )
character(len=*), intent(in)  :: grid_file
character(len=*), intent(in)  :: phis_file
type(cam_grid),   intent(out) :: grid

integer :: ncid

! put this in a subroutine that deals with the grid
call nc_check( nf90_open(grid_file, NF90_NOWRITE, ncid), &
               'read_grid_info', 'open '//trim(grid_file))

! Get the grid info
call get_grid_dimensions(ncid)

!#! call get_grid()

call nc_check( nf90_close(ncid), 'read_grid_info', 'close '//trim(grid_file))

end subroutine read_grid_info

!-----------------------------------------------------------------------
!>
!> 
!>   

subroutine get_grid_dimensions(ncid)
integer, intent(in) :: ncid

call fill_cam_1d_array(ncid, 'lon',  grid_data%lon)
call fill_cam_1d_array(ncid, 'lat',  grid_data%lat)
call fill_cam_1d_array(ncid, 'lev',  grid_data%lev)
call fill_cam_1d_array(ncid, 'ilev', grid_data%ilev) ! for staggered vertical grid
call fill_cam_1d_array(ncid, 'slon', grid_data%slon)
call fill_cam_1d_array(ncid, 'slat', grid_data%slat)
call fill_cam_1d_array(ncid, 'gw',   grid_data%gw)   ! gauss weights
call fill_cam_1d_array(ncid, 'hyai', grid_data%hyai)
call fill_cam_1d_array(ncid, 'hybi', grid_data%hybi)
call fill_cam_1d_array(ncid, 'hyam', grid_data%hyam)
call fill_cam_1d_array(ncid, 'hybm', grid_data%hybm)

! P0 is a scalar with no dimensionality
allocate(grid_data%P0%vals(1))
grid_data%P0%nsize = 1

end subroutine get_grid_dimensions


!-----------------------------------------------------------------------
!>
!> 
!>   


subroutine fill_cam_1d_array(ncid, varname, grid_array)
integer,            intent(in)    :: ncid
character(len=*),   intent(in)    :: varname
type(cam_1d_array), intent(inout) :: grid_array

!>@todo need to check that this exists
call nc_get_variable_size(ncid, varname, grid_array%nsize)
allocate(grid_array%vals(grid_array%nsize))

call nc_get_variable(ncid, varname, grid_array%vals)

if (debug > 10) then
   print*, 'variable name', trim(varname), grid_array%vals
endif

end subroutine fill_cam_1d_array

!-----------------------------------------------------------------------
!>
!> 
!>   


subroutine free_cam_1d_array(grid_array)
type(cam_1d_array), intent(inout) :: grid_array

deallocate(grid_array%vals)

grid_array%nsize = -1

end subroutine free_cam_1d_array


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
