! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE).

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, MISSING_R8, MISSING_R4, MISSING_I, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_time, get_date, &
                             print_time, print_date, set_calendar_type,         &
                             operator(*),  operator(+), operator(-),            &
                             operator(>),  operator(<), operator(/),            &
                             operator(/=), operator(<=), time_type

use     location_mod, only : location_type,      get_close_maxdist_init, &
                             get_close_obs_init, get_close_obs, set_location, &
                             set_location_missing

use    utilities_mod, only : register_module, error_handler, nc_check, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_MSG, logfileunit, file_exist


use     obs_kind_mod, only : KIND_SOIL_TEMPERATURE,   &
                             KIND_SOIL_MOISTURE,      &
                             KIND_LIQUID_WATER,       &
                             KIND_ICE,                &
                             KIND_SNOWCOVER_FRAC,     &
                             KIND_SNOW_THICKNESS,     &
                             KIND_LEAF_CARBON,        &
                             KIND_WATER_TABLE_DEPTH,  &
                             KIND_GEOPOTENTIAL_HEIGHT,&
                             paramname_length,        &
                             get_raw_obs_kind_index

use typesizes
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.
public :: cable_state_to_dart_vector, &
          dart_vector_to_model_file, &
          get_cable_restart_filename

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   '$URL$'
character(len=32 ), parameter :: revision = '$Revision$'
character(len=128), parameter :: revdate  = '$Date$'

type(location_type), allocatable :: state_loc(:)

!------------------------------------------------------------------
! namelist
! things which can/should be in the model_nml

integer :: nfields
integer, parameter :: max_state_variables = 40
integer, parameter :: num_state_table_columns = 2
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .true.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: cable_restart_filename = 'cable_restart.nc'
character(len=256) :: cable_history_filename = 'cable_history.nc'
character(len=obstypelength) :: cable_variables(max_state_variables*num_state_table_columns) = ' '


namelist /model_nml/            &
   cable_restart_filename,      &
   cable_history_filename,      &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   cable_variables

!------------------------------------------------------------------
! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: numdims
   integer :: maxlevels
   integer :: xtype
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   integer :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value      ! intended for future use
   logical  :: has_missing_value   ! intended for future use
   character(len=paramname_length) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

!------------------------------------------------------------------
! Module storage ... general purpose variables, metadata, that sort of thing

integer :: Nmland    ! number of land cells
integer :: Nmp_patch ! number of patches in a gridcell
integer :: Nsoil     ! number of soil layers
integer :: Nsnow     ! number of snow layers

real(r8), allocatable, dimension(:) :: snowd      ! Liquid water eqivalent snow depth
real(r8), allocatable, dimension(:) :: zse        ! depth of each soil layer
real(r8), allocatable, dimension(:) :: mland_lats
real(r8), allocatable, dimension(:) :: mland_lons
integer,  allocatable, dimension(:) :: nap  ! number of active patches
integer,  allocatable, dimension(:) :: patchfrac  ! fraction of vegetated grid cell area

! These are the lookup table equivalents of the location metadata
integer,  allocatable, dimension(:) :: lonixy ! longitude index of parent gridcell
integer,  allocatable, dimension(:) :: latjxy ! latitude  index of parent gridcell
real(r8), allocatable, dimension(:) :: levels ! depth
real(r8), allocatable, dimension(:) :: pfrac  ! fraction of vegetated grid cell area

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

integer         :: model_size      ! the state vector length
type(time_type) :: model_time      ! valid time of the model state
type(time_type) :: time_step       ! smallest time to adv model

real(r8), allocatable, dimension(:) :: ens_mean ! may be needed for forward ops

INTERFACE vector_to_prog_var
     MODULE PROCEDURE vector_to_1d_prog_var
     MODULE PROCEDURE vector_to_2d_prog_var
END INTERFACE

INTERFACE DART_get_var
     MODULE PROCEDURE get_var_1d
     MODULE PROCEDURE get_var_2d
END INTERFACE

INTERFACE get_state_time
     MODULE PROCEDURE get_state_time_ncid
     MODULE PROCEDURE get_state_time_fname
END INTERFACE

contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model.

real(r8) :: x_loc

character(len=NF90_MAX_NAME)          :: varname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, dimlen, varsize
integer :: iunit, io, ivar
integer :: i, j, xi, xj, index1, indexN, indx
integer :: ss, dd

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

if ( module_initialized ) return ! only need to do this once.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_output()) call error_handler(E_MSG,'static_init_model','model_nml values are')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

!---------------------------------------------------------------
! Set the time step ... causes cable namelists to be read.
! FIXME Ensure  time_step is multiple of 'model dynamics timestep' if possible

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time = get_state_time(cable_restart_filename)
time_step  = set_time(assimilation_period_seconds, assimilation_period_days)

call print_date(model_time,'model date is ')
call print_time(model_time,'model time is ')
call print_time(time_step,'model timestep is ')

!---------------------------------------------------------------
! Compile the list of cable variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the cable restart file
!
! Compute the offsets into the state vector for the start of each
! variable type. Requires reading shapes from the cable restart file.
! Record the extent of the data type in the state vector.

call verify_state_variables( cable_variables, cable_restart_filename, nfields )

!---------------------------------------------------------------
! Read all the model metadata

call read_metadata( cable_restart_filename, cable_history_filename )

model_size = progvar(nfields)%indexN

allocate(state_loc(model_size))
allocate(ens_mean(model_size))

! We need a fast way to determine the location and variable type given an index
! into the DART state vector. We are going to make a look-up table assuming
! storage is cheap.

allocate(lonixy(model_size), latjxy(model_size), levels(model_size), pfrac(model_size))

call fill_local_metadata( nfields )

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   ! must do one of these:
   !state_loc(i) =  set_location(x_loc)
   !state_loc(i) =  set_location(x_loc,y_loc,v_loc,v_type)
end do


end subroutine static_init_model




subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time


write(string1,*)'DART should not be trying to advance cable'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

! just so suppress compiler warnings. code unreachable
x(:) = MISSING_R8

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
obs_val = MISSING_R8

! The return code for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus = 1

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

! these should be set to the actual location and obs kind
location = set_location_missing()
if (present(var_type)) var_type = 0

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).

deallocate(state_loc, ens_mean)
deallocate(mland_lats, mland_lons, nap, patchfrac, zse, snowd)
deallocate(lonixy, latjxy, levels, pfrac)

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! As it stands, this routine will work for ANY model, with no modification.
!
! The simplest possible netCDF file would contain a 3D field
! containing the state of 'all' the ensemble members. This requires
! three coordinate variables -- one for each of the dimensions
! [model_size, ensemble_member, time]. A little metadata is useful,
! so we can also create some 'global' attributes.
! This is what is implemented here.
!
! Once the simplest case is working, this routine (and nc_write_model_vars)
! can be extended to create a more logical partitioning of the state vector,
! fundamentally creating a netCDF file with variables that are easily
! plotted. The bgrid model_mod is perhaps a good one to view, keeping
! in mind it is complicated by the fact it has two coordinate systems.
! There are stubs in this template, but they are only stubs.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
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

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

integer :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID       ! netCDF pointer to time dimension           (unlimited)
integer :: DimIDmland, DimIDmp_patch, DimIDsoil, DimIDsnow
integer :: myndims
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array
integer :: VarIDlatitude, VarIDlongitude, VarIDnap, VarIDpatchfrac, VarIDzse, VarID

! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=NF90_MAX_NAME) :: str1
character(len=NF90_MAX_NAME) :: varname

integer :: i, ivar

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file,
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     'nc_write_model_atts', 'inquire')
call nc_check(nf90_redef(ncFileID), 'nc_write_model_atts', 'redef')

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension.
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                            'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid= TimeDimID), &
                            'nc_write_model_atts', 'inq_dimid time')

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
                     ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable',  &
        len=model_size, dimid=StateVarDimID), 'nc_write_model_atts', 'def_dim state')

call nc_check(nf90_def_dim(ncid=ncFileID, name='mland',  &
        len=Nmland, dimid=DimIDmland), 'nc_write_model_atts', 'def_dim Nmland')

call nc_check(nf90_def_dim(ncid=ncFileID, name='mp_patch',  &
        len=Nmp_patch, dimid=DimIDmp_patch), 'nc_write_model_atts', 'def_dim Nmp_patch')

call nc_check(nf90_def_dim(ncid=ncFileID, name='soil',  &
        len=Nsoil, dimid=DimIDsoil), 'nc_write_model_atts', 'def_dim Nsoil')

call nc_check(nf90_def_dim(ncid=ncFileID, name='snow',  &
        len=Nsnow, dimid=DimIDsnow), 'nc_write_model_atts', 'def_dim Nsnow')

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1), &
                          'nc_write_model_atts', 'put_att creation_date')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source), &
                          'nc_write_model_atts', 'put_att model_source')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
                          'nc_write_model_atts', 'put_att model_revision')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate), &
                          'nc_write_model_atts', 'put_att model_revdate')
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model','template'), &
                          'nc_write_model_atts', 'put_att model')

call nc_check(nf90_def_var(ncid=ncFileID,name='latitude', xtype=NF90_FLOAT, &
        dimids=DimIDmland,varid=VarIDlatitude),'nc_write_model_atts','def_var latitude')
call nc_check(nf90_put_att(ncFileID, VarIDlatitude,'units','degrees_north'), &
        'nc_write_model_atts', 'put_att latitude long_name')

call nc_check(nf90_def_var(ncid=ncFileID,name='longitude', xtype=NF90_FLOAT, &
        dimids=DimIDmland,varid=VarIDlongitude),'nc_write_model_atts','def_var longitude')
call nc_check(nf90_put_att(ncFileID, VarIDlongitude,'units','degrees_east'), &
        'nc_write_model_atts', 'put_att longitude long_name')

call nc_check(nf90_def_var(ncid=ncFileID,name='nap', xtype=NF90_FLOAT, &
        dimids=DimIDmland,varid=VarIDnap),'nc_write_model_atts','def_var nap')
call nc_check(nf90_put_att(ncFileID, VarIDnap,'long_name','Number of active patches'), &
        'nc_write_model_atts', 'put_att nap long_name')

call nc_check(nf90_def_var(ncid=ncFileID,name='patchfrac', xtype=NF90_FLOAT, &
        dimids=DimIDmp_patch,varid=VarIDpatchfrac),'nc_write_model_atts','def_var patchfrac')
call nc_check(nf90_put_att(ncFileID, VarIDpatchfrac,'long_name', &
        'Fraction of vegetated grid cell area occupied by a vegetation/soil patch'), &
        'nc_write_model_atts', 'put_att patchfrac long_name')

call nc_check(nf90_def_var(ncid=ncFileID,name='zse', xtype=NF90_FLOAT, &
        dimids=DimIDsoil,varid=VarIDzse),'nc_write_model_atts','def_var zse')
call nc_check(nf90_put_att(ncFileID, VarIDzse,'long_name','Depth of each soil layer'), &
        'nc_write_model_atts', 'put_att zse long_name')

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
   call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=StateVarVarID), &
                             'nc_write_model_atts', 'def_var StateVariable')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID,'long_name','State Variable ID'), &
                             'nc_write_model_atts', 'put_att StateVariable long_name')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units',     'indexical'), &
                             'nc_write_model_atts', 'put_att StateVariable units')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'valid_range', (/ 1, model_size /)), &
                             'nc_write_model_atts', 'put_att StateVariable valid_range')

   ! Define the actual (3D) state vector, which gets filled as time goes on ...
   call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=StateVarID), 'nc_write_model_atts', 'def_var state')
   call nc_check(nf90_put_att(ncFileID, StateVarID, 'long_name', 'model state or fcopy'), &
                             'nc_write_model_atts', 'put_att state long_name')

   ! Leave define mode so we can fill variables.
   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'enddef')

   ! Fill the metadata variables
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)), &
                                    'nc_write_model_atts', 'put_var state')
else

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)

      ! match shape of the variable to the dimension IDs

      call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

      write(*,*)'TJH checking ',trim(varname), progvar(ivar)%xtype, mydimids(1:myndims)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
                    xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, &
              'long_name', trim(progvar(ivar)%long_name)), &
              'nc_write_model_atts', trim(string1)//' put_att long_name' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'DART_kind', trim(progvar(ivar)%kind_string)), &
              'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'units', trim(progvar(ivar)%units)), &
              'nc_write_model_atts', trim(string1)//' put_att units' )

      ! Preserve the DART missing_value/_FillValue code.

      if (  progvar(ivar)%xtype == NF90_INT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', MISSING_I), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  MISSING_I), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', MISSING_R4), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  MISSING_R4), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )

      elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 'missing_value', MISSING_R8), &
                 'nc_write_model_atts', trim(string1)//' put_att missing_value' )
         call nc_check(nf90_put_att(ncFileID, VarID, &
                 '_FillValue',  MISSING_R8), &
                 'nc_write_model_atts', trim(string1)//' put_att _FillValue' )
      endif

   enddo

   ! Leave define mode so we can fill variables.
   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts', 'enddef')

endif

call nc_check(nf90_put_var(ncFileID, VarIDlatitude, mland_lats), &
                                    'nc_write_model_atts', 'put_var mland_lats')
call nc_check(nf90_put_var(ncFileID, VarIDlongitude, mland_lons), &
                                    'nc_write_model_atts', 'put_var mland_lons')
call nc_check(nf90_put_var(ncFileID, VarIDnap, nap), &
                                    'nc_write_model_atts', 'put_var nap')
call nc_check(nf90_put_var(ncFileID, VarIDpatchfrac, patchfrac), &
                                    'nc_write_model_atts', 'put_var patchfrac')
call nc_check(nf90_put_var(ncFileID, VarIDzse, zse), &
                                    'nc_write_model_atts', 'put_var zse')

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!------------------------------------------------------------------
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

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, ncNdims, dimlen
integer :: TimeDimID, CopyDimID, VarID
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array

ierr = -1 ! assume things go poorly

!-------------------------------------------------------------------------------
! We will need to know the dimIDs for 'copy' and 'time' ...
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
        'nc_write_model_vars', 'inq_dimid copy')

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
        'nc_write_model_vars', 'inq_dimid time')

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
        'nc_write_model_vars', 'inquire')

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Blast out the DART state vector
   !----------------------------------------------------------------------------

   call nc_check(nf90_inq_varid(ncFileID, 'state', VarID), &
                               'nc_write_model_vars', 'inq_varid state' )
   call nc_check(nf90_put_var(ncFileID, VarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             'nc_write_model_vars', 'put_var state')

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields  ! Very similar to loop in sv_to_restart_file

      varname = trim(progvar(ivar)%varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(varname))

      call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(varname))

      mystart = 1   ! These are arrays, actually
      mycount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(varname)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*)trim(varname),' dim/dimlen ',i,dimlen, &
                            ' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ((debug > 7) .and. do_output()) then
         write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',mystart(1:ncNdims)
         write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',mycount(1:ncNdims)
      endif

      if (     progvar(ivar)%numdims == 1 ) then

         if ( ncNdims /= 3 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
         call vector_to_prog_var(statevec, ivar, data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(varname))
         deallocate(data_1d_array)

      elseif ( progvar(ivar)%numdims == 2 ) then

         if ( ncNdims /= 4 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                                 progvar(ivar)%dimlens(2) ))
         call vector_to_prog_var(statevec, ivar, data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(varname))
         deallocate(data_2d_array)

      else

         write(string1,*)'do not know how to handle CABLE variables with more than 2 dimensions'
         write(string2,*)trim(progvar(ivar)%varname),'has shape', &
                              progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
         call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)

      endif

   enddo
endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================


subroutine cable_state_to_dart_vector(filename, state_vector, state_time)
!------------------------------------------------------------------
! Reads the current time and state variables from the CABLE restart
! file and packs them into a dart state vector.
! Some of the variables are not exploiting the missing_value or _FillValue
! netCDF attributes and so special processing must be employed.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: state_time

! temp space to hold data while we are reading it
integer  :: i, j, ni, nj, ivar, indx
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer :: ncid, dimid, ncNdims, VarID, dimlen
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer :: ind1,ind2

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'cable_state_to_dart_vector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
              'cable_state_to_dart_vector','open '//trim(filename))

state_time = get_state_time(ncid)

if (do_output()) call print_time(state_time,'time in CABLE file '//trim(filename))
if (do_output()) call print_date(state_time,'date in CABLE file '//trim(filename))

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string3 = trim(filename)//' '//trim(varname)
   ind1    = index(progvar(ivar)%long_name,'snow ')
   ind2    = index(progvar(ivar)%long_name,'Snow ')

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'cable_state_to_dart_vector', 'inq_varid '//trim(string3))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'cable_state_to_dart_vector', 'inquire '//trim(string3))

   ! Check the rank of the variable against the local metadata

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'cable_state_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'cable_state_to_dart_vector', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string3),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'cable_state_to_dart_vector',string1,source,revision,revdate)
      endif

   enddo

   ! Pack the variable into the DART state vector
   ! Could/should fill metadata arrays at the same time ...
   ! Snow variables have indeterminate values when the snow depth is zero.
   ! DART requires the indeterminate values to be MISSING_R8
   ! FIXME ... there may be other variables that behave this way ...

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call DART_get_var(ncid, varname, data_1d_array)

      ! pack it into the DART state vector.

      if ((ind1 > 0) .or. (ind2 > 0)) then ! We are a snow variable, take action.
         do i = 1, ni
            if (snowd(i) == MISSING_R8) data_1d_array(i) = MISSING_R8
         enddo
      endif

      do i = 1, ni
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call DART_get_var(ncid, varname, data_2d_array)

      ! pack it into the DART state vector.

      if ((ind1 > 0) .or. (ind2 > 0)) then ! We are a snow variable, take action.
         do i = 1, ni
            if (snowd(i) == MISSING_R8) data_2d_array(i,:) = MISSING_R8
         enddo
      endif

      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'cable_state_to_dart_vector', string1, &
                        source,revision,revdate)
   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'cable_state_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

enddo

call nc_check(nf90_close(ncid),'cable_state_to_dart_vector','close '//trim(filename))
ncid = 0

end subroutine cable_state_to_dart_vector


subroutine dart_vector_to_model_file(state_vector, filename, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a CABLE netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statedate

! temp space to hold data while we are writing it
integer :: i, ni, nj, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer         :: VarID, ncNdims, dimlen
integer         :: ncid
type(time_type) :: file_time

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for writing.'
   call error_handler(E_ERR,'dart_vector_to_model_file',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncid), &
             'dart_vector_to_model_file','open '//trim(filename))

! make sure the time in the file is the same as the time on the data
! we are trying to insert.  we are only updating part of the contents
! of the cable restart file, and state vector contents from a different
! time won't be consistent with the rest of the file.

file_time = get_state_time(ncid)

if ( file_time /= model_time ) then
   call print_time(model_time,'DART  current time',logfileunit)
   call print_time(file_time, 'CABLE current time',logfileunit)
   call print_time(model_time,'DART  current time')
   call print_time(file_time, 'CABLE current time')
   write(string1,*)trim(filename),' current time /= model time. FATAL error.'
   call error_handler(E_ERR,'dart_vector_to_model_file',string1,source,revision,revdate)
endif

if (do_output()) call print_time(file_time,'time of restart file '//trim(filename))
if (do_output()) call print_date(file_time,'date of restart file '//trim(filename))

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'dart_vector_to_model_file', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'dart_vector_to_model_file', 'inquire '//trim(string2))

   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'dart_vector_to_model_file', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'dart_vector_to_model_file', string1, &
                         source, revision, revdate, text2=string2)
      endif

   enddo DimCheck

   ! When called with a 4th argument, vector_to_prog_var() replaces the DART
   ! missing code with the value in the corresponding variable in the original netCDF file.
   ! Any clamping to physically meaningful values occurrs in vector_to_prog_var.

   if (progvar(ivar)%numdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call vector_to_prog_var(state_vector, ivar, data_1d_array, ncid)

      call nc_check(nf90_put_var(ncid, VarID, data_1d_array), &
            'dart_vector_to_model_file', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call vector_to_prog_var(state_vector, ivar, data_2d_array, ncid)

      call nc_check(nf90_put_var(ncid, VarID, data_2d_array), &
            'dart_vector_to_model_file', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'dart_vector_to_model_file', string1, &
                        source,revision,revdate)
   endif

   ! TJH FIXME ... this works perfectly if it were not for a bug in netCDF.
   ! When they fix the bug, this will be a useful thing to restore.
   ! Make note that the variable has been updated by DART
   call nc_check(nf90_Redef(ncid), &
           'dart_vector_to_model_file','redef '//trim(filename))
   call nc_check(nf90_put_att(ncid, VarID,'DART','variable modified by DART'),&
           'dart_vector_to_model_file', 'modified '//trim(varname))
   call nc_check(nf90_enddef(ncid), &
           'dart_vector_to_model_file','state enddef '//trim(filename))

enddo UPDATE

call nc_check(nf90_close(ncid),'dart_vector_to_model_file','close '//trim(filename))
ncid = 0

write(string1,*)'dart_vector_to_model_file needs to be checked'
call error_handler(E_MSG,'dart_vector_to_model_file',string1,source,revision,revdate)

end subroutine dart_vector_to_model_file


!------------------------------------------------------------------


subroutine get_cable_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(cable_restart_filename)

end subroutine get_cable_restart_filename


!------------------------------------------------------------------


function get_state_time_ncid( ncid )
!------------------------------------------------------------------
! The restart netcdf files have the time of the state.

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid

integer :: VarID
integer :: model_seconds(1)
integer :: year, month, day, hour, minute, second
character(len=256) :: unitstring

type(time_type) :: base_time, curr_time


if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, 'time', VarID), 'get_state_time_ncid', &
                      &  'inq_varid time '//trim(cable_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID, model_seconds), 'get_state_time_ncid', &
                      &            'get_var rst_curr_ymd '//trim(cable_restart_filename))

if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
   call nc_check( nf90_get_att(ncid, VarID, 'units' , unitstring), &
           'get_state_time_ncid', 'get_att units '//trim(string2))
else
   write(string1,*) 'cannot get time base from  file <', trim(cable_restart_filename),'>.'
   call error_handler(E_ERR,'get_state_time_ncid',string1,source,revision,revdate)
endif

read(unitstring,'(14x,i4,5(1x,i2))'),year,month,day,hour,minute,second

base_time = set_date(year, month, day, hour, minute, second)
curr_time = set_time(model_seconds(1), 0)

get_state_time_ncid = curr_time + base_time

end function get_state_time_ncid


!------------------------------------------------------------------


function get_state_time_fname(filename)
!------------------------------------------------------------------
! the static_init_model ensures that the cable namelists are read.
!
type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer         :: ncid

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_state_time_fname',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_state_time_fname', 'open '//trim(filename))

get_state_time_fname = get_state_time_ncid(ncid)

call nc_check(nf90_close(ncid),'get_state_time_fname', 'close '//trim(filename))

end function get_state_time_fname


!------------------------------------------------------------------


subroutine verify_state_variables( state_variables, filename, ngood)

character(len=*), dimension(:),   intent(in)  :: state_variables
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood

integer :: index1, indexN, ivar, nrows, ncols
integer :: ncid, VarID, dimlen, varsize
integer :: i
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=obstypelength) :: kind_string
character(len=obstypelength) :: dimname
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs

integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: spvalR8

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, ncid), &
              'verify_state_variables','open '//trim(filename))

! Just make sure the variables exist and that the DART kinds are valid.

nrows = size(variable_table,1)
ncols = size(variable_table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname = trim(state_variables(2*i -1))
   dartstr = trim(state_variables(2*i   ))
   variable_table(i,1) = trim(varname)
   variable_table(i,2) = trim(dartstr)

   if ( variable_table(i,1) == ' ' .and. variable_table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( variable_table(i,1) == ' ' .or. variable_table(i,2) == ' ' ) then
      string1 = 'model_nml:cable_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file

   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 0) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', trim(variable_table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', trim(variable_table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

! Now that we know that the variables exist, read their shapes, etc.
! and fill the progvar vector.

index1  = 1;
indexN  = 0;
do ivar = 1, ngood

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%spvalINT    = MISSING_I
   progvar(ivar)%spvalR4     = MISSING_R4
   progvar(ivar)%spvalR8     = MISSING_R8
   progvar(ivar)%missingINT  = MISSING_I
   progvar(ivar)%missingR4   = MISSING_R4
   progvar(ivar)%missingR8   = MISSING_R8
   progvar(ivar)%has_fill_value    = .false.
   progvar(ivar)%has_missing_value = .false.
   progvar(ivar)%maxlevels   = 0

   string2 = trim(filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'verify_state_variables', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'verify_state_variables', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'verify_state_variables', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'verify_state_variables', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue, missing_value attributes so I can use it when I read and write ...
   ! CESM1_1_1 ... no attributes in the restart file for rank1 or greater
   ! variables.

   if (progvar(ivar)%xtype == NF90_INT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalINT) == NF90_NOERR) then
          progvar(ivar)%spvalINT       = spvalINT
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalINT) == NF90_NOERR) then
          progvar(ivar)%missingINT        = spvalINT
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR4) == NF90_NOERR) then
          progvar(ivar)%spvalR4        = spvalR4
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR4) == NF90_NOERR) then
          progvar(ivar)%missingR4         = spvalR4
          progvar(ivar)%has_missing_value = .true.
       endif

   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
       if (nf90_get_att(ncid, VarID, '_FillValue'    , spvalR8) == NF90_NOERR) then
          progvar(ivar)%spvalR8        = spvalR8
          progvar(ivar)%has_fill_value = .true.
       endif
       if (nf90_get_att(ncid, VarID, 'missing_value' , spvalR8) == NF90_NOERR) then
          progvar(ivar)%missingR8         = spvalR8
          progvar(ivar)%has_missing_value = .true.
       endif
   endif

   ! Now fill the shapes of the variables
   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'verify_state_variables', string1)
      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize = varsize
   progvar(ivar)%index1  = index1
   progvar(ivar)%indexN  = index1 + varsize - 1
   index1                = index1 + varsize      ! sets up for next variable

enddo

call state_report()

call nc_check(nf90_close(ncid),'verify_state_variables','close '//trim(filename))

end subroutine verify_state_variables





subroutine state_report()

integer :: ivar

do ivar = 1,nfields

      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string
      write(logfileunit,*) '  spvalINT    ',progvar(ivar)%spvalINT
      write(logfileunit,*) '  spvalR4     ',progvar(ivar)%spvalR4
      write(logfileunit,*) '  spvalR8     ',progvar(ivar)%spvalR8
      write(logfileunit,*) '  missingINT  ',progvar(ivar)%missingINT
      write(logfileunit,*) '  missingR4   ',progvar(ivar)%missingR4
      write(logfileunit,*) '  missingR8   ',progvar(ivar)%missingR8
      write(logfileunit,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
      write(logfileunit,*) '  has_missing_value ',progvar(ivar)%has_missing_value

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
      write(     *     ,*) '  spvalINT    ',progvar(ivar)%spvalINT
      write(     *     ,*) '  spvalR4     ',progvar(ivar)%spvalR4
      write(     *     ,*) '  spvalR8     ',progvar(ivar)%spvalR8
      write(     *     ,*) '  missingINT  ',progvar(ivar)%missingINT
      write(     *     ,*) '  missingR4   ',progvar(ivar)%missingR4
      write(     *     ,*) '  missingR8   ',progvar(ivar)%missingR8
      write(     *     ,*) '  has_fill_value    ',progvar(ivar)%has_fill_value
      write(     *     ,*) '  has_missing_value ',progvar(ivar)%has_missing_value

enddo

end subroutine state_report


subroutine read_metadata( restart_filename, history_filename )
!
! Some stuff can be read from the restart file, some stuff from the history file
! The restart file does not have the latitude or longitude arrays for the gridcells.
! We need that if we hope to have enough metadata to define a rectangular grid
! and populate it with values. (so we can plot a field of 'soil moisture' for example)

character(len=*), intent(in) :: restart_filename
character(len=*), intent(in) :: history_filename

integer :: i, ncid, dimid, VarID, numdims, dimlen, xtype
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=obstypelength) :: kind_string
character(len=obstypelength) :: dimname
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs

! These are module storage variables set by this routine
! integer :: Nmland    ! number of land cells
! integer :: Nmp_patch ! number of patches in a gridcell
! integer :: Nsoil     ! number of soil layers
! integer :: Nsnow     ! number of snow layers
! real(r8), allocatable, dimension(:) :: mland_lats
! real(r8), allocatable, dimension(:) :: mland_lons
! integer,  allocatable, dimension(:) :: nap  ! number of active patches
! integer,  allocatable, dimension(:) :: patchfrac  ! fraction of vegetated grid cell area
! integer,  allocatable, dimension(:) :: zse  ! depth of each soil layer
! integer,  allocatable, dimension(:) :: snowd

call nc_check(nf90_open(trim(restart_filename), NF90_NOWRITE, ncid), &
                  'read_metadata', 'open '//trim(restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'mland', dimid), &
         'read_metadata','inq_dimid mland '//trim(restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nmland), &
                  'get_history_dims','inquire_dimension mland '//trim(restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'mp_patch', dimid), &
         'read_metadata','inq_dimid mp_patch '//trim(restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nmp_patch), &
                  'get_history_dims','inquire_dimension mp_patch '//trim(restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'soil', dimid), &
         'read_metadata','inq_dimid soil '//trim(restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nsoil), &
                  'get_history_dims','inquire_dimension soil '//trim(restart_filename))

call nc_check(nf90_inq_dimid(ncid, 'snow', dimid), &
         'read_metadata','inq_dimid snow '//trim(restart_filename))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nsnow), &
                  'get_history_dims','inquire_dimension snow '//trim(restart_filename))

allocate(  mland_lats(Nmland))
allocate(  mland_lons(Nmland))
allocate(         nap(Nmland))
allocate(patchfrac(Nmp_patch))
allocate(    snowd(Nmp_patch))
allocate(         zse(Nsoil))

! FIXME ... check shapes of everything ...

! call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=numdims,xtype=xtype),&
!         'read_metadata', 'inquire   '//trim(string3))
! call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), name=dimname, len=Nmland), &
!         'read_metadata', 'inqure dimension '//trim(string3))

string3 = 'latitude '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'latitude', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, mland_lats), &
        'read_metadata', 'get_var'//trim(string3))

string3 = 'longitude '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'longitude', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, mland_lons), &
        'read_metadata', 'get_var'//trim(string3))

string3 = 'nap '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'nap', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, nap), &
        'read_metadata', 'get_var'//trim(string3))

string3 = 'patchfrac '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'patchfrac', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, patchfrac), &
        'read_metadata', 'get_var'//trim(string3))

string3 = 'snowd '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'snowd', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, snowd), &
        'read_metadata', 'get_var'//trim(string3))
where(snowd == 0.0_r8) snowd = MISSING_R8

string3 = 'zse '//trim(restart_filename)
call nc_check(nf90_inq_varid(ncid, 'zse', VarID), &
        'read_metadata', 'inq_varid '//trim(string3))
call nc_check(nf90_get_var(ncid, VarID, zse), &
        'read_metadata', 'get_var'//trim(string3))

end subroutine read_metadata



subroutine fill_local_metadata( ngood )
! Create the metadata arrays that are the same shape as the state vector.
! The metadata arrays will provide the ability to determine what grid cell is the parent
! of the state vector index in question ... as well as the actual surface area.
! This MUST stride through the state vector the same way the state vector is filled.

! integer,  allocatable, dimension(:) :: lonixy ! longitude index of parent gridcell
! integer,  allocatable, dimension(:) :: latjxy ! latitude  index of parent gridcell
! real(r8), allocatable, dimension(:) :: levels ! depth
! real(r8), allocatable, dimension(:) :: pfrac  ! fraction of vegetated grid cell area

integer, intent(in) :: ngood

integer :: ivar, i, j, indx
! Initialize all levels to surface. If there is a level, we will explicitly specify it.
levels(:) = 0.0_r8

do ivar=1, ngood

   ! All variables are at the surface until proven otherwise.
   progvar(ivar)%maxlevels = 1

   indx = progvar(ivar)%index1

   if (progvar(ivar)%numdims == 1) then

      if ((debug > 0) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
      endif

      SELECT CASE ( trim(progvar(ivar)%dimnames(1)) )
         CASE ("mp_patch")
            do i = 1, progvar(ivar)%dimlens(1)
!              xi             = grid1d_ixy(i)
!              xj             = grid1d_jxy(i) ! always unity if unstructured
!              if (unstructured) then
!                 lonixy(  indx) = xi
!                 latjxy(  indx) = xi
!                 landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi)
!              else
                  lonixy(indx) = i        ! FIXME
                  latjxy(indx) = i        ! FIXME
                  pfrac( indx) = patchfrac(i)
!              endif
               indx = indx + 1
            enddo

!        CASE ("plant_carbon_pools")
!           do i = 1, progvar(ivar)%dimlens(1)
!              xi             = cols1d_ixy(i)
!              xj             = cols1d_jxy(i) ! always unity if unstructured
!              if (unstructured) then
!                 lonixy(  indx) = xi
!                 latjxy(  indx) = xi
!                 landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * cols1d_wtxy(i)
!              else
!                 lonixy(  indx) = xi
!                 latjxy(  indx) = xj
!                 landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * cols1d_wtxy(i)
!              endif
!              indx = indx + 1
!           enddo

!        CASE ("soil_carbon_pools")
!           do i = 1, progvar(ivar)%dimlens(1)
!              xi             = pfts1d_ixy(i)
!              xj             = pfts1d_jxy(i) ! always unity if unstructured
!              if (unstructured) then
!                 lonixy(  indx) = xi
!                 latjxy(  indx) = xi
!                 landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi) * pfts1d_wtxy(i)
!              else
!                 lonixy(  indx) = xi
!                 latjxy(  indx) = xj
!                 landarea(indx) = AREA2D(xi,xj) * LANDFRAC2D(xi,xj) * pfts1d_wtxy(i)
!              endif
!              indx = indx + 1
!           enddo

         CASE DEFAULT
            write(string1,*)'(1d) unknown Dimension name '//trim(progvar(ivar)%dimnames(1))// &
             ' while trying to create metadata arrays.'
            call error_handler(E_ERR,'fill_local_metadata', string1, source, revision, revdate)

      END SELECT

   elseif (progvar(ivar)%numdims == 2) then

      if ((debug > 0) .and. do_output()) then
         write(*,*)
         write(*,*)'variable ',trim(progvar(ivar)%varname)
         write(*,*)'dimension 1 (i) ',progvar(ivar)%dimnames(1),progvar(ivar)%dimlens(1)
         write(*,*)'dimension 2 (j) ',progvar(ivar)%dimnames(2),progvar(ivar)%dimlens(2)
      endif

      ! Only dimension 1 matters for the weights, the other is the vertical

      SELECT CASE ( trim(progvar(ivar)%dimnames(1)) )
         CASE ("mp_patch")
!           if ((debug > 0) .and. do_output()) write(*,*)'length grid1d_ixy ',size(grid1d_ixy)
            do j = 1, progvar(ivar)%dimlens(2)
!              xi = grid1d_ixy(j)
!              xj = grid1d_jxy(j) ! always unity if unstructured
               do i = 1, progvar(ivar)%dimlens(1)
!                 if (unstructured) then
!                    lonixy(  indx) = xi
!                    latjxy(  indx) = xi
!                    landarea(indx) = AREA1D(xi) * LANDFRAC1D(xi)
!                 else
                     lonixy(indx) = i
                     latjxy(indx) = j
                     pfrac( indx) = patchfrac(j)
!                 endif
                  indx = indx + 1
               enddo
            enddo

         CASE DEFAULT
            write(string1,*)'(2d) unknown Dimension name '//trim(progvar(ivar)%dimnames(1))// &
             ' while trying to create metadata arrays.'
            call error_handler(E_ERR,'fill_local_metadata', string1, source, revision, revdate)

      END SELECT

   else

      write(string1,*)'variables of rank ',progvar(ivar)%numdims,' are unsupported.'
      write(string2,*)trim(progvar(ivar)%varname),' is dimensioned ',&
                           progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_ERR,'fill_local_metadata', string1, source, revision, &
                         revdate, text2=string2)

   endif

enddo

end subroutine fill_local_metadata



subroutine get_var_1d(ncid, varname, var1d)
! This function will return a R8 array with the netCDF attributes applied.
! scale_factor, offset will be applied,
! missing_value, _FillValue will be replaced by the DART missing value ...

! If _FillValue is defined then it should be scalar and of the same type as the variable.
! If the variable is packed using scale_factor and add_offset attributes (see below),
! the _FillValue attribute should have the data type of the packed data.
!
! missing_value
! When scale_factor and add_offset are used for packing, the value(s) of the missing_value
! attribute should be specified in the domain of the data in the file (the packed data),
! so that missing values can be detected before the scale_factor and add_offset are applied.
!
! scale_factor
! If present for a variable, the data are to be multiplied by this factor after the data
! are read by the application that accesses the data.  If valid values are specified using
! the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
! specified in the domain of the data in the file (the packed data), so that they can be
! interpreted before the scale_factor and add_offset are applied.
!
! add_offset
! If present for a variable, this number is to be added to the data after it is read by
! the application that accesses the data. If both scale_factor and add_offset attributes
! are present, the data are first scaled before the offset is added.

integer,                intent(in)  :: ncid
character(len=*),       intent(in)  :: varname
real(r8), dimension(:), intent(out) :: var1d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs,dimlens
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:) :: intarray
real(r4), allocatable, dimension(:) :: r4array

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_1d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_1d', 'inquire_variable')
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
              'get_var_1d', 'inquire_dimension')

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, intarray), 'get_var_1d', 'get_var')
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, r4array), 'get_var_1d', 'get_var')
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, var1d), 'get_var_1d', 'get_var')

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var1d /= MISSING_R8) var1d = var1d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

end subroutine get_var_1d



subroutine get_var_2d(ncid, varname, var2d)
! This function will return a R8 array with the netCDF attributes applied.
! scale_factor, offset will be applied,
! missing_value, _FillValue will be replaced by the DART missing value ...

! If _FillValue is defined then it should be scalar and of the same type as the variable.
! If the variable is packed using scale_factor and add_offset attributes (see below),
! the _FillValue attribute should have the data type of the packed data.
!
! missing_value
! When scale_factor and add_offset are used for packing, the value(s) of the missing_value
! attribute should be specified in the domain of the data in the file (the packed data),
! so that missing values can be detected before the scale_factor and add_offset are applied.
!
! scale_factor
! If present for a variable, the data are to be multiplied by this factor after the data
! are read by the application that accesses the data.  If valid values are specified using
! the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
! specified in the domain of the data in the file (the packed data), so that they can be
! interpreted before the scale_factor and add_offset are applied.
!
! add_offset
! If present for a variable, this number is to be added to the data after it is read by
! the application that accesses the data. If both scale_factor and add_offset attributes
! are present, the data are first scaled before the offset is added.

integer,                  intent(in)  :: ncid
character(len=*),         intent(in)  :: varname
real(r8), dimension(:,:), intent(out) :: var2d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:) :: intarray
real(r4), allocatable, dimension(:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_2d', 'inquire_variable')

if ( (numdims /= 2)  ) then
   write(string1,*) trim(varname)//' is not a 2D variable as expected.'
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
              'get_var_2d', 'inquire_dimension 1')
call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dimlens(2)), &
              'get_var_2d', 'inquire_dimension 2')

if ( (size(var2d,1) /= dimlens(1))  .or. (size(var2d,2) /= dimlens(2)) ) then
   write(string1,*) trim(varname)//' has shape ', dimlens(1), dimlens(2)
   write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, intarray), 'get_var_2d', 'get_var')
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, r4array), 'get_var_2d', 'get_var')
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, var2d), 'get_var_2d', 'get_var')

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1,source,revision,revdate)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var2d /= MISSING_R8) var2d = var2d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

end subroutine get_var_2d



subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, ncid)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset, into a 1d array.
!
! If the optional argument (ncid) is specified, some additional
! processing takes place. The variable in the netcdf is read.
! This must be the same shape as the intended output array.
! Anywhere the DART MISSING code is encountered in the input array,
! the corresponding (i.e. original) value from the netCDF file is
! used.

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
integer, optional,        intent(in)  :: ncid

integer :: i,ii, VarID
real(r8), allocatable, dimension(:) :: org_array

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do i = 1, progvar(ivar)%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

! Replace the DART fill value with the original value and apply any clamping.
! Get the 'original' variable from the netcdf file.

if (present(ncid)) then

   allocate(org_array(size(data_1d_array)))

   call nc_check(nf90_inq_varid(ncid, progvar(ivar)%varname, VarID), &
            'vector_to_1d_prog_var', 'inq_varid '//trim(progvar(ivar)%varname))

   call nc_check(nf90_get_var(ncid, VarID, org_array), &
            'vector_to_1d_prog_var', 'get_var '//trim(progvar(ivar)%varname))

   ! restoring the indeterminate original values

   where(data_1d_array == MISSING_R8) data_1d_array = org_array

   ! clamping the assimilated values to physically meaningful ranges.

!  if (trim(progvar(ivar)%varname) == 'SNOWDP') &
!     where((data_1d_array < 0.0_r8)) data_1d_array = org_array

!  if (trim(progvar(ivar)%varname) == 'H2OSNO') &
!     where((data_1d_array <= 0.0_r8)) data_1d_array = org_array

   deallocate(org_array)

endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------


subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, ncid)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
integer, optional,        intent(in)  :: ncid

integer :: i,j,ii, VarID
real(r8), allocatable, dimension(:,:) :: org_array

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

! Replace the DART fill value with the original value and apply any clamping.
! Get the 'original' variable from the netcdf file if need be.

if (present(ncid)) then

   allocate(org_array(size(data_2d_array,1),size(data_2d_array,2)))

   call nc_check(nf90_inq_varid(ncid, progvar(ivar)%varname, VarID), &
            'vector_to_2d_prog_var', 'inq_varid '//trim(progvar(ivar)%varname))

   call nc_check(nf90_get_var(ncid, VarID, org_array), &
            'vector_to_2d_prog_var', 'get_var '//trim(progvar(ivar)%varname))

   ! restoring the indeterminate original values

   where(data_2d_array == MISSING_R8 ) data_2d_array = org_array

   ! clamping the assimilated values to physically meaningful ranges.

!  if     (trim(progvar(ivar)%varname) == 'DZSNO') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'ZSNO') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'ZISNO') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'H2OSOI_LIQ') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'H2OSOI_ICE') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'T_SOISNO') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'T_LAKE') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = org_array
!  elseif (trim(progvar(ivar)%varname) == 'leafc') then
!     where((data_2d_array < 0.0_r8)) data_2d_array = 0.0_r8
!  endif

   deallocate(org_array)

endif

end subroutine vector_to_2d_prog_var



subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

integer :: i, mydimid

ndims = 0

do i = 1,progvar(ivar)%numdims

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid

enddo

ndims = ndims + 1
dimids(ndims) = memberdimid
ndims = ndims + 1
dimids(ndims) = unlimitedDimid

if ((debug > 8) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(logfileunit,*)' thus dimids ',dimids(1:ndims)
   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' has dimnames ', &
                       progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
   write(     *     ,*)' thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
