! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

! This is the interface between JULES and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, SECPERDAY, MISSING_R8,                    &
                             MISSING_I, MISSING_R4, rad2deg, deg2rad, PI,      &
                             obstypelength
use time_manager_mod, only : time_type, set_time, get_time, set_date, get_date,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+),  operator(-),          &
                             operator(>),  operator(<),  operator(/),          &
                             operator(/=), operator(<=), operator(==)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, get_close_obs, LocationDims

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             file_exist, find_textfile_dims, file_to_text,     &
                             open_file, close_file

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

use mpi_utilities_mod, only: my_task_id

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf

implicit none
private

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

public :: jules_to_dart_state_vector,   &
          dart_to_jules_restart,        &
          get_jules_restart_filename,   &
          get_state_time,               &
          get_grid_vertval,             &
          compute_gridcell_value,       &
          gridcell_components,          &
          DART_get_var,                 &
          get_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

!------------------------------------------------------------------
!
!  The DART state vector may consist of things like:
!
!  The variables in the jules restart file that are used to create the
!  DART state vector are specified in the input.nml:model_nml namelist.
!
!------------------------------------------------------------------

integer, parameter :: LAKE = 3

! A gridcell may consist of up to 9 tiles. 
! The tiles occupy a fraction of the gridcell.
! If the fractions do not add up to 100%, the
! rest is considered bare soil (tile type 8),
! WHETHER YOU SPECIFY IT OR NOT!
!
! tile types

integer, parameter ::  BROADLEAF_TREE   = 1
integer, parameter ::  NEEDLE_LEAF_TREE = 2
integer, parameter ::  C3_GRASS         = 3
integer, parameter ::  C4_GRASS         = 4 
integer, parameter ::  SHRUBS           = 5
integer, parameter ::  URBAN            = 6
integer, parameter ::  INLAND_WATER     = 7
integer, parameter ::  BARE_SOIL        = 8
integer, parameter ::  LAND_ICE         = 9

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer :: nfields
integer, parameter :: max_state_variables = 20
integer, parameter :: num_state_table_columns = 6
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variables
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_ORIGININDX   = 5 ! ... file of origin
integer, parameter :: VT_STATEINDX    = 6 ! ... update (state) or not

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: jules_restart_filename = 'jules_restart.nc'
character(len=256) :: jules_output_filename = 'jules_history.nc'

character(len=obstypelength) :: variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/            &
   jules_restart_filename,      &
   jules_output_filename,       &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug,                       &
   variables

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer  :: numdims
   integer  :: maxlevels
   integer  :: xtype
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   integer  :: spvalINT, missingINT
   real(r4) :: spvalR4, missingR4
   real(r8) :: spvalR8, missingR8
   logical  :: has_fill_value      ! intended for future use
   logical  :: has_missing_value   ! intended for future use
   character(len=paramname_length) :: kind_string
   character(len=512) :: origin    ! the file it came from
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

!----------------------------------------------------------------------
! how many and which columns are in each gridcell
!----------------------------------------------------------------------

type gridcellcolumns !  given a gridcell, which columns contribute
   private
   integer  :: ncols
   integer, pointer, dimension(:) :: columnids
   integer  :: Ntiles
   integer, pointer, dimension(:) :: tileIds
end type gridcellcolumns
type(gridcellcolumns), allocatable, dimension(:,:), target :: gridCellInfo

!------------------------------------------------------------------------------

integer :: Nlon    = -1   ! output file, AKA 'x'
integer :: Nlat    = -1   ! output file, AKA 'y'
integer :: Nsoil   = -1   ! output, restart, jules_soil.nml
integer :: Ntile   = -1   ! output, restart
integer :: Nland   = -1   ! Number of gridcells containing land
integer :: Nscpool = -1   ! Number of soil carbon pools

real(r8), allocatable :: LONGITUDE(:,:)    ! output file, grid cell centers
real(r8), allocatable ::  LATITUDE(:,:)    ! output file, grid cell centers
real(r8), allocatable :: SOILLEVEL(:)      ! jules_soil.nml, soil interfaces

real(r8), allocatable ::  AREA1D(:),   LANDFRAC1D(:)   ! masked grid
real(r8), allocatable ::  AREA2D(:,:), LANDFRAC2D(:,:) ! 2D grid

logical :: masked = .false.

!------------------------------------------------------------------------------
! These are the 'sparse' type arrays pertaining to the land gridcells.
!
! Unlike the soil layers, snow layer thickness, as well as snow layer depth,
! may change as time goes on (due to snow metamorphism, overburden and the like).
! So there's no uniform levsno as SOILLEVEL coordinate variable.

! The packing order is: Z-LONGITUDE-LATITUDE, Z moving fastest.
! all levels at a location, then
! scan along longitudes, then
! move to next latitude.

integer,  allocatable, dimension(:)  :: grid1d_ixy, grid1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: land1d_ixy, land1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: cols1d_ixy, cols1d_jxy ! 2D lon/lat index of corresponding gridcell
integer,  allocatable, dimension(:)  :: pfts1d_ixy, pfts1d_jxy ! 2D lon/lat index of corresponding gridcell
real(r8), allocatable, dimension(:)  :: land1d_wtxy    ! landunit weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: cols1d_wtxy    ! column   weight relative to corresponding gridcell
real(r8), allocatable, dimension(:)  :: pfts1d_wtxy    ! pft      weight relative to corresponding gridcell
integer,  allocatable, dimension(:)  :: cols1d_ityplun ! columntype ... lake, forest, city ...

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.

real(r8), allocatable, dimension(:) :: ens_mean     ! may be needed for forward ops
integer,  allocatable, dimension(:) :: lonixy       ! longitude index of parent gridcell
integer,  allocatable, dimension(:) :: latjxy       ! latitude  index of parent gridcell
real(r8), allocatable, dimension(:) :: levels       ! depth
real(r8), allocatable, dimension(:) :: landarea     ! land area ... 'support' ... 'weight'

!------------------------------------------------------------------------------
! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer         :: model_size      ! the state vector length
type(time_type) :: model_time      ! valid time of the model state
type(time_type) :: model_timestep  ! smallest time to adv model


INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
END INTERFACE

INTERFACE DART_get_var
      MODULE PROCEDURE get_var_1d
      MODULE PROCEDURE get_var_2d
      MODULE PROCEDURE get_var_3d
      MODULE PROCEDURE get_var_4d
END INTERFACE

INTERFACE get_state_time
      MODULE PROCEDURE get_state_time_ncid
      MODULE PROCEDURE get_state_time_fname
END INTERFACE


contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


function get_model_size()
!------------------------------------------------------------------
! Returns the size of the model as an integer.
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
!> adv_1step is responsible for advancing JULES as a subroutine call.
!> since JULES is not subroutine callable, this is a stub.
!> any 'async' value other than 0 will result in an error.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

call error_handler(E_MSG, 'adv_1step', 'FIXME RAFAEL routine not tested yet (needs long obs_seq file) ', source, revision, revdate)
 
if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*)'DART should not be trying to advance JULES as a subroutine.'
write(string2,*)'DART can only run JULES as a shell script advance.'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate,text2=string2)

! just so suppress compiler warnings. code unreachable
x(:) = MISSING_R8

end subroutine adv_1step


!------------------------------------------------------------------
!> get_state_meta_data is responsible for returning DART metadata.
!> Given an integer index into the state vector structure, returns the
!> associated array indices for lat, lon, and height, as well as the type.

subroutine get_state_meta_data(indx, location, var_type)

integer, intent(in)            :: indx
type(location_type)            :: location
integer, OPTIONAL, intent(out) :: var_type

! Local variables

integer  :: n

! Module variables

! LONGITUDE
! LATITUDE
! lonixy
! latjxy
! levels
! progvar

call error_handler(E_ERR, 'get_state_meta_data', 'FIXME TJH routine not written', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! TJH location = set_location( LONGITUDE(lonixy(indx)), LATITUDE(latjxy(indx)), levels(indx), VERTISHEIGHT)

if (present(var_type)) then

   var_type = MISSING_I

   FINDTYPE : do n = 1,nfields
      if((indx >= progvar(n)%index1) .and. &
         (indx <= progvar(n)%indexN) ) then
         var_type = progvar(n)%dart_kind
         exit FINDTYPE
      endif
   enddo FINDTYPE

   if( var_type == MISSING_I ) then
      write(string1,*) 'Problem, cannot find base_offset, indx is: ', indx
      call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
   endif

endif

return
end subroutine get_state_meta_data



!------------------------------------------------------------------
!> model_interpolate is the basis for all 'forward observation operator's 
!> For a given lat, lon, and height, interpolate the correct state value
!> to that location for the filter from the jules state vectors
!>
!> Reconstructing the vertical profile of the gridcell is complicated.
!> Each land unit/column can have a different number of vertical levels.
!> Do I just try to recontruct the whole profile and then interpolate?
!> Impossible to know which elements are 'above' and 'below' without
!> finding all the elements in the first place. The vertical information
!> is in the levels() array for each state vector component.

subroutine model_interpolate(x, location, obs_kind, interp_val, istatus)

real(r8),            intent(in)  :: x(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_kind
real(r8),            intent(out) :: interp_val
integer,             intent(out) :: istatus

! Local storage

real(r8), dimension(LocationDims) :: loc_array
real(r8) :: llon, llat, lheight
real(r8) :: interp_val_2
integer  :: istatus_2

call error_handler(E_ERR, 'model_interpolate', 'FIXME routine not written', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

interp_val   = MISSING_R8     ! the DART bad value flag
interp_val_2 = MISSING_R8     ! the DART bad value flag
istatus = 99                ! unknown error

! Get the individual locations values

loc_array = get_location(location)
llon      = loc_array(1)
llat      = loc_array(2)
lheight   = loc_array(3)

if ((debug > 6) .and. do_output()) print *, 'requesting interpolation at ', llon, llat, lheight

! FIXME may be better to check the %maxlevels and kick the interpolation to the
! appropriate routine based on that ... or check the dimnames for the
! vertical coordinate  ...

if (obs_kind == KIND_SOIL_TEMPERATURE) then
   call get_grid_vertval(x, location, 'T_SOISNO',  interp_val, istatus )

elseif (obs_kind == KIND_SOIL_MOISTURE) then
   ! TJH FIXME - actually ROLAND FIXME
   ! This is terrible ... the COSMOS operator wants m3/m3 ... jules is kg/m2
   call get_grid_vertval(x, location, 'H2OSOI_LIQ',interp_val  , istatus   )
   call get_grid_vertval(x, location, 'H2OSOI_ICE',interp_val_2, istatus_2 )
   if ((istatus == 0) .and. (istatus_2 == 0)) then
      interp_val = interp_val + interp_val_2
   else
      interp_val = MISSING_R8
      istatus = 6
   endif

elseif (obs_kind == KIND_LIQUID_WATER ) then
   call get_grid_vertval(x, location, 'H2OSOI_LIQ',interp_val, istatus )
elseif (obs_kind == KIND_ICE ) then
   call get_grid_vertval(x, location, 'H2OSOI_ICE',interp_val, istatus )
elseif (obs_kind == KIND_SNOWCOVER_FRAC ) then
   call compute_gridcell_value(x, location, 'frac_sno', interp_val, istatus)
elseif (obs_kind == KIND_LEAF_CARBON ) then
   call compute_gridcell_value(x, location, 'leafc',    interp_val, istatus)
elseif (obs_kind == KIND_WATER_TABLE_DEPTH ) then
   call compute_gridcell_value(x, location, 'ZWT',    interp_val, istatus)
elseif (obs_kind == KIND_SNOW_THICKNESS ) then
   write(string1,*)'model_interpolate for DZSNO not written yet.'
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
   istatus = 5
elseif ((obs_kind == KIND_GEOPOTENTIAL_HEIGHT) .and. vert_is_level(location)) then
   if (nint(lheight) > Nsoil) then
      interp_val = MISSING_R8
      istatus = 1
   else
      interp_val = SOILLEVEL(nint(lheight))
      istatus = 0
   endif
else
   write(string1,*)'model_interpolate not written for (integer) kind ',obs_kind
   call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)
   istatus = 5
endif

if ((debug > 6) .and. do_output()) write(*,*)'interp_val ',interp_val

end subroutine model_interpolate


!------------------------------------------------------------------
!> get_model_time_step
!> Returns the the time step of the model; the smallest increment
!> in time that the model is capable of advancing the state in a given
!> implementation. This interface is required for all applications.

function get_model_time_step()

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step


!------------------------------------------------------------------

!> static_init_model
!> Called to do one time initialization of the model.
!>
!> All the grid information comes from the initialization of
!> the dart_jules_mod module.

subroutine static_init_model()

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: dimname
integer :: ncid, TimeDimID, VarID, dimlen, varsize
integer :: ncidR ! handle to JULES restart file
integer :: ncidO ! handle to JULES output  file 
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

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_output()) call error_handler(E_MSG,'static_init_model','model_nml values are')
if (do_output()) write(logfileunit, nml=model_nml)
if (do_output()) write(     *     , nml=model_nml)

!---------------------------------------------------------------
! Set the time step ... causes jules namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_time     = get_state_time(jules_output_filename)
model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1)

!---------------------------------------------------------------
! The jules output file has the grid metadata.
! The jules restart files are intentionally lean and, in so doing,
! do not have the full lat/lon arrays nor any depth information.
!
call get_jules_output_dimensions( jules_output_filename )
call get_jules_restart_dimensions(jules_restart_filename)

allocate(LONGITUDE(Nlon,Nlat), LATITUDE(Nlon,Nlat),  SOILLEVEL(Nsoil))

! The soil level values are only available from the namelist, apparently.
call Read_Jules_Soil_Namelist()
call get_full_grid(jules_output_filename)

!---------------------------------------------------------------
! The jules grid in a restart file can take on two forms.
! If there is no masked (non-land) gridcells, it is a regular 2D grid.
! If a mask has been applied, the number of longitudes takes on the number
! of useful land gridcells and the number of latitudes is 1.

! TJH allocate(grid1d_ixy(Nland), grid1d_jxy(Nland))
! TJH allocate(cols1d_ixy(Nscpool),   cols1d_jxy(Nscpool))
! TJH allocate(cols1d_wtxy(Nscpool),  cols1d_ityplun(Nscpool))
! TJH allocate(pfts1d_ixy(Ntile),      pfts1d_jxy(Ntile)     , pfts1d_wtxy(Ntile))

! TJH FIXME 
! call get_sparse_geog(ncidO, jules_restart_filename, 'close')

!---------------------------------------------------------------
! Generate list of land units in each gridcell

allocate(gridCellInfo(Nlon,Nlat))
! TJH FIXME
! call SetLocatorArrays()

!---------------------------------------------------------------
! Compile the list of jules variables to use in the creation
! of the DART state vector. This just checks to see that the
! DART KIND is valid and that the variable was specified correctly.
! Whether or not the JULES variables exist is checked in the
! rather long loop that follows.

nfields = parse_variable_table()

! Compute the offsets into the state vector for the start of each
! variable type. Requires reading shapes from the jules restart file.
! Record the extent of the data type in the state vector.

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   ! convey the information in the variable_table to each progvar. 

   call SetVariableAttributes(ivar)

   ! Open the file for each variable and get dimensions, etc.
   ! Since the JULES restart files have no metadata (!?)
   ! we must get all the per-variable metadata from the 
   ! (hopefully identical) variables in the output file.

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'static_init_model','open '//trim(progvar(ivar)%origin))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   string2 = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_inq_varid(ncid, trim(progvar(ivar)%varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.
   ! FIXME ... at this point, we could check the same variable in the output file
   ! to get the metadata

   if( nf90_inquire_attribute(    ncid, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(ncid, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = progvar(ivar)%varname
   endif

   if( nf90_inquire_attribute(    ncid, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(ncid, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! Saving any FillValue, missing_value attributes so I can use it when I read and write ...

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

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimname, len=dimlen), &
                                          'static_init_model', string1)

      ! The 'land' dimension of this variable must match the number 
      ! of elements in the LONGITUDE, LATITUDE arrays.
      ! The 'x' and 'y' dimensions in the output file already match
      ! so there is no need to check those.
      if ( (dimname == 'land') .and. (dimlen /= Nlon*Nlat) )then
         write(string1,*)'dimension mismatch between restart and output'
         write(string2,*)trim(progvar(ivar)%varname),' land dimension is ',dimlen
         write(string3,*)'number of active land cells from output is ',Nlon*Nlat
         call error_handler(E_ERR, 'static_init_model', string1, &
                    source, revision, revdate, text2=string2, text3=string3)
      endif

      ! The 'soil' dimension of this variable must match the number 
      ! of soil layers from the namelist
      if ( (dimname == 'soil') .and. (dimlen /= Nsoil) )then
         write(string1,*)'dimension mismatch between file and namelist'
         write(string2,*)trim(progvar(ivar)%varname),' soil dimension is ',dimlen
         write(string3,*)'number of soil layers from namelist is ',Nsoil
         call error_handler(E_ERR, 'static_init_model', string1, &
                    source, revision, revdate, text2=string2, text3=string3)
      endif

      ! Only reserve space for a single time slice 
      if (dimIDs(i) == TimeDimID) dimlen = 1

      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ((debug > 0) .and. do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  filename    ',trim(progvar(ivar)%origin)
      write(logfileunit,*) '  update      ',progvar(ivar)%update
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims

      do i = 1,progvar(ivar)%numdims
         write(logfileunit,'(''   dimension ('',i1,'') length '',i10,'' name '',A)') &
                    i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimnames(i))
      enddo

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
      write(logfileunit,*)'   rangeRestricted   ',progvar(ivar)%rangeRestricted
      write(logfileunit,*)'   minvalue          ',progvar(ivar)%minvalue
      write(logfileunit,*)'   maxvalue          ',progvar(ivar)%maxvalue

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  filename    ',trim(progvar(ivar)%origin)
      write(     *     ,*) '  update      ',progvar(ivar)%update
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims

      do i = 1,progvar(ivar)%numdims
         write(  *,'(''   dimension ('',i1,'') length '',i10,'' name '',A)') &
                    i,progvar(ivar)%dimlens(i),trim(progvar(ivar)%dimnames(i))
      enddo

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
      write(     *     ,*)'   rangeRestricted   ',progvar(ivar)%rangeRestricted
      write(     *     ,*)'   minvalue          ',progvar(ivar)%minvalue
      write(     *     ,*)'   maxvalue          ',progvar(ivar)%maxvalue
   endif

   call nc_check(nf90_close(ncid),'static_init_model','close '//trim(string2))
   ncid = 0

enddo

model_size = progvar(nfields)%indexN

if ((debug > 99) .and. do_output()) then
  write(logfileunit, *)
  write(logfileunit,'("grid: Nlon, Nlat, Nsoil =",3(1x,i6))') Nlon, Nlat, Nsoil
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)
  write(     *     ,'("grid: Nlon, Nlat, Nsoil =",3(1x,i6))') Nlon, Nlat, Nsoil
  write(     *     , *)'model_size = ', model_size
endif

allocate(ens_mean(model_size))

end subroutine static_init_model


!------------------------------------------------------------------
!> end_model
!> Does any shutdown and clean-up needed for model.

subroutine end_model()

! if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR, 'end_model', 'FIXME routine not written', source, revision, revdate)

if (masked) then
   deallocate(AREA1D, LANDFRAC1D)
else
   deallocate(AREA2D, LANDFRAC2D)
endif

deallocate(LATITUDE, LONGITUDE, SOILLEVEL)
deallocate(grid1d_ixy, grid1d_jxy)
deallocate(land1d_ixy, land1d_jxy, land1d_wtxy)
deallocate(cols1d_ixy, cols1d_jxy, cols1d_wtxy, cols1d_ityplun)
deallocate(pfts1d_ixy, pfts1d_jxy, pfts1d_wtxy)

deallocate(ens_mean)
deallocate(lonixy, latjxy, landarea)

end subroutine end_model


!------------------------------------------------------------------
!> init_time 
!> Companion interface to init_conditions. Returns a time that is somehow
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'No good way to specify an arbitrary initial time for JULES.'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

! Just set to 0 to silence the compiler warnings.
time = set_time(0,0)

end subroutine init_time


!------------------------------------------------------------------
!> init_conditions
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.
!> If this option is not to be used in perfect_model_obs, or if no
!> synthetic data experiments using perfect_model_obs are planned,
!> this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'No good way to specify an arbitrary initial conditions for JULES.'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

! Just set to 0 to silence the compiler warnings.
x = 0.0_r8

end subroutine init_conditions


!------------------------------------------------------------------
!> nc_write_model_atts
!> Writes the model-specific attributes to a netCDF file.
!> This includes coordinate variables and some metadata, but NOT
!> the actual model state.
!
!> Typical sequence for adding new dimensions,variables,attributes:
!> NF90_OPEN             ! open existing netCDF dataset
!>    NF90_redef         ! put into define mode
!>    NF90_def_dim       ! define additional dimensions (if any)
!>    NF90_def_var       ! define variables: from name, type, and dims
!>    NF90_put_att       ! assign attribute values
!> NF90_ENDDEF           ! end definitions: leave define mode
!>    NF90_put_var       ! provide values for variable
!> NF90_CLOSE            ! close: save updated netCDF dataset

function nc_write_model_atts( ncFileID ) result (ierr)

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

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
integer ::   nlonDimID
integer ::   nlatDimID
integer ::  nsoilDimID
integer ::  nlandDimID
integer ::  ntileDimID
integer :: NscpoolDimID

! for the prognostic variables
integer :: ivar, VarID

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
character(len=NF90_MAX_NAME) :: varname
integer, dimension(NF90_MAX_VAR_DIMS) :: mydimids
integer :: i, myndims

character(len=128) :: filename

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

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                          'nc_write_model_atts', 'inq_dimid copy '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                          'nc_write_model_atts', 'inq_dimid time '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
        dimid = StateVarDimID),'nc_write_model_atts', 'def_dim StateVariable '//trim(filename))

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
           'nc_write_model_atts', 'put_att creation '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
           'nc_write_model_atts', 'put_att source '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
           'nc_write_model_atts', 'put_att revision '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
           'nc_write_model_atts', 'put_att revdate '//trim(filename))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'JULES' ), &
           'nc_write_model_atts', 'put_att model '//trim(filename))

!----------------------------------------------------------------------------
! We need to output the prognostic variables.
!----------------------------------------------------------------------------
! Define the new dimensions IDs
!----------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name='x' , len = Nlon, &
          dimid=nlonDimID),   'nc_write_model_atts', 'def_dim lon '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='y' , len = Nlat, &
          dimid=nlatDimID),   'nc_write_model_atts', 'def_dim lat '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='soil', len = Nsoil, &
          dimid=nsoilDimID),  'nc_write_model_atts', 'def_dim soil '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='land', len = Nland, &
          dimid=nlandDimID),  'nc_write_model_atts', 'def_dim land '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='tile', len = Ntile, &
          dimid=ntileDimID),  'nc_write_model_atts', 'def_dim tile '//trim(filename))

call nc_check(nf90_def_dim(ncid=ncFileID, name='scpool', len = Nscpool, &
          dimid=NscpoolDimID),'nc_write_model_atts', 'def_dim scpool '//trim(filename))

!----------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and the Attributes
!----------------------------------------------------------------------------

! Grid Longitudes
call nc_check(nf90_def_var(ncFileID,name='longitude', xtype=nf90_real, &
              dimids=(/ nlonDimID, nlatDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var longitude '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate longitude'), &
              'nc_write_model_atts', 'put_att longitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'X'),  &
              'nc_write_model_atts', 'put_att longitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
              'nc_write_model_atts', 'put_att longitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
              'nc_write_model_atts', 'put_att longitude valid_range '//trim(filename))

! Grid Latitudes
call nc_check(nf90_def_var(ncFileID,name='latitude', xtype=nf90_real, &
              dimids=(/ nlonDimID, nlatDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var latitude '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate latitude'), &
              'nc_write_model_atts', 'put_att latitude long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
              'nc_write_model_atts', 'put_att latitude cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
              'nc_write_model_atts', 'put_att latitude units '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
              'nc_write_model_atts', 'put_att latitude valid_range '//trim(filename))

! subsurface levels
call nc_check(nf90_def_var(ncFileID,name='soil', xtype=nf90_real, &
              dimids=(/ nsoilDimID /), varid=VarID),&
              'nc_write_model_atts', 'def_var soil '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'coordinate soil levels'), &
              'nc_write_model_atts', 'put_att soil long_name '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
              'nc_write_model_atts', 'put_att soil cartesian_axis '//trim(filename))
call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'm'),  &
              'nc_write_model_atts', 'put_att soil units '//trim(filename))

!----------------------------------------------------------------------------
! Create the (empty) Prognostic Variables and the Attributes
!----------------------------------------------------------------------------

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string1 = trim(filename)//' '//trim(varname)

   ! match shape of the variable to the dimension IDs

   call define_var_dims(ivar, ncFileID, MemberDimID, unlimitedDimID, myndims, mydimids)

   ! define the variable and set the attributes

   call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), &
                 xtype=progvar(ivar)%xtype, &
                 dimids = mydimids(1:myndims), varid=VarID),&
                 'nc_write_model_atts', 'def_var '//trim(string1))

   call nc_check(nf90_put_att(ncFileID, VarID, &
           'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', 'put_att long_name '//trim(string1))
   call nc_check(nf90_put_att(ncFileID, VarID, &
           'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', 'put_att DART_kind '//trim(string1))
   call nc_check(nf90_put_att(ncFileID, VarID, &
           'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', 'put_att units '//trim(string1))

   ! Preserve the original missing_value/_FillValue code.

   if (  progvar(ivar)%xtype == NF90_INT ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalINT), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalINT), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))

   elseif (  progvar(ivar)%xtype == NF90_FLOAT ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalR4), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalR4), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))

   elseif (  progvar(ivar)%xtype == NF90_DOUBLE ) then
      call nc_check(nf90_put_att(ncFileID, VarID, &
              'missing_value', progvar(ivar)%spvalR8), &
              'nc_write_model_atts', 'put_att missing_value '//trim(string1))
      call nc_check(nf90_put_att(ncFileID, VarID, &
              '_FillValue',  progvar(ivar)%spvalR8), &
              'nc_write_model_atts', 'put_att _FillValue '//trim(string1))
   endif

enddo

!----------------------------------------------------------------------------
! Finished with dimension/variable definitions, must end 'define' mode to fill.
!----------------------------------------------------------------------------

call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

!----------------------------------------------------------------------------
! Fill the coordinate variables
!----------------------------------------------------------------------------

call nc_check(nf90_inq_varid(ncFileID, 'longitude', VarID), &
             'nc_write_model_atts', 'inq_varid longitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, LONGITUDE ), &
             'nc_write_model_atts', 'put_var longitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'latitude', VarID), &
             'nc_write_model_atts', 'inq_varid latitude '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, LATITUDE ), &
             'nc_write_model_atts', 'put_var latitude '//trim(filename))

call nc_check(nf90_inq_varid(ncFileID, 'soil', VarID), &
             'nc_write_model_atts', 'inq_varid soil '//trim(filename))
call nc_check(nf90_put_var(ncFileID, VarID, SOILLEVEL ), &
             'nc_write_model_atts', 'put_var soil '//trim(filename))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts


!------------------------------------------------------------------
!> nc_write_model_vars
!> Writes the model variables to a netCDF file.
!
!> All errors are fatal, so the
!> return code is always '0 == normal', since the fatal errors stop execution.

function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME)          :: varname
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

character(len=128) :: filename

call error_handler(E_MSG, 'nc_write_model_vars', 'FIXME TJH routine not tested', &
           source, revision, revdate)

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

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

!----------------------------------------------------------------------------
! We need to process the prognostic variables.
!----------------------------------------------------------------------------

do ivar = 1,nfields  ! Very similar to loop in dart_to_jules_restart

   varname = trim(progvar(ivar)%varname)
   string2 = trim(filename)//' '//trim(varname)

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
         'nc_write_model_vars', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
         'nc_write_model_vars', 'inquire '//trim(string2))

   ncstart = 1   ! These are arrays, actually
   nccount = 1
   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'nc_write_model_vars', trim(string1))

      if (progvar(ivar)%dimnames(i) == 'time') cycle DimCheck

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*)trim(string2),' dim/dimlen ',i,dimlen, &
                         ' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                         source, revision, revdate, text2=trim(string2))
      endif

      nccount(i) = dimlen

   enddo DimCheck

   where(dimIDs == CopyDimID) ncstart = copyindex
   where(dimIDs == CopyDimID) nccount = 1
   where(dimIDs == TimeDimID) ncstart = timeindex
   where(dimIDs == TimeDimID) nccount = 1

   if ((debug > 9) .and. do_output()) then
      write(*,*)'nc_write_model_vars '//trim(varname)//' start is ',ncstart(1:ncNdims)
      write(*,*)'nc_write_model_vars '//trim(varname)//' count is ',nccount(1:ncNdims)
   endif

   ! Since 'time' is a singleton dimension, we can use the same logic
   ! as if it the variable had one less dimension.

   if (     (progvar(ivar)%numdims == 1) .or. &
           ((progvar(ivar)%numdims == 2) .and. &
            (progvar(ivar)%dimnames(2) == 'time')) )then

      if ( ncNdims /= 3 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 3 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_1d_array( progvar(ivar)%dimlens(1) ) )
      call vector_to_prog_var(state_vec, ivar, data_1d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_1d_array)

   elseif ( (progvar(ivar)%numdims == 2) .or. &
           ((progvar(ivar)%numdims == 3) .and. &
            (progvar(ivar)%dimnames(3) == 'time')) )then

      if ( ncNdims /= 4 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 4 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_2d_array( progvar(ivar)%dimlens(1),  &
                              progvar(ivar)%dimlens(2) ))
      call vector_to_prog_var(state_vec, ivar, data_2d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_2d_array)

   elseif ( (progvar(ivar)%numdims == 3) .or. &
           ((progvar(ivar)%numdims == 4) .and. &
            (progvar(ivar)%dimnames(4) == 'time')) )then

      if ( ncNdims /= 5 ) then
         write(string1,*)trim(varname),' no room for copy,time dimensions.'
         write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
         call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                         source, revision, revdate, text2=string2)
      endif

      allocate(data_3d_array( progvar(ivar)%dimlens(1),  &
                              progvar(ivar)%dimlens(2),  &
                              progvar(ivar)%dimlens(3) ))
      call vector_to_prog_var(state_vec, ivar, data_3d_array)
      call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
          start = ncstart(1:ncNdims), count=nccount(1:ncNdims)), &
                'nc_write_model_vars', 'put_var '//trim(string2))
      deallocate(data_3d_array)

   else

      write(string1,*)'do not know how to handle jules variables with more than 3 dimensions'
      write(string2,*)trim(progvar(ivar)%varname),'has shape', &
                           progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      call error_handler(E_ERR,'nc_write_model_vars',string1,source,revision,revdate)

   endif

enddo

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars


!------------------------------------------------------------------
!> pert_model_state
!> Perturbs a single model state for generating initial ensembles.
!> This (required interface) is unsupported in JULES and any attempt
!> to use it will cause DART to terminate. Initial ensemble members
!> are generated externally for JULES applications.

subroutine pert_model_state(state, pert_state, interf_provided)

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
logical, save :: random_seq_init = .false.

call error_handler(E_MSG, 'pert_model_state', 'Unsupported for JULES (TJH please check!!!)', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

call error_handler(E_ERR,'pert_model_state', &
                  'JULES cannot be started from a single vector', &
                  source, revision, revdate, &
                  text2='see comments in jules/model_mod.f90::pert_model_state()')

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! This does keep the compiler error messages down, but this
! section of code must never be reached.
do i=1,size(state)
   pert_state(i) = random_gaussian(random_seq, state(i), &
                                   model_perturbation_amplitude)
enddo

end subroutine pert_model_state


!------------------------------------------------------------------
!> ens_mean_for_model
!> If needed by the model interface, this is the current mean
!> for all state vector items across all ensembles.

subroutine ens_mean_for_model(filter_ens_mean)

real(r8), intent(in) :: filter_ens_mean(:)

call error_handler(E_ERR, 'ens_mean_for_model', 'FIXME routine not tested', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!==================================================================
! The remaining PUBLIC interfaces come next
!==================================================================


!------------------------------------------------------------------
!> jules_to_dart_state_vector
!> Reads the current time and state variables from a JULES restart
!> file and packs them into a dart state vector. This better happen
!> in the same fashion as the metadata arrays are built.

subroutine jules_to_dart_state_vector(state_vector, restart_time)

real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: restart_time

! temp space to hold data while we are reading it
integer  :: i, j, k, ni, nj, nk, ivar, indx, numsnowlevels
integer,  allocatable, dimension(:)         :: snlsno
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array
real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME) :: varname
integer :: io, TimeDimID, VarID, ncNdims, dimlen
integer :: ncid
character(len=256) :: myerrorstring

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

call nc_check(nf90_open(trim(jules_output_filename), NF90_NOWRITE, ncid), &
              'jules_to_dart_state_vector', 'open for time '//jules_output_filename)

restart_time = get_state_time(ncid)

if (do_output()) call print_time(restart_time,'time of model state '//jules_output_filename)
if (do_output()) call print_date(restart_time,'date of model state '//jules_output_filename)

call nc_check(nf90_close(ncid),'jules_to_dart_state_vector','close '//jules_output_filename)

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.
! The DART state vector may be composed of variables from either the restart or output file.
! Sometimes it is useful to have variables from the output file to use for forward ops.

do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   myerrorstring = trim(progvar(ivar)%origin)//' '//trim(progvar(ivar)%varname)

   call nc_check(nf90_open(trim(progvar(ivar)%origin), NF90_NOWRITE, ncid), &
              'jules_to_dart_state_vector','open '//trim(myerrorstring))

   ! File is not required to have a time dimension
   io = nf90_inq_dimid(ncid, 'time', TimeDimID)
   if (io /= NF90_NOERR) TimeDimID = MISSING_I

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'jules_to_dart_state_vector', 'inq_varid '//trim(myerrorstring))

   call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ncNdims), &
            'jules_to_dart_state_vector', 'inquire '//trim(myerrorstring))

   ! Check the rank of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(myerrorstring)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'jules_to_dart_state_vector', string1)

      ! Time dimension will be unity in progvar, but not necessarily
      ! in origin file. We only want a single matching time.
      ! static_init_model() only reserves space for a single time.
      
      if ( dimIDs(i) == TimeDimID ) dimlen = 1
          
      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(myerrorstring),' dim/dimlen ',i,dimlen, &
                              ' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'jules_to_dart_state_vector',string1,source,revision,revdate)
      endif

   enddo

   ! Pack the variable into the DART state vector

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call DART_get_var(ncid, varname, data_1d_array)

      do i = 1, ni
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      ! restart file variables are at most 2D, and never have time
      !   float canopy(tile, land) ;
      !   float t_soil(soil, land) ;
      !   float cs(scpool, land) ;
      !   float gs(land) ;

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call DART_get_var(ncid, varname, data_2d_array)

      do j = 1, nj
      do i = 1, ni
         state_vector(indx) = data_2d_array(i, j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   elseif (ncNdims == 3) then

      ! restart file variables never have 3 dimensions
      ! output  file variables may have 3 dimensions
      !   float precip(time, y, x) ;

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'x')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'y')   .and. &
               (trim(progvar(ivar)%dimnames(3)) == 'time') ) then

         ni = progvar(ivar)%dimlens(1)
         nj = progvar(ivar)%dimlens(2)
       ! nk = progvar(ivar)%dimlens(3) not needed ... time is always a singleton

         allocate(data_3d_array(ni, nj, 1))
         call DART_get_var(ncid, varname, data_3d_array)

         do j = 1, nj
         do i = 1, ni
            state_vector(indx) = data_3d_array(i, j, 1)
            indx = indx + 1
         enddo
         enddo
         deallocate(data_3d_array)
      else

         write(string1, *) '3D variable unexpected shape -- only support x, y, time(=1)'
         write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
         write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
         call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                           source, revision, revdate, text2=string2, text3=string3)

      endif

   elseif (ncNdims == 4) then

      ! output  file variables may have 4 dimensions
      !   float soil_wet(time, soil, y, x)
      !   float    esoil(time, tile, y, x)

      if     ( (trim(progvar(ivar)%dimnames(1)) == 'x')   .and. &
               (trim(progvar(ivar)%dimnames(2)) == 'y')   .and. &
               (trim(progvar(ivar)%dimnames(4)) == 'time') ) then

         ni = progvar(ivar)%dimlens(1)
         nj = progvar(ivar)%dimlens(2)
         nk = progvar(ivar)%dimlens(3)
       ! nl = progvar(ivar)%dimlens(3) not needed ... time is always a singleton

         allocate(data_4d_array(ni, nj, nk, 1))
         call DART_get_var(ncid, varname, data_4d_array)

         do k = 1, nk
         do j = 1, nj
         do i = 1, ni
            state_vector(indx) = data_4d_array(i, j, k, 1)
            indx = indx + 1
         enddo
         enddo
         enddo
         deallocate(data_4d_array)
      else

         write(string1, *) '4D variable unexpected shape -- only support x, y, xxx, time(=1)'
         write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
         write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
         call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                           source, revision, revdate, text2=string2, text3=string3)

      endif

   else

      write(string1, *) 'no support for data array of dimension ', ncNdims
      write(string2, *) 'variable [',trim(progvar(ivar)%varname),']'
      write(string3, *) 'file [',trim(progvar(ivar)%origin),']'
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source, revision, revdate, text2=string2, text3=string3)
   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'jules_to_dart_state_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   call nc_check(nf90_close(ncid),'jules_to_dart_state_vector','close '//progvar(ivar)%origin)
   ncid = 0

enddo

end subroutine jules_to_dart_state_vector


!------------------------------------------------------------------
!> dart_to_jules_restart
!> Writes the current time and state variables from a dart state
!> vector (1d array) into a JULES netcdf restart file.

subroutine dart_to_jules_restart( state_vector )

real(r8),         intent(in) :: state_vector(:)

! temp space to hold data while we are writing it
integer :: i, ni, nj, ivar
real(r8), allocatable, dimension(:)         :: data_1d_array
real(r8), allocatable, dimension(:,:)       :: data_2d_array

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
integer         :: VarID, ncNdims, dimlen
integer         :: ncFileID
type(time_type) :: file_time

call error_handler(E_MSG, 'dart_to_jules_restart', 'FIXME RAFAEL routine not tested', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! Check that the output file exists ...

if ( .not. file_exist(jules_restart_filename) ) then
   write(string1,*) 'cannot open file ', trim(jules_restart_filename),' for writing.'
   call error_handler(E_ERR,'dart_to_jules_restart',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(jules_restart_filename), NF90_WRITE, ncFileID), &
             'dart_to_jules_restart','open '//trim(jules_restart_filename))

! FIXME ... if/when the restart file gets metadata indicating the time,
! it would be a good idea to check the model_time against the time in
! the restart file.

! The DART prognostic variables are only defined for a single time.
! We already checked the assumption that variables are xy2d or xyz3d ...
! IF the netCDF variable has a TIME dimension, it must be the last dimension,
! and we need to read the LAST timestep and effectively squeeze out the
! singleton dimension when we stuff it into the DART state vector.

UPDATE : do ivar=1, nfields

   varname = trim(progvar(ivar)%varname)
   string2 = trim(progvar(ivar)%origin)//' '//trim(varname)

   if ( .not. progvar(ivar)%update ) then
      write(string1,*)'intentionally not updating '//trim(string2) 
      write(string3,*)'as per namelist control in model_nml:variables'
      call error_handler(E_MSG, 'dart_to_jules_restart', string1, text2=string3)
      cycle UPDATE
   endif

   ! Ensure netCDF variable is conformable with progvar quantity.
   ! The TIME and Copy dimensions are intentionally not queried
   ! by looping over the dimensions stored in the progvar type.

   call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'dart_to_jules_restart', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(ncFileID,VarID,dimids=dimIDs,ndims=ncNdims), &
            'dart_to_jules_restart', 'inquire '//trim(string2))

   DimCheck : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
            'dart_to_jules_restart', string1)

      if ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         write(string2,*)' but it should be.'
         call error_handler(E_ERR, 'dart_to_jules_restart', string1, &
                         source, revision, revdate, text2=string2)
      endif

   enddo DimCheck

   ! When called with a 4th argument, vector_to_prog_var() 
   ! clamps to physically meaningful values as specified in the model_mod namelist.

   if (progvar(ivar)%numdims == 1) then

      ni = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(ni))
      call vector_to_prog_var(state_vector, ivar, data_1d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array), &
            'dart_to_jules_restart', 'put_var '//trim(varname))
      deallocate(data_1d_array)

   elseif (progvar(ivar)%numdims == 2) then

      ni = progvar(ivar)%dimlens(1)
      nj = progvar(ivar)%dimlens(2)
      allocate(data_2d_array(ni, nj))
      call vector_to_prog_var(state_vector, ivar, data_2d_array, ncFileID)

      call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array), &
            'dart_to_jules_restart', 'put_var '//trim(varname))
      deallocate(data_2d_array)

   else
      write(string1, *) 'no support for data array of dimension ', ncNdims
      call error_handler(E_ERR,'dart_to_jules_restart', string1, &
                        source,revision,revdate)
   endif

   ! TJH FIXME ... this works perfectly if it were not for a bug in netCDF.
   ! When they fix the bug, this will be a useful thing to restore.
   ! Make note that the variable has been updated by DART
!  call nc_check(nf90_Redef(ncFileID),'dart_to_jules_restart', 'redef '//trim(jules_restart_filename))
!  call nc_check(nf90_put_att(ncFileID, VarID,'DART','variable modified by DART'),&
!                'dart_to_jules_restart', 'modified '//trim(varname))
!  call nc_check(nf90_enddef(ncfileID),'dart_to_jules_restart','state enddef '//trim(jules_restart_filename))

enddo UPDATE

call nc_check(nf90_close(ncFileID),'dart_to_jules_restart','close '//trim(jules_restart_filename))
ncFileID = 0

end subroutine dart_to_jules_restart


!------------------------------------------------------------------
!> get_jules_restart_filename
!> provides access to the filename normally in module storage
!> the filename originally comes from the dart namelist.

subroutine get_jules_restart_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(jules_restart_filename)

end subroutine get_jules_restart_filename


!------------------------------------------------------------------
!> get_jules_output_filename
!> provides access to the filename normally in module storage
!> the filename originally comes from the dart namelist.

subroutine get_jules_output_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(jules_output_filename)

end subroutine get_jules_output_filename


!------------------------------------------------------------------
!> (get_state_time): get_state_time_ncid
!> Since the JULES restart files do not have the time as part of the
!> metadata, the times must be read from the companion output file.
!> The last time in the output file is the time used as the valid time
!> of the model state.
!>
!>   float time(time) ;
!>           time:standard_name = "time" ;
!>           time:long_name = "Time of data" ;
!>           time:units = "seconds since 2014-01-01 03:00:00" ;
!>           time:bounds = "time_bounds" ;
!>           time:calendar = "standard" ;

function get_state_time_ncid( ncid )

type(time_type) :: get_state_time_ncid
integer, intent(in) :: ncid

integer :: VarID, numdims, xtype, ios, dimlen
integer :: rst_curr_ymd, rst_curr_tod, leftover
integer :: year, month, day, hour, minute, second

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames
integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=NF90_MAX_NAME) :: attvalue

real(r8), allocatable :: time_array(:)

type(time_type) :: base_time
type(time_type) :: forecast_length

if ( .not. module_initialized ) call static_init_model

call nc_check(nf90_Inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimID),&
                                   'get_state_time_ncid:', 'inquire '//trim(jules_output_filename))

call nc_check(nf90_inq_varid(ncid, 'time', VarID), 'get_state_time_ncid:', &
                      &  'inq_varid time '//trim(jules_output_filename))
call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_state_time_ncid:', 'inquire_variable time '//trim(jules_output_filename))

if (numdims /= 1) then
   write(string1,*)'"time" is supposed to be 1D.'
   write(string2,*)'"time" has ',numdims,' dimensions.'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2)
endif

if (dimIDs(1) /= unlimitedDimID) then
   write(string1,*)'"time" is supposed to be the unlimited dimension.'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate) 
endif

call nc_check(nf90_inquire_dimension(ncid, unlimitedDimID, len=dimlen), &
        'get_state_time_ncid', 'inquire_dimension time '//trim(jules_output_filename))

! Make sure 'time' is the unlimited dimension and just grab the whole thing,
! use the LAST one ...
! Get the time units attribute so we can add the offset to the base, etc.

call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'get_state_time_ncid:', 'time get_att units '//varname)

! Make sure the calendar is 'standard', which implies a gregorian calendar

! time:units = "seconds since 2014-01-01 03:00:00" ;
!               1234567890123

if (attvalue(1:13) /= 'seconds since') then
   write(string1,*)'expecting time units of [seconds since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'get_state_time_ncid:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

allocate(time_array(dimlen))

call nc_check(nf90_get_var(ncid, VarID, time_array), 'get_state_time_ncid', &
                      &            'get_var time '//trim(jules_output_filename))

forecast_length = set_time(int(time_array(dimlen)),0)
base_time = set_date(year, month, day, hour, minute, second)

get_state_time_ncid = base_time + forecast_length

end function get_state_time_ncid


!------------------------------------------------------------------
!> (get_state_time): get_state_time_fname
!> sometimes it is useful to use the netCDF file name

function get_state_time_fname(filename)

type(time_type) :: get_state_time_fname
character(len=*), intent(in) :: filename

integer :: ncid

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
!> get_grid_vertval
!>
!> Calculate the expected vertical value for the gridcell.
!> Each gridcell value is an area-weighted value of an unknown number of
!> column-based quantities.

subroutine get_grid_vertval(x, location, varstring, interp_val, istatus)

real(r8),            intent(in)  :: x(:)         ! state vector
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
character(len=*),    intent(in)  :: varstring    ! T_SOISNO, H2OSOI_LIQ, H2OSOI_ICE
real(r8),            intent(out) :: interp_val   ! area-weighted result
integer,             intent(out) :: istatus      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, indexi, counter1, counter2
integer  :: gridloni,gridlatj
real(r8), dimension(LocationDims) :: loc
real(r8) :: loc_lat, loc_lon, loc_lev
real(r8) :: value_below, value_above, total_area
real(r8) :: depthbelow, depthabove
real(r8) :: topwght, botwght
real(r8), dimension(1) :: loninds,latinds

real(r8), allocatable, dimension(:)   :: above, below
real(r8), allocatable, dimension(:,:) :: myarea

call error_handler(E_ERR, 'get_grid_vertval', 'FIXME routine not written', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! TJH ! Let's assume failure.  Set return val to missing, then the code can
! TJH ! just set istatus to something indicating why it failed, and return.
! TJH ! If the interpolation is good, the interp_val will be set to the
! TJH ! good value, and the last line here sets istatus to 0.
! TJH ! make any error codes set here be in the 10s
! TJH 
! TJH interp_val = MISSING_R8  ! the DART bad value flag
! TJH istatus    = 99          ! unknown error
! TJH 
! TJH loc        = get_location(location)  ! loc is in DEGREES
! TJH loc_lon    = loc(1)
! TJH loc_lat    = loc(2)
! TJH loc_lev    = loc(3)
! TJH 
! TJH if ( loc_lev < 0.0_r8 ) then
! TJH    write(string1,*)'Cannot support above-ground vertical interpolation.'
! TJH    write(string2,*)'requested a value at a depth of ',loc_lev
! TJH    write(string3,*)'jules has negative depths to indicate above-ground values.'
! TJH    call error_handler(E_ERR,'get_grid_vertval', string1, &
! TJH       source, revision, revdate, text2=string2, text3=string3)
! TJH endif
! TJH 
! TJH ! determine the portion of interest of the state vector
! TJH ivar   = findVarIndex(varstring, 'get_grid_vertval')
! TJH index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
! TJH indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here
! TJH 
! TJH ! BOMBPROOFING - check for a vertical dimension for this variable
! TJH if (progvar(ivar)%maxlevels < 2) then
! TJH    write(string1, *)'Variable '//trim(varstring)//' should not use this routine.'
! TJH    write(string2, *)'use compute_gridcell_value() instead.'
! TJH    call error_handler(E_ERR,'get_grid_vertval', string1, &
! TJH                   source, revision, revdate, text2=string2)
! TJH endif
! TJH 
! TJH ! determine the grid cell for the location
! TJH latinds  = minloc(abs(LATITUDE - loc_lat))   ! these return 'arrays' ...
! TJH loninds  = minloc(abs(LONGITUDE - loc_lon))   ! these return 'arrays' ...
! TJH gridlatj = latinds(1)
! TJH gridloni = loninds(1)
! TJH 
! TJH if ((debug > 4) .and. do_output()) then
! TJH    write(*,*)'get_grid_vertval:targetlon, lon, lon index, level is ', &
! TJH               loc_lon,LONGITUDE(gridloni),gridloni,loc_lev
! TJH    write(*,*)'get_grid_vertval:targetlat, lat, lat index, level is ', &
! TJH               loc_lat,LATITUDE(gridlatj),gridlatj,loc_lev
! TJH endif
! TJH 
! TJH ! Determine the level 'above' and 'below' the desired vertical
! TJH ! The above-ground 'depths' are calculated from ZISNO and are negative.
! TJH ! The 'depths' are all positive numbers, increasingly positive is deeper.
! TJH ! The variables currently supported use the subsurface definitions in
! TJH ! the module variable LEVNGRND.
! TJH 
! TJH if (loc_lev  <= SOILLEVEL(1)) then  ! the top level is so close to the surface
! TJH    depthabove = SOILLEVEL(1)        ! just use the top level
! TJH    depthbelow = SOILLEVEL(1)
! TJH elseif (loc_lev >= maxval(SOILLEVEL)) then  ! at depth, however ... do we
! TJH    depthabove    = maxval(SOILLEVEL)        ! fail or just use the deepest
! TJH    depthbelow    = maxval(SOILLEVEL)        ! I am using the deepest.
! TJH else
! TJH 
! TJH    LAYERS : do indexi = 2,size(SOILLEVEL)
! TJH       if (loc_lev < SOILLEVEL(indexi)) then
! TJH          depthabove = SOILLEVEL(indexi-1)
! TJH          depthbelow = SOILLEVEL(indexi  )
! TJH          exit LAYERS
! TJH       endif
! TJH    enddo LAYERS
! TJH 
! TJH endif
! TJH 
! TJH if ((debug > 4) .and. do_output()) then
! TJH    write(*,*)'get_grid_vertval:depthbelow ',depthbelow,'>= loc_lev', &
! TJH                    loc_lev,'>= depthabove',depthabove
! TJH endif
! TJH 
! TJH ! Determine how many elements can contribute to the gridcell value.
! TJH ! There are multiple column-based contributors, each column has a
! TJH ! separate area-based weight. There are multiple levels.
! TJH ! I believe I have to keep track of all of them to sort out how to
! TJH ! calculate the gridcell value at a particular depth.
! TJH 
! TJH counter1 = 0
! TJH counter2 = 0
! TJH GRIDCELL : do indexi = index1, indexN
! TJH    if ( lonixy(indexi) /=  gridloni )  cycle GRIDCELL
! TJH    if ( latjxy(indexi) /=  gridlatj )  cycle GRIDCELL
! TJH    if (      x(indexi) == MISSING_R8)  cycle GRIDCELL
! TJH 
! TJH    if (levels(indexi) == depthabove) counter1 = counter1 + 1
! TJH    if (levels(indexi) == depthbelow) counter2 = counter2 + 1
! TJH 
! TJH enddo GRIDCELL
! TJH 
! TJH if ( (counter1+counter2) == 0 ) then
! TJH    if ((debug > 0) .and. do_output()) then
! TJH       write(string1, *)'statevector variable '//trim(varstring)//' had no viable data'
! TJH       write(string2, *)'at gridcell lon/lat = (',gridloni,',',gridlatj,')'
! TJH       write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
! TJH       call error_handler(E_MSG,'get_grid_vertval', string1, &
! TJH                      text2=string2,text3=string3)
! TJH    endif
! TJH    return
! TJH endif
! TJH 
! TJH allocate(above(counter1),below(counter2),myarea(max(counter1,counter2),2))
! TJH above  = MISSING_R8
! TJH below  = MISSING_R8
! TJH myarea = 0.0_r8
! TJH 
! TJH counter1 = 0
! TJH counter2 = 0
! TJH ELEMENTS : do indexi = index1, indexN
! TJH 
! TJH    if ( lonixy(indexi) /=  gridloni )  cycle ELEMENTS
! TJH    if ( latjxy(indexi) /=  gridlatj )  cycle ELEMENTS
! TJH    if (      x(indexi) == MISSING_R8)  cycle ELEMENTS
! TJH 
! TJH !  write(*,*)'level ',indexi,' is ',levels(indexi),' location depth is ',loc_lev
! TJH 
! TJH    if (levels(indexi)     == depthabove) then
! TJH       counter1            = counter1 + 1
! TJH       above( counter1)    =        x(indexi)
! TJH       myarea(counter1,1)  = landarea(indexi)
! TJH    endif
! TJH    if (levels(indexi)     == depthbelow) then
! TJH       counter2            = counter2 + 1
! TJH       below( counter2)    =        x(indexi)
! TJH       myarea(counter2,2)  = landarea(indexi)
! TJH    endif
! TJH 
! TJH    if ((levels(indexi) /= depthabove) .and. &
! TJH        (levels(indexi) /= depthbelow)) then
! TJH       cycle ELEMENTS
! TJH    endif
! TJH 
! TJH    if ((debug > 4) .and. do_output()) then
! TJH    write(*,*)
! TJH    write(*,*)'gridcell location match at statevector index',indexi
! TJH    write(*,*)'statevector value is (',x(indexi),')'
! TJH    write(*,*)'area is          (',landarea(indexi),')'
! TJH    write(*,*)'LONGITUDE index is     (',lonixy(indexi),')'
! TJH    write(*,*)'LATITUDE index is     (',latjxy(indexi),')'
! TJH    write(*,*)'gridcell LONGITUDE is  (',LONGITUDE(gridloni),')'
! TJH    write(*,*)'gridcell LATITUDE is  (',LATITUDE(gridlatj),')'
! TJH    write(*,*)'depth        is  (',levels(indexi),')'
! TJH    endif
! TJH 
! TJH enddo ELEMENTS
! TJH 
! TJH ! could arise if the above or below was 'missing' ... but the mate was not.
! TJH 
! TJH if ( counter1 /= counter2 ) then
! TJH    write(string1, *)'Variable '//trim(varstring)//' has peculiar interpolation problems.'
! TJH    write(string2, *)'uneven number of values "above" and "below"'
! TJH    write(string3, *)'counter1 == ',counter1,' /= ',counter2,' == counter2'
! TJH    call error_handler(E_MSG,'get_grid_vertval', string1, &
! TJH                   text2=string2,text3=string3)
! TJH    return
! TJH endif
! TJH 
! TJH ! Determine the value for the level above the depth of interest.
! TJH 
! TJH total_area = sum(myarea(1:counter1,1))
! TJH 
! TJH if ( total_area /= 0.0_r8 ) then
! TJH    ! normalize the area-based weights
! TJH    myarea(1:counter1,1) = myarea(1:counter1,1) / total_area
! TJH    value_above = sum(above(1:counter1) * myarea(1:counter1,1))
! TJH else
! TJH    write(string1, *)'Variable '//trim(varstring)//' had no viable data above'
! TJH    write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthabove,')'
! TJH    write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
! TJH    call error_handler(E_ERR,'get_grid_vertval', string1, &
! TJH                   source, revision, revdate, text2=string2,text3=string3)
! TJH endif
! TJH 
! TJH ! Determine the value for the level below the depth of interest.
! TJH 
! TJH total_area = sum(myarea(1:counter2,2))
! TJH 
! TJH if ( total_area /= 0.0_r8 ) then
! TJH    ! normalize the area-based weights
! TJH    myarea(1:counter2,2) = myarea(1:counter2,2) / total_area
! TJH    value_below = sum(below(1:counter2) * myarea(1:counter2,2))
! TJH else
! TJH    write(string1, *)'Variable '//trim(varstring)//' had no viable data below'
! TJH    write(string2, *)'at gridcell lon/lat/lev = (',gridloni,',',gridlatj,',',depthbelow,')'
! TJH    write(string3, *)'obs lon/lat/lev (',loc_lon,',',loc_lat,',',loc_lev,')'
! TJH    call error_handler(E_ERR,'get_grid_vertval', string1, &
! TJH                   source, revision, revdate, text2=string2,text3=string3)
! TJH endif
! TJH 
! TJH if (depthbelow == depthabove) then
! TJH    topwght = 1.0_r8
! TJH    botwght = 0.0_r8
! TJH else
! TJH    topwght = (depthbelow - loc_lev) / (depthbelow - depthabove)
! TJH    botwght = (loc_lev - depthabove) / (depthbelow - depthabove)
! TJH endif
! TJH 
! TJH interp_val = value_above*topwght + value_below*botwght
! TJH istatus    = 0
! TJH 
! TJH deallocate(above, below, myarea)

end subroutine get_grid_vertval


!------------------------------------------------------------------
!> compute_gridcell_value
!>
!> Each gridcell may contain values for several land units, each land unit may contain
!> several columns, each column may contain several pft's. BUT this routine never
!> aggregates across multiple pft's. So, each gridcell value
!> is an area-weighted value of an unknown number of column-based quantities.

subroutine compute_gridcell_value(x, location, varstring, interp_val, istatus)

real(r8),            intent(in)  :: x(:)         ! state vector
type(location_type), intent(in)  :: location     ! location somewhere in a grid cell
character(len=*),    intent(in)  :: varstring    ! frac_sno, leafc
real(r8),            intent(out) :: interp_val   ! area-weighted result
integer,             intent(out) :: istatus      ! error code (0 == good)

! Local storage

integer  :: ivar, index1, indexN, indexi, counter
integer  :: gridloni,gridlatj
real(r8) :: loc_lat, loc_lon
real(r8) :: total, total_area
real(r8), dimension(1) :: loninds,latinds
real(r8), dimension(LocationDims) :: loc

call error_handler(E_ERR, 'compute_gridcell_value', 'FIXME routine not written', source, revision, revdate)
if ( .not. module_initialized ) call static_init_model

! TJH ! Let's assume failure.  Set return val to missing, then the code can
! TJH ! just set istatus to something indicating why it failed, and return.
! TJH ! If the interpolation is good, the interp_val will be set to the
! TJH ! good value, and the last line here sets istatus to 0.
! TJH ! make any error codes set here be in the 10s
! TJH 
! TJH interp_val = MISSING_R8  ! the DART bad value flag
! TJH istatus    = 99          ! unknown error
! TJH 
! TJH loc        = get_location(location)  ! loc is in DEGREES
! TJH loc_lon    = loc(1)
! TJH loc_lat    = loc(2)
! TJH 
! TJH ! determine the portion of interest of the state vector
! TJH ivar   = findVarIndex(varstring, 'compute_gridcell_value')
! TJH index1 = progvar(ivar)%index1 ! in the DART state vector, start looking here
! TJH indexN = progvar(ivar)%indexN ! in the DART state vector, stop  looking here
! TJH 
! TJH ! BOMBPROOFING - check for a vertical dimension for this variable
! TJH if (progvar(ivar)%maxlevels > 1) then
! TJH    write(string1, *)'Variable '//trim(varstring)//' cannot use this routine.'
! TJH    write(string2, *)'use get_grid_vertval() instead.'
! TJH    call error_handler(E_ERR,'compute_gridcell_value', string1, &
! TJH                   source, revision, revdate, text2=string2)
! TJH endif
! TJH 
! TJH ! determine the grid cell for the location
! TJH latinds  = minloc(abs(LATITUDE - loc_lat))   ! these return 'arrays' ...
! TJH loninds  = minloc(abs(LONGITUDE - loc_lon))   ! these return 'arrays' ...
! TJH gridlatj = latinds(1)
! TJH gridloni = loninds(1)
! TJH 
! TJH if ((debug > 5) .and. do_output()) then
! TJH    write(*,*)'compute_gridcell_value:targetlon, lon, lon index is ',&
! TJH                   loc_lon,LONGITUDE(gridloni),gridloni
! TJH    write(*,*)'compute_gridcell_value:targetlat, lat, lat index is ',&
! TJH                   loc_lat,LATITUDE(gridlatj),gridlatj
! TJH endif
! TJH 
! TJH ! If there is no vertical component, the problem is greatly simplified.
! TJH ! Simply area-weight an average of all pieces in the grid cell.
! TJH ! FIXME ... this is the loop that can exploit the knowledge of what 
! TJH ! columnids or tileIds are needed for any particular gridcell.
! TJH ! gridCellInfo%tileIds, gridCellInfo%columnids
! TJH 
! TJH counter    = 0
! TJH total      = 0.0_r8      ! temp storage for state vector
! TJH total_area = 0.0_r8      ! temp storage for area
! TJH ELEMENTS : do indexi = index1, indexN
! TJH 
! TJH    if (   lonixy(indexi) /=  gridloni ) cycle ELEMENTS
! TJH    if (   latjxy(indexi) /=  gridlatj ) cycle ELEMENTS
! TJH    if (        x(indexi) == MISSING_R8) cycle ELEMENTS
! TJH    if ( landarea(indexi) ==   0.0_r8  ) cycle ELEMENTS
! TJH 
! TJH    counter    = counter    + 1
! TJH    total      = total      + x(indexi)*landarea(indexi)
! TJH    total_area = total_area +           landarea(indexi)
! TJH 
! TJH    if ((debug > 5) .and. do_output()) then
! TJH       write(*,*)
! TJH       write(*,*)'gridcell location match',counter,'at statevector index',indexi
! TJH       write(*,*)'statevector value is (',x(indexi),')'
! TJH       write(*,*)'area is              (',landarea(indexi),')'
! TJH       write(*,*)'LONGITUDE index is         (',lonixy(indexi),')'
! TJH       write(*,*)'LATITUDE index is         (',latjxy(indexi),')'
! TJH       write(*,*)'closest LONGITUDE is       (',LONGITUDE(gridloni),')'
! TJH       write(*,*)'closest LATITUDE is       (',LATITUDE(gridlatj),')'
! TJH       write(*,*)'closest lev is       (',levels(indexi),')'
! TJH    endif
! TJH 
! TJH enddo ELEMENTS
! TJH 
! TJH if (total_area /= 0.0_r8) then ! All good.
! TJH    interp_val = total/total_area
! TJH    istatus    = 0
! TJH else
! TJH    if ((debug > 4) .and. do_output()) then
! TJH       write(string1, *)'Variable '//trim(varstring)//' had no viable data'
! TJH       write(string2, *)'at gridcell ilon/jlat = (',gridloni,',',gridlatj,')'
! TJH       write(string3, *)'obs lon/lat = (',loc_lon,',',loc_lat,')'
! TJH       call error_handler(E_MSG,'compute_gridcell_value', string1, &
! TJH                      text2=string2,text3=string3)
! TJH    endif
! TJH endif
! TJH 
! TJH ! Print more information for the really curious
! TJH if ((debug > 5) .and. do_output()) then
! TJH    write(string1,*)'counter, total, total_area', counter, total, total_area
! TJH    write(string2,*)'interp_val, istatus', interp_val, istatus
! TJH    call error_handler(E_MSG,'compute_gridcell_value', string1, text2=string2)
! TJH endif

end subroutine compute_gridcell_value


!------------------------------------------------------------------
!> gridcell_components
!>
!> In order to exercise some of the routines, it is necessary to know
!> which gridcells have multiple land units
!> which gridcells have multiple columns
!> which gridcells have multiple PFTs
!>
!> This routine simply tells me which gridcells are 'interesting'.
!> Each level counts separately. 1 column with 20 levels ... yields a count of 20
!>
!> It is very similar to SetLocatorArrays(), but only does the landunit,column,pft
!> that the variable uses. It is also public - currently only used by
!> model_mod_check.

subroutine gridcell_components(varstring)

character(len=*), intent(in) :: varstring    ! T_SOISNO, H2OSOI_LIQ

! Local storage

integer :: ivar, indexi, i, j
integer, allocatable, dimension(:,:) :: countmat

call error_handler(E_ERR, 'gridcell_components', 'FIXME routine not written', source, revision, revdate)

if ( .not. module_initialized ) call static_init_model

! Empty ... possibly not needed for jules.

end subroutine gridcell_components


!------------------------------------------------------------------
!> (DART_get_var) get_var_1d
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.

subroutine get_var_1d(ncid, varname, var1d)

integer,                intent(in)  :: ncid
character(len=*),       intent(in)  :: varname
real(r8), dimension(:), intent(out) :: var1d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:) :: intarray
real(r4), allocatable, dimension(:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_1d', 'inq_varid '//varname)
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_1d', 'inquire_variable '//varname)
call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1)), &
              'get_var_1d', 'inquire_dimension '//varname)

if ((numdims /= 1) .or. (size(var1d) /= dimlens(1)) ) then
   write(string1,*) trim(varname)//' is not the expected shape/length of ', size(var1d)
   call error_handler(E_ERR,'get_var_1d',string1,source,revision,revdate)
endif

ncstart = 1
nccount = dimlens(1)

if (dimIDs(1) == TimeDimID) then
   call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_1d', 'inquire_dimension time '//trim(varname))
   timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
   ncstart(1) = timeindex
   nccount(1) = 1
endif

if (do_output() .and. (debug > 8)) then
   write(*,*)'get_var_1d: variable ['//trim(varname)//']'
   write(*,*)'get_var_1d: start ',ncstart(1:numdims)
   write(*,*)'get_var_1d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)
   var1d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
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
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)
   var1d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
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

   call nc_check(nf90_get_var(ncid, VarID, values=var1d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_1d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var1d == spvalR8) var1d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_1d',string1)
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


!------------------------------------------------------------------
!> (DART_get_var) get_var_2d
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.

subroutine get_var_2d(ncid, varname, var2d)

integer,                  intent(in)  :: ncid
character(len=*),         intent(in)  :: varname
real(r8), dimension(:,:), intent(out) :: var2d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: VarID, numdims, xtype, io1, io2, i
integer  :: TimeDimID, time_dimlen, timeindex
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:) :: intarray
real(r4), allocatable, dimension(:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

io1 = nf90_inq_dimid(ncid, 'time', TimeDimID)
if (io1 /= NF90_NOERR) TimeDimID = MISSING_I

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_2d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_2d', 'inquire_variable')

if ( (numdims /= 2)  ) then
   write(string1,*) trim(varname)//' is not a 2D variable as expected.'
   call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate)
endif

ncstart(:) = 1
nccount(:) = 1

DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
              'get_var_2d', string1)

   if ( dimIDs(i) == TimeDimID ) then
      call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen), &
            'get_var_2d', 'inquire_dimension time '//trim(varname))
      timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
      ncstart(i) = timeindex
      dimlens(i) = 1

   elseif ( size(var2d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var2d,1),size(var2d,2)
      call error_handler(E_ERR,'get_var_2d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 8)) then
   write(*,*)'get_var_2d: variable ['//trim(varname)//']'
   write(*,*)'get_var_2d: start ',ncstart(1:numdims)
   write(*,*)'get_var_2d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)
   var2d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
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
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)
   var2d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
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

   call nc_check(nf90_get_var(ncid, VarID, values=var2d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_2d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var2d == spvalR8) var2d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_2d',string1)
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


!------------------------------------------------------------------
!> (DART_get_var) get_var_3d
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.

subroutine get_var_3d(ncid, varname, var3d)

integer,                    intent(in)  :: ncid
character(len=*),           intent(in)  :: varname
real(r8), dimension(:,:,:), intent(out) :: var3d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:,:) :: intarray
real(r4), allocatable, dimension(:,:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

! 3D fields must have a time dimension.
! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'get_var_3d', 'inq_dimid time '//varname)
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_3d', 'inquire_dimension time '//varname)

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_3d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_3d', 'inquire_variable '//varname)

if ( (numdims /= 3)  ) then
   write(string1,*) trim(varname)//' is not a 3D variable as expected.'
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

! only expecting [Nlon,Nlat,time]

ncstart(:) = 1
nccount(:) = 1
DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                 'get_var_3d', trim(string1))

   if ( dimIDs(i) == TimeDimID ) then
       timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
       ncstart(i) = timeindex
       dimlens(i) = 1

   elseif ( size(var3d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var3d,i)
      call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 8)) then
   write(*,*)'get_var_3d: variable ['//trim(varname)//']'
   write(*,*)'get_var_3d: start ',ncstart(1:numdims)
   write(*,*)'get_var_3d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)
   var3d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3)))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)
   var3d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var3d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_3d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var3d == spvalR8) var3d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_3d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var3d /= MISSING_R8) var3d = var3d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_3d',string1,source,revision,revdate)
endif

end subroutine get_var_3d


!------------------------------------------------------------------
!> (DART_get_var) get_var_4d
!>
!> This function will return a R8 array with the netCDF attributes applied.
!> scale_factor, offset will be applied,
!> missing_value, _FillValue will be replaced by the DART missing value ...
!>
!> If _FillValue is defined then it should be scalar and of the same type as the variable.
!> If the variable is packed using scale_factor and add_offset attributes (see below),
!> the _FillValue attribute should have the data type of the packed data.
!>
!> missing_value
!> When scale_factor and add_offset are used for packing, the value(s) of the missing_value
!> attribute should be specified in the domain of the data in the file (the packed data),
!> so that missing values can be detected before the scale_factor and add_offset are applied.
!>
!> scale_factor
!> If present for a variable, the data are to be multiplied by this factor after the data
!> are read by the application that accesses the data.  If valid values are specified using
!> the valid_min, valid_max, valid_range, or _FillValue attributes, those values should be
!> specified in the domain of the data in the file (the packed data), so that they can be
!> interpreted before the scale_factor and add_offset are applied.
!>
!> add_offset
!> If present for a variable, this number is to be added to the data after it is read by
!> the application that accesses the data. If both scale_factor and add_offset attributes
!> are present, the data are first scaled before the offset is added.

subroutine get_var_4d(ncid, varname, var4d)

integer,                      intent(in)  :: ncid
character(len=*),             intent(in)  :: varname
real(r8), dimension(:,:,:,:), intent(out) :: var4d

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, dimlens, ncstart, nccount
integer  :: i, TimeDimID, time_dimlen, timeindex
integer  :: VarID, numdims, xtype, io1, io2
integer  :: spvalINT
real(r4) :: spvalR4
real(r8) :: scale_factor, add_offset
real(r8) :: spvalR8

integer,  allocatable, dimension(:,:,:,:) :: intarray
real(r4), allocatable, dimension(:,:,:,:) :: r4array

if ( .not. module_initialized ) call static_init_model

! a little whitespace makes this a lot more readable
if (do_output() .and. (debug > 1)) then
   write(*,*)
   write(logfileunit,*)
endif

! 4D fields must have a time dimension.
! Need to know the Time Dimension ID and length

call nc_check(nf90_inq_dimid(ncid, 'time', TimeDimID), &
         'get_var_4d', 'inq_dimid time '//varname)
call nc_check(nf90_inquire_dimension(ncid, TimeDimID, len=time_dimlen ), &
         'get_var_4d', 'inquire_dimension time '//varname)

call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), 'get_var_4d', 'inq_varid')
call nc_check(nf90_inquire_variable( ncid, VarID, dimids=dimIDs, ndims=numdims, &
              xtype=xtype), 'get_var_4d', 'inquire_variable '//varname)

if ( (numdims /= 4)  ) then
   write(string1,*) trim(varname)//' is not a 4D variable as expected.'
   call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate)
endif

! only expecting [x,y,[soil,tile],time]

ncstart(:) = 1
nccount(:) = 1
DimCheck : do i = 1,numdims

   write(string1,'(a,i2,1x,A)') 'inquire dimension ',i,trim(varname)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                 'get_var_4d', trim(string1))

   if ( dimIDs(i) == TimeDimID ) then
       timeindex  = FindDesiredTimeIndx(ncid, time_dimlen, varname)
       ncstart(i) = timeindex
       dimlens(i) = 1

   elseif ( size(var4d,i) /= dimlens(i) ) then
      write(string1,*) trim(varname)//' has shape ', dimlens(1:numdims)
      write(string2,*) 'which is not the expected shape of ', size(var4d,i)
      call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate,text2=string2)
   endif

   nccount(i) = dimlens(i)

enddo DimCheck

if (do_output() .and. (debug > 8)) then
   write(*,*)'get_var_4d: variable ['//trim(varname)//']'
   write(*,*)'get_var_4d: start ',ncstart(1:numdims)
   write(*,*)'get_var_4d: count ',nccount(1:numdims)
endif

if (xtype == NF90_INT) then

   allocate(intarray(dimlens(1),dimlens(2),dimlens(3),1))
   call nc_check(nf90_get_var(ncid, VarID, values=intarray, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//varname)
   var4d = intarray  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalINT)
   if (  io1 == NF90_NOERR) where (intarray == spvalINT) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalINT)
   if (  io2 == NF90_NOERR) where (intarray == spvalINT) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalINT,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

   deallocate(intarray)

elseif (xtype == NF90_FLOAT) then

   allocate(r4array(dimlens(1),dimlens(2),dimlens(3),1))
   call nc_check(nf90_get_var(ncid, VarID, values=r4array, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//varname)
   var4d = r4array  ! perform type conversion to desired output type

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR4)
   if (  io1 == NF90_NOERR) where (r4array == spvalR4) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR4)
   if (  io2 == NF90_NOERR) where (r4array == spvalR4) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR4,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

   deallocate(r4array)

elseif (xtype == NF90_DOUBLE) then

   call nc_check(nf90_get_var(ncid, VarID, values=var4d, &
           start=ncstart(1:numdims), count=nccount(1:numdims)), &
           'get_var_4d', 'get_var '//varname)

   io1 = nf90_get_att(ncid, VarID, '_FillValue' , spvalR8)
   if (  io1 == NF90_NOERR) where (var4d == spvalR8) var4d = MISSING_R8
   if ( (io1 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing _FillValue ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io2 = nf90_get_att(ncid, VarID, 'missing_value' , spvalR8)
   if (  io2 == NF90_NOERR) where (var4d == spvalR8) var4d = MISSING_R8
   if ( (io2 == NF90_NOERR) .and. do_output() ) then
      write(string1,*)trim(varname)//': replacing missing_value ',spvalR8,' with ',MISSING_R8
      call error_handler(E_MSG,'get_var_4d',string1)
   endif

   io1 = nf90_get_att(ncid, VarID, 'scale_factor', scale_factor)
   io2 = nf90_get_att(ncid, VarID, 'add_offset'  , add_offset)

   if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor + add_offset
   elseif (io1 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d * scale_factor
   elseif (io2 == NF90_NOERR) then
      where (var4d /= MISSING_R8) var4d = var4d + add_offset
   endif

else
   write(string1,*) trim(varname)//' has unsupported (by DART) xtype of', xtype
   call error_handler(E_ERR,'get_var_4d',string1,source,revision,revdate)
endif

end subroutine get_var_4d


!------------------------------------------------------------------
!>  This routine 

function get_model_time()
type(time_type) :: get_model_time

if ( .not. module_initialized ) call static_init_model

get_model_time = model_time

end function get_model_time



!==================================================================
! The remaining (private) interfaces come last
!==================================================================

!------------------------------------------------------------------
!> (vector_to_prog_var) vector_to_1d_prog_var 
!>
!> convert the values from a 1d array, starting at an offset, into a 1d array.
!>
!> If the optional argument (ncid) is specified, some additional
!> processing takes place. 

subroutine vector_to_1d_prog_var(x, ivar, data_1d_array, ncid)

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:),   intent(out) :: data_1d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,ii, VarID
real(r8), allocatable, dimension(:) :: org_array

call error_handler(E_MSG, 'vector_to_1d_prog_var', 'FIXME routine not tested', source, revision, revdate)

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

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the jules restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array > progvar(ivar)%maxvalue)) &
              data_1d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_1d_array /= MISSING_R8) .and. &
             (data_1d_array < progvar(ivar)%minvalue)) &
              data_1d_array = progvar(ivar)%minvalue
   endif

   ! Replace the DART missing value flag with the one JULES uses.
   ! FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_1d_array == MISSING_I) data_1d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_1d_array == MISSING_R4) data_1d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_1d_array == MISSING_R8) data_1d_array = progvar(ivar)%spvalR8
   endif

endif

end subroutine vector_to_1d_prog_var


!------------------------------------------------------------------
!> (vector_to_prog_var) vector_to_2d_prog_var 
!>
!> convert the values from a 1d array, starting at an offset,
!> into a 2d array.

subroutine vector_to_2d_prog_var(x, ivar, data_2d_array, ncid)

real(r8), dimension(:),   intent(in)  :: x
integer,                  intent(in)  :: ivar
real(r8), dimension(:,:), intent(out) :: data_2d_array
integer, OPTIONAL,        intent(in)  :: ncid

integer :: i,j,ii, VarID
real(r8), allocatable, dimension(:,:) :: org_array

call error_handler(E_MSG, 'vector_to_2d_prog_var', 'FIXME routine not tested', source, revision, revdate)

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

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the jules restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array > progvar(ivar)%maxvalue)) &
              data_2d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_2d_array /= MISSING_R8) .and. &
             (data_2d_array < progvar(ivar)%minvalue)) &
              data_2d_array = progvar(ivar)%minvalue
   endif

   ! replace the missing values with the original missing values.
   ! FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_2d_array == MISSING_I) data_2d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_2d_array == MISSING_R4) data_2d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_2d_array == MISSING_R8) data_2d_array = progvar(ivar)%spvalR8
   endif

endif

end subroutine vector_to_2d_prog_var


!------------------------------------------------------------------
!> (vector_to_prog_var) vector_to_3d_prog_var 
!>
!> convert the values from a 1d array, starting at an offset,
!> into a 3d array.

subroutine vector_to_3d_prog_var(x, ivar, data_3d_array, ncid)

real(r8), dimension(:),     intent(in)  :: x
integer,                    intent(in)  :: ivar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array
integer, OPTIONAL,          intent(in)  :: ncid

integer :: i,j,k,ii, VarID
real(r8), allocatable, dimension(:,:,:) :: org_array

call error_handler(E_MSG, 'vector_to_3d_prog_var', 'FIXME routine not tested', source, revision, revdate)

! unpack the right part of the DART state vector into a 1D array.

ii = progvar(ivar)%index1

do k = 1,progvar(ivar)%dimlens(3)
do j = 1,progvar(ivar)%dimlens(2)
do i = 1,progvar(ivar)%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

ii = ii - 1
if ( ii /= progvar(ivar)%indexN ) then
   write(string1, *)'Variable '//trim(progvar(ivar)%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

! Apply the min/max values, if applicable
! This should only be true when converting to a variable that will
! be reinserted into the jules restart file. This is indicated
! by the presence of the ncid variable.

if (present(ncid)) then

   if ((progvar(ivar)%rangeRestricted == BOUNDED_ABOVE ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array > progvar(ivar)%maxvalue)) &
              data_3d_array = progvar(ivar)%maxvalue
   endif

   if ((progvar(ivar)%rangeRestricted == BOUNDED_BELOW ) .or. &
       (progvar(ivar)%rangeRestricted == BOUNDED_BOTH )) then
      where ((data_3d_array /= MISSING_R8) .and. &
             (data_3d_array < progvar(ivar)%minvalue)) &
              data_3d_array = progvar(ivar)%minvalue
   endif

   ! replace the missing values with the original missing values.
   ! FIXME ... I am not sure if there would ever be any missing values
   ! in a JULES restart file.

   if     (progvar(ivar)%xtype == NF90_INT) then
      where(data_3d_array == MISSING_I) data_3d_array = progvar(ivar)%spvalINT
   elseif (progvar(ivar)%xtype == NF90_FLOAT) then
      where(data_3d_array == MISSING_R4) data_3d_array = progvar(ivar)%spvalR4
   elseif (progvar(ivar)%xtype == NF90_DOUBLE) then
      where(data_3d_array == MISSING_R8) data_3d_array = progvar(ivar)%spvalR8
   endif

endif

end subroutine vector_to_3d_prog_var


!------------------------------------------------------------------
!> get_jules_output_dimensions
!>
!> Read the dimensions from the history netcdf file.
!> The file name comes from module storage ... namelist.

subroutine get_jules_output_dimensions( fname )

character(len=*), intent(in) :: fname

integer :: dimid
integer :: ncid

! sets module variables Nlon, Nlat, Nsoil, Ntile

call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
            'get_jules_output_dimensions','open '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'x', dimid), &
            'get_jules_output_dimensions','inq_dimid x '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlon), &
            'get_jules_output_dimensions','inquire_dimension x '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'y', dimid), &
            'get_jules_output_dimensions','inq_dimid y '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nlat), &
            'get_jules_output_dimensions','inquire_dimension y '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'soil', dimid), &
            'get_jules_output_dimensions','inq_dimid soil '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nsoil), &
            'get_jules_output_dimensions','inquire_dimension soil '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'tile', dimid), &
            'get_jules_output_dimensions','inq_dimid tile '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Ntile), &
            'get_jules_output_dimensions','inquire_dimension tile '//trim(fname))

if ((debug > 9) .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'get_jules_output_dimensions output follows:'
   write(logfileunit,*)'Nlon  = ',Nlon
   write(logfileunit,*)'Nlat  = ',Nlat
   write(logfileunit,*)'Nsoil = ',Nsoil
   write(logfileunit,*)'Ntile = ',Ntile
   write(     *     ,*)
   write(     *     ,*)'get_jules_output_dimensions output follows:'
   write(     *     ,*)'Nlon  = ',Nlon
   write(     *     ,*)'Nlat  = ',Nlat
   write(     *     ,*)'Nsoil = ',Nsoil
   write(     *     ,*)'Ntile = ',Ntile
endif

call nc_check(nf90_close(ncid),'get_jules_output_dimensions','close '//trim(fname) )

end subroutine get_jules_output_dimensions


!------------------------------------------------------------------
!> get_jules_restart_dimensions
!>
!> Read the dimensions from the history netcdf file.
!> The file name comes from module storage ... namelist.

subroutine get_jules_restart_dimensions( fname )

character(len=*), intent(in) :: fname

integer :: dimid
integer :: ncid
integer :: myNsoil, myNtile

! sets module variables Nland, Nscpool, Nsoil, Ntile

call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
            'get_jules_restart_dimensions','open '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'land', dimid), &
            'get_jules_restart_dimensions','inq_dimid land '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nland), &
            'get_jules_restart_dimensions','inquire_dimension land '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'scpool', dimid), &
            'get_jules_restart_dimensions','inq_dimid scpool '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=Nscpool), &
            'get_jules_restart_dimensions','inquire_dimension scpool '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'soil', dimid), &
            'get_jules_restart_dimensions','inq_dimid soil '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=myNsoil), &
            'get_jules_restart_dimensions','inquire_dimension soil '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'tile', dimid), &
            'get_jules_restart_dimensions','inq_dimid tile '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=myNtile), &
            'get_jules_restart_dimensions','inquire_dimension tile '//trim(fname))

if (myNsoil /= Nsoil) then
   write(string1,*)'Confusion about number of soil levels.'
   write(string2,*)'number of soil levels in output  file is ',Nsoil
   write(string3,*)'number of soil levels in restart file is ',myNsoil
   call error_handler(E_ERR,'get_jules_restart_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if (myNtile /= Ntile) then
   write(string1,*)'Confusion about number of tiles.'
   write(string2,*)'number of tiles in output  file is ',Ntile
   write(string3,*)'number of tiles in restart file is ',myNtile
   call error_handler(E_ERR,'get_jules_restart_dimensions',string1, &
              source, revision, revdate, text2=string2, text3=string3)
endif

if ((debug > 9) .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'get_jules_restart_dimensions output follows:'
   write(logfileunit,*)'Nland   = ',Nland
   write(logfileunit,*)'Nscpool = ',Nscpool
   write(logfileunit,*)'Nsoil   = ',myNsoil
   write(logfileunit,*)'Ntile   = ',myNtile
   write(     *     ,*)
   write(     *     ,*)'get_jules_restart_dimensions output follows:'
   write(     *     ,*)'Nland   = ',Nland
   write(     *     ,*)'Nscpool = ',Nscpool
   write(     *     ,*)'Nsoil   = ',myNsoil
   write(     *     ,*)'Ntile   = ',Ntile
endif

call nc_check(nf90_close(ncid),'get_jules_restart_dimensions','close '//trim(fname) )

end subroutine get_jules_restart_dimensions


!------------------------------------------------------------------
!> get_full_grid
!>
!> Read the grid dimensions from the jules output netcdf file.
!> LONGITUDE, LATITUDE, AREA, LANDFRAC, ... all have module scope

subroutine get_full_grid(fname)

character(len=*), intent(in)    :: fname

integer :: ncid

! Make sure the variables are the right size ...
! at some point in the future ...

call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
        'get_full_grid','open '//trim(fname))

! The lat/lon matrices in the history file have been masked by
! the land values such that the wet cells are 'missing' values.
! This makes it less than useful for us.

call DART_get_var(ncid, 'longitude', LONGITUDE)
call DART_get_var(ncid, 'latitude',  LATITUDE)

! just to make sure we are [0,360] and [-90,90]

where (LONGITUDE <   0.0_r8) LONGITUDE = LONGITUDE + 360.0_r8
where (LONGITUDE > 360.0_r8) LONGITUDE = LONGITUDE - 360.0_r8

if (any(LONGITUDE < 0.0_r8)) then
   write(string1,*)'longitudes in history file variable "lon" still negative.'
   call error_handler(E_ERR,'get_full_grid',string1,source,revision,revdate)
endif

where (LATITUDE < -90.0_r8) LATITUDE = -90.0_r8
where (LATITUDE >  90.0_r8) LATITUDE =  90.0_r8

call nc_check(nf90_close(ncid),'get_full_grid','close '//trim(fname) )

! A little sanity check

if ((debug > 7) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'history_file grid information as interpreted ...'
   write(logfileunit,*)'longitude  range ', minval(LONGITUDE), maxval(LONGITUDE)
   write(logfileunit,*)'latitude   range ', minval( LATITUDE), maxval( LATITUDE)
   write(logfileunit,*)'soillevel  is ', SOILLEVEL
   write(     *     ,*)
   write(     *     ,*)'history_file grid information as interpreted ...'
   write(     *     ,*)'longitude  range ', minval(LONGITUDE), maxval(LONGITUDE)
   write(     *     ,*)'longitude  range ', minval( LATITUDE), maxval(LATITUDE)
   write(     *     ,*)'SOILLEVEL  is ', SOILLEVEL

endif

return
end subroutine get_full_grid


!------------------------------------------------------------------
!> get_sparse_geog
!>
!> Read the geography information from from the restart netcdf file.

subroutine get_sparse_geog(ncid, fname, cstat)

integer,          intent(inout) :: ncid
character(len=*), intent(in)    :: fname
character(len=*), intent(in)    :: cstat

integer :: VarID

call error_handler(E_ERR, 'get_sparse_geog', 'FIXME TJH routine not tested', source, revision, revdate)

if (ncid == 0) then
   call nc_check(nf90_open(trim(fname), nf90_nowrite, ncid), &
               'get_sparse_geog','open '//trim(fname))
endif

! Make sure the variables are the right size ...
! by comparing agains the size of the variable ...

if ( Nland < 0 ) then
   write(string1,*)'Unable to read the number of land units.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( Nscpool < 0 ) then
   write(string1,*)'Unable to read the number of columns.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

if ( Ntile < 0 ) then
   write(string1,*)'Unable to read the number of pfts.'
   call error_handler(E_ERR,'get_sparse_geog',string1,source,revision,revdate)
endif

! Read the netcdf file data

call nc_check(nf90_inq_varid(ncid, 'grid1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid grid1d_ixy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_ixy),     'get_sparse_geog', &
                                   'get_var grid1d_ixy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'grid1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid grid1d_jxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   grid1d_jxy),     'get_sparse_geog', &
                                   'get_var grid1d_jxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid land1d_ixy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_ixy),     'get_sparse_geog', &
                                   'get_var land1d_ixy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid land1d_jxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_jxy),     'get_sparse_geog', &
                                   'get_var land1d_jxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'land1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid land1d_wtxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   land1d_wtxy),    'get_sparse_geog', &
                                   'get_var land1d_wtxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid cols1d_ixy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_ixy),     'get_sparse_geog', &
                                   'get_var cols1d_ixy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid cols1d_jxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_jxy),     'get_sparse_geog', &
                                   'get_var cols1d_jxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid cols1d_wtxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_wtxy),    'get_sparse_geog', &
                                   'get_var cols1d_wtxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'cols1d_ityplun', VarID), 'get_sparse_geog', &
                         'inq_varid cols1d_ityplun '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   cols1d_ityplun), 'get_sparse_geog', &
                                   'get_var cols1d_ityplun '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_ixy', VarID),     'get_sparse_geog', &
                         'inq_varid pfts1d_ixy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_ixy),     'get_sparse_geog', &
                                   'get_var pfts1d_ixy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_jxy', VarID),     'get_sparse_geog', &
                         'inq_varid pfts1d_jxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_jxy),     'get_sparse_geog', &
                                   'get_var pfts1d_jxy '//trim(jules_restart_filename))

call nc_check(nf90_inq_varid(ncid, 'pfts1d_wtxy', VarID),    'get_sparse_geog', &
                         'inq_varid pfts1d_wtxy '//trim(jules_restart_filename))
call nc_check(nf90_get_var(  ncid, VarID,   pfts1d_wtxy),    'get_sparse_geog', &
                                   'get_var pfts1d_wtxy '//trim(jules_restart_filename))

if (cstat == 'close') then
   call nc_check(nf90_close(ncid),'get_sparse_geog','close '//trim(fname) )
   ncid = 0
endif

! A little sanity check

if ((debug > 7) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'Raw lat/lon information as read ...'
   write(logfileunit,*)'grid1d_ixy     range ',minval(grid1d_ixy),    maxval(grid1d_ixy)
   write(logfileunit,*)'grid1d_jxy     range ',minval(grid1d_jxy),    maxval(grid1d_jxy)

   write(logfileunit,*)'land1d_ixy     range ',minval(land1d_ixy),    maxval(land1d_ixy)
   write(logfileunit,*)'land1d_jxy     range ',minval(land1d_jxy),    maxval(land1d_jxy)
   write(logfileunit,*)'land1d_wtxy    range ',minval(land1d_wtxy),   maxval(land1d_wtxy)

   write(logfileunit,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(logfileunit,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(logfileunit,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(logfileunit,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(logfileunit,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(logfileunit,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(logfileunit,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)

   write(     *     ,*)
   write(     *     ,*)'Raw lat/lon information as read ...'
   write(     *     ,*)'grid1d_ixy     range ',minval(grid1d_ixy),    maxval(grid1d_ixy)
   write(     *     ,*)'grid1d_jxy     range ',minval(grid1d_jxy),    maxval(grid1d_jxy)

   write(     *     ,*)'land1d_ixy     range ',minval(land1d_ixy),    maxval(land1d_ixy)
   write(     *     ,*)'land1d_jxy     range ',minval(land1d_jxy),    maxval(land1d_jxy)
   write(     *     ,*)'land1d_wtxy    range ',minval(land1d_wtxy),   maxval(land1d_wtxy)

   write(     *     ,*)'cols1d_ixy     range ',minval(cols1d_ixy),    maxval(cols1d_ixy)
   write(     *     ,*)'cols1d_jxy     range ',minval(cols1d_jxy),    maxval(cols1d_jxy)
   write(     *     ,*)'cols1d_wtxy    range ',minval(cols1d_wtxy),   maxval(cols1d_wtxy)
   write(     *     ,*)'cols1d_ityplun range ',minval(cols1d_ityplun),maxval(cols1d_ityplun)

   write(     *     ,*)'pfts1d_ixy     range ',minval(pfts1d_ixy),    maxval(pfts1d_ixy)
   write(     *     ,*)'pfts1d_jxy     range ',minval(pfts1d_jxy),    maxval(pfts1d_jxy)
   write(     *     ,*)'pfts1d_wtxy    range ',minval(pfts1d_wtxy),   maxval(pfts1d_wtxy)

endif

return
end subroutine get_sparse_geog


!------------------------------------------------------------------
!> set_model_time_step
!> This defines the window used for assimilation.
!> all observations +/- half this timestep are assimilated.

function set_model_time_step()
!------------------------------------------------------------------

type(time_type) :: set_model_time_step

call error_handler(E_MSG, 'set_model_time_step', 'FIXME SHAMS routine is not tested')

! FIXME ... should check to see that time step is attainable given the JULES namelist values.

set_model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

end function set_model_time_step


!------------------------------------------------------------------
!> parse_variable_table
!>
!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:variables  variable.
!>  Each variable must have 6 entries.
!>  1: variable name
!>  2: DART KIND
!>  3: minimum value - as a character string - if none, use 'NA'
!>  4: maximum value - as a character string - if none, use 'NA'
!>  5: what file contains the variable - '.r. => restart', '.h0. => h0 history file'
!>  6: does the variable get updated in the restart file or not ...
!>     only variables from restart files may be updated.
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file
!>     all these variables will be updated INTERNALLY IN DART
!>     only variables marked '.r', 'UPDATE' will be modified for jules.
!>
!>  The calling code should check to see if the variable exists.

function parse_variable_table() result(ngood)

integer :: ngood
! character variables(:)        is module scope
! character variable_table(:,:) is module scope

integer :: i
character(len=NF90_MAX_NAME) :: varname       ! column 1
character(len=NF90_MAX_NAME) :: dartstr       ! column 2
character(len=NF90_MAX_NAME) :: minvalstring  ! column 3
character(len=NF90_MAX_NAME) :: maxvalstring  ! column 4
character(len=NF90_MAX_NAME) :: origin_file   ! column 5
character(len=NF90_MAX_NAME) :: state_or_aux  ! column 6

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do i = 1, max_state_variables

   varname      = trim(variables(num_state_table_columns*i - 5))
   dartstr      = trim(variables(num_state_table_columns*i - 4))
   minvalstring = trim(variables(num_state_table_columns*i - 3))
   maxvalstring = trim(variables(num_state_table_columns*i - 2))
   origin_file  = trim(variables(num_state_table_columns*i - 1))
   state_or_aux = trim(variables(num_state_table_columns*i    ))

   call to_upper(origin_file)
   call to_upper(state_or_aux)

   variable_table(i,VT_VARNAMEINDX) = trim(varname)
   variable_table(i,VT_KINDINDX)    = trim(dartstr)
   variable_table(i,VT_MINVALINDX)  = trim(minvalstring)
   variable_table(i,VT_MAXVALINDX)  = trim(maxvalstring)
   variable_table(i,VT_ORIGININDX)  = trim(origin_file)
   variable_table(i,VT_STATEINDX)   = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( variable_table(i,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(variable_table(i,:) == ' ') ) then
      string1 = 'input.nml &model_nml:clm_variables not fully specified'
      string2 = 'must be 6 entries per variable. Last known variable name is'
      string3 = '['//trim(variable_table(i,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 8) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', &
                                               trim(variable_table(i,2)), ' ', &
                                               trim(variable_table(i,3)), ' ', &
                                               trim(variable_table(i,4)), ' ', &
                                               trim(variable_table(i,5)), ' ', &
                                               trim(variable_table(i,6))
      write(     *     ,*)'variable ',i,' is ',trim(variable_table(i,1)), ' ', &
                                               trim(variable_table(i,2)), ' ', &
                                               trim(variable_table(i,3)), ' ', &
                                               trim(variable_table(i,4)), ' ', &
                                               trim(variable_table(i,5)), ' ', &
                                               trim(variable_table(i,6))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == max_state_variables) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_table',string1,text2=string2)
endif

end function parse_variable_table


!------------------------------------------------------------------
!> define_var_dims()
!> takes the N-dimensional variable and appends the DART
!> dimensions of 'copy' and 'time'. If the variable initially had a 'time'
!> dimension, it is ignored because (by construction) it is a singleton
!> dimension.

subroutine define_var_dims(ivar, ncid, memberdimid, unlimiteddimid, ndims, dimids)

integer,               intent(in)  :: ivar, ncid, memberdimid, unlimiteddimid
integer,               intent(out) :: ndims
integer, dimension(:), intent(out) :: dimids

character(len=NF90_MAX_NAME),dimension(NF90_MAX_VAR_DIMS) :: dimnames

integer :: i, mydimid

ndims = 0

DIMLOOP : do i = 1,progvar(ivar)%numdims

   if (progvar(ivar)%dimnames(i) == 'time') cycle DIMLOOP

   call nc_check(nf90_inq_dimid(ncid=ncid, name=progvar(ivar)%dimnames(i), dimid=mydimid), &
                           'define_var_dims','inq_dimid '//trim(progvar(ivar)%dimnames(i)))

   ndims         = ndims + 1
   dimids(ndims) = mydimid
   dimnames(ndims) = progvar(ivar)%dimnames(i)

enddo DIMLOOP

! The last two dimensions are always 'copy' and 'time'
ndims           = ndims + 1
dimids(ndims)   = memberdimid
dimnames(ndims) = 'copy'
ndims           = ndims + 1
dimids(ndims)   = unlimitedDimid
dimnames(ndims) = 'time'

if ((debug > 9) .and. do_output()) then

   write(logfileunit,*)
   write(logfileunit,*)'define_var_dims knowledge'
   write(logfileunit,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(logfileunit,*)'thus dimids ',dimids(1:ndims)

   write(     *     ,*)
   write(     *     ,*)'define_var_dims knowledge'
   write(     *     ,*)trim(progvar(ivar)%varname),' repackaging into dimnames: ', &
                       (/ (trim(dimnames(i))//' ',i=1,ndims) /)
   write(     *     ,*)'thus dimids ',dimids(1:ndims)

endif

return
end subroutine define_var_dims


!------------------------------------------------------------------
!> findVarIndex
!>
!> finds the index into the progvar structure for the named variable


function findVarIndex(varstring, caller)
character(len=*), intent(in) :: varstring
character(len=*), intent(in) :: caller
integer                      :: findVarIndex

integer :: i

findVarIndex = -1

! Skip to the right variable
VARTYPES : do i = 1,nfields
    findVarIndex = i
    if ( trim(progvar(i)%varname) == varstring) exit VARTYPES
enddo VARTYPES

if (findVarIndex < 1) then
   write(string1,*) trim(caller)//' cannot find "'//trim(varstring)//'" in list of DART state variables.'
   call error_handler(E_ERR,'findVarIndex',string1,source,revision,revdate)
endif

end function findVarIndex


!------------------------------------------------------------------
!> SetLocatorArrays
!>
!> This function will create the relational table that will indicate how many
!> and which columns pertain to the gridcells. A companion function will
!> return the column indices that are needed to recreate the gridcell value.
!>
!> This fills the gridCellInfo(:,:) structure.
!> given a gridcell, the gridCellInfo(:,:) structure will indicate how many and
!> which columns are part of the gridcell.

subroutine SetLocatorArrays()

integer :: ilon, ilat, ij, iunit
integer :: icol, currenticol(Nlon,Nlat)
integer :: ipft, currentipft(Nlon,Nlat)

call error_handler(E_ERR, 'SetLocatorArrays', 'FIXME routine not written', source, revision, revdate)

! TJH gridCellInfo(:,:)%ncols = 0
! TJH gridCellInfo(:,:)%Ntiles = 0
! TJH 
! TJH ! Count up how many columns are in each gridcell
! TJH 
! TJH do ij = 1,Ncolumn
! TJH    ilon = cols1d_ixy(ij)
! TJH    ilat = cols1d_jxy(ij)
! TJH    gridCellInfo(ilon,ilat)%ncols = gridCellInfo(ilon,ilat)%ncols + 1
! TJH enddo
! TJH 
! TJH ! Count up how many pfts are in each gridcell
! TJH 
! TJH do ij = 1,Ntile
! TJH    ilon = pfts1d_ixy(ij)
! TJH    ilat = pfts1d_jxy(ij)
! TJH    gridCellInfo(ilon,ilat)%Ntiles = gridCellInfo(ilon,ilat)%Ntiles + 1
! TJH enddo
! TJH 
! TJH ! Create storage for the list of column,pft indices
! TJH 
! TJH do ilon = 1,Nlon
! TJH do ilat = 1,Nlat
! TJH    if ( gridCellInfo(ilon,ilat)%ncols > 0 ) then
! TJH       allocate( gridCellInfo(ilon,ilat)%columnids( gridCellInfo(ilon,ilat)%ncols ))
! TJH    endif
! TJH    if ( gridCellInfo(ilon,ilat)%Ntiles > 0 ) then
! TJH       allocate( gridCellInfo(ilon,ilat)%tileIds( gridCellInfo(ilon,ilat)%Ntiles ))
! TJH    endif
! TJH enddo
! TJH enddo
! TJH 
! TJH ! Fill the column pointer arrays
! TJH 
! TJH currenticol(:,:) = 0
! TJH do ij = 1,Ncolumn
! TJH 
! TJH    ilon = cols1d_ixy(ij)
! TJH    ilat = cols1d_jxy(ij)
! TJH 
! TJH    currenticol(ilon,ilat) = currenticol(ilon,ilat) + 1
! TJH    icol = currenticol(ilon,ilat)
! TJH 
! TJH    if ( icol <= gridCellInfo(ilon,ilat)%ncols ) then
! TJH       gridCellInfo(ilon,ilat)%columnids(icol) = ij
! TJH    else
! TJH       write(string1,'(''gridcell('',i4,'','',i4,'') has at most '',i4,'' columns.'')') &
! TJH          ilon, ilat, gridCellInfo(ilon,ilat)%ncols
! TJH       write(string2,'(''Found '',i8,'' at dart index '',i12)') icol, ij
! TJH       call error_handler(E_ERR, 'SetLocatorArrays', string1, &
! TJH                    source, revision, revdate, text2=string2)
! TJH    endif
! TJH enddo
! TJH 
! TJH ! Fill the pft pointer arrays
! TJH 
! TJH currentipft(:,:) = 0
! TJH do ij = 1,Ntile
! TJH 
! TJH    ilon = pfts1d_ixy(ij)
! TJH    ilat = pfts1d_jxy(ij)
! TJH 
! TJH    currentipft(ilon,ilat) = currentipft(ilon,ilat) + 1
! TJH    ipft = currentipft(ilon,ilat)
! TJH 
! TJH    if ( ipft <= gridCellInfo(ilon,ilat)%Ntiles ) then
! TJH       gridCellInfo(ilon,ilat)%tileIds(ipft) = ij
! TJH    else
! TJH       write(string1,'(''gridcell('',i4,'','',i4,'') has at most '',i4,'' pfts.'')') &
! TJH          ilon, ilat, gridCellInfo(ilon,ilat)%Ntiles
! TJH       write(string2,'(''Found '',i8,'' at dart index '',i12)') ipft, ij
! TJH       call error_handler(E_ERR, 'SetLocatorArrays', string1, &
! TJH                    source, revision, revdate, text2=string2)
! TJH    endif
! TJH enddo
! TJH 
! TJH ! Check block
! TJH 
! TJH if ((debug > 99) .and. do_output()) then
! TJH 
! TJH    iunit = open_file('gridcell_column_table.txt',form='formatted',action='write')
! TJH 
! TJH    do ilon = 1,Nlon
! TJH    do ilat = 1,Nlat
! TJH       if (gridCellInfo(ilon,ilat)%ncols > 0) then
! TJH          write(iunit,'(''gridcell'',i8,1x,i8,'' has '', i6, '' columns:'')') &
! TJH                    ilon,ilat,gridCellInfo(ilon,ilat)%ncols
! TJH          write(iunit,*)gridCellInfo(ilon,ilat)%columnids
! TJH       endif
! TJH    enddo
! TJH    enddo
! TJH 
! TJH    call close_file(iunit)
! TJH 
! TJH    iunit = open_file('gridcell_pft_table.txt',form='formatted',action='write')
! TJH 
! TJH    do ilon = 1,Nlon
! TJH    do ilat = 1,Nlat
! TJH       if (gridCellInfo(ilon,ilat)%Ntiles > 0) then
! TJH          write(iunit,'(''gridcell'',i8,1x,i8,'' has '', i6, '' pfts : '')') &
! TJH                    ilon,ilat,gridCellInfo(ilon,ilat)%Ntiles
! TJH          write(iunit,*)gridCellInfo(ilon,ilat)%tileIds
! TJH       endif
! TJH    enddo
! TJH    enddo
! TJH 
! TJH    call close_file(iunit)
! TJH 
! TJH endif

end subroutine SetLocatorArrays


!------------------------------------------------------------------
!> SetVariableAttributes
!> converts the information in the variable_table 
!> to the progvar structure for each variable.
!> If the numerical limit does not apply, it is set to MISSING_R8, even if
!> it is the maximum that does not apply.

subroutine SetVariableAttributes(ivar)

integer, intent(in) :: ivar

integer  :: ios
real(r8) :: minvalue, maxvalue

progvar(ivar)%varname     = trim(variable_table(ivar,VT_VARNAMEINDX))
progvar(ivar)%kind_string = trim(variable_table(ivar,VT_KINDINDX))
progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
progvar(ivar)%maxlevels   = 0
progvar(ivar)%dimlens     = 0
progvar(ivar)%dimnames    = ' '
progvar(ivar)%spvalINT    = -9999        ! from CESM/jules jules_varcon.F90
progvar(ivar)%spvalR4     = 1.e36_r4     ! from CESM/jules jules_varcon.F90
progvar(ivar)%spvalR8     = 1.e36_r8     ! from CESM/jules jules_varcon.F90
progvar(ivar)%missingINT  = MISSING_I
progvar(ivar)%missingR4   = MISSING_R4
progvar(ivar)%missingR8   = MISSING_R8
progvar(ivar)%rangeRestricted   = BOUNDED_NONE
progvar(ivar)%minvalue          = MISSING_R8
progvar(ivar)%maxvalue          = MISSING_R8
progvar(ivar)%has_fill_value    = .true.
progvar(ivar)%has_missing_value = .true.
progvar(ivar)%update            = .false.

if (variable_table(ivar,VT_ORIGININDX) == 'RESTART') then
   progvar(ivar)%origin = trim(jules_restart_filename)
else
   variable_table(ivar,VT_ORIGININDX) = 'OUTPUT'
   progvar(ivar)%origin = trim(jules_output_filename)
endif

if ((variable_table(ivar,VT_STATEINDX)  == 'UPDATE') .and. &
    (variable_table(ivar,VT_ORIGININDX) == 'RESTART')) progvar(ivar)%update = .true.

! set the default values

minvalue = MISSING_R8
maxvalue = MISSING_R8
progvar(ivar)%minvalue = MISSING_R8
progvar(ivar)%maxvalue = MISSING_R8

! If the character string can be interpreted as an r8, great.
! If not, there is no value to be used.

read(variable_table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
if (ios == 0) progvar(ivar)%minvalue = minvalue

read(variable_table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
if (ios == 0) progvar(ivar)%maxvalue = maxvalue

! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
       (progvar(ivar)%maxvalue /= MISSING_R8) ) then
   progvar(ivar)%rangeRestricted = BOUNDED_BOTH

elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

elseif (progvar(ivar)%minvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_BELOW

else
   progvar(ivar)%rangeRestricted = BOUNDED_NONE

endif

! Check to make sure min is less than max if both are specified.

if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
   if (maxvalue < minvalue) then
      write(string1,*)'&model_nml state_variable input error for ',trim(progvar(ivar)%varname)
      write(string2,*)'minimum value (',minvalue,') must be less than '
      write(string3,*)'maximum value (',maxvalue,')'
      call error_handler(E_ERR,'SetVariableAttributes',string1, &
         source,revision,revdate,text2=trim(string2),text3=trim(string3))
   endif
endif

end subroutine SetVariableAttributes


!------------------------------------------------------------------
!> FindDesiredTimeIndx
!> returns the index into the time array that matches
!> the model_time from the jules restart file.

function FindDesiredTimeIndx(ncid, ntimes, varname)

integer,          intent(in) :: ncid
integer,          intent(in) :: ntimes
character(len=*), intent(in) :: varname
integer                      :: FindDesiredTimeIndx

integer :: VarID
real(r8),  dimension(ntimes) :: mytimes
character(len=NF90_MAX_NAME) :: attvalue

type(time_type) :: thistime, basetime
integer :: ios, itime, basedays, baseseconds
integer :: iyear, imonth, iday, ihour, imin, isec

FindDesiredTimeIndx = MISSING_I   ! initialize to failure setting

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'FindDesiredTimeIndx:', 'inq_varid time '//varname)
call nc_check(nf90_get_var(  ncid, VarID, mytimes), &
        'FindDesiredTimeIndx:', 'get_var   time '//varname)
call nc_check(nf90_get_att(ncid, VarID, 'units', attvalue), &
        'FindDesiredTimeIndx:', 'time get_att units '//varname)

! time:units = "seconds since 2004-01-01 00:00:00" ;
!               12345678901234

if (attvalue(1:13) /= 'seconds since') then
   write(string1,*)'expecting time units of [seconds since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2)
endif

read(attvalue,'(14x,i4,5(1x,i2))',iostat=ios)iyear,imonth,iday,ihour,imin,isec
if (ios /= 0) then
   write(string1,*)'Unable to read time units. Error status was ',ios
   write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

basetime = set_date(iyear, imonth, iday, ihour, imin, isec)

! convert each time to a DART time and compare to desired

TIMELOOP : do itime = 1,ntimes

   isec     = int(mytimes(itime))
   thistime = set_time(isec, 0) + basetime

   if (thistime == model_time) then
      FindDesiredTimeIndx = itime
      exit TIMELOOP
   endif

enddo TIMELOOP

! FIXME ... do we actually need a perfect match ... or do we just use the last one
if ( FindDesiredTimeIndx == MISSING_I ) then
   call print_time(model_time,str='model time is ',iunit=logfileunit)
   call print_time(model_time,str='model time is ')
   call print_date(model_time,str='model date is ',iunit=logfileunit)
   call print_date(model_time,str='model date is ')
   write(string1,*)'No time matching model_time found for '//trim(varname)
   call error_handler(E_ERR, 'FindDesiredTimeIndx:', string1, &
          source, revision, revdate )
endif

if ((debug > 0) .and. do_output()) then
   write(string1,*)trim(varname)//' matching time index is ',FindDesiredTimeIndx
   call error_handler(E_MSG, 'FindDesiredTimeIndx:', string1)
endif

end function FindDesiredTimeIndx


!------------------------------------------------------------------
!> Read_Jules_Soil_Namelist
!> reads the namelist specifying the number of soil layers and the
!> associated layer thicknesses. DART wants the soil layer depth,
!> so the thicknesses are converted to depths.

subroutine Read_Jules_Soil_Namelist()

! integer,  intent(out) :: Nsoil     ! were it not for module scope
! real(r8), intent(out) :: SOILLEVEL ! were it not for module scope

integer :: iunit, io, i
integer, PARAMETER :: MAX_SOIL_LEVELS = 200

real(r8) :: confrac
real(r8) :: dzsoil_io(MAX_SOIL_LEVELS) 
logical :: l_bedrock
logical :: l_dpsids_dsdz
logical :: l_soil_sat_down
logical :: l_vg_soil
integer :: sm_levels
integer :: soilhc_method
real(r8) :: zsmc
real(r8) :: zst

namelist /jules_soil/ confrac, dzsoil_io, l_bedrock, l_dpsids_dsdz, &
    l_soil_sat_down, l_vg_soil, sm_levels, soilhc_method, zsmc, zst

call find_namelist_in_file('jules_soil.nml', 'jules_soil', iunit)
read(iunit, nml = jules_soil, iostat = io)
call check_namelist_read(iunit, io, 'jules_soil')

if (Nsoil /= sm_levels) then
   write(string1,*)'Number of soil layers unknown.'
   write(string2,*)'Number of soil layers from namelist is ',sm_levels
   write(string3,*)'Number of soil layers from netCDF   is ',Nsoil
   call error_handler(E_ERR, 'Read_Jules_Soil_Namelist:', string1, &
          source, revision, revdate, text2=string2, text3=string3 )
endif

! Add up the thicknesses 

SOILLEVEL(1) = dzsoil_io(1)
do i = 2,Nsoil
   SOILLEVEL(i) = SOILLEVEL(i-1) + dzsoil_io(i)
enddo

end subroutine Read_Jules_Soil_Namelist


!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
