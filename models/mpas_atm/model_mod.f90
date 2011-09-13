! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is the interface between the model model and DART.

! Modules that are absolutely required for use are listed
use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI
use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      & 
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, loc_get_close_obs => get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text, close_file

use     obs_kind_mod, only : paramname_length,        &
                             get_raw_obs_kind_index,  &
                             get_raw_obs_kind_name

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

public :: get_model_analysis_filename,  &
          analysis_file_to_statevector, &
          statevector_to_analysis_file, &
          get_base_time,                &
          get_state_time

! version controlled file description for error handling, do not edit

character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

character(len=256) :: string1, string2
logical, save :: module_initialized = .false.

! Storage for a random sequence for perturbing a single initial state

type(random_seq_type) :: random_seq

! things which can/should be in the model_nml

integer            :: assimilation_period_days = 0
integer            :: assimilation_period_seconds = 60
real(r8)           :: model_perturbation_amplitude = 0.2
logical            :: output_state_vector = .true.
integer            :: debug = 0   ! turn up for more and more debug messages
character(len=32)  :: calendar = 'Gregorian'
character(len=256) :: model_analysis_filename = 'mpas_restart.nc'

namelist /model_nml/  &
   model_analysis_filename,        &
   output_state_vector,         &
   assimilation_period_days,    &  ! for now, this is the timestep
   assimilation_period_seconds, &
   model_perturbation_amplitude,&
   calendar,                    &
   debug

!------------------------------------------------------------------
! DART state vector are specified in the input.nml:mpas_vars_nml namelist.
!------------------------------------------------------------------

integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 2
character(len=NF90_MAX_NAME) :: mpas_state_variables(max_state_variables * num_state_table_columns ) = ' '
character(len=NF90_MAX_NAME) :: variable_table(max_state_variables, num_state_table_columns )

namelist /mpas_vars_nml/ mpas_state_variables

integer :: nfields

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: xtype
   integer :: numdims
   integer :: numvertical
   integer :: numcells
   logical :: ZonHalf     ! vertical coordinate has dimension nVertLevels
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=paramname_length) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! Grid parameters - the values will be read from an mpas analysis file. 

integer :: nCells=-1, nVertices=-1, nEdges=-1   ! horizontal grid dimensions
integer :: nVertLevels=-1, nVertLevelsP1=-1     ! vertical grid dimensions
integer :: maxEdges=-1, vertexDegree=-1         ! unstructured mesh
integer :: nSoilLevels=-1                       ! soil layers

! scalar grid positions

real(r8), allocatable :: lonCell(:) ! cell center longitudes (degrees)
real(r8), allocatable :: latCell(:) ! cell center latitudes  (degrees)
real(r8), allocatable :: zgrid(:,:) ! cell center geometric height at cell centers (ncells,nvert)
integer,  allocatable :: CellsOnVertex(:,:) ! list of cell centers defining a triangle

integer               :: model_size      ! the state vector length
type(time_type)       :: model_time      ! valid time of the model state
type(time_type)       :: model_timestep  ! smallest time to adv model
real(r8), allocatable :: ens_mean(:)     ! may be needed for forward ops

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.

logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

!------------------------------------------------------------------
! The model analysis manager namelist variables
!------------------------------------------------------------------

character(len= 64) :: ew_boundary_type, ns_boundary_type

INTERFACE vector_to_prog_var
      MODULE PROCEDURE vector_to_1d_prog_var
      MODULE PROCEDURE vector_to_2d_prog_var
      MODULE PROCEDURE vector_to_3d_prog_var
      MODULE PROCEDURE vector_to_4d_prog_var
END INTERFACE

INTERFACE get_base_time
      MODULE PROCEDURE get_base_time_ncid
      MODULE PROCEDURE get_base_time_fname
END INTERFACE

contains

!==================================================================
! All the REQUIRED interfaces come first - just by convention.
!==================================================================


function get_model_size()
!------------------------------------------------------------------
! Done - TJH.
! Returns the size of the model as an integer. 
! Required for all applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! Done - TJH.
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called IF the namelist parameter
! async is set to 0 in perfect_model_obs or filter -OR- if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If none of these options
! are used (the model will only be advanced as a separate 
! model-specific executable), this can be a NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
   call print_time(time,'NULL interface adv_1step (no advance) DART time is')
   call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

! FIXME: put an error handler call here - we cannot advance the model
! this way and it would be an error if filter called it.

end subroutine adv_1step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
! FIXME ... just needs checking ....
! given an index into the state vector, return its location and
! if given, the var kind.   despite the name, var_type is a generic
! kind, like those in obs_kind/obs_kind_mod.f90, starting with KIND_

integer, intent(in)            :: index_in
type(location_type)            :: location
integer, optional, intent(out) :: var_type
  
! Local variables

integer  :: nxp, nzp, iloc, kloc, nf, n
integer  :: myindx, cell_index
real(r8) :: height

if ( .not. module_initialized ) call static_init_model

myindx = -1
nf     = -1

! Determine the right variable
FindIndex : DO n = 1,nfields
    IF( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) THEN
      nf = n
      myindx = index_in - progvar(n)%index1 + 1
      EXIT FindIndex
    ENDIF
ENDDO FindIndex

IF( myindx == -1 ) THEN
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
ENDIF

! Now that we know the variable, find the cell 

nxp = progvar(nf)%numcells
nzp = progvar(nf)%numvertical

iloc   = 1 + (myindx-1) / nzp  ! cell index
kloc = myindx - (iloc-1)*nzp   ! vertical level index

! If full sigma level

if ( progvar(nf)%ZonHalf ) then
   height = (zgrid(iloc,kloc) + zgrid(iloc,kloc+1))*0.5
else
   height = zgrid(iloc,kloc)
endif

location = set_location(lonCell(iloc),latCell(iloc), height, VERTISHEIGHT)

if (debug > 5) then

    write(*,'("INDEX_IN / myindx / IVAR / NX, NZ: ",2(i10,2x),3(i5,2x))') index_in, myindx, nf, nxp, nzp
    write(*,'("                       ILOC, KLOC: ",2(i5,2x))') iloc, kloc
    write(*,'("                      LON/LAT/HGT: ",3(f12.3,2x))') lonCell(iloc), latCell(iloc), height

endif

if (present(var_type)) then
     var_type = progvar(nf)%dart_kind
endif

end subroutine get_state_meta_data



subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

!############################################################################
!
!     PURPOSE:
!
!     For a given lat, lon, and height, interpolate the correct state value
!     to that location for the filter from the model state vectors
!
!       ERROR codes:
!
!       ISTATUS = 99:  general error in case something terrible goes wrong...
!       ISTATUS = 15:  dont know what to do with vertical coord supplied
!       ISTATUS = 16:  lonCell illegal
!       ISTATUS = 17:  latCell illegal
!       ISTATUS = 18:  altitude illegal
!
!############################################################################

! Passed variables

  real(r8),            intent(in)  :: x(:)
  type(location_type), intent(in)  :: location
  integer,             intent(in)  :: obs_type
  real(r8),            intent(out) :: interp_val
  integer,             intent(out) :: istatus

! Local storage

  real(r8)         :: loc_array(3), llon, llat, lheight, lon_fract, lat_fract, alt_fract
  integer          :: base_offset, end_offset, blon(2), blat(2), balt(2), i, j, k, ier
  real(r8)         :: cube(2, 2, 2), square(2, 2), line(2)

  IF ( .not. module_initialized ) call static_init_model

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the interp_val will be set to the 
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s

  interp_val = MISSING_R8     ! the DART bad value flag
  istatus = 99                ! unknown error

! Get the individual locations values

  loc_array = get_location(location)
  llon      = loc_array(1)
  llat      = loc_array(2)
  lheight   = loc_array(3)

  IF (debug > 2) print *, 'requesting interpolation at ', llon, llat, lheight

! Only height for vertical location type is supported at this point
  IF(.not. vert_is_height(location) ) THEN 
     istatus = 15
     return
  ENDIF

! Find the start and end offsets for this field in the state vector x(:)


! FIXME ... Jeffs homework starts here ...


end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_timestep

end function get_model_time_step



subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model.
! 
! All the grid information comes from the initialization of
! the dart_model_mod module.

! Local variables - all the important ones have module scope


integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=paramname_length)       :: kind_string
integer :: ncid, VarID, numdims, varsize, dimlen
integer :: iunit, io, ivar, i, index1, indexN
integer :: ss, dd
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID

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

! Read the MPAS variable list to populate DART state vector
! Once parsed, the values will be recorded for posterity
call find_namelist_in_file('input.nml', 'mpas_vars_nml', iunit)
read(iunit, nml = mpas_vars_nml, iostat = io)
call check_namelist_read(iunit, io, 'mpas_vars_nml')

!---------------------------------------------------------------
! Set the time step ... causes mpas namelists to be read.
! Ensures model_timestep is multiple of 'dynamics_timestep'

call set_calendar_type( calendar )   ! comes from model_mod_nml

model_timestep = set_model_time_step()

call get_time(model_timestep,ss,dd) ! set_time() assures the seconds [0,86400)

write(string1,*)'assimilation period is ',dd,' days ',ss,' seconds'
call error_handler(E_MSG,'static_init_model',string1,source,revision,revdate)

!---------------------------------------------------------------
! 1) get grid dimensions
! 2) allocate space for the grids 
! 3) read them from the analysis file

! FIXME: Do we need to define nVertLevelsP1 in get_grid_dims for staggered grids?

call get_grid_dims(nCells, nVertices, nVertLevels, nVertLevelsP1, vertexDegree, nSoilLevels)

allocate(latCell(nCells), lonCell(nCells)) 
allocate(zgrid(nCells, nVertLevelsP1))
allocate(CellsOnVertex(nVertices, vertexDegree))

call get_grid(nCells, nVertices, nVertLevels, vertexDegree, &
              latCell, lonCell, zgrid, CellsOnVertex)
              
!---------------------------------------------------------------
! Compile the list of model variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the model analysis file
!
! Compute the offsets into the state vector for the start of each
! different variable type. Requires reading shapes from the model
! analysis file. As long as TIME is the LAST dimension, we're OK.
!
! Record the extent of the data type in the state vector.


call nc_check( nf90_open(trim(model_analysis_filename), NF90_NOWRITE, ncid), &
                  'static_init_model', 'open '//trim(model_analysis_filename))

call verify_state_variables( mpas_state_variables, ncid, model_analysis_filename, &
                             nfields, variable_table )

TimeDimID = FindTimeDimension( ncid )

if (TimeDimID < 0 ) then
   write(string1,*)'unable to find a dimension named xtime.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

call nc_check(nf90_Inquire(ncid,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                    'static_init_model', 'inquire '//trim(model_analysis_filename))

if ( (TimeDimID > 0) .and. (unlimitedDimID > 0) .and. (TimeDimID /= unlimitedDimID)) then
   write(string1,*)'IF Time is not the unlimited dimension, I am lost.'
   call error_handler(E_MSG,'static_init_model', string1, source, revision, revdate)
endif

index1  = 1;
indexN  = 0;
do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0

   string2 = trim(model_analysis_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'long_name' , progvar(ivar)%long_name), &
            'static_init_model', 'get_att long_name '//trim(string2))

   call nc_check( nf90_get_att(ncid, VarId, 'units' , progvar(ivar)%units), &
            'static_init_model', 'get_att units '//trim(string2))

   call nc_check(nf90_inquire_variable(ncid, VarId, dimids=dimIDs, ndims=progvar(ivar)%numdims), &
            'static_init_model', 'inquire '//trim(string2))

   ! TIME is the unlimited dimension, and in Fortran, this is always the
   ! LAST dimension in this loop. Since we are not concerned with it, we
   ! need to skip it. 
   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      if (dimIDs(i) == TimeDimID) cycle DimensionLoop

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
                                          'static_init_model', string1)
      progvar(ivar)%dimlens(i) = dimlen
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ( debug > 0 ) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens(1:progvar(ivar)%numdims)
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

enddo

call nc_check( nf90_close(ncid), &
                  'static_init_model', 'close '//trim(model_analysis_filename))

model_size = progvar(nfields)%indexN

if ( debug > 0 ) then
  write(logfileunit,'("grid: nCells, nVertices, nVertLevels =",3(1x,i6))') nCells, nVertices, nVertLevels
  write(     *     ,'("grid: nCells, nVertices, nVertLevels =",3(1x,i6))') nCells, nVertices, nVertLevels
  write(logfileunit, *)'model_size = ', model_size
  write(     *     , *)'model_size = ', model_size
endif

allocate( ens_mean(model_size) )


end subroutine static_init_model



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

if (allocated(latCell))        deallocate(latCell)
if (allocated(lonCell))        deallocate(lonCell)
if (allocated(zgrid))          deallocate(zgrid)
if (allocated(CellsOnVertex))  deallocate(CellsOnVertex)

end subroutine end_model



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

if ( .not. module_initialized ) call static_init_model

! for now, just set to 0
time = set_time(0,0)

! FIXME: put an error handler call here - we cannot initialize the model
! this way and it would be an error if filter called it.

end subroutine init_time



subroutine init_conditions(x)
!------------------------------------------------------------------
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model
 
x = 0.0_r8

! FIXME: put an error handler call here - we cannot initialize the model
! this way and it would be an error if filter called it.

end subroutine init_conditions



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
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
integer :: NCellsDimID
integer :: NVertLevelsDimID
integer :: NVertLevelsP1DimID

! for the prognostic variables
integer :: ivar, VarID

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_model_namelist

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

call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
                           'nc_write_model_atts','inq_dimid NMLlinelen')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
                           'nc_write_model_atts', 'copy dimid '//trim(filename))
call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
                           'nc_write_model_atts', 'time dimid '//trim(filename))

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)'Time Dimension ID ',TimeDimID, &
             ' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
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
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'model' ), &
           'nc_write_model_atts', 'model put '//trim(filename))

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
   call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
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
   
   call nc_check(nf90_def_dim(ncid=ncFileID, name='nCells', &
          len = nCells, dimid = NCellsDimID),'nc_write_model_atts', 'nCells def_dim '//trim(filename))

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nVertLevelsP1', &
          len = nVertLevelsP1, dimid = NVertLevelsP1DimID),'nc_write_model_atts', &
                                      'nVertLevelsP1 def_dim '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Coordinate Variables and the Attributes
   !----------------------------------------------------------------------------

   ! Grid Longitudes
   call nc_check(nf90_def_var(ncFileID,name='lonCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'lonCell def_var '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, VarID, 'type', 'x1d'),  &
!                'nc_write_model_atts', 'lonCell type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center longitudes'), &
                 'nc_write_model_atts', 'lonCell long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_east'), &
                 'nc_write_model_atts', 'lonCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'valid_range', (/ 0.0_r8, 360.0_r8 /)), &
                 'nc_write_model_atts', 'lonCell valid_range '//trim(filename))

   ! Grid Latitudes
   call nc_check(nf90_def_var(ncFileID,name='latCell', xtype=nf90_double, &
                 dimids=nCellsDimID, varid=VarID),&
                 'nc_write_model_atts', 'latCell def_var '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID, VarID, 'type', 'y1d'),  &
!                'nc_write_model_atts', 'latCell type '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'cell center latitudes'), &
                 'nc_write_model_atts', 'latCell long_name '//trim(filename))
!  call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Y'),   &
!                'nc_write_model_atts', 'latCell cartesian_axis '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'degrees_north'),  &
                 'nc_write_model_atts', 'latCell units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID,'valid_range',(/ -90.0_r8, 90.0_r8 /)), &
                 'nc_write_model_atts', 'latCell valid_range '//trim(filename))

   ! Grid vertical information
   call nc_check(nf90_def_var(ncFileID,name='zgrid',xtype=nf90_real, &
                 dimids=(/ nVertLevelsP1DimID, nCellsDimID /) ,varid=VarID), &
                 'nc_write_model_atts', 'zgrid def_var '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'long_name', 'grid zgrid'), &
                 'nc_write_model_atts', 'zgrid long_name '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'units', 'meters'),  &
                 'nc_write_model_atts', 'zgrid units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'positive', 'up'),  &
                 'nc_write_model_atts', 'zgrid units '//trim(filename))
   call nc_check(nf90_put_att(ncFileID,  VarID, 'cartesian_axis', 'Z'),   &
                 'nc_write_model_atts', 'zgrid cartesian_axis '//trim(filename))

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and the Attributes
   !----------------------------------------------------------------------------

   do ivar=1, nfields

      varname = trim(progvar(ivar)%varname)
      string1 = trim(filename)//' '//trim(varname)

      ! match shape of the variable to the dimension IDs

 !    call define_var_dims(progvar(ivar), myndims, mydimids, MemberDimID, unlimitedDimID, &
 !                    NLONDimID, NLATDimID, NALTDimID, NgridLon, NgridLat, NgridAlt)

      ! define the variable and set the attributes

      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(varname), xtype=progvar(ivar)%xtype, &
                    dimids = mydimids(1:myndims), varid=VarID),&
                    'nc_write_model_atts', trim(string1)//' def_var' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(progvar(ivar)%long_name)), &
           'nc_write_model_atts', trim(string1)//' put_att long_name' )

      call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(progvar(ivar)%kind_string)), &
           'nc_write_model_atts', trim(string1)//' put_att dart_kind' )
      call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(progvar(ivar)%units)), &
           'nc_write_model_atts', trim(string1)//' put_att units' )

   enddo

   !----------------------------------------------------------------------------
   ! Finished with dimension/variable definitions, must end 'define' mode to fill.
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncfileID), 'prognostic enddef '//trim(filename))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(NF90_inq_varid(ncFileID, 'lonCell', VarID), &
                 'nc_write_model_atts', 'lonCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, lonCell ), &
                'nc_write_model_atts', 'lonCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'latCell', VarID), &
                 'nc_write_model_atts', 'latCell inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, latCell ), &
                'nc_write_model_atts', 'latCell put_var '//trim(filename))

   call nc_check(NF90_inq_varid(ncFileID, 'zgrid', VarID), &
                 'nc_write_model_atts', 'zgrid inq_varid '//trim(filename))
   call nc_check(nf90_put_var(ncFileID, VarID, zgrid ), &
                'nc_write_model_atts', 'zgrid put_var '//trim(filename))

endif

!-------------------------------------------------------------------------------
! Fill the variables we can
!-------------------------------------------------------------------------------

if (has_model_namelist) then
   call file_to_text('model_vars.nml', textblock)
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



function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
!
! TJH 29 Aug 2011 -- all errors are fatal, so the
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

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: state_vec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
character(len=NF90_MAX_NAME)          :: varname 
integer :: i, ivar, VarID, ncNdims, dimlen
integer :: TimeDimID, CopyDimID

real(r8), allocatable, dimension(:)       :: data_1d_array
real(r8), allocatable, dimension(:,:)     :: data_2d_array
real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array

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
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
            'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
            'nc_write_model_vars', 'inq_dimid time '//trim(filename))

if ( output_state_vector ) then

   call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
                 'nc_write_model_vars', 'state inq_varid '//trim(filename))
   call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
                 'nc_write_model_vars', 'state put_var '//trim(filename))

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   do ivar = 1,nfields  

      varname = trim(progvar(ivar)%varname)
      string2 = trim(filename)//' '//trim(varname)

      ! Ensure netCDF variable is conformable with progvar quantity.
      ! The TIME and Copy dimensions are intentionally not queried
      ! by looping over the dimensions stored in the progvar type.

      call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
            'nc_write_model_vars', 'inq_varid '//trim(string2))

      call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
            'nc_write_model_vars', 'inquire '//trim(string2))

      mystart = 1   ! These are arrays, actually
      mycount = 1
      DimCheck : do i = 1,progvar(ivar)%numdims

         write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(string2)
         call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
               'nc_write_model_vars', trim(string1))

         if ( dimlen /= progvar(ivar)%dimlens(i) ) then
            write(string1,*) trim(string2),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
            write(string2,*)' but it should be.'
            call error_handler(E_ERR, 'nc_write_model_vars', trim(string1), &
                            source, revision, revdate, text2=trim(string2))
         endif

         mycount(i) = dimlen

      enddo DimCheck

     ! FIXME - wouldn't hurt to make sure each of these match something.
     !         could then eliminate the if ncndims /= xxx checks below.

      where(dimIDs == CopyDimID) mystart = copyindex
      where(dimIDs == CopyDimID) mycount = 1
      where(dimIDs == TimeDimID) mystart = timeindex
      where(dimIDs == TimeDimID) mycount = 1

      if ( debug > 1 ) then
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
         call vector_to_prog_var(state_vec, progvar(ivar), data_1d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_1d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
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
         call vector_to_prog_var(state_vec, progvar(ivar), data_2d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_2d_array)

      elseif ( progvar(ivar)%numdims == 3) then

         if ( ncNdims /= 5 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 5 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_3d_array( progvar(ivar)%dimlens(1), &
                                 progvar(ivar)%dimlens(2), &
                                 progvar(ivar)%dimlens(3)))
         call vector_to_prog_var(state_vec, progvar(ivar), data_3d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_3d_array)

      elseif ( progvar(ivar)%numdims == 4) then

         if ( ncNdims /= 6 ) then
            write(string1,*)trim(varname),' no room for copy,time dimensions.'
            write(string2,*)'netcdf file should have 6 dimensions, has ',ncNdims
            call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                            source, revision, revdate, text2=string2)
         endif

         allocate(data_4d_array( progvar(ivar)%dimlens(1), &
                                 progvar(ivar)%dimlens(2), &
                                 progvar(ivar)%dimlens(3), &
                                 progvar(ivar)%dimlens(4)))
         call vector_to_prog_var(state_vec, progvar(ivar), data_4d_array)
         call nc_check(nf90_put_var(ncFileID, VarID, data_4d_array, &
             start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                   'nc_write_model_vars', 'put_var '//trim(string2))
         deallocate(data_4d_array)

      else

         ! FIXME put an error message here

      endif

   enddo


endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), 'nc_write_model_vars', 'sync '//trim(filename))

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
! may do so by adding a perturbation to each model state 
! variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer :: i
logical, save :: random_seq_init = .false.

if ( .not. module_initialized ) call static_init_model

interf_provided = .true.

! Initialize my random number sequence
if(.not. random_seq_init) then
   call init_random_seq(random_seq, my_task_id())
   random_seq_init = .true.
endif

! add some uncertainty to each ...
do i=1,size(state)
   pert_state(i) = random_gaussian(random_seq, state(i), &
                                   model_perturbation_amplitude)
enddo

end subroutine pert_model_state



subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, &
                         obs_loc, obs_kind, num_close, close_ind, dist)
!------------------------------------------------------------------

! Given a DART location (referred to as "base") and a set of candidate
! locations & kinds (obs, obs_kind), returns the subset close to the
! "base", their indices, and their distances to the "base" ...

! For vertical distance computations, general philosophy is to convert all
! vertical coordinates to a common coordinate. This coordinate type is defined
! in the namelist with the variable "vert_localization_coord".

type(get_close_type),              intent(in)    :: gc
type(location_type),               intent(inout) :: base_obs_loc
integer,                           intent(in)    :: base_obs_kind
type(location_type), dimension(:), intent(inout) :: obs_loc
integer,             dimension(:), intent(in)    :: obs_kind
integer,                           intent(out)   :: num_close
integer,             dimension(:), intent(out)   :: close_ind
real(r8),            dimension(:), intent(out)   :: dist

integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc

! Initialize variables to missing status

num_close = 0
close_ind = -99
dist      = 1.0e9   !something big and positive (far away)
istatus1  = 0
istatus2  = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_obs_loc)
base_which = nint(query_location(base_obs_loc))

! fixme ... 
if (.not. horiz_dist_only) then
!  if (base_which /= wrf%dom(1)%vert_coord) then
!     call vert_interpolate(ens_mean, base_obs_loc, base_obs_kind, istatus1)
!  elseif (base_array(3) == MISSING_R8) then
!     istatus1 = 1
!  endif
endif

if (istatus1 == 0) then

   ! Loop over potentially close subset of obs priors or state variables
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for obs_loc).
   call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind)

   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = obs_loc(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if (.not. horiz_dist_only) then
 !fixme       if (local_obs_which /= wrf%dom(1)%vert_coord) then
 !fixme           call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            obs_loc(t_ind) = local_obs_loc
 !fixme        endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is missing
      ! or vert_interpolate returned error (istatus2=1)
      local_obs_array = get_location(local_obs_loc)
      if (( (.not. horiz_dist_only)             .and. &
            (local_obs_array(3) == MISSING_R8)) .or.  &
            (istatus2 == 1)                   ) then
            dist(k) = 1.0e9
      else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      endif

   enddo
endif

end subroutine get_close_obs



subroutine ens_mean_for_model(filter_ens_mean)
!------------------------------------------------------------------
! If needed by the model interface, this is the current mean
! for all state vector items across all ensembles.

real(r8), intent(in) :: filter_ens_mean(:)

if ( .not. module_initialized ) call static_init_model

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!==================================================================
! The (model-specific) optional interfaces come next
!==================================================================


subroutine get_model_analysis_filename( filename )

character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call static_init_model

filename = trim(model_analysis_filename)

end subroutine get_model_analysis_filename



subroutine analysis_file_to_statevector(dirname, state_vector, model_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a model analysis
! file and packs them into a dart state vector.

character(len=*), intent(in)  :: dirname 
real(r8),         intent(out) :: state_vector(:)
type(time_type),  intent(out) :: model_time

integer :: ivar, i
character(len=NF90_MAX_NAME) :: varname

if ( .not. module_initialized ) call static_init_model

state_vector = MISSING_R8

! this is going to have to loop over all the blocks, both to get
! the data values and to get the full grid spacings.

model_time = get_state_time(dirname)

if (do_output()) &
    call print_time(model_time,'time in analysis file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in analysis file '//trim(dirname)//'/header.rst')

! sort the required fields into the order they exist in the
! binary analysis files and fill in the state vector as you
! read each field.  when this routine returns all the data has
! been read.

call get_data(dirname, state_vector)

end subroutine analysis_file_to_statevector



subroutine statevector_to_analysis_file(state_vector, dirname, statedate)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a model netcdf analysis file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: dirname 
type(time_type),  intent(in) :: statedate


integer :: ivar 
character(len=NF90_MAX_NAME) :: varname
character(len=128) :: dirnameout

if ( .not. module_initialized ) call static_init_model

print *, 'in statevector_to_analysis_file, debug, nfields = ', debug, nfields

! sort the required fields into the order they exist in the
! binary analysis files and write out the state vector data
! field by field.  when this routine returns all the data has
! been written.

dirnameout = trim(dirname) // '.out'
call put_data(dirname, dirnameout, state_vector)

! FIXME:
! write out model_time to a text file?
   call print_time(model_time)

 
if (do_output()) &
    call print_time(model_time,'time in analysis file '//trim(dirname)//'/header.rst')
if (do_output()) &
    call print_date(model_time,'date in analysis file '//trim(dirname)//'/header.rst')

end subroutine statevector_to_analysis_file



function get_base_time_ncid( ncid )
!------------------------------------------------------------------
! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_ncid

integer, intent(in) :: ncid

integer :: year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')

get_base_time_ncid = set_date(year, month, day, hour, minute, second)

end function get_base_time_ncid



function get_base_time_fname(filename)
!------------------------------------------------------------------
! The analysis netcdf files have the start time of the experiment.
! The time array contains the time trajectory since then.
! This routine returns the start time of the experiment.

type(time_type) :: get_base_time_fname

character(len=*), intent(in) :: filename

integer :: ncid, year, month, day, hour, minute, second

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR,'get_base_time',string1,source,revision,revdate)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'get_base_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'YEAR'  , year), &
                  'get_base_time', 'get_att year')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MONTH' , month), &
                  'get_base_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'DAY'   , day), &
                  'get_base_time', 'get_att day')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'HOUR'  , hour), &
                  'get_base_time', 'get_att hour')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'MINUTE', minute), &
                  'get_base_time', 'get_att minute')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'SECOND', second), &
                  'get_base_time', 'get_att second')
call nc_check(nf90_close(ncid), 'get_base_time', 'close '//trim(filename))

get_base_time_fname = set_date(year, month, day, hour, minute, second)

end function get_base_time_fname



function get_state_time( dirname )
!------------------------------------------------------------------
! the static_init_model ensures that the model namelists are read.
!
type(time_type)              :: get_state_time
character(len=*), intent(in) :: dirname

type(time_type) :: model_offset, base_time

integer  :: iunit, i, ios
integer  :: istep
real(r8) :: tsimulation
integer  :: iyear, imonth, iday, ihour, imin, isec
integer  :: ndays,nsec

character(len=256) :: filename
character(len=100) :: cLine

if ( .not. module_initialized ) call static_init_model

tsimulation = MISSING_R8
istep       = -1
iyear       = -1
imonth      = -1
iday        = -1
ihour       = -1
imin        = -1
isec        = -1

write(filename,'(a,''/header.rst'')') trim(dirname)

iunit = open_file(trim(filename), action='read')

FILEREAD : do i = 1, 100

   read(iunit,'(a)',iostat=ios) cLine

   if (ios < 0) exit FILEREAD  ! end of file

   if (ios /= 0) then
      write(string1,*) 'cannot read ',trim(filename)
      call error_handler(E_ERR,'get_grid_info',string1,source,revision,revdate)
   endif

   select case( cLine(1:6) ) 
      case('#ISTEP')
         read(iunit,*)istep
      case('#TSIMU')
         read(iunit,*)tsimulation
      case('#TIMES')
         read(iunit,*)iyear
         read(iunit,*)imonth
         read(iunit,*)iday
         read(iunit,*)ihour
         read(iunit,*)imin
         read(iunit,*)isec
      case default
   end select

enddo FILEREAD

call close_file(iunit)

base_time      = set_date(iyear, imonth, iday, ihour, imin, isec)
ndays          = tsimulation/86400
nsec           = tsimulation - ndays*86400
model_offset   = set_time(nsec,ndays)
get_state_time = base_time + model_offset

if (debug > 8) then
   write(*,*)'get_state_time : iyear       ',iyear
   write(*,*)'get_state_time : imonth      ',imonth
   write(*,*)'get_state_time : iday        ',iday
   write(*,*)'get_state_time : ihour       ',ihour
   write(*,*)'get_state_time : imin        ',imin
   write(*,*)'get_state_time : isec        ',isec
   write(*,*)'get_state_time : tsimulation ',tsimulation
   write(*,*)'get_state_time : ndays       ',ndays
   write(*,*)'get_state_time : nsec        ',nsec

   call print_date(     base_time, 'get_state_time:model base date')
   call print_time(     base_time, 'get_state_time:model base time')
   call print_time(  model_offset, 'get_state_time:model offset')
   call print_date(get_state_time, 'get_state_time:model date')
   call print_time(get_state_time, 'get_state_time:model time')
endif

end function get_state_time



!==================================================================
! The (model-specific) private interfaces come last
!==================================================================


subroutine get_grid_dims(nCells, nVertices, nVertLevels, nVertLevelsP1, vertexDegree, nSoilLevels)
!------------------------------------------------------------------
!
! Read the grid dimensions from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: nCells         ! Number of cells
integer, intent(out) :: nVertices      ! Number of vertices
integer, intent(out) :: nVertLevels    ! Number of levels
integer, intent(out) :: nVertLevelsP1  ! Number of levels
integer, intent(out) :: vertexDegree   ! Number of edges that touch a vertex
integer, intent(out) :: nSoilLevels    ! Number of soil layers

integer :: grid_id, dimid

if ( .not. module_initialized ) call static_init_model

! get the ball rolling ...

call nc_check( nf90_open(trim(model_analysis_filename), NF90_NOWRITE, grid_id), &
              'get_grid_dims', 'open '//trim(model_analysis_filename))

! nCells : get dimid for 'nCells' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nCells', dimid), &
              'get_grid_dims','inq_dimid nCells '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nCells), &
            'get_grid_dims','inquire_dimension nCells '//trim(model_analysis_filename))

! nVertices : get dimid for 'nVertices' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertices', dimid), &
              'get_grid_dims','inq_dimid nVertices '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertices), &
            'get_grid_dims','inquire_dimension nVertices '//trim(model_analysis_filename))

! nVertLevels : get dimid for 'nVertLevels' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertLevels', dimid), &
              'get_grid_dims','inq_dimid nVertLevels '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevels), &
            'get_grid_dims','inquire_dimension nVertLevels '//trim(model_analysis_filename))

! nVertLevelsP1 : get dimid for 'nVertLevelsP1' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nVertLevelsP1', dimid), &
              'get_grid_dims','inq_dimid nVertLevelsP1 '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nVertLevelsP1), &
            'get_grid_dims','inquire_dimension nVertLevelsP1 '//trim(model_analysis_filename))

! vertexDegree : get dimid for 'vertexDegree' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'vertexDegree', dimid), &
              'get_grid_dims','inq_dimid vertexDegree '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=vertexDegree), &
            'get_grid_dims','inquire_dimension vertexDegree '//trim(model_analysis_filename))

! nSoilLevels : get dimid for 'nSoilLevels' and then get value

call nc_check(nf90_inq_dimid(grid_id, 'nSoilLevels', dimid), &
              'get_grid_dims','inq_dimid nSoilLevels '//trim(model_analysis_filename))
call nc_check(nf90_inquire_dimension(grid_id, dimid, len=nSoilLevels), &
            'get_grid_dims','inquire_dimension nSoilLevels '//trim(model_analysis_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_grid_dims','close '//trim(model_analysis_filename) )

end subroutine get_grid_dims



subroutine get_grid(nCells, nVertices, nVertLevels, vertexDegree, &
                     latCell, lonCell, zgrid, CellsOnVertex) 
!------------------------------------------------------------------
!
! Read the grid dimensions from the MPAS netcdf file.
!
! The file name comes from module storage ... namelist.

integer, intent(in) :: nCells          ! Number of cells
integer, intent(in) :: nVertices       ! Number of vertices
integer, intent(in) :: nVertLevels     ! Number of levels
integer, intent(in) :: vertexDegree    ! Degree of vertices - No. of edges that touch the vertex

real(r8), dimension(:),   intent(out) :: latCell, lonCell  
real(r8), dimension(:,:), intent(out) :: zgrid          ! geometric height - FIXME: do we need it here?
integer,  dimension(:,:), intent(out) :: CellsOnVertex

integer  :: ncid, VarID

! Read the netcdf file data

call nc_check(nf90_open(trim(model_analysis_filename), nf90_nowrite, ncid), 'get_grid', 'open '//trim(model_analysis_filename))

! Read the variables

call nc_check(nf90_inq_varid(ncid, 'latCell', VarID), &
      'get_grid', 'inq_varid latCell '//trim(model_analysis_filename))
call nc_check(nf90_get_var( ncid, VarID, latCell), &
      'get_grid', 'get_var latCell '//trim(model_analysis_filename))

call nc_check(nf90_inq_varid(ncid, 'lonCell', VarID), &
      'get_grid', 'inq_varid lonCell '//trim(model_analysis_filename))
call nc_check(nf90_get_var( ncid, VarID, lonCell), &
      'get_grid', 'get_var lonCell '//trim(model_analysis_filename))

call nc_check(nf90_inq_varid(ncid, 'zgrid', VarID), &
      'get_grid', 'inq_varid zgrid '//trim(model_analysis_filename))
call nc_check(nf90_get_var( ncid, VarID, zgrid), &
      'get_grid', 'get_var zgrid '//trim(model_analysis_filename))

call nc_check(nf90_inq_varid(ncid, 'CellsOnVertex', VarID), &
      'get_grid', 'inq_varid CellsOnVertex '//trim(model_analysis_filename))
call nc_check(nf90_get_var( ncid, VarID, CellsOnVertex), &
      'get_grid', 'get_var CellsOnVertex '//trim(model_analysis_filename))

call nc_check(nf90_close(ncid), 'get_grid','close '//trim(model_analysis_filename) )

! MPAS analysis files are in radians - at this point DART needs degrees.

latCell = latCell * rad2deg
lonCell = lonCell * rad2deg

! A little sanity check

if ( debug > 1 ) then

   write(*,*)'latCell       range ',minval(latCell),      maxval(latCell)
   write(*,*)'lonCell       range ',minval(lonCell),      maxval(lonCell)
   write(*,*)'zgrid         range ',minval(zgrid  ),      maxval(zgrid  )
   write(*,*)'CellsOnVertex range ',minval(CellsOnVertex),maxval(CellsOnVertex)

endif

end subroutine get_grid



subroutine vector_to_1d_prog_var(x, progvar, data_1d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 1d array.
!
real(r8), dimension(:),   intent(in)  :: x
type(progvartype),        intent(in)  :: progvar
real(r8), dimension(:),   intent(out) :: data_1d_array

integer :: i,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do i = 1, progvar%dimlens(1)
   data_1d_array(i) = x(ii)
   ii = ii + 1
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_1d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_1d_prog_var



subroutine vector_to_2d_prog_var(x, progvar, data_2d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 2d array.
!
real(r8), dimension(:),   intent(in)  :: x
type(progvartype),        intent(in)  :: progvar
real(r8), dimension(:,:), intent(out) :: data_2d_array

integer :: i,j,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_2d_array(i,j) = x(ii)
   ii = ii + 1
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_2d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_2d_prog_var



subroutine vector_to_3d_prog_var(x, progvar, data_3d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 3d array.
!
real(r8), dimension(:),     intent(in)  :: x
type(progvartype),          intent(in)  :: progvar
real(r8), dimension(:,:,:), intent(out) :: data_3d_array

integer :: i,j,k,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do k = 1,progvar%dimlens(3)
do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_3d_array(i,j,k) = x(ii)
   ii = ii + 1
enddo
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_3d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_3d_prog_var



subroutine vector_to_4d_prog_var(x, progvar, data_4d_array)
!------------------------------------------------------------------
! convert the values from a 1d array, starting at an offset,
! into a 4d array.
!
real(r8), dimension(:),       intent(in)  :: x
type(progvartype),            intent(in)  :: progvar
real(r8), dimension(:,:,:,:), intent(out) :: data_4d_array

integer :: i,j,k,m,ii

if ( .not. module_initialized ) call static_init_model

ii = progvar%index1

do m = 1,progvar%dimlens(4)
do k = 1,progvar%dimlens(3)
do j = 1,progvar%dimlens(2)
do i = 1,progvar%dimlens(1)
   data_4d_array(i,j,k,m) = x(ii)
   ii = ii + 1
enddo
enddo
enddo
enddo

ii = ii - 1
if ( ii /= progvar%indexN ) then
   write(string1, *)'Variable '//trim(progvar%varname)//' filled wrong.'
   write(string2, *)'Should have ended at ',progvar%indexN,' actually ended at ',ii
   call error_handler(E_ERR,'vector_to_4d_prog_var', string1, &
                    source, revision, revdate, text2=string2)
endif

end subroutine vector_to_4d_prog_var




function set_model_time_step()
!------------------------------------------------------------------
! the static_init_model ensures that the model namelists are read.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call static_init_model

! FIXME - determine when we can stop the model

   set_model_time_step = set_time(0, 1) ! (seconds, days)

end function set_model_time_step



subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )
!------------------------------------------------------------------

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, varid
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)
ncols = size(table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      string1 = 'model_vars_nml:model state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in model analysis variable list

!%!   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
!%!   call nc_check(NF90_inq_varid(ncid, trim(varname), varid), &
!%!                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector 

   if ( debug > 0 ) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), ' ', trim(table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables



function FindTimeDimension(ncid) result(timedimid)
!------------------------------------------------------------------

! Find the Time Dimension ID in a netCDF file.
! If there is none - (spelled the obvious way) - the routine
! returns a negative number. You don't HAVE to have a TIME dimension.

integer                      :: timedimid
integer,          intent(in) :: ncid

integer :: nc_rc

TimeDimID = -1 ! same as the netCDF library routines. 
nc_rc = nf90_inq_dimid(ncid,'TIME',dimid=TimeDimID)
if ( nc_rc /= NF90_NOERR ) then ! did not find it - try another spelling
   nc_rc = nf90_inq_dimid(ncid,'Time',dimid=TimeDimID)
   if ( nc_rc /= NF90_NOERR ) then ! did not find it - try another spelling
      nc_rc = nf90_inq_dimid(ncid,'time',dimid=TimeDimID)
   endif
endif

end function FindTimeDimension


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
