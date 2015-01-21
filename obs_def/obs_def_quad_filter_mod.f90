! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! An example of a simple forward operator that involves more than
! just interpolating directly from a state vector in a model.
!
! This section defines a specific type in the left column and
! can be any string you want to use for an observation.  The
! right column must be a generic kind that already exists in
! the obs_kind/DEFAULT_obs_kind_mod.F90 file.

! BEGIN DART PREPROCESS KIND LIST
! QUAD_FILTER_SQUARED_ERROR,    KIND_QUAD_FILTER_SQUARED_ERROR
! END DART PREPROCESS KIND LIST

! This section will be added to the main obs_def_mod.f90 that
! is going to be generated, to allow it to call the code we
! are defining here.

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_quad_filter_mod, only : get_expected_quad_sq_err
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! This section will be dropped into a large case statement in the
! main obs_def_mod.f90 code to control what happens with each
! observation type that is processed.

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!   case(QUAD_FILTER_SQUARED_ERROR)
!        call get_expected_quad_sq_err(state, location, obs_val, istatus) 
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! The next few sections do nothing because there is no additional
! data to read, write, or prompt for.  But there still needs to be a
! case statement in the large select, so they must be here.

! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(QUAD_FILTER_SQUARED_ERROR)
!     continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(QUAD_FILTER_SQUARED_ERROR)
!     continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(QUAD_FILTER_SQUARED_ERROR)
!     continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! This is the code that implements the forward operator.
! Define a module, and make public anything that will be called
! from the main obs_def_mod.f90 file.  Here it is just the
! get_expected routine.  There isn't any initialization needed
! but the stub is there; it could read a namelist if there are
! any run-time options to be set.

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_quad_filter_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module
use     location_mod, only : location_type

implicit none
private

public :: get_expected_quad_sq_err

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

contains

! ---------------------------------------------------

subroutine initialize_module

! Handle any module initialization tasks

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

! ---------------------------------------------------

subroutine get_expected_quad_sq_err(state_vector, location, sqerr, istatus)  
 real(r8),            intent(in)  :: state_vector(:)
 type(location_type), intent(in)  :: location
 real(r8),            intent(out) :: sqerr
 integer,             intent(out) :: istatus

! Pseudo forward operator for the paired observations.
! Always returns error and missing r8, intentionally.

if ( .not. module_initialized ) call initialize_module

istatus = 1
sqerr = MISSING_R8

end subroutine get_expected_quad_sq_err

! ---------------------------------------------------

end module obs_def_quad_filter_mod
! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
