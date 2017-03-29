! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This module is used to wrap around the basic portions of existing dynamical models to
!> add capabilities needed by the standard assimilation methods.
!>
!>@todo FIXME explore redundant routines , especially 'diagnostic' ones

module assim_model_mod

use    types_mod, only : r8

use location_mod, only : location_type, &
                         get_close_obs_init   => get_close_init,     &
                         get_close_state_init => get_close_init,     &
                         get_close_obs_destroy   => get_close_destroy, &
                         get_close_state_destroy => get_close_destroy

use time_manager_mod, only : time_type,                                            &
                             operator(<), operator(>), operator(+), operator(-),   &
                             operator(/), operator(*), operator(==), operator(/=)

use utilities_mod, only : register_module

use ensemble_manager_mod, only : ensemble_type

! these are the required model_mod interfaces:
use     model_mod, only : get_model_size, static_init_model, get_state_meta_data,  &
                          get_model_time_step => shortest_time_between_assimilations, &
                          shortest_time_between_assimilations,                     &
                          init_conditions, init_time, adv_1step, end_model,        &
                          pert_model_copies, get_close_obs, get_close_state,       &
                          convert_vertical_obs, convert_vertical_state,            &
                          interpolate => model_interpolate,                        &
                          read_model_time, write_model_time

implicit none
private

public :: static_init_assim_model, &
          end_assim_model, &
          get_initial_condition, &
          get_closest_state_time_to, &
          get_state_meta_data, &
          get_model_time_step, &
          get_model_size,  &
          adv_1step, &
          interpolate, &
          pert_model_copies, &
          get_close_obs_init, &
          get_close_obs, &
          get_close_obs_destroy, &
          get_close_state_init, &
          get_close_state, &
          get_close_state_destroy, &
          convert_vertical_obs, &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


!>@todo FIXME:
! DEPRECATE THIS. i believe this is unused
! in the current system.

!! Type to keep model state and time together
!type assim_model_type
!   private
!   real(r8), pointer :: state_vector(:)
!   type(time_type)   :: time
!   integer           :: model_size
!   integer           :: copyID
!end type assim_model_type

! Permanent class storage for model_size
integer :: model_size

! Ensure init code is called exactly once
logical :: module_initialized = .false.

!-------------------------------------------------------------

contains

!======================================================================


  subroutine static_init_assim_model()
!----------------------------------------------------------------------
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
! is initializing the position of the state variables as location types.

! only execute this code once, even if called multiple times.
if (module_initialized) return

! First thing to do is echo info to logfile ... 
call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine static_init_assim_model


function get_closest_state_time_to(model_time, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

type(time_type) :: get_closest_state_time_to
type(time_type), intent(in) :: model_time, time

type(time_type) :: time_step

! Get the model time step capabilities - this is the
! shortest amount of time you can/want to ask the model
! to advance the state
time_step = shortest_time_between_assimilations()

if(model_time > time) then
   ! If model_time is past start of obs window, don't advance it
   get_closest_state_time_to = model_time
   return
endif

get_closest_state_time_to = model_time

do while((time_step + 2*get_closest_state_time_to) < 2*time)
   get_closest_state_time_to = get_closest_state_time_to + time_step
enddo

end function get_closest_state_time_to


subroutine get_initial_condition(time, x)
!----------------------------------------------------------------------
!

type(time_type), intent(out) :: time
real(r8),        intent(out) :: x(:)

call init_conditions(x)

call init_time(time)

end subroutine get_initial_condition


!-------------------------------------------------------------------

subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model. For now, only thing to do is tell model to end.

call end_model()

end subroutine end_assim_model

!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
