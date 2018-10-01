! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program stdin_set_def_state

!> This creates the text input for create_obs_sequence that will
!> create an observation at every state variable location.
!> The initial version used explicit knowledge of how the model
!> state is constructed so that only SIF observations are created.
!>
!> It may be useful to run model_mod_check with verbose=.true.
!> to determine which indices in the model state are consistent
!> with the observations you want to create.

use     types_mod, only : r8, i8
use  location_mod, only : location_type, get_location
use utilities_mod, only : get_unit
use  obs_kind_mod, only : QTY_SURFACE_PRESSURE
use     model_mod, only : static_init_model, get_model_size, &
                          get_state_meta_data

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: i, model_size, var_type, iunit
type(location_type) :: location
real(r8) :: loc3d(3)
real(r8) :: lon, lat, vert

! Write to file
iunit  = get_unit()
open(unit = iunit, file = 'stdin_set_def_state.out')

! Get the model size
call static_init_model()
model_size = get_model_size()

! Set the number of state variables, all observed
! write(iunit, *) model_size    ! Input upper bound on number of observations in sequence
write(iunit, *) 55296    ! Input upper bound on number of observations in sequence

! No values or qc
write(iunit, *) 0 ! Input number of copies of data (0 for just a definition)
write(iunit, *) 0 ! Input number of quality control values per field (0 or greater)

! Loop through all the state variables, set the obs variance accordingly
!do i = 1, model_size
do i = 199476, model_size

   call get_state_meta_data(int(i,i8), location, var_type)
   loc3d = get_location(location)
   lon  = loc3d(1)
   lat  = loc3d(2)
   vert = loc3d(3)

   ! There are more obs
   write(iunit, *) 0  ! input a -1 if there are no more obs

   ! tricky bit ...  what kind of observation
   write(iunit, *) 31  ! 31 OCO2_SIF  - maybe !!!

                        ! Vertical coordinate options
   write(iunit, *) 3    ! 3  --> height
   write(iunit, *) 0.0  !  Vertical coordinate height (in meters)
   write(iunit, *) lon  ! Input longitude: value 0 to 360.0
   write(iunit, *) lat  ! Input latitude: value -90.0 to 90.0

   write(iunit, *) 2000, 1, 1, 0, 0, 0 ! Time is 0 days 0 seconds for create obs sequence

   write(iunit, *) 0.25 ! Input the error variance for this observation definition

end do

! Output the default set_def.out file name
write(iunit, '(''set_def.out'')')

end program stdin_set_def_state

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
