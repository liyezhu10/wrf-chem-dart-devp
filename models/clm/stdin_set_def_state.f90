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

use           types_mod, only : r8, i8

use        location_mod, only : location_type, get_location

use       utilities_mod, only : open_file, close_file, &
                                E_ERR, error_handler

use        obs_kind_mod, only : QTY_SURFACE_PRESSURE, &
                                QTY_LANDMASK

use state_structure_mod, only : get_num_domains, &
                                get_num_variables, &
                                get_variable_name, &
                                get_index_start, &
                                get_index_end

use           model_mod, only : static_init_model,  &
                                get_model_size,     &
                                get_state_meta_data

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

integer :: var_type, iunit, idomain, var_domain, ivar, var_id
integer(i8) :: i, model_size, indx1, indxN
type(location_type) :: location
real(r8) :: loc3d(3)
real(r8) :: lon, lat, vert
character(len=20) :: varname

! Write to file
iunit = open_file('stdin_set_def_state.out',form='formatted',action='write')

! Get the model size
call static_init_model()
model_size = get_model_size()

! Figure out where the variable of interest starts/stop in the DART state.
var_domain = -1
DOMAIN: do idomain = 1,get_num_domains()
   VARIABLES : do ivar = 1, get_num_variables(idomain)
      varname = get_variable_name(idomain,ivar)
      if (varname == 'FSIF') then
         var_domain = idomain
         var_id     = ivar
         exit DOMAIN
      endif
   enddo VARIABLES 
enddo DOMAIN

if (var_domain < 0) call error_handler(E_ERR,'bob','woe',source,revision,revdate)

indx1 = get_index_start(var_domain, var_id)
indxN = get_index_end(  var_domain, var_id)
write(*,*)'domain ', var_domain, var_id, indx1, indxN

! Set the maximum number of intended observations
write(iunit, *) model_size ! Input upper bound on number of observations in sequence

! No values or qc
write(iunit, *) 0 ! Input number of copies of data (0 for just a definition)
write(iunit, *) 0 ! Input number of quality control values per field (0 or greater)

! Loop through all the state variables, set the obs variance accordingly
MODEL : do i = indx1, indxN

   call get_state_meta_data(i, location, var_type)

   ! skip locations masked as water
   ! We don't have a QTY_WATERMASK, so I am abusing QTY_LANDMASK
   if (var_type == QTY_LANDMASK) cycle MODEL

   loc3d = get_location(location)
   lon  = loc3d(1)
   lat  = loc3d(2)
   vert = loc3d(3)

   ! There are more obs
   write(iunit, *) 0  ! input a -1 if there are no more obs

   ! tricky bit ...  what kind of observation
   write(iunit, *) 31   ! 31 OCO2_SIF  - maybe !!!

                        ! Vertical coordinate options
   write(iunit, *) 3    ! 3  --> height
   write(iunit, *) 0.0  !  Vertical coordinate height (in meters)
   write(iunit, *) lon  ! Input longitude: value 0 to 360.0
   write(iunit, *) lat  ! Input latitude: value -90.0 to 90.0

   write(iunit, *) 2000, 1, 1, 0, 0, 0 ! Time is 0 days 0 seconds for create obs sequence

   write(iunit, *) 0.25 ! Input the error variance for this observation definition

enddo MODEL

! There are more obs
write(iunit, *) -1  ! input a -1 if there are no more obs

! Output the default set_def.out file name
write(iunit, '(''set_def.out'')')

call close_file(iunit)

end program stdin_set_def_state

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
