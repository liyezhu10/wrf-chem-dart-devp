! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program paste_inflation_restarts

!----------------------------------------------------------------------
! purpose: interface between DART and the CESM model
!
! method: Read DART state vector and overwrite values in a CESM restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called cesm_in.DART is created with a time_manager_nml namelist
!         appropriate to advance CESM to the requested time.
!
!         The paste_inflation_restarts_nml namelist setting for advance_time_present
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             nmlfileunit, error_handler, E_MSG, do_nml_file, do_nml_term
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart, &
                             open_restart_write, awrite_state_restart
use time_manager_mod, only : time_type, print_time, print_date
use        model_mod, only : static_init_model, get_model_size, set_start_end

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=265) :: cam_inflation_restart = ''
character(len=265) :: pop_inflation_restart = ''
real(r8)           :: cam_damping = 1.0_r8
real(r8)           :: pop_damping = 1.0_r8
character(len=265) :: combined_inflation_output = ''

namelist /paste_inflation_restarts_nml/ &
   cam_inflation_restart, &
   pop_inflation_restart, &
   cam_damping,           &
   pop_damping,           &
   combined_inflation_output

!----------------------------------------------------------------------

integer               :: iunit1, iunit2, x_start, x_end, x_size, io
type(time_type)       :: model_time1, model_time2
real(r8), allocatable :: inf_mean(:), inf_sd(:)
character(len=512)    :: msgstring

!----------------------------------------------------------------------

call initialize_utilities(progname='paste_inflation_restarts')


!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the CESM namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

! Read the namelist to get the input filename.

call find_namelist_in_file("input.nml", "paste_inflation_restarts_nml", iunit1)
read(iunit1, nml = paste_inflation_restarts_nml, iostat = io)
call check_namelist_read(iunit1, io, "paste_inflation_restarts_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=paste_inflation_restarts_nml)
if (do_nml_term()) write(     *     , nml=paste_inflation_restarts_nml)


write(msgstring, '(A)') 'converting CAM and POP inflation files into a single CESM inflation file'
call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)

write(msgstring, '(A)') 'reading CAM file '//trim(cam_inflation_restart)
call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)

write(msgstring, '(A)') 'reading POP file '//trim(pop_inflation_restart)
call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)

if (cam_damping /= 1.0_r8) then
   write(msgstring, '(A,F9.3)') 'applying damping to CAM inflation of ', cam_damping
   call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)
endif

if (pop_damping /= 1.0_r8) then
   write(msgstring, '(A,F9.3)') 'applying damping to POP inflation of ', pop_damping
   call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)
endif

write(msgstring, '(A)') 'writing CESM file '//trim(combined_inflation_output)
call error_handler(E_MSG, 'paste_inflation_restarts: ', msgstring)


! combined state vector size (computed by static_init_model() above)

x_size = get_model_size()
allocate(inf_mean(x_size), inf_sd(x_size))

! open both restart files.
iunit1 = open_restart_read(cam_inflation_restart)
iunit2 = open_restart_read(pop_inflation_restart)

! get the inflation means
call set_start_end('CAM', x_start, x_end)
call aread_state_restart(model_time1, inf_mean(x_start:x_end), iunit1)
call set_start_end('POP', x_start, x_end)
call aread_state_restart(model_time2, inf_mean(x_start:x_end), iunit2)

! now the sd copies
call set_start_end('CAM', x_start, x_end)
call aread_state_restart(model_time1, inf_sd(x_start:x_end), iunit1)
call set_start_end('POP', x_start, x_end)
call aread_state_restart(model_time2, inf_sd(x_start:x_end), iunit2)

! now if the damping is /= 1, apply the damping
if (cam_damping /= 1.0_r8 .or. pop_damping /= 1.0_r8) then
   call set_start_end('CAM', x_start, x_end)
   inf_mean(x_start:x_end) = inf_mean(x_start:x_end) * cam_damping
   call set_start_end('POP', x_start, x_end)
   inf_mean(x_start:x_end) = inf_mean(x_start:x_end) * pop_damping
endif

call close_restart(iunit1)
call close_restart(iunit2)

! and now write them out as a single vector
iunit1 = open_restart_write(combined_inflation_output)

call awrite_state_restart(model_time1, inf_mean, iunit1)
call awrite_state_restart(model_time1, inf_sd,   iunit1)

call close_restart(iunit1)


end program paste_inflation_restarts

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
