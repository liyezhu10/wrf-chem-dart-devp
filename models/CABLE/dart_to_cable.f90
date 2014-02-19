! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cable

!----------------------------------------------------------------------
! purpose: interface between DART and the CABLE model
!
! method: Read DART state vector and overwrite values in a CABLE restart file.
!         If the DART state vector has an 'advance_to_time' present, 
!         it is read ... but nothing happens with it at this time.
!
!         The dart_to_cable_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, dart_vector_to_model_file, &
                             get_model_size, get_cable_restart_filename, &
                             get_time_origin

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_cable_input_file = 'dart_restart'
logical               :: advance_time_present   = .false.

namelist /dart_to_cable_nml/ dart_to_cable_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: cable_restart_filename
integer               :: iunit, io, x_size, newunit
type(time_type)       :: model_time, adv_to_time
type(time_type)       :: cable_origin, temp_time
real(r8), allocatable :: statevector(:)
integer               :: days, seconds, currentseconds, futureseconds

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_cable')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the cable namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_cable_nml", iunit)
read(iunit, nml = dart_to_cable_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cable_nml")

call get_cable_restart_filename( cable_restart_filename )

write(*,*)
write(*,'(''dart_to_cable:converting DART file '',A, &
      &'' to cable restart file '',A)') &
     trim(dart_to_cable_input_file), trim(cable_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_cable_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)

   ! FIXME convey the advance_to_time to CABLE 
   ! There is more to do. Convert the model_time and adv_to_time to
   ! the units of the "time" variable in the restart file (same as the
   ! atmospheric forcing file). Those will be used to construct a call
   ! to ncks to subset the gigantic forcing files in to ones that are
   ! appropriate for this particular model advance.

   cable_origin = get_time_origin()   ! from the restart file, essentially

   ! figure out the current model time 
   temp_time = model_time - cable_origin
   call get_time(temp_time, seconds, days)
   currentseconds = seconds + days*86400

   ! figure out the advance-to-time time 
   temp_time = adv_to_time - cable_origin
   call get_time(temp_time, seconds, days)
   futureseconds = seconds + days*86400

   newunit = open_file('time_control.txt',action='rewind')
   write(newunit,*)currentseconds
   write(newunit,*)futureseconds
   call print_date( model_time,'dart_to_cable:DART   model date',newunit)
   call print_date(adv_to_time,'dart_to_cable:DART desired date',newunit)

else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current cable state vector
!----------------------------------------------------------------------

call dart_vector_to_model_file(statevector, cable_restart_filename, model_time)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_cable:cable model date')
call print_time( model_time,'dart_to_cable:DART  model time')
call print_date( model_time,'dart_to_cable:cable model date',logfileunit)
call print_time( model_time,'dart_to_cable:DART  model time',logfileunit)

if ( advance_time_present ) then
   call print_time(adv_to_time,'dart_to_cable:advance_to time')
   call print_date(adv_to_time,'dart_to_cable:advance_to date')
   call print_time(adv_to_time,'dart_to_cable:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_cable:advance_to date',logfileunit)
endif

call finalize_utilities('dart_to_cable')

end program dart_to_cable

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
