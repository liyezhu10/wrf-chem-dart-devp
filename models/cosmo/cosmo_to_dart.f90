! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program cosmo_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between ncommas and DART
!
! method: Read ncommas "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The cosmo filename is read from the cosmo_in namelist
!         <edit cosmo_to_dart_output_file in input.nml:cosmo_to_dart_nml>
!         cosmo_to_dart
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size,get_state_vector,get_state_time
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date,&
                             set_calendar_type,JULIAN

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: cosmo_to_dart_output_file  = 'dart.ud'

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
real(r8),allocatable  :: x(:)
type(time_type)       :: model_time
character(len=256)    :: cosmo_restart_filename

!======================================================================

call initialize_utilities(progname='cosmo_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

model_time = get_state_time()
x_size = get_model_size()

print*,x_size

allocate(x(1:x_size))
x(:)=get_state_vector()

iunit = open_restart_write(cosmo_to_dart_output_file)

call awrite_state_restart(model_time, x, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call set_calendar_type(JULIAN)

call print_date(model_time, str='cosmo_to_dart:cosmo model date')
call print_time(model_time, str='cosmo_to_dart:DART model time')
call timestamp(string1=source, pos='end')

end program cosmo_to_dart
