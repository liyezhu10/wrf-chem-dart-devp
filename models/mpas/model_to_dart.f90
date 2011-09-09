! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between model and DART
!
! method: Read model "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The model dirname is read from the model_in namelist
!         <edit model_to_dart_output_file in input.nml:model_to_dart_nml>
!         model_to_dart
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, restart_file_to_sv, &
                             get_model_restart_dirname
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: model_to_dart_output_file  = 'dart.ud'
character(len=256) :: model_restart_dirname = 'model_restartdir'

namelist /model_to_dart_nml/    &
     model_to_dart_output_file, &
     model_restart_dirname

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)

!======================================================================

call initialize_utilities(progname='model_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output dirname.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "model_to_dart_nml", iunit)
read(iunit, nml = model_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "model_to_dart_nml") ! closes, too.

write(*,*)
write(*,*) 'model_to_dart: converting model restart files in directory ', &
           "'"//trim(model_restart_dirname)//"'" 
write(*,*) ' to DART file ', "'"//trim(model_to_dart_output_file)//"'"

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call get_model_restart_dirname( model_restart_dirname )

call restart_file_to_sv(model_restart_dirname, statevector, model_time) 

iunit = open_restart_write(model_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call print_date(model_time, str='model_to_dart:model model date')
call print_time(model_time, str='model_to_dart:DART model time')
call timestamp(string1=source, pos='end')

end program model_to_dart

