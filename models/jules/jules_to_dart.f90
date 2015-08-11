! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program jules_to_dart

!----------------------------------------------------------------------
! purpose: interface between jules and DART
!
! method: Read jules "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The jules filename is read from the jules_in namelist
!         <edit jules_to_dart_output_file in input.nml:jules_to_dart_nml>
!         jules_to_dart
!
! author: Tim Hoar 10 August 2015
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, jules_to_dart_state_vector
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=512) :: jules_to_dart_output_file  = 'dart_ics'

namelist /jules_to_dart_nml/ jules_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)

!======================================================================

call initialize_utilities(progname='jules_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "jules_to_dart_nml", iunit)
read(iunit, nml = jules_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "jules_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

! Each variable specifies its own file of origin.
call jules_to_dart_state_vector(statevector, model_time) 

iunit = open_restart_write(jules_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='jules_to_dart: model date')
call print_time(model_time, str='jules_to_dart: model time')

call finalize_utilities('jules_to_dart')

end program jules_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
