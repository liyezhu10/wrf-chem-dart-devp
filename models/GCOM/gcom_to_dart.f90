! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program gcom_to_dart

!----------------------------------------------------------------------
! purpose: interface between GCOM and DART
!
! method: Read GCOM "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The GCOM filename is read from the gcom_in namelist
!         <edit gcom_to_dart_output_file in input.nml:gcom_to_dart_nml>
!         gcom_to_dart
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : restart_file_to_sv, static_init_model, &
                             get_model_size, get_gcom_restart_filename
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

use netcdf
implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: gcom_to_dart_output_file  = 'dart_ics'

namelist /gcom_to_dart_nml/ gcom_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character (len = 128) :: gcom_restart_filename = 'no_gcom_restart_filename' 

!----------------------------------------------------------------------

call initialize_utilities(progname='gcom_to_dart')

!----------------------------------------------------------------------
! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, starting date, deltaT, etc.
!----------------------------------------------------------------------

call static_init_model()

!----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "gcom_to_dart_nml", iunit)
read(iunit, nml = gcom_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "gcom_to_dart_nml") ! closes, too.

call get_gcom_restart_filename( gcom_restart_filename )

write(*,*)
write(*,'(''gcom_to_dart:converting GCOM restart file '',A, &
      &'' to DART file '',A)') &
       trim(gcom_restart_filename), trim(gcom_to_dart_output_file)

!----------------------------------------------------------------------
! Now that we know the names, get to work.
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))
call restart_file_to_sv(gcom_restart_filename, statevector, model_time) 

iunit = open_restart_write(gcom_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='gcom_to_dart:GCOM  model date')
call print_time(model_time, str='gcom_to_dart:DART model time')
call finalize_utilities('gcom_to_dart')

end program gcom_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
