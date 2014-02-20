! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program cable_to_dart

!----------------------------------------------------------------------
! purpose: interface between CABLE and DART
!
! method: Read CABLE "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The CABLE filename is read from the cable_in namelist
!         <edit cable_to_dart_output_file in input.nml:cable_to_dart_nml>
!         cable_to_dart
!
! author: Tim Hoar 20 February 2014
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, error_handler, nc_check, file_exist, logfileunit
use        model_mod, only : get_model_size, cable_state_to_dart_vector, &
                             get_cable_restart_filename, static_init_model
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date, set_date, &
                             get_time, operator(-)

use typesizes
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

logical               :: replace_cable_time         = .false.
integer, dimension(6) :: cable_state_yyyymmddhhmmss = 0
character(len=128)    :: cable_to_dart_output_file  = 'dart_ics'

namelist /cable_to_dart_nml/ cable_to_dart_output_file, &
              replace_cable_time, cable_state_yyyymmddhhmmss

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: cable_restart_filename

!======================================================================

call initialize_utilities(progname='cable_to_dart')

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the cable namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "cable_to_dart_nml", iunit)
read(iunit, nml = cable_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "cable_to_dart_nml") ! closes, too.

call get_cable_restart_filename( cable_restart_filename )

write(*,*)
write(*,'(''cable_to_dart:converting cable restart file '',A, &
      &'' to DART file '',A)') &
       trim(cable_restart_filename), trim(cable_to_dart_output_file)

write(logfileunit,*)
write(logfileunit,'(''cable_to_dart:converting cable restart file '',A, &
      &'' to DART file '',A)') &
       trim(cable_restart_filename), trim(cable_to_dart_output_file)

if ( replace_cable_time ) call replace_timestamp( cable_restart_filename )

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call cable_state_to_dart_vector(cable_restart_filename, statevector, model_time) 

iunit = open_restart_write(cable_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call print_date(model_time, str='cable_to_dart:CABLE model date')
call print_time(model_time, str='cable_to_dart:DART  model time')
call print_date(model_time, str='cable_to_dart:CABLE model date',iunit=logfileunit)
call print_time(model_time, str='cable_to_dart:DART  model time',iunit=logfileunit)

call finalize_utilities('cable_to_dart')


contains


subroutine replace_timestamp( filename )

character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer               :: ncid, VarID, ncNdims, dimlen
integer               :: iyear, imonth, iday, ihour, imin, isec
character(len=128)    :: unitstring, string1, string2
type(time_type)       :: time_base, time_desired, time_offset
integer, dimension(1) :: seconds

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*)'Cannot open file ', trim(filename),' for WRITING.'
   call error_handler(E_ERR,'replace_timestamp',string1,source,revision,revdate)
endif

call nc_check(nf90_open(trim(filename), NF90_WRITE, ncid), &
        'replace_timestamp','open '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'time', VarID), &
        'replace_timestamp', 'inq_varid time'//trim(filename))

call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
        'replace_timestamp', 'inquire time'//trim(filename))

if (ncNdims /= 1) then
   write(string1,*)'time variable has unexpected rank.'
   write(string2,*)'expected 1, got ',ncNdims
   call error_handler(E_ERR, &
        'replace_timestamp', string1, source, revision, revdate, text2=string2)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
        'replace_timestamp', trim(string1))

if (dimlen /= 1) then
   write(string1,*)'time variable has unexpected length.'
   write(string2,*)'expected 1, got ',dimlen
   call error_handler(E_ERR, &
        'replace_timestamp', string1, source, revision, revdate, text2=string2)
endif

! FIXME check for reasonable input values - although set_date does that.

iyear  = cable_state_yyyymmddhhmmss(1)
imonth = cable_state_yyyymmddhhmmss(2)
iday   = cable_state_yyyymmddhhmmss(3)
ihour  = cable_state_yyyymmddhhmmss(4)
imin   = cable_state_yyyymmddhhmmss(5)
isec   = cable_state_yyyymmddhhmmss(6)

write(unitstring,'(''seconds since '',i4.4,''-01-01 00:00:00'')') iyear

! replace units string
call nc_check(nf90_Redef(ncid), 'replace_timestamp','redef '//trim(filename))
call nc_check(nf90_put_att(ncid, VarID,'units',trim(unitstring)),&
        'replace_timestamp','put_att units '//trim(filename))
call nc_check(nf90_put_att(ncid, VarID,'DART','redefined by DART'),&
        'replace_timestamp', 'put att DART '//trim(filename))

call nc_check(nf90_enddef(ncid), 'replace_timestamp', 'enddef')

! calculate and replace value

time_base    = set_date(iyear, 1, 1, 0, 0, 0)
time_desired = set_date(iyear, imonth, iday, ihour, imin, isec)
time_offset  = time_desired - time_base

call get_time(time_offset,isec,iday)

seconds(1) = iday * 86400 + isec

call nc_check(nf90_put_var(ncid, VarID, seconds), &
        'replace_timestamp', 'put_var time'//trim(filename))

call nc_check(nf90_close(ncid),'replace_timestamp','close '//trim(filename))

call print_date(time_desired, str='cable_to_dart: OVERRRIDE : NEW CABLE model date')
call print_time(time_desired, str='cable_to_dart: OVERRRIDE : NEW CABLE model time')

end subroutine replace_timestamp

end program cable_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
