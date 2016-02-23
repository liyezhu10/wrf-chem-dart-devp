! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program openggcm_to_netcdf

!----------------------------------------------------------------------
! purpose: interface between openggcm and DART
!
! method: Read openggcm "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The openggcm filename is read from the openggcm_in namelist
!         <edit openggcm_to_netcdf_output_file in input.nml:openggcm_to_netcdf_nml>
!         openggcm_to_netcdf
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r4, r8, digits12
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, nc_check
use time_manager_mod, only : time_type, print_time, print_date, set_date, &
                             set_calendar_type, operator(-), get_time

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

character(len=256) :: openggcm_to_netcdf_input_file  = '../data/da0002.dart.pot'
character(len=256) :: openggcm_to_netcdf_output_file = 'openggcm.nc'

namelist /openggcm_to_netcdf_nml/ openggcm_to_netcdf_output_file, &
                                  openggcm_to_netcdf_input_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit
type(time_type)       :: model_time, model_time_base, model_time_offset
real(r4), allocatable :: statevector(:,:)
real(r8), allocatable :: latitude(:)
real(r8), allocatable :: longitude(:)
real(digits12) :: time_offset_seconds

integer :: iyear, imonth, iday, ihour, iminute, isecond
integer :: ilat, ilon
integer :: nthe, nphi
real(r8) :: dlat, dlon

integer :: ncid, NtheDimID, NphiDimID, VarID, LonVarID, LatVarID
integer :: TimeVarID

!----------------------------------------------------------------------

call initialize_utilities(progname='openggcm_to_netcdf')

!----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "openggcm_to_netcdf_nml", iunit)
read(iunit, nml = openggcm_to_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "openggcm_to_netcdf_nml") ! closes, too.

write(*,*)
write(*,'(''openggcm_to_netcdf:converting openggcm restart file '',A, &
      &'' to DART file '',A)') &
       trim(openggcm_to_netcdf_input_file), trim(openggcm_to_netcdf_output_file)

iunit = open_file(openggcm_to_netcdf_input_file, form='unformatted',action='read')
read(iunit) nphi, nthe
read(iunit) iyear, imonth, iday, ihour, iminute, isecond
allocate(statevector(nphi,nthe))
allocate(latitude(nthe), longitude(nphi))
read(iunit) statevector
call close_file(iunit)

iyear = iyear + 1900

dlat = 180.0_r8/real(nthe-1,r8)
do ilat = 1,nthe  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   latitude(ilat) = dlat*real(ilat-1,r8)
enddo

dlon = 360.0_r8/real(nphi-1,r8)
do ilon = 1,nphi  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   longitude(ilon) = dlon*real(ilon-1,r8)
enddo

write(*,*)'TJH iyear, imonth, iday, ihour, iminute, isecond', iyear, imonth, iday, ihour, iminute, isecond
write(*,*)'TJH nphi, nthe is ',nphi, nthe
write(*,*)'file size should be ',4 + 4+4 + 4 + 4 + nphi*nthe*4 + 4

call set_calendar_type('Gregorian')
model_time      = set_date(iyear, imonth, iday, ihour, iminute, isecond)
model_time_base = set_date(1900, 1, 1, 0, 0, 0)
model_time_offset = model_time - model_time_base
call get_time(model_time_offset, isecond, iday)

time_offset_seconds = real(iday,digits12) + real(isecond,digits12)/86400.0_digits12

call print_date(model_time, str='openggcm_to_netcdf:openggcm  model date')
call print_time(model_time, str='openggcm_to_netcdf:DART model time')

call nc_check(nf90_create(openggcm_to_netcdf_output_file, NF90_CLOBBER, ncid),'openggcm_to_netcdf')

! define dimensions

call nc_check(nf90_def_dim(ncid, 'nphi', nphi, NphiDimID),'openggcm_to_netcdf')
call nc_check(nf90_def_dim(ncid, 'nthe', nthe, NtheDimID),'openggcm_to_netcdf')

call nc_check(nf90_def_var(ncid, 'time', nf90_double, TimeVarID),'openggcm_to_netcdf')
call nc_check(nf90_put_att(ncid, TimeVarID, 'units', 'days since 1900-01-01'), &
                  'put_att time:units: openggcm_to_netcdf')

call nc_check(nf90_def_var(ncid, 'nphi', nf90_real, (/ NphiDimID /), LonVarID),'openggcm_to_netcdf')
call nc_check(nf90_put_att(ncid, LonVarID, 'short_name', 'geomagnetic longitude'), &
                  'put_att nphi:short_name: openggcm_to_netcdf')

call nc_check(nf90_def_var(ncid, 'nthe', nf90_real, (/ NtheDimID /), LatVarID),'openggcm_to_netcdf')
call nc_check(nf90_put_att(ncid, LatVarID, 'short_name', 'geomagnetic latitude'), &
                  'put_att nthe:short_name: openggcm_to_netcdf')

call nc_check(nf90_def_var(ncid, 'pot' , nf90_real, (/ NphiDimID, NtheDimID /), VarID),'openggcm_to_netcdf')

call nc_check(nf90_put_att(ncid, VarID, 'short_name', 'electric potential'), &
                  'put_att pot:short_name: openggcm_to_netcdf')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'ionosphere electric potential'), &
                  'put_att pot:short_name: openggcm_to_netcdf')
call nc_check(nf90_put_att(ncid, VarID, 'units', 'volts'), &
                  'put_att pot:short_name: openggcm_to_netcdf')

call nc_check(nf90_enddef(ncid),'openggcm_to_netcdf')

call nc_check(nf90_put_var(ncid,  LonVarID,           longitude),'openggcm_to_netcdf')
call nc_check(nf90_put_var(ncid,  LatVarID,            latitude),'openggcm_to_netcdf')
call nc_check(nf90_put_var(ncid,     VarID,         statevector),'openggcm_to_netcdf')
call nc_check(nf90_put_var(ncid, TimeVarID, time_offset_seconds),'openggcm_to_netcdf')

call nc_check(nf90_close(ncid),'openggcm_to_netcdf')

call finalize_utilities('openggcm_to_netcdf')

end program openggcm_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
