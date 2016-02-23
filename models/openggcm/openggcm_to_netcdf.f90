! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program openggcm_to_netcdf

!-----------------------------------------------------------------------
! purpose: create a netCDF file from an openggcm unformatted binary
!
! USAGE:  The filenames are read from the openggcm_to_netcdf_nml namelist
!
! author: Tim Hoar 2/23/16
!-----------------------------------------------------------------------

use        types_mod, only : r4, r8, digits12
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, nc_check, &
                             error_handler, E_ERR, E_MSG
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
logical            :: verbose = .false.

namelist /openggcm_to_netcdf_nml/ openggcm_to_netcdf_output_file, &
                                  openggcm_to_netcdf_input_file, &
                                  verbose

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

integer               :: io, iunit
type(time_type)       :: model_time, model_time_base, model_time_offset
real(r4), allocatable :: statevector(:,:)
real(r8), allocatable :: latitude(:)
real(r8), allocatable :: longitude(:)
real(digits12) :: time_offset_seconds

integer  :: ncid, NtheDimID, NphiDimID, VarID, LonVarID, LatVarID
integer  :: TimeVarID

integer  :: iyear, imonth, iday, ihour, iminute, isecond
integer  :: ilat, ilon
integer  :: nthe, nphi
real(r8) :: dlat, dlon

character(len=512) :: string1, string2

!=======================================================================

call initialize_utilities(progname='openggcm_to_netcdf')

!-----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.

call find_namelist_in_file("input.nml", "openggcm_to_netcdf_nml", iunit)
read(iunit, nml = openggcm_to_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "openggcm_to_netcdf_nml")

if (verbose) then
   write(string1,*)'..  converting openggcm binary file >'//trim(openggcm_to_netcdf_input_file)//'<'
   write(string2,*)'to netcdf file >'//trim(openggcm_to_netcdf_output_file)//'<'
   call error_handler(E_MSG,'openggcm_to_netcdf',string1, text2=string2)
endif

!-----------------------------------------------------------------------

iunit = open_file(openggcm_to_netcdf_input_file, form='unformatted',action='read')

read(iunit,iostat=io) nphi, nthe
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read nphi, nthe'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
elseif (verbose) then
   write(string1,'(''nphi, nthe'',2(1x,i6))'), nphi, nthe
   call error_handler(E_MSG,'openggcm_to_netcdf',string1)
endif

read(iunit,iostat=io) iyear, imonth, iday, ihour, iminute, isecond
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read time record of 6 integers'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
elseif (verbose) then
   write(string1,'(''iyear, imonth, iday, ihour, iminute, isecond '',i4,4(1x,i2),1x,i5)') &
                    iyear, imonth, iday, ihour, iminute, isecond
   call error_handler(E_MSG,'openggcm_to_netcdf',string1)
endif

allocate(statevector(nphi,nthe))
allocate(latitude(nthe), longitude(nphi))
read(iunit,iostat=io) statevector
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read the "pot" variable.'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
endif

call close_file(iunit)

! The test file has year since 1900 ... 
iyear = iyear + 1900

dlat = 180.0_r8/real(nthe-1,r8)
do ilat = 1,nthe  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   latitude(ilat) = dlat*real(ilat-1,r8)
enddo

dlon = 360.0_r8/real(nphi-1,r8)
do ilon = 1,nphi  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   longitude(ilon) = dlon*real(ilon-1,r8)
enddo

call set_calendar_type('Gregorian')
model_time      = set_date(iyear, imonth, iday, ihour, iminute, isecond)
model_time_base = set_date(1900, 1, 1, 0, 0, 0)
model_time_offset = model_time - model_time_base
call get_time(model_time_offset, isecond, iday)

time_offset_seconds = real(iday,digits12) + real(isecond,digits12)/86400.0_digits12

call print_date(model_time, str='openggcm_to_netcdf:openggcm  model date')
call print_time(model_time, str='openggcm_to_netcdf:in DART timeframe')

!-----------------------------------------------------------------------
! Write the netCDF file.

call nc_check(nf90_create(openggcm_to_netcdf_output_file, NF90_CLOBBER, ncid),'openggcm_to_netcdf')

call nc_check(nf90_def_dim(ncid, 'nphi', nphi, NphiDimID),'def_dim nphi')
call nc_check(nf90_def_dim(ncid, 'nthe', nthe, NtheDimID),'def_dim nthe')

call nc_check(nf90_def_var(ncid, 'time', nf90_double, TimeVarID),'def_var time')
call nc_check(nf90_put_att(ncid, TimeVarID, 'units', 'days since 1900-01-01'), 'put_att:time:units')

call nc_check(nf90_def_var(ncid, 'nphi', nf90_real, (/ NphiDimID /), LonVarID),'def_var:nphi')
call nc_check(nf90_put_att(ncid, LonVarID, 'short_name', 'geomagnetic longitude'), &
                  'put_att nphi:short_name')

call nc_check(nf90_def_var(ncid, 'nthe', nf90_real, (/ NtheDimID /), LatVarID),'def_var:nthe')
call nc_check(nf90_put_att(ncid, LatVarID, 'short_name', 'geomagnetic latitude'), &
                  'put_att nthe:short_name')

call nc_check(nf90_def_var(ncid, 'pot' , nf90_real, (/ NphiDimID, NtheDimID /), VarID),'def_var:pot')

call nc_check(nf90_put_att(ncid, VarID, 'short_name', 'electric potential'), &
                  'put_att pot:short_name')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'ionosphere electric potential'), &
                  'put_att pot:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units', 'volts'), &
                  'put_att pot:units')

call nc_check(nf90_enddef(ncid),'enddef')

call nc_check(nf90_put_var(ncid,  LonVarID,           longitude),'put_var: lon')
call nc_check(nf90_put_var(ncid,  LatVarID,            latitude),'put_var: lat')
call nc_check(nf90_put_var(ncid,     VarID,         statevector),'put_var: pot')
call nc_check(nf90_put_var(ncid, TimeVarID, time_offset_seconds),'put_var: time')

call nc_check(nf90_close(ncid),'openggcm_to_netcdf')

call finalize_utilities('openggcm_to_netcdf')

end program openggcm_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
