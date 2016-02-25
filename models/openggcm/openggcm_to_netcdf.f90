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

use       netcdf_mod, only :  wr_netcdf_model_time, &
                              wr_netcdf_ctim_grid, &
                              wr_netcdf_interface_grid, &
                              wr_netcdf_r4_2D, &
                              wr_netcdf_r8_3D

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

namelist /openggcm_to_netcdf_nml/ openggcm_to_netcdf_input_file, &
                                  openggcm_to_netcdf_output_file, &
                                  verbose

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

integer               :: io, iunit
type(time_type)       :: model_time, model_time_base, model_time_offset
real(r8), allocatable :: tensor3D(:,:,:)
real(r4), allocatable :: tensor2D(:,:)
real(digits12) :: time_offset_seconds

integer  :: ncid

integer  :: iyear, imonth, iday, ihour, iminute, isecond
integer  :: nthe, nphi, nz

character(len=512) :: string1, string2

!=======================================================================

call initialize_utilities(progname='openggcm_to_netcdf')

!-----------------------------------------------------------------------
! Read the namelist to get the input and output filenames.

call find_namelist_in_file("input.nml", "openggcm_to_netcdf_nml", iunit)
read(iunit, nml = openggcm_to_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "openggcm_to_netcdf_nml")

if (verbose) then
   write(string1,*)"..  converting openggcm binary file '"//trim(openggcm_to_netcdf_input_file)//"'"
   write(string2,*)"to netcdf file '"//trim(openggcm_to_netcdf_output_file)//"'"
   call error_handler(E_MSG,'openggcm_to_netcdf',string1, text2=string2)
endif

call set_calendar_type('Gregorian')

!-----------------------------------------------------------------------
! Write the netCDF file.

call nc_check(nf90_create(openggcm_to_netcdf_output_file, NF90_CLOBBER, ncid),'openggcm_to_netcdf')

!-----------------------------------------------------------------------

iunit = open_file(openggcm_to_netcdf_input_file, form='unformatted',action='read')

read(iunit,iostat=io) iyear, imonth, iday, ihour, iminute, isecond
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read time record of 6 integers'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
elseif (verbose) then
   write(string1,'(''iyear, imonth, iday, ihour, iminute, isecond '',i4,4(1x,i2),1x,i5)') &
                    iyear, imonth, iday, ihour, iminute, isecond
   call error_handler(E_MSG,'openggcm_to_netcdf',string1)
endif
iyear = iyear + 1900
model_time      = set_date(iyear, imonth, iday, ihour, iminute, isecond)
model_time_base = set_date(1966, 1, 1, 0, 0, 0)
model_time_offset = model_time - model_time_base
call get_time(model_time_offset, isecond, iday)
time_offset_seconds = real(iday,digits12)*86400.0_digits12 + real(isecond,digits12)

call wr_netcdf_model_time(ncid, time_offset_seconds)

read(iunit,iostat=io) nphi, nthe, nz
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read nphi, nthe, nz'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
elseif (verbose) then
   write(string1,'(''nphi, nthe, nz'',3(1x,i6))'), nphi, nthe, nz
   call error_handler(E_MSG,'openggcm_to_netcdf',string1)
endif

allocate(tensor3D(nphi,nthe,nz))

read(iunit,iostat=io) tensor3D
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read the "pot" variable.'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
endif

call close_file(iunit)

call wr_netcdf_ctim_grid(ncid,nphi,nthe,nz)
call wr_netcdf_r8_3D(ncid,'ctim',nphi,nthe,nz,tensor3D,'ion','degrees kelvin','ionospheric miracle')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

iunit = open_file('../data/da0002.dart.pot', form='unformatted',action='read')

read(iunit,iostat=io) nphi, nthe
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read nphi, nthe, nz'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
elseif (verbose) then
   write(string1,'(''nphi, nthe, nz'',3(1x,i6))'), nphi, nthe, nz
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

allocate(tensor2D(nphi,nthe))

read(iunit,iostat=io) tensor2D
if (io /= 0) then
   write(string1,*)'read error was ',io,' when trying to read the "pot" variable.'
   call error_handler(E_ERR,'real_obs_sequence', string1, source, revision, revdate)
endif
call close_file(iunit)

call wr_netcdf_interface_grid(ncid,nphi,nthe,7)
call wr_netcdf_r4_2D(ncid,'interface',nphi,nthe,tensor2D,'pot','degrees trout','potential miracle')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

call nc_check(nf90_close(ncid),'openggcm_to_netcdf')

call finalize_utilities('openggcm_to_netcdf')

end program openggcm_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
