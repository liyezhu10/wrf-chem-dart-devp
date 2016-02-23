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

use        types_mod, only : r4
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, nc_check
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

character(len=256) :: openggcm_to_netcdf_input_file  = '../data/da0002.dart.pot'
character(len=256) :: openggcm_to_netcdf_output_file = 'openggcm.nc'

namelist /openggcm_to_netcdf_nml/ openggcm_to_netcdf_output_file, &
                                  openggcm_to_netcdf_input_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit
type(time_type)       :: model_time
real(r4), allocatable :: statevector(:)
real(r4), allocatable :: datmat(:,:)

integer :: nthe, nphi

integer :: ncid, NtheDimID, NphiDimID, VarID

!----------------------------------------------------------------------

call initialize_utilities(progname='openggcm_to_netcdf')

!----------------------------------------------------------------------
! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, starting date, deltaT, etc.
!----------------------------------------------------------------------

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
write(*,*)'TJH nphi, nthe is ',nphi, nthe
allocate(statevector(nphi*nthe))
read(iunit) statevector
call close_file(iunit)

write(*,*)'file size should be ',4 + 4+4 + 4 + 4 + nphi*nthe*4 + 4

call nc_check(nf90_create(openggcm_to_netcdf_output_file, NF90_CLOBBER, ncid),'openggcm_to_netcdf')

! define dimensions

call nc_check(nf90_def_dim(ncid, 'nphi', nphi, NphiDimID),'openggcm_to_netcdf')
call nc_check(nf90_def_dim(ncid, 'nthe', nthe, NtheDimID),'openggcm_to_netcdf')
call nc_check(nf90_def_var(ncid, 'pot', nf90_real, (/ NphiDimID, NtheDimID /), VarID),'openggcm_to_netcdf')
call nc_check(nf90_enddef(ncid),'openggcm_to_netcdf')

datmat = reshape(statevector, (/ nphi, nthe /))
call nc_check(nf90_put_var(ncid,  VarID,  datmat),'openggcm_to_netcdf')

call nc_check(nf90_close(ncid),'openggcm_to_netcdf')

!call print_date(model_time, str='openggcm_to_netcdf:openggcm  model date')
!call print_time(model_time, str='openggcm_to_netcdf:DART model time')
call finalize_utilities('openggcm_to_netcdf')

end program openggcm_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
