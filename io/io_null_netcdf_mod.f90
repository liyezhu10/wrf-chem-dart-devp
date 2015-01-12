! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: io_filenames_mod.f90 7314 2014-12-23 20:11:42Z hkershaw $

module io_netcdf_mod

!> \defgroup io_filenames io_filenames
!> Aim is to store the io filenames
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Any module can set the filenames here, then state_vector_io_mod
!> can read from this module to get the filenames.
!> Maybe this is a bit lazy, but I just want a way for the 
!> different modules to set the filenames
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use utilities_mod,    only : do_nml_file, nmlfileunit, do_nml_term, check_namelist_read, &
                             find_namelist_in_file, nc_check
use time_manager_mod, only : time_type, set_date !> idealy would not have these here
use netcdf

implicit none

private :: get_date

! These should probably be set and get functions rather than 
! direct access

public :: get_variable_list,    &
          get_info_file_name,   &
          get_model_time,       &
          get_variables_domains

contains

!--------------------------------------------------------------------
!> read namelist and set up filename arrays
function get_variable_list(num_variables_in_state)

integer             :: num_variables_in_state
character(len=256)  :: get_variable_list(num_variables_in_state)

end function get_variable_list

!--------------------------------------------------------------------
!> get the input file name
function get_info_file_name(domain)

integer, intent(in) :: domain
character(len=265)  :: get_info_file_name

end function get_info_file_name

!--------------------------------------------------------------------
!> read the time from the input file
function get_model_time(filename) 

character(len=1024), intent(in) :: filename

end function get_model_time

!--------------------------------------------------------------------
!> pass number of variables in the state out to filter 
subroutine get_variables_domains(num_variables_in_state, num_doms)

integer, intent(out) :: num_variables_in_state
integer, intent(out) :: num_doms !< number of domains

end subroutine get_variables_domains

!--------------------------------------------------------------------
!>  get the date from a tstring
subroutine get_date(tstring, year, month,  day, hour, minute, second)

integer,           intent(out) :: year, month, day, hour, minute, second
character(len=19), intent(in)  :: tstring

read(tstring( 1: 4),'(i4)') year
read(tstring( 6: 7),'(i2)') month
read(tstring( 9:10),'(i2)') day
read(tstring(12:13),'(i2)') hour
read(tstring(15:16),'(i2)') minute
read(tstring(18:19),'(i2)') second

return

end subroutine get_date

end module io_netcdf_mod
