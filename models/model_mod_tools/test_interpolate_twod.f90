! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for twod locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check, nc_add_global_creation_time, &
                                  nc_create_file, nc_close_file, nc_begin_define_mode, &
                                  nc_define_dimension, nc_put_variable, &
                                  nc_add_attribute_to_variable, &
                                  nc_define_double_variable, nc_end_define_mode 

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist

use          obs_kind_mod, only : get_name_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : model_interpolate

use netcdf

implicit none

public :: test_interpolate_range, &
          test_interpolate_single, &
          find_closest_state_item

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!-------------------------------------------------------------------------------
! Do a interpolation on a range of x and y values.  Returns the number of failures.
!-------------------------------------------------------------------------------
function test_interpolate_range( ens_handle,            &
                                 ens_size,              &
                                 interp_test_dx,        &
                                 interp_test_dy,        &
                                 interp_test_dz,        &
                                 interp_test_vertcoord, &
                                 interp_test_xrange,    &
                                 interp_test_yrange,    &
                                 interp_test_zrange,    &
                                 quantity_string,       &
                                 verbose )

type(ensemble_type),    intent(inout) :: ens_handle
integer,                intent(in)    :: ens_size
real(r8),               intent(in)    :: interp_test_dx
real(r8),               intent(in)    :: interp_test_dy
real(r8),               intent(in)    :: interp_test_dz
character(len=*),       intent(in)    :: interp_test_vertcoord
real(r8), dimension(2), intent(in)    :: interp_test_yrange
real(r8), dimension(2), intent(in)    :: interp_test_xrange
real(r8), dimension(2), intent(in)    :: interp_test_zrange
character(len=*),       intent(in)    :: quantity_string
logical,                intent(in)    :: verbose

! function to exercise the model_mod:model_interpolate() function
! This will result in a netCDF file with all salient metadata
integer :: test_interpolate_range

character(len=metadatalength) :: kind_of_interest

! Local variables

character(len=*), parameter :: routine = 'test_interpolate_range'

real(r8), allocatable :: X(:), Y(:)
real(r8), allocatable :: field(:,:,:)
integer :: nx, ny
integer :: i, j, nfailed
character(len=128) :: ncfilename,txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nxDimID, nyDimID
integer :: VarID(ens_size), XVarID, YVarID

character(len=256) :: output_file = 'check_me'

! for message strings
character(len=512) :: string1, string2

character(len=32)  :: field_name
type(location_type) :: loc
integer :: iunit, ios_out(ens_size), imem
integer, allocatable :: all_ios_out(:,:)

test_interpolate_range = 0

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! round down to avoid exceeding the specified range
ny  = aint(( interp_test_yrange(2) -  interp_test_yrange(1))/interp_test_dy) + 1
nx  = aint(( interp_test_xrange(2) -  interp_test_xrange(1))/interp_test_dx) + 1

allocate(X(nx), Y(ny), field(nx,ny,ens_size))
allocate(all_ios_out(nx*ny*nz,ens_size))
nfailed = 0

! matlab format
iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nx = '',i8,'';'')')nx
write(iunit,'(''ny = '',i8,'';'')')ny
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

do i = 1, nx
   X(i) = interp_test_xrange(1) + real(i-1,r8) * interp_test_dx
   do j = 1, ny
      Y(j) = interp_test_yrange(1) + real(j-1,r8) * interp_test_dy

      loc = set_location(X(i), Y(j))

      call model_interpolate(ens_handle, ens_size, loc, mykindindex, field(i,j,:), ios_out)

      write(iunit,*) field(i,j,:)
      if (any(ios_out /= 0)) then
        if (verbose) then
           write(string2,'(''i,j,X,Y'',2(1x,i6),2(1x,f14.6))') &
                       i,j,X(i),Y(j)
           write(string1,*) 'interpolation return code was', ios_out
           call error_handler(E_MSG,'test_interpolate_range',string1,source,revision,revdate,text2=string2)
        endif
        nfailed = nfailed + 1
        all_ios_out(nfailed,:) = ios_out
      endif

   end do
end do

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,ny,nx,nens);'')')
write(iunit,'(''datmat = permute(datmat,[3,1,2]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total interpolations  : ', nx
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out, nfailed)


! Write out the netCDF file for easy exploration.

ncid = nc_create_file(ncfilename, routine)
call nc_add_global_creation_time(ncid, 'test_interpolate_range', 'creation put '//trim(ncfilename))

! Define dimensions

call nc_define_dimension(ncid, 'X',  nx, routine)
call nc_define_dimension(ncid, 'Y',  ny, routine)

! Define variables

call nc_define_double_variable(   ncid, 'X', 'X', routine)
call nc_add_attribute_to_variable(ncid, 'X', 'range', interp_test_xrange, routine)
call nc_add_attribute_to_variable(ncid, 'X', 'cartesian_axis', 'X', routine)

call nc_define_double_variable(   ncid, 'Y', 'Y', routine)
call nc_add_attribute_to_variable(ncid, 'Y', 'range', interp_test_yrange, routine)
call nc_add_attribute_to_variable(ncid, 'Y', 'cartesian_axis', 'Y', routine)

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_define_double_variable(ncid, field_name, (/ 'X', 'Y' /), routine)

   call nc_add_attribute_to_variable(ncid, field_name, 'long_name', quantity_string, routine)
   call nc_add_attribute_to_variable(ncid, field_name, '_FillValue', MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'missing_value', MISSING_R8, routine)
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_put_variable(ncid, 'X', X, routine)
call nc_put_variable(ncid, 'Y', Y, routine)

do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_put_variable(ncid, field_name, field(:,:,imem), routine)
enddo

call nc_close_file(ncid, routine)


deallocate(X, Y, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range

!-------------------------------------------------------------------------------
! Do a single interpolation on a given location and kind.  Returns the
! interpolated values and ios_out. Returns the number of ensemble members that
! passed
!-------------------------------------------------------------------------------
function test_interpolate_single( ens_handle,       &
                                  ens_size,         &
                                  vertcoord_string, &
                                  xval,             &
                                  yval,             &
                                  zval,             &
                                  mykindindex,      &
                                  interp_vals,      &
                                  ios_out)

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
character(len=*)      , intent(in)    :: vertcoord_string
real(r8)              , intent(in)    :: xval
real(r8)              , intent(in)    :: yval
real(r8)              , intent(in)    :: vertval
integer               , intent(in)    :: mykindindex
real(r8)              , intent(out)   :: interp_vals(ens_size)
integer               , intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: loc
integer :: imem, num_passed
character(len=128) :: my_location

num_passed = 0

loc = set_location(xval, yval)

if ( do_output() ) then
   call write_location(0, loc, charstring=my_location)
   write(*,'(A)') ''
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'("interpolating at ",A)') trim(my_location)
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(A)') ''
endif

call model_interpolate(ens_handle, ens_size, loc, mykindindex, interp_vals, ios_out)

do imem = 1, ens_size
   if (ios_out(imem) == 0 ) then
      if (do_output() ) then
         write(*,'(A)') '-------------------------------------------------------------'
         write(*,'("member ",I3,", model_interpolate SUCCESS with value   :: ",F10.3)') imem, interp_vals(imem)
         write(*,'(A)') '-------------------------------------------------------------'
         num_passed = num_passed + 1
      endif
   else
      if (do_output()) then
         write(*,'(A)') '-------------------------------------------------------------'
         write(*,'("member ",I3,", model_interpolate ERROR with error code :: ",I2  )') imem, ios_out(imem)
         write(*,'(A)') '-------------------------------------------------------------'
      endif
   endif
enddo

test_interpolate_single = num_passed

end function test_interpolate_single

!-------------------------------------------------------------------------------
! Count the number of different error codes and output the results.  This
! is just a helper function for test_interpolate_range. Only sums error codes
! for the first ensemble member
!-------------------------------------------------------------------------------
subroutine count_error_codes(error_codes, num_failed)

integer, intent(in) :: error_codes(:,:)
integer, intent(in) :: num_failed

integer :: i, count_errors, results

count_errors = 1

i = 1
do while (count_errors < num_failed)
   results = count(error_codes(:,1) == i)
   if (results /= 0) then
      if ( do_output() ) &
         write(*,'(i10, a, i3)') results + 1, " failed with ios_out ", i
      count_errors = count_errors + results
   endif
   i = i+1
enddo

end subroutine count_error_codes

!-------------------------------------------------------------------------------

subroutine find_closest_state_item(loc_of_interest, vert_string,
quantity_string)

real(r8),         intent(in) :: loc_of_interest(:)
character(len=*), intent(in) :: vert_string
character(len=*), intent(in) :: quantity_string

character(len=*), parameter :: routine = 'find_closest_state_item'

!>@todo FIXME  write me

end subroutine find_closest_state_item

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
