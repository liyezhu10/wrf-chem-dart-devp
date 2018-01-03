! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed sphere locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check, nc_add_global_creation_time,        & 
                                  nc_create_file, nc_close_file, nc_begin_define_mode, &
                                  nc_define_dimension, nc_put_variable,         &
                                  nc_add_attribute_to_variable,                 &
                                  nc_define_double_variable, nc_end_define_mode

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location, query_location, &
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate


implicit none

public :: test_interpolate_single, &
          test_interpolate_range, &
          find_closest_state_item

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
! for messages
character(len=512) :: string1, string2, string3

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of lat, lon, and vert values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

! Do a interpolation on a range of lat, lon, vert values.  Returns the
! number of failures.
!>@todo FIXME: this is totally unaware of MPI.  for now i'm going to change
! the code to only run on task 0.  but maybe we should run this on all tasks
! and have each write out to a different file with the task number in the name?
!-------------------------------------------------------------------------------
function test_interpolate_range( ens_handle,            &
                                 ens_size,              &
                                 interp_test_dlon,      &
                                 interp_test_dlat,      &
                                 interp_test_dvert,     &
                                 interp_test_vertcoord, &
                                 interp_test_lonrange,  &
                                 interp_test_latrange,  &
                                 interp_test_vertrange, &
                                 quantity_string,       &
                                 verbose )

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
real(r8)              , intent(in)    :: interp_test_dlon
real(r8)              , intent(in)    :: interp_test_dlat
real(r8)              , intent(in)    :: interp_test_dvert
character(len=*)      , intent(in)    :: interp_test_vertcoord
real(r8), dimension(2), intent(in)    :: interp_test_latrange
real(r8), dimension(2), intent(in)    :: interp_test_lonrange
real(r8), dimension(2), intent(in)    :: interp_test_vertrange
character(len=*),       intent(in)    :: quantity_string
logical               , intent(in)    :: verbose

integer :: test_interpolate_range

! Local variables

character(len=*), parameter :: routine = 'test_interpolate_range'

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
real(r8) :: lonrange_top
integer :: nlon, nlat, nvert
integer :: ilon, jlat, kvert, nfailed
character(len=128) :: ncfilename, txtfilename
character(len=256)  :: output_file = 'check_me'
character(len=32)   :: field_name
type(location_type) :: loc
integer :: ncid
integer :: iunit, ios_out(ens_size), imem
integer :: quantity_index, vertcoord
integer,  allocatable :: all_ios_out(:,:)

test_interpolate_range = 0

! only do this on a single task.  this should be fixed.
if (.not. do_output()) return

if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   if ( do_output() ) then
      write(*,'(A)')    'Skipping the rigorous interpolation test because one of'
      write(*,'(A)')    'interp_test_dlon,interp_test_dlat are < 0.0'
      write(*,'(A,I2)') 'interp_test_dlon  = ',interp_test_dlon
      write(*,'(A,I2)') 'interp_test_dlat  = ',interp_test_dlat
      write(*,'(A,I2)') 'interp_test_dvert = ',interp_test_dvert
   endif
   return
endif

vertcoord = convert_string_to_index(interp_test_vertcoord)
quantity_index = get_index_for_quantity(quantity_string)

!>@todo FIXME add the task number to the names?
write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! for longitude, allow wrap.
lonrange_top = interp_test_lonrange(2)
if (interp_test_lonrange(2) < interp_test_lonrange(1)) &
   lonrange_top = interp_test_lonrange(2) + 360.0_r8

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) - interp_test_latrange(1))  / interp_test_dlat) + 1
nlon  = aint((            lonrange_top - interp_test_lonrange(1))  / interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1)) / interp_test_dvert) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''nvert = '',i8,'';'')')nvert
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,ens_size))
allocate(all_ios_out(nlon*nlat*nvert,ens_size))
nfailed = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   if (lon(ilon) >= 360.0_r8) lon(ilon) = lon(ilon) - 360.0_r8
   if (lon(ilon) < 0.0_r8)   lon(ilon) = lon(ilon) + 360.0_r8
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat
      do kvert = 1, nvert
         vert(kvert) = interp_test_vertrange(1) + real(kvert-1,r8) * interp_test_dvert

         loc = set_location(lon(ilon), lat(jlat), vert(kvert), vertcoord)

         call model_interpolate(ens_handle, ens_size, loc, quantity_index, &
                                field(ilon,jlat,kvert,:), ios_out)
         write(iunit,*) field(ilon,jlat,kvert,:)
         if (any(ios_out /= 0)) then
            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                                 ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
               call error_handler(E_MSG, routine, string1, &
                                  source, revision, revdate, text2=string2)
            endif
            nfailed = nfailed + 1
            all_ios_out(nfailed,:) = ios_out
         endif

      enddo
   enddo
enddo

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon,nens);'')')
write(iunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nlon*nlat*nvert
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out, nfailed)


! Write out the netCDF file for easy exploration.

ncid = nc_create_file(ncfilename, routine)
call nc_add_global_creation_time(ncid, 'test_interpolate_range', 'creation put '//trim(ncfilename))

! Define dimensions

call nc_define_dimension(ncid, 'lon',  nlon,  routine)
call nc_define_dimension(ncid, 'lat',  nlat,  routine)
call nc_define_dimension(ncid, 'vert', nvert, routine)

! Define variables

call nc_define_double_variable(   ncid, 'lon', 'lon', routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'range', interp_test_lonrange, routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'cartesian_axis', 'X', routine)

call nc_define_double_variable(   ncid, 'lat', 'lat', routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'range', interp_test_latrange, routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'cartesian_axis', 'Y', routine)

call nc_define_double_variable(   ncid, 'vert', 'vert', routine)
call nc_add_attribute_to_variable(ncid, 'vert', 'range', interp_test_vertcoord, routine)
call nc_add_attribute_to_variable(ncid, 'vert', 'cartesian_axis', 'Z', routine)

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_define_double_variable(ncid, field_name, (/ 'lon ', 'lat ', 'vert' /), routine)

   call nc_add_attribute_to_variable(ncid, field_name, 'long_name', quantity_string, routine)
   call nc_add_attribute_to_variable(ncid, field_name, '_FillValue',    MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'missing_value', MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'interp_test_vertcoord', interp_test_vertcoord, routine)
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_put_variable(ncid, 'lon',  lon,  routine)
call nc_put_variable(ncid, 'lat',  lat,  routine)
call nc_put_variable(ncid, 'vert', vert, routine)

do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_put_variable(ncid, field_name, field(:,:,:,imem), routine)
enddo


call nc_close_file(ncid, routine)

deallocate(lon, lat, vert)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
!> Do a single interpolation on a given location and kind.
!> Returns the interpolated values and ios_out.
!> Returns the number of ensemble members that passed.

function test_interpolate_single( ens_handle,       &
                                  ens_size,         &
                                  vertcoord_string, &
                                  lonval,           &
                                  latval,           &
                                  vertval,          &
                                  quantity_string,  &
                                  interp_vals,      &
                                  ios_out)

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
character(len=*)      , intent(in)    :: vertcoord_string
real(r8)              , intent(in)    :: lonval
real(r8)              , intent(in)    :: latval
real(r8)              , intent(in)    :: vertval
character(len=*)      , intent(in)    :: quantity_string
real(r8)              , intent(out)   :: interp_vals(ens_size)
integer               , intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: loc
integer :: imem, num_passed, vertcoord
character(len=128) :: my_location
integer :: quantity_index

quantity_index = get_index_for_quantity(quantity_string)
vertcoord = convert_string_to_index(vertcoord_string)

loc = set_location(lonval, latval, vertval, vertcoord)

if ( do_output() ) then
   call write_location(0, loc, charstring=my_location)
   write(*,'(A)') ''
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'("interpolating at ",A)') trim(my_location)//' for "'//trim(quantity_string)//'"'
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(A)') ''
endif

call model_interpolate(ens_handle, ens_size, loc, quantity_index, interp_vals, ios_out)

num_passed = 0
do imem = 1, ens_size
   if (ios_out(imem) == 0 ) then
      if (do_output()) then
         write(string1,*)'model_interpolate SUCCESS with value    :: ', interp_vals(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
         num_passed = num_passed + 1
      endif
   else
      if (do_output()) then
         write(string1,*)'model_interpolate ERROR with error code :: ', ios_out(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
      endif
   endif
enddo

if ( do_output() ) write(*,'(A)') ''

test_interpolate_single = num_passed

end function test_interpolate_single


!-------------------------------------------------------------------------------
!> Count the number of different error codes and output the results.
!> This is just a helper function for test_interpolate_range.
!> Only sums error codes for the first ensemble member.

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
!> need to convert the character string for the test vertical coordinate into
!> the corresponding dart index.

function  convert_string_to_index(test_vertcoord)

character(len=*) , intent(in) :: test_vertcoord

integer :: convert_string_to_index

select case (test_vertcoord)
   case ('VERTISUNDEF')
      convert_string_to_index = VERTISUNDEF
   case ('VERTISSURFACE')
      convert_string_to_index = VERTISSURFACE
   case ('VERTISLEVEL')
      convert_string_to_index = VERTISLEVEL
   case ('VERTISPRESSURE')
      convert_string_to_index = VERTISPRESSURE
   case ('VERTISHEIGHT')
      convert_string_to_index = VERTISHEIGHT
   case ('VERTISSCALEHEIGHT')
      convert_string_to_index = VERTISSCALEHEIGHT
   case default
      convert_string_to_index = VERTISUNDEF
end select

end function  convert_string_to_index

!-----------------------------------------------------------------------
!> Expensive exhaustive search to find the indices into the
!> state vector of a particular lon/lat/vert. At present, only for a
!> single variable - could be extended to identify the closest location
!> for every variable in each domain. This could help ensure grid
!> staggering is being handled correctly.

subroutine find_closest_state_item(loc_of_interest, vert_string, quantity_string)

real(r8),         intent(in) :: loc_of_interest(:)
character(len=*), intent(in) :: vert_string
character(len=*), intent(in) :: quantity_string

character(len=*), parameter :: routine = 'find_closest_state_item'

type(location_type)   :: loc0, loc1
integer(i8)           :: i
integer               :: quantity_index, var_type, state_vtype, vert_type
real(r8)              :: closest, rlon, rlat, rvert
logical               :: matched
real(r8), allocatable :: thisdist(:)
real(r8),   parameter :: FARAWAY = huge(r8)
character(len=metadatalength) :: myquantity

!>@todo there should be arrays of length state_structure_mod:get_num_variables(domid)
!>      get_num_domains(), get_num_variables() ...

allocate( thisdist(get_model_size()) )
thisdist  = FARAWAY
matched   = .false.

! Trying to support the ability to specify matching a particular QUANTITY.
! With staggered grids, the closest state_item might not be of the quantity
! of interest.

quantity_index = get_index_for_quantity(quantity_string)

rlon  = loc_of_interest(1)
rlat  = loc_of_interest(2)
rvert = loc_of_interest(3)

vert_type = convert_string_to_index(vert_string)
loc0 = set_location(rlon, rlat, rvert, vert_type)

write(string1,'("Computing indices into the state vector that are closest to")')
call write_location(0, loc0, charstring=string2)
write(string3,'("for (",A,") variables.")')trim(quantity_string)
call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through
! the array and come back to find all the 'identical' values.

DISTANCE : do i = 1,get_model_size()

   call get_state_meta_data(i, loc1, var_type)

   ! this doesn't support ALL
   if (var_type .ne. quantity_index) cycle DISTANCE

   ! Grab the vert_type from the grid and
   ! set out target location to have the same.
   ! Compute the distance.

   !>@todo FIXME you can't do this - what if the vertical in the namelist
   ! is in meters and the model uses pressure or level?

   state_vtype   = nint(query_location(loc1))
   if (state_vtype == vert_type) then
      thisdist(i) = get_dist( loc1, loc0, no_vert=.false.)
   else
      thisdist(i) = get_dist( loc1, loc0, no_vert=.true.)
   endif
   matched     = .true.

enddo DISTANCE

if (.not. matched) then
   write(string1,*)'No state vector elements of type "'//trim(quantity_string)//'"'
   call error_handler(E_MSG, routine, string1)
   deallocate( thisdist )
   return
endif

closest = minval(thisdist)

! Now that we know the distances ... report
! If more than one quantity has the same distance, report all.
! Be aware that if 'approximate_distance' is .true., everything
! in the box has a single location.

REPORT: do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1, var_type)
      myquantity = get_name_for_quantity(var_type)

      call write_location(0, loc1, charstring=string2)
      write(string1,'(A,I12,A)')trim(string2)//' is index ',i,' ('//trim(myquantity)//')'
      call error_handler(E_MSG, routine, string1)
   endif

enddo REPORT

deallocate( thisdist )

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
