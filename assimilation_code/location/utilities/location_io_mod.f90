! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> common routines which can write netcdf arrays no matter which
!> location module is compiled in.

!>@todo FIXME  define ALL the VERTISXXX here, add a call to has_vertical_choice
!>to the xxx/location_mod, and if no, return without needing routines
!>in each of the specific location routines.  if yes, call query and 
!>match it with the type to return true/false.  remove vert_is_xxx from 
!>all the other location routines!!!  :)
!>will have to replicate the VERTISxxx in the location modules because
!>fortran doesn't allow circular 'use's between modules.  ugh.

module location_io_mod

use            types_mod, only : r8, MISSING_I
use netcdf_utilities_mod, only : nc_check

!>@todo FIXME: should there be accessor functions for 5 LocationXXX variables below?
use        location_mod, only : location_type, get_location, &
                                LocationDims, LocationName, LocationLName, &
                                LocationStorageOrder, LocationUnits, &
                                has_vertical_choice, query_location


use typeSizes
use netcdf

implicit none
private

public :: nc_write_location_atts, nc_get_location_varids, &
          nc_write_location, nc_write_location_vert

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! should import these but they don't exist in low order locations mods
!>@todo define them all here and replicate in the ones which use them.

integer, parameter :: VERTISUNDEF       = -2  ! has no specific vertical location (undefined)
integer, parameter :: VERTISSURFACE     = -1  ! surface value (value is surface elevation in m)
integer, parameter :: VERTISLEVEL       =  1  ! by level
integer, parameter :: VERTISPRESSURE    =  2  ! by pressure (in pascals)
integer, parameter :: VERTISHEIGHT      =  3  ! by height (in meters)
integer, parameter :: VERTISSCALEHEIGHT =  4  ! by scale height (unitless)

interface nc_write_location
   module procedure nc_write_single_location
   module procedure nc_write_multiple_locations
end interface

contains

!----------------------------------------------------------------------------
!> Create and add attributes to a 'location' dimension and variable.

!>@todo FIXME does the last arg need to be an optional actual dim?
!>or an additional dim?  check obs_seq_verify for usage

subroutine nc_write_location_atts(ncFileID, dimlen, use_dimID, fname) 
 
integer,                    intent(in) :: ncFileID    ! handle to the netcdf file
integer,                    intent(in) :: dimlen      ! number of locations to be created
integer,          optional, intent(in) :: use_dimID   ! if other than locations dim, use this
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocDimID, LDimID, VarID
integer :: rc


! define the rank/dimension of the location information
rc = nf90_def_dim(ncid=ncFileID, name='location', len=dimlen, dimid=LocDimID)
call checkit(rc, 'nc_write_location_atts', 'def_dim:location', fname)

if (LocationDims > 1) then
   rc = nf90_def_dim(ncid=ncFileID, name='locdim', len=LocationDims, dimid=LDimID)
   call checkit(rc, 'nc_write_location_atts', 'def_dim:locdim', fname)
endif

! Define the location variable and attributes

if (LocationDims > 1) then
   if (present(use_dimID)) then
      call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                dimids=(/ LDimID, use_dimID /), varid=VarID), &
               'nc_write_location_atts', 'location:def_var')
   else
      call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                dimids=(/ LDimID, LocDimID /), varid=VarID), &
               'nc_write_location_atts', 'location:def_var')
   endif
else
   if (present(use_dimID)) then
      call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                dimids=(/ use_dimID /), varid=VarID), &
               'nc_write_location_atts', 'location:def_var')
   else
      call nc_check(nf90_def_var(ncFileID, 'location', xtype=nf90_double, &
                dimids=(/ LocDimID /), varid=VarID), &
               'nc_write_location_atts', 'location:def_var')
   endif
endif

call nc_check(nf90_put_att(ncFileID, VarID, 'description', 'location coordinates'), &
              'nc_write_location_atts', 'location:description')
call nc_check(nf90_put_att(ncFileID, VarID, 'location_type', trim(LocationName)), &
              'nc_write_location_atts', 'location:location_type')
call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(LocationLName)), &
              'nc_write_location_atts', 'location:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'storage_order', trim(LocationStorageOrder)),  &
              'nc_write_location_atts', 'location:storage_order')
call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(LocationUnits)),   &
              'nc_write_location_atts', 'location:units')

end subroutine nc_write_location_atts

!----------------------------------------------------------------------------
!> Define the ancillary vertical array and attributes

subroutine nc_write_location_vert( ncFileID, fname )

integer,           intent(in) :: ncFileID    ! handle to the netcdf file
character(len=*),  intent(in) :: fname       ! file name (for printing purposes)

integer :: VarID

call nc_check(nf90_def_var(ncid=ncFileID, name='which_vert', xtype=nf90_int, &
          dimids=(/ nf90_unlimited /), varid=VarID), &
            'nc_write_location_vert', 'which_vert:def_var')

call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', 'vertical coordinate system code'), &
           'nc_write_location_vert', 'which_vert:long_name')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISUNDEF', VERTISUNDEF), &
           'nc_write_location_vert', 'which_vert:VERTISUNDEF')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSURFACE', VERTISSURFACE), &
           'nc_write_location_vert', 'which_vert:VERTISSURFACE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISLEVEL', VERTISLEVEL), &
           'nc_write_location_vert', 'which_vert:VERTISLEVEL')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISPRESSURE', VERTISPRESSURE), &
           'nc_write_location_vert', 'which_vert:VERTISPRESSURE')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISHEIGHT', VERTISHEIGHT), &
           'nc_write_location_vert', 'which_vert:VERTISHEIGHT')
call nc_check(nf90_put_att(ncFileID, VarID, 'VERTISSCALEHEIGHT', VERTISSCALEHEIGHT), &
           'nc_write_location_vert', 'which_vert:VERTISSCALEHEIGHT')

end subroutine nc_write_location_vert

!----------------------------------------------------------------------------
!> Return the LocationVarID and WhichVertVarID variables from a given netCDF file.
!>@todo FIXME why do we need this?
!>
!> ncFileId         the netcdf file descriptor
!> fname            the name of the netcdf file (for error messages only)
!> LocationVarID    the integer ID of the 'location' variable in the netCDF file
!> WhichVertVarID   the integer ID of the 'which_vert' variable in the netCDF file

subroutine nc_get_location_varids(ncFileID, LocationVarID, WhichVertVarID, fname)

integer,                    intent(in)  :: ncFileID   ! handle to the netcdf file
integer,                    intent(out) :: LocationVarID
integer,          optional, intent(out) :: WhichVertVarID
character(len=*), optional, intent(in)  :: fname      ! file name (for printing purposes)

integer :: rc

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call checkit(rc, 'nc_get_location_varids', 'inq_varid:location ', fname)

if (present(WhichVertVarID)) then
  rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
  if (rc /= NF90_NOERR) then
      WhichVertVarID = MISSING_I
  endif
endif

end subroutine nc_get_location_varids

!----------------------------------------------------------------------------
!> Writes a SINGLE location to the specified netCDF variable and file.
!> The LocationVarID and WhichVertVarID must be the values returned from
!> the nc_get_location_varids call.

subroutine nc_write_single_location(ncFileID, loc, locindex, do_vert, fname)
 
integer,             intent(in) :: ncFileID
type(location_type), intent(in) :: loc
integer,             intent(in) :: locindex
logical, optional,   intent(in) :: do_vert
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocationVarID
integer :: WhichVertVarID
real(r8), dimension(LocationDims) :: locations
integer,  dimension(1) :: intval
logical :: write_vert
integer :: rc

write_vert = .false.
if (present(do_vert)) write_vert = do_vert

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call checkit(rc, 'nc_write_single_location', 'inq_varid:location ', fname)

locations = get_location( loc ) 

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
              start=(/ 1, locindex /), count=(/ LocationDims, 1 /) ), &
              'nc_write_single_location', 'put_var:location')

if (write_vert) then
   rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
   if (rc /= NF90_NOERR) then
     intval = query_location(loc, 'WHICH_VERT')
     call nc_check(nf90_put_var(ncFileID, WhichVertVarID, intval, &
                   start=(/ locindex /), count=(/ 1 /) ), &
                   'nc_write_single_location','put_var:vert' )
   endif
endif

end subroutine nc_write_single_location

!----------------------------------------------------------------------------
!> Writes an array of locations to the specified netCDF variable and file.
!> The LocationVarID and WhichVertVarID must be the values returned from
!> the nc_get_location_varids call.

subroutine nc_write_multiple_locations(ncFileID, loc, loccount, startlocindex, do_vert, fname)

integer,             intent(in) :: ncFileID
type(location_type), intent(in) :: loc(:)
integer,             intent(in) :: loccount
integer, optional,   intent(in) :: startlocindex
logical, optional,   intent(in) :: do_vert
character(len=*), optional, intent(in) :: fname       ! file name (for error printing purposes)

integer :: LocationVarID
integer :: WhichVertVarID
real(r8), allocatable :: locations(:,:)
integer,  allocatable :: intvals(:)
logical :: write_vert
integer :: rc, i, starthere

write_vert = .false.
if (present(do_vert)) write_vert = do_vert

starthere = 1
if (present(startlocindex)) starthere = startlocindex

rc = nf90_inq_varid(ncFileID, 'location', varid=LocationVarID)
call checkit(rc, 'nc_write_multiple_locations', 'inq_varid:location ', fname)

allocate(locations(LocationDims,loccount))
if (write_vert) allocate(intvals(loccount))

do i=1, loccount
   locations(:,i) = get_location( loc(i) ) 
   if (write_vert) intvals(i) = query_location(loc(i), 'WHICH_VERT')
enddo

call nc_check(nf90_put_var(ncFileID, LocationVarId, locations, &
              start=(/ 1, starthere /), count=(/ LocationDims, loccount /) ), &
              'nc_write_multiple_locations', 'put_var:location')

if (write_vert) then
   rc = nf90_inq_varid(ncFileID, 'which_vert', varid=WhichVertVarID)
   if (rc /= NF90_NOERR) then
     call nc_check(nf90_put_var(ncFileID, WhichVertVarID, intvals, &
                   start=(/ starthere /), count=(/ loccount /) ), &
                   'nc_write_multiple_locations','put_var:vert' )
   endif
endif

deallocate(locations)
if (write_vert) deallocate(intvals)

end subroutine nc_write_multiple_locations

!----------------------------------------------------------------------------

subroutine checkit(rc, subname, action, fname)

integer,                    intent(in) :: rc
character(len=*),           intent(in) :: subname
character(len=*),           intent(in) :: action
character(len=*), optional, intent(in) :: fname

if (present(fname)) then
   call nc_check(rc, subname, action//trim(fname))
else
   call nc_check(rc, subname, action)
endif

end subroutine checkit

!----------------------------------------------------------------------------

end module location_io_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
