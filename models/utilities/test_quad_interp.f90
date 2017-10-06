! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program test_quad_interp

! intended to show how the state structure and quad code can be used
! together.  start with a simple regional grid and work out from there.

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8, deg2rad, rad2deg
use    utilities_mod, only : error_handler, open_file, close_file
use   random_seq_mod, only : init_random_seq, random_seq_type, &
                             random_uniform, random_gaussian

use quad_utils_mod, only : quad_interp_handle, init_quad_interp, finalize_quad_interp, set_quad_coords,       &
                           quad_lon_lat_locate, quad_lon_lat_evaluate, GRID_QUAD_FULLY_REGULAR,               &
                           GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR, GRID_QUAD_UNKNOWN_TYPE, &
                           QUAD_LOCATED_UNKNOWN, QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES,           &
                           QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS


implicit none


type(quad_interp_handle) :: h

! data grid size
integer, parameter :: nx = 9
integer, parameter :: ny = 5

! locations of data grid corners
real(r8) :: data_lons(nx, ny) = MISSING_R8
real(r8) :: data_lats(nx, ny) = MISSING_R8

! extents of the data grid (these mimic a regional model's grid)
real(r8) :: start_lon = 100.0_r8
real(r8) :: end_lon   = 150.5_r8
real(r8) :: start_lat = -11.4_r8
real(r8) :: end_lat   =  34.1_r8

! angle to rotate data grid in degrees
! positive is counterclockwise; will rotate
! around lower left grid point (start lon/lat).
!real(r8) :: angle = 45.0_r8
!real(r8) :: angle = 30.0_r8
!real(r8) :: angle =  90.0_r8
!real(r8) :: angle = -30.0_r8
!real(r8) :: angle = -10.0_r8
 real(r8) :: angle = 0.0_r8

! deform grid by this fraction of the deltas
real(r8) :: lon_def = 0.3_r8
real(r8) :: lat_def = 0.3_r8

! data values on the grid
real(r8) :: grid_data(nx, ny) = MISSING_R8


! sampling grid size
integer, parameter :: nrx = 210
integer, parameter :: nry = 150

! locations of sampling grid
real(r8) :: reg_lons(nrx) = MISSING_R8
real(r8) :: reg_lats(nrx) = MISSING_R8

! extents of the sampling grid
real(r8) :: reg_start_lon = 110.0_r8
real(r8) :: reg_end_lon   = 140.0_r8
real(r8) :: reg_start_lat = -20.0_r8
real(r8) :: reg_end_lat   =  30.0_r8

! where interpolated values are stored on reg grid
real(r8) :: interp_data(nrx, nry) = MISSING_R8


type(random_seq_type) :: ran

integer  :: i, j
real(r8) :: lon_del, lat_del, reg_lon_del, reg_lat_del
integer  :: lon_bot, lat_bot, lon_top, lat_top
integer  :: istatus
real(r8) :: invals(4), outval
integer  :: iunit_orig, iunit_interp

call init_random_seq(ran)

iunit_orig = open_file('original_data.txt', action='write')
iunit_interp = open_file('interp_data.txt', action='write')

lon_del = (end_lon - start_lon) / (nx-1)
lat_del = (end_lat - start_lat) / (ny-1)

! "data grid" corners and data vals
do i=1, nx
   do j=1, ny
      data_lons(i, j) = start_lon + (i-1)*lon_del + deform(lon_del, lon_def, ran)
      data_lats(i, j) = start_lat + (j-1)*lat_del + deform(lat_del, lat_def, ran)
      if (angle /= 0.0_r8) &
         call rotate(data_lons(i, j), data_lats(i, j), angle, start_lon, start_lat)
      ! pick one:
      ! constant by row
      !grid_data(i, j) = (j-1)*nx
      ! based on lon only
      !grid_data(i, j) = data_lons(i, j)
      ! based on lat only
      !grid_data(i, j) = data_lats(i, j)
      ! increasing monotonically
      !grid_data(i, j) = (j-1)*nx + i
      ! random between (0-10]
      grid_data(i, j) = random_uniform(ran) * 10.0_r8
      write(iunit_orig, *) i, j, data_lons(i,j), data_lats(i, j), grid_data(i, j)
   enddo
enddo

reg_lon_del = (reg_end_lon - reg_start_lon) / nrx
reg_lat_del = (reg_end_lat - reg_start_lat) / nry

! "sampled grid" spacing along each axis
do i=1, nrx
   reg_lons(i) = reg_start_lon + (i-1)*reg_lon_del
enddo
do j=1, nry
   reg_lats(j) = reg_start_lat + (j-1)*reg_lat_del
enddo

! end of data setup - now call interp routines

!print *, 'calling init_quad_interp: nx,ny = ', nx, ny
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, QUAD_LOCATED_CELL_CENTERS, .false., .false., .false., h)
!print *, 'calling set_quad_coords'
call set_quad_coords(h, data_lons, data_lats)
!print *, 'done with setup'

! later, do this in a loop
! "sampled grid" spacing along each axis
do i=1, nrx
   do j=1, nry
      call quad_lon_lat_locate(h, reg_lons(i), reg_lats(j), lon_bot, lat_bot, lon_top, lat_top, istatus)
      if (istatus /= 0) then
         !print *, 'location outside of grid: ', reg_lons(i), reg_lats(j)
         interp_data(i, j) = MISSING_R8 
         cycle
      endif
      print *, i, j, lon_bot, lat_bot, lon_top, lat_top, reg_lons(i), reg_lats(j)

      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      invals(1) = grid_data(lon_bot, lat_bot)
      invals(2) = grid_data(lon_top, lat_bot)
      invals(3) = grid_data(lon_top, lat_top)
      invals(4) = grid_data(lon_bot, lat_top)

      call quad_lon_lat_evaluate(h, reg_lons(i), reg_lats(j), lon_bot, lat_bot, lon_top, lat_top, &
                                 invals, outval, istatus)

      interp_data(i, j) = outval

      write(iunit_interp, *) i, j, reg_lons(i), reg_lats(j), interp_data(i,j)

   enddo
enddo

! this program doesn't currently have any missing locations - but i'll test that next.
print *, 'number of missing values in  input data: ', count(grid_data(:,:) == MISSING_R8)
print *, 'number of missing values in output data: ', count(interp_data(:,:) == MISSING_R8)

call close_file(iunit_orig)
call close_file(iunit_interp)

call finalize_quad_interp(h)

print *, 'closed files and finalized interp handle'

contains

!------------------------------------------------------------
! rotate vector a counterclockwise by angle theta, relative
! to the given origin point.

subroutine rotate(x, y, theta, x0, y0)
 real(r8), intent(inout) :: x, y
 real(r8), intent(in)    :: theta
 real(r8), intent(in)    :: x0, y0

real(r8) :: a(2), b(2)
real(r8) :: r(2,2)
real(r8) :: rads

a(1) = x - x0
a(2) = y - y0

rads = theta * deg2rad

r(1,1) = cos(rads)
r(1,2) = sin(rads)
r(2,1) = sin(-rads)
r(2,2) = cos(rads)

b(1) = r(1,1)*a(1) + r(1,2)*a(2)
b(2) = r(2,1)*a(1) + r(2,2)*a(2)

x = b(1) + x0
y = b(2) + y0

end subroutine rotate

!------------------------------------------------------------
! compute +/- a random value based on a width and percentage
! of that width

function deform(width, fraction, seq)

use random_seq_mod

 real(r8), intent(in) :: width
 real(r8), intent(in) :: fraction
 type(random_seq_type), intent(inout) :: seq
 real(r8)             :: deform

real(r8) :: val

! random val between -1 and 1
val = (random_uniform(seq) * 2.0_r8) - 1.0_r8

deform = val * width * fraction

end function deform

!------------------------------------------------------------

end program test_quad_interp

! <next few lines under version control, do not edit>
! $URL$
! $Revision$
! $Date$
