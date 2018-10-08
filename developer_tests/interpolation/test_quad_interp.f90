! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id: test_quad_interp.f90 12852 2018-09-25 20:32:09Z thoar@ucar.edu $

! test cases for the quad interpolation routines.
! 
! eventually, want to support all of these combinations, but
! not there yet.
!
! data grids:
! 
! global grid, fully regular
! global grid, semi regular
! global grid, irregular
! 
! regional grid, fully regular
! regional grid, semi regular
! regional grid, irregular
! 
! with regional grids that:
! 
! cover poles
! cross greenwich line
! are cyclic in longitude
! are rotated
! 
! and data:
! 
! sin/cos pattern
! moving average (rows/cols)
! lon or lat vals (lon+lat?)
! 
! 
! notes: 
! 
! generating output suitable for plotting:
! 
! generate a sampling grid that's fully regular
! for ease in plotting.  make it higher density
! than the data grid to ensure smooth interpolation
! across the original cells.
! 
! automating tests more:
! 
! interpolate from data grid to sample grid, and
! then interpolate back to original data grid and
! quantify the residuals.  no plot output but print
! high and low values.
! 
! 
! make sampling grids that:
! 
! cover poles (regional)
! cross greenwich line (regional)
! are cyclic in longitude (regional)
! are global
! 

program test_quad_interp


use      types_mod, only : r8, i8, MISSING_R8, deg2rad, rad2deg
use  utilities_mod, only : error_handler, initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read, open_file, close_file, &
                           nmlfileunit, do_output, do_nml_file, do_nml_term
use random_seq_mod, only : init_random_seq, random_seq_type, random_uniform, random_gaussian
use parse_args_mod, only : get_args_from_string
use quad_utils_mod, only : quad_interp_handle, init_quad_interp, finalize_quad_interp, set_quad_coords,       &
                           quad_lon_lat_locate, quad_lon_lat_evaluate, GRID_QUAD_FULLY_REGULAR,               &
                           GRID_QUAD_IRREG_SPACED_REGULAR, GRID_QUAD_FULLY_IRREGULAR, GRID_QUAD_UNKNOWN_TYPE, &
                           QUAD_LOCATED_UNKNOWN, QUAD_LOCATED_CELL_CENTERS, QUAD_LOCATED_LON_EDGES,           &
                           QUAD_LOCATED_LAT_EDGES, QUAD_LOCATED_CELL_CORNERS


implicit none

! make a derived type to hold configs for a list of test cases
! always sampled on regular grid (for ease in plotting)
!
! for fully regular, need delta lon,lat and rep count or extents
! for partially regular, need full 1d lon,lat arrays 
! for fully irregular, need full 2d lon,lat arrays

! source data grid to be sampled
type src_testtype
 integer  :: grid_type
 integer  :: relative_cell_location
 integer  :: lon_count, lat_count
 real(r8) :: lon_start, lon_end
 real(r8) :: lat_start, lat_end
 real(r8), allocatable :: lon_1d_array(:),   lat_1d_array(:)
 real(r8), allocatable :: lon_2d_array(:,:), lat_2d_array(:,:)
 real(r8), allocatable :: grid_2d_array(:,:)
 integer  :: data_pattern
 real(r8) :: miss_percent
 logical  :: global, cyclic, polar
 real(r8) :: angle, lon_def, lat_def
end type src_testtype

! resulting sampling grid to be plotted - always fully regular
type smp_testtype
 integer  :: lon_count, lat_count
 real(r8) :: lon_start, lon_end
 real(r8) :: lat_start, lat_end
 real(r8), allocatable :: lon_1d_array(:), lat_1d_array(:)
 real(r8), allocatable :: data_2d_array(:,:)
 character(len=512) :: grid_lon_file
 character(len=512) :: grid_lat_file
 character(len=512) :: data_file
end type smp_testtype

integer :: debug = 0
type(src_testtype) :: dtt   ! source data
type(smp_testtype) :: stt   ! for now always regular, global?
type(quad_interp_handle) :: h
type(random_seq_type) :: ran
integer :: ntest

! needed?
integer  :: istatus, iunit, io

! variables changing the test; added to a namelist.

! interpolation test grid.  we construct a different grid
! and call the interpolation code on each corner of this
! other grid.  called 'sampling grid' to differentiate it
! from the 'data grid'.  usually much denser so we can look
! for discontinuties or errors in the interp code.

logical :: done

character(len=512) :: input_config_filename = "quad_test_input.txt"

!>@todo change this to include the input test filename (.txt)

namelist /test_quad_interp_nml/ input_config_filename, debug
  

! start of executable code

call initialize_utilities("test_quad_interp")
call init_random_seq(ran)

! Read the namelist entry
call find_namelist_in_file("input.nml", "test_quad_interp_nml", iunit)
read(iunit, nml = test_quad_interp_nml, iostat = io)
call check_namelist_read(iunit, io, "test_quad_interp_nml")

! Record the namelist values
if (do_nml_file()) write(nmlfileunit, nml=test_quad_interp_nml)
if (do_nml_term()) write(    *      , nml=test_quad_interp_nml)

call sample_grid_setup(stt)

! loop here over multiple tests
iunit = open_file(input_config_filename, action="read")

do 

   call data_grid_setup(dtt)

   call read_next_test(iunit, dtt, ntest, done)
   if (done) exit
   
   call fill_data_grid(dtt)
   
   call fill_sample_coords(stt)
   
   if (debug > 5) then
      call dump_src_grid(dtt)
      call dump_smp_grid(stt)
   endif
   
   call interpolate(dtt, stt)
   
   call write_original_grid(dtt)
   call write_sample_grid(stt, ntest)
   
   call takedown_data_grid(dtt)
enddo

! end test loop

call finalize_quad_interp(h)

if (debug > 0) print *, "closed files and finalized interp handle"

call finalize_utilities("test_quad_interp")


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

real(r8),              intent(in)    :: width
real(r8),              intent(in)    :: fraction
type(random_seq_type), intent(inout) :: seq
real(r8)                             :: deform

real(r8) :: val

! random val between -1 and 1
val = (random_uniform(seq) * 2.0_r8) - 1.0_r8

deform = val * width * fraction

end function deform

!------------------------------------------------------------

subroutine write_original_grid(tt)
type(src_testtype), intent(in) :: tt

! integer  :: lon_count, lat_count
! real(r8) :: lon_start, lon_end
! real(r8) :: lat_start, lat_end
! real(r8), allocatable :: lon_1d_array(:),   lat_1d_array(:)
! real(r8), allocatable :: lon_2d_array(:,:), lat_2d_array(:,:)
! real(r8), allocatable :: grid_2d_array(:,:)

select case(tt%grid_type)
 case (GRID_QUAD_FULLY_REGULAR)
   call writeit_1d("data_lons_1d_full_reg_test.txt", tt%lon_count, tt%lon_1d_array)
   call writeit_1d("data_lats_1d_full_reg_test.txt", tt%lat_count, tt%lat_1d_array)
   call writeit_2d("data_data_2d_full_reg_test.txt", tt%lon_count, tt%lat_count, tt%grid_2d_array)
 case (GRID_QUAD_IRREG_SPACED_REGULAR)
   call writeit_1d("data_lons_1d_part_reg_test.txt", tt%lon_count, tt%lon_1d_array)
   call writeit_1d("data_lats_1d_part_reg_test.txt", tt%lat_count, tt%lat_1d_array)
   call writeit_2d("data_data_2d_part_reg_test.txt", tt%lon_count, tt%lat_count, tt%grid_2d_array)
 case (GRID_QUAD_FULLY_IRREGULAR)
   call writeit_2d("data_lons_2d_irreg_test.txt", tt%lon_count, tt%lat_count, tt%lon_2d_array)
   call writeit_2d("data_lats_2d_irreg_test.txt", tt%lon_count, tt%lat_count, tt%lat_2d_array)
   call writeit_2d("data_data_2d_irreg_test.txt", tt%lon_count, tt%lat_count, tt%grid_2d_array)
 case default
   print *, "unexpected error: unrecognized grid_type in write_original_grid()"
   stop
end select

end subroutine write_original_grid

!------------------------------------------------------------

subroutine write_sample_grid(tt, tnum)
type(smp_testtype), intent(in) :: tt
integer,            intent(in) :: tnum

character(len=512) :: fname
character(len=16)  :: form

if (tnum < 10) then
  form='(A,I1,A)'
else
  form='(A,I2,A)'
endif

write(fname, form) trim(tt%grid_lon_file)//'_', tnum, '.txt'
call writeit_1d(fname, tt%lon_count, tt%lon_1d_array)

write(fname, form) trim(tt%grid_lat_file)//'_', tnum, '.txt'
call writeit_1d(fname, tt%lat_count, tt%lat_1d_array)

write(fname, form) trim(tt%data_file)//'_', tnum, '.txt'
call writeit_2d(fname, tt%lon_count, tt%lat_count, tt%data_2d_array)

end subroutine write_sample_grid

!------------------------------------------------------------
!> should we output directly with matlab formats?

subroutine writeit_1d(fname, nx, dataarray)
 character(len=*), intent(in) :: fname
 integer, intent(in) :: nx
 real(r8), intent(in) :: dataarray(:)

integer :: i, iunit

iunit = open_file(fname, action="write")

do i=1, nx
   write(iunit, *) dataarray(i)
enddo

call close_file(iunit)

end subroutine writeit_1d

!------------------------------------------------------------

subroutine writeit_2d(fname, nx, ny, dataarray)
 character(len=*), intent(in) :: fname
 integer, intent(in) :: nx, ny
 real(r8), intent(in) :: dataarray(nx, ny)

integer :: i, j, iunit

iunit = open_file(fname, action="write")

do j=1, ny
   do i=1, nx
      write(iunit, *) dataarray(i, j)
   enddo
enddo

call close_file(iunit)

end subroutine writeit_2d

!------------------------------------------------------------
!> for now, hardcode the sampling grid (global, regular)

subroutine sample_grid_setup(tt)
type(smp_testtype), intent(out) :: tt

! integer  :: lon_count, lat_count
! real(r8) :: lon_start, lon_end
! real(r8) :: lat_start, lat_end
! real(r8), allocatable :: data_2d_array(:,:)
! character(len=512) :: grid_lon_file, grid_lat_file
! character(len=512) :: data_file

tt%lon_count = 360
tt%lat_count = 180

tt%lon_start =   0.0_r8
tt%lon_end   = 360.0_r8
tt%lat_start = -90.0_r8
tt%lat_end   =  90.0_r8

allocate(tt%lon_1d_array(tt%lon_count), tt%lat_1d_array(tt%lat_count))
allocate(tt%data_2d_array(tt%lon_count, tt%lat_count))

tt%grid_lon_file = "sampled_lon"
tt%grid_lat_file = "sampled_lat"
tt%data_file     = "sampled_data"

end subroutine sample_grid_setup

!------------------------------------------------------------
! give all the items in the derived type initial values
! before reading anything in.

subroutine data_grid_setup(tt)
type(src_testtype), intent(inout) :: tt

tt%grid_type = -1
tt%relative_cell_location = -1
tt%lon_count = 0
tt%lat_count = 0
tt%lon_start = 0
tt%lon_end = 0
tt%lat_start = 0
tt%lat_end = 0
tt%data_pattern = -1
tt%miss_percent = 0.0_r8
tt%global = .false.
tt%cyclic = .false.
tt%polar = .false.
tt%angle = 0.0_r8
tt%lon_def = 0.0_r8
tt%lat_def = 0.0_r8

end subroutine data_grid_setup

!------------------------------------------------------------
!>@todo FIXME: too restrictive.  would like to do something
!>like:  lons min max count
!>on a single line.  use get_args_from_string() instead of
!>parsing a limited format.
!>
!> make each line contain at most 2 tokens.
!> read into a string buffer and parse from there.
!> # in column 1 is comment and entire line ignored
!> blank lines ignored
!> parse first word; interpretation of second word (string, int, real)
!>  depends on which word.

!   read(unitnum, "(A256)") line
!   call get_args_from_string(line, wordcount, words)


subroutine read_next_test(iunit, tt, ntest, done)
integer,            intent(in)    :: iunit
type(src_testtype), intent(inout) :: tt
integer,            intent(out)   :: ntest
logical,            intent(out)   :: done

character(len=512) :: next_line
character(len=128) :: words(25)  ! what's a reasonable limit here?
integer :: ios, ival, wordcount, i
real(r8) :: rval

done = .false.
do
   read(iunit, '(A512)', iostat=ios) next_line
   if (ios /= 0) then
      print *, 'error reading next line'
      exit
   endif
!print *, 'next line: ', trim(next_line)
   call get_args_from_string(next_line, wordcount, words)
!print *, 'wordcount = ', wordcount
!do i=1, wordcount
!  print *, ' word ', i, ' value "', trim(words(i)), '"'
!enddo

   if (wordcount == 0) cycle
   !if (next_line(1:1) == "#") cycle
   if (words(1)(1:1) == "#") cycle
   if (len_trim(next_line) == 0) cycle

   select case (words(1))
    case ('quad_test') 
       continue
    case ('version') 
       if (words(2) /= '1.0') print *, 'bad version number'
    case ('grid_type')
      select case (words(2))
       case ('GRID_QUAD_FULLY_REGULAR')
        tt%grid_type = GRID_QUAD_FULLY_REGULAR
       case ('GRID_QUAD_IRREG_SPACED_REGULAR')
        tt%grid_type = GRID_QUAD_IRREG_SPACED_REGULAR
       case ('GRID_QUAD_FULLY_IRREGULAR')
        tt%grid_type = GRID_QUAD_FULLY_IRREGULAR
       case default
        print *, 'unrecognized grid type'
      end select
    case ('data_loc') 
      select case (words(2))
       case('QUAD_LOCATED_CELL_CENTERS')
         tt%relative_cell_location = QUAD_LOCATED_CELL_CENTERS
       case('QUAD_LOCATED_LON_EDGES')
         tt%relative_cell_location = QUAD_LOCATED_LON_EDGES
       case('QUAD_LOCATED_LAT_EDGES')
         tt%relative_cell_location = QUAD_LOCATED_LAT_EDGES
       case('QUAD_LOCATED_CELL_CORNERS')
         tt%relative_cell_location = QUAD_LOCATED_CELL_CORNERS
       case default
        print *, 'unrecognized grid type'
      end select
    case ('lon_count') 
      read(words(2), *) ival
      tt%lon_count = ival
    case ('lat_count') 
      read(words(2), *) ival
      tt%lat_count = ival
    case ('lon_start') 
      read(words(2), *) rval
      tt%lon_start = rval
    case ('lon_end') 
      read(words(2), *) rval
      tt%lon_end = rval
    case ('lat_start') 
      read(words(2), *) rval
      tt%lat_start = rval
    case ('lat_end') 
      read(words(2), *) rval
      tt%lat_end = rval
    case ('data_pattern') 
      read(words(2), *) ival
      tt%data_pattern = ival
    case ('angle') 
      read(words(2), *) rval
      tt%angle = rval
    case ('lon_def') 
      read(words(2), *) rval
      tt%lon_def = rval
    case ('lat_def') 
      read(words(2), *) rval
      tt%lat_def = rval
    case ('start_test') 
      read(words(2), *) ntest
      print *, 'test number ', ntest
    case ('end_test')
      print *, 'end test ', ntest
      return
    case ('end_file')
      done = .true.
      exit
    case default
      print *, 'unrecognized token'
   end select

enddo

end subroutine read_next_test

!------------------------------------------------------------

subroutine fill_data_grid(tt)
type(src_testtype), intent(inout) :: tt

integer :: i, j
integer :: ndx, ndy
real(r8) :: del_lon, del_lat
real(r8) :: next_data, prev_data

! just allocate everything, may not be used but
! makes it easier to rotate a regular grid into
! an irregular one.
allocate(tt%lon_1d_array(tt%lon_count), tt%lat_1d_array(tt%lat_count))
allocate(tt%lon_2d_array(tt%lon_count, tt%lat_count))
allocate(tt%lat_2d_array(tt%lon_count, tt%lat_count))
allocate(tt%grid_2d_array(tt%lon_count, tt%lat_count))

ndx = tt%lon_count
ndy = tt%lat_count

del_lon = (tt%lon_end - tt%lon_start) / ndx
del_lat = (tt%lat_end - tt%lat_start) / ndy

do j=1, tt%lat_count
   tt%lat_1d_array(j) = tt%lat_start + (j-1)*del_lat
   do i=1, tt%lon_count
      tt%lon_1d_array(i) = tt%lon_start + (i-1)*del_lon

      tt%lon_2d_array(i, j) = tt%lon_1d_array(i)
      tt%lat_2d_array(i, j) = tt%lat_1d_array(j)

      if (tt%angle /= 0.0_r8) &
         call rotate(tt%lon_2d_array(i, j), tt%lat_2d_array(i, j), tt%angle, tt%lon_start, tt%lat_start)

      if (tt%lon_def /= 0.0_r8) &
         tt%lon_2d_array(i, j) = tt%lon_start + (i-1)*del_lon + deform(del_lon, tt%lon_def, ran)

      if (tt%lat_def /= 0.0_r8) &
         tt%lat_2d_array(i, j) = tt%lat_start + (j-1)*del_lat + deform(del_lat, tt%lat_def, ran)
   enddo
enddo

if (tt%grid_type /= GRID_QUAD_FULLY_IRREGULAR) then
   if (tt%angle /= 0.0_r8 .or. tt%lon_def /= 0.0_r8 .or. tt%lat_def /= 0.0_r8) &
      tt%grid_type = GRID_QUAD_FULLY_IRREGULAR
endif

print *, 'fill data grid: ', ndx, ndy
do j=1, tt%lat_count
   do i=1, tt%lon_count
      ! generate the data values on the corners.  pick one:
      select case (tt%data_pattern) 
      case (1)
         ! increasing monotonically 
         tt%grid_2d_array(i, j) = (j-1)*tt%lon_count + i
      case (2)
         ! constant by row
         tt%grid_2d_array(i, j) = j
      case (3)
         ! constant by column
         tt%grid_2d_array(i, j) = i
      case (4)
         ! based on lon only
         tt%grid_2d_array(i, j) = tt%lon_start + (del_lon * (i-1))
      case (5) 
         ! based on lat only
         tt%grid_2d_array(i, j) = tt%lat_start + (del_lat * (j-1))
      case (6)
         ! random between (0-10)
         tt%grid_2d_array(i, j) = random_uniform(ran) * 10.0_r8
      case (7)
         ! running average with gaussian noise added
         next_data = random_gaussian(ran, 0.0_r8, 1.0_r8)
         next_data = prev_data + (0.1 * next_data)
         tt%grid_2d_array(i, j) = next_data
         prev_data = next_data
      case default
         ! gaussian with mean 0 and stddev 1
         tt%grid_2d_array(i, j) = random_gaussian(ran, 0.0_r8, 1.0_r8)
      end select

      if (tt%miss_percent > 0.0_r8) then
        if (random_uniform(ran) * 100.0_r8 < tt%miss_percent) tt%grid_2d_array(i, j) = MISSING_R8 
      endif

   enddo
enddo

end subroutine fill_data_grid

!------------------------------------------------------------

subroutine fill_sample_coords(tt)
type(smp_testtype), intent(inout) :: tt

integer  :: i
real(r8) :: sample_del_lon, sample_del_lat

! "sampled grid" spacing along each axis
!do i=1, nsx
!   do j=1, nsy
!      sample_lons(i, j) = sample_start_lon + (i-1)*sample_del_lon + deform(sample_del_lon, lon_def, ran)
!      sample_lats(i, j) = sample_start_lat + (j-1)*sample_del_lat + deform(sample_del_lat, lat_def, ran)
!      ! generate locations of the corners of all the quads
!      if (angle /= 0.0_r8) &
!         call rotate(sample_lons(i, j), sample_lats(i, j), angle, sample_start_lon, sample_start_lat)
!    enddo
!enddo

sample_del_lon = (tt%lon_end - tt%lon_start) / tt%lon_count
sample_del_lat = (tt%lat_end - tt%lat_start) / tt%lat_count

do i=1, tt%lon_count
   tt%lon_1d_array(i) = tt%lon_start + (i-1)*sample_del_lon
enddo
do i=1, tt%lat_count
   tt%lat_1d_array(i) = tt%lat_start + (i-1)*sample_del_lat
enddo

end subroutine fill_sample_coords

!------------------------------------------------------------

subroutine interpolate(dtt, stt)
type(src_testtype), intent(in)    :: dtt
type(smp_testtype), intent(inout) :: stt

integer  :: i, j, k
integer  :: lon_indices(4), lat_indices(4)
real(r8) :: delta_lon, delta_lat
real(r8) :: invals(4)
real(r8) :: lon_fract, lat_fract
type(quad_interp_handle) :: h


delta_lon = (dtt%lon_end - dtt%lon_start) / dtt%lon_count
delta_lat = (dtt%lat_end - dtt%lat_start) / dtt%lat_count

stt%data_2d_array(:,:) = MISSING_R8

!subroutine init_quad_interp(grid_type, num_lons, num_lats, cell_relative, &
!                            global, spans_lon_zero, pole_wrap, interp_handle)


call init_quad_interp(dtt%grid_type, dtt%lon_count, dtt%lat_count, dtt%relative_cell_location, &
                      dtt%global, dtt%cyclic, dtt%polar, h)

if (dtt%grid_type == GRID_QUAD_FULLY_IRREGULAR) then
   call set_quad_coords(h, dtt%lon_2d_array, dtt%lat_2d_array)
else
! FIXME add semi-regular case here

   call set_quad_coords(h, dtt%lon_start, delta_lon, dtt%lat_start, delta_lat)
endif


do i=1, stt%lon_count
   do j=1, stt%lat_count

      ! find the quad corners
      if (dtt%grid_type == GRID_QUAD_FULLY_IRREGULAR) then
         call quad_lon_lat_locate(h, stt%lon_1d_array(i), stt%lat_1d_array(j), lon_indices, lat_indices, &
                                  istatus)
         if(debug > 0)print *, i, j, lon_indices, lat_indices, stt%lon_1d_array(i), stt%lat_1d_array(j)
      else
         call quad_lon_lat_locate(h, stt%lon_1d_array(i), stt%lat_1d_array(j), lon_indices, lat_indices, &
                                  lon_fract, lat_fract, istatus)
         if(debug > 0)print *, i, j, lon_indices, lat_indices, lon_fract, lat_fract, &
                               stt%lon_1d_array(i), stt%lat_1d_array(j)
  
      endif
      if (istatus /= 0) then
!print *, 'cannot locate ', sst%lon_1d_array(i), sst%lat_1d_array(j)
         cycle
      endif


      ! get values of data at lon/lat bot/top indices, counterclockwise around quad
      do k=1, 4
         invals(k) = dtt%grid_2d_array(lon_indices(k), lat_indices(k))
      enddo

      if (dtt%grid_type == GRID_QUAD_FULLY_IRREGULAR) then
         call quad_lon_lat_evaluate(h, stt%lon_1d_array(i), stt%lat_1d_array(j), lon_indices, lat_indices, &
                                    invals, stt%data_2d_array(i,j), istatus)
      else
         call quad_lon_lat_evaluate(h, lon_fract, lat_fract, invals, stt%data_2d_array(i,j), istatus)
      endif
      if (istatus /= 0) then
print *, 'cannot evaluate ', stt%lon_1d_array(i), stt%lat_1d_array(j)
         cycle
      endif

   enddo
enddo

end subroutine interpolate

!------------------------------------------------------------

subroutine takedown_data_grid(tt)
type(src_testtype), intent(inout) :: tt

if (allocated(tt%lon_1d_array)) deallocate(tt%lon_1d_array)
if (allocated(tt%lat_1d_array)) deallocate(tt%lat_1d_array)

if (allocated(tt%lon_2d_array)) deallocate(tt%lon_2d_array)
if (allocated(tt%lat_2d_array)) deallocate(tt%lat_2d_array)

deallocate(tt%grid_2d_array)

end subroutine takedown_data_grid

!------------------------------------------------------------

subroutine dump_src_grid(tt)
type(src_testtype), intent(in) :: tt

print *, 'data grid info:'
print *, 'grid type: ', tt%grid_type
print *, 'relative cell loc: ', tt%relative_cell_location
print *, 'lon count: ', tt%lon_count
print *, 'lon start/end: ', tt%lon_start, tt%lon_end
print *, 'lat count: ', tt%lat_count
print *, 'lat start/end: ', tt%lat_start, tt%lat_end
print *, 'allocated(lon_1d_array): ', allocated(tt%lon_1d_array)
print *, 'allocated(lat_1d_array): ', allocated(tt%lat_1d_array)
print *, 'allocated(lon_2d_array): ', allocated(tt%lon_2d_array)
print *, 'allocated(lat_2d_array): ', allocated(tt%lat_2d_array)
print *, 'allocated(grid_2d_array): ', allocated(tt%grid_2d_array)
print *, 'data pattern: ', tt%data_pattern
print *, 'miss percent: ', tt%miss_percent
print *, 'global, : ', tt%global
print *, 'cyclic, : ', tt%cyclic 
print *, 'polar: ', tt%polar
print *, 'angle : ', tt%angle 
print *, 'lon_def : ', tt%lon_def 
print *, 'lat_def: ', tt%lat_def

end subroutine dump_src_grid

!------------------------------------------------------------

subroutine dump_smp_grid(tt)
type(smp_testtype), intent(in) :: tt

print *, 'sampling grid info:'
print *, 'lon count: ', tt%lon_count
print *, 'lon start/end: ', tt%lon_start, tt%lon_end
print *, 'lat count: ', tt%lat_count
print *, 'lat start/end: ', tt%lat_start, tt%lat_end
print *, 'allocated(lon_1d_array): ', allocated(tt%lon_1d_array)
print *, 'allocated(lat_1d_array): ', allocated(tt%lat_1d_array)
print *, 'allocated(data_2d_array): ', allocated(tt%data_2d_array)
print *, 'lon file: ', trim(tt%grid_lon_file)
print *, 'lat file: ', trim(tt%grid_lat_file)
print *, 'data file: ', trim(tt%data_file)

end subroutine dump_smp_grid

!------------------------------------------------------------

end program test_quad_interp

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/roms_interpolation/developer_tests/interpolation/test_quad_interp.f90 $
! $Revision: 12852 $
! $Date: 2018-09-25 14:32:09 -0600 (Tue, 25 Sep 2018) $
