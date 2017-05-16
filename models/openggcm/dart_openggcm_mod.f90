! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_openggcm_mod

use        types_mod, only : r8, rad2deg, PI, SECPERDAY
use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, get_calendar_string, &
                             print_date, print_time, operator(==), operator(-)
use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, nc_check, &
                             find_namelist_in_file, check_namelist_read, &
                             E_ERR, E_MSG, find_textfile_dims, &
                             logfileunit

use typesizes
use netcdf

implicit none
private

public :: get_openggcm_calendar, set_model_time_step, &
          get_horiz_grid_dims, get_vert_grid_dim, &
          read_horiz_grid, read_topography, read_vert_grid, &
          write_openggcm_namelist, get_openggcm_restart_filename

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring
logical, save :: module_initialized = .false.

character(len=256) :: ic_filename      = 'openggcm.r.nc'
!character(len=256) :: restart_filename = 'dart_openggcm_mod_restart_filename_not_set'

! set this to true if you want to print out the current time
! after each N observations are processed, for benchmarking.
logical :: print_timestamps = .false.
integer :: print_every_Nth  = 10000

!------------------------------------------------------------------
! The openggcm time manager namelist variables
!------------------------------------------------------------------

character(len=100) :: accel_file ! length consistent with openggcm
character(len= 64) :: stop_option, runid, dt_option, time_mix_opt
character(len=  1) :: date_separator
logical  :: impcor, laccel, allow_leapyear
real(r8) :: dtuxcel, dt_count
integer  :: iyear0, imonth0, iday0, ihour0, iminute0, isecond0
integer  :: stop_count, fit_freq, time_mix_freq

namelist /time_manager_nml/ runid, time_mix_opt, time_mix_freq, &
    impcor, laccel, accel_file, dtuxcel, iyear0, imonth0, &
    iday0, ihour0, iminute0, isecond0, dt_option, dt_count, &
    stop_option, stop_count, date_separator, allow_leapyear, fit_freq

!------------------------------------------------------------------
! The openggcm I/O namelist variables
!------------------------------------------------------------------

character(len=100) :: log_filename, pointer_filename
logical :: lredirect_stdout, luse_pointer_files
logical :: luse_nf_64bit_offset
integer :: num_iotasks

namelist /io_nml/ num_iotasks, lredirect_stdout, log_filename, &
    luse_pointer_files, pointer_filename, luse_nf_64bit_offset

!------------------------------------------------------------------
! The openggcm restart manager namelist variables
!------------------------------------------------------------------

character(len=100) :: restart_outfile ! length consistent with openggcm
character(len= 64) :: restart_freq_opt, restart_start_opt, restart_fmt
logical :: leven_odd_on, pressure_correction
integer :: restart_freq, restart_start, even_odd_freq

namelist /restart_nml/ restart_freq_opt, restart_freq, &
     restart_start_opt, restart_start, restart_outfile, &
    restart_fmt, leven_odd_on, even_odd_freq, pressure_correction

!------------------------------------------------------------------
! The openggcm initial temperature and salinity namelist 
!------------------------------------------------------------------

character(len=100) :: init_ts_file ! length consistent with openggcm
character(len=100) :: init_ts_outfile
character(len= 64) :: init_ts_option, init_ts_suboption
character(len= 64) :: init_ts_file_fmt, init_ts_outfile_fmt

namelist /init_ts_nml/ init_ts_option, init_ts_suboption, &
                       init_ts_file, init_ts_file_fmt, & 
                       init_ts_outfile, init_ts_outfile_fmt

!------------------------------------------------------------------
! The openggcm domain namelist 
!------------------------------------------------------------------

character(len= 64) :: clinic_distribution_type, tropic_distribution_type
character(len= 64) :: ew_boundary_type, ns_boundary_type
integer :: nprocs_clinic, nprocs_tropic

namelist /domain_nml/ clinic_distribution_type, nprocs_clinic, &
                      tropic_distribution_type, nprocs_tropic, &
                      ew_boundary_type, ns_boundary_type

!------------------------------------------------------------------
! The openggcm grid info namelist 
!------------------------------------------------------------------
!
! openggcm grid information comes in several files:
!   horizontal grid lat/lons in one, 
!   topography (lowest valid vert level) in another, and 
!   the vertical grid spacing in a third.
!
!------------------------------------------------------------------
!
! Here is what we can get from the (binary) horiz grid file:
!   real (r8), dimension(:,:), allocatable :: &
!      ULAT,            &! latitude  (radians) of U points
!      ULON,            &! longitude (radians) of U points
!      HTN ,            &! length (cm) of north edge of T box
!      HTE ,            &! length (cm) of east  edge of T box
!      HUS ,            &! length (cm) of south edge of U box
!      HUW ,            &! length (cm) of west  edge of U box
!      ANGLE             ! angle
!
! Here is what we can get from the topography file:
!   integer, dimension(:,:), allocatable :: &
!      KMT               ! k index of deepest grid cell on T grid
!
! These must be derived or come from someplace else ...
!      KMU               ! k index of deepest grid cell on U grid
!      HT                ! real(r8) value of deepest valid T depth (in cm)
!      HU                ! real(r8) value of deepest valid U depth (in cm)
!
! The vert grid file is ascii, with 3 columns/line:
!    cell thickness(in cm)   cell center(in m)   cell bottom(in m)
!
!------------------------------------------------------------------

character(len=100) :: horiz_grid_file, vert_grid_file, &
                      topography_file, topography_outfile, &
                      bathymetry_file, region_info_file, &
                      bottom_cell_file, region_mask_file
character(len= 64) :: horiz_grid_opt, sfc_layer_opt, vert_grid_opt, &
                      topography_opt
logical :: partial_bottom_cells, topo_smooth, flat_bottom, lremove_points
integer :: kmt_kmin, n_topo_smooth

namelist /grid_nml/ horiz_grid_opt, horiz_grid_file, sfc_layer_opt, &
    vert_grid_opt, vert_grid_file, topography_opt, kmt_kmin, &
    topography_file, topography_outfile, bathymetry_file, &
    partial_bottom_cells, bottom_cell_file, n_topo_smooth, &
    region_mask_file, topo_smooth, flat_bottom, lremove_points, &
    region_info_file

!======================================================================
contains
!======================================================================


subroutine initialize_module
!------------------------------------------------------------------
integer :: iunit, io

! Read openggcm calendar information
! In 'restart' mode, this is primarily the calendar type and 'stop'
! information. The time attributes of the restart file override 
! the namelist time information. 

call find_namelist_in_file('openggcm_in', 'time_manager_nml', iunit)
read(iunit, nml = time_manager_nml, iostat = io)
call check_namelist_read(iunit, io, 'time_manager_nml')

! FIXME : Real observations are always GREGORIAN dates ...
! but stomping on that here gets in the way of running
! a perfect_model experiment for pre-1601 AD cases.

! STOMP if ( allow_leapyear ) then
   call set_calendar_type('gregorian')
! STOMP else
! STOMP    call set_calendar_type('noleap')
! STOMP endif

! Read openggcm I/O information (for restart file ... grid dimensions)
! Read openggcm initial information (for input/restart filename)

call find_namelist_in_file('openggcm_in', 'io_nml', iunit)
read(iunit, nml = io_nml, iostat = io)
call check_namelist_read(iunit, io, 'io_nml')

call find_namelist_in_file('openggcm_in', 'init_ts_nml', iunit)
read(iunit, nml = init_ts_nml, iostat = io)
call check_namelist_read(iunit, io, 'init_ts_nml')

! Is it a pointer file or not ...
!if ( luse_pointer_files ) then
!
!   restart_filename = trim(pointer_filename)//'.restart'
!
!   if ( .not. file_exist(restart_filename) ) then
!      msgstring = 'openggcm_in:pointer file '//trim(restart_filename)//' not found'
!      call error_handler(E_ERR,'initialize_module', &
!             msgstring, source, revision, revdate)
!   endif
!
!   iunit = open_file(restart_filename,'formatted')
!   read(iunit,'(A)')ic_filename
!
!   restart_filename = ' '  
!   write(*,*)'DEBUG ... pointer filename dereferenced to ',trim(ic_filename )
!
!else
!   ic_filename = trim(init_ts_file)//'.'//trim(init_ts_file_fmt)
!endif

! Make sure we have a openggcm restart file (for grid dims)
if ( .not. file_exist(ic_filename) ) then
   msgstring = 'openggcm_in:init_ts_file '//trim(ic_filename)//' not found'
   call error_handler(E_ERR,'initialize_module', &
          msgstring, source, revision, revdate)
endif

! Read openggcm restart information (for model timestepping)
call find_namelist_in_file('openggcm_in', 'restart_nml', iunit)
read(iunit, nml = restart_nml, iostat = io)
call check_namelist_read(iunit, io, 'restart_nml')

! Read openggcm domain information (for lon wrapping or not)
call find_namelist_in_file('openggcm_in', 'domain_nml', iunit)
read(iunit, nml = domain_nml, iostat = io)
call check_namelist_read(iunit, io, 'domain_nml')

! Read openggcm grid information (for grid filenames)
call find_namelist_in_file('openggcm_in', 'grid_nml', iunit)
read(iunit, nml = grid_nml, iostat = io)
call check_namelist_read(iunit, io, 'grid_nml')

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

end subroutine initialize_module



subroutine get_horiz_grid_dims(Nx, Ny)
!------------------------------------------------------------------
! subroutine get_horiz_grid_dims(Nx, Ny)
!
! Read the lon, lat grid size from the restart netcdf file.
! The actual grid file is a binary file with no header information.
!
! The file name comes from module storage ... namelist.

integer, intent(out) :: Nx   ! Number of Longitudes
integer, intent(out) :: Ny   ! Number of Latitudes

integer :: grid_id, dimid, nc_rc

if ( .not. module_initialized ) call initialize_module

! get the ball rolling ...

call nc_check(nf90_open(trim(ic_filename), nf90_nowrite, grid_id), &
         'get_horiz_grid_dims','open '//trim(ic_filename))

! Longitudes : get dimid for 'i' or 'nlon', and then get value
nc_rc = nf90_inq_dimid(grid_id, 'i', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlon', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'i' or 'nlon' in file "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate) 
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Nx), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(ic_filename))

! Latitudes : get dimid for 'j' or 'nlat' ... and then get value
nc_rc = nf90_inq_dimid(grid_id, 'j', dimid)
if (nc_rc /= nf90_noerr) then
   nc_rc = nf90_inq_dimid(grid_id, 'nlat', dimid)
   if (nc_rc /= nf90_noerr) then
      msgstring = "unable to find either 'j' or 'nlat' in "//trim(ic_filename)
      call error_handler(E_ERR, 'get_horiz_grid_dims', msgstring, &
                         source,revision,revdate)
   endif
endif

call nc_check(nf90_inquire_dimension(grid_id, dimid, len=Ny), &
         'get_horiz_grid_dims','inquire_dimension i '//trim(ic_filename))

! tidy up

call nc_check(nf90_close(grid_id), &
         'get_horiz_grid_dims','close '//trim(ic_filename) )

end subroutine get_horiz_grid_dims



  subroutine get_vert_grid_dim(Nz)
!------------------------------------------------------------------
! subroutine get_vert_grid_dim(Nz)
!
! count the number of lines in the ascii file to figure out max
! number of vert blocks.

integer, intent(out) :: Nz

integer :: linelen ! disposable

if ( .not. module_initialized ) call initialize_module

call find_textfile_dims(vert_grid_file, Nz, linelen)

end subroutine get_vert_grid_dim


   
subroutine get_openggcm_calendar(calstring)
!------------------------------------------------------------------
! the initialize_module ensures that the openggcm namelists are read and 
! the DART time manager gets the openggcm calendar setting.
!
! Then, the DART time manager is queried to return what it knows ...
!
character(len=*), INTENT(OUT) :: calstring

if ( .not. module_initialized ) call initialize_module

call get_calendar_string(calstring)

end subroutine get_openggcm_calendar



function set_model_time_step()
!------------------------------------------------------------------
! the initialize_module ensures that the openggcm namelists are read.
! The restart times in the openggcm_in&restart_nml are used to define
! appropriate assimilation timesteps.
!
type(time_type) :: set_model_time_step

if ( .not. module_initialized ) call initialize_module

! Check the 'restart_freq_opt' and 'restart_freq' to determine
! when we can stop the model

if ( trim(restart_freq_opt) == 'nday' ) then
   set_model_time_step = set_time(0, restart_freq) ! (seconds, days)
else if ( trim(restart_freq_opt) == 'nyear' ) then
   ! FIXME ... CCSM_openggcm uses a bogus value for this
   set_model_time_step = set_time(0, 1) ! (seconds, days)
else
   call error_handler(E_ERR,'set_model_time_step', &
              'restart_freq_opt must be nday', source, revision, revdate)
endif

end function set_model_time_step




subroutine write_openggcm_namelist(model_time, adv_to_time)
!------------------------------------------------------------------
!
type(time_type), INTENT(IN) :: model_time, adv_to_time
type(time_type) :: offset

integer :: iunit, secs, days

if ( .not. module_initialized ) call initialize_module

offset = adv_to_time - model_time
call get_time(offset, secs, days)

if (secs /= 0 ) then
   write(msgstring,*)'adv_to_time has seconds == ',secs,' must be zero'
   call error_handler(E_ERR,'write_openggcm_namelist', msgstring, source, revision, revdate)
endif

! call print_date( model_time,'write_openggcm_namelist:dart model date')
! call print_date(adv_to_time,'write_openggcm_namelist:advance_to date')
! call print_time( model_time,'write_openggcm_namelist:dart model time')
! call print_time(adv_to_time,'write_openggcm_namelist:advance_to time')
! call print_time(     offset,'write_openggcm_namelist:a distance of')
! write( *,'(''write_openggcm_namelist:TIME_MANAGER_NML   STOP_COUNT '',i10,'' days'')') days

!Convey the information to the namelist 'stop option' and 'stop count'

if ( trim(stop_option) == 'nday' ) then
   stop_count = days
else
   call error_handler(E_ERR,'write_openggcm_namelist', &
              'stop_option must be "nday"', source, revision, revdate)
endif

iunit = open_file('openggcm_in.DART',form='formatted',action='rewind')
write(iunit, nml=time_manager_nml)
write(iunit, '('' '')')
close(iunit)

end subroutine write_openggcm_namelist



  subroutine read_horiz_grid(nx, ny, ULAT, ULON, TLAT, TLON)
!------------------------------------------------------------------
! subroutine read_horiz_grid(nx, ny, ULAT, ULON, TLAT, TLON)
!
! Open and read the binary grid file

integer,                    intent(in)  :: nx, ny
real(r8), dimension(nx,ny), intent(out) :: ULAT, ULON, TLAT, TLON

!real(r8), dimension(nx,ny) :: &
!     HTN ,  &! length (cm) of north edge of T box
!     HTE ,  &! length (cm) of east  edge of T box
!     HUS ,  &! length (cm) of south edge of U box
!     HUW ,  &! length (cm) of west  edge of U box
!     ANGLE  ! angle

integer :: grid_unit, reclength

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(horiz_grid_file) ) then
   msgstring = 'openggcm_in:horiz_grid_file '//trim(horiz_grid_file)//' not found'
   call error_handler(E_ERR,'read_horiz_grid', &
          msgstring, source, revision, revdate)
endif

! Open it and read them in the EXPECTED order.
! Actually, we only need the first two, so I'm skipping the rest.

grid_unit = get_unit()
INQUIRE(iolength=reclength) ULAT

open(grid_unit, file=trim(horiz_grid_file), form='unformatted', &
            access='direct', recl=reclength, status='old', action='read' )
read(grid_unit, rec=1) ULAT
read(grid_unit, rec=2) ULON
!read(grid_unit, rec=3) HTN
!read(grid_unit, rec=4) HTE
!read(grid_unit, rec=5) HUS
!read(grid_unit, rec=6) HUW
!read(grid_unit, rec=7) ANGLE
close(grid_unit)

call calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)

! convert from radians to degrees

ULAT = ULAT * rad2deg
ULON = ULON * rad2deg
TLAT = TLAT * rad2deg
TLON = TLON * rad2deg

! ensure [0,360) [-90,90]

where (ULON <   0.0_r8) ULON = ULON + 360.0_r8
where (ULON > 360.0_r8) ULON = ULON - 360.0_r8
where (TLON <   0.0_r8) TLON = TLON + 360.0_r8
where (TLON > 360.0_r8) TLON = TLON - 360.0_r8

where (ULAT < -90.0_r8) ULAT = -90.0_r8
where (ULAT >  90.0_r8) ULAT =  90.0_r8
where (TLAT < -90.0_r8) TLAT = -90.0_r8
where (TLAT >  90.0_r8) TLAT =  90.0_r8

end subroutine read_horiz_grid


  subroutine calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)
!------------------------------------------------------------------
! subroutine calc_tpoints(nx, ny, ULAT, ULON, TLAT, TLON)
!
! mimic openggcm grid.F90:calc_tpoints(), but for one big block.

integer,                    intent( in) :: nx, ny
real(r8), dimension(nx,ny), intent( in) :: ULAT, ULON
real(r8), dimension(nx,ny), intent(out) :: TLAT, TLON

integer  :: i, j
real(r8) :: xc,yc,zc,xs,ys,zs,xw,yw,zw   ! Cartesian coordinates for
real(r8) :: xsw,ysw,zsw,tx,ty,tz,da      ! nbr points

real(r8), parameter ::  c0 = 0.000_r8, c1 = 1.000_r8
real(r8), parameter ::  c2 = 2.000_r8, c4 = 4.000_r8
real(r8), parameter :: p25 = 0.250_r8, p5 = 0.500_r8
real(r8)            :: pi, pi2, pih, radian

if ( .not. module_initialized ) call initialize_module

! Define some constants as in openggcm

pi     = c4*atan(c1)
pi2    = c2*pi
pih    = p5*pi
radian = 180.0_r8/pi

! initialize these arrays to 0. in the code below there is
! a column that is referenced by a where() construct before 
! those values are set.  make sure that it doesn't cause a
! floating point exception from random memory bits which aren't
! valid floating point numbers.
TLAT(:,:) = c0
TLON(:,:) = c0

do j=2,ny
do i=2,nx

   !*** convert neighbor U-cell coordinates to 3-d Cartesian coordinates 
   !*** to prevent problems with averaging near the pole

   zsw = cos(ULAT(i-1,j-1))
   xsw = cos(ULON(i-1,j-1))*zsw
   ysw = sin(ULON(i-1,j-1))*zsw
   zsw = sin(ULAT(i-1,j-1))

   zs  = cos(ULAT(i  ,j-1))
   xs  = cos(ULON(i  ,j-1))*zs
   ys  = sin(ULON(i  ,j-1))*zs
   zs  = sin(ULAT(i  ,j-1))

   zw  = cos(ULAT(i-1,j  ))
   xw  = cos(ULON(i-1,j  ))*zw
   yw  = sin(ULON(i-1,j  ))*zw
   zw  = sin(ULAT(i-1,j  ))

   zc  = cos(ULAT(i  ,j  ))
   xc  = cos(ULON(i  ,j  ))*zc
   yc  = sin(ULON(i  ,j  ))*zc
   zc  = sin(ULAT(i  ,j  ))

   !*** straight 4-point average to T-cell Cartesian coords

   tx = p25*(xc + xs + xw + xsw)
   ty = p25*(yc + ys + yw + ysw)
   tz = p25*(zc + zs + zw + zsw)

   !*** convert to lat/lon in radians

   da = sqrt(tx**2 + ty**2 + tz**2)

   TLAT(i,j) = asin(tz/da)

   if (tx /= c0 .or. ty /= c0) then
      TLON(i,j) = atan2(ty,tx)
   else
      TLON(i,j) = c0
   endif

end do
end do

!*** for bottom row of domain where sw 4pt average is not valid,
!*** extrapolate from interior
!*** NOTE: THIS ASSUMES A CLOSED SOUTH BOUNDARY - WILL NOT
!***       WORK CORRECTLY FOR CYCLIC OPTION

do i=1,nx
   TLON(i,1) =    TLON(i,1+1)
   TLAT(i,1) = c2*TLAT(i,1+1) - TLAT(i,1+2)
end do

where (TLON(:,:) > pi2) TLON(:,:) = TLON(:,:) - pi2
where (TLON(:,:) < c0 ) TLON(:,:) = TLON(:,:) + pi2

!*** this leaves the leftmost/western edge to be filled 
!*** if the longitudes wrap, this is easy.
!*** the gx3v5 grid TLON(:,2) and TLON(:,nx) are both about 2pi,
!*** so taking the average is reasonable.
!*** averaging the latitudes is always reasonable.

if ( trim(ew_boundary_type) == 'cyclic' ) then

   TLAT(1,:) = (TLAT(2,:) + TLAT(nx,:))/c2
   TLON(1,:) = (TLON(2,:) + TLON(nx,:))/c2

else
   write(msgstring,'(''openggcm_in&domain_nml:ew_boundary_type '',a,'' unknown.'')') &
                                    trim(ew_boundary_type)
   call error_handler(E_ERR,'calc_tpoints',msgstring,source,revision,revdate)
endif

end subroutine calc_tpoints



  subroutine read_topography(nx, ny, KMT, KMU)
!------------------------------------------------------------------
! subroutine read_topography(nx, ny, KMT, KMU)
!
! Open and read the binary topography file

integer,                   intent(in)  :: nx, ny
integer, dimension(nx,ny), intent(out) :: KMT, KMU

integer  :: i, j, topo_unit, reclength

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(topography_file) ) then
   msgstring = 'openggcm_in:topography_file '//trim(topography_file)//' not found'
   call error_handler(E_ERR,'read_topography', &
          msgstring, source, revision, revdate)
endif

! read the binary file

topo_unit = get_unit()
INQUIRE(iolength=reclength) KMT

open( topo_unit, file=trim(topography_file), form='unformatted', &
            access='direct', recl=reclength, status='old', action='read' )
read( topo_unit, rec=1) KMT
close(topo_unit)

! the equation numbered 3.2 on page 15 of this document:
!  http://www.cesm.ucar.edu/models/cesm1.0/openggcm2/doc/sci/openggcmRefManual.pdf
! is WRONG.  (WRONG == inconsistent with the current openggcm source code.)
!
! for any U(i,j), the T(i,j) point with the same index values is located
! south and west. so the T points which surround any U(i,j) point are
! in fact at indices i,i+1, and j,j+1 .
!
!  NO: KMU(i,j) = min(KMT(i, j), KMT(i-1, j), KMT(i, j-1), KMT(i-1, j-1)) 
! YES: KMU(i,j) = min(KMT(i, j), KMT(i+1, j), KMT(i, j+1), KMT(i+1, j+1))
!
! the latter matches the openggcm source code, on yellowstone, lines 908 and 909 in:
!  /glade/p/cesm/releases/cesm1_2_2/models/ocn/openggcm2/source/grid.F90
!
! wrap around longitude boundary at i == nx.  set the topmost (last) latitude
! U row to the same value in all cases. in the shifted pole grid currently in 
! use all these points are on land and so are 0.  in the original unshifted
! lat/lon grid these last row U points are above the final T row and are believed
! to be unused.  for completeness we set all values in the last U row to the 
! minimum of the all T row values immediately below it, for all longitudes.

do j=1,ny-1
   do i=1,nx-1
      KMU(i,j) = min(KMT(i, j), KMT(i+1, j), KMT(i, j+1), KMT(i+1, j+1))
   enddo
   KMU(nx,j) = min(KMT(nx, j), KMT(1, j), KMT(nx, j+1), KMT(1, j+1))
enddo
KMU(:,ny) = minval(KMT(:,ny))

end subroutine read_topography



  subroutine read_vert_grid(nz, ZC, ZG)
!------------------------------------------------------------------
! subroutine read_vert_grid(nz, ZC, ZG)
!
! Open and read the ASCII vertical grid information
!
! The vert grid file is in ascii, with either 3 columns/line
!    cell thickness(in cm)   cell center(in m)   cell bottom(in m)
! or it can contain 2 columns/line
!    cell thickness(in cm)   1.0
! in which case we compute the cell center (ZC) and cell bottom (ZG)
! and ignore the second column.

integer,  intent(in)  :: nz
real(r8), intent(out) :: ZC(nz), ZG(nz)

integer  :: iunit, i, ios
real(r8) :: depth

logical :: three_columns
character(len=256)  :: line

real(r8), parameter :: centemeters_to_meters = 0.01_r8

if ( .not. module_initialized ) call initialize_module

! Check to see that the file exists.

if ( .not. file_exist(vert_grid_file) ) then
   msgstring = 'openggcm_in:vert_grid_file '//trim(vert_grid_file)//' not found'
   call error_handler(E_ERR,'read_vert_grid', &
          msgstring, source, revision, revdate)
endif

! read the ASCII file

iunit = open_file(trim(vert_grid_file), action = 'read')

! determine the number of columns
read(iunit,'(A)') line

read(line,*,iostat=ios) depth, ZC(1), ZG(1)

if(ios == 0) then
   three_columns = .true.
else
   three_columns = .false.
   
   ! read depth and calculate center and bottom of cells
   read(line,*,iostat=ios) depth 

   ZC(1) = depth*centemeters_to_meters*0.5 
   ZG(1) = depth*centemeters_to_meters
endif

do i=2, nz
   
   if(three_columns) then
      read(iunit,*,iostat=ios) depth, ZC(i), ZG(i)
   else
      read(iunit,*,iostat=ios) depth

      ZC(i) = ZG(i-1) + depth*centemeters_to_meters*0.5
      ZG(i) = ZG(i-1) + depth*centemeters_to_meters
   endif

   if ( ios /= 0 ) then ! error
      write(msgstring,*)'error reading depths, line ',i
      call error_handler(E_ERR,'read_vert_grid',msgstring,source,revision,revdate)
   endif

enddo

end subroutine read_vert_grid


!------------------------------------------------------------------


subroutine get_openggcm_restart_filename( filename )
character(len=*), intent(OUT) :: filename

if ( .not. module_initialized ) call initialize_module

filename   = trim(ic_filename)

end subroutine get_openggcm_restart_filename


end module dart_openggcm_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
