! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_gridded_sif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_gridded_sif - program that reads a sif netCDF profiler
!                          wind observation file and writes a DART
!                          obs_seq file using the DART library routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              check_namelist_read, find_namelist_in_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN, &
                              print_time, print_date
use      location_mod, only : VERTISHEIGHT
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : SOLAR_INDUCED_FLUORESCENCE, PARNORM_SIF
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname

use           netcdf

implicit none

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: ncid, nlon, nlat, ilon, ilat, i, oday, osec, nused, io, iunit
logical  :: file_exist, first_obs

real(r8) :: oerr, qc

real(r8), allocatable :: lat(:), lon(:), sif(:,:), normsif(:,:), sif_std(:,:), normsif_std(:,:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time

character(len=256) :: netcdf_file = 'sif_input.nc'
character(len=256) :: out_file    = 'obs_seq.out'

namelist /convert_gridded_sif_nml/ netcdf_file, out_file

!------------
! start of executable code
!------------

call initialize_utilities('convert_gridded_sif')

! Read the DART namelist for this model
call find_namelist_in_file('input.nml', 'convert_gridded_sif_nml', iunit)
read(iunit, nml = convert_gridded_sif_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_gridded_sif_nml')

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

first_obs = .true.

call nc_check( nf90_open(netcdf_file, nf90_nowrite, ncid), &
               'convert_gridded_sif', 'opening file "'//trim(netcdf_file)//'"')

call getdimlen(ncid, "longitude", nlon)
call getdimlen(ncid, "latitude" , nlat)

allocate(lat(nlat))
allocate(lon(nlon))

allocate(sif(        nlon,nlat))
allocate(sif_std(    nlon,nlat))
allocate(normsif(    nlon,nlat))
allocate(normsif_std(nlon,nlat))

! read in the data arrays

call getvar_real(ncid, 'latitude',   lat)
call getvar_real(ncid, 'longitude',  lon)
where(lon < 0.0_r8 )  lon = lon + 360.0_r8

call getvar_real_2d(ncid, "SIF_740",     sif    )
call getvar_real_2d(ncid, "SIF_740_std", sif_std)
call getvar_real_2d(ncid, "Par_normalized_SIF_740",     normsif)
call getvar_real_2d(ncid, "Par_normalized_SIF_740_std", normsif_std)

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=out_file, exist=file_exist)

if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(out_file, 0, 0, 2*nlon*nlat, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, 2*nlon*nlat)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

time_obs = get_time_from_filename(netcdf_file)
call get_time(time_obs, osec, oday)

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

nused = 0
latloop: do ilat = 1, nlat
   if ( lat(ilat) >  90.0_r8 .or. lat(ilat) <  -90.0_r8 ) cycle latloop

   lonloop: do ilon = 1, nlon
   
     ! check the lat/lon values to see if they are ok
     if ( lon(ilon) < 0.0_r8 .or. lon(ilon) > 360.0_r8 ) cycle lonloop
 
     !>@todo MARYIA ... how to get a land/sea mask ... 
     if (sif(ilon,ilat) <= 0.0_r8)  cycle lonloop

     if (sif(ilon,ilat) /= sif(ilon,ilat))  cycle lonloop  ! Test for NAN-ness
     if (sif_std(ilon,ilat) /= sif_std(ilon,ilat))  cycle lonloop  ! Test for NAN-ness

     oerr = sif_std(ilon,ilat)*2.0_r8  !>@todo MARYIA ... what about representativeness error
   
     call create_3d_obs(lat(ilat), lon(ilon), 0.0_r8, VERTISHEIGHT, sif(ilon,ilat), &
                        SOLAR_INDUCED_FLUORESCENCE, oerr, oday, osec, qc, obs)
     call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
     if (normsif(ilon,ilat) /= normsif(ilon,ilat))  cycle lonloop  ! Test for NAN-ness
     if (normsif_std(ilon,ilat) /= normsif_std(ilon,ilat))  cycle lonloop  ! Test for NAN-ness

     oerr = normsif_std(ilon,ilat)*2.0_r8  !>@todo MARYIA ... what about representativeness error

     call create_3d_obs(lat(ilat), lon(ilon), 0.0_r8, VERTISHEIGHT, normsif(ilon,ilat), &
                        PARNORM_SIF, oerr, oday, osec, qc, obs)
     call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
     nused = nused + 1
   
   enddo lonloop
enddo latloop

! need to wait to close file because in the loop it queries the
! report types.
call nc_check( nf90_close(ncid) , &
               'convert_gridded_sif', 'closing file '//trim(netcdf_file))

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, out_file)

! end of main program
call finalize_utilities()

contains

function get_time_from_filename(fname) result (mytime)
character(len=*), intent(in) :: fname
type(time_type) :: mytime

character(len=512) :: mystring
integer :: strlen, year,month,day1,dayN

strlen = len_trim(fname)

mystring = fname((strlen-13):strlen)

write(*,*)' I think the date is "'//trim(mystring)//'"'

! ret_f_nr5_nsvd12_v26_waves734_nolog.grid_SIF_gome2__0.5x0.5_20080101_07.nc
read(mystring,'(i4,i2,i2,1x,i2,3x)')year,month,day1,dayn

write(*,*)'year is ',year
write(*,*)'month is ',month
write(*,*)'day1 is ',day1
write(*,*)'dayN is ',dayN

mytime = set_date(year,month,day1)

call print_time(mytime,'just checking')
call print_date(mytime,'just checking')

end function get_time_from_filename

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
