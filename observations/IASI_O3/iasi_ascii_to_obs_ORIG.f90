! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program iasi_ascii_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   iasi_ascii_to_obs - a program that only needs minor customization to read
!      in a iasi_ascii-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     this is work in progress. IASI dataset are in HDF format. I do not
!     have HDF libraries for now, so Gabi Pfister reads the hdf file in IDL and
!     did some processing before she dumped the data in ascii. what you are
!     reading here is a 'processed' dataset of IASI O3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r4, r8, digits12
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), GREGORIAN, operator(+), print_date
use      location_mod, only : VERTISUNDEF, VERTISHEIGHT,VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data
use      obs_kind_mod, only : IASI_O3_RETRIEVAL
use obs_utilities_mod, only : add_obs_to_seq

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, parameter :: debug = .true.  ! set to .true. to print info

character(len=64), parameter :: iasi_ascii_input_file = 'iasi_asciidata.input'
character(len=64), parameter :: obs_out_file          = 'iasi_obs_seq.out'
character(len=84) :: input_line
character(len=10) :: otype_char

logical  :: first_obs, file_exist

integer, parameter :: nlevels = 40

integer :: i, oday, osec, rcio, iunit, ilev
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs

integer  :: qc_count
real(r8) :: iasi_col, iasi_vmr, iasi_err, sza, cloud, dfs100, dfs300
real(r8) :: lat, lon, retlev2, qc, rnlev_use
real(r8) :: altretlev2(nlevels) = 0.0_r8
real(r8) :: akcol(nlevels) = 0.0_r8
real(r8) :: apcol(nlevels) = 0.0_r8
real(r8) :: iasi_col_prof(nlevels) = 0.0_r8
real(r8) :: iasi_vmr_prof(nlevels) = 0.0_r8
real(r8) :: aircol_val(nlevels) = 0.0_r8
real(r8) :: apcol_val
real(r8) :: seconds
integer  :: sellev, nlev_use
!real(r8) :: aircol_val

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

! Start of executable code
call initialize_utilities('iasi_ascii_to_obs')

! Time setup
call set_calendar_type(GREGORIAN)

!! some times are supplied as number of seconds since some reference
!! date.  This is an example of how to support that.
!! put the reference date into DART format
!comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

! open input iasi_ascii file

iunit = open_file(iasi_ascii_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(iasi_ascii_input_file)

! Each observation in this series will have a single observation value
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
max_obs    = 10000000
num_copies = 1
num_qc     = 1

! Call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! Create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! The first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! if you want to append to existing files (e.g. you have a lot of
! small iasi_ascii files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files
! once they are in DART obs_seq format.

!  ! existing file found, append to it
!  inquire(file=obs_out_file, exist=file_exist)
!  if ( file_exist ) then
!     call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)
!  endif

! Set the DART data quality control.   0 is good data.
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8
qc_count = 0

! Loop for reading IASI ascii data file
obsloop: do    ! No end limit - have the loop break when input ends

   ! Read in a line from the iasi_ascii file.   What you need to create an obs:
   !  location: lat, lon, and height in pressure or meters
   !  time: when the observation was taken
   !  type: from the DART list of obs types
   !  error: very important - the instrument error plus representativeness error
   !        (see html file for more info)
   !  averaging kernel and a priori profile
   !  for now, we chose 2 'retrieval levels' corresponding to highest sensitivity
   !  assume to be independent from each other
   !  assume here a line is a type (1/2), location, time, value, obs error

   ! read in entire iasi_ascii line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   !print *, input_line
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code from input file, rcio = ', rcio
      exit obsloop
   endif

   read(input_line(1:10), *, iostat=rcio) otype_char
   if (rcio /= 0) then
      if (debug) print *, 'got bad read code trying to get obs type, rcio = ', rcio
      exit obsloop
   endif

   ! otype is fixed to 1 for IASI O3
   if (debug) print *, 'next observation type = ', otype_char
   read(input_line(12:15), *, iostat=rcio) year
   read(input_line(16:17), *, iostat=rcio) month
   read(input_line(18:19), *, iostat=rcio) day
   ! next line
   read(iunit,"(9F14.4)", iostat=rcio) seconds, lat, lon, sza, cloud, retlev2, &
                            rnlev_use, dfs100, dfs300
   !write(*,"(9F14.4)", iostat=rcio) seconds, lat, lon, sza, cloud, retlev2, &
   !                         rnlev_use, dfs100, dfs300

   nlev_use = int(rnlev_use)
   retlev2  = retlev2*1000.0_r8


   hour =  seconds/3600
   minute = (seconds - hour*3600)/60
   second = (seconds - hour*3600 - minute*60)
   ! next line
   read(iunit, *, iostat=rcio) (altretlev2(i),i=1,nlev_use)

   altretlev2(1:nlev_use) = altretlev2(1:nlev_use)*1000.0_r8

   ! find sellev
   sellev = 0
   do ilev = 1,nlev_use
      if ( altretlev2(ilev) <= retlev2 ) then
           sellev = ilev
      endif
   enddo

   ! next line
   read(iunit,*, iostat=rcio) iasi_col, iasi_vmr
   !write(*, iostat=rcio) iasi_col, iasi_vmr

   ! next line
   read(iunit,*, iostat=rcio) iasi_err
   iasi_err = iasi_err*iasi_col

   ! next line
   read(iunit,*, iostat=rcio) (akcol(i),i=1,nlev_use)

   ! next line
   read(iunit,*, iostat=rcio) (apcol(i),i=1,nlev_use)

   apcol_val = 0.0_r8
   ! transform to (I-A)xa
   do ilev = 1, nlev_use
        apcol_val = apcol_val + akcol(ilev)*apcol(ilev)
   enddo
   apcol_val = apcol(sellev) - apcol_val

   read(iunit,*, iostat=rcio) (iasi_col_prof(i), i=1,nlev_use)
   !print *, iasi_col_prof
   read(iunit,*, iostat=rcio) (iasi_vmr_prof(i), i=1,nlev_use)
   !print *, iasi_vmr_prof
   !print *, year, month, day, hour, minute, second
   !print *, seconds, lat, lon, sza, cloud, retlev2, nlev_use, dfs100, dfs300
   !print *, 'alt'
   !print *, altretlev2
   !print *, 'col vmr'
   !print *, iasi_col, iasi_vmr
   !print *, 'err'
   !print *, iasi_err
   !print *, 'akcol'
   !print *, akcol
   !print *, 'apcol'
   !print *, apcol
   !print *, 'debug'
   !exit
   !endif
   !if (nlev_use < nlevels ) then
   !   cycle obsloop
   !else
   !if (iasi_vmr > 0.0_r8) then
   do ilev = 1, nlevels
      !print *, iasi_vmr_prof(i), iasi_col_prof(i)
      if (iasi_vmr_prof(ilev)>0.0_r8) then
      aircol_val(ilev) = iasi_col_prof(ilev)/iasi_vmr_prof(ilev)
      else
      aircol_val(ilev) = 0.0_r8
      endif
   enddo
   !else
   !   print *,'AFAJ DEBUG '
   !   exit
   !endif

    ! temp correction
    !day = day + 1
    ! play with the error for now
    if (debug) print *, 'next observation located at lat, lon = ', lat, lon
    if (rcio /= 0) then
       if (debug) print *, 'got bad read code getting rest of iasi o3 obs, rcio = ', rcio
       exit obsloop
    endif

   ! check the lat/lon values to see if they are ok
   !if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle obsloop
   !if ( lon <   0.0_r8 .or. lon >  360.0_r8 ) cycle obsloop


   ! if lon comes in between -180 and 180, use these lines instead:
   if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
   if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360

   ! put date into a dart time format
   time_obs = set_date(year, month, day, hour, minute, second)

   if (debug) call print_date(time_obs, 'next obs time is')

   !! if time is given in seconds since 1/1/1970, here's how to add it.
   !time_obs = comp_day0 + time_obs

   ! extract time of observation into gregorian day, sec.
   call get_time(time_obs, osec, oday)

!   if (retlev2 < 10000.0_r8) then
   print *, oday, osec
   qc_count = qc_count + 1
   ! this example assumes there is an obs type, where otype=1 is
   ! a temperature measured in height, and if otype=2, there's a wind
   ! speed and direction and height is pressure.  any kind of observation
   ! can use any of the vertical types; this is just an example.

   ! fixed to otype=1 for IASI O3
      ! no height since it is a column integrated quantity

      ! make an obs derived type, and then add it to the sequence
      call create_3d_obs(lat, lon, retlev2, VERTISHEIGHT, iasi_col, &
                         IASI_O3_RETRIEVAL, iasi_err, oday, osec, qc, obs, &
                         akcol, apcol_val, altretlev2, aircol_val, qc_count, nlev_use, apcol)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      if (debug) print *, 'added iasi obs to output seq'
!   endif
end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! End of main program

call finalize_utilities()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.
!
!       NOTE: assumes the code is using the threed_sphere locations module,
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!    metadata (AFAJ)
!    akcol - averaging kernel
!    apcol_val - a priori sub-column (I-A)xa
!    nlevels - number of levels
!    aircol_val - air sub column
!    qc_count - obs count (key)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, &
                         obs, akcol, apcol_val, altretlev, aircol_val, qc_count, nlev_use, apcol)
use        types_mod, only : r8
use      obs_def_mod, only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, set_obs_def_key
use obs_def_iasi_mod, only : set_obs_def_iasi_o3
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

integer,        intent(in)    :: okind, vkind, day, sec
real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
type(obs_type), intent(inout) :: obs
real(r8),       intent(in)    :: akcol(nlevels)
real(r8),       intent(in)    :: altretlev(nlevels)
real(r8),       intent(in)    :: apcol(nlevels)
real(r8),       intent(in)    :: aircol_val(nlevels)
real(r8),       intent(in)    :: apcol_val ! aircol_val
integer,        intent(in)    :: qc_count
integer,        intent(in)    :: nlev_use

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def_iasi_o3(qc_count, akcol(1:nlev_use), apcol_val, altretlev(1:nlev_use), &
                         aircol_val(1:nlev_use), nlev_use, apcol(1:nlev_use))
call set_obs_def_key(obs_def, qc_count)
!call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)
call set_obs_def(obs, obs_def)

end subroutine create_3d_obs


end program iasi_ascii_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

