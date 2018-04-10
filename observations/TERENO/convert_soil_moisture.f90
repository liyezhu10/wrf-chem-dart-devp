! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_soil_moisture

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_soil_moisture - program that reads a TERENO netCDF file and 
!     writes a DART obs_seq file. The TERENO soil moisture observations
!     are taken by two instruments about 8cm apart to increase 
!     representativeness - we are simply going to average the two.
!
! Apr 2018, Tim Hoar, NCAR/DAReS
!
! SoilWaterContent_0.2mSensor2_QC:comment = "Processing status (qcdim[0]) and quality flag (qcdim[1]) 
! is mapped to variable processing_status and quality_flag" ;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8, digits12

use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              open_file, close_file, find_namelist_in_file, &
                              check_namelist_read, nmlfileunit, logfileunit, &
                              do_nml_file, do_nml_term, file_exist, &
                              find_textfile_dims, get_next_filename, &
                              error_handler, E_ERR, E_MSG

use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              get_time, set_time, operator(+), GREGORIAN, &
                              print_time, print_date

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : SOIL_MOISTURE

use obs_utilities_mod, only : getvar_real, getvar_real_2d, getvar_int_3d, &
                              getdimlen, create_3d_obs, add_obs_to_seq

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = "convert_soil_moisture"

character(len=512) :: string1, string2, string3

!-----------------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file    = 'obs_seq.out'
logical            :: verbose         = .false.

namelist /convert_soil_moisture_nml/ &
         input_file_list, obs_out_file, verbose

!-----------------------------------------------------------------------
! MAX_NUM_INPUT_FILES : max number of input files to be processed

integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile
character(len=256), dimension(MAX_NUM_INPUT_FILES) :: filename_seq_list
character(len=256) :: filename

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: io, iunit, ncid, i, oday, osec
logical  :: first_obs
integer  :: idepth, istation, nstations, itime, ntimes, nqcdims
integer  :: counts, max_obs
real(r8) :: obs_val, err_std, qc

real(r8), allocatable :: latitude(:), longitude(:), altitude(:) 
real(r8), allocatable ::   sensor1(:,:,:),   sensor2(:,:,:)
integer,  allocatable :: qcsensor1(:,:,:,:), qcsensor2(:,:,:,:)
type(time_type), allocatable :: dart_time(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: prev_time

integer, parameter :: NDEPTHS = 3
real(r8) :: depths(NDEPTHS) = (/ 0.05, 0.2, 0.5 /) ! meters
character(len=*), parameter :: varname1(NDEPTHS) = &
                (/'SoilWaterContent_0.05mSensor1', &
                  'SoilWaterContent_0.2mSensor1 ', &
                  'SoilWaterContent_0.5mSensor1 '/)
character(len=*), parameter :: varname2(NDEPTHS) = &
                (/'SoilWaterContent_0.05mSensor2', &
                  'SoilWaterContent_0.2mSensor2 ', &
                  'SoilWaterContent_0.5mSensor2 '/)

character(len=*), parameter :: qcvarname1(NDEPTHS) = &
                (/'SoilWaterContent_0.05mSensor1_QC', &
                  'SoilWaterContent_0.2mSensor1_QC ', &
                  'SoilWaterContent_0.5mSensor1_QC '/)
character(len=*), parameter :: qcvarname2(NDEPTHS) = &
                (/'SoilWaterContent_0.05mSensor2_QC', &
                  'SoilWaterContent_0.2mSensor2_QC ', &
                  'SoilWaterContent_0.5mSensor2_QC '/)

real(r8) ::   sensor1_miss(NDEPTHS),   sensor2_miss(NDEPTHS)
integer  :: qcsensor1_miss(NDEPTHS), qcsensor2_miss(NDEPTHS)

! These values came from a dump of 'processing_status,quality_flag'
!>@todo determine these from the string instead of the presumed index
integer, parameter :: processing_status = 4 ! processing_status == 'approved upload'
integer, parameter :: quality_flag = 2      ! quality_flag      == 'ok: ok'

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities(routine)

call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "convert_soil_moisture_nml", iunit)
read(iunit, nml=convert_soil_moisture_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_soil_moisture_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_soil_moisture_nml)
if (do_nml_term()) write(     *     , nml=convert_soil_moisture_nml)

num_input_files = Check_Input_Files(input_file_list, filename_seq_list)

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.

counts = 2 * 3 * 2890 * 300  ! FIXME ... nsensors * ndepths * nT * nStations

max_obs = counts*num_input_files ! overkill

! either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

if ( file_exist(obs_out_file) ) then ! existing file found, append to it

  call read_obs_seq(obs_out_file, 0, 0, max_obs, obs_seq)

else ! create a new one

  call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'observation')
  enddo
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  enddo

endif

!-----------------------------------------------------------------------
! Loop over all the input data files.

FileLoop: do ifile = 1,num_input_files

   filename = filename_seq_list(ifile)

   io = nf90_open(filename, nf90_nowrite, ncid)
   call nc_check(io, routine, 'opening "'//trim(filename)//'"')

   call getdimlen(ncid, "station", nstations)
   call getdimlen(ncid, "time", ntimes)
   call getdimlen(ncid, "qcdim", nqcdims)

   allocate(longitude(nstations), latitude(nstations), altitude(nstations))
   allocate(dart_time(ntimes))

   call getvar_real(ncid, "lon",       longitude)
   call getvar_real(ncid, "lat",       latitude )
   call getvar_real(ncid, "altitude",  altitude )

   if (verbose) then
      !>@todo better verbose ... what file, for example
      write(string1,*)'..  range of longitude ',minval(longitude),maxval(longitude)
      write(string2,*)    'range of  latitude ',minval( latitude),maxval( latitude)
      write(string3,*)    'range of  altitude ',minval( altitude),maxval( altitude)
      call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)
   endif

   ! ensure the lat/longitude values are in range
   where( latitude <  -90.0_r8 ) latitude  = -90.0_r8
   where( latitude >   90.0_r8 ) latitude  =  90.0_r8
   where( longitude <   0.0_r8 ) longitude = longitude + 360.0_r8
   where( longitude > 360.0_r8 ) longitude = longitude - 360.0_r8

   call convert_to_dart_time(ncid, dart_time, filename)

   allocate(  sensor1(NDEPTHS,ntimes,nstations),   sensor2(NDEPTHS,ntimes,nstations))
   allocate(qcsensor1(NDEPTHS,nqcdims,ntimes,nstations), qcsensor2(NDEPTHS,nqcdims,ntimes,nstations))

   ! read in the observation data arrays and their qcs

   !>@todo what is the second qcdim?

   do idepth = 1,NDEPTHS
      call getvar_real_2d(ncid,   varname1(idepth), sensor1(  idepth,:,:), sensor1_miss(idepth))
      call getvar_real_2d(ncid,   varname2(idepth), sensor2(  idepth,:,:), sensor2_miss(idepth))
      call getvar_int_3d( ncid, qcvarname1(idepth), qcsensor1(idepth,:,:,:), qcsensor1_miss(idepth))
      call getvar_int_3d( ncid, qcvarname2(idepth), qcsensor2(idepth,:,:,:), qcsensor2_miss(idepth))
   enddo

   call nc_check(nf90_close(ncid),routine,'closing "'//trim(filename)//'"')

   ! Having a NaN as the _FillValue is just the dumbest thing ever.
   ! By definition they fail every test ... so you cannot test for equality
   where(.not. (sensor1 <= 0.0) .and. .not. (sensor1 > 0.0)) sensor1 = MISSING_R8
   where(.not. (sensor2 <= 0.0) .and. .not. (sensor2 > 0.0)) sensor2 = MISSING_R8

   if (verbose) then
      write(string1,*)'..  range of sensor1 ',minval(  sensor1),maxval(  sensor1), &
                                              minval(qcsensor1),maxval(qcsensor1)
      write(string2,*)    'range of sensor2 ',minval(  sensor2),maxval(  sensor2), &
                                              minval(qcsensor2),maxval(qcsensor2)
      call error_handler(E_MSG,routine,string1,text2=string2)
   endif

   ! I want to add the all the observations at a particular time FIRST
   ! before going on to the next time.

   TIMELOOP:    do itime   =1,ntimes
      call get_time(dart_time(itime), osec, oday)

      STATIONLOOP: do istation=1,nstations
      DEPTHLOOP:   do idepth  =1,NDEPTHS

         if ( sensor1(idepth,itime,istation) == MISSING_R8 ) cycle DEPTHLOOP
         if ( sensor2(idepth,itime,istation) == MISSING_R8 ) cycle DEPTHLOOP

         ! If both sensors have good quality, average them.
         ! If only one sensor has good quality, only use it.

         ! processing_status and quality_flag came from the netCDF metadata
         if ( qcsensor1(idepth,1,itime,istation) == processing_status .and. &
              qcsensor1(idepth,2,itime,istation) == quality_flag .and. &
              qcsensor2(idepth,1,itime,istation) == processing_status .and. &
              qcsensor2(idepth,2,itime,istation) == quality_flag ) then
            obs_val = 0.5_r8*sensor1(idepth,itime,istation) + &
                      0.5_r8*sensor2(idepth,itime,istation)
            qc      = 0
         elseif ( qcsensor1(idepth,1,itime,istation) == processing_status .and. &
                  qcsensor1(idepth,2,itime,istation) == quality_flag .and. &
            obs_val = sensor1(idepth,itime,istation)
            qc      = 0
         elseif ( qcsensor2(idepth,1,itime,istation) == processing_status .and. &
                  qcsensor2(idepth,2,itime,istation) == quality_flag ) then
            obs_val = sensor2(idepth,itime,istation)
            qc      = 0
         else
            cycle DEPTHLOOP
         endif
 
         write(*,*)obs_val, qcsensor1(idepth,:,itime,istation), qcsensor2(idepth,:,itime,istation)
 
         ! The observation error is either 20% or 0.02 
         err_std = max(obs_val/20.0_r8, 0.02_r8)
  
         !>@todo be more descriptive than SOIL_MOISTURE ... what kind of instrument 

         call create_3d_obs(latitude(istation), longitude(istation), depths(idepth), &
                 VERTISHEIGHT, obs_val, SOIL_MOISTURE, err_std, oday, osec, qc, obs)
   
         call add_obs_to_seq(obs_seq, obs, dart_time(itime), prev_obs, prev_time, first_obs)

      enddo DEPTHLOOP
      enddo STATIONLOOP
   enddo TIMELOOP

   deallocate(longitude, latitude, altitude, dart_time)
   deallocate(sensor1, sensor2, qcsensor1, qcsensor2)

enddo FileLoop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   write(string1,*) 'writing observations: obs_count = ', get_num_obs(obs_seq)
   call error_handler(E_MSG, routine, string1)
   call write_obs_seq(obs_seq, obs_out_file)
else
   call error_handler(E_MSG, routine, 'no observations to write out.')
endif

! end of main program
call finalize_utilities()

contains 

!-----------------------------------------------------------------------
!> Read a list of files to process.
!> Make sure each of the files exists.

function Check_Input_Files(input_list, output_list)
character(len=*), intent(in)  :: input_list      !> filename containing list
character(len=*), intent(out) :: output_list(:)
integer                       :: Check_Input_Files

character(len=256) :: filename
character(len=256) :: ladjusted
integer :: iunit, iline, nlines

character(len=*), parameter :: routine='Check_Input_Files'

Check_Input_files = -1

call find_textfile_dims(input_list, nlines)

iunit = open_file(trim(input_list), 'formatted', 'read')

if (nlines >= MAX_NUM_INPUT_FILES ) then
   write(string1,*)'Too many files to process. Increase MAX_NUM_INPUT_FILES, recompile, and try again.'
   write(string2,*)'MAX_NUM_INPUT_FILES currently set to ',MAX_NUM_INPUT_FILES
   write(string3,*)'There were ',nlines,' files specified in ',trim(input_list)
   call error_handler(E_ERR,routine,string1,source,revision,revdate, text2=string2, text3=string3)
endif

Check_Input_Files = 0
FileNameLoop: do iline = 1,nlines ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=io) filename
   if (io > 0) then
      write(string1,*) 'While reading ', trim(input_list)
      write(string2,*) 'got read code (iostat) = ', io,' around line ',iline
      call error_handler(E_ERR, routine, string1, &
                    source, revision, revdate, text2=string2)
   elseif (io < 0) then
      ! Normal end of file
      exit FileNameLoop
   else
      Check_Input_Files = Check_Input_Files + 1

      ladjusted = adjustl(filename)
      if ( file_exist(ladjusted) ) then
         output_list(Check_Input_Files) = ladjusted
      else
         write(string1,*)'following file does not exist:'
         call error_handler(E_ERR, routine, string1, &
             source, revision, revdate, text2='"'//trim(ladjusted)//'"')
      endif
   endif

enddo FileNameLoop

call close_file(iunit)

end function Check_Input_Files


!-----------------------------------------------------------------------
!> read the observation time array and convert to an array of DART times
!>  double time(time) ;
!>         time:long_name = "time" ;
!>         time:units = "milliseconds since 1970-01-01T00:00:00" ;

subroutine convert_to_dart_time(ncid, dart_times, filename) 

integer,                    intent(in)  :: ncid
type(time_type),            intent(out) :: dart_times(:)
character(len=*), optional, intent(in)  :: filename

character(len=*), parameter :: routine = 'convert_to_dart_time'

type(time_type) :: origin, offset
real(digits12)  :: timearray(ntimes)
real(digits12)  :: seconds(ntimes)
integer :: io, VarID
integer :: iyear,imonth,iday,ihour,imin,isec
character(len=NF90_MAX_NAME) :: attvalue
character(len=512) :: mystring

if (present(filename)) then
   mystring = '"'//trim(filename)//'"'
else
   mystring = ''
endif

io = nf90_inq_varid(ncid, 'time', VarID)
call nc_check(io, routine, 'inq_varid time '//mystring)

io = nf90_get_var(ncid, VarID, timearray)
call nc_check(io, routine, 'get_var   time '//mystring)

io = nf90_get_att(ncid, VarID, 'units', attvalue)
call nc_check(io, routine, 'time get_att units '//mystring)

! time:units = "milliseconds since 1970-01-01T00:00:00" ;
!               1234567890123456789

if (attvalue(1:18) /= 'milliseconds since') then
   write(string1,*)'expecting time units of [milliseconds since ... ]'
   write(string2,*)'read time units of ['//trim(attvalue)//']'
   call error_handler(E_ERR, routine, string1, &
          source, revision, revdate, text2=string2,text3=mystring)
endif

read(attvalue,'(19x,i4,5(1x,i2))',iostat=io)iyear,imonth,iday,ihour,imin,isec
if (io /= 0) then
   write(string1,*)'Unable to read time units ',trim(mystring)
   write(string2,*)'expected "milliseconds since YYYY-MM-DDTHH:MM:SS"'
   write(string3,*)'was      "'//trim(attvalue)//'"'
   call error_handler(E_ERR, routine, string1, &
          source, revision, revdate, text2=string2, text3=string3)
endif

origin = set_date(iyear, imonth, iday, ihour, imin, isec)

! convert each time to a DART time and compare to desired
! first, convert msec to seconds (without overflowing?)
! DART only uses integer seconds as lowest unit of time ...

seconds = timearray / 1000.0_digits12

TIMELOOP : do itime = 1,ntimes

   iday   = int(seconds(itime)/86400.0_digits12)
   isec   = int(seconds(itime) - real(iday,digits12)*86400.0_digits12)
   offset = set_time(isec, iday)

   dart_time(itime) = offset + origin

enddo TIMELOOP

if (verbose) then
   call print_time(dart_time(1),     str='first obs time is ',iunit=logfileunit)
   call print_time(dart_time(1),     str='first obs time is ')
   call print_time(dart_time(ntimes),str='last  obs time is ',iunit=logfileunit)
   call print_time(dart_time(ntimes),str='last  obs time is ')

   call print_date(dart_time(1),     str='first obs date is ',iunit=logfileunit)
   call print_date(dart_time(1),     str='first obs date is ')
   call print_date(dart_time(ntimes),str='last  obs date is ',iunit=logfileunit)
   call print_date(dart_time(ntimes),str='last  obs date is ')
endif

end subroutine convert_to_dart_time

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
