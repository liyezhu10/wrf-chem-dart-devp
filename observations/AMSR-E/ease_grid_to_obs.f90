! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ease_grid_to_obs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   program to convert a series of flat binary files 
!
!   created 7 Nov 2013   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use          types_mod, only : r4, r8, PI, DEG2RAD

use          netcdf

use      utilities_mod, only : initialize_utilities, finalize_utilities, &
                               open_file, close_file, find_namelist_in_file, &
                               check_namelist_read, nmlfileunit, get_unit, &
                               do_nml_file, do_nml_term, get_next_filename, &
                               error_handler, E_ERR, E_MSG, file_exist, nc_check

use   time_manager_mod, only : time_type, set_calendar_type, set_date, get_date, &
                               operator(>=), increment_time, set_time, get_time, &
                               operator(-), operator(+), GREGORIAN, &
                               print_time, print_date

use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                               static_init_obs_sequence, init_obs, write_obs_seq, & 
                               init_obs_sequence, get_num_obs, &
                               set_copy_meta_data, set_qc_meta_data

use            location_mod, only : VERTISSURFACE, set_location
use       obs_utilities_mod, only : add_obs_to_seq, create_3d_obs
use            obs_kind_mod, only : AMSRE_BRIGHTNESS_T
use obs_def_brightnessT_mod, only : set_amsre_metadata
use      EASE_utilities_mod, only : get_grid_dims, ezlh_inverse, read_ease_Tb, &
                                    deconstruct_filename

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!----------------------------------------------------------------
! Namelist input with default values

character(len=256) :: input_file_list = 'file_list.txt'
character(len=256) :: obs_out_file = 'obs_seq.out'
logical            :: verbose = .false.

namelist /ease_grid_to_obs_nml/ &
         input_file_list, obs_out_file, verbose

!----------------------------------------------------------------

! max_num_input_files : max number of input files to be processed
integer, parameter :: max_num_input_files = 500  ! Long; change from 500
integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile, istatus
character(len=256), dimension(max_num_input_files) :: filename_seq_list

! information gleaned from filenaming convention
integer          :: iyear, idoy
character(len=2) :: gridarea
character(len=1) :: passdir
character(len=1) :: polarization
logical          :: is_time_file
real(r8)         :: footprint
real             :: frequency   ! real type to match EASE type

character(len=256) :: input_line
character(len=256) :: msgstring1,msgstring2,msgstring3
! character(len=256) :: outfilelog='/scratch/02714/zhaol/data/AMSR-E/&
!                 AMSR-E_interpolated_090_125/dart_obs_seq/latlonlog.txt' !==============Long

integer :: oday, osec, iocode, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
integer :: landcode = 0  ! FIXME ... totally bogus for now
           
logical  :: first_obs

! The interpolated EASE grid is a double array of data
real(r8), allocatable, dimension(:) :: Tb
real(r8), allocatable, dimension(:) :: RI6v
real(r8), allocatable, dimension(:) :: LON
real(r8), allocatable, dimension(:) :: LAT
integer :: key, icount, counts, ncid, VarID

real(r8) :: temp, terr, qc
real     :: rlat, rlon

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: cal_day0, time_obs, prev_time

!----------------------------------------------------------------
! start of executable code

call initialize_utilities('ease_grid_to_obs')

! time setup
call set_calendar_type(GREGORIAN)

! Read the namelist entry
call find_namelist_in_file("input.nml", "ease_grid_to_obs_nml", iunit)
read(iunit, nml = ease_grid_to_obs_nml, iostat = iocode)
call check_namelist_read(iunit, iocode, "ease_grid_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=ease_grid_to_obs_nml)
if (do_nml_term()) write(     *     , nml=ease_grid_to_obs_nml)

num_input_files = Check_Input_Files(input_file_list, filename_seq_list) 
write(*,*)' There are ',num_input_files,' input files.'

! need some basic information from the first file
istatus = deconstruct_filename( filename_seq_list(1), &
             gridarea, iyear, idoy, passdir, frequency, polarization, is_time_file)
if (istatus /= 0) then
   write(msgstring2,*) 'filename nonconforming'
   call error_handler(E_ERR, 'main', trim(filename_seq_list(1)), &
                 source, revision, revdate, text2=msgstring2)
endif

counts = 17199     !====192x288====================FIXME Long

allocate( Tb(counts) ) 
allocate( RI6v(counts) )  
allocate( LON(counts))
allocate( LAT(counts))

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
! There is a lot of observations per day. Should only do 1 day at at time.
max_obs    = counts*num_input_files ! overkill
num_copies = 1
num_qc     = 1

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call   set_qc_meta_data(obs_seq, 1,     'Data QC')

! --------------------------------------------------------------------------------
! Loop over all the input data files.

! open(77,file=outfilelog, status='unknown')      !=====================Long

iunit = get_unit()  ! Get free unit for AMSR-E data file
FileLoop: do ifile = 1,num_input_files

   ! A little helpful logging
   write(msgstring1,*)'.. Converting file',ifile,' of ',num_input_files
   call error_handler(E_MSG, 'main', msgstring1, text2 = trim(filename_seq_list(ifile)))

   istatus = deconstruct_filename( filename_seq_list(ifile), &
             gridarea, iyear, idoy, passdir, frequency, polarization, is_time_file)
   if (istatus /= 0) then
      call error_handler(E_ERR, 'main', 'filename nonconforming', &
             source, revision, revdate, text2= trim(filename_seq_list(ifile)))
   endif

   if (gridarea(2:2) == 'L') then
!      footprint = 25.0_r8  ! EASE grid is 25km resolution
      footprint = 100.0_r8  ! interpolated AMSR-E data is about 50km resolution
   else
      call error_handler(E_ERR, 'main', 'unknown footprint', &
             source, revision, revdate, text2= trim(filename_seq_list(ifile)))
   endif

   ! FIXME - crude time information at this point
   ! put date into a dart time format
   ! calculate first day of the year, then add in day-of-year
   time_obs = set_date(iyear,1,1,0,0,0)
   call get_time(time_obs, osec, oday)
   time_obs = set_time(osec, oday+idoy-1)
   call get_time(time_obs, osec, oday)

   if (verbose) then
      write(*,*)trim(filename_seq_list(ifile))
      call print_date(time_obs, str='file date is')
      call print_time(time_obs, str='file time is')
   endif

   ! Ignoring time files for now
   if ( is_time_file ) cycle FileLoop

   ! Fills up matrix of brightness temperatures

   call nc_check(nf90_open(trim(filename_seq_list(ifile)), nf90_nowrite, ncid), &
       'ease_grid_to_dart','open '//trim(filename_seq_list(ifile)))

   call nc_check(nf90_inq_varid(ncid, 'lat', VarID), &
            'ease_grid_to_dart', 'inq_varid lat')
   call nc_check(nf90_get_var(ncid, VarID, LAT), &
            'ease_grid_to_dart', 'get_var LAT')

   call nc_check(nf90_inq_varid(ncid, 'lon', VarID), &
            'ease_grid_to_dart', 'inq_varid lon')
   call nc_check(nf90_get_var(ncid, VarID, LON), &
            'ease_grid_to_dart', 'get_var LON')

   call nc_check(nf90_inq_varid(ncid, 'tb', VarID), &
            'ease_grid_to_dart', 'inq_varid tb')
   call nc_check(nf90_get_var(ncid, VarID, Tb), &
            'ease_grid_to_dart', 'get_var Tb')

   call nc_check(nf90_inq_varid(ncid, 'ri', VarID), &
            'ease_grid_to_dart', 'inq_varid ri')
   call nc_check(nf90_get_var(ncid, VarID, RI6v), &
            'ease_grid_to_dart', 'get_var RI6v')


   call nc_check(nf90_close(ncid), 'ease_grid_to_dart','close '//trim(filename_seq_list(ifile)))

   COUNTLOOP: do icount=1,counts

      if ( Tb(icount) ==  0.0_r8 ) cycle COUNTLOOP

      rlat = LAT(icount)
      rlon = LON(icount)

      ! ensure the lat/lon values are in range
      if ( rlat >  90.0_r8 .or. rlat <  -90.0_r8 ) cycle COUNTLOOP
      if ( rlon > 360.0_r8 .or. rlon <    0.0_r8 ) cycle COUNTLOOP

!==================================================================================Long

      if ( rlat < -60.0_r8 .or. rlat >   66.5_r8 ) then 
         if (frequency==6.9 .or. frequency==10.7) cycle COUNTLOOP 
      endif                                                              !=========Long
      
!==========for globe
      if ((rlon > 235.0_r8 .and. rlon < 300.0_r8 .and. rlat < 49.0_r8 .and. rlat > 25.0_r8) .or. &
!          (rlon > 35.0_r8  .and. rlon < 75.0_r8  .and. rlat < 35.0_r8 .and. rlat > 15.0_r8)) then
          (rlon > 35.0_r8  .and. rlon < 90.0_r8  .and. rlat < 35.0_r8 .and. rlat > 15.0_r8) .or. &
          RI6v(icount) >= 1.0_r8 ) then
         if (frequency==6.9) cycle COUNTLOOP
      else 
         if (frequency==10.7) cycle COUNTLOOP 
      endif

      if ( rlat < 25.0_r8 ) then 
         if (frequency==18.7 .or. frequency==23.8) cycle COUNTLOOP 
      endif                                                              !=========Long
      

!==================================================================================Long
  
      ! make an obs derived type, and then add it to the sequence
      temp = Tb(icount)

!     if ( verbose) then
!        write(77,'(A20,i6,f10.4,f10.4,f8.1)')'icount,lat,lon,Tb',icount,rlat,rlon,temp
!     endif
 
      terr = 2.0_r8   ! temps are [200,300] so this is about one percent
!      terr = 5.0_r8  ! set terr to a larger value
!      terr =Tbstd(icount) ! use std of Tbs within this CLM grids when upscaling
      
      call set_amsre_metadata(key, real(frequency,r8), footprint, polarization, landcode)

      call create_3d_obs(real(rlat,r8), real(rlon,r8), 0.0_r8, VERTISSURFACE, temp, &
                            AMSRE_BRIGHTNESS_T, terr, oday, osec, qc, obs, key)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
    

   enddo COUNTLOOP

enddo FileLoop

! close(77)

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

deallocate(Tb)
deallocate(RI6v)
deallocate(LON)
deallocate(LAT)

! end of main program
call finalize_utilities()


contains


function Check_Input_Files(input_list, output_list)
! Read a list of files to process 
character(len=*),               intent(in)  :: input_list
character(len=*), dimension(:), intent(out) :: output_list
integer                                     :: Check_Input_Files

character(len=256) :: ladjusted
integer, parameter :: MAXLINES = 500             ! Long; change from 1000
integer :: iline

Check_Input_files = -1

iunit = open_file(trim(input_list), 'formatted', 'read')

Check_Input_Files = 0
FileNameLoop: do iline = 1,MAXLINES ! a lot of lines 

   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=iocode) input_line
   if (iocode > 0) then 
      write(msgstring1,*) 'While reading ', trim(input_list)
      write(msgstring2,*) 'got read code (iostat) = ', iocode,' around line ',iline
      call error_handler(E_ERR, 'Check_Input_Files', msgstring1, &
                    source, revision, revdate, text2=msgstring2)
   elseif (iocode < 0) then 
      ! Normal end of file
      exit FileNameLoop
   else
      Check_Input_Files = Check_Input_Files + 1

      ladjusted = adjustl(input_line)
      if ( file_exist(trim(ladjusted)) ) then
         output_list(Check_Input_Files) = trim(ladjusted)
      else
         write(msgstring1,*)'file does not exist.'
         call error_handler(E_ERR,'Check_Input_Files',&
          msgstring1,source,revision,revdate,text2=trim(ladjusted))
      endif
   endif

enddo FileNameLoop

if (Check_Input_Files >= MAXLINES-1 ) then
   write(msgstring1,*)'Too many files to process. Increase MAXLINES and try again.'
   call error_handler(E_ERR,'Check_Input_Files',msgstring1,source,revision,revdate)
endif

end function Check_Input_Files


end program ease_grid_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
