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

use      utilities_mod, only : initialize_utilities, finalize_utilities, &
                               open_file, close_file, find_namelist_in_file, &
                               check_namelist_read, nmlfileunit, get_unit, &
                               do_nml_file, do_nml_term, get_next_filename, &
                               error_handler, E_ERR, E_MSG

use   time_manager_mod, only : time_type, set_calendar_type, set_date, get_date, &
                               operator(>=), increment_time, set_time, get_time, &
                               operator(-), operator(+), GREGORIAN, &
                               print_time, print_date

use       location_mod, only : VERTISUNDEF

use   obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                               static_init_obs_sequence, init_obs, write_obs_seq, & 
                               init_obs_sequence, get_num_obs, & 
                               set_copy_meta_data, set_qc_meta_data

use  obs_utilities_mod, only : create_3d_obs, add_obs_to_seq

use       obs_kind_mod, only : KIND_BRIGHTNESS_TEMPERATURE

use EASE_utilities_mod, only : get_grid_dims, ezlh_inverse, deconstruct_filename, &
                               read_ease_Tb

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
integer, parameter :: max_num_input_files = 500
integer            :: num_input_files = 0  ! actual number of files
integer            :: ifile, istatus
character(len=256), dimension(max_num_input_files) :: filename_seq_list

! information gleaned from filenaming convention
integer          :: iyear, idoy, channel
character(len=2) :: gridarea
character(len=1) :: passdir
character(len=1) :: polarization
logical          :: is_time_file

character(len=256) :: input_line
character(len=256) :: msgstring1,msgstring2,msgstring3

integer :: oday, osec, iocode, iunit, otype
integer :: year, month, day, hour, minute, second
integer :: num_copies, num_qc, max_obs
           
logical  :: file_exist, first_obs

! The EASE grid is a 2D array of 16 bit unsigned integers
integer, allocatable, dimension(:,:) :: Tb
integer, dimension(2) :: ndims
integer :: nrows, ncols, irow, icol

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


! TJH FIXME ... using 1 input file
! TJH FIXME ... using 1 input file
num_input_files = 1
! TJH FIXME ... using 1 input file
! TJH FIXME ... using 1 input file


! need some basic information from the first file
istatus = deconstruct_filename( filename_seq_list(1), &
             gridarea, iyear, idoy, passdir, channel, polarization, is_time_file)
if (istatus /= 0) then
   write(msgstring2,*) 'filename nonconforming'
   call error_handler(E_ERR, 'main', trim(filename_seq_list(ifile)), &
                 source, revision, revdate, text2=msgstring2)
endif

ndims = get_grid_dims(gridarea)
nrows = ndims(1)
ncols = ndims(2)

allocate( Tb(nrows,ncols) )

! each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
! There is a lot of observations per day. Should only do 1 day at at time.
max_obs    = nrows*ncols*num_input_files ! overkill
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

iunit = get_unit()  ! Get free unit for AMSR-E data file
FileLoop: do ifile = 1,num_input_files

   ! A little helpful logging
   write(msgstring1,*)'.. Converting file',ifile,' of ',num_input_files
   call error_handler(E_MSG, 'main', msgstring1, text2 = trim(filename_seq_list(ifile)))

   istatus = deconstruct_filename( filename_seq_list(ifile), &
             gridarea, iyear, idoy, passdir, channel, polarization, is_time_file)
   if (istatus /= 0) then
      call error_handler(E_ERR, 'main', 'filename nonconforming', &
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
   iocode = read_ease_Tb(filename_seq_list(ifile), iunit, Tb)

   ! FIXME ... check row/col vs. col/row
   ROWLOOP: do irow=1,nrows
   COLLOOP: do icol=1,ncols

      if ( Tb(irow,icol) == 0 ) cycle COLLOOP

      ! Convert icol,irow to lat/lon using EASE routine
      iocode = ezlh_inverse(gridarea, real(icol), real(irow), rlat, rlon)
      if (iocode /= 0) cycle COLLOOP

      if ( verbose) then
         write(*,*)'icol,irow,lat,lon',icol,irow,rlat,rlon
      endif

      ! check the lat/lon values to see if they are ok
      if ( rlat >  90.0 .or. rlat <  -90.0 ) cycle COLLOOP
      if ( rlon <   0.0 .or. rlon >  360.0 ) cycle COLLOOP
   
      ! if lon comes in between -180 and 180, use these lines instead:
      !if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle COLLOOP
      !if ( lon < 0.0_r8 )  lon = lon + 360.0_r8 ! changes into 0-360
   
      ! make an obs derived type, and then add it to the sequence
      ! TJH FIXME ... polarization, frequency ... metadata
      temp = real(Tb(irow,icol),r8) / 10.0_r8
      call create_3d_obs(real(rlat,r8), real(rlon,r8), 0.0_r8, VERTISUNDEF, temp, &
                            KIND_BRIGHTNESS_TEMPERATURE, terr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   
   enddo COLLOOP
   enddo ROWLOOP

enddo FileLoop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

deallocate(Tb)

! end of main program
call finalize_utilities()


contains


function Check_Input_Files(input_list, output_list)
! Read a list of files to process 
character(len=*),               intent(in)  :: input_list
character(len=*), dimension(:), intent(out) :: output_list
integer                                     :: Check_Input_Files

integer, parameter :: MAXLINES = 1000
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
      output_list(Check_Input_Files) = adjustl(input_line)
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
