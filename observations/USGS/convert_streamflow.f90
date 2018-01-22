! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program convert_streamflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_streamflow - reads streamflow data and writes a DART 
!                      obs_seq file using the DART library routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8

use      location_mod, only : VERTISHEIGHT

use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities, &
                              nmlfileunit, do_nml_file, do_nml_term, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG

use  time_manager_mod, only : time_type, set_calendar_type, set_date, operator(>=), &
                              increment_time, get_time, operator(-), GREGORIAN

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : STREAM_FLOW

use obs_utilities_mod, only : getvar_real, getvar_int, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, set_missing_name

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = 'convert_streamflow:'

character(len=512) :: string1, string2 ! strings for messages

! input file has data quality control fields, whether to use or ignore them.

integer,  parameter :: NUM_COPIES      = 1      ! number of copies in sequence
integer,  parameter :: NUM_QC          = 1      ! number of QC entries
real(r8), parameter :: MIN_OBS_ERR_STD = 0.5_r8 ! m^3/sec

integer :: nobs, n, i, nlinks, iunit, io, indx
integer :: ncid_d, ncid_l, varid
logical :: file_exist, first_obs

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: prev_obs
type(time_type)         :: time_obs, prev_time

! These should/must match what is in the netCDF files ... should error check
! or at least provide a decent error message.

integer, parameter :: IDLength = 15
integer, parameter :: stationIdStrLen = 15
integer, parameter :: timeStrLen = 19

character(len=timeStrLen),      allocatable ::    time_string(:)
character(len=stationIdStrLen), allocatable :: station_string(:)
character(len=IDLength),        allocatable ::    gage_string(:)
integer,                        allocatable :: discharge_quality(:)
character(len=stationIdStrLen), allocatable :: gages_list(:)

! variables needed for create_3d_obs 

real(r8), allocatable :: lat(:)
real(r8), allocatable :: lon(:)
real(r8), allocatable :: altitude(:)
real(r8), allocatable :: discharge(:)
real(r8)              :: oerr, qc
integer               :: oday, osec
type(obs_type)        :: obs

! namelist variables
character(len=256) :: input_files     = 'inputs.txt' ! This exists in the namelist ...
! which is the address to a file having all the files we want to process 
character(len=256) :: input_file      = 'input.nc'  ! This would be a line of ...
! input_files which is one time slice file 
character(len=256) :: location_file   = 'location.nc'
character(len=256) :: output_file     = 'obs_seq.out'
character(len=256) :: gages_list_file = 'gage_lists.txt'
real(r8)           :: obs_fraction_for_error = 0.01
integer            :: verbose = 0

! local variabls
integer            :: fileEnd = 0
integer            :: nlines  = 0

namelist / convert_streamflow_nml / input_files, output_file, location_file, &
                                    obs_fraction_for_error, verbose, &
                                    gages_list_file

! Get going

call initialize_utilities(routine)

! Read the DART namelist
call find_namelist_in_file('input.nml', 'convert_streamflow_nml', iunit)
read(iunit, nml = convert_streamflow_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_streamflow_nml')

! print the content of the namelist for clarification 
print*, "input_files : ", input_files
print*, "output_file : ", output_file
print*, "location_file : ", location_file
print*, "obs_fraction_for_error : ", obs_fraction_for_error
print*, "verbose : ", verbose
print*, "gages_list_file : ", gages_list_file

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_streamflow_nml)
if (do_nml_term()) write(     *     , nml=convert_streamflow_nml)

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

!*****************************************************************************************
! READ the file with the location information 
!*****************************************************************************************
call nc_check( nf90_open(location_file, nf90_nowrite, ncid_l), &
             routine, 'opening file '//trim(location_file) )
call getdimlen(ncid_l, "linkDim", nlinks)

write(string1,*)'number of links is ',nlinks
call error_handler(E_MSG, routine, string1)

allocate( lat(nlinks))
allocate( lon(nlinks))
allocate(altitude(nlinks))
allocate(gage_string(nlinks))

! read in the data arrays
call getvar_real(     ncid_l, "lat",  lat       )  ! latitudes
call getvar_real(     ncid_l, "lon",  lon       )  ! longitudes
call getvar_real(     ncid_l, "alt",  altitude  )  ! elevation
call get_gage_strings(ncid_l, "gages"           ) ! character string ID

! convert all lons to [0,360], replace any bad lat lims
! check the lat/lon values to see if they are ok
! change lon from -180 to 180 into 0-360

where (lat >  90.0_r8) lat =  90.0_r8
where (lat < -90.0_r8) lat = -90.0_r8
where (lon <   0.0_r8) lon = lon + 360.0_r8

call nc_check(nf90_close(ncid_l), routine, 'closing file '//trim(location_file))

call static_init_obs_sequence()
call init_obs(obs,      num_copies=NUM_COPIES, num_qc=NUM_QC)
call init_obs(prev_obs, num_copies=NUM_COPIES, num_qc=NUM_QC)

!*****************************************************************************************
!    Read the list of the gages to proprocess
!*****************************************************************************************

!>@todo FIXME ... remove hardcoded units, use function to count lines, etc.

! find out the number of lines 
if (gages_list_file /= "") then 
        nlines = 0
        OPEN(222,FILE=gages_list_file)
        do
         READ(222, *, IOSTAT = fileEnd) 
         if (fileEnd == 0) then
         nlines = nlines + 1
         else if (fileEnd < 0 ) then
                exit
         end if
        enddo
        CLOSE(222)

        print*, "number of gages to read from the time slices :", nlines

        ! allocate ...
        allocate(gages_list(nlines))
        OPEN(222,FILE=gages_list_file)
        do i = 1, nlines
                READ(222, '(A15)' , IOSTAT = fileEnd), gages_list(i)
        enddo
        print*, "List of the gages to process", gages_list
endif

!*****************************************************************************************
!   Loop through the time slices and adding the data to obs_seq.out 
!   if the gage is in the gage_list
!*****************************************************************************************

OPEN(111,FILE=input_files)
do
 READ(111, '(A256)', IOSTAT = fileEnd) input_file
 print*, input_file
 if (fileEnd > 0 ) then 
     Print*, "There was a probelm reading : ", input_files
     stop
 else if (fileEnd == 0 ) then 
        
        first_obs = .true.
        call nc_check( nf90_open(input_file, nf90_nowrite, ncid_d), &
                       routine, 'opening file '//trim(input_file) )

        call getdimlen(ncid_d, "stationIdInd", nobs)

        write(string1,*)'number of obs in the time slice is ',nobs
        call error_handler(E_MSG, routine, string1)

        allocate(discharge(nobs))
        allocate(discharge_quality(nobs))
        allocate(time_string(nobs))
        allocate(station_string(nobs))

        ! read in the data arrays
        call getvar_real(ncid_d, "discharge", discharge )  ! streamflow
        call getvar_int( ncid_d, "discharge_quality", discharge_quality)

        call get_station_strings(ncid_d, "stationId" ) ! Name of the gage
        call get_time_strings(   ncid_d, "time"      ) ! observation time

        call nc_check(nf90_close(ncid_d), routine, 'closing file '//trim(input_file) )

        ! Done with reading the data ... 

        !  either read existing obs_seq or create a new one

        inquire(file=output_file, exist=file_exist)
        if ( file_exist ) then

          ! existing file found, append to it
          call read_obs_seq(output_file, 0, 0, nobs, obs_seq)

        else

          ! create a new one ... 
          call init_obs_sequence(obs_seq, NUM_COPIES, NUM_QC, nobs)  !!!!Arezoo: I want to reduce the number of observation .... to what is in the gage list 
        
          do i=1,NUM_COPIES ! kinda silly ... only 1 type of observation
             call set_copy_meta_data(obs_seq, i, 'observation')
          enddo
          do i=1,NUM_QC ! kinda silly ... only 1 type of qc
             call set_qc_meta_data(obs_seq, i, 'Data QC')
          enddo

        endif
        ! Set the DART data quality control.  Be consistent with NCEP codes;
        ! 0 is 'must use', 1 is good, no reason not to use it.

        qc = 1.0_r8   ! modify based on discharge_quality ... perhaps

        obsloop1: do n = 1, nobs

           if (.not.allocated(gages_list) .or. ANY(gages_list == adjustl(station_string(n)))) then 
               if ( discharge(n) < 0.0_r8 ) cycle

               ! oerr is the observation error standard deviation in this application.
               ! The observation error variance encoded in the observation file
               ! will be oerr*oerr
               oerr = max(discharge(n)*obs_fraction_for_error, MIN_OBS_ERR_STD)

               call convert_time_string(time_string(n),oday,osec)

               ! find the lat/lon for the station,link that match this obs.
               indx = find_matching_gage_index(n)

               call create_3d_obs(lat(indx), lon(indx), altitude(indx), VERTISHEIGHT, discharge(n), &
                                      STREAM_FLOW, oerr, oday, osec, qc, obs)
               call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
           endif
        end do obsloop1
    

        ! if we added any obs to the sequence, write it now.
        if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, output_file)

        deallocate(discharge, discharge_quality)
        deallocate(time_string, station_string)

 else if (fileEnd < 0 ) then
        print*, "Reached End of the file"
        exit
 end if
        
end do
CLOSE(111)
stop

! end of main program

deallocate(lat,lon, altitude, discharge, discharge_quality)
deallocate(gage_string, time_string, station_string)
if (allocated(gages_list)) deallocate(gages_list)
call finalize_utilities()

contains

!-----------------------------------------------------------------------
!> Match the stationId string to the gage string to determine which gage
!> has the location matching the stationId.

function find_matching_gage_index(station_id) result(id)

integer, intent(in) :: station_id
integer :: id

integer :: i

i = 0

LINKS : do i = 1,nlinks
   if (station_string(station_id) == gage_string(i)) then
      if (verbose > 0) then
         write(string1,*)'identified station "',station_string(station_id),'"'
         call error_handler(E_MSG, 'find_matching_gage_index:', string1)
      endif
      id = i
      exit LINKS
   endif
enddo LINKS

if (i == 0) then
   write(string1,*)'Unable to match station id for obs #',station_id,' "',station_string(station_id),'"'
   call error_handler(E_ERR, 'find_matching_gage_index', string1, source, revision, revdate)
endif

end function find_matching_gage_index


!-----------------------------------------------------------------------
!> read the character matrix from the netCDF file and parse into
!> useable strings


subroutine convert_time_string(string,days,seconds)

character(len=*), intent(in) :: string
integer, intent(out) :: days, seconds

integer :: year, month, day, hour, minute, second
type(time_type) :: darttime

read(string,'(i4,5(1x,i2.2))') year, month, day, hour, minute, second

if (verbose > 1) then
   write(string1,*) ' ..  read ',trim(string)
   write(string2,*)' interpreted as ', year, month, day, hour, minute, second
   call error_handler(E_MSG,'convert_time_string:',string1,text2=string2)
endif

darttime = set_date(year, month, day, hour, minute, second)
call get_time(darttime, seconds, days)

end subroutine convert_time_string


!-----------------------------------------------------------------------
!> read the character matrix of UTC times from the data netCDF file and
!> parse into an array of strings.
!> dimensions:
!>         stationIdInd = UNLIMITED ; // (2 currently)
!>         timeStrLen = 19 ;
!> variables:
!>         char time(stationIdInd, timeStrLen) ;
!>              time:units = "UTC" ;
!>              time:long_name = "YYYY-MM-DD_HH:mm:ss UTC" ;

subroutine  get_time_strings(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname

integer :: io, dim1

call getdimlen(ncid, "timeStrLen", dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_time_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,time_string)
call nc_check(io, 'get_time_strings', 'get_var "'//varname//'"')

if (verbose > 2) then
   do i = 1,nobs
      write(string1,*) 'time ',i,' is "'//time_string(i)//'"'
      call error_handler(E_MSG, 'get_time_strings:', string1)
   enddo
endif

end subroutine  get_time_strings


!-----------------------------------------------------------------------
!> read the character matrix of USGS station identifiers from the data
!> netCDF file and parse into an array of strings.
!> dimensions:
!>         stationIdInd = UNLIMITED ; // (2 currently)
!>         timeStrLen = 19 ;
!> variables:
!>         char stationId(stationIdInd, stationIdStrLen) ;
!>              stationId:long_name = "USGS station identifer of length 15" ;


subroutine  get_station_strings(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname

integer :: io, dim1

call getdimlen(ncid, "stationIdStrLen", dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_station_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,station_string)
call nc_check(io, 'get_station_strings', 'get_var "'//varname//'"')

if (verbose > 2) then
   do i = 1,nobs
      write(string1,*) 'station ',i,' is "'//station_string(i)//'"'
      call error_handler(E_MSG, 'get_station_strings:', string1)
   enddo
endif

end subroutine  get_station_strings


!-----------------------------------------------------------------------
!> read the character matrix of NHD Gage Event IDs  from the metadata
!> netCDF file and parse into an array of strings.
!> dimensions:
!>        linkDim = 157 ;
!>        IDLength = 15 ;
!> variables:
!>        char gages(linkDim, IDLength) ;
!>             gages:long_name = "NHD Gage Event ID from SOURCE_FEA field in Gages feature class" ;
!>             gages:coordinates = "lat lon" ;


subroutine  get_gage_strings(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname

integer :: io, dim1

call getdimlen(ncid, "linkDim", dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_gage_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,gage_string)
call nc_check(io, 'get_gage_strings', 'get_var "'//varname//'"')

if (verbose > 2) then
   do i = 1,nlinks
      write(string1,*) 'gage_string ',i,' is "'//gage_string(i)//'"'
      call error_handler(E_MSG, 'get_gage_strings:', string1)
   enddo
endif

end subroutine  get_gage_strings


end program convert_streamflow

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
