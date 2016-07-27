! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!----------------------------------------------------------------------
! This module provides support for observations from MOPITT.
! Each MOPITT installation has soil properties that have been measured
! and are needed by the MOPITT algorithm that converts moisture to
! counts of neutron intensity. These properties constitute additional
! metadata for each observation. The routines in this module read and
! write that metadata and provide the 'get_expected_neutron_intensity()'
! routine that is the 'observation operator'
!
!  OBS            1
!           -1           2          -1
! obdef
! loc3d
! 3.687678362484394        0.2270109638878783         90000.00000000000  -2
! kind
!          114
!          10
!   0.519130249717886
!    101188.200000000
!  -9.971600253259284E-002 -0.179213126683249      -0.186301389096441
!  -0.251874414846400      -0.347876389743044      -0.429032618622858
!  -0.471976941876570      -0.481665813523456      -0.305828023764864
!  -3.167302937284040E-002
!        11680
!  73129     148814
!
!----------------------------------------------------------------------

! BEGIN DART PREPROCESS KIND LIST
! MOPITT_CO_RETRIEVAL,           KIND_MOPITT_CO_RETRIEVAL
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_MOPITT_mod, only : read_MOPITT_metadata, &
!                                 write_MOPITT_metadata
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!      case(MOPITT_CO_RETRIEVAL)
!         continue
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!   case(MOPITT_CO_RETRIEVAL)
!      call read_MOPITT_metadata(obs_def%key, key, ifile, fform)
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!   case(MOPITT_CO_RETRIEVAL)
!      call write_MOPITT_metadata(obs_def%key, ifile, fform)
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!   case(MOPITT_CO_RETRIEVAL)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_MOPITT_mod

! <next few lines under version control, do not edit>
! $URL: $
! $Id: $
! $Revision: $
! $Date: $

use        types_mod, only : r8, PI, metadatalength, MISSING_R8, MISSING_I
use    utilities_mod, only : register_module, error_handler, E_ERR, E_WARN, E_MSG, &
                             ascii_file_format
use     obs_kind_mod, only : KIND_MOPITT_CO_RETRIEVAL

use typesizes
use netcdf

implicit none
private

public ::            set_MOPITT_metadata, &
                     get_MOPITT_metadata, &
                    read_MOPITT_metadata, &
                   write_MOPITT_metadata

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: $", &
   revision = "$Revision: $", &
   revdate  = "$Date: $"

character(len=256) :: string1, string2
logical, save      :: module_initialized = .false.

! Metadata for MOPITT observations.
! There are soil parameters for each site that must be added to each
! observation in the sequence. Also MOPITT parameters ...

integer, parameter :: MAXNLEVELS = 20

type site_metadata
   private
   integer  :: nlevels
   real(r8) :: mystery
   real(r8) :: psurf
   real(r8),dimension(MAXNLEVELS) :: profile
   integer  :: retrieval_number
end type site_metadata

type(site_metadata), allocatable, dimension(:) :: observation_metadata
type(site_metadata) :: missing_metadata

integer :: i
logical :: debug = .FALSE.
integer :: MAXMOPITTkey = 1000000  ! one year of hourly data - to start
integer ::    MOPITTkey = 0       ! useful length of metadata arrays

!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------


  subroutine initialize_module
!----------------------------------------------------------------------------
! subroutine initialize_module
!

call register_module(source, revision, revdate)

module_initialized = .true.

missing_metadata%nlevels          = MISSING_I
missing_metadata%mystery          = MISSING_R8
missing_metadata%psurf            = MISSING_R8
missing_metadata%profile          = MISSING_R8
missing_metadata%retrieval_number = MISSING_I

allocate(observation_metadata(MAXMOPITTkey))

observation_metadata(:) = missing_metadata

end subroutine initialize_module



 subroutine set_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )
!----------------------------------------------------------------------
!subroutine set_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )
!
! Fill the module storage metadata for a particular observation.

integer,  intent(out) :: key
integer,  intent(in)  :: nlevels
real(r8), intent(in)  :: mystery
real(r8), intent(in)  :: psurf
real(r8), intent(in)  :: profile(:)
integer,  intent(in)  :: retrieval_number

if ( .not. module_initialized ) call initialize_module

MOPITTkey = MOPITTkey + 1  ! increase module storage used counter

! Make sure the new key is within the length of the metadata arrays.
call grow_metadata(MOPITTkey,'set_MOPITT_metadata')

key = MOPITTkey ! now that we know its legal

! FIXME : make sure nlevels is less than length(profile)

observation_metadata(key)%nlevels            = nlevels
observation_metadata(key)%mystery            = mystery
observation_metadata(key)%psurf              = psurf
observation_metadata(key)%profile(1:nlevels) = profile(1:nlevels)
observation_metadata(key)%retrieval_number   = retrieval_number

end subroutine set_MOPITT_metadata



 subroutine get_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )
!----------------------------------------------------------------------
!subroutine get_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )
!
! Query the metadata in module storage for a particular observation.
! This can be useful for post-processing routines, etc.

integer,                intent(in)  :: key
integer,                intent(out) :: nlevels
real(r8),               intent(out) :: mystery
real(r8),               intent(out) :: psurf
real(r8), dimension(:), intent(out) :: profile
integer,                intent(out) :: retrieval_number

if ( .not. module_initialized ) call initialize_module

! Make sure the desired key is within the length of the metadata arrays.
call key_within_range(key,'get_MOPITT_metadata')

profile(:) = MISSING_R8

nlevels            = observation_metadata(key)%nlevels
mystery            = observation_metadata(key)%mystery
psurf              = observation_metadata(key)%psurf
profile(1:nlevels) = observation_metadata(key)%profile(1:nlevels)
retrieval_number   = observation_metadata(key)%retrieval_number

end subroutine get_MOPITT_metadata



 subroutine read_MOPITT_metadata(key,       obsID, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_MOPITT_metadata(obs_def%key, key, ifile, fform)
!
! This routine reads the metadata for neutron intensity observations.
!
integer,          intent(out)          :: key    ! index into local metadata
integer,          intent(in)           :: obsID
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

! temp variables
logical           :: is_asciifile
integer           :: ierr
character(len=6)  :: header
integer           :: oldkey

integer  :: nlevels
real(r8) :: mystery
real(r8) :: psurf
real(r8), dimension(MAXNLEVELS) :: profile
integer  :: retrieval_number

if ( .not. module_initialized ) call initialize_module

is_asciifile = ascii_file_format(fform)

write(string2,*)'observation #',obsID

if ( is_asciifile ) then

   read(ifile, *, iostat=ierr) nlevels
   call check_iostat(ierr,'read_MOPITT_metadata','nlevels ',string2)
   read(ifile, *, iostat=ierr) mystery
   call check_iostat(ierr,'read_MOPITT_metadata','mystery ',string2)
   read(ifile, *, iostat=ierr) psurf
   call check_iostat(ierr,'read_MOPITT_metadata','psurf ',string2)

!   read(ifile, *, iostat=ierr) profile(1), profile(2), profile(3)
!   call check_iostat(ierr,'read_MOPITT_metadata','profile(1:3) ',string2)
!   read(ifile, *, iostat=ierr) profile(4), profile(5), profile(6)
!   call check_iostat(ierr,'read_MOPITT_metadata','profile(4:6) ',string2)
!   read(ifile, *, iostat=ierr) profile(7), profile(8), profile(9)
!   call check_iostat(ierr,'read_MOPITT_metadata','profile(7:9) ',string2)
!   read(ifile, *, iostat=ierr) profile(10)
!   call check_iostat(ierr,'read_MOPITT_metadata','profile(10) ',string2)

   read(ifile, *, iostat=ierr) (profile(i),i=1,nlevels)
   call check_iostat(ierr,'read_MOPITT_metadata','profile',string2)

   read(ifile, *, iostat=ierr) retrieval_number
   call check_iostat(ierr,'read_MOPITT_metadata','retrieval_number ',string2)

else

   read(ifile, iostat=ierr) nlevels
   call  check_iostat(ierr,'read_MOPITT_metadata','nlevels',string2)
   read(ifile, iostat=ierr) mystery
   call  check_iostat(ierr,'read_MOPITT_metadata','mystery ',string2)
   read(ifile, iostat=ierr) psurf
   call  check_iostat(ierr,'read_MOPITT_metadata','psurf ',string2)

!   read(ifile, iostat=ierr) profile(1), profile(2), profile(3)
!   call  check_iostat(ierr,'read_MOPITT_metadata','profile(1:3)',string2)
!   read(ifile, iostat=ierr) profile(4), profile(5), profile(6)
!   call  check_iostat(ierr,'read_MOPITT_metadata','profile(4:6)',string2)
!   read(ifile, iostat=ierr) profile(7), profile(8), profile(9)
!   call  check_iostat(ierr,'read_MOPITT_metadata','profile(7:9)',string2)
!   read(ifile, iostat=ierr) profile(10)
!   call  check_iostat(ierr,'read_MOPITT_metadata','profile(10)',string2)

   read(ifile, iostat=ierr) (profile(i),i=1,nlevels)
   call check_iostat(ierr,'read_MOPITT_metadata','profile',string2)

   read(ifile, iostat=ierr) retrieval_number
   call  check_iostat(ierr,'read_MOPITT_metadata','retrieval_number ',string2)
endif

! The oldkey is thrown away.

! Store the metadata in module storage and record the new length of the metadata arrays.
call set_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number)

! The new 'key' is returned.

end subroutine read_MOPITT_metadata



 subroutine write_MOPITT_metadata(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_MOPITT_metadata(key, ifile, fform)
!
! writes the metadata for neutron intensity observations.

integer,           intent(in)           :: key
integer,           intent(in)           :: ifile
character(len=*),  intent(in), optional :: fform

logical  :: is_asciifile
integer  :: nlevels
real(r8) :: mystery
real(r8) :: psurf
real(r8), dimension(MAXNLEVELS) :: profile
integer  :: retrieval_number

if ( .not. module_initialized ) call initialize_module

! given the index into the local metadata arrays - retrieve
! the metadata for this particular observation.

call get_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )

is_asciifile = ascii_file_format(fform)

if (is_asciifile) then
   write(ifile, *) nlevels
   write(ifile, *) mystery
   write(ifile, *) psurf
!   write(ifile, *) profile(1), profile(2), profile(3)
!   write(ifile, *) profile(4), profile(5), profile(6)
!   write(ifile, *) profile(7), profile(8), profile(9)
!   write(ifile, *) profile(10)
   write(ifile, *) (profile(i),i=1,nlevels)
   write(ifile, *) retrieval_number
else
   write(ifile   ) nlevels
   write(ifile   ) mystery
   write(ifile   ) psurf
!   write(ifile   ) profile(1), profile(2), profile(3)
!   write(ifile   ) profile(4), profile(5), profile(6)
!   write(ifile   ) profile(7), profile(8), profile(9)
!   write(ifile   ) profile(10)
   write(ifile   ) (profile(i),i=1,nlevels)
   write(ifile   ) retrieval_number
endif

end subroutine write_MOPITT_metadata



!  subroutine interactive_MOPITT_metadata(key)
! !----------------------------------------------------------------------
! !subroutine interactive_MOPITT_metadata(key)
! !
! integer, intent(out) :: key
!
! real(r8) :: nlevels, mystery, psurf, profile, retrieval_number
!
! if ( .not. module_initialized ) call initialize_module
!
! ! Prompt for input for the required metadata
!
! nlevels          = interactive('"nlevels"   number of levels '                ,minvalue = 0.0_r8)
! mystery          = interactive('"mystery"   dunno [m^3/m^3]'               ,minvalue = 0.0_r8)
! psurf            = interactive('"psurf"     surface pressure ]'         ,minvalue = 0.0_r8)
! profile          = interactive('"profile"   10 numbers, please ',minvalue=0.0_r8)
! retrieval_number = interactive('"retrieval_number"  inteer   ]'  ,minvalue = 0.0_r8)
!
! call set_MOPITT_metadata(key, nlevels, mystery, psurf, profile, retrieval_number )
!
! end subroutine interactive_MOPITT_metadata



function interactive(str1,minvalue,maxvalue)
real(r8)                       :: interactive
character(len=*),   intent(in) :: str1
real(r8), optional, intent(in) :: minvalue
real(r8), optional, intent(in) :: maxvalue

integer :: i

interactive = MISSING_R8

! Prompt with a minimum amount of error checking

if     (present(minvalue) .and. present(maxvalue)) then

   interactive = minvalue - 1.0_r8
   MINMAXLOOP : do i = 1,10
      if ((interactive >= minvalue) .and. (interactive <= maxvalue)) exit MINMAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MINMAXLOOP

elseif (present(minvalue)) then

   interactive = minvalue - 1.0_r8
   MINLOOP : do i=1,10
      if (interactive >= minvalue) exit MINLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MINLOOP

elseif (present(maxvalue)) then

   interactive = maxvalue + 1.0_r8
   MAXLOOP : do i=1,10
      if (interactive <= maxvalue) exit MAXLOOP
      write(*, *) 'Enter '//str1
      read( *, *) interactive
   end do MAXLOOP

else ! anything goes ... cannot check
      write(*, *) 'Enter '//str1
      read( *, *) interactive
endif

end function interactive




 subroutine check_iostat(istat, routine, varname, msgstring)
!----------------------------------------------------------------------
!

integer,          intent(in) :: istat
character(len=*), intent(in) :: routine
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: msgstring

if ( istat /= 0 ) then
   write(string1,*)'istat should be 0 but is ',istat,' for '//varname
   call error_handler(E_ERR, routine, string1, source, revision, revdate, text2=msgstring)
end if

end subroutine check_iostat



subroutine key_within_range(key, routine)
!----------------------------------------------------------------------
! Make sure we are addressing within the metadata arrays

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

! fine -- no problem.
if ((key > 0) .and. (key <= MOPITTkey)) return

! Bad news. Tell the user.
write(string1, *) 'key (',key,') not within known range ( 1,', MOPITTkey,')'
call error_handler(E_ERR,routine,string1,source,revision,revdate)

end subroutine key_within_range



subroutine grow_metadata(key, routine)
!----------------------------------------------------------------------
! If the allocatable metadata arrays are not big enough ... try again

integer,          intent(in) :: key
character(len=*), intent(in) :: routine

integer :: orglength
type(site_metadata), allocatable, dimension(:) :: safe_metadata

! fine -- no problem.
if ((key > 0) .and. (key <= MAXMOPITTkey)) return

orglength    =     MAXMOPITTkey
MAXMOPITTkey = 2 * orglength

! Check for some error conditions.
if (key < 1) then
   write(string1, *) 'key (',key,') must be >= 1'
   call error_handler(E_ERR,routine,string1,source,revision,revdate)
elseif (key >= 2*MAXMOPITTkey) then
   write(string1, *) 'key (',key,') really unexpected.'
   write(string2, *) 'doubling storage will not help.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate, &
                      text2=string2)
endif

! News. Tell the user we are increasing storage.
write(string1, *) 'key (',key,') exceeds Nmax_neutron_intensity (',orglength,')'
write(string2, *) 'Increasing Nmax_neutron_intensity to ',MAXMOPITTkey
call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)

allocate(safe_metadata(orglength))
safe_metadata(:) = observation_metadata(:)

deallocate(observation_metadata)
  allocate(observation_metadata(MAXMOPITTkey))

observation_metadata(1:orglength)              = safe_metadata(:)
observation_metadata(orglength+1:MAXMOPITTkey) = missing_metadata

deallocate(safe_metadata)

end subroutine grow_metadata



end module obs_def_MOPITT_mod

! END DART PREPROCESS MODULE CODE

