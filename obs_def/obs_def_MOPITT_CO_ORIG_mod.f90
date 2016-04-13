! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! MOPITT_CO_RETRIEVAL, KIND_CO
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_mopitt_mod, only : write_mopitt_co, read_mopitt_co, &
!                                  interactive_mopitt_co, get_expected_mopitt_co, &
!                                  set_obs_def_mopitt_co
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(MOPITT_CO_RETRIEVAL)
!            call get_expected_mopitt_co(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL)
!         call read_mopitt_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL)
!         call write_mopitt_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(MOPITT_CO_RETRIEVAL)
!         call interactive_mopitt_co(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_MOPITT_CO
!      case(MOPITT_CO_RETRIEVAL)
!         call set_obs_def_mopitt_co(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_MOPITT_CO

! BEGIN DART PREPROCESS MODULE CODE

module obs_def_mopitt_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISLEVEL

use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_CO

implicit none
private

public :: write_mopitt_co, &
          read_mopitt_co, &
          interactive_mopitt_co, &
          get_expected_mopitt_co, &
          set_obs_def_mopitt_co

! Storage for the special information required for observations of this type
integer, parameter :: MAX_MOPITT_CO_OBS = 10000000
integer, parameter :: MOPITT_DIM = 10
integer            :: num_mopitt_co_obs = 0
integer            :: counts1 = 0

real(r8) :: mopitt_prior(  MAX_MOPITT_CO_OBS)
real(r8) :: mopitt_psurf(  MAX_MOPITT_CO_OBS)
integer  :: mopitt_nlevels(MAX_MOPITT_CO_OBS)
real(r8) :: avg_kernel(    MAX_MOPITT_CO_OBS,MOPITT_DIM)

real(r8) :: mopitt_pressure(MOPITT_DIM) = (/ &
                  95000.,90000.,80000.,70000.,60000.,50000.,40000.,30000.,20000.,10000. /)

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
logical, save      :: module_initialized = .false.

contains

!----------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------
!>

subroutine read_mopitt_co(key, ifile, fform)

integer,                    intent(out) :: key
integer,                    intent(in)  :: ifile
character(len=*), optional, intent(in)  :: fform

character(len=32) :: fileformat

integer :: mopitt_nlevels_1
real(r8):: mopitt_prior_1
real(r8):: mopitt_psurf_1
real(r8), dimension(MOPITT_DIM):: avg_kernels_1
integer :: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = adjustl(fform)

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

avg_kernels_1(:) = 0.0_r8

SELECT CASE (trim(fileformat))

   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   mopitt_nlevels_1 = read_mopitt_nlevels(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlevels_1)  = read_mopitt_avg_kernels(ifile, mopitt_nlevels_1, fileformat)
   read(ifile) keyin

   CASE DEFAULT
   mopitt_nlevels_1 = read_mopitt_nlevels(ifile, fileformat)
   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   mopitt_psurf_1 = read_mopitt_psurf(ifile, fileformat)
   avg_kernels_1(1:mopitt_nlevels_1)  = read_mopitt_avg_kernels(ifile, mopitt_nlevels_1, fileformat)
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_mopitt_co(key, avg_kernels_1, mopitt_prior_1, mopitt_psurf_1, &
                           mopitt_nlevels_1)

end subroutine read_mopitt_co

!----------------------------------------------------------------------
!>

subroutine write_mopitt_co(key, ifile, fform)

integer,                    intent(in) :: key
integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform

character(len=32) :: fileformat

real(r8) :: avg_kernels_temp(MOPITT_DIM)

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = adjustl(fform)

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

avg_kernels_temp=avg_kernel(key,:)

SELECT CASE (trim(fileformat))

   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_mopitt_nlevels(ifile, mopitt_nlevels(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlevels(key), fileformat)
   write(ifile) key

   CASE DEFAULT
   call write_mopitt_nlevels(ifile, mopitt_nlevels(key), fileformat)
   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_mopitt_psurf(ifile, mopitt_psurf(key), fileformat)
   call write_mopitt_avg_kernels(ifile, avg_kernels_temp, mopitt_nlevels(key), fileformat)
   write(ifile, *) key
END SELECT

end subroutine write_mopitt_co

!----------------------------------------------------------------------
!>

subroutine interactive_mopitt_co(key)

! Initializes the specialized part of a MOPITT observation
! Passes back up the key for this one

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   write(string1, *)'Not enough space for a mopitt CO obs.'
   write(string2, *)'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'interactive_mopitt_co',string1,source,revision,revdate,text2=string2)
endif

! Increment the index
num_mopitt_co_obs = num_mopitt_co_obs + 1
key = num_mopitt_co_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_mopitt_co observation'
write(*, *) 'Input the MOPITT Prior '
read(*, *) mopitt_prior
write(*, *) 'Input MOPITT Surface Pressure '
read(*, *) mopitt_psurf(num_mopitt_co_obs)
write(*, *) 'Input the 10 Averaging Kernel Weights '
read(*, *) avg_kernel(num_mopitt_co_obs,:)

end subroutine interactive_mopitt_co

!----------------------------------------------------------------------
!>

subroutine get_expected_mopitt_co(state, location, key, val, istatus)

real(r8),            intent(in)  :: state(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: val
integer,             intent(out) :: istatus

integer :: i
type(location_type) :: loc2
real(r8) :: mloc(3)
real(r8) :: obs_val, level

integer  :: nlevels

if ( .not. module_initialized ) call initialize_module

mloc = get_location(location)
! Apply MOPITT Averaging kernel A and MOPITT Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector

val = 0.0_r8
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif

mopitt_pressure(1)=mopitt_psurf(key)
nlevels = mopitt_nlevels(key)
level   = 1.0_r8

do i=1,nlevels
   if (i == 1) then
   loc2 = set_location(mloc(1),mloc(2),level, VERTISLEVEL)
   else
   loc2 = set_location(mloc(1),mloc(2),mopitt_pressure(i), VERTISPRESSURE)
   endif
   obs_val = 0.0_r8
   istatus = 0

   call interpolate(state, loc2, KIND_CO, obs_val, istatus)

   !print *, 'AFAJ ',istatus, obs_val
   if (istatus /= 0) then
      val = 0
      return
   endif
!   if (avg_kernel(key,i)<0d0) then
!      avg_kernel(key,i)=0d0
!   endif
   val = val + avg_kernel(key,i) * (obs_val)
enddo
val = val + mopitt_prior(key)
!print *, val
!print *,'AFAJ DEBUG ', val
!stop

end subroutine get_expected_mopitt_co

!----------------------------------------------------------------------
!>

subroutine set_obs_def_mopitt_co(key, co_avgker, co_prior, co_psurf, co_nlevels)

! Allows passing of obs_def special information

integer,  intent(in):: key, co_nlevels
real(r8), intent(in):: co_avgker(MOPITT_DIM)
real(r8), intent(in):: co_prior
real(r8), intent(in):: co_psurf

if ( .not. module_initialized ) call initialize_module

if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   write(string1,*) 'Not enough space for a mopitt CO obs.'
   write(string2,*) 'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'set_obs_def_mopitt_co',string1,source,revision,revdate, text2=string2)
endif

avg_kernel(key,:)   = co_avgker(:)
mopitt_prior(key)   = co_prior
mopitt_psurf(key)   = co_psurf
mopitt_nlevels(key) = co_nlevels

end subroutine set_obs_def_mopitt_co

!----------------------------------------------------------------------
!>

function read_mopitt_prior(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_mopitt_prior

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_prior
   CASE DEFAULT
      read(ifile, *) read_mopitt_prior
END SELECT

end function read_mopitt_prior

!----------------------------------------------------------------------
!>

function read_mopitt_nlevels(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
integer                                :: read_mopitt_nlevels

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_nlevels
   CASE DEFAULT
      read(ifile, *) read_mopitt_nlevels
END SELECT

end function read_mopitt_nlevels

!----------------------------------------------------------------------
!>

subroutine write_mopitt_prior(ifile, mopitt_prior_temp, fform)

integer,           intent(in) :: ifile
real(r8),          intent(in) :: mopitt_prior_temp
character(len=32), intent(in) :: fform

character(len=32) :: fileformat

fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_prior_temp
   CASE DEFAULT
      write(ifile, *) mopitt_prior_temp
END SELECT

end subroutine write_mopitt_prior

!----------------------------------------------------------------------
!>

subroutine write_mopitt_nlevels(ifile, mopitt_nlevels_temp, fform)

integer,          intent(in) :: ifile
integer,          intent(in) :: mopitt_nlevels_temp
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_nlevels_temp
   CASE DEFAULT
      write(ifile, *) mopitt_nlevels_temp
END SELECT

end subroutine write_mopitt_nlevels

!----------------------------------------------------------------------
!>

function read_mopitt_psurf(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_mopitt_psurf

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_psurf
   CASE DEFAULT
      read(ifile, *) read_mopitt_psurf
END SELECT

end function read_mopitt_psurf

!----------------------------------------------------------------------
!>

subroutine write_mopitt_psurf(ifile, mopitt_psurf_temp, fform)

integer,          intent(in) :: ifile
real(r8),         intent(in) :: mopitt_psurf_temp
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_psurf_temp
   CASE DEFAULT
      write(ifile, *) mopitt_psurf_temp
END SELECT

end subroutine write_mopitt_psurf

!----------------------------------------------------------------------
!>

function read_mopitt_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_mopitt_avg_kernels(MOPITT_DIM)

character(len=32)  :: fileformat

read_mopitt_avg_kernels(:) = 0.0_r8

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_mopitt_avg_kernels(1:nlevels)
END SELECT

end function read_mopitt_avg_kernels

!----------------------------------------------------------------------
!>

subroutine write_mopitt_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)

integer,          intent(in) :: ifile, nlevels_temp
real(r8),         intent(in) :: avg_kernels_temp(MOPITT_DIM)
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat

fileformat = adjustl(fform)

SELECT CASE (trim(fileformat))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels_temp)
END SELECT

end subroutine write_mopitt_avg_kernels

!----------------------------------------------------------------------
!>

end module obs_def_mopitt_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
