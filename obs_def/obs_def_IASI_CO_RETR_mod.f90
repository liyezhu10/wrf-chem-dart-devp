! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! IASI_CO_RETRIEVAL, KIND_CO
! END DART PREPROCESS KIND LIST
!
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_iasi_co_mod, only : write_iasi_co, read_iasi_co, &
!                                  interactive_iasi_co, get_expected_iasi_co, &
!                                  set_obs_def_iasi_co
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(IASI_CO_RETRIEVAL)                                                           
!            call get_expected_iasi_co(state, location, obs_def%key, obs_val, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!
! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call read_iasi_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF
!
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call write_iasi_co(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF
!
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_CO_RETRIEVAL)
!         call interactive_iasi_co(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!
! BEGIN DART PREPROCESS SET_OBS_DEF_IASI_CO
!      case(IASI_CO_RETRIEVAL)
!         call set_obs_def_iasi_co(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_IASI_CO
!
! BEGIN DART PREPROCESS MODULE CODE

module obs_def_iasi_CO_mod

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISSURFACE
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_CO, KIND_SURFACE_PRESSURE

implicit none
private

public :: write_iasi_co, &
          read_iasi_co, &
          interactive_iasi_co, &
          get_expected_iasi_co, &
          set_obs_def_iasi_co

! Storage for the special information required for observations of this type
integer, parameter :: MAX_IASI_CO_OBS = 10000000
integer, parameter :: IASI_DIM = 19
integer, parameter :: IASI_DIMP = 20
integer            :: num_iasi_co_obs = 0
integer            :: counts1 = 0

real(r8), dimension(MAX_IASI_CO_OBS,IASI_DIM)  :: avg_kernel
real(r8), dimension(MAX_IASI_CO_OBS,IASI_DIMP) :: pressure
real(r8), dimension(MAX_IASI_CO_OBS)           :: iasi_prior
real(r8), dimension(MAX_IASI_CO_OBS)           :: iasi_psurf
integer,  dimension(MAX_IASI_CO_OBS)           :: iasi_nlevels
integer,  dimension(MAX_IASI_CO_OBS)           :: iasi_nlevelsp
!
! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
logical, save :: module_initialized = .false.

contains

!----------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)

module_initialized = .true.

end subroutine initialize_module

!----------------------------------------------------------------------
!>

subroutine read_iasi_co(key, ifile, fform)

integer,                    intent(out) :: key
integer,                    intent(in)  :: ifile
character(len=*), optional, intent(in)  :: fform

character(len=32) :: fileformat

integer  :: iasi_nlevels_1
integer  :: iasi_nlevelsp_1
real(r8) :: iasi_prior_1
real(r8) :: iasi_psurf_1
real(r8), dimension(IASI_DIM)  :: avg_kernels_1
real(r8), dimension(IASI_DIMP) :: pressure_1
integer :: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_1(:) = 0.0_r8
pressure_1(:) = 0.0_r8
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   iasi_nlevels_1 = read_iasi_nlevels(ifile, fileformat)
   iasi_nlevelsp_1 = iasi_nlevels_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(1:iasi_nlevels_1) = read_iasi_avg_kernels(ifile, iasi_nlevels_1, fileformat)
   pressure_1(1:iasi_nlevelsp_1) = read_iasi_pressure(ifile, iasi_nlevelsp_1, fileformat)
   read(ifile) keyin
   CASE DEFAULT
   iasi_nlevels_1 = read_iasi_nlevels(ifile, fileformat)
   iasi_nlevelsp_1 = iasi_nlevels_1+1
   iasi_prior_1 = read_iasi_prior(ifile, fileformat)
   iasi_psurf_1 = read_iasi_psurf(ifile, fileformat)
   avg_kernels_1(1:iasi_nlevels_1) = read_iasi_avg_kernels(ifile, iasi_nlevels_1, fileformat)
   pressure_1(1:iasi_nlevelsp_1) = read_iasi_pressure(ifile, iasi_nlevelsp_1, fileformat)
   read(ifile, *) keyin
END SELECT
counts1 = counts1 + 1
key = counts1
call set_obs_def_iasi_co(key, avg_kernels_1, pressure_1, iasi_prior_1, iasi_psurf_1, &
                           iasi_nlevels_1, iasi_nlevelsp_1)
end subroutine read_iasi_co

!----------------------------------------------------------------------
!>

subroutine write_iasi_co(key, ifile, fform)

integer,                    intent(in) :: key
integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform

character(len=32) :: fileformat
real(r8) :: avg_kernels_temp(IASI_DIM)
real(r8) :: pressure_temp(IASI_DIMP)

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
avg_kernels_temp=avg_kernel(key,:)
pressure_temp=pressure(key,:)
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_iasi_nlevels(ifile, iasi_nlevels(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlevels(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevelsp(key), fileformat)
   write(ifile) key
   CASE DEFAULT
   call write_iasi_nlevels(ifile, iasi_nlevels(key), fileformat)
   call write_iasi_prior(ifile, iasi_prior(key), fileformat)
   call write_iasi_psurf(ifile, iasi_psurf(key), fileformat)
   call write_iasi_avg_kernels(ifile, avg_kernels_temp, iasi_nlevels(key), fileformat)
   call write_iasi_pressure(ifile, pressure_temp, iasi_nlevelsp(key), fileformat)
   write(ifile, *) key
END SELECT 

end subroutine write_iasi_co

!----------------------------------------------------------------------
!>

subroutine interactive_iasi_co(key)

! Initializes the specialized part of a IASI observation
! Passes back up the key for this one
integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module
!
! Make sure there's enough space, if not die for now (clean later)
if(num_iasi_co_obs >= MAX_IASI_CO_OBS) then
   ! PUT IN ERROR HANDLER CALL
   write(string1, *)'Not enough space for a iasi CO obs.'
   write(string2, *)'Can only have max_iasi_co obs (currently ',MAX_IASI_CO_OBS,')'
   call error_handler(E_ERR,'interactive_iasi_co',string1,source,revision,revdate,text2=string2)
endif
!
! Increment the index
num_iasi_co_obs = num_iasi_co_obs + 1
key = num_iasi_co_obs
!
! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_iasi_co observation'
write(*, *) 'Input the IASI Prior '
read(*, *) iasi_prior
write(*, *) 'Input IASI Surface Pressure '
read(*, *) iasi_psurf(num_iasi_co_obs)
write(*, *) 'Input the 19 Averaging Kernel Weights '
read(*, *) avg_kernel(num_iasi_co_obs,:)
write(*, *) 'Input the 20 Averaging Pressure Levels '
read(*, *) pressure(num_iasi_co_obs,:)

end subroutine interactive_iasi_co

!----------------------------------------------------------------------
!>

subroutine get_expected_iasi_co(state, location, key, val, istatus)

real(r8),            intent(in)  :: state(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: val
integer,             intent(out) :: istatus

integer :: i,kstr
type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)            :: obs_val,wrf_psf,level
real(r8)            :: co_min,obs_val_min,iasi_prs_mid,iasi_psf,iasi_psf_save,iasi_prs
integer             :: nlevels,nlevelsp,nnlevels
integer             :: iflg

real(r8)            :: vert_mode_filt

! Initialize DART
if ( .not. module_initialized ) call initialize_module

! Initialize variables
co_min=1.e-2
obs_val_min=4.e-2

! Get iasi data
nlevels       = iasi_nlevels(key)
nlevelsp      = iasi_nlevelsp(key)
iasi_psf      = iasi_psurf(key)
iasi_psf_save = iasi_psurf(key)
!
! Get location information
mloc = get_location(location)
if (mloc(2)>90.0_r8) then
    mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
    mloc(2)=-90.0_r8
endif

! Get wrf surface pressure
wrf_psf = 0.0_r8
istatus = 0
loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)  

! Correct iasi surface pressure
if(istatus/=0) then
   write(string1, *)'APM NOTICE: IASI CO WRF psf is bad ',wrf_psf,istatus
   call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
   obs_val=missing_r8
   return
endif

if(iasi_psf.gt.wrf_psf) then
   if((iasi_psf-wrf_psf).gt.10000.) then
      write(string1, *)'APM: NOTICE - reject IASI CO - WRF PSF too large ', &
      iasi_psf,wrf_psf
      call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
      obs_val=missing_r8
      istatus=2
      return
   else
!      write(string1, *)'APM: NOTICE correct IASI CO psf with WRF psf ',iasi_psf,wrf_psf
!      call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
      iasi_psf=wrf_psf
   endif
endif

! Find kstr - the surface level index
kstr=0
do i=1,iasi_dim
   if (i.eq.1 .and. iasi_psf.gt.pressure(key,2)) then
      kstr=i
      exit
   endif
   if (i.ne.1 .and. i.ne.iasi_dim .and. pressure(key,i).ge.iasi_psf.and. &
      iasi_psf.gt.pressure(key,i+1)) then
      kstr=i
      exit   
   endif
enddo
if (kstr.eq.0) then
   write(string1, *)'APM: ERROR in IASI CO obs def kstr=0: iasi_psf=',iasi_psf
   call error_handler(E_ERR,'set_obs_def_iasi_co',string1,source,revision,revdate)
elseif (kstr.gt.7) then
   write(string1, *)'APM: ERROR IASI CO surface pressure is unrealistic: iasi_psf, wrf_psf= ',iasi_psf,wrf_psf
   call error_handler(E_ERR,'set_obs_def_iasi_co',string1,source,revision,revdate)
endif

! Reject ob when number of IASI levels from WRF cannot equal actual number of IASI levels
nnlevels=iasi_dim-kstr+1
if(nnlevels.ne.nlevels) then
   write(string1, *)'APM: NOTICE reject IASI CO ob - WRF IASI lvls .ne. IASI levels, nnlvls,nlvls ', &
   nnlevels,nlevels
   call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
   obs_val=missing_r8
   istatus=2
   return
endif   

! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
val = 0.0_r8
do i=1,nlevels

   ! APM: remove the if test to use layer average data
   if (i.eq.1) then
      iasi_prs_mid=(iasi_psf+pressure(key,kstr+i))/2.
      loc2 = set_location(mloc(1),mloc(2),iasi_prs_mid, VERTISPRESSURE)
   else
      iasi_prs_mid=(pressure(key,kstr+i-1)+pressure(key,kstr+i))/2.
      loc2 = set_location(mloc(1),mloc(2),iasi_prs_mid, VERTISPRESSURE)
   endif

   ! Interpolate WRF CO data to IASI pressure level midpoint
   obs_val = 0.0_r8
   istatus = 0
   call interpolate(state, loc2, KIND_CO, obs_val, istatus)  
   if (istatus /= 0) then
      write(string1, *)'APM NOTICE: WRF extrapolation needed reject IASI CO ob '
      call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
      obs_val = missing_r8
      return
   endif

   ! Check for min values
   if (obs_val.lt.co_min) then
      write(string1, *)'APM NOTICE: resetting minimum IASI CO value ', &
      obs_val,co_min,obs_val_min
      call error_handler(E_MSG,'set_obs_def_iasi_co',string1,source,revision,revdate)
      obs_val=obs_val_min
   endif

   ! apply averaging kernel
   val = val + avg_kernel(key,i) * obs_val*1.e3  

enddo

val = val + iasi_prior(key)

end subroutine get_expected_iasi_co

!----------------------------------------------------------------------
!>

subroutine set_obs_def_iasi_co(key, co_avgker, co_press, co_prior, co_psurf, co_nlevels, co_nlevelsp)

! Allows passing of obs_def special information 
integer,  intent(in) :: key, co_nlevels, co_nlevelsp
real(r8), intent(in) :: co_avgker(IASI_DIM)
real(r8), intent(in) :: co_press(IASI_DIMP)
real(r8), intent(in) :: co_prior
real(r8), intent(in) :: co_psurf

if ( .not. module_initialized ) call initialize_module

if(num_iasi_co_obs >= MAX_IASI_CO_OBS) then
   write(string1, *)'Not enough space for a iasi CO obs.'
   write(string2, *)'Can only have MAX_IASI_CO_OBS (currently ',MAX_IASI_CO_OBS,')'
   call error_handler(E_ERR,'set_obs_def_iasi_co',string1,source,revision,revdate,text2=string2)
endif

avg_kernel(key,:) = co_avgker(:)
pressure(  key,:) = co_press(:)
iasi_prior(key)   = co_prior
iasi_psurf(key)   = co_psurf
iasi_nlevels(key) = co_nlevels
iasi_nlevelsp(key)= co_nlevelsp

end subroutine set_obs_def_iasi_co

!----------------------------------------------------------------------
!>

function read_iasi_prior(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_prior

character(len=32)  :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior
   CASE DEFAULT
      read(ifile, *) read_iasi_prior
END SELECT

end function read_iasi_prior

!----------------------------------------------------------------------
!>

function read_iasi_nlevels(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
integer                                :: read_iasi_nlevels

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlevels
   CASE DEFAULT
      read(ifile, *) read_iasi_nlevels
END SELECT

end function read_iasi_nlevels

!----------------------------------------------------------------------
!>

function read_iasi_nlevelsp(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
integer                                :: read_iasi_nlevelsp

character(len=32)  :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_nlevelsp
   CASE DEFAULT
      read(ifile, *) read_iasi_nlevelsp
END SELECT

end function read_iasi_nlevelsp

!----------------------------------------------------------------------
!>

subroutine write_iasi_prior(ifile, iasi_prior_temp, fform)

integer,          intent(in) :: ifile
real(r8),         intent(in) :: iasi_prior_temp
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_prior_temp
   CASE DEFAULT
      write(ifile, *) iasi_prior_temp
END SELECT

end subroutine write_iasi_prior

!----------------------------------------------------------------------
!>

subroutine write_iasi_nlevels(ifile, iasi_nlevels_temp, fform)

integer,          intent(in) :: ifile
integer,          intent(in) :: iasi_nlevels_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlevels_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlevels_temp
END SELECT

end subroutine write_iasi_nlevels

!----------------------------------------------------------------------
!>

subroutine write_iasi_nlevelsp(ifile, iasi_nlevelsp_temp, fform)

integer,          intent(in) :: ifile
integer,          intent(in) :: iasi_nlevelsp_temp
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_nlevelsp_temp
   CASE DEFAULT
      write(ifile, *) iasi_nlevelsp_temp
END SELECT

end subroutine write_iasi_nlevelsp

!----------------------------------------------------------------------
!>

function read_iasi_psurf(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_psurf

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_psurf
   CASE DEFAULT
      read(ifile, *) read_iasi_psurf
END SELECT

end function read_iasi_psurf

!----------------------------------------------------------------------
!>

subroutine write_iasi_psurf(ifile, iasi_psurf_temp, fform)

integer,          intent(in) :: ifile
real(r8),         intent(in) :: iasi_psurf_temp
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) iasi_psurf_temp
   CASE DEFAULT
      write(ifile, *) iasi_psurf_temp
END SELECT

end subroutine write_iasi_psurf

!----------------------------------------------------------------------
!>

function read_iasi_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_avg_kernels(IASI_DIM)

character(len=32) :: fileformat

read_iasi_avg_kernels(:) = 0.0_r8

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_avg_kernels(1:nlevels)
END SELECT

end function read_iasi_avg_kernels

!----------------------------------------------------------------------
!>

subroutine write_iasi_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)

integer,          intent(in) :: ifile, nlevels_temp
real(r8),         intent(in) :: avg_kernels_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels_temp)
END SELECT

end subroutine write_iasi_avg_kernels

!----------------------------------------------------------------------
!>

function read_iasi_pressure(ifile, nlevelsp, fform)

integer,                    intent(in) :: ifile, nlevelsp
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_pressure(IASI_DIMP)

character(len=32) :: fileformat

read_iasi_pressure(:) = 0.0_r8

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_pressure(1:nlevelsp)
   CASE DEFAULT
      read(ifile, *) read_iasi_pressure(1:nlevelsp)
END SELECT

end function read_iasi_pressure

!----------------------------------------------------------------------
!>

subroutine write_iasi_pressure(ifile, pressure_temp, nlevelsp_temp, fform)

integer,          intent(in) :: ifile, nlevelsp_temp
real(r8),         intent(in) :: pressure_temp(IASI_DIMP)
character(len=*), intent(in) :: fform

character(len=32)  :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) pressure_temp(1:nlevelsp_temp)
   CASE DEFAULT
      write(ifile, *) pressure_temp(1:nlevelsp_temp)
END SELECT

end subroutine write_iasi_pressure

!----------------------------------------------------------------------
!>

end module obs_def_iasi_CO_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
