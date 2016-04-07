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

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISSURFACE

use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_CO, KIND_SURFACE_PRESSURE

implicit none
private
public :: write_mopitt_co, read_mopitt_co, interactive_mopitt_co, &
          get_expected_mopitt_co, set_obs_def_mopitt_co

! Storage for the special information required for observations of this type
integer, parameter               :: max_mopitt_co_obs = 10000000
integer, parameter               :: mopitt_dim = 10
integer                          :: num_mopitt_co_obs = 0
real(r8), dimension(max_mopitt_co_obs,10) :: avg_kernel
real(r8), dimension(max_mopitt_co_obs)	 :: mopitt_prior
real(r8)   :: mopitt_pressure(mopitt_dim) =(/ &
                              100000.,90000.,80000.,70000.,60000.,50000.,40000.,30000.,20000.,1000. /)
real(r8)   :: mopitt_pressure_mid(mopitt_dim) =(/ &
                              100000.,85000.,75000.,65000.,55000.,45000.,35000.,25000.,15000.,7500. /)
real(r8), dimension(max_mopitt_co_obs)	 :: mopitt_psurf	
integer,  dimension(max_mopitt_co_obs)   :: mopitt_nlevels

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer  :: counts1 = 0

contains

!----------------------------------------------------------------------

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module



 subroutine read_mopitt_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_mopitt_co(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32) 		:: fileformat

integer			:: mopitt_nlevels_1
real(r8)			:: mopitt_prior_1
real(r8)			:: mopitt_psurf_1
real(r8), dimension(mopitt_dim)	:: avg_kernels_1
integer 			:: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

avg_kernels_1(:) = 0.0_r8

SELECT CASE (fileformat)

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


 subroutine write_mopitt_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_mopitt_co(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional 	:: fform

character(len=32) 		:: fileformat
real(r8), dimension(mopitt_dim) :: avg_kernels_temp

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
   
avg_kernels_temp=avg_kernel(key,:)

SELECT CASE (fileformat)
   
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


 subroutine interactive_mopitt_co(key)
!----------------------------------------------------------------------
!subroutine interactive_mopitt_co(key)
!
! Initializes the specialized part of a MOPITT observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_mopitt_co_obs >= max_mopitt_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'interactive_mopitt_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_mopitt_co_obs (currently ',max_mopitt_co_obs,')'
   call error_handler(E_ERR,'interactive_mopitt_co',msgstring,source,revision,revdate)
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


 subroutine get_expected_mopitt_co(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_mopitt_co(state, location, key, val, istatus)

    real(r8), intent(in)            :: state(:)
    type(location_type), intent(in) :: location
    integer, intent(in)             :: key
    real(r8), intent(out)           :: val
    integer, intent(out)            :: istatus
!
    integer :: i,kstr
    type(location_type) :: loc2
    real(r8)            :: mloc(3)
    real(r8)	        :: obs_val,wrf_psf,level,missing
    real(r8)            :: co_min,mopitt_prs_mid,mopitt_psf
    real(r8)            :: vert_mode_filt
!
    integer             :: nlevels,nnlevels
    integer             :: iflg
    character(len=129)  :: msgstring
!
! Initialize DART
    if ( .not. module_initialized ) call initialize_module
!
! Initialize variables
    co_min=1.e-4
    missing=-888888.0_r8
!
! Get mopitt data
    nlevels = mopitt_nlevels(key)
    mopitt_psf = mopitt_psurf(key)
!
! Get location infomation
    mloc = get_location(location)
    if (mloc(2)>90.0_r8) then
        mloc(2)=90.0_r8
    elseif (mloc(2)<-90.0_r8) then
        mloc(2)=-90.0_r8
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Try phase space localization 
!    vert_mode_filt=10000.
!    if(mloc(3).le.vert_mode_filt) then
!       istatus=2
!       obs_val=missing
!       return
!    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get wrf surface pressure
    wrf_psf = 0.0_r8
    istatus = 0
    loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
    call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)  
!
! Correct mopitt surface pressure
    if(istatus/=0) then
       write(msgstring, *)'APM NOTICE: MOPITT CO WRF psf is bad ',wrf_psf,istatus
       call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
       obs_val=missing
       return
    endif              
!
    if(mopitt_psf.gt.wrf_psf) then
       if((mopitt_psf-wrf_psf).gt.10000.) then
          write(msgstring, *)'APM: NOTICE - reject MOPITT CO - WRF PSF too large ',mopitt_psf,wrf_psf
          call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
          obs_val=missing
          istatus=2
          return
       else
          write(msgstring, *)'APM: NOTICE correct MOPITT CO psf with WRF psf ',mopitt_psf,wrf_psf
          call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
          mopitt_psf=wrf_psf
       endif
    endif
!
! Find kstr - the surface level index
    kstr=0
    do i=1,mopitt_dim
       if (i.eq.1 .and. mopitt_psf.gt.mopitt_pressure(2)) then
          kstr=i
          exit
       endif
       if (i.ne.1 .and. i.ne.mopitt_dim .and. mopitt_pressure(i).ge.mopitt_psf .and. &
       mopitt_psf.gt.mopitt_pressure(i+1)) then
          kstr=i
          exit   
       endif
    enddo
    if (kstr.eq.0) then
       write(msgstring, *)'APM: ERROR in MOPITT CO obs def kstr=0: mopitt_psf= ',mopitt_psf
       call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
       stop
    elseif (kstr.gt.6) then
       write(msgstring, *)'APM: ERROR MOPITT CO surface pressure is unrealistic: mopitt_psf, wrf_psf= ',mopitt_psf,wrf_psf
       call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
       stop
    endif
!
! Reject ob when number of MOPITT levels from WRF cannot equal actual number of MOPITT levels
    nnlevels=mopitt_dim-kstr+1
    if(nnlevels.ne.nlevels) then
       write(msgstring, *)'APM: NOTICE reject MOPITT CO ob - WRF MOP  levels .ne. MOP levels, nnlvls,nlvls ',nnlevels,nlevels
       call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
       obs_val=missing
       istatus=2
       return
    endif   
!
! Apply MOPITT Averaging kernel A and MOPITT Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
    val = 0.0_r8
    do i=1,nlevels
!
! APM: remove the if test to use layer average data
       if (i.eq.1) then
          mopitt_prs_mid=(mopitt_psf+mopitt_pressure(kstr+i))/2.
          loc2 = set_location(mloc(1),mloc(2),mopitt_prs_mid, VERTISPRESSURE)
       else
          mopitt_prs_mid=mopitt_pressure_mid(kstr+i-1)
          loc2 = set_location(mloc(1),mloc(2),mopitt_prs_mid, VERTISPRESSURE)
       endif
!
! Interpolate WRF CO data to MOPITT pressure level midpoint
       obs_val = 0.0_r8
       istatus = 0
       call interpolate(state, loc2, KIND_CO, obs_val, istatus)  
       if (istatus /= 0) then
          write(msgstring, *)'APM NOTICE: WRF extrapolation needed reject MOPITT CO ob '
          call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
          obs_val = missing 
          return
       endif
!
! Check for WRF CO lower bound
       if (obs_val.lt.co_min) then
          write(msgstring, *)'APM NOTICE: resetting minimum MOPITT CO value '
          call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
          obs_val=co_min
       endif
!
! apply averaging kernel
       val = val + avg_kernel(key,i) * log10(obs_val*1.e-6)  
    enddo
!   val = val + mopitt_prior(key)
!
end subroutine get_expected_mopitt_co
!
!----------------------------------------------------------------------

 subroutine set_obs_def_mopitt_co(key, co_avgker, co_prior, co_psurf, co_nlevels)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,	 	intent(in)	:: key, co_nlevels
real*8,dimension(10),	intent(in)	:: co_avgker	
real*8,			intent(in)	:: co_prior
real*8,			intent(in)	:: co_psurf
character(len=129) 			:: msgstring

if ( .not. module_initialized ) call initialize_module

if(num_mopitt_co_obs >= max_mopitt_co_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_mopitt_co_obs (currently ',max_mopitt_co_obs,')'
   call error_handler(E_ERR,'set_obs_def_mopitt_co',msgstring,source,revision,revdate)
endif

avg_kernel(key,:) 	= co_avgker(:)
mopitt_prior(key)	= co_prior
mopitt_psurf(key)	= co_psurf
mopitt_nlevels(key)     = co_nlevels

end subroutine set_obs_def_mopitt_co


function read_mopitt_prior(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_prior
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_prior
   CASE DEFAULT
      read(ifile, *) read_mopitt_prior
END SELECT

end function read_mopitt_prior

function read_mopitt_nlevels(ifile, fform)

integer,                    intent(in) :: ifile
integer                               :: read_mopitt_nlevels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_nlevels
   CASE DEFAULT
      read(ifile, *) read_mopitt_nlevels
END SELECT

end function read_mopitt_nlevels



subroutine write_mopitt_prior(ifile, mopitt_prior_temp, fform)

integer,                    intent(in) :: ifile
real(r8), 		    intent(in) :: mopitt_prior_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_prior_temp
   CASE DEFAULT
      write(ifile, *) mopitt_prior_temp
END SELECT

end subroutine write_mopitt_prior

subroutine write_mopitt_nlevels(ifile, mopitt_nlevels_temp, fform)

integer,                    intent(in) :: ifile
integer,                    intent(in) :: mopitt_nlevels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_nlevels_temp
   CASE DEFAULT
      write(ifile, *) mopitt_nlevels_temp
END SELECT

end subroutine write_mopitt_nlevels



function read_mopitt_psurf(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_psurf
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_psurf
   CASE DEFAULT
      read(ifile, *) read_mopitt_psurf
END SELECT

end function read_mopitt_psurf

subroutine write_mopitt_psurf(ifile, mopitt_psurf_temp, fform)

integer,                    intent(in) :: ifile
real(r8),		    intent(in) :: mopitt_psurf_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) mopitt_psurf_temp
   CASE DEFAULT
      write(ifile, *) mopitt_psurf_temp
END SELECT

end subroutine write_mopitt_psurf

function read_mopitt_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(10)        :: read_mopitt_avg_kernels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

read_mopitt_avg_kernels(:) = 0.0_r8

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_mopitt_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_mopitt_avg_kernels(1:nlevels)
END SELECT

end function read_mopitt_avg_kernels

subroutine write_mopitt_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)

integer,                    intent(in) :: ifile, nlevels_temp
real(r8), dimension(10), intent(in)  :: avg_kernels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels_temp)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels_temp)
END SELECT

end subroutine write_mopitt_avg_kernels



end module obs_def_mopitt_mod
! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
