! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! Oct 16, 2014 LXL: 
! This is the version fixing the "wrf has no level higher than 6000Pa" error by 
! only sum the subcolumn under 60000Pa level/starting at 120000Pa.
! Due to get_model_pressure_profile returns top press(>120000) somehow, decide
! to sum starting 200000 Pa by changing "nlevels-9" to "nlevels-10"

! BEGIN DART PREPROCESS KIND LIST
! OMI_NO2_COLUMN, KIND_NO2
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_omi_mod, only : write_omi_no2, read_omi_no2, &
!                               interactive_omi_no2, get_expected_omi_no2, &
!                               set_obs_def_omi_no2
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(OMI_NO2_COLUMN)                                                           
!            call get_expected_omi_no2(state, location, obs_def%key, obs_val, istatus)  
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(OMI_NO2_COLUMN)
!         call read_omi_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(OMI_NO2_COLUMN)
!         call write_omi_no2(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(OMI_NO2_COLUMN)
!         call interactive_omi_no2(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_OMI_NO2
!      case(OMI_NO2_COLUMN)
!         call set_obs_def_omi_no2(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF


! BEGIN DART PREPROCESS MODULE CODE
module obs_def_omi_mod

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISSURFACE, VERTISLEVEL

use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_NO2, KIND_SURFACE_PRESSURE

implicit none

public :: write_omi_no2, read_omi_no2, interactive_omi_no2, &
          get_expected_omi_no2, set_obs_def_omi_no2

! Storage for the special information required for observations of this type
integer, parameter               :: max_omi_no2_obs = 10000000
integer, parameter               :: omi_dim = 35
integer                          :: num_omi_no2_obs = 0
real(r8), dimension(max_omi_no2_obs,35) :: avg_kernel
! lxl:real(r8), dimension(max_omi_no2_obs)	 :: mopitt_prior
real(r8)   :: omi_pressure(omi_dim) =(/ &
        102000., 101000., 100000., 99000., 97500., 96000., 94500., &
         92500.,  90000.,  87500., 85000., 82500., 80000., 77000., &
         74000.,  70000.,  66000., 61000., 56000., 50000., 45000., &
         40000.,  35000.,  28000., 20000., 12000.,  6000.,  3500., &
          2000.,   1200.,    800.,   500.,   300.,   150.,     80. /)    
real(r8), dimension(max_omi_no2_obs)   :: omi_psurf	
real(r8), dimension(max_omi_no2_obs)   :: omi_ptrop	
integer,  dimension(max_omi_no2_obs)   :: omi_nlevels

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
!----------------------------------------------------------------------------
! subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

subroutine read_omi_no2(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_omi_no2(key, ifile, fform)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional    :: fform
character(len=32) 		:: fileformat

integer		        	:: omi_nlevels_1
! lxl:real(r8)			:: mopitt_prior_1
real(r8)			:: omi_psurf_1
real(r8)			:: omi_ptrop_1
real(r8), dimension(omi_dim)	:: avg_kernels_1
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
   omi_nlevels_1 = read_omi_nlevels(ifile, fileformat)
!lxl   omi_prior_1 = read_mopitt_prior(ifile, fileformat)
   omi_psurf_1 = read_omi_psurf(ifile, fileformat)
   omi_ptrop_1 = read_omi_ptrop(ifile, fileformat)
   avg_kernels_1(1:omi_nlevels_1)  = read_omi_avg_kernels(ifile, omi_nlevels_1, fileformat)
   read(ifile) keyin

   CASE DEFAULT
   omi_nlevels_1 = read_omi_nlevels(ifile, fileformat)
!lxl   mopitt_prior_1 = read_mopitt_prior(ifile, fileformat)
   omi_psurf_1 = read_omi_psurf(ifile, fileformat)
   omi_ptrop_1 = read_omi_ptrop(ifile, fileformat)
   avg_kernels_1(1:omi_nlevels_1)  = read_omi_avg_kernels(ifile, omi_nlevels_1, fileformat)
   read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key = counts1
call set_obs_def_omi_no2(key, avg_kernels_1, omi_psurf_1, omi_ptrop_1, omi_nlevels_1)

end subroutine read_omi_no2

subroutine write_omi_no2(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine write_omi_no2(key, ifile, fform)

integer, intent(in)             :: key
integer, intent(in)             :: ifile
character(len=*), intent(in), optional 	:: fform

character(len=32) 		:: fileformat
real(r8), dimension(omi_dim) :: avg_kernels_temp

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
   
avg_kernels_temp=avg_kernel(key,:)

SELECT CASE (fileformat)
   
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
   call write_omi_nlevels(ifile, omi_nlevels(key), fileformat)
!lxl   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_omi_psurf(ifile, omi_psurf(key), fileformat)
   call write_omi_ptrop(ifile, omi_ptrop(key), fileformat)
   call write_omi_avg_kernels(ifile, avg_kernels_temp, omi_nlevels(key), fileformat)
   write(ifile) key

   CASE DEFAULT
   call write_omi_nlevels(ifile, omi_nlevels(key), fileformat)
!lxl   call write_mopitt_prior(ifile, mopitt_prior(key), fileformat)
   call write_omi_psurf(ifile, omi_psurf(key), fileformat)
   call write_omi_ptrop(ifile, omi_ptrop(key), fileformat)
   call write_omi_avg_kernels(ifile, avg_kernels_temp, omi_nlevels(key), fileformat)
   write(ifile, *) key
END SELECT 

end subroutine write_omi_no2

subroutine interactive_omi_no2(key)
!----------------------------------------------------------------------
!subroutine interactive_omi_no2(key)
!
! Initializes the specialized part of a MOPITT observation
! Passes back up the key for this one

integer, intent(out) :: key

character(len=129) :: msgstring

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_omi_no2_obs >= max_omi_no2_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a omi NO2 obs.'
   call error_handler(E_MSG,'interactive_omi_no2',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_omi_no2_obs (currently ',max_omi_no2_obs,')'
   call error_handler(E_ERR,'interactive_omi_no2',msgstring,source,revision,revdate)
endif

! Increment the index
num_omi_no2_obs = num_omi_no2_obs + 1
key = num_omi_no2_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_omi_no2 observation'
!write(*, *) 'Input the MOPITT Prior '
!read(*, *) mopitt_prior
write(*, *) 'Input OMI Surface Pressure '
read(*, *) omi_psurf(num_omi_no2_obs)
write(*, *) 'Input OMI Tropopause Pressure '
read(*, *) omi_ptrop(num_omi_no2_obs)
write(*, *) 'Input the 35 Averaging Kernel Weights '
read(*, *) avg_kernel(num_omi_no2_obs,:)

end subroutine interactive_omi_no2

subroutine get_expected_omi_no2(state, location, key, val, istatus)
!----------------------------------------------------------------------
!subroutine get_expected_omi_no2(state, location, key, val, istatus)
!
   real(r8), intent(in)            :: state(:)
   type(location_type), intent(in) :: location
   integer, intent(in)             :: key
   real(r8), intent(out)           :: val
   integer, intent(out)            :: istatus
!
   integer :: i,kstr,kend
   type(location_type) :: loc2
   real(r8)            :: mloc(3)
   real(r8)	       :: obs_val,wrf_psf,level,missing
   real(r8)            :: no2_min,omi_psf,omi_ptrp,omi_psf_save,mg !,mopitt_prs_mid
   real(r8), dimension(omi_dim) :: no2_vmr
!
   integer             :: nlevels,nnlevels
   integer             :: iflg
   character(len=129)  :: msgstring
!
! Initialize DART
   if ( .not. module_initialized ) call initialize_module
!
! 1ppt unit:mixing ratio !1.e-4 LXL: i THINK UNIT IS ppmv
! (mg  = (28.97/6.02E23)*1E-3    *   9.8      *     1E4;  
   no2_min=1.e-6
   missing=-1.2676506e30
   mg=4.716046511627907e-21  
   level   = 1.0_r8
!
! Get omi data
   nlevels = omi_nlevels(key)
   omi_psf = omi_psurf(key)
   omi_ptrp = omi_ptrop(key)
!
! Get location infomation
   mloc = get_location(location)
   if (mloc(2)>90.0_r8) then
       mloc(2)=90.0_r8
   elseif (mloc(2)<-90.0_r8) then
       mloc(2)=-90.0_r8
   endif
!
! Get wrf surface pressure
   wrf_psf = 0.0_r8
   istatus = 0
   loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
   call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)  
!
! Correct omi surface pressure
   print *, 'APM: istatus ',istatus
   print *, 'APM: omi_psf ',omi_psf
   print *, 'APM: wrf_psf ',wrf_psf
   if (omi_psf .gt. wrf_psf) then
      omi_psf=wrf_psf
   endif
   omi_psf_save=omi_psf
!
! Find kstr - the surface level index
   kstr=0
   do i=1,omi_dim
      if (i .eq. 1 .and. omi_psf .gt. omi_pressure(2)) then
         kstr=i
         exit
      endif
      if (i .ne. 1 .and. i .ne. omi_dim .and. omi_psf .le. omi_pressure(i) .and. &
      omi_psf .gt. omi_pressure(i+1)) then
         kstr=i
         exit   
      endif
   enddo
!
! Find kend - index for the first OMI level above omi_ptrop
   kend=0
   do i=1,omi_dim-1
      if (omi_ptrp .lt. omi_pressure(i) .and. omi_ptrp .ge. omi_pressure(i+1)) then
         kend=i+1
         exit   
      endif
   enddo
!
   if (kstr .eq. 0) then
      write(msgstring, *)'APM: ERROR in OMI obs def kstr=0: omi_psf=',omi_psf
      call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      print *, 'APM: omi_psf ',omi_psf
      print *, 'APM: wrf_psf ',wrf_psf
      print *, 'APM: omi_pressure ',omi_pressure


      call abort
   elseif (kstr .gt. 20) then
      write(msgstring, *)'APM: ERROR surface pressure is unrealistic: omi_psf=',omi_psf
      call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      call abort
   endif
   if (kend .eq. 0) then
      write(msgstring, *)'APM: ERROR in OMI obs def kend=0: omi_ptrp=',omi_ptrp
      call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      call abort
   endif
!
! Reject ob when number of OMI levels from WRF cannot equal actual number of OMI levels
   nnlevels=omi_dim-kstr+1-(omi_dim-kend)
   if (nnlevels .ne. nlevels) then
      istatus=2
      obs_val=missing
      write(msgstring, *)'APM: NOTICE reject ob - # of WRF OMI levels .ne. # of OMI levels  '
      call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      return
   endif   
!
! Find the lowest pressure level midpoint
! lxl: omi_prs=(mopitt_psf+omi_pressure(kstr+1))/2.
!
! Migliorini forward operators assimilation A*x_t
! Apply MOPITT Averaging kernel A and MOPITT Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
!
   no2_vmr(:)=0.
   do i=1,nlevels
!
! APM: remove the if test to use layer average data
      if (i .eq. 1) then
         loc2 = set_location(mloc(1),mloc(2),level, VERTISLEVEL)
!         write(msgstring, *)'APM NOTICE: Location for surface '
!         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      else if (i .eq. nlevels) then
         omi_psf=omi_ptrp
         loc2 = set_location(mloc(1),mloc(2),omi_psf, VERTISPRESSURE)
!         write(msgstring, *)'APM NOTICE: Location for tropopause '
!         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      else
         omi_psf=omi_pressure(kstr+i-1)
         loc2 = set_location(mloc(1),mloc(2),omi_psf, VERTISPRESSURE)
!         write(msgstring, *)'APM NOTICE: Location for free atm '
!         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      endif
!
! APM: check whether OMI pressure is less than omi_ptrop
! APM: Note = omi_pressure(nlevels) is first OMI pressure level above omi_ptrp
      if (i .ne. nlevels .and. omi_pressure(kstr+i-1) .lt. omi_ptrp) then
         write(msgstring, *)'APM ERROR: OMI pressure is less than ptrop ',omi_pressure(kstr+i-1),omi_ptrp 
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      endif
!
! Interpolate WRF NO2 data to OMI pressure level midpoint
      obs_val = 0.0_r8
      istatus = 0
      call interpolate(state, loc2, KIND_NO2, obs_val, istatus)  
      if (istatus .ne. 0 .and. istatus .ne. 2) then
         write(msgstring, *)'APM ERROR: istatus,kstr,obs_val ',istatus,kstr,obs_val 
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
         write(msgstring, *)'APM ERROR: wrf_psf,omi_psurf,omi_psf ', &
         wrf_psf,omi_psurf(key),omi_psf
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
         write(msgstring, *)'APM ERROR: i, nlevels ',i,nlevels
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
         call abort
      endif
      if (istatus .eq. 2 .and. kstr+i-1 .ge. nlevels-2) then
         istatus=0
         obs_val=no2_min
!         write(msgstring, *)'APM NOTICE: obs_def_omi - NEED NO2 ABOVE MODEL TOP lev_idx ',kstr+i-1
!         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
      endif           
      if (istatus .eq. 2 .and. i .lt. 3) then
         write(msgstring, *)'APM NOTICE: NO2 MODEL SURF - reject ob ',kstr,kstr+i-1,omi_psf_save,omi_pressure(kstr+i-1)
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
         istatus=2
         obs_val=missing
         return
      endif           
!
! Check for WRF CO lower bound
      if (obs_val .lt. no2_min) then
         write(msgstring, *)'APM NOTICE: RESETTING NO2 ',istatus,kstr+i-1,omi_pressure(kstr+i-1),obs_val
         call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
         obs_val=no2_min
      endif
!
! Convert from ppmv to no unit
      no2_vmr(i)=obs_val*1.e-6
   enddo
   val = 0.0_r8
   do i=1,nlevels-1
!
! apply averaging kernel onto vmr to calculate subcol
      if (i .eq. 1) then 
         val = val + 0.5*(avg_kernel(key,i)*no2_vmr(i)+avg_kernel(key,i+1)*no2_vmr(i+1))*(omi_psf_save-omi_pressure(kstr+i))/mg
      else
         val = val + 0.5*(avg_kernel(key,i)*no2_vmr(i)+avg_kernel(key,i+1)*no2_vmr(i+1))*(omi_pressure(kstr+i-1)-omi_pressure(kstr+i))/mg
      endif
   enddo
!
end subroutine get_expected_omi_no2
!
!----------------------------------------------------------------------

 subroutine set_obs_def_omi_no2(key, no2_avgker, no2_psurf, no2_ptrop, no2_nlevels)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,	 	intent(in)	:: key, no2_nlevels
real*8,dimension(35),	intent(in)	:: no2_avgker	
!real*8,			intent(in)	:: co_prior
real*8,			intent(in)	:: no2_psurf
real*8,			intent(in)	:: no2_ptrop
character(len=129) 			:: msgstring

if ( .not. module_initialized ) call initialize_module

if(num_omi_no2_obs >= max_omi_no2_obs) then
   ! PUT IN ERROR HANDLER CALL
   write(msgstring, *)'Not enough space for a omi NO2 obs.'
   call error_handler(E_MSG,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
   write(msgstring, *)'Can only have max_omi_no2_obs (currently ',max_omi_no2_obs,')'
   call error_handler(E_ERR,'set_obs_def_omi_no2',msgstring,source,revision,revdate)
endif

avg_kernel(key,:) 	= no2_avgker(:)
!mopitt_prior(key)	
omi_psurf(key)	= no2_psurf
omi_ptrop(key)	= no2_ptrop
omi_nlevels(key)     = no2_nlevels

end subroutine set_obs_def_omi_no2

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

function read_omi_nlevels(ifile, fform)

integer,                    intent(in) :: ifile
integer                                :: read_omi_nlevels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_omi_nlevels
   CASE DEFAULT
      read(ifile, *) read_omi_nlevels
END SELECT

end function read_omi_nlevels

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

subroutine write_omi_nlevels(ifile, omi_nlevels_temp, fform)

integer,                    intent(in) :: ifile
integer,                    intent(in) :: omi_nlevels_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) omi_nlevels_temp
   CASE DEFAULT
      write(ifile, *) omi_nlevels_temp
END SELECT

end subroutine write_omi_nlevels

function read_omi_psurf(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_omi_psurf
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_omi_psurf
   CASE DEFAULT
      read(ifile, *) read_omi_psurf
END SELECT

end function read_omi_psurf

subroutine write_omi_psurf(ifile, omi_psurf_temp, fform)

integer,                    intent(in) :: ifile
real(r8),		    intent(in) :: omi_psurf_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) omi_psurf_temp
   CASE DEFAULT
      write(ifile, *) omi_psurf_temp
END SELECT

end subroutine write_omi_psurf

function read_omi_ptrop(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_omi_ptrop
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_omi_ptrop
   CASE DEFAULT
      read(ifile, *) read_omi_ptrop
END SELECT

end function read_omi_ptrop

subroutine write_omi_ptrop(ifile, omi_ptrop_temp, fform)

integer,                    intent(in) :: ifile
real(r8),		    intent(in) :: omi_ptrop_temp
character(len=32),          intent(in) :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

if ( .not. module_initialized ) call initialize_module

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) omi_ptrop_temp
   CASE DEFAULT
      write(ifile, *) omi_ptrop_temp
END SELECT

end subroutine write_omi_ptrop

function read_omi_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
real(r8), dimension(35)        :: read_omi_avg_kernels
character(len=*), intent(in), optional :: fform

character(len=5)   :: header
character(len=129) :: errstring
character(len=32)  :: fileformat

read_omi_avg_kernels(:) = 0.0_r8

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_omi_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_omi_avg_kernels(1:nlevels)
END SELECT

end function read_omi_avg_kernels

subroutine write_omi_avg_kernels(ifile, avg_kernels_temp, nlevels_temp, fform)

integer,                    intent(in) :: ifile, nlevels_temp
real(r8), dimension(35), intent(in)  :: avg_kernels_temp
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

end subroutine write_omi_avg_kernels

end module obs_def_omi_mod
! END DART PREPROCESS MODULE CODE
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$