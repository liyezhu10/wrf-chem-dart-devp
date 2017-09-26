! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

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

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                             nmlfileunit, check_namelist_read, &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use     location_mod, only : location_type, set_location, get_location, VERTISPRESSURE, VERTISSURFACE, &
        VERTISUNDEF

use  assim_model_mod, only : interpolate
use    obs_kind_mod, only  : KIND_CO, KIND_SURFACE_PRESSURE, KIND_LANDMASK

implicit none
private

public :: write_mopitt_co, &
          read_mopitt_co, &
          interactive_mopitt_co, &
          get_expected_mopitt_co, &
          set_obs_def_mopitt_co

! Storage for the special information required for observations of this type
integer, parameter               :: MAX_MOPITT_CO_OBS = 10000000
integer, parameter               :: MOPITT_DIM = 10
integer                          :: num_mopitt_co_obs = 0
real(r8)   :: mopitt_pressure(MOPITT_DIM) =(/ &
                              100000.,90000.,80000.,70000.,60000.,50000.,40000.,30000.,20000.,1000. /)
real(r8)   :: mopitt_pressure_mid(MOPITT_DIM) =(/ &
                              100000.,85000.,75000.,65000.,55000.,45000.,35000.,25000.,15000.,7500. /)

real(r8), allocatable, dimension(:,:) :: avg_kernel
real(r8), allocatable, dimension(:) :: mopitt_prior
real(r8), allocatable, dimension(:) :: mopitt_psurf
integer,  allocatable, dimension(:) :: mopitt_nlevels

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2

logical, save :: module_initialized = .false.
integer  :: counts1 = 0

character(len=129)  :: MOPITT_CO_retrieval_type
logical             :: use_log_co
!
! MOPITT_CO_retrieval_type:
!     RAWR - retrievals in VMR (ppb) units
!     RETR - retrievals in log10(VMR ([ ])) units
!     QOR  - quasi-optimal retrievals
!     CPSR - compact phase space retrievals
    namelist /obs_def_MOPITT_CO_nml/ MOPITT_CO_retrieval_type, use_log_co

contains

!----------------------------------------------------------------------

subroutine initialize_module

integer :: iunit, rc

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

allocate (avg_kernel(    MAX_MOPITT_CO_OBS,MOPITT_DIM))
allocate (mopitt_prior(  MAX_MOPITT_CO_OBS))
allocate (mopitt_psurf(  MAX_MOPITT_CO_OBS))
allocate (mopitt_nlevels(MAX_MOPITT_CO_OBS))

! Read the namelist entry.
MOPITT_CO_retrieval_type='RETR'
use_log_co=.false.
call find_namelist_in_file("input.nml", "obs_def_MOPITT_CO_nml", iunit)
read(iunit, nml = obs_def_MOPITT_CO_nml, iostat = rc)
call check_namelist_read(iunit, rc, "obs_def_MOPITT_CO_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_MOPITT_CO_nml)
if (do_nml_term()) write(     *     , nml=obs_def_MOPITT_CO_nml)

end subroutine initialize_module

subroutine read_mopitt_co(key, ifile, fform)
!----------------------------------------------------------------------
!subroutine read_mopitt_co(key, ifile, fform)

integer,          intent(out)          :: key
integer,          intent(in)           :: ifile
character(len=*), intent(in), optional :: fform

character(len=32)                   :: fileformat

integer  :: mopitt_nlevels_1
real(r8) :: mopitt_prior_1
real(r8) :: mopitt_psurf_1
real(r8), dimension(MOPITT_DIM) :: avg_kernels_1
integer  :: keyin

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
character(len=*), intent(in), optional :: fform

character(len=32) :: fileformat
real(r8), dimension(MOPITT_DIM) :: avg_kernels_temp

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

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)
if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   ! PUT IN ERROR HANDLER CALL
   write(string1, *)'Not enough space for a mopitt CO obs.'
   write(string2, *)'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'interactive_mopitt_co',string1,source,revision,revdate, text2=string2)
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

    integer,parameter               :: nlvls=33
    real(r8),dimension(nlvls)       :: prs_prf_ocn,co_prf_ocn,prs_prf_lnd,co_prf_lnd
    real(r8), intent(in)            :: state(:)
    type(location_type), intent(in) :: location
    integer, intent(in)             :: key
    real(r8), intent(out)           :: val
    integer, intent(out)            :: istatus
!
    integer :: i,ii,kstr
    type(location_type) :: loc2
    real(r8)            :: mloc(3),fac,obs_val_min,obs_val_min_log
    real(r8)            :: obs_val,wrf_psf,level,wrf_landmask
    real(r8)            :: co_min,co_min_log,mopitt_prs_mid,mopitt_psf
    real(r8)            :: wt_lnd_up,wt_lnd_dw
    real(r8)            :: wt_ocn_up,wt_ocn_dw
    real(r8)            :: val_lnd,val_ocn
!
    integer             :: nlevels,nnlevels,check_min
    integer             :: ilv_ocn,ilv_lnd

!
! Initialize DART
    if ( .not. module_initialized ) call initialize_module
!
! Initialize variables
    co_min=1.e-2
    co_min_log=-8.
    obs_val_min=4.e-2
    obs_val_min_log=-8.
    if ( use_log_co ) then
       co_min=co_min_log
       obs_val_min=obs_val_min_log
    endif
!
! Set background profiles
    prs_prf_ocn(1:nlvls)=(/101282.5,100419.1,99251.30,97776.66, &
    95951.41,93671.55,90937.23,87145.57,82497.19,77853.89,73266.28,66969.46,59372.49,52488.63, &
    46259.94,40637.37,35567.76,31051.95,26937.60,23274.84,20014.34,17104.80,14545.98,12288.55, &
    10292.10,8526.316,6990.605,5656.344,4497.387,3498.531,2640.769,1910.597,1288.596/)
!
    co_prf_ocn(1:nlvls)=(/0.1053032,0.1052274,0.1050205, &
    0.1049154,0.1044448,0.1032934,0.1016622,9.9698491E-02,9.7816505E-02,9.6885450E-02, &
    9.7135067E-02,9.8520398E-02,0.1009575,0.1030745,0.1053934,0.1080828,0.1108766,0.1125486, &
    0.1116490,0.1075780,0.1003375,8.8890545E-02,7.5979680E-02,6.2236544E-02,5.1421385E-02, &
    4.3278556E-02,3.5766449E-02,2.5966400E-02,1.5804987E-02,1.7072648E-02,1.8331194E-02, &   
    2.1269552E-02,2.4853116E-02/)
!
    prs_prf_lnd(1:nlvls)=(/92520.25,91734.84,90671.42,89330.86, &
    87666.75,85588.16,83095.80,79636.74,75398.21,71163.90,66979.16,61235.73,54303.09,48017.58,  &
    42331.53,37197.62,32569.12,28444.63,24687.65,21343.15,18365.25,15708.11,13371.58,11309.93,  &
    9486.382,7873.569,6471.439,5252.575,4194.045,3282.142,2498.509,1831.732,1263.482/)
!
    co_prf_lnd(1:nlvls)=(/0.1451281,0.1406095,0.1350518,0.1303771, &
    0.1271011,0.1245595,0.1224129,0.1199306,0.1172802,0.1150899,0.1134757,0.1111645,0.1083032,  &
    0.1057542,0.1039881,0.1032084,0.1028377,0.1018478,9.8806620E-02,9.2105977E-02,8.2429104E-02, &
    7.1644865E-02,6.0651097E-02,5.1242340E-02,4.3855492E-02,3.6937837E-02,2.9593525E-02,2.1433521E-02, &
    1.6749103E-02,1.8461850E-02,1.9997887E-02,2.1560470E-02,2.3468491E-02/)
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
! Get wrf surface pressure
    wrf_psf = 0.0_r8
    istatus = 0
    loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
    call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)  
!
! Correct mopitt surface pressure
    if(istatus/=0) then
!       write(string1, *)'APM NOTICE: MOPITT CO WRF psf is bad ',wrf_psf,istatus
!       call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
       obs_val=missing_r8
       return
    endif              
!
    if(mopitt_psf.gt.wrf_psf) then
       if((mopitt_psf-wrf_psf).gt.10000.) then
          write(string1, *)'APM: NOTICE - reject MOPITT CO - WRF PSF too large ', &
          mopitt_psf,wrf_psf
          call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
          obs_val=missing_r8
          istatus=2
          return
       else
!          write(string1, *)'APM: NOTICE correct MOPITT CO psf with WRF psf ',mopitt_psf,wrf_psf
!          call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
          mopitt_psf=wrf_psf
       endif
    endif
!
! Find kstr - the surface level index
    kstr=0
    do i=1,MOPITT_DIM
       if (i.eq.1 .and. mopitt_psf.gt.mopitt_pressure(2)) then
          kstr=i
          exit
       endif
       if (i.ne.1 .and. i.ne.MOPITT_DIM .and. mopitt_pressure(i).ge.mopitt_psf .and. &
       mopitt_psf.gt.mopitt_pressure(i+1)) then
          kstr=i
          exit   
       endif
    enddo
    if (kstr.eq.0) then
       write(string1, *)'APM: ERROR in MOPITT CO obs def kstr=0: mopitt_psf= ',mopitt_psf
       call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
       call exit_all()
    elseif (kstr.gt.6) then
       write(string1, *)'APM: ERROR MOPITT CO surface pressure is unrealistic: mopitt_psf, wrf_psf= ', &
       mopitt_psf,wrf_psf
       call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
       call exit_all()
    endif
!
! Reject ob when number of MOPITT levels from WRF cannot equal actual number of MOPITT levels
    nnlevels=MOPITT_DIM-kstr+1
    if(nnlevels.ne.nlevels) then
       write(string1, *)'APM: NOTICE reject MOPITT CO ob - WRF MOP  levels .ne. MOP levels, nnlvls,nlvls ', &
       nnlevels,nlevels
       call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
       obs_val=missing_r8
       istatus=2
       return
    endif   
!
! Get landmask
    istatus = 0
    wrf_landmask=-9999
    loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISUNDEF)
    call interpolate(state, loc2, KIND_LANDMASK, wrf_landmask, istatus)  
    if (istatus /= 0) then
       write(string1, *)'APM ERROR: WRF interpolate wrf landmask error ', &
       wrf_landmask,istatus
       call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
       wrf_landmask = missing_r8 
       obs_val = missing_r8
       istatus = 2
       call exit_all()
    endif
!
! Apply MOPITT Averaging kernel A and MOPITT Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 10 element vector 
    val = 0.0_r8
    do i=1,nlevels
!       print *, ' '
!       if(i.eq.1) print *, 'APM: Printing START '
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
          write(string1, *)'APM NOTICE: WRF extrapolation needed reject MOPITT CO ob '
          call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
          obs_val = missing_r8 
          return
       endif
!       print *, 'level, mdl_val,prs ',i,obs_val,mopitt_prs_mid
!       print *, 'kstr,psf,nlevels ',kstr,mopitt_psf,nlevels
!       print *, 'loc_1, loc_2, loc_3 ',mloc(1),mloc(2),mloc(3)
!       print *, 'avg_k ',i,avg_kernel(key,:)
!
! Check for min values
       check_min=0 
       if (check_min.eq.0) then
          if (obs_val.lt.co_min) then
             write(string1, *)'APM NOTICE: resetting minimum MOPITT CO value ', &
             obs_val,co_min,obs_val_min
             call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source, &
             revision,revdate)
             obs_val=obs_val_min
          endif
       else
!
! Check for min values
          ilv_ocn=-9999
          do ii=1,nlvls-1
             if(mopitt_prs_mid.ge.prs_prf_ocn(1)) then
                 ilv_ocn=1
                 exit
             elseif(mopitt_prs_mid.lt.prs_prf_ocn(nlvls)) then
                 ilv_ocn=nlvls
                 exit
             elseif(mopitt_prs_mid.lt.prs_prf_ocn(ii) .and. &
             mopitt_prs_mid.ge.prs_prf_ocn(ii+1)) then
                 ilv_ocn=ii
                 exit
             endif
          enddo     
          ilv_lnd=-9999
          do ii=1,nlvls-1
             if(mopitt_prs_mid.ge.prs_prf_lnd(1)) then
                 ilv_lnd=1
                 exit
             elseif(mopitt_prs_mid.lt.prs_prf_lnd(nlvls)) then
                 ilv_lnd=nlvls
                 exit
             elseif(mopitt_prs_mid.lt.prs_prf_lnd(ii) .and. &
             mopitt_prs_mid.ge.prs_prf_lnd(ii+1)) then
                 ilv_lnd=ii
                 exit
             endif
          enddo
          wt_ocn_up=-9999.
          wt_ocn_dw=-9999.
          wt_lnd_up=-9999.
          wt_lnd_dw=-9999.
          if(ilv_ocn.eq.1) then
             val_ocn=co_prf_ocn(1)
          elseif(ilv_ocn.eq.nlvls) then
             val_ocn=co_prf_ocn(nlvls)
          else       
             wt_ocn_up=(mopitt_prs_mid-prs_prf_ocn(ilv_ocn+1))
             wt_ocn_dw=(prs_prf_ocn(ilv_ocn)-mopitt_prs_mid)
             val_ocn=(co_prf_ocn(ilv_ocn)*wt_ocn_up + &
             co_prf_ocn(ilv_ocn+1)*wt_ocn_dw)/(wt_ocn_up+wt_ocn_dw)
          endif
          if(ilv_lnd.eq.1) then
             val_lnd=co_prf_lnd(1)
          elseif(ilv_lnd.eq.nlvls) then
             val_lnd=co_prf_lnd(nlvls)
          else       
             wt_lnd_up=(mopitt_prs_mid-prs_prf_lnd(ilv_lnd+1))
             wt_lnd_dw=(prs_prf_lnd(ilv_lnd)-mopitt_prs_mid)
             val_lnd=(co_prf_lnd(ilv_lnd)*wt_lnd_up + &
             co_prf_lnd(ilv_lnd+1)*wt_lnd_dw)/(wt_lnd_up+wt_lnd_dw)
          endif
          obs_val_min=(wrf_landmask-1.)*val_ocn + (2.-wrf_landmask)*val_lnd
!          print *,'ilv_ocn, ilv_lnd ',ilv_ocn,ilv_lnd
!          print *,'wt_ocn_up, wt_ocn_dw ',wt_ocn_up,wt_ocn_dw
!          print *,'wt_lnd_up, wt_lnd_dw ',wt_lnd_up,wt_lnd_dw
!          print *,'val_ocn, val_lnd ',val_ocn,val_lnd
!          print *,'landmask ',wrf_landmask
!          print *,'obs_val_min,obs_val ',obs_val_min,obs_val
          if(obs_val_min.ge..100) then
             fac=.5
          elseif(obs_val_min.lt..100 .and. obs_val_min.ge..040) then
             fac=.75
          elseif(obs_val_min.lt..040) then
             fac=1.
          endif
          if (obs_val.lt.fac*obs_val_min) then
             write(string1, *)'APM NOTICE: resetting minimum MOPITT CO value ', &
             obs_val,fac*obs_val_min
             call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
             obs_val=fac*obs_val_min
          endif
       endif
!
! apply averaging kernel
!
! Use this form for RAWR, RETR, QOR, and CPSR
       if(use_log_co) then
          val = val + avg_kernel(key,i) * obs_val  
       else
          val = val + avg_kernel(key,i) * log10(obs_val*1.e-6)  
       endif
!       print *, 'val_itr ',i,val
!       print *, 'avg_ker, obs_val, val ',avg_kernel(key,i),obs_val, &
!       log10(obs_val*1.e-6)
!       print *, 'level, mdl_val,prs ',i,obs_val,mopitt_prs_mid
!       print *, 'kstr,psf,nlevels ',kstr,mopitt_psf,nlevels
!       print *, 'loc_1, loc_2, loc_3 ',mloc(1),mloc(2),mloc(3)
!       print *, 'avg_k ',i,avg_kernel(key,:)
    enddo
    val = val + mopitt_prior(key)
! 
! Use this correction for raw retrievals
    if(trim(MOPITT_CO_retrieval_type) .eq. 'RAWR') then
       val = (10.**val)*1.e6
    endif
!
end subroutine get_expected_mopitt_co
!
!----------------------------------------------------------------------

 subroutine set_obs_def_mopitt_co(key, co_avgker, co_prior, co_psurf, co_nlevels)
!----------------------------------------------------------------------
! Allows passing of obs_def special information 

integer,                 intent(in) :: key, co_nlevels
real(r8), dimension(10), intent(in) :: co_avgker
real(r8),                intent(in) :: co_prior
real(r8),                intent(in) :: co_psurf

if ( .not. module_initialized ) call initialize_module

if(num_mopitt_co_obs >= MAX_MOPITT_CO_OBS) then
   
   write(string1, *)'Not enough space for a mopitt CO obs.'
   call error_handler(E_MSG,'set_obs_def_mopitt_co',string1,source,revision,revdate)
   write(string1, *)'Can only have MAX_MOPITT_CO_OBS (currently ',MAX_MOPITT_CO_OBS,')'
   call error_handler(E_ERR,'set_obs_def_mopitt_co',string1,source,revision,revdate)
endif

avg_kernel(key,:)   = co_avgker(:)
mopitt_prior(key)   = co_prior
mopitt_psurf(key)   = co_psurf
mopitt_nlevels(key) = co_nlevels

end subroutine set_obs_def_mopitt_co


function read_mopitt_prior(ifile, fform)

integer,                    intent(in) :: ifile
real(r8)                               :: read_mopitt_prior
character(len=*), intent(in), optional :: fform

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

integer,           intent(in) :: ifile
real(r8),          intent(in) :: mopitt_prior_temp
character(len=32), intent(in) :: fform

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

integer,           intent(in) :: ifile
real(r8),          intent(in) :: mopitt_psurf_temp
character(len=32), intent(in) :: fform

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

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
