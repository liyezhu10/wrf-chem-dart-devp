! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! BEGIN DART PREPROCESS KIND LIST
! IASI_O3_RETRIEVAL, KIND_O3
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!   use obs_def_iasi_O3_mod, only : write_iasi_o3, read_iasi_o3, &
!                                  interactive_iasi_o3, get_expected_iasi_o3, &
!                                  set_obs_def_iasi_o3
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(IASI_O3_RETRIEVAL)
!            call get_expected_iasi_o3(state, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call read_iasi_o3(obs_def%key, ifile, fileformat)
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call write_iasi_o3(obs_def%key, ifile, fileformat)
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!      case(IASI_O3_RETRIEVAL)
!         call interactive_iasi_o3(obs_def%key)
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS SET_OBS_DEF_IASI_O3
!      case(IASI_O3_RETRIEVAL)
!         call set_obs_def_iasi_o3(obs_def%key)
! END DART PREPROCESS SET_OBS_DEF_IASI_O3

! BEGIN DART PREPROCESS MODULE CODE

module obs_def_iasi_O3_mod

use         types_mod, only : r4, r8
use     utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use      location_mod, only : location_type, set_location, get_location, VERTISHEIGHT,&
                              VERTISPRESSURE, VERTISLEVEL, VERTISSURFACE
use   assim_model_mod, only : interpolate
use      obs_kind_mod, only : KIND_O3, KIND_SURFACE_PRESSURE, KIND_PRESSURE
use mpi_utilities_mod, only : my_task_id

implicit none
private

public :: write_iasi_o3, &
          read_iasi_o3, &
          interactive_iasi_o3, &
          get_expected_iasi_o3, &
          set_obs_def_iasi_o3

! Storage for the special information required for observations of this type
integer, parameter          :: MAX_IASI_O3_OBS = 6000000
integer, parameter          :: IASI_DIM = 40
integer                     :: num_iasi_o3_obs = 0
integer                     :: counts1 = 0

real(r8) :: iasi_o3_prior(  MAX_IASI_O3_OBS)          ! prior term of x=Ax + (I-A)xa + Gey
integer  :: iasi_nlevels(   MAX_IASI_O3_OBS)          ! number of iasi levels used
real(r8) :: avg_kernel(     MAX_IASI_O3_OBS,IASI_DIM) ! iasii averaging kernel
real(r8) :: iasi_air_column(MAX_IASI_O3_OBS,IASI_DIM) ! iasi air column profile
real(r8) :: iasi_heights(   MAX_IASI_O3_OBS,IASI_DIM) ! iasi retrieval heights
real(r8) :: iasi_pressure(  MAX_IASI_O3_OBS,IASI_DIM)

! nominal iasi height levels in m
real(r8)                    :: iasi_altitude(IASI_DIM) =(/ &
                               500.,1500.,2500.,3500.,4500., &
                               5500.,6500.,7500.,8500.,9500., &
                               10500.,11500.,12500.,13500.,14500., &
                               15500.,16500.,17500.,18500.,19500., &
                               20500.,21500.,22500.,23500.,24500., &
                               25500.,26500.,27500.,28500.,29500., &
                               30500.,31500.,32500.,33500.,34500., &
                               35500.,36500.,37500.,38500.,39500. /)

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3
logical, save :: module_initialized = .false.

contains

!-----------------------------------------------------------------------
!>

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!-----------------------------------------------------------------------
!>

subroutine read_iasi_o3(key, ifile, fform)
integer,                    intent(out) :: key
integer,                    intent(in)  :: ifile
character(len=*), optional, intent(in)  :: fform

! temp variables

character(len=32) :: fileformat

integer  :: nlevel_1
real(r8) :: prior_1
real(r8) :: altitude_1(  IASI_DIM)
real(r8) :: pressure_1(  IASI_DIM)
real(r8) :: avg_kernel_1(IASI_DIM)
real(r8) :: aircol_1(    IASI_DIM)
integer  :: keyin

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

nlevel_1     = read_iasi_num_levels(  ifile, fileformat)
prior_1      = read_iasi_prior_column(ifile, fileformat)
altitude_1   = read_iasi_heights(     ifile, nlevel_1, fileformat)
pressure_1   = read_iasi_pressure(    ifile, nlevel_1, fileformat)
avg_kernel_1 = read_iasi_avg_kernels( ifile, nlevel_1, fileformat)
aircol_1     = read_iasi_air_column(  ifile, nlevel_1, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) keyin

   CASE DEFAULT
      read(ifile, *) keyin
END SELECT

counts1 = counts1 + 1
key     = counts1

call set_obs_def_iasi_o3(key,prior_1,altitude_1(1:nlevel_1),pressure_1(1:nlevel_1), &
              avg_kernel_1(1:nlevel_1),aircol_1(1:nlevel_1),nlevel_1)

end subroutine read_iasi_o3

!-----------------------------------------------------------------------
!>

subroutine write_iasi_o3(key, ifile, fform)

integer,                    intent(in) :: key
integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform

! temp variables

character(len=32) :: fileformat
integer  :: nlevel_1
real(r8) :: iasi_o3_prior_1
real(r8) :: altitude_1(       IASI_DIM)
real(r8) :: pressure_1(       IASI_DIM)
real(r8) :: avg_kernel_1(     IASI_DIM)
real(r8) :: iasi_air_column_1(IASI_DIM)

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

nlevel_1                      = iasi_nlevels(   key)
iasi_o3_prior_1               = iasi_o3_prior(  key)
altitude_1(       1:nlevel_1) = iasi_heights(   key,:)
pressure_1(       1:nlevel_1) = iasi_pressure(  key,:)
avg_kernel_1(     1:nlevel_1) = avg_kernel(     key,:)
iasi_air_column_1(1:nlevel_1) = iasi_air_column(key,:)

call write_iasi_num_levels(  ifile, nlevel_1, fileformat)
call write_iasi_prior_column(ifile, iasi_o3_prior_1, fileformat)
call write_iasi_heights(     ifile, altitude_1(       1:nlevel_1), nlevel_1, fileformat)
call write_iasi_pressure(    ifile, pressure_1(       1:nlevel_1), nlevel_1, fileformat)
call write_iasi_avg_kernels( ifile, avg_kernel_1(     1:nlevel_1), nlevel_1, fileformat)
call write_iasi_air_column(  ifile, iasi_air_column_1(1:nlevel_1), nlevel_1, fileformat)

SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) key
   CASE DEFAULT
      write(ifile, *) key
END SELECT

end subroutine write_iasi_o3

!-----------------------------------------------------------------------
!> Initializes the specialized part of a IASI observation

subroutine interactive_iasi_o3(key)

! Passes back up the key for this one

integer, intent(out) :: key

if ( .not. module_initialized ) call initialize_module

! Make sure there's enough space, if not die for now (clean later)

if(num_iasi_o3_obs >= MAX_IASI_O3_OBS) then
   write(string1, *)'Not enough space for a iasi O3 obs.'
   write(string2, *)'Can only have MAX_IASI_O3_OBS (currently ',MAX_IASI_O3_OBS,')'
   call error_handler(E_ERR, 'interactive_iasi_o3', string1, &
              source, revision, revdate, text2=string2)
endif

! Increment the index

num_iasi_o3_obs = num_iasi_o3_obs + 1
key = num_iasi_o3_obs

! Otherwise, prompt for input for the three required beasts

write(*, *) 'Creating an interactive_iasi_o3 observation'

end subroutine interactive_iasi_o3

!-----------------------------------------------------------------------
!>

subroutine get_expected_iasi_o3(state, location, key, val, istatus)

real(r8),            intent(in)  :: state(:)
type(location_type), intent(in)  :: location
integer,             intent(in)  :: key
real(r8),            intent(out) :: val
integer,             intent(out) :: istatus

integer,parameter   :: nlev_wrf=32
integer             :: i, ilev, kstr, apm_dom, apm_mm, nnlevels
type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)            :: wrf_prs(nlev_wrf)
real(r8)            :: obs_val,ubv_obs_val,level,missing
real(r8)            :: iasi_psf,iasi_psf_sv,wrf_psf,ubv_delt_prs
real(r8)            :: iasi_prs,iasi_o3_min,ylat,ylon,ylev
integer             :: nlevels,iubv,icnt
character(len=20)   :: apm_spec
real(r8)            :: vert_mode_filt

if ( .not. module_initialized ) call initialize_module

mloc = get_location(location)

! Initialize variables
iasi_o3_min = 1.e-6
missing     = -9999.9_r8
icnt=0

! Get IASI data
nlevels     = iasi_nlevels(key)
iasi_psf    = iasi_pressure(key,1)
iasi_psf_sv = iasi_pressure(key,1)

! Get location information
if (mloc(2) .gt. 90.0_r8) then
   mloc(2)=90.0_r8
elseif (mloc(2) .lt. -90.0_r8) then
   mloc(2)=-90.0_r8
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Try phase space localization
vert_mode_filt=10000.
if(mloc(3).le.vert_mode_filt) then
   istatus=2
   obs_val=missing
   return
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get wrf surface pressure
istatus = 0
wrf_psf = 0.0_r8
loc2 = set_location(mloc(1), mloc(2), 0.0_r8, VERTISSURFACE)
call interpolate(state, loc2, KIND_SURFACE_PRESSURE, wrf_psf, istatus)

! Correct iasi surface pressure
if(istatus/=0) then
   write(string1, *)'APM NOTICE: IASI O3 WRF psf is bad ',wrf_psf,istatus
   call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
   obs_val=missing
   return
endif

if(iasi_psf.gt.wrf_psf) then
   if((iasi_psf-wrf_psf).gt.10000.) then
      write(string1, *)'APM: NOTICE - reject IASI O3 - WRF PSF too large ',iasi_psf,wrf_psf
      call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
      istatus=2
      obs_val=missing
      return
   endif
else
!         write(string1, *)'APM: NOTICE correct IASI O3 psf with WRF psf ',iasi_psf,wrf_psf
!         call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
   iasi_psf=wrf_psf
endif

! Find kstr - the surface level index
kstr=0
do i=1,nlevels
   if (i.eq.1 .and. iasi_psf.gt.iasi_pressure(key,2)) then
      kstr=i
      exit
   endif
   if (i.ne.1 .and. i.ne.nlevels .and. iasi_pressure(key,i).ge.iasi_psf .and. &
   iasi_psf.gt.iasi_pressure(key,i+1)) then
      kstr=i
      exit
   endif
enddo
if (kstr.eq.0) then
   write(string1, *)'APM: ERROR in IASI O3 obs def kstr=0: iasi_psf= ',iasi_psf
   call error_handler(E_ERR,'set_obs_def_iasi_o3',string1,source,revision,revdate)
elseif (kstr.gt.6) then
   write(string1, *)'APM: ERROR IASI O3 psf is unrealistic: iasi_psf, wrf_psf= ',iasi_psf,wrf_psf
   call error_handler(E_MSG,'set_obs_def_iasi_03',string1,source,revision,revdate)
endif

! Reject ob when number of IASI levels from WRF cannot equal actual number of IASI levels
nnlevels=nlevels-kstr+1
if(nnlevels.ne.nlevels) then
   write(string1, *)'APM: NOTICE reject IASI O3 ob - WRF IASI levels .ne. IASI levels, nnlvls,nlvls ',nnlevels,nlevels
   call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
   istatus=2
   obs_val=missing
   return
endif

! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 40-element vector !

! Get wrf pressue profile
wrf_prs(1:nlev_wrf)=0.0_r8
do ilev = 1, nlev_wrf
   level = real(ilev)
   loc2 = set_location(mloc(1),mloc(2), level, VERTISLEVEL)
   call interpolate(state, loc2, KIND_PRESSURE, wrf_prs(ilev), istatus)
enddo

! loop through IASI levels
val = 0.0_r8
do ilev = 1, nlevels

! get location of obs
   if (ilev .eq. 1) then
      iasi_prs = (iasi_psf+iasi_pressure(key,ilev))/2.
      loc2 = set_location(mloc(1),mloc(2), iasi_prs, VERTISPRESSURE)
   else
      iasi_prs = (iasi_pressure(key,ilev-1)+iasi_pressure(key,ilev))/2.
      loc2 = set_location(mloc(1),mloc(2), iasi_prs, VERTISPRESSURE)
   endif
   istatus = 0
   obs_val = 0.0_r8

! interpolate to obs location
   call interpolate(state, loc2, KIND_O3, obs_val, istatus)

! fix for failed near surface extrapolation
   if (istatus .eq. 2 .and. ilev .eq. 1) then
      loc2 = set_location(mloc(1),mloc(2), wrf_prs(1), VERTISPRESSURE)
      call interpolate(state, loc2, KIND_O3, obs_val, istatus)
!            write(string1, *)'APM: NOTICE IASI O3 sfc extrap fix ',ilev,istatus,obs_val
!            call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
   endif

! check for problems with the interpolation
   if (istatus .eq. 2 .and. iasi_prs .le. wrf_prs(nlev_wrf)) then
      istatus = 0
      apm_dom=01
      apm_spec='OX'
      apm_mm=6
      ylon=mloc(1)
      if(mloc(1) >= 180.) ylon=ylon-360.
      ylat=mloc(2)
      call wrf_dart_ubval_interp(ubv_obs_val,ubv_delt_prs, &
      apm_dom,apm_spec,ylon,ylat,iasi_prs,apm_mm,istatus)

! convert uvb to ppm
      obs_val=ubv_obs_val*1.e6
   endif

! interpolation failed
   if (istatus /= 0) then
      write(string1, *)'APM: NOTICE reject IASI O3 ob - WRF interpolation failed ',ilev,istatus,obs_val
      call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
      obs_val = missing
      return
   endif

! Check for WRF O3 lower bound
   if (obs_val.lt.iasi_o3_min) then
      write(string1, *)'APM: NOTICE resetting minimum IASI O3 value  '
      call error_handler(E_MSG,'set_obs_def_iasi_o3',string1,source,revision,revdate)
      obs_val = iasi_o3_min
   endif
   obs_val = obs_val * 1000.0_r8

! apply averaging kernel
   val = val + avg_kernel(key,ilev) * obs_val
enddo
!   val = val + iasi_prior(key)

end subroutine get_expected_iasi_o3

!-----------------------------------------------------------------------
!>

subroutine wrf_dart_ubval_interp(obs_val,del_prs,domain,species,lon,lat,lev,im2,istatus)

use netcdf
implicit none

integer,parameter                                 :: nx1=2,ny1=96,nz1=38,nm1=12,nmspc1=8,nchr1=20
integer,parameter                                 :: nx2=100,ny2=40,nz2=66,nm2=12
integer                                           :: fid1,fid2,domain,rc,im2
integer                                           :: i,j,k,imn1,istatus
integer,dimension(nm1)                            :: month
character(len=20)                                 :: species
character(len=180)                                :: file_nam,file_in1,file_in2,path
character(len=nchr1),dimension(nmspc1)            :: spec_nam
real(r8)                                          :: lon,lat,lev,del_lon1
real(r8)                                          :: obs_val,del_prs
real(r4),dimension(nx1,ny1)                       :: xlon1,xlat1
real(r4),dimension(ny1)                           :: xlat_tmp1
real(r4),dimension(nz1)                           :: xlev1,prs_tmp1
real(r4),dimension(ny1,nz1)                       :: fld1,fld_tmp1
real(r4),dimension(nx1,ny1,nz1)                   :: fldd1
real(r4),dimension(nx2,ny2)                       :: xlat2,xlon2
real(r4),dimension(nz2)                           :: xlev2,prs_tmp2
real(r4),dimension(nm2)                           :: ddoyr
real(r4),dimension(nx2,ny2,nz2,nm2)               :: o3_col_dens
real(r4),dimension(nx2,ny2,nz2)                   :: fld_tmp2
logical                                           :: use_interp_1

! decide with upper boundary data to use
use_interp_1=.true.
del_lon1=2.

! assign upper boundary profile files
path='./'

! this file has OX NOX HNO3 CH4 CO N2X N2O5 H20
file_in1='ubvals_b40.20th.track1_1996-2005.nc'

! this file has o3 only
file_in2='exo_coldens_d01'
if(domain.eq.2) then
   file_in2='exo_coldens_d02'
endif

! open upper boundary profile files
if (use_interp_1) then
   file_nam=trim(path)//trim(file_in1)
   rc = nf90_open(trim(file_nam),NF90_NOWRITE,fid1)
   if(rc.ne.0) then
      print *, 'APM: nc_open error file=',trim(file_nam)
      call abort
   endif
!      print *, 'opened ',trim(file_nam)
else
   file_nam=trim(path)//trim(file_in2)
   rc = nf90_open(trim(file_nam),NF90_NOWRITE,fid2)
   if(rc.ne.0) then
      print *, 'APM: nc_open error file=',trim(file_nam)
      call abort
   endif
!      print *, 'opened ',trim(file_nam)
endif

! select upper boundary data from ubvals_b40.20th.track1_1996-2005.nc
if (use_interp_1) then
   imn1=6
   call apm_get_ubvals(fid1,species,imn1,fld1,xlat_tmp1,xlev1)
   rc=nf90_close(fid1)
else

! select upper boundary data from exo_coldens_dxx
   call apm_get_exo_coldens(fid2,'XLAT',xlat2,nx2,ny2,1,1)
!      print *, 'XLAT',xlat2(1,1),xlat2(nx2,ny2)
   call apm_get_exo_coldens(fid2,'XLONG',xlon2,nx2,ny2,1,1)
   !      print *, 'XLON',xlon2(1,1),xlon2(nx2,ny2)
   call apm_get_exo_coldens(fid2,'coldens_levs',xlev2,nz2,1,1,1)
!      print *, 'coldens_levs',xlev2(:)
   call apm_get_exo_coldens(fid2,'days_of_year',ddoyr,nm2,1,1,1)
!      print *, 'ddoyr',ddoyr(1),ddoyr(nm2)
   call apm_get_exo_coldens(fid2,'o3_column_density',o3_col_dens,nx2,ny2,nz2,nm2)
!      print *, 'o3_coldens',o3_col_dens(1,1,1,1),o3_col_dens(nx2,ny2,nz2,nm2)
   rc=nf90_close(fid2)
endif
!   print *, 'ny1,nz1 ',ny1,nz1
!   print *, 'fld1 ',fld1
!   print *, 'xlat1 ',xlat1
!   print *, 'xlev1 ',xlev1

! convert longitude to 0 - 360
if (.not.  use_interp_1) then
   do i=1,nx2
      do j=1,ny2
         if(xlon2(i,j).lt.0.) then
            xlon2(i,j)=xlon2(i,j)+360.
         endif
      enddo
   enddo
endif

! invert the pressure grid and data
if (use_interp_1) then
   do k=1,nz1
      prs_tmp1(nz1-k+1)=xlev1(k)
      do j=1,ny1
         fld_tmp1(j,nz1-k+1)=fld1(j,k)
      enddo
   enddo
   xlev1(1:nz1)=prs_tmp1(1:nz1)*100.
   fldd1(1,1:ny1,1:nz1)=fld_tmp1(1:ny1,1:nz1)
   fldd1(2,1:ny1,1:nz1)=fld_tmp1(1:ny1,1:nz1)

! interpolate data1 to (lat,lev) point
   do j=1,ny1
      xlon1(1,j)=lon-del_lon1
      xlon1(2,j)=lon+del_lon1
      if(lon.lt.0.) then
         xlon1(1,j)=lon+360.-del_lon1
         xlon1(2,j)=lon+360.+del_lon1
      endif
      do i=1,nx1
         xlat1(i,j)=xlat_tmp1(j)
      enddo
   enddo
!      print *, 'IN UBVAL SUB: lon,lat,lev ',lon,lat,lev
!      print *, 'IN UBVAL SUB: xlon,xlat,xlev ',xlon1(1,48),xlat1(1,48)
!      do j=1,nz1
!        print *, 'IN UBVAL SUB: fldd1 ',j,xlev1(j),fldd1(1,48,j)
!      enddo
   call apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon1,xlat1,xlev1, &
   fldd1,nx1,ny1,nz1,istatus)
!      print *, 'IN UBVAL SUB: obs_val,del_prs ',obs_val,del_prs
else
   do k=1,nz2
      prs_tmp2(nz2-k+1)=xlev2(k)
      do i=1,nx2
         do j=1,ny2
            fld_tmp2(i,j,nz2-k+1)=o3_col_dens(i,j,k,im2)
         enddo
      enddo
   enddo
   xlev2(1:nz2)=prs_tmp2(1:nz2)
   o3_col_dens(1:nx2,1:ny2,1:nz2,im2)=fld_tmp2(1:nx2,1:ny2,1:nz2)

! interpolate data2 to (lat,lon,lev) point
   call apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon2,xlat2,xlev2, &
   o3_col_dens(1,1,1,im2),nx2,ny2,nz2,istatus)
endif

end subroutine wrf_dart_ubval_interp

!-----------------------------------------------------------------------
!>

subroutine apm_get_exo_coldens(fid,fldname,dataf,nx,ny,nz,nm)

use netcdf
implicit none

integer,parameter                      :: maxdim=4
integer                                :: nx,ny,nz,nm
integer                                :: i,rc,v_ndim,natts,domain,fid
integer                                :: v_id,typ
integer,dimension(maxdim)              :: v_dimid,v_dim,one
character(len=*)                       :: fldname
character(len=180)                     :: vnam
real,dimension(nx,ny,nz,nm)            :: dataf

! get variables identifiers
rc = nf90_inq_varid(fid,trim(fldname),v_id)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error'
   call abort
endif

! get dimension identifiers
v_dimid=0
rc = nf90_inquire_variable(fid,v_id,vnam,typ,v_ndim,v_dimid,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error'
   call abort
endif
if(maxdim.lt.v_ndim) then
   print *, 'ERROR: maxdim is too small ',maxdim,v_ndim
   call abort
endif

! get dimensions
v_dim(:)=1
do i=1,v_ndim
   rc = nf90_inquire_dimension(fid,v_dimid(i),len=v_dim(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error'
      call abort
   endif
enddo

! check dimensions
if(nx.ne.v_dim(1)) then
   print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
   call abort
else if(ny.ne.v_dim(2)) then
   print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
   call abort
else if(nz.ne.v_dim(3)) then
   print *, 'ERROR: nz dimension conflict ',nz,v_dim(3)
   call abort
else if(nm.ne.v_dim(4)) then
   print *, 'ERROR: nm dimension conflict ',nm,v_dim(4)
   call abort
endif

! get data
one(:)=1
rc = nf90_get_var(fid,v_id,dataf,one,v_dim)

end subroutine apm_get_exo_coldens

!-----------------------------------------------------------------------
!>

subroutine apm_get_ubvals(fid,species,imn,dataf,lats,levs)

use netcdf
implicit none

integer,parameter                                 :: maxdim=4
integer,parameter                                 :: ny1=96,nz1=38,nm1=12,nmspc1=8,nmchr1=20
integer                                           :: i,j,idx,fid,rc,typ,natts,imn
integer                                           :: vid1,vid2,vid3,vid4,vid5
integer                                           :: vndim1,vndim2,vndim3,vndim4,vndim5
integer,dimension(maxdim)                         :: vdimid1,vdimid2,vdimid3,vdimid4,vdimid5
integer,dimension(maxdim)                         :: one,vdim1,vdim2,vdim3,vdim4,vdim5
integer,dimension(nm1)                            :: mths
character(len=20)                                 :: species
character(len=nmchr1),dimension(nmchr1,nmspc1)    :: spcs
character(len=180)                                :: vnam1,vnam2,vnam3,vnam4,vnam5
real(r4),dimension(ny1)                           :: lats
real(r4),dimension(nz1)                           :: levs
real(r4),dimension(ny1,nmspc1,nm1,nz1)            :: vmrs
real(r4),dimension(ny1,nz1)                       :: dataf

! get variables identifiers
rc = nf90_inq_varid(fid,'lat',vid1)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error_1'
   call abort
endif
rc = nf90_inq_varid(fid,'lev',vid2)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error_2'
   call abort
endif
rc = nf90_inq_varid(fid,'month',vid3)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error_3'
   call abort
endif
rc = nf90_inq_varid(fid,'specname',vid4)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error_4'
   call abort
endif
rc = nf90_inq_varid(fid,'vmr',vid5)
if(rc.ne.0) then
   print *, 'APM: nf_inq_varid error_5'
   call abort
endif

! get dimension identifiers
vdimid1=0
rc = nf90_inquire_variable(fid,vid1,vnam1,typ,vndim1,vdimid1,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error_1'
   call abort
endif
vdimid2=0
rc = nf90_inquire_variable(fid,vid2,vnam2,typ,vndim2,vdimid2,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error_2'
   call abort
endif
vdimid3=0
rc = nf90_inquire_variable(fid,vid3,vnam3,typ,vndim3,vdimid3,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error_3'
   call abort
endif
vdimid4=0
rc = nf90_inquire_variable(fid,vid4,vnam4,typ,vndim4,vdimid4,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error_4'
   call abort
endif
vdimid5=0
rc = nf90_inquire_variable(fid,vid5,vnam5,typ,vndim5,vdimid5,natts)
if(rc.ne.0) then
   print *, 'APM: nc_inq_var error_5'
   call abort
endif

! test the number of dimensions
if(1.lt.vndim1) then
   print *, 'ERROR: maxdim is too small 1 ',1,vndim1
   call abort
endif
if(1.lt.vndim2) then
   print *, 'ERROR: maxdim is too small 2 ',1,vndim2
   call abort
endif
if(1.lt.vndim3) then
   print *, 'ERROR: maxdim is too small 3 ',1,vndim3
   call abort
endif
if(2.lt.vndim4) then
   print *, 'ERROR: maxdim is too small 4',1,vndim4
   call abort
endif
if(4.lt.vndim5) then
   print *, 'ERROR: maxdim is too small 5',1,vndim5
   call abort
endif

! get dimensions
vdim1(:)=1
do i=1,vndim1
   rc = nf90_inquire_dimension(fid,vdimid1(i),len=vdim1(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error_1'
      call abort
   endif
enddo
vdim2(:)=1
do i=1,vndim2
   rc = nf90_inquire_dimension(fid,vdimid2(i),len=vdim2(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error_2'
      call abort
   endif
enddo
vdim3(:)=1
do i=1,vndim3
   rc = nf90_inquire_dimension(fid,vdimid3(i),len=vdim3(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error_3'
      call abort
   endif
enddo
vdim4(:)=1
do i=1,vndim4
   rc = nf90_inquire_dimension(fid,vdimid4(i),len=vdim4(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error_4'
      call abort
   endif
enddo
vdim5(:)=1
do i=1,vndim5
   rc = nf90_inquire_dimension(fid,vdimid5(i),len=vdim5(i))
   if(rc.ne.0) then
      print *, 'APM: nf_inq_dimlen error_5'
      call abort
   endif
enddo

! check dimensions
if(ny1.ne.vdim1(1)) then
   print *, 'ERROR: ny1 dimension conflict 1 ',ny1,vdim1(1)
   call abort
else if(nz1.ne.vdim2(1)) then
   print *, 'ERROR: nz1 dimension conflict 2 ',nz1,vdim2(1)
   call abort
else if(nm1.ne.vdim3(1)) then
   print *, 'ERROR: nm1 dimension conflict 3 ',nm1,vdim3(1)
   call abort
endif
if(nmchr1.ne.vdim4(1)) then
   print *, 'ERROR: nmchr1 dimension conflict 4 ',nmchr1,vdim4(1)
   call abort
else if(nmspc1.ne.vdim4(2)) then
   print *, 'ERROR: nmspc1 dimension conflict 4 ',nmspc1,vdim4(2)
   call abort
endif
if(ny1.ne.vdim5(1)) then
   print *, 'ERROR: ny1 dimension conflict 5 ',ny1,vdim5(1)
   call abort
else if(nmspc1.ne.vdim5(2)) then
   print *, 'ERROR: nmspc1 dimension conflict 5 ',nmspc1,vdim5(2)
   call abort
else if(nm1.ne.vdim5(3)) then
   print *, 'ERROR: nm1 dimension conflict 5 ',nm1,vdim5(3)
   call abort
else if(nz1.ne.vdim5(4)) then
   print *, 'ERROR: nz1 dimension conflict 5 ',nz1,vdim5(4)
   call abort
endif

! get data
one(:)=1
rc = nf90_get_var(fid,vid1,lats,one,vdim1)
if(rc.ne.0) then
   print *, 'APM: get_var error_1'
   call abort
endif
!   print *, 'lats ',lats
one(:)=1
rc = nf90_get_var(fid,vid2,levs,one,vdim2)
if(rc.ne.0) then
   print *, 'APM: get_var error_2'
   call abort
endif
!   print *, 'levs ',levs
one(:)=1
rc = nf90_get_var(fid,vid3,mths,one,vdim3)
if(rc.ne.0) then
   print *, 'APM: get_var error_3'
   call abort
endif
!   print *, 'mths ',mths
one(:)=1
rc = nf90_get_var(fid,vid4,spcs,one,vdim4)
if(rc.ne.0) then
   print *, 'APM: get_var error_4'
   call abort
endif
!   print *, 'spcs ',spcs
one(:)=1
rc = nf90_get_var(fid,vid5,vmrs,one,vdim5)
if(rc.ne.0) then
   print *, 'APM: get_var error_5'
   call abort
endif
!   print *, 'vmrs ',vmrs

! locate requested field
  do i=1,nmspc1
  if(trim(species).eq.trim(spcs(i,1))) then
     idx=i
     exit
  endif
  enddo
do i=1,ny1
   do j=1,nz1
      dataf(i,j)=vmrs(i,idx,imn,j)
   enddo
enddo

end subroutine apm_get_ubvals

!-----------------------------------------------------------------------
!>

subroutine apm_interpolate(obs_val,del_prs,lon,lat,lev,xlon,xlat,xlev,dataf,nx,ny,nz,istatus)

! longitude and latitude must be in degrees
! pressure grid must be in hPa and go from bottom to top

implicit none

integer                                :: nx,ny,nz,nzm,istatus
integer                                :: i,j,k,im,ip,jm,jp,quad
integer                                :: k_lw,k_up,i_min,j_min
real(r8)                               :: obs_val,del_prs
real(r8)                               :: lon,lat,lev
real(r4)                               :: l_lon,l_lat,l_lev
real(r4)                               :: fld_lw,fld_up
real(r4)                               :: xlnp_lw,xlnp_up,xlnp_pt
real(r4)                               :: dz_lw,dz_up
real(r4)                               :: mop_x,mop_y
real(r4)                               :: re,pi,rad2deg
real(r4)                               :: rad,rad_crit,rad_min,mod_x,mod_y
real(r4)                               :: dx_dis,dy_dis
real(r4)                               :: w_q1,w_q2,w_q3,w_q4,wt
real(r4),dimension(nz)                 :: xlev
real(r4),dimension(nx,ny)              :: xlon,xlat
real(r4),dimension(nx,ny,nz)           :: dataf

! set constants
pi=4.*atan(1.)
rad2deg=360./(2.*pi)
re=6371000.
rad_crit=200000.
quad=0

! find the closest point
rad_min=1.e10
l_lon=lon
l_lat=lat
l_lev=lev
if(l_lon.lt.0.) l_lon=l_lon+360.
!   print *, 'lon,lat,lev ',l_lon,l_lat,l_lev

do i=1,nx
   do j=1,ny
      mod_x=(xlon(i,j))/rad2deg
      if(xlon(i,j).lt.0.) mod_x=(360.+xlon(i,j))/rad2deg
      mod_y=xlat(i,j)/rad2deg
      mop_x=l_lon/rad2deg
      mop_y=l_lat/rad2deg
      dx_dis=abs(mop_x-mod_x)*cos((mop_y+mod_y)/2.)*re
      dy_dis=abs(mop_y-mod_y)*re
      rad=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
      rad_min=min(rad_min,rad)
      if(rad.eq.rad_min) then
         i_min=i
         j_min=j
      endif
   enddo
enddo
if(rad_min.gt.rad_crit) then
   print *, 'APM: ERROR in intrp - min dist exceeds threshold ',rad_min, rad_crit
   print *, 'grid ',i_min,j_min,xlon(i_min,j_min),xlat(i_min,j_min)
   print *, 'point ',l_lon,l_lat
   istatus=2
   return
!      call abort
endif

! do interpolation
im=i_min-1
if(im.eq.0) im=1
ip=i_min+1
if(ip.eq.nx+1) ip=nx
jm=j_min-1
if(jm.eq.0) jm=1
jp=j_min+1
if(jp.eq.ny+1) jp=ny

! find quadrant and interpolation weights
quad=0
mod_x=xlon(i_min,j_min)
if(xlon(i_min,j_min).lt.0.) mod_x=xlon(i_min,j_min)+360.
mod_y=xlat(i_min,j_min)
if(mod_x.ge.l_lon.and.mod_y.ge.l_lat) quad=1
if(mod_x.le.l_lon.and.mod_y.ge.l_lat) quad=2
if(mod_x.le.l_lon.and.mod_y.le.l_lat) quad=3
if(mod_x.ge.l_lon.and.mod_y.le.l_lat) quad=4
if(quad.eq.0) then
   print *, 'APM: ERROR IN INTERPOLATE quad = 0 '
   call abort
endif

! Quad 1
if (quad.eq.1) then
   mod_x=xlon(i_min,j_min)
   if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
   w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(im,j_min)
   if(xlon(im,j_min).lt.0.) mod_x=360.+xlon(im,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(im,j_min))/rad2deg*re
   w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(im,jm)
   if(xlon(im,jm).lt.0.) mod_x=360.+xlon(im,jm)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jm))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(im,jm))/rad2deg*re
   w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,jm)
   if(xlon(i_min,jm).lt.0.) mod_x=360.+xlon(i_min,jm)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jm))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,jm))/rad2deg*re
   w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 2
else if (quad.eq.2) then
   mod_x=xlon(ip,j_min)
   if(xlon(ip,j_min).lt.0.) mod_x=360.+xlon(ip,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(ip,j_min))/rad2deg*re
   w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,j_min)
   if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
   w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,jm)
   if(xlon(i_min,jm).lt.0.) mod_x=360.+xlon(i_min,jm)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jm))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,jm))/rad2deg*re
   w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(ip,jm)
   if(xlon(ip,jm).lt.0.) mod_x=360.+xlon(ip,jm)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,jm))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(ip,jm))/rad2deg*re
   w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 3
else if (quad.eq.3) then
   mod_x=xlon(ip,jp)
   if(xlon(ip,jp).lt.0.) mod_x=360.+xlon(ip,jp)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,jp))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(ip,jp))/rad2deg*re
   w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,jp)
   if(xlon(i_min,jp).lt.0.) mod_x=360.+xlon(i_min,jp)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jp))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,jp))/rad2deg*re
   w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,j_min)
   if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
   w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(ip,j_min)
   if(xlon(ip,j_min).lt.0.) mod_x=360.+xlon(ip,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(ip,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(ip,j_min))/rad2deg*re
   w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
! Quad 4
else if (quad.eq.4) then
   mod_x=xlon(i_min,jp)
   if(xlon(i_min,jp).lt.0.) mod_x=360.+xlon(i_min,jp)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,jp))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,jp))/rad2deg*re
   w_q1=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(im,jp)
   if(xlon(im,jp).lt.0.) mod_x=360.+xlon(im,jp)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jp))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(im,jp))/rad2deg*re
   w_q2=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(im,jm)
   if(xlon(im,jm).lt.0.) mod_x=360.+xlon(im,jm)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(im,jm))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(im,jm))/rad2deg*re
   w_q3=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
   mod_x=xlon(i_min,j_min)
   if(xlon(i_min,j_min).lt.0.) mod_x=360.+xlon(i_min,j_min)
   dx_dis=abs(l_lon-mod_x)/rad2deg*cos((l_lat+xlat(i_min,j_min))/rad2deg/2.)*re
   dy_dis=abs(l_lat-xlat(i_min,j_min))/rad2deg*re
   w_q4=sqrt(dx_dis*dx_dis + dy_dis*dy_dis)
endif
if(l_lon.ne.xlon(i_min,j_min).or.l_lat.ne.xlat(i_min,j_min)) then
   wt=1./w_q1+1./w_q2+1./w_q3+1./w_q4
endif

! find vertical indexes
nzm=nz-1
k_lw=-1
k_up=-1
do k=1,nzm
   if(k.eq.1 .and. l_lev.gt.xlev(k)) then
      k_lw=k
      k_up=k
      exit
   endif
   if(l_lev.le.xlev(k) .and. l_lev.gt.xlev(k+1)) then
      k_lw=k
      k_up=k+1
      exit
   endif
   if(k.eq.nzm .and. l_lev.ge.xlev(k+1)) then
      k_lw=k+1
      k_up=k+1
      exit
   endif
enddo
if(k_lw.le.0 .or. k_up.le.0) then
   print *, 'APM: ERROR IN K_LW OR K_UP ',k_lw,k_up
   call abort
endif

! horizontal interpolation
fld_lw=0.
fld_up=0.
if(l_lon.eq.xlon(i_min,j_min).and.l_lat.eq.xlat(i_min,j_min)) then
   fld_lw=dataf(i_min,j_min,k_lw)
   fld_up=dataf(i_min,j_min,k_up)
else if(quad.eq.1) then
   fld_lw=(1./w_q1*dataf(i_min,j_min,k_lw)+1./w_q2*dataf(im,j_min,k_lw)+ &
   1./w_q3*dataf(im,jm,k_lw)+1./w_q4*dataf(i_min,jm,k_lw))/wt
   fld_up=(1./w_q1*dataf(i_min,j_min,k_up)+1./w_q2*dataf(im,j_min,k_up)+ &
   1./w_q3*dataf(im,jm,k_up)+1./w_q4*dataf(i_min,jm,k_up))/wt
else if(quad.eq.2) then
   fld_lw=(1./w_q1*dataf(ip,j_min,k_lw)+1./w_q2*dataf(i_min,j_min,k_lw)+ &
   1./w_q3*dataf(i_min,jm,k_lw)+1./w_q4*dataf(ip,jm,k_lw))/wt
   fld_up=(1./w_q1*dataf(ip,j_min,k_up)+1./w_q2*dataf(i_min,j_min,k_up)+ &
   1./w_q3*dataf(i_min,jm,k_up)+1./w_q4*dataf(ip,jm,k_up))/wt
else if(quad.eq.3) then
   fld_lw=(1./w_q1*dataf(ip,jp,k_lw)+1./w_q2*dataf(i_min,jp,k_lw)+ &
   1./w_q3*dataf(i_min,j_min,k_lw)+1./w_q4*dataf(ip,j_min,k_lw))/wt
   fld_up=(1./w_q1*dataf(ip,jp,k_up)+1./w_q2*dataf(i_min,jp,k_up)+ &
   1./w_q3*dataf(i_min,j_min,k_up)+1./w_q4*dataf(ip,j_min,k_up))/wt
else if(quad.eq.4) then
   fld_lw=(1./w_q1*dataf(i_min,jp,k_lw)+1./w_q2*dataf(im,jp,k_lw)+ &
   1./w_q3*dataf(im,j_min,k_lw)+1./w_q4*dataf(i_min,j_min,k_lw))/wt
   fld_up=(1./w_q1*dataf(i_min,jp,k_up)+1./w_q2*dataf(im,jp,k_up)+ &
   1./w_q3*dataf(im,j_min,k_up)+1./w_q4*dataf(i_min,j_min,k_up))/wt
endif
!   print *,'fld_lw ',fld_lw
!   print *,'fld_up ',fld_up

! vertical interpolation
!   print *,'p_lw,p_up,p ',xlev(k_lw),xlev(k_up),l_lev

xlnp_lw=log(xlev(k_lw))
xlnp_up=log(xlev(k_up))
xlnp_pt=log(l_lev)
dz_lw=xlnp_lw-xlnp_pt
dz_up=xlnp_pt-xlnp_up
if(dz_lw.eq.0.) then
   obs_val=fld_lw
else if(dz_up.eq.0.) then
   obs_val=fld_up
else if(dz_lw.ne.0. .and. dz_up.ne.0.) then
   obs_val=(1./dz_lw*fld_lw+1./dz_up*fld_up)/(1./dz_lw+1./dz_up)
endif
del_prs=xlev(k_lw)-xlev(k_up)

end subroutine apm_interpolate

!-----------------------------------------------------------------------
!> Allows passing of obs_def special information

subroutine set_obs_def_iasi_o3(key, apcol_val, altretlev, pressure, akcol, aircol_val, nlev_use)

integer,  intent(in) :: key
integer,  intent(in) :: nlev_use
real(r8), intent(in) :: apcol_val
real(r8), intent(in) :: altretlev( IASI_DIM)
real(r8), intent(in) :: akcol(     IASI_DIM)
real(r8), intent(in) :: aircol_val(IASI_DIM)
real(r8), intent(in) :: pressure(  IASI_DIM)

if ( .not. module_initialized ) call initialize_module

! Check for sufficient space
if(num_iasi_o3_obs >= MAX_IASI_O3_OBS) then
   write(string1, *)'Not enough space for a iasi O3 obs.'
   write(string2, *)'Can only have MAX_IASI_O3_OBS (currently ',MAX_IASI_O3_OBS,')'
   call error_handler(E_ERR, 'set_obs_def_iasi_o3', string1, &
              source, revision, revdate, text2=string2)
endif

iasi_nlevels(   key)             = nlev_use
iasi_o3_prior(  key)             = apcol_val
iasi_heights(   key, 1:nlev_use) = altretlev( 1:nlev_use)
avg_kernel(     key, 1:nlev_use) = akcol(     1:nlev_use)
iasi_air_column(key, 1:nlev_use) = aircol_val(1:nlev_use)
iasi_pressure(  key, 1:nlev_use) = pressure(  1:nlev_use)

end subroutine set_obs_def_iasi_o3

!=================================
! other functions and subroutines
!=================================

function read_iasi_prior_column(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_prior_column

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)

   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior_column
   CASE DEFAULT
      read(ifile, *) read_iasi_prior_column
END SELECT

end function read_iasi_prior_column

!-----------------------------------------------------------------------
!>

subroutine write_iasi_prior_column(ifile, iasi_prior_temp, fform)

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

end subroutine write_iasi_prior_column

!-----------------------------------------------------------------------
!>

function read_iasi_num_levels(ifile, fform)

integer,                    intent(in) :: ifile
character(len=*), optional, intent(in) :: fform
integer                                :: read_iasi_num_levels

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_num_levels
   CASE DEFAULT
      read(ifile, *) read_iasi_num_levels
END SELECT

end function read_iasi_num_levels

!-----------------------------------------------------------------------
!>

subroutine write_iasi_num_levels(ifile, number_of_levels_temp, fform)

integer,          intent(in) :: ifile
integer,          intent(in) :: number_of_levels_temp
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) number_of_levels_temp
   CASE DEFAULT
      write(ifile, *) number_of_levels_temp
END SELECT

end subroutine write_iasi_num_levels


!-----------------------------------------------------------------------
!>

function read_iasi_avg_kernels(ifile, nlevels, fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_avg_kernels(IASI_DIM)

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_avg_kernels(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_avg_kernels(1:nlevels)
END SELECT

end function read_iasi_avg_kernels

!-----------------------------------------------------------------------
!>

function read_iasi_heights(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_heights(IASI_DIM)

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_heights(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_heights(1:nlevels)
END SELECT

end function read_iasi_heights

!-----------------------------------------------------------------------
!>

function read_iasi_pressure(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_pressure(IASI_DIM)

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))
SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_pressure(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_pressure(1:nlevels)
END SELECT

end function read_iasi_pressure

!-----------------------------------------------------------------------
!>

function read_iasi_prior_prof(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_prior_prof(IASI_DIM)

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_prior_prof(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_prior_prof(1:nlevels)
END SELECT

end function read_iasi_prior_prof

!-----------------------------------------------------------------------
!>

function read_iasi_air_column(ifile, nlevels,fform)

integer,                    intent(in) :: ifile, nlevels
character(len=*), optional, intent(in) :: fform
real(r8)                               :: read_iasi_air_column(IASI_DIM)

character(len=32) :: fileformat

fileformat = "ascii"    ! supply default

if(present(fform)) fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) read_iasi_air_column(1:nlevels)
   CASE DEFAULT
      read(ifile, *) read_iasi_air_column(1:nlevels)
END SELECT

end function read_iasi_air_column

!-----------------------------------------------------------------------
!>

subroutine write_iasi_avg_kernels(ifile, avg_kernels_temp, nlevels, fform)

integer,          intent(in) :: ifile, nlevels
real(r8),         intent(in) :: avg_kernels_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) avg_kernels_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) avg_kernels_temp(1:nlevels)
END SELECT

end subroutine write_iasi_avg_kernels

!-----------------------------------------------------------------------
!>

subroutine write_iasi_heights(ifile, height_temp, nlevels, fform)

integer,          intent(in) :: ifile, nlevels
real(r8),         intent(in) :: height_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) height_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) height_temp(1:nlevels)
END SELECT

end subroutine write_iasi_heights

!-----------------------------------------------------------------------
!>

subroutine write_iasi_pressure(ifile, pressure_temp, nlevels, fform)

integer,          intent(in) :: ifile, nlevels
real(r8),         intent(in)  :: pressure_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) pressure_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) pressure_temp(1:nlevels)
END SELECT

end subroutine write_iasi_pressure

!-----------------------------------------------------------------------
!>

subroutine write_iasi_prior_prof(ifile, prior_prof_temp, nlevels, fform)

integer,          intent(in) :: ifile, nlevels
real(r8),         intent(in) :: prior_prof_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) prior_prof_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) prior_prof_temp(1:nlevels)
END SELECT

end subroutine write_iasi_prior_prof

!-----------------------------------------------------------------------
!>

subroutine write_iasi_air_column(ifile, aircol_prof_temp, nlevels, fform)

integer,          intent(in) :: ifile, nlevels
real(r8),         intent(in) :: aircol_prof_temp(IASI_DIM)
character(len=*), intent(in) :: fform

character(len=32) :: fileformat

fileformat = trim(adjustl(fform))

SELECT CASE (fileformat)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) aircol_prof_temp(1:nlevels)
   CASE DEFAULT
      write(ifile, *) aircol_prof_temp(1:nlevels)
END SELECT

end subroutine write_iasi_air_column

end module obs_def_iasi_O3_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
