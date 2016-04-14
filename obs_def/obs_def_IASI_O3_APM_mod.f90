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

call set_obs_def_iasi_o3(key,prior_1,altitude_1(1:nlevel_1),avg_kernel_1(1:nlevel_1), &
                         aircol_1(1:nlevel_1),nlevel_1)

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
real(r8) :: avg_kernel_1(     IASI_DIM)
real(r8) :: iasi_air_column_1(IASI_DIM)

if ( .not. module_initialized ) call initialize_module

fileformat = "ascii"   ! supply default
if(present(fform)) fileformat = trim(adjustl(fform))

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading

nlevel_1                      = iasi_nlevels(   key)
altitude_1(       1:nlevel_1) = iasi_heights(   key,:)
avg_kernel_1(     1:nlevel_1) = avg_kernel(     key,:)
iasi_air_column_1(1:nlevel_1) = iasi_air_column(key,:)

call write_iasi_num_levels(  ifile, nlevel_1, fileformat)
call write_iasi_prior_column(ifile, iasi_o3_prior(key), fileformat)
call write_iasi_heights(     ifile, altitude_1(       1:nlevel_1), nlevel_1, fileformat)
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

type(location_type) :: loc2
real(r8)            :: mloc(3)
real(r8)            :: obs_val, tmp_val, level
integer             :: i, ilev, valid_lev, nlevels

if ( .not. module_initialized ) call initialize_module

mloc = get_location(location)

val = 0.0_r8
if (mloc(2)>90.0_r8) then
   mloc(2)=90.0_r8
elseif (mloc(2)<-90.0_r8) then
   mloc(2)=-90.0_r8
endif

! Apply IASI Averaging kernel A and IASI Prior (I-A)xa
! x = Axm + (I-A)xa , where x is a 40-element vector

nlevels   = iasi_nlevels(key)
obs_val   = 0.0_r8
istatus   = 0
valid_lev = 0
level     = 1.0_r8

! loop through all levels
do ilev = 1, nlevels

   !get location of obs
   if (ilev == 1) then
      loc2 = set_location(mloc(1),mloc(2),level, VERTISLEVEL)
   else
      loc2 = set_location(mloc(1),mloc(2),iasi_heights(key,ilev), VERTISHEIGHT)
   endif

   !initialize obs_val
   obs_val = 0.0_r8
   istatus = 0

   ! interpolate to obs location
   call interpolate(state, loc2, KIND_O3, obs_val, istatus)

   ! check for problems with the interpolation
   if (istatus /= 0) then
      if (istatus == 2) then

         ! it looks like it's a problem with vertical level interpolation
         ! use the prior iasi sub column instead
         ! APM: need to remove this because we don't keep the prior and not realistic
         tmp_val = 0.0_r8
         istatus = 0
      else

         ! interpolation failed (outside of domain?)
         tmp_val = 0.0_r8
         istatus = istatus
         return
      endif
   else

      ! assign interpolated subcolumn (multiply by air column and 1000)
      if ( obs_val > 0.0_r8 ) then
         tmp_val = obs_val * iasi_air_column(key,ilev) * 1000.0_r8
         istatus = 0
      else

         ! the interpolated value cannot be negative
         tmp_val = 0.0_r8
         istatus = 25
         return
      endif
   endif

   ! update expected obs with the averaging kernel
   val = val + avg_kernel(key,ilev) * tmp_val

enddo

end subroutine get_expected_iasi_o3

!-----------------------------------------------------------------------
!>

subroutine set_obs_def_iasi_o3(key, apcol_val, altretlev, akcol, aircol_val, nlev_use)

integer,  intent(in) :: key
integer,  intent(in) :: nlev_use
real(r8), intent(in) :: apcol_val
real(r8), intent(in) :: altretlev( IASI_DIM)
real(r8), intent(in) :: akcol(     IASI_DIM)
real(r8), intent(in) :: aircol_val(IASI_DIM)

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
