! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module cov_cutoff_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                          do_output, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read, &
                          file_exist, get_unit
use location_mod,  only : location_type, get_location
use  obs_kind_mod, only : get_obs_kind_name, get_raw_obs_kind_name,         &
                          KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT,     &
                          KIND_SURFACE_PRESSURE, KIND_TEMPERATURE

implicit none
private

public :: comp_cov_factor

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"


!============================================================================

!---- namelist with default values
logical :: namelist_initialized = .false.

integer :: select_localization = 1
! Value 1 selects default Gaspari-Cohn cutoff
! Value 2 selects boxcar
! Value 3 selects ramped boxcar

namelist / cov_cutoff_nml / select_localization

! Lili (2013-04-18)
logical :: first_get_ELF_horizontal = .true.
logical :: first_get_ELF_vertical   = .true.
logical :: first_get_precip_loc     = .true.
logical :: first_print_T_OP = .true.
logical :: first_print_T_ON = .true.
logical :: first_print_W_OP = .true.
logical :: first_print_W_ON = .true.
character(len = 129) :: elf_file_horizontal_name = "elf_horizontal.in"
character(len = 129) :: elf_file_vertical_name   = "elf_vertical.in"
character(len = 129) :: precip_loc_name = 'precip_loc.dat'
character(len = 129) :: errstring
real(r8), allocatable   :: elf_horizontal(:,:)
real(r8), allocatable   :: elf_vertical(:,:)
real(r8), allocatable   :: precip_loc(:,:)
integer   :: num_precip_loc

!============================================================================

contains

!======================================================================



function comp_cov_factor(z_in, c_in, obs_loc, obs_kind, target_loc, target_kind, &
   localization_override)
!----------------------------------------------------------------------
! function comp_cov_factor(z_in, c)
!
! Computes a covariance cutoff function from Gaspari and Cohn
! QJRMS, 125, 723-757.  (their eqn. 4.10)
!
! z_in is the distance while c is the cutoff distance. 
! For distances greater than 2c, the cov_factor returned goes to 0.

! Other ramping shapes are also available and can be selected by a namelist
! parameter. At present, these include a boxcar with the given halfwidth
! and a ramp in which the weight is set at 1.0 up to the half-width 
! distance and then decreases linearly to 0 at twice the half-width 
! distance.

! Additional information is passed in about the location and kind of the
! observation and the location and kind of the variable being targeted for
! increments. These can be used for more refined algorithms that want to 
! make the cutoff a function of these additional arguments. 

implicit none

real(r8),                      intent(in) :: z_in, c_in
type(location_type), optional, intent(in) :: obs_loc, target_loc
integer,             optional, intent(in) :: obs_kind, target_kind
integer,             optional, intent(in) :: localization_override
!integer,             optional, intent(in) :: ilatindex
real(r8)                                  :: comp_cov_factor

real(r8) :: z, r
integer  :: iunit, io
integer  :: localization_type

integer  :: i, j, k, pindex
real(r8) :: k1, k2, k3, k4, k5
real(r8) :: b1, b2, b3, b4, b5
real(r8) :: c_vert, z_joint, z_adjust, height_norm
real(r8) :: c
real(r8), dimension(3) :: obs_array, var_array

logical  :: obs_with_precip

!--------------------------------------------------------
! Initialize namelist if not already done
if(.not. namelist_initialized) then

   call register_module(source, revision, revdate)

   namelist_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "cov_cutoff_nml", iunit)
   read(iunit, nml = cov_cutoff_nml, iostat = io)
   call check_namelist_read(iunit, io, "cov_cutoff_nml")

   if (do_nml_file()) write(nmlfileunit,nml=cov_cutoff_nml)
   if (do_nml_term()) write(     *     ,nml=cov_cutoff_nml)


   if (do_output()) then
      select case (select_localization)
         case (1)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Standard Gaspari Cohn localization selected')
         case (2)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Boxcar localization selected')
         case (3)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Ramped localization selected')
         case default
            call error_handler(E_ERR,'comp_cov_factor', &
               'Illegal value of "select_localization" in cov_cutoff_mod namelist', &
                source, revision, revdate )
      end select
   endif

endif
!---------------------------------------------------------

if(present(localization_override)) then
   localization_type = localization_override
else
   localization_type = select_localization
endif

z = abs(z_in)

c = c_in

!----------------------------------------------------------

!----------------------------------------------------------
! read in locations (lat, lon) with precipitation
if ( first_get_precip_loc ) then
   if ( file_exist(precip_loc_name) ) then
      write(errstring,*) 'precip_loc file exist ', trim(elf_file_vertical_name)
      call error_handler(E_MSG,'comp_cov_factor:', errstring,source,revision,revdate)

      iunit = get_unit()
      open(unit = iunit, file = precip_loc_name )
      read(iunit, *) num_precip_loc
      allocate(precip_loc(num_precip_loc,2))
      do i = 1, num_precip_loc
         read(iunit,*) precip_loc(i,1:2)
      enddo

print *, 'num_precip_loc = ', num_precip_loc
do i = 1, num_precip_loc
   print *, precip_loc(i,1:2)
enddo
      close(iunit)
   else
      write(errstring,*) 'CANNOT find precip_loc file ', trim(elf_file_vertical_name)
      call error_handler(E_ERR,'comp_cov_factor:',errstring,source,revision,revdate)
   endif
   first_get_precip_loc = .false.
endif

!----------------------------------------------------------

!----------------------------------------------------------
! whether the obs has precipitation
obs_with_precip = .false.

obs_array = get_location(obs_loc)   !lon,lat,vloc
!print *, 'obs_array = ', obs_array
do i = 1, num_precip_loc
   if ( abs(obs_array(1)-precip_loc(i,1))<0.01_r8 .and. abs(obs_array(2)-precip_loc(i,2))<0.01_r8 ) then
       obs_with_precip = .true.
   endif
enddo 

!print *, 'obs with precipitation ? ', obs_with_precip

if ( obs_with_precip == .true. ) then
    if ( obs_kind == KIND_TEMPERATURE ) then
       c = 0.07_r8
       c_vert  = 0.10_r8
       z_joint = 3000.0_r8
       height_norm = 80000.0_r8
       k1 = 3500.0
       k2 = 7000.0
       k3 = 10500.0
       k4 = 14000.0
       b1 = 0.0941e-11
       b2 = 0.5622e-11
       b3 = -0.6231e-11
       b4 = 0.2198e-11
       if ( first_print_T_OP ) then
          print *, 'T_OP, c = ', c, 'betas = ', b1, b2, b3, b4
          first_print_T_OP = .false.
       endif
    else
       c = 0.06_r8
       c_vert  = 0.08_r8
       z_joint = 3000.0_r8
       height_norm = 80000.0_r8
       k1 = 3500.0
       k2 = 7000.0
       k3 = 10500.0
       k4 = 14000.0
       b1 = -0.0140e-10
       b2 = 0.1038e-10
       b3 = -0.0930e-10
       b4 = 0.0291e-10
       if ( first_print_W_OP ) then
          print *, 'W_OP, c = ', c, 'betas = ', b1, b2, b3, b4
          first_print_W_OP = .false.
       endif
    endif
else
    if ( obs_kind == KIND_TEMPERATURE ) then
       c = 0.1_r8
       c_vert  = 0.06_r8
       z_joint = 2000.0_r8
       height_norm = 80000.0_r8
       k1 = 3500.0
       k2 = 7000.0
       k3 = 10500.0
       k4 = 14000.0
       b1 = -0.7995e-11
       b2 = 0.9136e-11
       b3 = -0.7145e-11
       b4 = 0.2277e-11
       if ( first_print_T_ON ) then
          print *, 'T_ON, c = ', c, 'betas = ', b1, b2, b3, b4
          first_print_T_ON = .false.
       endif
    else
       c = 0.09_r8
       c_vert  = 0.06_r8
       z_joint = 2000.0_r8
       height_norm = 80000.0_r8
       k1 = 3500.0
       k2 = 7000.0
       k3 = 10500.0
       k4 = 14000.0
       b1 = -0.6003e-11
       b2 = 0.7613e-11
       b3 = -0.5701e-11
       b4 = 0.1827e-11
       if ( first_print_W_ON ) then
          print *, 'W_ON, c = ', c, 'betas = ', b1, b2, b3, b4
          first_print_W_ON = .false.
       endif
    endif
endif

!print *, 'cutoff = ', c, 'betas = ', b1, b2, b3, b4 

!----------------------------------------------------------

if(localization_type == 1) then ! Standard Gaspari Cohn localization

   if( z >= c*2.0_r8 ) then

      comp_cov_factor = 0.0_r8

   else if( z <= c ) then
      r = z / c
      comp_cov_factor = &
           ( ( ( -0.25_r8*r +0.5_r8 )*r +0.625_r8 )*r -5.0_r8/3.0_r8 )*r**2 + 1.0_r8
!!$           r**5 * (-0.25_r8 ) + &
!!$           r**4 / 2.0_r8 +              &
!!$           r**3 * 5.0_r8/8.0_r8 -       &
!!$           r**2 * 5.0_r8/3.0_r8 + 1.0_r8
   else

      r = z / c
      comp_cov_factor = &
           ( ( ( ( r/12.0_r8 -0.5_r8 )*r +0.625_r8 )*r +5.0_r8/3.0_r8 )*r -5.0_r8 )*r &
!!$           r**5 / 12.0_r8  -  &
!!$           r**4 / 2.0_r8   +  &
!!$           r**3 * 5.0_r8 / 8.0_r8 + &
!!$           r**2 * 5.0_r8 / 3.0_r8 - 5.0_r8*r &
           + 4.0_r8 - 2.0_r8 / (3.0_r8 * r) 
   endif

else if(localization_type == 2) then ! BOXCAR localization

   if(z < 2.0_r8 * c) then
      comp_cov_factor = 1.0_r8
   else
      comp_cov_factor = 0.0_r8
   endif

else if(localization_type == 3) then ! Ramped localization

   if(z >= 2.0_r8 * c) then
      comp_cov_factor = 0.0_r8
   else if(z >= c .and. z < 2.0_r8 * c) then
      comp_cov_factor = (2.0_r8 * c - z) / c
   else
      comp_cov_factor = 1.0_r8
   endif

else if(localization_type == 4) then ! horizontal cubic spline fitted function

!    ! f(z)=beta1*(knots1-z)^3+beta2*(knots2-z)^3+beta3*(knots3-z)^3... for each
!    ! term, z <= knots
!    if ( first_get_elf_horizontal ) then
!       if ( file_exist(elf_file_horizontal_name) ) then
!          write(errstring,*) 'ELF file of cubicsplinefit exist ', trim(elf_file_horizontal_name)
!          call error_handler(E_MSG,'comp_cov_factor:', errstring,source,revision,revdate)
! 
!          allocate(elf_horizontal(6,3))
!          iunit = get_unit()
!          open(unit = iunit, file = elf_file_horizontal_name)
!          ! first row: beta, SH
!          read(iunit, *) elf_horizontal(1,1:3)
!          ! second row: knots, SH
!          read(iunit, *) elf_horizontal(2,1:3)
!          read(iunit, *) elf_horizontal(3,1:3)  ! beta, TP
!          read(iunit, *) elf_horizontal(4,1:3)  ! knot, TP
!          read(iunit, *) elf_horizontal(5,1:3)  ! beta, NH
!          read(iunit, *) elf_horizontal(6,1:3)  ! knot, NH
!do i = 1, 6
!   print *, elf_horizontal(i,1:3)
!enddo
!          close(iunit)
!       else
!          write(errstring,*) 'CANNOT find ELF file of cubicsplinefit ', trim(elf_file_horizontal_name)
!          call error_handler(E_ERR,'comp_cov_factor:',errstring,source,revision,revdate)
!       endif
!       first_get_elf_horizontal = .false.
!    endif
!    b1 = elf_horizontal((ilatindex-1)*2+1,1)
!    b2 = elf_horizontal((ilatindex-1)*2+1,2)
!    b3 = elf_horizontal((ilatindex-1)*2+1,3)
!    k1 = elf_horizontal((ilatindex-1)*2+2,1)
!    k2 = elf_horizontal((ilatindex-1)*2+2,2)
!    k3 = elf_horizontal((ilatindex-1)*2+2,3)
!
!    if ( z <= k1 ) then
!       comp_cov_factor = b1*(k1-z)**3 + b2*(k2-z)**3 + b3*(k3-z)**3
!    else if ( z > k1 .and. z <= k2 ) then
!       comp_cov_factor = b2*(k2-z)**3 + b3*(k3-z)**3
!    else if ( z > k2 .and. z <= k3 ) then
!       comp_cov_factor = b3*(k3-z)**3
!    else
!       comp_cov_factor = 0.0_r8
!    endif
!    comp_cov_factor = min(comp_cov_factor,1.0_r8)

else if(localization_type == 5) then ! vertical cubic spline fitted function

    ! f(z)=beta1*(knots1-z)^3+beta2*(knots2-z)^3+beta3*(knots3-z)^3... for each
    ! term, z <= knots

    if ( z > z_joint ) then
       z_adjust = z - z_joint
       if ( z_adjust <= k1 ) then
          comp_cov_factor = b1*(k1-z_adjust)**3 + b2*(k2-z_adjust)**3 + b3*(k3-z_adjust)**3 + b4*(k4-z_adjust)**3
       else if ( z_adjust > k1 .and. z_adjust <= k2 ) then
          comp_cov_factor = b2*(k2-z_adjust)**3 + b3*(k3-z_adjust)**3 + b4*(k4-z_adjust)**3
       else if ( z_adjust > k2 .and. z_adjust <= k3 ) then
          comp_cov_factor = b3*(k3-z_adjust)**3 + b4*(k4-z_adjust)**3
       else if ( z_adjust > k3 .and. z_adjust <= k4 ) then
          comp_cov_factor = b4*(k4-z_adjust)**3
       else
          comp_cov_factor = 0.0_r8
       endif
!       comp_cov_factor = min(comp_cov_factor,1.0_r8)
    else
       z_adjust = z/height_norm
       if( z_adjust >= c_vert*2.0_r8 ) then
          comp_cov_factor = 0.0_r8
       else if( z_adjust <= c_vert ) then
          r = z_adjust / c_vert
          comp_cov_factor = &
               ( ( ( -0.25_r8*r +0.5_r8 )*r +0.625_r8 )*r -5.0_r8/3.0_r8 )*r**2 + 1.0_r8
       else
          r = z_adjust / c_vert
          comp_cov_factor = &
               ( ( ( ( r/12.0_r8 -0.5_r8 )*r +0.625_r8 )*r +5.0_r8/3.0_r8 )*r-5.0_r8 )*r &
               + 4.0_r8 - 2.0_r8 / (3.0_r8 * r)
       endif
    endif

    comp_cov_factor = min(comp_cov_factor,1.0_r8)

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(E_ERR,'comp_cov_factor', &
              'Illegal value of "localization" in cov_cutoff_mod namelist', &
               source, revision, revdate )

endif

end function comp_cov_factor

end module cov_cutoff_mod
