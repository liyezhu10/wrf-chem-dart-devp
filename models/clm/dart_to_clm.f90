! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program dart_to_clm

!----------------------------------------------------------------------
! purpose: interface between DART and the CLM model
!
! method: Read DART state vector and overwrite values in a CLM restart file.
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG
use        model_mod, only : static_init_model

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

integer            :: iunit, io
character(len=512) :: string1, string2, string3

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256) :: input_dart_file = 'dart_posterior.nc'
character(len=256) :: clm_file_to_update = 'clm_restart_file.nc'

namelist /dart_to_clm_nml/ input_dart_file, clm_file_to_update

!======================================================================

call initialize_utilities(progname='dart_to_clm')

call static_init_model()

call find_namelist_in_file("input.nml", "dart_to_clm_nml", iunit)
read(iunit, nml = dart_to_clm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_clm_nml")

write(string1,*)'converting DART file "'//trim(input_dart_file)//'"'
write(string2,*)'to clm restart file "'//trim(clm_file_to_update)//'"'
call error_handler(E_MSG,'dart_to_clm',string1,text2=string2)

! this is where something useful might happen ... like rebalancing the
! SWE into snow layers or making the soil temperature and moisture compatible ...

call finalize_utilities('dart_to_clm')

end program dart_to_clm

