! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module options_mod

! Aim: 
! 
! 

use types_mod, only : r8

implicit none

public :: get_missing_ok_status, set_missing_ok_status

logical :: allow_missing_r8 = .FALSE.

contains

!------------------------------------------------------------------------

subroutine set_missing_ok_status(allow_missing)
logical, intent(in) :: allow_missing

allow_missing_r8 = allow_missing

end subroutine set_missing_ok_status

!------------------------------------------------------------------------

function get_missing_ok_status()
logical :: get_missing_ok_status

get_missing_ok_status = allow_missing_r8

end function get_missing_ok_status

!------------------------------------------------------------------------


!-------------------------------------------------------
!> Null version
!> Check whether you need to error out, clamp, or
!> do nothing depending on the variable bounds

end module options_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

