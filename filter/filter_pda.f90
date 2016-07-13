! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: filter.f90 7492 2015-01-27 22:36:58Z hkershaw $

!> \dir filter  Main program contained here
!> \file filter.f90 Main program

program filter

use filter_mod

implicit none

!----------------------------------------------------------------

!call filter_main()
call pda_main()

end program filter

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/pda/filter/filter.f90 $
! $Id: filter.f90 7492 2015-01-27 22:36:58Z hkershaw $
! $Revision: 7492 $
! $Date: 2015-01-27 15:36:58 -0700 (Tue, 27 Jan 2015) $
