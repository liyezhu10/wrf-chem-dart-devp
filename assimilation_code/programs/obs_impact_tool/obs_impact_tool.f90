! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> This program assists in constructing a table which can be read
!> by filter at run-time to disable or alter how the assimilation
!> of different types of observations impact state vector values
!> based on their quantity.  This tool allows users to group related
!> collections of observation types and state vector quantities by name
!> and then express the relationship of the named groups to each
!> other in a concise way.
!>
!> Normally filter determines the impact of an observation
!> on the state based on the covariance and localization values.
!> The impact factor is an additional multiplier in this process.
!> For example, if the model does not
!> accurately represent the relationship between different 
!> parts of the model state (e.g part of the state is
!> computed offline, or two different models are run
!> which do not exchange information otherwise), then
!> this tool can be used to prevent all impact of
!> an observation on those parts of the state.
!>
!> We recommend initially only using values of 0.0 or 1.0,
!> although other values can be used after careful analysis
!> of the results.
!>
!> All the listed observation types and state vector quantities
!> must be known by the system.  If they are not, look at the
!> &preprocess_nml :: input_items namelist which specifies
!> which obs_def_xxx_mod.f90 files are included, which is
!> where observation types are defined.  Quantities are defined
!> in the assimilation_code/modules/observations/DEFAULT_obs_kinds_mod.F90 file.
!> (Note you must add new quantities in 2 places 
!> if you do alter this file.)
!>



! program to read an ascii file with directions for which state and observation
! QTYS should impact which other state and observation QTYS.

! the format of the ascii input file is:
!
! # rest of line is comment after hash mark
! GROUP groupname1
!  QTY_xxx  QTY_xxx  QTY_xxx
!  QTY_xxx
! END GROUP
!
! GROUP groupname2
!  QTY_xxx  QTY_xxx  QTY_xxx
!  QTY_xxx
! END GROUP
!
! GROUP groupnameM
!  ALL EXCEPT QTY_xxx QTY_xxx
!  QTY_xxx
! END GROUP
!
! # to choose all quantities except a select few
! GROUP groupnameN
!  ALL EXCEPT groupnameY
! END GROUP
!
! also ALLTYPES, ALLQTYS, as well as ALL
!
! IMPACT
!  QTY_xxx     QTY_xxx      0.0
!  QTY_xxx     groupname1   0.0
!  groupname1  QTY_xxx      0.0
!  groupname1  groupname2   0.0
! END IMPACT

! GROUP groupnameX
!  different_groupname  # no circular dependencies allowed
! END GROUP

! the output of this tool is a ascii file containing lines:
!  QTY1_string   QTY2_string    0.0
!

program obs_impact_tool

use      types_mod, only : r8
use  utilities_mod, only : register_module, initialize_utilities, finalize_utilities, &
                           find_namelist_in_file, check_namelist_read, E_MSG,         &
                           do_nml_file, do_nml_term, nmlfileunit, error_handler
use obs_impact_mod, only : create_impact_table

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"
character(len=128), parameter :: id  = "$Id$"

integer :: funit, ios

! namelist: input/output names, values, etc
character(len=512) :: input_filename  = ''
character(len=512) :: output_filename = ''
logical :: debug = .false.  ! .true. for more output

! namelist
namelist /obs_impact_tool_nml/  &
   input_filename,  &
   output_filename, &
   debug


! initialization and setup

call initialize_utilities('obs_impact_tool')
call register_module(source,revision,revdate)

call find_namelist_in_file("input.nml", "obs_impact_tool_nml", funit)
read(funit, nml = obs_impact_tool_nml, iostat = ios)
call check_namelist_read(funit, ios, "obs_impact_tool_nml")

if (do_nml_file()) write(nmlfileunit, nml=obs_impact_tool_nml)
if (do_nml_term()) write(     *     , nml=obs_impact_tool_nml)

if (debug) call error_handler(E_MSG, 'obs_impact_tool', ' debug on')


! build and output impact_table
call create_impact_table(input_filename, output_filename, debug) 

! clean up
call finalize_utilities('obs_impact_tool')

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
