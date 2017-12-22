! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module netcdf_utilities_mod

!>@todo FIXME: NEEDS MUCH MORE WORK.
!>
!> start of a module to help streamline the code that
!> calls netcdf routines.  usually adds error checking.
!> 
!> also isolates as much access to the netcdf libs to routines
!> inside this module.
!>
!> routines in this file are prefixed with nc_
!>
!> the intent is that these are simple routines that aren't trying to
!> handle all possible options.  if you have complicated needs, write your
!> own specialized routine, either here or in the calling code.
!> for example, the 'put_var' routine here doesn't have start or count, 
!> intentionally, for array entries.


use types_mod, only : r4, r8, digits12, i4, i8, PI, MISSING_R8, MISSING_I
use utilities_mod, only : error_handler, E_DBG, E_MSG, E_ALLMSG, E_WARN, E_ERR
           
use netcdf

implicit none
private

public :: nc_check,                       &
          nc_add_global_attribute,        &
          nc_add_attribute_to_variable,   &
          nc_define_dimension,            &
          nc_define_integer_variable,     &
          nc_define_real_variable,        &
          nc_define_double_variable,      &
          nc_global_attribute_exists,     &
          nc_dimension_exists,            &
          nc_variable_exists,             &
          nc_put_variable,                &
          nc_get_variable,                &
          nc_add_global_creation_time,    &
          nc_get_variable_num_dimensions, &
          nc_get_variable_size,           &
          nc_open_file_readonly,          &
          nc_open_file_readwrite,         &
          nc_create_file,                 &
          nc_close_file,                  &
          nc_begin_define_mode,           &
          nc_end_define_mode,             &
          nc_synchronize_file


interface nc_add_global_attribute
   module procedure nc_add_global_char_att
   module procedure nc_add_global_int_att
   module procedure nc_add_global_real_att
   module procedure nc_add_global_real_array_att
end interface

interface nc_add_attribute_to_variable
   module procedure nc_add_char_att_to_var
   module procedure nc_add_int_array_att_to_var
   module procedure nc_add_int_att_to_var
   module procedure nc_add_real_att_to_var
   module procedure nc_add_real_array_att_to_var
end interface

interface nc_define_integer_variable
   module procedure nc_define_var_int_scalar
   module procedure nc_define_var_int_1d
   module procedure nc_define_var_int_Nd
end interface

interface nc_define_real_variable
   module procedure nc_define_var_real_scalar
   module procedure nc_define_var_real_1d
   module procedure nc_define_var_real_Nd
end interface

interface nc_define_double_variable
   module procedure nc_define_var_double_scalar
   module procedure nc_define_var_double_1d
   module procedure nc_define_var_double_Nd
end interface

interface nc_put_variable
   module procedure nc_put_single_int_1d
   module procedure nc_put_int_1d
   module procedure nc_put_single_real_1d
   module procedure nc_put_real_1d
   module procedure nc_put_int_2d
   module procedure nc_put_real_2d
   module procedure nc_put_int_3d
   module procedure nc_put_real_3d
end interface

interface nc_get_variable
   module procedure nc_get_single_int_1d
   module procedure nc_get_int_1d
   module procedure nc_get_single_real_1d
   module procedure nc_get_real_1d
   module procedure nc_get_int_2d
   module procedure nc_get_real_2d
   module procedure nc_get_int_3d
   module procedure nc_get_real_3d
end interface

interface nc_get_variable_size
   module procedure nc_get_variable_size_1d
   module procedure nc_get_variable_size_Nd
end interface

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring1, msgstring2, msgstring3

!> make a derived type that is (ncid, filename) 
!> store filename on file open, delete it on file close. cache the
!> last N filenames - look them up on error and stop
!> having to keep the filename around.

!> NOTE!!
!> this assumes that you are doing the open/read/write/close
!> operations on the same task.  which i believe is true for all
!> our current code.

integer, parameter :: MAX_NCFILES = 10
integer, parameter :: FH_EMPTY = -1

type ncinfo_type
  integer :: file_handle = FH_EMPTY
  character(len=256) :: file_name = ''
end type

! for now hardcode max size.  could make it allocatable
! and extend it if you need more open slots

type(ncinfo_type) :: ncinfo(MAX_NCFILES)

! do we need one of these?
!namelist /netcdf_utilities_nml/ 

contains

!------------------------------------------------------------------
!> check return code from previous call. on error, print and stop.
!> if you want to continue after an error don't use this call. 

subroutine nc_check(istatus, subr_name, context, context2, filename, ncid)

integer,          intent(in)           :: istatus
character(len=*), intent(in)           :: subr_name
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: context2
character(len=*), intent(in), optional :: filename
integer,          intent(in), optional :: ncid
  
character(len=256) :: saved_filename

if (istatus == nf90_noerr) return

! something wrong.  construct an error string and call the handler.
msgstring1 = nf90_strerror(istatus)

! context is optional, but is very useful if specified.
if (present(context)) then
  msgstring1 = trim(context) // ': ' // trim(msgstring1)
endif

! you can give this routine a file handle and it will try
! to extract the filename from it.
if (.not. present(filename) .and. present(ncid)) then
   call find_name_from_fh(ncid, saved_filename)
else if (present(filename)) then
   saved_filename = filename
else
   saved_filename = ''
endif

! this does not return 
call error_handler(E_ERR, subr_name, msgstring1, source, revision, revdate, &
                   text2=context2, text3=saved_filename)
  

end subroutine nc_check

!------------------------------------------------------------------
!--------------------------------------------------------------------
! global attribute section

subroutine nc_add_global_char_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_char_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_char_att

!--------------------------------------------------------------------

subroutine nc_add_global_int_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
integer,          intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_int_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_int_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_real_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_real_att

!--------------------------------------------------------------------

subroutine nc_add_global_real_array_att(ncid, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_global_real_array_att'
integer :: ret

ret = nf90_put_att(ncid, NF90_GLOBAL, attname, val)
call nc_check(ret, routine, 'adding the global attribute: '//trim(attname), context, filename, ncid)

end subroutine nc_add_global_real_array_att

!------------------------------------------------------------------
!--------------------------------------------------------------------
! adding attributes to variables section

subroutine nc_add_char_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
character(len=*), intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_char_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_char_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_int_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
integer,          intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_int_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_int_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_int_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
integer,          intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_int_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_int_array_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_real_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_real_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_real_att_to_var

!--------------------------------------------------------------------

subroutine nc_add_real_array_att_to_var(ncid, varname, attname, val, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: attname
real(r8),         intent(in) :: val(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_add_real_array_att_to_var'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_att(ncid, varid, attname, val)
call nc_check(ret, routine, 'adding the attribute: '//trim(attname)//' to variable: '//trim(varname), context, filename, ncid)

end subroutine nc_add_real_array_att_to_var

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! defining dimensions section

subroutine nc_define_dimension(ncid, dimname, dimlen, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
integer,          intent(in) :: dimlen
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_dimension'
integer :: ret, dimid

ret = nf90_def_dim(ncid, dimname, dimlen, dimid)
call nc_check(ret, routine, 'define dimension '//trim(dimname), context, filename, ncid)

end subroutine nc_define_dimension

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! defining variables section

!>@todo FIXME: ndim is not required or used here, but without it
!>the signatures for scalar, 1d, Nd cannot be distinguished.
!> same for define real, double scalar.  bother.
!> (the optional string arguments make certain combinations of
!> arguments identical)

subroutine nc_define_var_int_scalar(ncid, varname, ndim, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: ndim ! ignore for now
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_int, varid=varid)
call nc_check(ret, routine, 'define int variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_int_scalar

!--------------------------------------------------------------------

subroutine nc_define_var_int_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_int, dimid, varid)
call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename, ncid)

print *, trim(varname), varid, trim(dimname), dimid, ret

end subroutine nc_define_var_int_1d

!--------------------------------------------------------------------

subroutine nc_define_var_int_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_int_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename, ncid)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename, ncid)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename, ncid)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d integer variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_int, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define integer variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_int_Nd

!--------------------------------------------------------------------

subroutine nc_define_var_real_scalar(ncid, varname, ndim, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: ndim ! ignore for now
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_real, varid=varid)
call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_real_scalar

!--------------------------------------------------------------------

subroutine nc_define_var_real_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_real, dimid, varid)
call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_real_1d

!--------------------------------------------------------------------

subroutine nc_define_var_real_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_real_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename, ncid)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename, ncid)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename, ncid)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d real variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_real, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define real variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_real_Nd

!--------------------------------------------------------------------

subroutine nc_define_var_double_scalar(ncid, varname, ndim, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: ndim ! ignore for now
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_scalar'
integer :: ret, varid

ret = nf90_def_var(ncid, varname, nf90_double, varid=varid)
call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_double_scalar

!--------------------------------------------------------------------

subroutine nc_define_var_double_1d(ncid, varname, dimname, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimname
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_1d'
integer :: ret, dimid, varid

ret = nf90_inq_dimid(ncid, dimname, dimid)
call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimname), context, filename, ncid)

ret = nf90_def_var(ncid, varname, nf90_double, dimid, varid)
call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_double_1d

!--------------------------------------------------------------------

subroutine nc_define_var_double_Nd(ncid, varname, dimnames, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: dimnames(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_define_var_double_Nd'
integer :: ret, dimid1, dimid2, dimid3, varid

if (size(dimnames) >= 1) then
   ret = nf90_inq_dimid(ncid, dimnames(1), dimid1)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(1)), context, filename, ncid)
endif

if (size(dimnames) >= 2) then
   ret = nf90_inq_dimid(ncid, dimnames(2), dimid2)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(2)), context, filename, ncid)
endif

if (size(dimnames) >= 3) then
   ret = nf90_inq_dimid(ncid, dimnames(3), dimid3)
   call nc_check(ret, routine, 'inquire dimension id for dim '//trim(dimnames(3)), context, filename, ncid)
endif

if (size(dimnames) >= 4) then
   call error_handler(E_ERR, routine, 'only 1d, 2d and 3d double variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

if (size(dimnames) == 1) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimid1, varid=varid)
else if (size(dimnames) == 2) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimids=(/ dimid1, dimid2 /), varid=varid)
else if (size(dimnames) == 3) then
   ret = nf90_def_var(ncid, varname, nf90_double, dimids=(/ dimid1, dimid2, dimid3 /), varid=varid)
endif

call nc_check(ret, routine, 'define double variable '//trim(varname), context, filename, ncid)

end subroutine nc_define_var_double_Nd

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! check if vars, dims, or global atts exist (without error if not)
! these are functions, unlike the rest of these routines.

function nc_global_attribute_exists(ncid, attname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: attname
logical                      :: nc_global_attribute_exists

character(len=*), parameter :: routine = 'nc_global_attribute_exists'
integer :: ret

ret = nf90_inquire_attribute(ncid, NF90_GLOBAL, attname)
if (ret == NF90_NOERR) then
   nc_global_attribute_exists = .true.
else
   nc_global_attribute_exists = .false.
endif

end function nc_global_attribute_exists

!--------------------------------------------------------------------

function nc_dimension_exists(ncid, dimname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: dimname
logical                      :: nc_dimension_exists

character(len=*), parameter :: routine = 'nc_dimension_exists'
integer :: ret, dimid

ret = nf90_inq_dimid(ncid, dimname, dimid)
if (ret == NF90_NOERR) then
   nc_dimension_exists = .true.
else
   nc_dimension_exists = .false.
endif

end function nc_dimension_exists

!--------------------------------------------------------------------

function nc_variable_exists(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
logical                      :: nc_variable_exists

character(len=*), parameter :: routine = 'nc_variable_exists'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
if (ret == NF90_NOERR) then
   nc_variable_exists = .true.
else
   nc_variable_exists = .false.
endif

end function nc_variable_exists

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! put values into variables

subroutine nc_put_single_int_1d(ncid, varname, varindex, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varindex
integer,          intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varval, start=(/varindex/))
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename, ncid)

end subroutine nc_put_single_int_1d

!--------------------------------------------------------------------

subroutine nc_put_int_1d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_1d

!--------------------------------------------------------------------

subroutine nc_put_single_real_1d(ncid, varname, varindex, varval, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varindex
real(r8),         intent(in) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varval, start=(/varindex/))
call nc_check(ret, routine, 'put value for '//trim(varname), context, filename, ncid)

end subroutine nc_put_single_real_1d

!--------------------------------------------------------------------

subroutine nc_put_real_1d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_1d

!--------------------------------------------------------------------

subroutine nc_put_int_2d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_2d

!--------------------------------------------------------------------

subroutine nc_put_real_2d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_2d

!--------------------------------------------------------------------

subroutine nc_put_int_3d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
integer,          intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_int_3d

!--------------------------------------------------------------------

subroutine nc_put_real_3d(ncid, varname, varvals, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname
real(r8),         intent(in) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_put_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_put_var(ncid, varid, varvals)
call nc_check(ret, routine, 'put values for '//trim(varname), context, filename, ncid)

end subroutine nc_put_real_3d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! get values from variables

subroutine nc_get_single_int_1d(ncid, varname, varval, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_single_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varval)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_single_int_1d

!--------------------------------------------------------------------

subroutine nc_get_int_1d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_1d

!--------------------------------------------------------------------

subroutine nc_get_single_real_1d(ncid, varname, varval, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varval
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_single_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varval)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_single_real_1d

!--------------------------------------------------------------------

subroutine nc_get_real_1d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_1d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_1d

!--------------------------------------------------------------------

subroutine nc_get_int_2d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_2d

!--------------------------------------------------------------------

subroutine nc_get_real_2d(ncid, varname, varvals, context, filename)
integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_2d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_2d

!--------------------------------------------------------------------

subroutine nc_get_int_3d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_int_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_int_3d

!--------------------------------------------------------------------

subroutine nc_get_real_3d(ncid, varname, varvals, context, filename)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
real(r8),         intent(out) :: varvals(:,:,:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_real_3d'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_get_var(ncid, varid, varvals)
call nc_check(ret, routine, 'get values for '//trim(varname), context, filename, ncid)

end subroutine nc_get_real_3d

!------------------------------------------------------------------
!--------------------------------------------------------------------
! inquire variable info

subroutine nc_get_variable_size_1d(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_1d'
integer :: ret, varid, dimids(NF90_MAX_VAR_DIMS)

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

ret = nf90_inquire_dimension(ncid, dimids(1), len=varsize)
call nc_check(ret, routine, 'inquire dimension length for dimension 1', context, filename, ncid)

end subroutine nc_get_variable_size_1d

!--------------------------------------------------------------------

subroutine nc_get_variable_size_Nd(ncid, varname, varsize, context, filename)      

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
integer,          intent(out) :: varsize(:)
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_size_Nd'
integer :: ret, varid, dimids(NF90_MAX_VAR_DIMS)

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, dimids=dimids)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

if (size(varsize) >= 1) then
   ret = nf90_inquire_dimension(ncid, dimids(1), len=varsize(1))
   call nc_check(ret, routine, 'inquire dimension 1 length for variable '//trim(varname), context, filename, ncid)
endif

if (size(varsize) >= 2) then
   ret = nf90_inquire_dimension(ncid, dimids(2), len=varsize(2))
   call nc_check(ret, routine, 'inquire dimension 2 length for variable '//trim(varname), context, filename, ncid)
endif

if (size(varsize) >= 3) then
   ret = nf90_inquire_dimension(ncid, dimids(3), len=varsize(3))
   call nc_check(ret, routine, 'inquire dimension 3 length for variable '//trim(varname), context, filename, ncid)
endif

if (size(varsize) >= 4) then
   ret = nf90_inquire_dimension(ncid, dimids(4), len=varsize(4))
   call nc_check(ret, routine, 'inquire dimension 4 length for variable '//trim(varname), context, filename, ncid)
endif

if (size(varsize) >= 5) then
   ret = nf90_inquire_dimension(ncid, dimids(5), len=varsize(5))
   call nc_check(ret, routine, 'inquire dimension 5 length for variable '//trim(varname), context, filename, ncid)
endif

if (size(varsize) >= 6) then
   call error_handler(E_ERR, routine, 'only 1d to 5d variables supported', &
                      source, revision, revdate, text2='variable '//trim(varname))
endif

end subroutine nc_get_variable_size_Nd

!------------------------------------------------------------------

subroutine nc_get_variable_num_dimensions(ncid, varname, numdims, context, filename) 

integer, intent(in) :: ncid
character(len=*), intent(in):: varname
integer, intent(out) :: numdims
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_get_variable_num_dimensions'
integer :: ret, varid

ret = nf90_inq_varid(ncid, varname, varid)
call nc_check(ret, routine, 'inquire variable id for '//trim(varname), context, filename, ncid)

ret = nf90_inquire_variable(ncid, varid, ndims=numdims)
call nc_check(ret, routine, 'inquire dimensions for variable '//trim(varname), context, filename, ncid)

end subroutine nc_get_variable_num_dimensions 

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! misc section: file operations, standard timestamp routine

subroutine nc_add_global_creation_time(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=128) :: str1

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_add_global_char_att(ncid, "creation_date", str1, context, filename)

end subroutine nc_add_global_creation_time

!--------------------------------------------------------------------

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.
! it is also a function that returns an integer.  (caller doesn't need netcdf)

function nc_open_file_readonly(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readonly

character(len=*), parameter :: routine = 'nc_open_file_readonly'
integer :: ret, ncid

ret = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(ret, routine, 'open '//trim(filename)//' read only', context)

call add_fh_to_list(ncid, filename)
nc_open_file_readonly = ncid

end function nc_open_file_readonly

!--------------------------------------------------------------------

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.

function nc_open_file_readwrite(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_open_file_readwrite

character(len=*), parameter :: routine = 'nc_open_file_readwrite'
integer :: ret, ncid, oldmode

ret = nf90_open(filename, NF90_WRITE, ncid)
call nc_check(ret, routine, 'open '//trim(filename)//' read/write', context)

call add_fh_to_list(ncid, filename)
nc_open_file_readwrite = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)
call nc_check(ret, routine, 'set nofill mode', context, ncid=ncid)

end function nc_open_file_readwrite

!--------------------------------------------------------------------

! the opens are the only routines in which filename is not the last argument.
! all other start with ncid.  this one starts with filename and it's required.

function nc_create_file(filename, context)

character(len=*), intent(in)  :: filename
character(len=*), intent(in), optional :: context
integer                       :: nc_create_file

character(len=*), parameter :: routine = 'nc_create_file'
integer :: ret, ncid, oldmode

ret = nf90_create(filename, NF90_CLOBBER, ncid)
call nc_check(ret, routine, 'create '//trim(filename)//' read/write', context)

call add_fh_to_list(ncid, filename)
nc_create_file = ncid

! faster if we don't fill the vars first with 'fill' value.
! this works if we are planning to write all vars.  (which we are.)

ret = nf90_set_fill(ncid, NF90_NOFILL, oldmode)
call nc_check(ret, routine, 'set nofill mode', context, ncid=ncid)

end function nc_create_file

!--------------------------------------------------------------------

subroutine nc_close_file(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_close_file'
integer :: ret

ret = nf90_close(ncid)
call nc_check(ret, routine, 'close file', context, filename, ncid)

call del_fh_from_list(ncid)

end subroutine nc_close_file

!--------------------------------------------------------------------

subroutine nc_begin_define_mode(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_begin_define_mode'
integer :: ret

ret = nf90_redef(ncid)
call nc_check(ret, routine, 'begin file define mode', context, filename, ncid)

end subroutine nc_begin_define_mode

!--------------------------------------------------------------------

subroutine nc_end_define_mode(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_end_define_mode'
integer :: ret

ret = nf90_enddef(ncid)
call nc_check(ret, routine, 'end file define mode', context, filename, ncid)

end subroutine nc_end_define_mode

!--------------------------------------------------------------------

subroutine nc_synchronize_file(ncid, context, filename)

integer,          intent(in) :: ncid
character(len=*), intent(in), optional :: context
character(len=*), intent(in), optional :: filename

character(len=*), parameter :: routine = 'nc_synchronize_file'
integer :: ret

ret = nf90_sync(ncid)
call nc_check(ret, routine, 'synchronize file contents', context, filename, ncid)

end subroutine nc_synchronize_file

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! internal bookkeeping of mapping between ncids and filenames

!------------------------------------------------------------------
!> add file handle and filename to an available slot

subroutine add_fh_to_list(ncid, filename)
integer, intent(in) :: ncid
character(len=*), intent(in) :: filename

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= FH_EMPTY) cycle

   ncinfo(i)%file_handle = ncid
   ncinfo(i)%file_name = filename
   return 
enddo

end subroutine add_fh_to_list

!------------------------------------------------------------------
!> remove an entry when file is closed

subroutine del_fh_from_list(ncid)
integer, intent(in) :: ncid

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= ncid) cycle

   ncinfo(i)%file_handle = FH_EMPTY
   ncinfo(i)%file_name = ''
   return 
enddo

end subroutine del_fh_from_list

!------------------------------------------------------------------
!> look up and return the filename for this handle

subroutine find_name_from_fh(ncid, filename)
integer, intent(in) :: ncid
character(len=*), intent(out) :: filename

integer :: i

do i=1, MAX_NCFILES
   if (ncinfo(i)%file_handle /= ncid) cycle

   filename = ncinfo(i)%file_name 
   return 
enddo

filename = ''

end subroutine find_name_from_fh

!------------------------------------------------------------------
!------------------------------------------------------------------

end module netcdf_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

