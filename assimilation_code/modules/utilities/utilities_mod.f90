! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module utilities_mod

use types_mod, only : r4, r8, digits12, i4, i8, PI, MISSING_R8, MISSING_I

!>@todo FIXME  netcdf is only needed for the nc_check() routine.
!>when it moves to the netcdf_utilities_mod module remove the 'use'
!>from here to simplify things.
use netcdf

implicit none
private

! module local data

integer, parameter :: E_DBG = -2, E_MSG = -1, E_ALLMSG = 0, E_WARN = 1, E_ERR = 2
integer, parameter :: NML_NONE = 0, NML_FILE = 1, NML_TERMINAL = 2, NML_BOTH = 3

real(r8), parameter :: TWOPI = PI * 2.0_r8

logical :: do_output_flag = .false.
integer :: nml_flag       = NML_FILE
logical :: single_task    = .true.
integer :: task_number    = 0
logical :: module_initialized = .false.
integer :: logfileunit = -1
integer :: nmlfileunit = -1


public :: get_unit, &
          open_file, &
          close_file, &
          file_exist, &
          error_handler, &
          to_upper, &
          squeeze_out_blanks, &
          nc_check, &
          next_file, &   ! deprecate this
          logfileunit, &
          nmlfileunit, &
          find_textfile_dims, &
          file_to_text, &
          timestamp, &
          set_tasknum, &
          set_output, &
          do_output, &
          set_nml_output, &
          E_DBG, &
          E_MSG, &
          E_ALLMSG, &
          E_WARN, &
          E_ERR, &
          is_longitude_between, &
          get_next_filename, &
          ascii_file_format, &
          set_filename_list, &
          set_multiple_filename_lists, &
          scalar, &
          string_to_real, &
          string_to_integer, &
          string_to_logical, &
          find_enclosing_indices, &
          array_dump, &
          dump_unit_attributes, &
          ! lowest level routines
          initialize_utilities, &
          finalize_utilities, &
          register_module, &
          find_namelist_in_file, &
          check_namelist_read, &
          do_nml_file, &
          do_nml_term, &
          log_it

! this routine is either in the null_mpi_utilities_mod.f90, or in
! the mpi_utilities_mod.f90 file, but it is not a module subroutine.
! the mpi files use this module, and you cannot have circular module
! references.  the point of this call is that in the mpi multi-task
! case, you want to call MPI_Abort() to kill the other tasks associated
! with this job when you exit.  in the non-mpi case, it just calls exit.
interface
 subroutine exit_all(exitval)
  integer, intent(in) :: exitval
 end subroutine exit_all
end interface

! make a converter from a 1 element 1D array to a scalar.
! (will this have problems compiling if "integer" is coerced to I8 by
! compiler flags?  making to_scalar_int and to_scalar_int8 the same?
! if so, make to_scalar_int explicitly I4, i guess.)
interface scalar
   module procedure to_scalar_real
   module procedure to_scalar_int
   module procedure to_scalar_int8
end interface

interface array_dump
   module procedure array_1d_dump
   module procedure array_2d_dump
   module procedure array_3d_dump
end interface

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: msgstring1, msgstring2, msgstring3

!----------------------------------------------------------------
! Namelist input with default values

! E_ERR All warnings/errors are assumed fatal.
integer            :: TERMLEVEL      = E_ERR   

! default log and namelist output filenames
character(len=256) :: logfilename    = 'dart_log.out'
character(len=256) :: nmlfilename    = 'dart_log.nml'

! output each module subversion details
logical            :: module_details = .true.  

! print messages labeled DBG
logical            :: print_debug    = .false. 

! where to write namelist values.
! valid strings:  'none', 'file', 'terminal', 'both'
character(len=32)  :: write_nml      = 'file'  

namelist /utilities_nml/ TERMLEVEL, logfilename, module_details, &
                         nmlfilename, print_debug, write_nml

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! base (lowest level) routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>

subroutine initialize_utilities(progname, alternatename, output_flag)
character(len=*), intent(in), optional :: progname
character(len=*), intent(in), optional :: alternatename
logical, intent(in), optional          :: output_flag
integer :: iunit, io

character(len=256) :: lname
character(len=512) :: string1,string2,string3

if ( module_initialized ) return

module_initialized = .true.



! now default to false, and only turn on if i'm task 0
! or the caller tells me to turn it on.
if (present(output_flag)) then
   do_output_flag = output_flag
else
   if (single_task .or. task_number == 0) do_output_flag = .true.
endif

! Since the logfile is not open yet, the error terminations
! must be handled differently than all other cases.
! The routines that normally write to the logfile cannot
! be used just yet. If we cannot open a logfile, we
! always abort execution at this step.

if ( present(progname) ) then
   if (do_output_flag) write(*,*)'Starting program ',trim(progname)
endif

if (do_output_flag) write(*,*)'Initializing the utilities module.'

! Read the namelist entry first before opening logfile, because
! you can rename the logfile via a utilities namelist item.

call find_namelist_in_file("input.nml", "utilities_nml", iunit)
read(iunit, nml = utilities_nml, iostat = io)
call check_namelist_read(iunit, io, "utilities_nml")

! Open the log file with the name from the namelist 
! does not return here on failure.
logfileunit = get_unit()

if (present(alternatename)) then
   lname = alternatename
else
   lname = logfilename
endif

if (do_output_flag) write(*,*)'Trying to log to unit ', logfileunit
if (do_output_flag) write(*,*)'Trying to open file ', trim(lname)

open(logfileunit, file=lname, form='formatted', &
                  action='write', position='append', iostat = io )
if ( io /= 0 ) call fatal_opening_log('initialize_utilities', lname)

! Log the run-time 

if (do_output_flag) then
   if ( present(progname) ) then
      call write_time (logfileunit, label='Starting ', &
                       string1='Program '//trim(progname))
      call write_time (             label='Starting ', &
                       string1='Program '//trim(progname))
   else
      call write_time (logfileunit, label='Starting ')
      call write_time (             label='Starting ')
   endif 
endif

! Check to make sure termlevel is set to a reasonable value
call checkTermLevel(TERMLEVEL)

! Echo the module information using normal mechanism
call register_module(source, revision, revdate)

! Set the defaults for logging the namelist values
call set_nml_output(write_nml)

! If nmlfilename != logfilename, open it.  otherwise set nmlfileunit
! to be same as logunit.
if (do_nml_file()) then
   if (nmlfilename /= lname) then
      if (do_output_flag) &
       write(*,*)'Trying to open namelist log ', trim(nmlfilename)
 
      nmlfileunit = get_unit()

      open(nmlfileunit, file=nmlfilename, form='formatted', &
           position='append', iostat = io )
      if ( io /= 0 ) then
         call error_handler(E_ERR,'initialize_utilities', &
             'Cannot open nm log file', source, revision, revdate)
      endif
 
   else
     nmlfileunit = logfileunit
   endif
endif

! Echo the namelist values for this module using normal mechanism
! including a separator line for this run.
if (do_output_flag) then
   if (do_nml_file() .and. (nmlfileunit /= logfileunit)) then
      if ( present(progname) ) then
         write(nmlfileunit, *) '!Starting Program '//trim(progname)
      else
         write(nmlfileunit, *) '!Starting Program '
      endif 
   endif
   if (do_nml_file()) write(nmlfileunit, nml=utilities_nml)
   if (do_nml_term()) write(     *     , nml=utilities_nml)
endif

! Record the values used for variable types:
if (do_output_flag .and. print_debug) then
  
   call log_it('') ! a little whitespace is nice

   write(string1,*)'..  digits12 is ',digits12
   write(string2,*)'r8       is ',r8
   write(string3,*)'r4       is ',r4
   call error_handler(E_DBG, 'initialize_utilities', string1, &
                      source, revision, revdate, text2=string2, text3=string3)

   write(string1,*)'..  integer  is ',kind(iunit) ! any integer variable will do
   write(string2,*)'i8       is ',i8
   write(string3,*)'i4       is ',i4
   call error_handler(E_DBG, 'initialize_utilities', string1, &
                      source, revision, revdate, text2=string2, text3=string3)
endif

end subroutine initialize_utilities

!-----------------------------------------------------------------------
!>

subroutine finalize_utilities(progname)
character(len=*), intent(in), optional :: progname

! if called multiple times, just return
if (.not. module_initialized) return

if (do_output_flag) then
   if ( present(progname) ) then
      call write_time (logfileunit, label='Finished ', &
                       string1='Program '//trim(progname))
      call write_time (             label='Finished ', &
                       string1='Program '//trim(progname))
      else
         call write_time (logfileunit, label='Finished ')
         call write_time (             label='Finished ')
      endif 

      if (do_nml_file() .and. (nmlfileunit /= logfileunit)) then
         if ( present(progname) ) then
            write(nmlfileunit, *) '!Finished Program '//trim(progname)
         else
            write(nmlfileunit, *) '!Finished Program '
         endif 
      endif 
   endif

   call close_file(logfileunit)
   if ((nmlfileunit /= logfileunit) .and. (nmlfileunit /= -1)) then
      call close_file(nmlfileunit)
   endif

   module_initialized = .false.

end subroutine finalize_utilities

!-----------------------------------------------------------------------
!> log the subversion information about the different source modules
!> being used in this run.

subroutine register_module(src, rev, rdate)
character(len=*), intent(in) :: src, rev, rdate

if ( .not. do_output_flag) return
if ( .not. module_details) return

! you cannot have this routine call init because it calls
! back into register module.  this is an error if this
! routine is called before initialize_utilities().
! AND you cannot use the error handler because it hasn't
! been initialized yet.

if ( .not. module_initialized ) call fatal_not_initialized('register_module')


call log_it('')
call log_it('Registering module :')
call log_it(src)
call log_it(rev)
call log_it(rdate)
call log_it('Registration complete.')
call log_it('')

end subroutine register_module

!-----------------------------------------------------------------------
!> unfortunately you can't pass a namelist as an argument, so all modules
!> with namelists have to call write() themselves.  so this routine and
!> the next are logicals to say whether (and to where) they should write.
!>
!> return whether nml should be written to the nml file

function do_nml_file ()

logical :: do_nml_file

if ( .not. module_initialized ) call initialize_utilities

if ( .not. do_output()) then
   do_nml_file = .false.
else
   do_nml_file = (nml_flag == NML_FILE .or. nml_flag == NML_BOTH)
endif

end function do_nml_file

!-----------------------------------------------------------------------
!> return whether nml should be written to terminal
!> for more details see above.

function do_nml_term ()

logical :: do_nml_term

if ( .not. module_initialized ) call initialize_utilities

if ( .not. do_output()) then
   do_nml_term = .false.
else
   do_nml_term = (nml_flag == NML_TERMINAL .or. nml_flag == NML_BOTH)
endif

end function do_nml_term


!-----------------------------------------------------------------------
!> Opens namelist_file_name if it exists on unit iunit, error if it
!> doesn't exist.
!> Searches file for a line containing ONLY the string  &nml_name, 
!> for instance &filter_nml. If found, backs up one record and
!> returns true. Otherwise, error message and terminates
!>

subroutine find_namelist_in_file(namelist_file_name, nml_name, iunit)

character(len=*),  intent(in)  :: namelist_file_name
character(len=*),  intent(in)  :: nml_name
integer,           intent(out) :: iunit

character(len=256) :: next_nml_string, test_string, string1
integer            :: io

if (.not. module_initialized) call fatal_not_initialized('find_namelist_in_file')

! Decide if there is a logfile or not by looking at the logfileunit.
! if >= 0, ok to write there.

! Check for namelist file existence; no file is an error
if(.not. file_exist(trim(namelist_file_name))) then

   write(msgstring1, *) 'Namelist input file: ', namelist_file_name, ' must exist.'
   if(logfileunit >= 0) then
      call error_handler(E_ERR, 'find_namelist_in_file', msgstring1, &
                         source, revision, revdate)
   else
      call fatal_error_w_no_log('find_namelist_in_file', msgstring1)
   endif

endif

iunit = open_file(namelist_file_name, action = 'read')

! Read each line until end of file is found
! Look for the start of a namelist with &nml_name
! Convert test string to all uppercase ... since that is
! what happens if Fortran writes a namelist.

string1 = adjustl(nml_name)
call to_upper(string1)             ! works in-place
test_string = '&' // trim(string1)

do
   read(iunit, '(A)', iostat = io) next_nml_string
   if(io /= 0) then
      ! Reached end of file and didn't find this namelist
      write(msgstring1, *) 'Namelist entry &', trim(nml_name), &
                           ' must exist in file ', trim(namelist_file_name)
      if(logfileunit >= 0) then
         call error_handler(E_ERR, 'find_namelist_in_file', msgstring1, &
                            source, revision, revdate)
      else
         call fatal_error_w_no_log('find_namelist_in_file', msgstring1)
      endif
   else
      ! see if this line starts the namelist we are asking for
      string1 = adjustl(next_nml_string)
      call to_upper(string1)

      if(string1 == test_string) then
         backspace(iunit)
         return
      endif
   endif
end do

end subroutine find_namelist_in_file


!-----------------------------------------------------------------------

!> Confirms that a namelist read was successful. If it failed
!> produces an error message and stops execution.

subroutine check_namelist_read(iunit, iostat_in, nml_name)

integer,          intent(in) :: iunit, iostat_in
character(len=*), intent(in) :: nml_name

character(len=256) :: nml_string
integer            :: io

! If the namelist read was successful, close the namelist file and we're done.
if(iostat_in == 0) then
   call close_file(iunit)
   return
endif

! If it wasn't successful, print the line on which it failed  
backspace(iunit)
read(iunit, '(A)', iostat = io) nml_string

! A failure in this read means that the namelist started but never terminated
! Result was falling off the end, so backspace followed by read fails
if(io /= 0) then
   write(msgstring1, *) 'Namelist ', trim(nml_name), ' started but never terminated'
else
   ! Didn't fall off end so bad entry in the middle of namelist
   write(msgstring1, *) 'INVALID NAMELIST ENTRY: ', trim(nml_string), ' in namelist ', trim(nml_name)
endif

if(logfileunit >= 0) then
   call error_handler(E_ERR, 'check_namelist_read', msgstring1, &
                      source, revision, revdate)
else
   call fatal_error_w_no_log('check_namelist_read', msgstring1)
endif

end subroutine check_namelist_read

!-----------------------------------------------------------------------
!> if trying to write an unformatted string, like "write(*,*)"
!> to both standard output and the logfile, call this routine instead.
!> it prevents you from having to maintain two copies of the same
!> output message.

subroutine log_it(message)
character(len=*), intent(in) :: message

                      write(     *     , *) trim(message)
if (logfileunit >= 0) write(logfileunit, *) trim(message)

end subroutine log_it


!-----------------------------------------------------------------------
!> call this routine if you cannot open the log file.

subroutine fatal_opening_log(from_routine, lname)

character(len=*), intent(in) :: from_routine
character(len=*), intent(in) :: lname

write(*,*)'FATAL ERROR in '//trim(from_routine)
write(*,*)'   unable to open the logfile for writing.'
write(*,*)'   the logfile name is "',trim(lname),'"'
write(*,*)'  ',trim(source)
write(*,*)'  ',trim(revision)
write(*,*)'  ',trim(revdate)
write(*,*)'   stopping.'
call exit_all(66)

end subroutine fatal_opening_log

!-----------------------------------------------------------------------
!> call this routine if you end up in a function or subroutine here
!> that CANNOT be called before the initialize_utilities() or
!> initialize_mpi_utilities() routine is called.

subroutine fatal_not_initialized(from_routine)

character(len=*), intent(in) :: from_routine

write(*,*)'FATAL ERROR in '//trim(from_routine)
write(*,*)'   initialize_utilities() or initialize_mpi_utilities()'
write(*,*)'   must be called before calling '//trim(from_routine)//'().'
write(*,*)'  ',trim(source)
write(*,*)'  ',trim(revision)
write(*,*)'  ',trim(revdate)
write(*,*)'   stopping.'
call exit_all(77)

end subroutine fatal_not_initialized

!-----------------------------------------------------------------------
!> call this routine if you find an error before the logfile
!> has been opened.

subroutine fatal_error_w_no_log(from_routine, msgstring1, msgstring2)

character(len=*), intent(in) :: from_routine
character(len=*), intent(in), optional :: msgstring1
character(len=*), intent(in), optional :: msgstring2

write(*,*)'FATAL ERROR in '//trim(from_routine)
if (present(msgstring1)) write(*,*) trim(msgstring1)
if (present(msgstring2)) write(*,*) trim(msgstring2)
write(*,*)'  ',trim(source)
write(*,*)'  ',trim(revision)
write(*,*)'  ',trim(revdate)
write(*,*)'   stopping.'
call exit_all(88)

end subroutine fatal_error_w_no_log

!-----------------------------------------------------------------------

subroutine checkTermLevel(level)

integer, intent(in) :: level

select case (level)
  case (E_MSG, E_ALLMSG, E_WARN, E_ERR, E_DBG)
    ! ok, do nothing
  case default
    write(msgstring1, *) 'bad integer value for "termlevel", must be one of'
    write(msgstring2, *) '-1 (E_MSG), 0 (E_ALLMSG), 1 (E_WARN), 2 (E_ERR), -2 (E_DBG)'
    call error_handler(E_ERR,'checkTermLevel', msgstring1, &
                       source, revision, revdate, text2=msgstring2)

  end select

end subroutine checkTermLevel


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! file and error handling
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>

function file_exist (file_name)

character(len=*), intent(in) :: file_name
logical :: file_exist

if ( .not. module_initialized ) call initialize_utilities

inquire (file=file_name, exist=file_exist)

end function file_exist


!-----------------------------------------------------------------------
!>

! get available file unit number
function get_unit() 
integer :: get_unit

integer :: i, iunit
logical :: already_open

if ( .not. module_initialized ) call fatal_not_initialized('get_unit')

iunit = -1
do i = 10, 80
   inquire (i, opened=already_open)
   if (.not. already_open) then
      get_unit = i
      return
   endif
enddo

! if you get here it is an error
write(msgstring1, *) 'Unable to find an available unit number between 10 and 80'
if (logfileunit >= 0) then
   call error_handler(E_ERR,'get_unit', msgstring1, source, revision, revdate)
else
   call fatal_error_w_no_log('get_unit', msgstring1)
endif

end function get_unit


!-----------------------------------------------------------------------
!> write to log and/or standard output, messages, warnings, debug, errors
!>


subroutine error_handler(level, routine, text, src, rev, rdate, aut, text2, text3 )

integer, intent(in) :: level
character(len=*), intent(in) :: routine, text
character(len=*), intent(in), optional :: src, rev, rdate, aut, text2, text3

character(len=16) :: taskstr, msgtype
character(len=246) :: wherefrom, wherecont

if ( .not. module_initialized ) call initialize_utilities

! early returns:

! messages only print if the 'do_output_flag' is on, which by default
! is only task 0.  debug messages only print if enabled.

if (level == E_MSG .and. .not. do_output_flag) return
if (level == E_DBG .and. .not. print_debug)    return

! if we get here, we're printing something.  set up some strings
! to make the code below simpler.

      if ( single_task ) then
   taskstr = ''
      else
        if (task_number == 0) then
      write(taskstr, '(a)' ) "PE 0: "
        else
      write(taskstr, '(a,i5,a)' ) "PE ", task_number, ": "
        endif
      endif

! these are going to get used a lot below.  make them
! single strings so the code is easier to parse.  but can't
! add trailing blanks here because trim will strip them below.

wherefrom = trim(taskstr)//' '//trim(routine)
wherecont = trim(taskstr)//' '//trim(routine)//' ...'

if (level == E_ERR)  msgtype = 'ERROR FROM:'
if (level == E_WARN) msgtype = 'WARNING FROM:'
if (level == E_DBG)  msgtype = 'DEBUG FROM:'

! current choice is to log all errors and warnings regardless
! of setting of output flag.  messages only print from those
! tasks which are allowed to print.

select case(level)
   case (E_MSG)
                          call log_it(trim(wherefrom)//' '//trim(text))
      if (present(text2)) call log_it(trim(wherecont)//' '//trim(text2))
      if (present(text3)) call log_it(trim(wherecont)//' '//trim(text3))

   case (E_ALLMSG)

         if ( single_task ) then
                             call log_it(trim(wherefrom)//' '//trim(text))
         if (present(text2)) call log_it(trim(wherecont)//' '//trim(text2))
         if (present(text3)) call log_it(trim(wherecont)//' '//trim(text3))
           else
        ! this has a problem that multiple tasks are writing to the same logfile.
        ! it's overwriting existing content.  short fix is to NOT write ALLMSGs
        ! to the log file, only stdout.
                            write(*,*) trim(trim(wherefrom)//' '//trim(text))
        if (present(text2)) write(*,*) trim(trim(wherecont)//' '//trim(text2))
        if (present(text3)) write(*,*) trim(trim(wherecont)//' '//trim(text3))
      endif

   case (E_DBG, E_WARN, E_ERR)

      call log_it(msgtype)
      call log_it(wherefrom)
      call log_it(' routine: '//trim(routine))
      call log_it(' message: '//trim(text))
      if (present(text2)) call log_it(' message: ... '//trim(text2))
      if (present(text3)) call log_it(' message: ... '//trim(text3))
      call log_it('')
      call log_it(' source file: '//trim(src))
      call log_it(' file revision: '//trim(rev))
      call log_it(' revision date: '//trim(rdate))
      if(present(aut)) call log_it(' last editor: '//trim(aut))

end select

! TERMLEVEL gets set in the namelist
if( level >= TERMLEVEL ) call exit_all( 99 ) 

end subroutine error_handler


!-----------------------------------------------------------------------
!>

function open_file (fname, form, action, access, convert, delim, reclen, return_rc) result (iunit)

character(len=*), intent(in)            :: fname
character(len=*), intent(in),  optional :: form, action, access, convert, delim
integer,          intent(in),  optional :: reclen
integer,          intent(out), optional :: return_rc
integer  :: iunit

integer           :: rc, rlen
logical           :: open, use_recl
character(len=32) :: format, pos, act, stat, acc, conversion, del

if ( .not. module_initialized ) call initialize_utilities

! if file already open, set iunit and return
inquire (file=trim(fname), opened=open, number=iunit, iostat=rc)
if (open) then
   if (present(return_rc)) return_rc = rc
   return
endif

! not already open, so open it.
      
! set defaults, and then modify depending on what user requests
! via the arguments.  this combination of settings either creates
! a new file or overwrites an existing file from the beginning.

format     = 'formatted'
act        = 'readwrite'
pos        = 'rewind'
stat       = 'unknown'
acc        = 'sequential'
rlen       = 1
del        = 'apostrophe'
conversion = 'native'

if (present(form)) format = form
call to_upper(format)  

! change defaults based on intended action.
if (present(action)) then
    select case(action)

       case ('read', 'READ')
          ! open existing file.  fail if not found.  read from start.
          act  = 'read'
          stat = 'old'

       case ('write', 'WRITE')
          ! create new file/replace existing file.  write at start.
          act  = 'write'
          stat = 'replace'

       case ('append', 'APPEND')
          ! create new/open existing file.  write at end if existing.
          act  = 'readwrite'
          pos  = 'append'

       case default
          ! if the user specifies an action, make sure it is a valid one.
          write(msgstring1,*) 'opening file "'//trim(fname)//'"'
          write(msgstring2,*) 'unrecognized action, "'//trim(action)//'"; valid values: "read", "write", "append"'
          call error_handler(E_ERR, 'open_file', msgstring1, source, revision, revdate, text2=msgstring2)
    end select
endif

! from the ibm help pages:
!   valid values for access: SEQUENTIAL, DIRECT or STREAM. 
!   If ACCESS= is DIRECT, RECL= must be specified. 
!   If ACCESS= is STREAM, RECL= must not be specified.
!   SEQUENTIAL is the default, for which RECL= is optional
! i can't see how to specify all the options in any kind of reasonable way.
! but i need to be able to specify 'stream'... so here's a stab at it.

if (present(access)) then
   acc = access
   call to_upper(acc)
endif

! recl can't apply to stream files, is required for direct,
! and is optional for sequential.  ugh.
if (present(reclen)) then
   rlen = reclen
   use_recl = .true.
else if (acc == 'DIRECT') then
   use_recl = .true.
else
   use_recl = .false.
endif

! endian-conversion only applies to binary files
! valid values seem to be:  'native', 'big-endian', 'little-endian', and possibly 'cray'
! depending on the compiler.
if (present(convert)) then 
   if (format == 'FORMATTED') then
      write(msgstring1,*) 'opening file "'//trim(fname)//'"'
      write(msgstring2,*) 'cannot specify binary conversion on a formatted file'
      call error_handler(E_ERR, 'open_file ', msgstring1, source, revision, revdate, text2=msgstring2)
   endif
   conversion = convert
endif

! string delimiters only apply to ascii files
if (present(delim)) then
   if (format /= 'FORMATTED') then
      write(msgstring1,*) 'opening file "'//trim(fname)//'"'
      write(msgstring2,*) 'cannot specify a delimiter on an unformatted file'
      call error_handler(E_ERR, 'open_file ', msgstring1, source, revision, revdate, text2=msgstring2)
   endif
   del = delim
endif

! ok, now actually open the file

iunit = get_unit()

if (format == 'FORMATTED') then
   ! formatted file: only pass in recl if required
   if (use_recl) then
      open (iunit, file=trim(fname), form=format, access=acc, recl=rlen, &
            delim=del, position=pos, action=act, status=stat, iostat=rc)
   else
      open (iunit, file=trim(fname), form=format, access=acc,            &
            delim=del, position=pos, action=act, status=stat, iostat=rc)
   endif
else  
   ! unformatted file - again, only pass in recl if required 
   if (use_recl) then
      open (iunit, file=trim(fname), form=format, access=acc, recl=rlen, &
            convert=conversion, position=pos, action=act, status=stat, iostat=rc)
   else
      open (iunit, file=trim(fname), form=format, access=acc,            &
            convert=conversion, position=pos, action=act, status=stat, iostat=rc)
   endif
endif
if (rc /= 0 .and. print_debug) call dump_unit_attributes(iunit) 

if (present(return_rc)) then
   return_rc = rc
   return
endif

if (rc /= 0) then
   write(msgstring1, *)'Cannot open file "'//trim(fname)//'" for '//trim(act)
   write(msgstring2,*)'File may not exist or permissions may prevent the requested operation'
   write(msgstring3,*)'Error code was ', rc
   call error_handler(E_ERR, 'open_file: ', msgstring1, source, revision, revdate, &
                      text2=msgstring2, text3=msgstring3)
endif

end function open_file

!-----------------------------------------------------------------------

!> Closes the given unit_number if that unit is open.
!> Not an error to call on an already closed unit.
!> Will print a message if the status of the unit cannot be determined.

subroutine close_file(iunit)

integer, intent(in) :: iunit

integer :: ios
logical :: open

if ( .not. module_initialized ) call initialize_utilities

inquire (unit=iunit, opened=open, iostat=ios)
if ( ios /= 0 ) then
write(msgstring1,*)'Unable to determine status of file unit ', iunit
call error_handler(E_MSG, 'close_file: ', msgstring1, source, revision, revdate)
endif

if (open) close(iunit)

end subroutine close_file

!-----------------------------------------------------------------------
!> Function that returns .true. if this unit number refers to an open file.

function is_file_open(iunit)

integer, intent(in) :: iunit
logical :: is_file_open

integer :: ios
logical :: open

if ( .not. module_initialized ) call initialize_utilities

inquire (unit=iunit, opened=open, iostat=ios)
if ( ios /= 0 ) then
   write(msgstring1,*)'Unable to determine status of file unit ', iunit
   call error_handler(E_MSG, 'is_file_open: ', msgstring1, source, revision, revdate)
endif

is_file_open = open

end function is_file_open

!-----------------------------------------------------------------------
!> Common routine for decoding read/write file format string.
!> Returns .true. for formatted/ascii file, .false. is unformatted/binary
!> Defaults (if fform not specified) to formatted/ascii.

function ascii_file_format(fform)

character(len=*), intent(in), optional :: fform
logical                                :: ascii_file_format


if ( .not. module_initialized ) call initialize_utilities

! Default to formatted/ascii.
if ( .not. present(fform)) then
   ascii_file_format = .true.
   return
endif

SELECT CASE (fform)
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      ascii_file_format = .false.
   CASE DEFAULT
      ascii_file_format = .true.
END SELECT

end function ascii_file_format

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! time and text/text file handling
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>@todo FIXME:  nsc opinion:
!> 1. this routine should NOT support 'end' anymore.  the calling code
!> should call finalize_utilities() directly.  
!> 2. 'brief' format should be the default (easier to grep for in output,
!> or sed for postprocessing)
!> 3. write_time() should be able to take a string to write into, and there should
!> be an easy way to write to both the log and standard out in a single call.

subroutine timestamp(string1,string2,string3,pos)

   character(len=*), optional, intent(in) :: string1
   character(len=*), optional, intent(in) :: string2
   character(len=*), optional, intent(in) :: string3
   character(len=*),           intent(in) :: pos

   if ( .not. module_initialized ) call initialize_utilities
   if ( .not. do_output_flag) return

!>@todo remove this option
   if (trim(adjustl(pos)) == 'end') then
      call finalize_utilities()
   else if (trim(adjustl(pos)) == 'brief') then
      call write_time (logfileunit, brief=.true., & 
                       string1=string1, string2=string2, string3=string3)
      call write_time (             brief=.true., &
                       string1=string1, string2=string2, string3=string3)
       
   else
      call write_time (logfileunit, & 
                       string1=string1, string2=string2, string3=string3)
      call write_time (string1=string1, string2=string2, string3=string3)
       
   endif

end subroutine timestamp

!-----------------------------------------------------------------------
!> write time to the given unit.  should have option to write to a string
!> or both log and stdout.  and brief should be the default.  (my opinion. nsc)
!>
!>   in: unit number (default is * if not specified)
!>   in: label (default is  "Time is" if not specified)
!>   in: string1,2,3 (no defaults)
!>
!>  default output is a block of 3-4 lines, with dashed line separators
!>  and up to 3 descriptive text strings.
!>  if brief specified as true, only string1 printed if given,
!>  and time printed on same line in YYYY/MM/DD HH:MM:SS format
!>  with the tag 'TIME:' before it.  should be easier to postprocess.

subroutine write_time (unit, label, string1, string2, string3, tz, brief)

integer,          optional, intent(in) :: unit
character(len=*), optional, intent(in) :: label
character(len=*), optional, intent(in) :: string1
character(len=*), optional, intent(in) :: string2
character(len=*), optional, intent(in) :: string3
logical,          optional, intent(in) :: tz
logical,          optional, intent(in) :: brief


integer :: lunit
character(len= 8) :: cdate
character(len=10) :: ctime
character(len= 5) :: zone
integer, dimension(8) :: values
logical :: oneline

if (present(unit)) then
   lunit = unit
else
   lunit = 6   ! this should be *
endif

call DATE_AND_TIME(cdate, ctime, zone, values)

! give up if no good values were returned
if (.not. any(values /= -HUGE(0)) ) return 

oneline = .false.
if (present(brief)) oneline = brief

!>@todo write into a string to avoid replicating complex lines
!> add on the label if it's there separately.

if (oneline) then
   if (present(string1)) then
      write(lunit,'(A,1X,I4,5(A1,I2.2))') string1//' TIME:', &
                     values(1), '/', values(2), '/', values(3), &
                     ' ', values(5), ':', values(6), ':', values(7)
   else
      write(lunit,'(A,1X,I4,5(A1,I2.2))') 'TIME: ', &
                     values(1), '/', values(2), '/', values(3), &
                     ' ', values(5), ':', values(6), ':', values(7)
   endif
else
   write(lunit,*)
   write(lunit,*)'--------------------------------------'
   if ( present(label) ) then
      write(lunit,*) label // '... at YYYY MM DD HH MM SS = '
   else
      write(lunit,*) 'Time is  ... at YYYY MM DD HH MM SS = '
   endif 
   write(lunit,'(17x,i4,5(1x,i2))') values(1), values(2), &
                     values(3),  values(5), values(6), values(7)

   if(present(string1)) write(lunit,*)trim(string1)
   if(present(string2)) write(lunit,*)trim(string2)
   if(present(string3)) write(lunit,*)trim(string3)

   if (present(tz)) then
      if ( values(4) /= -HUGE(0) .and. tz) &
         write(lunit,*)'time zone offset is ',values(4),' minutes.'
   endif

   write(lunit,*)'--------------------------------------'
   write(lunit,*)
endif

end subroutine write_time

!-----------------------------------------------------------------------
!> set whether output is written to a log file or simply ignored 
!>
!>   in:  doflag  = whether to output log information or not


subroutine set_output (doflag)

logical, intent(in) :: doflag

!! THIS ONE IS DIFFERENT.  Set the flag FIRST before doing the
!! standard initialization, so if you are turning off writing
!! for some tasks you do not get output you are trying to avoid.

do_output_flag = doflag

if ( .not. module_initialized ) call initialize_utilities

end subroutine set_output

!-----------------------------------------------------------------------
!> return whether output should be written from this task 

function do_output ()

logical :: do_output

if ( .not. module_initialized ) call initialize_utilities

do_output = do_output_flag

end function do_output

!-----------------------------------------------------------------------
!> set whether nml output is written to stdout file or only nml file
!>
!>    in:  doflag  = whether to output nml information to stdout 

subroutine set_nml_output (nmlstring)

character(len=*), intent(in) :: nmlstring

! initialize_utilities calls this routine, so you cannot call
! the init routine from here.

if ( .not. module_initialized ) call fatal_not_initialized('set_nml_output')

select case (nmlstring)
   case ('NONE', 'none')
      nml_flag = NML_NONE
      call error_handler(E_MSG, 'set_nml_output', &
                         'No echo of NML values')

   case ('FILE', 'file')
      nml_flag = NML_FILE
      call error_handler(E_MSG, 'set_nml_output', &
                         'Echo NML values to log file only')
  
   case ('TERMINAL', 'terminal')
      nml_flag = NML_TERMINAL
      call error_handler(E_MSG, 'set_nml_output', &
                         'Echo NML values to terminal output only')

   case ('BOTH', 'both')
      nml_flag = NML_BOTH
      call error_handler(E_MSG, 'set_nml_output', &
                         'Echo NML values to both log file and terminal')

   case default
      call error_handler(E_ERR, 'set_nml_output', &
                        'unrecognized input string: '//trim(nmlstring), &
                        source, revision, revdate)
 
end select

end subroutine set_nml_output

!-----------------------------------------------------------------------
!> for multiple-task jobs, set the task number for error msgs 
!>
!>    in:  tasknum  = task number, 0 to N-1

subroutine set_tasknum (tasknum)

integer, intent(in) :: tasknum

if ( .not. module_initialized ) call initialize_utilities

single_task = .false. 
task_number = tasknum

end subroutine set_tasknum

!-----------------------------------------------------------------------

!>@todo FIXME: move this to the netcdf utilities module

subroutine nc_check(istatus, subr_name, context)
   integer, intent (in)                   :: istatus
   character(len=*), intent(in)           :: subr_name
   character(len=*), intent(in), optional :: context
  
   character(len=512) :: error_msg
  
   ! if no error, nothing to do here.  we are done.
   if( istatus == nf90_noerr) return


   ! something wrong.  construct an error string and call the handler.

   ! context is optional, but is very useful if specified.
   ! if context + error code > 512, the assignment will truncate.
   if (present(context) ) then
       error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))
   else
       error_msg = nf90_strerror(istatus)
   endif

   ! this does not return 
   call error_handler(E_ERR, 'nc_check', error_msg, source, revision, revdate, &
                      text2=subr_name)
  

end subroutine nc_check


!-----------------------------------------------------------------------
!> convert a string to upper case *in place*

subroutine to_upper( string )

character(len=*), intent(inout) :: string

integer :: ismalla, ibiga, i

ismalla = ichar('a')
ibiga   = ichar('A')

do i = 1,len(string)
   if ((string(i:i) >= 'a') .and. (string(i:i) <= 'z')) then
        string(i:i)  = char( ichar(string(i:i)) + ibiga - ismalla)
   endif
enddo

end subroutine to_upper


!-----------------------------------------------------------------------
!> copy instring to outstring, omitting all internal blanks
!> outstring must be at least as long as instring

subroutine squeeze_out_blanks(instring, outstring)
character(len=*), intent(in)  :: instring
character(len=*), intent(out) :: outstring

integer :: i, o

outstring = ''

o = 1
do i = 1,len_trim(instring)
   if (instring(i:i) == ' ') cycle
   outstring(o:o) = instring(i:i)
   o = o + 1
enddo

end subroutine squeeze_out_blanks

!-----------------------------------------------------------------------
!> Determines the number of lines and maximum line length
!> of an ascii file.

subroutine find_textfile_dims( fname, nlines, linelen )

character(len=*),  intent(in)  :: fname
integer,           intent(out) :: nlines
integer, optional, intent(out) :: linelen

integer :: i, maxlen, mylen, ios, funit

character(len=1024) :: oneline
character(len=512)  :: error_msg

! if there is no file, return -1 for both counts
if (.not. file_exist(fname)) then
  nlines = -1
  if (present(linelen)) linelen = -1
  return
endif

! the file exists, go count things up.
nlines  = 0
maxlen  = 0
funit   = open_file(fname, form="FORMATTED", action="READ")

READLOOP : do i = 1,100000

   read(funit, '(A)', iostat=ios) oneline
   if (ios < 0) exit READLOOP  ! end of file
   if (ios > 0) then
      write(error_msg,'(A,'' read around line '',i8)')trim(fname),nlines
      call error_handler(E_ERR,'find_textfile_dims', error_msg, &
                         source, revision, revdate)
   endif

   nlines = nlines + 1
   mylen  = len_trim(oneline)

   if (mylen > maxlen) maxlen = mylen

enddo READLOOP

call close_file(funit)

if (present(linelen)) linelen = maxlen

end subroutine find_textfile_dims


!-----------------------------------------------------------------------
!> Reads a text file into a character variable.
!> Initially needed to read a namelist file into a variable that could 
!> then be inserted into a netCDF file. Due to a quirk in the way Fortran
!> and netCDF play together, I have not figured out how to dynamically
!> create the minimal character length ... so any line longer than
!> the declared length of the textblock variable is truncated.

subroutine file_to_text( fname, textblock )

character(len=*),               intent(in)  :: fname
character(len=*), dimension(:), intent(out) :: textblock

integer :: i, ios, funit
integer :: mynlines, mylinelen, strlen

character(len=512)  :: string

call find_textfile_dims(fname, mynlines, mylinelen)

strlen = len(textblock)

if ( ( mynlines /= size(textblock) ) .or. &
     (mylinelen >      strlen    ) ) then
   write(string,'(A, '' file shape is '',i6,'' by '',i4, &
                    &'' truncating to '',i6,'' by '',i4)') &
   trim(fname),mynlines,mylinelen,size(textblock),strlen
   call error_handler(E_MSG,'file_to_text', trim(string))
endif

funit   = open_file(fname, form="FORMATTED", action="READ")

strlen  = min(mylinelen, strlen)

do i = 1,mynlines

   read(funit, '(A)', iostat=ios) string

   write(textblock(i),'(A)') string(1:strlen)

   if ( ios /= 0 ) then
      write(string,'(A,'' read around line '',i8)')trim(fname),i
      call error_handler(E_ERR,'file_to_text', trim(string), &
                         source, revision, revdate)
   endif

enddo 

call close_file(funit)

end subroutine file_to_text

!-----------------------------------------------------------------------
!> Arguments are the name of a file which contains a list of filenames.
!> This routine opens the listfile, and returns the index-th one.
!>@todo FIXME this should return 512 chars, not 256.  will that 
!>break some of the calling code?

function get_next_filename( listname, lineindex )

character(len=*),  intent(in) :: listname
integer,           intent(in) :: lineindex
character(len=256)            :: get_next_filename

integer :: i, ios, funit

character(len=512)  :: string

funit   = open_file(listname, form="FORMATTED", action="READ")

do i=1, lineindex

   read(funit, '(A)', iostat=ios) string

   ! reached end of file, return '' as indicator.
   if ( ios /= 0 ) then
      get_next_filename = ''
      call close_file(funit)
      return
   endif

enddo

! check for length problems
if (len_trim(string) > len(get_next_filename)) then
   call error_handler(E_ERR, 'get_next_filename', &
                      'maximum filename length of 256 exceeded', &
                      source, revision, revdate)   
endif

! is this overly complicated?  either just adjustl() or nothing?
get_next_filename = adjustl(string(1:len(get_next_filename)))
call close_file(funit)

end function get_next_filename


!-----------------------------------------------------------------------
!> this function is intended to be used when there are 2 ways to specify 
!> an unknown number of input files, most likely in a namelist.
!>
!> e.g. to specify 3 input files named 'file1', 'file2', 'file3' you can
!> either give:
!>   input_files = 'file1', 'file2', 'file3' 
!>     OR
!>   input_file_list = 'flist'  
!>
!>   and the contents of text file 'flist' are (e.g. 'cat flist' gives)
!>          file1
!>          file2
!>          file3
!>   each filename is on a separate line with no quotes
!>
!> you pass both those variables into this function and when it returns
!> the 'name_array' will have all the filenames as an array of strings.
!> it will have opened the text files, if given, and read in the contents.
!> it returns the number of filenames it found.
!>
!> contrast this with the following function where you tell it how many
!> names and/or indirect files you expect, and how many filenames you
!> expect to have at the end.

function set_filename_list(name_array, listname, caller_name)

! return the count of names specified by either the name_array()
! or the inside the listname but not both.  caller_name is used
! for error messages.  verify that if a listname is used that
! it does not contain more than the allowed number of input names
! (specified by the length of the name_array).  the listname,
! if specified, must be the name of an ascii input file with
! a list of names, one per line.

character(len=*), intent(inout) :: name_array(:)
character(len=*), intent(in)    :: listname
character(len=*), intent(in)    :: caller_name
integer                         :: set_filename_list

integer :: fileindex, max_num_input_files
logical :: from_file
character(len=64) :: fsource

! here's the logic:
! if the user specifies neither name_array nor listname, error
! if the user specifies both, error.
! if the user gives a filelist, we make sure the length is not more
!   than maxfiles and read it into the explicit list and continue.
! when this routine returns, the function return val is the count
! and the names are in name_array()

if (name_array(1) == '' .and. listname == '') then
   call error_handler(E_ERR, caller_name, &
          'must specify either filenames in the namelist, or a filename containing a list of names', &
          source,revision,revdate)
endif
   
! make sure the namelist specifies one or the other but not both
if (name_array(1) /= '' .and. listname /= '') then
   call error_handler(E_ERR, caller_name, &
       'cannot specify both filenames in the namelist and a filename containing a list of names', &
       source,revision,revdate)
endif

! if they have specified a file which contains a list, read it into
! the name_array array and set the count.
if (listname /= '') then
   fsource = 'filenames contained in a list file'
   from_file = .true.
else
   fsource = 'filenames in the namelist'
   from_file = .false.
endif

! the max number of names allowed in a list file is the 
! size of the name_array passed in by the user.
max_num_input_files = size(name_array)

! loop over the inputs.  if the names were already specified in the
! name_array, just look for the '' to indicate the end of the list.
! if the names were specified in the listname file, read them in and
! fill in the name_array and then look for ''.
do fileindex = 1, max_num_input_files
   if (from_file) &
      name_array(fileindex) = get_next_filename(listname, fileindex)

   if (name_array(fileindex) == '') then
      if (fileindex == 1) then
         write(msgstring2,*)'reading file # ',fileindex
         write(msgstring3,*)'reading file name "'//trim(name_array(fileindex))//'"'
         call error_handler(E_ERR, caller_name, 'found no '//trim(fsource), &
                            source,revision,revdate,text2=msgstring2,text3=msgstring3)
      endif

      ! at the end of the list. return how many filenames were found, 
      ! whether the source was the name_array or the listname.
      set_filename_list = fileindex - 1
      return
   endif
enddo

! if you get here, you read in all max_num_input_files without
! seeing an empty string.  if the input names were already in the
! array, you're done - set the count and return.   but if you're
! reading names from a file it is possible to specify more names
! than fit in the list.  test for that and give an error if you
! aren't at the end of the list.

if (from_file) then
   if (get_next_filename(listname, max_num_input_files+1) /= '') then
      write(msgstring1, *) 'cannot specify more than ',max_num_input_files,' filenames in the list file'
      call error_handler(E_ERR, caller_name, msgstring1, source,revision,revdate)
   endif
endif

set_filename_list = max_num_input_files

end function set_filename_list

!-----------------------------------------------------------------------
!> this function is intended to be used when there are 2 ways to specify 
!> a KNOWN number of input files, most likely in a namelist.
!>
!> e.g. to specify 6 input files named 'file1a', 'file2a', 'file3a', 
!>                                     'file1b', 'file2b', 'file3b', 
!> either give:
!>   input_files = 'file1a', 'file2a', 'file3a', 'file1b', 'file2b', 'file3b', 
!>     OR
!>   input_file_list = 'flistA', 'flistB'
!>
!>   and the contents of text file 'flistA' are (e.g. 'cat flistA' gives)
!>          file1a
!>          file2a
!>          file3a
!>   and the contents of text file 'flistB' are (e.g. 'cat flistB' gives)
!>          file1b
!>          file2b
!>          file3b
!>
!>   each filename is on a separate line with no quotes, and the result
!>   is all of the contents of list1 followed by list2, or the files in
!>   order as they're given explicitly in the first string array.
!>
!> you pass both those variables into this function and when it returns
!> the 'name_array' will have all the filenames as an array of strings.
!> it will have opened the text files, if given, and read in the contents.
!> it dies with a fatal error if it can't construct a list of the right length.
!> (the 'right length' is nlists * nentries)
!>
!> contrast this with the previous function where you don't know (or care)
!> how many filenames are specified, and there's only a single listlist option.

subroutine set_multiple_filename_lists(name_array, listname, nlists, nentries, &
                                       caller_name, origin, origin_list)

! when this routine returns, name_array() contains (nlists * nentries) of names,
! either because they started out there or because we've opened up 'nlists'
! listname files and read in 'nentries' from each.  it's a fatal error not to
! have enough files. caller_name is used for error messages.

character(len=*), intent(inout) :: name_array(:)
character(len=*), intent(in)    :: listname(:)
integer,          intent(in)    :: nlists
integer,          intent(in)    :: nentries
character(len=*), intent(in)    :: caller_name
character(len=*), intent(in)    :: origin
character(len=*), intent(in)    :: origin_list

integer :: fileindex, max_num_input_files
logical :: from_file
character(len=64) :: fsource
integer :: nl, ne, num_lists

! here's the logic:
! if the user specifies neither name_array nor listname, error
! if the user specifies both, error.
! if the user gives a listname, make sure there are nlists of them
!   and each contains at least nentries
! when this routine returns the names are in name_array()

if (name_array(1) == '' .and. listname(1) == '') then
   call error_handler(E_ERR, caller_name, &
          'missing filenames',source,revision,revdate, &
          text2='must specify either "'//trim(origin)//'" in the namelist,', &
          text3='or a "'//trim(origin_list)//'" file containing a list of names')
endif
   
! make sure the namelist specifies one or the other but not both
if (name_array(1) /= '' .and. listname(1) /= '') then
   call error_handler(E_ERR, caller_name, &
          'can not specify both an array of files and a list of files', &
          source,revision,revdate, &
          text2='must specify either "'//trim(origin)//'" in the namelist,', &
          text3='or a "'//trim(origin_list)//'" file containing a list of names')
endif

! if they have specified a file which contains a list, read it into
! the name_array array and set the count.
if (listname(1) /= '') then
   fsource = ' contained in a list file'
   from_file = .true.
else
   fsource = ' in the namelist'
   from_file = .false.
endif

! this version of the code knows how many files should be in the listname
if (from_file) then
   num_lists = size(listname)
   if (num_lists < nlists) then
      write(msgstring1, *) 'expecting ', nlists, ' filenames in "'//trim(origin_list)//'", got ', num_lists
      call error_handler(E_ERR, caller_name, msgstring1, source,revision,revdate)
   endif
endif
   
! the max number of names allowed in a list file is the 
! size of the name_array passed in by the user.
max_num_input_files = size(name_array)
if (max_num_input_files < nlists * nentries) then
   write(msgstring1, *) 'list length = ', max_num_input_files, '  needs room for ', nlists * nentries
   call error_handler(E_ERR, caller_name, 'internal error: name_array not long enough to hold lists', &
                      source,revision,revdate, text2=msgstring1)
endif

! loop over the inputs.  if the names were already specified in the
! name_array leave them there.  if they were in the listname file,
! read them into the name_array.

do nl = 1, nlists
   do ne = 1, nentries
      fileindex = (nl-1) * nentries + ne

      if (from_file) &
         name_array(fileindex) = get_next_filename(listname(nl), ne)
   
      if (name_array(fileindex) == '') then
         write(msgstring1, *) 'Missing filename'

         if (from_file) then
            write(msgstring2,*)'reading entry # ', nl, ' from "'//trim(origin_list)//'"'
            write(msgstring3,*)'expecting ', nentries, ' files, have ', ne-1
         else
            write(msgstring2,*)'required entry # ', fileindex, ' from "'//trim(origin)//'"'
            write(msgstring3,*)'expecting ', nlists*nentries, ' filenames, have ', fileindex-1
         endif

         call error_handler(E_ERR, caller_name, trim(msgstring1)//trim(fsource), &
                            source,revision,revdate,text2=msgstring2,text3=msgstring3)
   
      endif
   enddo
enddo

end subroutine set_multiple_filename_lists

!-----------------------------------------------------------------------

function next_file(fname,ifile)

! FIXME: THIS FUNCTION IS DEPRECATED AND SHOULD BE REMOVED.
! FIXME: THIS FUNCTION IS DEPRECATED AND SHOULD BE REMOVED.
! FIXME: THIS FUNCTION IS DEPRECATED AND SHOULD BE REMOVED.

!----------------------------------------------------------------------
! The file name can take one of three forms:
! /absolute/path/to/nirvana/obs_001/obs_seq.final   (absolute path)
! obs_0001/obs_seq.final    (relative path)
! obs_seq.final      (no path ... local)
!
! If there is a '/' in the file name, we grab the portion before the
! slash and look for an underscore. Anything following the underscore
! is presumed to be the portion to increment.
!
! If there is no slash AND ifile is > 1 ... we have already read
! the 'one and only' obs_seq.final file and we return 'done'.
!----------------------------------------------------------------------

character(len=*), intent(in) :: fname
integer,          intent(in) :: ifile

character(len=len(fname)) :: next_file
character(len=len(fname)) :: dir_name

integer,            SAVE :: filenum = 0
integer,            SAVE :: dir_prec = 0
character(len=256), SAVE :: dir_base
character(len=256), SAVE :: filename
character(len=129), SAVE :: dir_ext

integer :: slashindex, splitindex, i, strlen, ios

character(len=512) :: string1, string2, string3

if (len(fname) > len(dir_base) ) then
   write(string1,*)'input filename not guaranteed to fit in local variables'
   write(string2,'('' input filename (len='',i3,'') tempvars are (len='',i3,'')'')') &
   len(fname),len(dir_base)
   write(string3,*)'increase len of dir_base, filename, dir_ext and recompile'
   call error_handler(E_MSG, 'next_file', string1, source, revision, revdate, &
              text2=string2, text3=string3)
endif

if (ifile == 1) then ! First time through ... find things.

   ! Start looking (right-to-left) for the 'slash'.
   ! Anything to the right of it must be a filename.
   ! Anything to the left must be the part that gets incremented.

   filename   = adjustl(fname)
   next_file  = trim(filename)
   strlen     = len_trim(filename)
   slashindex = 0

   SlashLoop : do i = strlen,1,-1
   if ( next_file(i:i) == '/' ) then
      slashindex = i
      exit SlashLoop
   endif
   enddo SlashLoop

   if (slashindex > 0) then ! we have a directory structure

      dir_name   = trim(fname(1:slashindex-1))
      filename   = trim(fname(slashindex+1:129))
      strlen     = len_trim(dir_name)
      splitindex = 0

      SplitLoop : do i = strlen,1,-1
      if ( dir_name(i:i) == '_' ) then
         splitindex = i
         exit SplitLoop
      elseif (dir_name(i:i) == '/') then
         ! there is no underscore in the directory node
         ! immediately preceeding the filename
         exit SplitLoop
      endif
      enddo SplitLoop

      if (splitindex <= 0) then
         filenum  = -1 ! indicates no next file
      else
         dir_base   = dir_name(1:splitindex-1)
         dir_ext    = dir_name(splitindex+1:slashindex-1)
         dir_prec   = slashindex - splitindex - 1
    
         read(dir_ext,*,iostat=ios) filenum
         if(ios /= 0) then
            ! Directory has an '_' separating two alphabetic parts
            ! Nothing to increment.
            filenum = -1
         endif
      endif

   else ! we have one single file - on the first trip through

      filenum  = -1 ! indicates no next file

   endif 

else

   if (filenum < 0) then
      next_file = 'doneDONEdoneDONE'
   else

      filenum = filenum + 1
      if (dir_prec == 1) then
      write(next_file,'(a,''_'',i1.1,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 2) then
      write(next_file,'(a,''_'',i2.2,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 3) then
      write(next_file,'(a,''_'',i3.3,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 4) then
      write(next_file,'(a,''_'',i4.4,''/'',a)') trim(dir_base),filenum,trim(filename)
      else
      write(next_file,'(a,''_'',i5.5,''/'',a)') trim(dir_base),filenum,trim(filename)
      endif

      write(string1,*)'WARNING: This feature is deprecated and will be removed in the next release.'
      write(string2,*)'to use multiple input files, use the "list" construct.'
      call error_handler(E_MSG,'next_file',string1,source,revision,revdate,text2=string2)

   endif

endif

end function next_file

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! generic routines needed by more than one part of the code
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!>  uniform way to treat longitude ranges, in degrees, on a globe.
!>  returns true if lon is between min and max, starting at min
!>  and going EAST until reaching max.  wraps across 0 longitude.
!>  if min == max, all points are inside.  includes edges.
!>  if optional arg doradians is true, do computation in radians 
!>  between 0 and 2*PI instead of 360.   if given, return the
!>  'lon' value possibly + 360 (or 2PI) which can be used for averaging
!>  or computing on a consistent set of longitude values.  after the
!>  computation is done if the answer is > 360 (or 2PI), subtract that
!>  value to get back into the 0 to 360 (or 2PI) range.

function is_longitude_between (lon, minlon, maxlon, doradians, newlon)

real(r8), intent(in)            :: lon, minlon, maxlon
logical,  intent(in),  optional :: doradians
real(r8), intent(out), optional :: newlon
logical :: is_longitude_between

real(r8) :: minl, maxl, lon2, circumf

circumf = 360.0_r8
if (present(doradians)) then
  if (doradians) circumf = TWOPI
endif

! ensure the valid region boundaries are between 0 and one circumference
! (must use modulo() and not mod() so negative vals are handled ok)
minl = modulo(minlon, circumf)
maxl = modulo(maxlon, circumf)

! boundary points are included in the valid region so if min=max 
! the 'region' is the entire globe and you can return early.
if (minl == maxl) then
   is_longitude_between = .true. 
   if (present(newlon)) newlon = lon
   return
endif

! ensure the test point is between 0 and one circumference
lon2  = modulo(lon, circumf)

! here's where the magic happens:
! minl will be bigger than maxl if the region of interest crosses the prime 
! meridian (longitude = 0).  in this case add one circumference to the 
! eastern boundary so maxl is guarenteed to be larger than minl (and valid 
! values are now between 0 and 2 circumferences).  
!
! if the test point longitude is west of the minl boundary add one circumference
! to it as well before testing against the bounds.  values that were east of 
! longitude 0 but west of maxl will now be shifted so they are again correctly 
! within the new range; values that were west of the prime meridian but east 
! of minl will stay in range; values west of minl and east of maxl will be 
! correctly shifted out of range.

if (minl > maxl) then
   maxl = maxl + circumf
   if (lon2 < minl) lon2 = lon2 + circumf
endif

is_longitude_between = ((lon2 >= minl) .and. (lon2 <= maxl))

! if requested, return the value that was tested against the bounds, which 
! will always be between 0 and 2 circumferences and monotonically increasing
! from minl to maxl.  if the region of interest doesn't cross longitude 0
! this value will be the same as the input value.  if the region does
! cross longitude 0 this value will be between 0 and 2 circumferences.
! it's appropriate for averaging values together or comparing them against
! other values returned from this routine with a simple greater than or less
! than without further computation for longitude 0.  to convert the values
! back into the range from 0 to one circumference, compare it to the
! circumference and if larger, subtract one circumference from the value.

if (present(newlon)) newlon = lon2

end function is_longitude_between 


!-----------------------------------------------------------------------
!>

pure function to_scalar_real(x)
 real(r8), intent(in) :: x(1)
 real(r8) :: to_scalar_real
 
to_scalar_real = x(1)

end function to_scalar_real

!-----------------------------------------------------------------------
!>

pure function to_scalar_int(x)
 integer, intent(in) :: x(1)
 integer :: to_scalar_int
 
to_scalar_int = x(1)

end function to_scalar_int

!-----------------------------------------------------------------------
!>

pure function to_scalar_int8(x)
 integer(i8), intent(in) :: x(1)
 integer(i8) :: to_scalar_int8
 
to_scalar_int8 = x(1)

end function to_scalar_int8

!-----------------------------------------------------------------------
!>

function string_to_real(inputstring)

character(len=*), intent(in) :: inputstring
real(r8)                     :: string_to_real

integer :: io

! if the string converts - great, if not, it's MISSING_R8
read(inputstring,*,iostat=io)string_to_real
if (io /= 0) string_to_real = MISSING_R8

end function string_to_real

!-----------------------------------------------------------------------
!>

function string_to_integer(inputstring)

character(len=*), intent(in) :: inputstring
integer                      :: string_to_integer

integer :: io

! if the string converts - great, if not, it's MISSING_I
read(inputstring,*,iostat=io)string_to_integer
if (io /= 0) string_to_integer = MISSING_I

end function string_to_integer


!-----------------------------------------------------------------------
!>  if matching true or false, uppercase the string so any
!>  combination of upper and lower case matches.  match with
!>  and without the enclosing periods (e.g. .true. and true
!>  both match, as would .TrUe.).  also allow plain T and F.
!>  if 'string to match' is provided, no upper casing is done
!>  and the string must match exactly.

function string_to_logical(inputstring, string_to_match)

character(len=*), intent(in)           :: inputstring
character(len=*), intent(in), optional :: string_to_match
logical                                :: string_to_logical

! to_upper() works in place, so if looking for true/false
! make a copy first so the input string can stay intent(in).

character(len=len_trim(inputstring)) :: ucase_instring

! see if the input matches the string given by the caller
if (present(string_to_match)) then
   string_to_logical = (inputstring == string_to_match)
   return
endif

! see if the input matches true or false
ucase_instring = trim(inputstring)
call to_upper(ucase_instring)

select case (ucase_instring)
   case ("TRUE", ".TRUE.", "T")
      string_to_logical = .true.
   case ("FALSE", ".FALSE.", "F")
      string_to_logical = .false.
   case default 
      ! we can't give any context here for where it was
      ! being called, but if it isn't true or false, error out.
      msgstring1 = '.TRUE., TRUE, T or .FALSE., FALSE, F are valid values'         
      call error_handler(E_ERR,'string_to_logical', &
                 'Cannot parse true or false value from string: "'//trim(inputstring)//'"', &
                  source, revision, revdate, text2=msgstring1)
end select

end function string_to_logical

!-----------------------------------------------------------------------
!> dump the contents of a 1d array with a max of N items per line.
!> optional arguments allow the caller to restrict the output to 
!> no more than X items, to write to an open file unit, and to
!> write a text label before the numerical dump.
!>

subroutine array_1d_dump(array, nper_line, max_items, funit, label)
real(r8),         intent(in)           :: array(:)
integer,          intent(in), optional :: nper_line
integer,          intent(in), optional :: max_items
integer,          intent(in), optional :: funit
character(len=*), intent(in), optional :: label

integer :: i, per_line, ounit, asize_i
logical :: has_label

! set defaults and override if arguments are present

per_line = 4
if (present(nper_line)) per_line = nper_line

asize_i = size(array)
if (present(max_items)) asize_i = min(asize_i, max_items)

ounit = 0
if (present(funit)) ounit = funit

has_label = .false.
if (present(label)) has_label = .true.

! output section

if (has_label) write(ounit, *) trim(label)

do i=1, asize_i, per_line
   write(ounit, *) i, ' : ', array(i:min(asize_i,i+per_line-1))
enddo

end subroutine array_1d_dump

!-----------------------------------------------------------------------
!> dump the contents of a 2d array with a max of N items per line.
!> optional arguments allow the caller to restrict the output to 
!> no more than X items, to write to an open file unit, and to
!> write a text label before the numerical dump.
!>

subroutine array_2d_dump(array, nper_line, max_i_items, max_j_items, funit, label)
real(r8),         intent(in)           :: array(:,:)
integer,          intent(in), optional :: nper_line
integer,          intent(in), optional :: max_i_items
integer,          intent(in), optional :: max_j_items
integer,          intent(in), optional :: funit
character(len=*), intent(in), optional :: label

integer :: i, j, per_line, ounit, asize_i, asize_j
logical :: has_label

! set defaults and override if arguments are present

per_line = 4
if (present(nper_line)) per_line = nper_line

asize_i = size(array, 1)
asize_j = size(array, 2)
if (present(max_i_items)) asize_i = min(asize_i, max_i_items)
if (present(max_j_items)) asize_j = min(asize_j, max_j_items)

ounit = 0
if (present(funit)) ounit = funit

has_label = .false.
if (present(label)) has_label = .true.

! output section

if (has_label) write(ounit, *) trim(label)

do j=1, asize_j
   do i=1, asize_i, per_line
      write(ounit, *) i, j, ' : ', array(i:min(asize_i,i+per_line-1), j)
   enddo
enddo

end subroutine array_2d_dump

!-----------------------------------------------------------------------
!> dump the contents of a 3d array with a max of N items per line.
!> optional arguments allow the caller to restrict the output to 
!> no more than X items, to write to an open file unit, and to
!> write a text label before the numerical dump.
!>

subroutine array_3d_dump(array, nper_line, max_i_items, max_j_items, max_k_items, funit, label)
real(r8),         intent(in)           :: array(:,:,:)
integer,          intent(in), optional :: nper_line
integer,          intent(in), optional :: max_i_items
integer,          intent(in), optional :: max_j_items
integer,          intent(in), optional :: max_k_items
integer,          intent(in), optional :: funit
character(len=*), intent(in), optional :: label

integer :: i, j, k, per_line, ounit, asize_i, asize_j, asize_k
logical :: has_label

! set defaults and override if arguments are present

per_line = 4
if (present(nper_line)) per_line = nper_line

asize_i = size(array, 1)
asize_j = size(array, 2)
asize_k = size(array, 3)
if (present(max_i_items)) asize_i = min(asize_i, max_i_items)
if (present(max_j_items)) asize_j = min(asize_j, max_j_items)
if (present(max_k_items)) asize_k = min(asize_k, max_k_items)

ounit = 0
if (present(funit)) ounit = funit

has_label = .false.
if (present(label)) has_label = .true.

! output section

if (has_label) write(ounit, *) trim(label)

do k=1, asize_k
   do j=1, asize_j
      do i=1, asize_i, per_line
         write(ounit, *) i, j, k, ' : ', array(i:min(asize_i,i+per_line-1), j, k)
      enddo
   enddo
enddo

end subroutine array_3d_dump

!-----------------------------------------------------------------------
!> given an array of sorted values and a value to find, return the
!> two indices that enclose that value, and the fraction between.
!>
!> fraction_across = 0.0 is the 100% the smaller index value, 
!>                   1.0 is the 100% the larger index value.
!>
!> if the array values are inverted (e.g. index 1 is the largest value),
!> set inverted = .true.  the interpretation in the calling code for
!> smaller index, larger index and fraction_across remain the same as the default case.
!>
!> if the fraction_across the enclosing level should be computed using a
!> log scale, set log_scale = .true.
!>
!> my_status values:
!>   0 = good return
!>  -1 = value_to_find is below smallest value
!>   1 = value_to_find is above largest value
!>  96 = cannot use log scale with negative data values
!>  97 = array only has a single value
!>  98 = interval has 0 width or values are inverted
!>  99 = unknown error
!>
!> bad output values use MISSING_I and MISSING_R8
!>
!> FIXME:
!> added to the utilities module, but this module should be split into
!> smaller modules because right now it's a dumping ground for every
!> random routine that is useful to more than one module.  (my fault
!> as much as anyones - nsc)

subroutine find_enclosing_indices(nitems, data_array, value_to_find,     &
                                  smaller_index, larger_index, fraction_across, my_status, &
                                  inverted, log_scale)

integer,  intent(in)  :: nitems
real(r8), intent(in)  :: data_array(nitems)
real(r8), intent(in)  :: value_to_find
integer,  intent(out) :: smaller_index
integer,  intent(out) :: larger_index
real(r8), intent(out) :: fraction_across
integer,  intent(out) :: my_status
logical,  intent(in), optional :: inverted
logical,  intent(in), optional :: log_scale

integer :: i, j, k
integer :: smaller_val_i, larger_val_i
logical :: one_is_smallest, linear_interp    ! the normal defaults

! set defaults and initialize intent(out) items
! so we can return immediately on error.

one_is_smallest = .true.
if (present(inverted)) one_is_smallest = .not. inverted

linear_interp = .true.
if (present(log_scale)) linear_interp = .not. log_scale

smaller_index = MISSING_I
larger_index  = MISSING_I
fraction_across = MISSING_R8
my_status = -99


! exclude malformed call cases
if (nitems <= 1) then
   my_status = 97
   return
endif

! set these indices so we can simplify the tests below
if (one_is_smallest) then
   smaller_val_i = 1
   larger_val_i  = nitems
else
   smaller_val_i = nitems
   larger_val_i  = 1
endif
   
! discard out of range values
if (value_to_find < data_array(smaller_val_i)) then
   my_status = -1
   return
endif

if (value_to_find > data_array(larger_val_i)) then
   my_status = 1
   return
endif

! bisection section (get it?)
i = 1
j = nitems

do
   k=(i+j)/2
   if ((value_to_find < data_array(k) .and.       one_is_smallest) .or. &
       (value_to_find > data_array(k) .and. .not. one_is_smallest)) then
      j=k
   else
      i=k
   endif
   if (i+1 >= j) exit
enddo

! return vals
smaller_index = i
larger_index  = i+1

! again, set these indices so we can simplify the tests below
! to compute fraction_across
if (one_is_smallest) then
   smaller_val_i = i
   larger_val_i  = i+1
else
   smaller_val_i = i+1
   larger_val_i  = i
endif
   
! avoid divide by 0.  at this point we'll be returning valid index values
! but returning a bad fraction and status code.  might help identify where
! in the incoming data array things are bad.
if ((data_array(larger_val_i) - data_array(smaller_val_i)) <= 0.0_r8) then
   my_status = 98
   return
endif

! no log computations if any data values are negative
! (do this on 2 lines to avoid testing the data value
!  unless we are planning to take the log.)
if (.not. linear_interp) then
   if (data_array(smaller_val_i) <= 0.0) then
      my_status = 96
      return
   endif
endif

! compute fraction here
if (linear_interp) then
   fraction_across = (value_to_find            - data_array(smaller_val_i)) / &
                     (data_array(larger_val_i) - data_array(smaller_val_i))
else
   fraction_across = (log(value_to_find)            - log(data_array(smaller_val_i))) / &
                     (log(data_array(larger_val_i)) - log(data_array(smaller_val_i)))

endif

! good return
my_status = 0

end subroutine find_enclosing_indices


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! debug code section
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 
!-----------------------------------------------------------------------
!>  Useful for dumping all the attributes for a file 'unit'
!>  A debugging routine, really. TJH Oct 2004

   subroutine dump_unit_attributes(iunit) 

      integer, intent(in) :: iunit

      logical :: exists, connected, named_file
      character(len=256) :: file_name
      character(len=512) :: str1
      character(len=32)  :: ynu     ! YES, NO, UNDEFINED ... among others
      integer :: ios, reclen, nextrecnum
character(len=*), parameter :: routine = "dump_unit_attributes"

      if ( .not. module_initialized ) call initialize_utilities

      write(str1,*)'for unit ',iunit 
call error_handler(E_MSG, routine, str1, source, revision, revdate)

      inquire(iunit, opened = connected, iostat=ios)
      if ( connected .and. (ios == 0) ) &
   call error_handler(E_MSG, routine, ' connected', source, revision, revdate)

      inquire(iunit, named = named_file, iostat=ios)
      if ( named_file .and. (ios == 0) ) &
   call error_handler(E_MSG, routine, ' file is named.', source, revision, revdate)

      inquire(iunit, name = file_name, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'file name is ' // trim(adjustl(file_name))
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, exist = exists, iostat=ios)
      if ( exists .and. (ios == 0) ) &
   call error_handler(E_MSG, routine, ' file exists', source, revision, revdate)

      inquire(iunit, recl = reclen, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'record length is ', reclen
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, nextrec = nextrecnum, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'next record is ', nextrecnum
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, access = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'access_type is ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, sequential = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'is file sequential ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, direct = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'is file direct ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, form = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'file format ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, action = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'action ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, read = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'read ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, write = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'write ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, readwrite = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'readwrite ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, blank = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'blank ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, position = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'position ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, delim = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'delim ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

      inquire(iunit, pad = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'pad ', ynu
   call error_handler(E_MSG, routine, str1, source, revision, revdate)
      endif

   end subroutine dump_unit_attributes


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!=======================================================================
! End of utilities_mod
!=======================================================================

end module utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
