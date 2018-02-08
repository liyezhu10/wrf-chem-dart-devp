! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program forcing_check

! program to take a netCDF file ...

use     types_mod, only : r4, r8
use utilities_mod, only : register_module, error_handler, E_ERR, E_MSG,       &
                          open_file, close_file, get_next_filename, &
                          find_namelist_in_file, check_namelist_read,         &
                          do_nml_file, do_nml_term, nmlfileunit,              &
                          initialize_utilities, finalize_utilities
use parse_args_mod, only : get_args_from_string

use netcdf_utilities_mod, only : nc_check

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'

character(len=*), parameter :: routine = 'forcing_check'

! variables used to read the netcdf info
integer, parameter :: maxd = 7
integer :: i, j, ndims, odims, ncrc, etype, nitems, nvars, xtype
integer :: ncinid1
integer :: invarid1
integer :: dimid(maxd), dimlen(maxd), odimid(maxd), odimlen(maxd)
character(128) :: dimname(maxd), odimname(maxd)
integer :: nin1Dimensions, nin1Variables, nin1Attributes, in1unlimitedDimID

! arrays for all possible dimensions, real and int
real(r4)              ::  r4missing
real(r4)              ::  zerod1
real(r4), allocatable ::   oned1(:)
real(r4), allocatable ::   twod1(:,:)
real(r4), allocatable :: threed1(:,:,:)
real(r4), allocatable ::  fourd1(:,:,:,:)
real(r4), allocatable ::  fived1(:,:,:,:,:)
real(r4), allocatable ::   sixd1(:,:,:,:,:,:)
real(r4), allocatable :: sevend1(:,:,:,:,:,:,:)

integer               :: imissing
integer               ::  izerod1
integer,  allocatable ::   ioned1(:)
integer,  allocatable ::   itwod1(:,:)
integer,  allocatable :: ithreed1(:,:,:)
integer,  allocatable ::  ifourd1(:,:,:,:)
integer,  allocatable ::  ifived1(:,:,:,:,:)
integer,  allocatable ::   isixd1(:,:,:,:,:,:)
integer,  allocatable :: isevend1(:,:,:,:,:,:,:)

logical, save :: module_initialized = .false.

! arg parsing code
character(len=256) :: argline
integer :: argcount = 1
character(len=NF90_MAX_NAME) :: argwords(3)

character(len=NF90_MAX_NAME) :: infile1
character(len=NF90_MAX_NAME) :: nextfield
logical :: from_file

character(len=512) :: string1, string2, string3

integer :: iunit, io

interface process
   procedure process_int_1d
   procedure process_int_2d
   procedure process_int_3d
   procedure process_r4_1d
   procedure process_r4_2d
   procedure process_r4_3d
end interface

logical :: debug = .false.                   ! or .true.
logical :: fail_on_missing_field = .true.    ! or .false.
logical :: do_all_numeric_fields = .true.    ! or .false.
logical :: only_report_differences = .true.  ! or .false.
character(len=NF90_MAX_NAME) :: fieldnames(1000) = ''  ! something large
character(len=256) :: fieldlist_file = ''

! fieldnames here?
namelist /forcing_check_nml/ debug, fail_on_missing_field, &
                             do_all_numeric_fields,        &
                             fieldnames, fieldlist_file,   &
                             only_report_differences

! flow:
!   initialization
call initialize_utilities(routine)
call initialize_module()

! Read the namelist entry
call find_namelist_in_file('input.nml', 'forcing_check_nml', iunit)
read(iunit, nml = forcing_check_nml, iostat = io)
call check_namelist_read(iunit, io, 'forcing_check_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=forcing_check_nml)
if (do_nml_term()) write(     *     , nml=forcing_check_nml)

if (debug) then
   call error_handler(E_MSG, routine, ' debug on')
endif

! whether to fail or just warn if a field is not found
if (fail_on_missing_field) then
   etype = E_ERR
else
   etype = E_MSG
endif

!   check inputs - get a string from stdin
!   eg:  echo infile1.nc | ./forcing_check
read(*, '(A)') argline
call get_args_from_string(argline, argcount, argwords)

if (argcount /= 1) then
   string1 = 'Usage: echo infile1.nc | ./forcing_check'
   call error_handler(E_ERR, routine, string1, source, revision, revdate) 
endif

infile1 = argwords(1)

call error_handler(E_MSG, routine, ' reading file: '//trim(infile1))

if (do_all_numeric_fields) then
   call error_handler(E_MSG, routine, ' doing all numeric fields')
else
   ! make sure the namelist specifies one or the other but not both
   ! only if we aren't trying to do all fields in the file.
   if (fieldnames(1) /= '' .and. fieldlist_file /= '') then
      call error_handler(E_ERR, routine, &
          'cannot specify both fieldnames and fieldlist_file', &
          source,revision,revdate)
   endif
   
   if (fieldlist_file /= '') then
      call error_handler(E_MSG, routine, 'list of fields file: '//trim(fieldlist_file))
      from_file = .true.
   else
      call error_handler(E_MSG, routine, 'field names specified in namelist.')
      from_file = .false.
   endif
endif

! open the files
ncrc = nf90_open(infile1, NF90_NOWRITE,   ncinid1)
call nc_check(ncrc, routine, 'nf90_open', infile1)

ncrc = nf90_inquire(ncinid1, nin1Dimensions, nin1Variables, &
                            nin1Attributes, in1unlimitedDimID)
call nc_check(ncrc, routine, 'nf90_inquire', infile1)

! for now, loop over the number of vars in file 1.  at some point we
! should print out vars that are in 1 but not 2, and in 2 but not 1.
! but for that we need more logic to track which vars are processed.
! this code loops over the names in file 1 and finds them (or not)
! in file 2 but doesn't complain about unused vars in file 2.
nvars = nin1Variables

if (debug) then
   write(string1, *) 'infile1 ndim, nvar, nattr:', nin1Dimensions, &
                      nin1Variables, nin1Attributes
   call error_handler(E_MSG, routine, string1)
endif

!    input files to get data from
!    list of netcdf fields to compare

fieldloop : do i=1, nvars

   ! get the variable name of interest
   if (do_all_numeric_fields) then
      ncrc = nf90_inquire_variable(ncinid1, i, nextfield)
      call nc_check(ncrc, routine, 'nf90_inquire_variable', nextfield, infile1)
   else
      if (from_file) then
         nextfield = get_next_filename(fieldlist_file, i)
      else
         nextfield = fieldnames(i)
      endif
      if (nextfield == '') exit fieldloop
   endif

   ! check if variable exists, get the ID 
   ncrc = nf90_inq_varid(ncinid1, nextfield, invarid1)
   if (ncrc /= NF90_NOERR) then
      string2 = ' not found in input file "'//trim(infile1)//'"'
      if (etype == E_ERR) then
         string1 = 'variable "'//trim(nextfield)//'"'
      else
         string1 = 'skipping variable '//trim(nextfield)//'"'
      endif
      call error_handler(etype, routine, string1, source, revision, revdate, text2=string2) 
      cycle fieldloop
   endif
   
   ncrc = nf90_inquire_variable(ncinid1, invarid1, xtype=xtype)
   call nc_check(ncrc, routine, 'inquire for xtype', nextfield, infile1)

   if (xtype /= NF90_INT .and. xtype /= NF90_FLOAT .and. xtype /= NF90_DOUBLE) then
      string1 = 'skipping variable "'//trim(nextfield)//'"'
      string2 = 'not integer, float, or double'
      call error_handler(E_MSG, routine, string1, &
                         source, revision, revdate, text2=string2) 
      cycle fieldloop
   endif

   if (debug) then
      write(string1, *) 'invarid1: ',invarid1, trim(nextfield)//' '//trim(infile1)
      call error_handler(E_MSG, routine, string1)
   endif

   ! get dimensions

   ncrc = nf90_inquire_variable(ncinid1, invarid1, ndims=ndims,  dimids=dimid)
   call nc_check(ncrc, routine, 'nf90_inquire_variable', nextfield, infile1)

   do j=1,ndims

      write(string3,*)j, trim(dimname(j)), trim(nextfield)

      ncrc = nf90_inquire_dimension(ncinid1,  dimid(j),  dimname(j),  dimlen(j))
      call nc_check(ncrc, routine, 'nf90_inquire_dimension', string3, infile1) 

      if (debug) then
         write(string1, '(2A,I5,A,I8,2A)') trim(infile1), ' dim: ', j, ' len: ', dimlen(j), ' name: ', trim(dimname(j))
         call error_handler(E_MSG, routine, string1)
      endif
      
   enddo

   select case(ndims)
      case (0)
         write(string3,*)       trim(nextfield), ' [scalar value]'
      case (1)
         continue
      case (2)
         continue
      case (3)
         continue
      case (4)
         write(string3,'(A,"(",4I6,")")') trim(nextfield), dimlen(1:4)
      case (5)
         write(string3,'(A,"(",5I6,")")') trim(nextfield), dimlen(1:5)
      case (6)
         write(string3,'(A,"(",6I6,")")') trim(nextfield), dimlen(1:6)
      case (7)
         write(string3,'(A,"(",7I6,")")') trim(nextfield), dimlen(1:7)
      case default
         ! "can't happen"
         write(string1, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
   end select

   ! allocate right dim array
   ! read/write and then deallocate
   ! this section has lot of replicated code, but with the strongly typed
   ! arrays it is hard to do much else.  we could overload a subroutine
   ! with each dimension but we'd end up with the same amount of replicated code

   select case(xtype)
    case(NF90_INT)
     select case(ndims)
      case (0)
         ncrc = nf90_get_var(ncinid1, invarid1, izerod1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
      case (1)
         allocate(ioned1(dimlen(1)))
         ncrc = nf90_get_var(ncinid1, invarid1, ioned1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, imissing, ioned1, nitems)
         deallocate(ioned1)
      case (2)
         allocate(itwod1(dimlen(1),dimlen(2)))
         ncrc = nf90_get_var(ncinid1, invarid1, itwod1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, imissing, itwod1, nitems)
         deallocate(itwod1)
      case (3)
         allocate(ithreed1(dimlen(1),dimlen(2),dimlen(3)))
         ncrc = nf90_get_var(ncinid1, invarid1, ithreed1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, imissing, ithreed1, nitems)
         deallocate(ithreed1)
      case (4)
         allocate(ifourd1(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
         ncrc = nf90_get_var(ncinid1, invarid1, ifourd1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(ifourd1)
      case (5)
         allocate(ifived1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
         ncrc = nf90_get_var(ncinid1, invarid1, ifived1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(ifived1)
      case (6)
         allocate(isixd1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
         ncrc = nf90_get_var(ncinid1, invarid1, isixd1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(isixd1)
      case (7)
         allocate(isevend1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6),dimlen(7)))
         ncrc = nf90_get_var(ncinid1, invarid1, isevend1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(isevend1)
      case default
         ! "really can't happen"
         write(string1, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
     end select

     ! common reporting code for integers
     if (nitems > 0) then
        write(string1, *) 'arrays differ in ', nitems, ' places'
        call error_handler(E_MSG, routine, string1, source, revision, revdate)
     else if (.not. only_report_differences) then
        write(string1, *) 'arrays same'
        call error_handler(E_MSG, routine, string1, source, revision, revdate)
     endif

    case default
     select case(ndims)
      case (0)
         ncrc = nf90_get_var(ncinid1, invarid1, zerod1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
      case (1)
         allocate(oned1(dimlen(1)))
         ncrc = nf90_get_var(ncinid1, invarid1, oned1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, r4missing, oned1, nitems)
         deallocate(oned1)
      case (2)
         allocate(twod1(dimlen(1),dimlen(2)))
         ncrc = nf90_get_var(ncinid1, invarid1, twod1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, r4missing, twod1, nitems)
         deallocate(twod1)
      case (3)
         allocate(threed1(dimlen(1),dimlen(2),dimlen(3)))
         ncrc = nf90_get_var(ncinid1, invarid1, threed1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         call process(nextfield, r4missing, threed1, nitems)
         deallocate(threed1)
      case (4)
         allocate(fourd1(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
         ncrc = nf90_get_var(ncinid1, invarid1, fourd1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(fourd1)
      case (5)
         allocate(fived1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
         ncrc = nf90_get_var(ncinid1, invarid1, fived1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(fived1)
      case (6)
         allocate(sixd1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
         ncrc = nf90_get_var(ncinid1, invarid1, sixd1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(sixd1)
      case (7)
         allocate(sevend1(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6),dimlen(7)))
         ncrc = nf90_get_var(ncinid1, invarid1, sevend1)
         call nc_check(ncrc, routine, 'nf90_get_var', nextfield, infile1)
         deallocate(sevend1)
      case default
         ! "really can't happen"
         write(string1, *) 'array dimension is illegal value: ', ndims
         call error_handler(E_ERR, routine, string1, source, revision, revdate)
     end select

     ! common reporting code for reals
     if (nitems > 0) then
        write(string1, *) 'arrays differ in ', nitems, ' places'
        call error_handler(E_MSG, routine, string1, source, revision, revdate)
     else if (.not. only_report_differences) then
        write(string1, *) 'arrays same'
        call error_handler(E_MSG, routine, string1, source, revision, revdate)
     endif
   end select

enddo fieldloop

!  close up
call nc_check(nf90_close(ncinid1), 'nf90_close', infile1)

if (debug) then
   write(string1, *) 'closing files ',  trim(infile1)
   call error_handler(E_MSG, routine, string1)
endif

call finalize_utilities(routine)

! end of main code


contains

!----------------------------------------------------------------------

subroutine initialize_module

  call register_module(source, revision, revdate)
  module_initialized = .true.

end subroutine initialize_module


!----------------------------------------------------------------------

subroutine process_int_1d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
integer,          intent(in)    :: missing
integer,          intent(inout) :: datmat(:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), minimum, maximum

end subroutine process_int_1d

!----------------------------------------------------------------------

subroutine process_int_2d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
integer,          intent(in)    :: missing
integer,          intent(inout) :: datmat(:,:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), minimum, maximum

end subroutine process_int_2d

!----------------------------------------------------------------------

subroutine process_int_3d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
integer,          intent(in)    :: missing
integer,          intent(inout) :: datmat(:,:,:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), minimum, maximum

end subroutine process_int_3d

!----------------------------------------------------------------------

subroutine process_r4_1d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
real(r4),         intent(in)    :: missing
real(r4),         intent(inout) :: datmat(:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), shape(datmat), 'minmax:', minimum, maximum

end subroutine process_r4_1d

!----------------------------------------------------------------------

subroutine process_r4_2d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
real(r4),         intent(in)    :: missing
real(r4),         intent(inout) :: datmat(:,:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), shape(datmat), 'minmmax:', minimum, maximum

end subroutine process_r4_2d

!----------------------------------------------------------------------

subroutine process_r4_3d(varname, missing, datmat, nitems)

character(len=*), intent(in)    :: varname
real(r4),         intent(in)    :: missing
real(r4),         intent(inout) :: datmat(:,:,:)
integer,          intent(out)   :: nitems

integer :: minimum, maximum

nitems  = 0
minimum = minval(datmat)
maximum = maxval(datmat)

write(*,*)trim(varname), shape(datmat), 'minmmax:', minimum, maximum

end subroutine process_r4_3d

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
