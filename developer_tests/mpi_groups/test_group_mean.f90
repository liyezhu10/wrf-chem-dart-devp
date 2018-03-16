! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_group_mean

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                  task_count, my_task_id, task_sync, set_group_size,&
                                  get_group_size

use      time_manager_mod, only : time_type, set_time, print_date, operator(-), &
                                  NO_CALENDAR

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, print_ens_handle

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains, get_model_variable_indices, &
                                  state_structure_info, add_domain

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename, set_file_metadata, &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use distributed_state_mod, only : create_state_window, free_state_window, &
                                  create_mean_window, free_mean_window

use             model_mod, only : static_init_model

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, parameter :: MAX_TESTS = 7

! this is max number of domains times number of ensemble members
! if you have more than one domain and your ensemble members are
! in separate files, the names should be listed in this order:
!  all filenames for ensemble members for domain 1
!  all filenames for ensemble members for domain 2, etc

integer, parameter :: MAX_FILES = 1000

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

logical                       :: single_file = .false.
integer                       :: num_ens = 1
integer                       :: dtype = 1
integer                       :: ltype = 1
integer                       :: ttype = 1
integer                       :: group_size
character(len=256)            :: input_state_files(MAX_FILES)  = 'null'
character(len=256)            :: output_state_files(MAX_FILES) = 'null'


namelist /test_group_mean_nml/ num_ens, single_file, input_state_files, output_state_files, &
                               dtype, ltype, ttype, group_size

! io variables
integer                   :: iunit, io
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.

! model state variables
type(ensemble_type)   :: ens_handle

type(time_type)       :: model_time
integer(i8)           :: model_size

! misc. variables
integer :: i, idom, imem, domid, num_domains

! message strings
character(len=512) :: my_base, my_desc
character(len=512) :: string1

character(len=256), allocatable :: file_array_input(:,:)
character(len=256), dimension(1) :: var_names = (/'temp'/)
integer,parameter :: one_domain = 1

!======================================================================
! start of executable code
!======================================================================

call initialize_modules_used()

call find_namelist_in_file("input.nml", "test_group_mean_nml", iunit)
read(iunit, nml = test_group_mean_nml, iostat = io)
call check_namelist_read(iunit, io, "test_group_mean_nml")

!----------------------------------------------------------------------
! read/write restart files
!----------------------------------------------------------------------

model_size = 60

if(my_task_id()==0) then
   print*, 'dtype ', dtype
   print*, 'ltype ', ltype
   print*, 'ttype ', ttype
endif

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle,                &
                           num_copies           = 1,  &
                           num_vars             = model_size, &
                           distribution_type_in = dtype, & ! round robin, pair round robin, block
                           layout_type          = ltype, & ! no vars, transposable, transpose and duplicate
                           transpose_type_in    = ttype)   ! no vars, transposable, transpose and duplicate

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)

num_domains = 1

allocate(file_array_input( num_ens, num_domains))
file_array_input  = RESHAPE(input_state_files,  (/num_ens,  num_domains/))

domid = add_domain('simple1.nc', 1, var_names)

! Test the read portion.
call io_filenames_init(file_info_input,              &
                       ncopies       = num_ens,      &
                       cycling       = single_file,  &
                       single_file   = single_file,  &
                       restart_files = file_array_input)

do imem = 1, num_ens
   write(my_base,'(A,I2)') 'inens_',    imem
   write(my_desc,'(A,I2)') 'input ens', imem
   call set_file_metadata(file_info_input,                      &
                          cnum     = imem,                      &
                          fnames   = file_array_input(imem,:),  &
                          basename = my_base,                   &
                          desc     = my_desc)

   call set_io_copy_flag(file_info_input,    &
                         cnum    = imem,     &
                         io_flag = READ_COPY)
enddo

input_restart_files = get_stage_metadata(file_info_input)

do imem = 1, num_ens
   write(string1, *) '- Reading File : ', trim(get_restart_filename(input_restart_files, imem, domain=one_domain))
enddo

call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

call print_ens_handle(ens_handle, force=.true., label='test_mean', contents=.false., limit=5)

deallocate(file_array_input)

! do i = 0, task_count()-1
!    call task_sync()
!    if(my_task_id() == i) then 
!       call print_ens_handle(ens_handle)
!    endif
! enddo

call task_sync()

!----------------------------------------------------------------------
! Check window
!----------------------------------------------------------------------

call create_mean_window(ens_handle, mean_copy=1, distribute_mean=.true.)

call free_mean_window()

!----------------------------------------------------------------------
! finalize test_group_mean
!----------------------------------------------------------------------

call finalize_modules_used()

!======================================================================
contains
!======================================================================


!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('test_group_mean')

call register_module(source,revision,revdate)

call state_vector_io_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------
!> clean up before exiting

subroutine finalize_modules_used()

! this must be last, and you can't print/write anything
! after this is called.
call finalize_mpi_utilities()

end subroutine finalize_modules_used

!----------------------------------------------------------------------
!> print the vars array

subroutine print_ens(ens_handle, lim, rank)
type(ensemble_type), intent(in) :: ens_handle
integer,             intent(in) :: rank
integer,             intent(in) :: lim
integer :: i

print*, 'rank :: ', rank
do i = 1, ens_handle%my_num_vars
   write(*,'(A,I3,A,i4)'), 'my_vars(', i,') = ', ens_handle%my_vars(i)
enddo
write(*,*) ''

end subroutine print_ens

end program test_group_mean

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
