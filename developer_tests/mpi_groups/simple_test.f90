! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program simple_test

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                  task_count, my_task_id, task_sync, datasize, &
                                  get_dart_mpi_comm, create_groups

use      time_manager_mod, only : time_type, set_time, print_date, operator(-), &
                                  NO_CALENDAR

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, print_ens_handle, &
                                  all_vars_to_all_copies, all_copies_to_all_vars, &
                                  my_vars_to_group_copies

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains, get_model_variable_indices, &
                                  state_structure_info, add_domain

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename, set_file_metadata, &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use distributed_state_mod, only : create_state_window, free_state_window, &
                                  create_mean_window, free_mean_window,   &
                                  get_state

use             model_mod, only : static_init_model

use netcdf

use mpi

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
logical                       :: debug      = .false.
integer                       :: group_size = 2
integer(KIND=MPI_OFFSET_KIND) :: NX         = 8 !< lengths of dimensions
integer                       :: max_iter   = 1000
integer                       :: dtype = 1
integer                       :: ltype = 1
integer                       :: ttype = 1

namelist /simple_test_nml/ NX, group_size, max_iter, dtype, ltype, ttype, debug

! io variables
integer                   :: iunit, io
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.

! model state variables
type(ensemble_type)   :: state_handle
type(ensemble_type)   :: mean_handle
type(ensemble_type)   :: group_mean_handle

! misc. variables
integer :: i

! message strings
character(len=512) :: my_base, my_desc
character(len=512) :: string1

character(len=256), allocatable  :: file_array_input(:,:)
character(len=256), dimension(1) :: var_names = (/'temp'/)
integer,parameter :: one_domain = 1

real(r8) ::  u, my_val

integer, allocatable :: group_members(:)
integer :: dart_group
integer :: sub_group
integer :: group_comm
integer :: group_rank

! grid window
integer               :: my_window
real(r8), allocatable :: my_array(:) !< local get_my_val info
real(r8)              :: duplicate_array(*)
pointer(aa, duplicate_array)

! index variables
integer :: ii, jj

! timing variables
real(r8) :: t1, t2, max_time, min_time

! MPI variables
integer :: ierr
integer :: my_rank

!======================================================================
! start of executable code
!======================================================================
call initialize_modules_used()

call find_namelist_in_file("input.nml", "simple_test_nml", iunit)
read(iunit, nml = simple_test_nml, iostat = io)
call check_namelist_read(iunit, io, "simple_test_nml")

if (group_size == -1) group_size = task_count()

!----------------------------------------------------------------------
! create groups
!----------------------------------------------------------------------
t1 = MPI_WTIME()
call create_groups(group_size, group_comm)
t2 = MPI_WTIME() - t1

call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, get_dart_mpi_comm(), ierr)
call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, get_dart_mpi_comm(), ierr)

if (my_task_id() == 0) then
   print*, 'group_size = ', group_size, ', NX = ', NX
   print*, 'create_groups : Max Time = ', max_time
   print*, 'create_groups : Min Time = ', min_time
endif

!----------------------------------------------------------------------
! create data array on task 0
!----------------------------------------------------------------------

! Set up the ensemble storage for mean
call init_ensemble_manager(group_mean_handle,               &
                           num_copies           = 1,  &
                           num_vars             = NX, &
                           distribution_type_in = dtype, & ! 1 - round robin, 2 - pair round robin, 3 - block
                           layout_type          = ltype, & ! 1 - standard,    2 - round-robbin,     3 - distribute mean in groups
                           transpose_type_in    = ttype, & ! 1 - no vars,     2 - transposable,     3 - transpose and duplicate
                           mpi_comm             = group_comm)  

call init_ensemble_manager(state_handle,              &
                           num_copies           = 1,  &
                           num_vars             = NX, &
                           distribution_type_in = 1,  & ! 1 - round robin, 2 - pair round robin, 3 - block
                           layout_type          = 1,  & ! 1 - standard,    2 - round-robbin,     3 - distribute mean in groups
                           transpose_type_in    = 2,  & ! 1 - no vars,     2 - transposable,     3 - transpose and duplicate
                           mpi_comm             = get_dart_mpi_comm())  
                           

! first task has the state_handle vars then sends copies to 
! other tasks doing a transpose (i.e all_vars_to_all_copies)
if (my_task_id() == 0) then
   do i = 1, NX
      ! Dimensioned (num_vars, my_num_copies)
      state_handle%vars(i, 1) = i
   enddo
endif

!print*, my_task_id(), 'initial transpose state_handle%copies  (:,:)',state_handle%copies(:,:)
!print*, my_task_id(), 'initial transpose state_handle%vars    (:,:)',state_handle%vars(:,:)

! let all tasks have a copy of the mean
call all_vars_to_all_copies(state_handle, label='state_handle%vars(:,:) -> state_handle%copies(:,:)')

call print_ens_handle(state_handle,             &
                      force    = .true.,        &
                      label    = 'state_handle', &
                      contents = .true.)

! mean_handle is optionally returned when creating mean window.
call create_mean_window(state_handle, mean_copy=1, distribute_mean=.false., state_mean_ens_handle=mean_handle)

call my_vars_to_group_copies(mean_handle, group_mean_handle, label='my_vars_to_group_copies')

call print_ens_handle(mean_handle,             &
                      force    = .true.,        &
                      label    = 'mean_handle', &
                      contents = .true.)

call print_ens_handle(group_mean_handle,             &
                      force    = .true.,        &
                      label    = 'group_mean_handle', &
                      contents = .true.)

!print*, 'actav state_handle%copies(:,:)',state_handle%copies(:,:)
!print*, 'actav state_handle%vars  (:,:)',state_handle%vars(:,:)

!call all_vars_to_all_copies(mean_handle)

!print*, 'avtac state_handle%copies(:,:)',state_handle%copies(:,:)
!print*, 'avtac state_handle%vars  (:,:)',state_handle%vars(:,:)

call free_mean_window()

! call print_ens_handle(mean_handle,              &
!                       force    = .true.,        &
!                       label    = 'mean_handle', &
!                       contents = .true.)
! print*, my_task_id(), 'AFTER PRINT'

! allocate(my_array(NX/group_size))
! 
! do ii = 1, NX/group_size
!    my_array(ii) = group_rank*NX/group_size + ii
! enddo
! 
! !----------------------------------------------------------------------
! ! create window
! !----------------------------------------------------------------------
! call create_window()
! 
! !----------------------------------------------------------------------
! ! timing test
! !----------------------------------------------------------------------
! t2 = 0.0_r8
! do ii = 1, max_iter
!    call random_number(u)
!    jj = FLOOR(NX*u)+1
!    t1 = MPI_WTIME()
!    my_val = get_my_val(jj)
!    t2 = t2 + (MPI_WTIME() - t1)
!    if (my_val /= jj) then
!       print*, 'jj /= my_val', jj, my_val
!    endif
! enddo
! 
! call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, get_dart_mpi_comm(), ierr)
! call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, get_dart_mpi_comm(), ierr)
! if (my_task_id() == 0) then
!    print*, 'get_value     : Max Time = ', max_time
!    print*, 'get_value     : Min Time = ', min_time
! endif
! 
! !----------------------------------------------------------------------
! ! print array
! !----------------------------------------------------------------------
! if (debug) then
!    do jj = 0, task_count()-1
!       if (my_task_id() == jj) then
!          print*, my_task_id(), '::', my_array(:)
!       else
!          call task_sync()
!       endif
!    enddo
! endif
! 
! !----------------------------------------------------------------------
! ! debug get_my_val
! !----------------------------------------------------------------------
! do jj = 0, task_count()-1
!    if (my_task_id() == jj .and. debug) then
!       do ii = 1, NX
!          print*, 'my_task_id() = ', my_task_id(), 'get_owner(ii)', ii, get_owner(ii, my_task_id()), get_my_val(ii)
!       enddo
!    else
!       call task_sync()
!    endif
! enddo
! 
! call task_sync()
! 
! !----------------------------------------------------------------------
! ! finalize simple_test
! !----------------------------------------------------------------------
! call free_window()

call task_sync()

call finalize_modules_used()

!======================================================================
contains
!======================================================================

!#! !-----------------------------------------------------------
!#! !> create window for get_my_val
!#! !> Need to add a loop for domain.  Just single domain for now.
!#! !> Window is filled with columns contiguous in memory, possibly
!#! !> we will grab a whole column at once.
!#! !> 
!#! !> If you want to change the order in the window, be sure to 
!#! !> change the subroutine who_has_grid_info
!#! subroutine create_window
!#! 
!#! integer ii, jj, kk, sizedouble, count
!#! integer(KIND=MPI_ADDRESS_KIND) :: window_size
!#! 
!#! ! datasize comes from mpi_utilities_mod
!#! call mpi_type_size(datasize, sizedouble, ierr)
!#! 
!#! window_size = (NX/group_size)*sizedouble
!#! aa          = malloc(NX)
!#! 
!#! call MPI_ALLOC_MEM(window_size, mpi_info_null, aa, ierr)
!#! 
!#! ! can't do my_array assignment with a cray pointer, so you need to loop
!#! !#! count = 1
!#! !#! do ii = 1, NX
!#! !#!    duplicate_array(count) = my_array(ii)
!#! !#!    count = count + 1
!#! !#! enddo
!#! 
!#! call mpi_win_create(my_array,        window_size,   &
!#!                     sizedouble,      MPI_INFO_NULL, &
!#!                     group_comm,   my_window,     ierr)
!#! 
!#! end subroutine create_window
!#! 
!#! !---------------------------------------------------------
!#! !> Free the mpi windows you created
!#! subroutine free_window
!#! integer :: ierr
!#! 
!#! call mpi_win_free(my_window, ierr)
!#! call MPI_FREE_MEM(duplicate_array, ierr)
!#! 
!#! deallocate(group_members, my_array)
!#! 
!#! end subroutine free_window

!---------------------------------------------------------
!> Function to get get_my_val
function get_my_val(i)

real(r8) :: get_my_val
integer, intent(in) :: i

integer                          :: owner !< which task has the part of get_my_val we need
integer(KIND=MPI_ADDRESS_KIND)   :: target_disp !< displacement

! caluclate who has the info
owner = get_owner(i, my_task_id())
target_disp = mod(i,NX/group_size)

if (debug) then
   print*, 'my_task_id(), owner, my_window', my_task_id(), owner, my_window
endif

! grab the info
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, my_window, ierr)
call mpi_get(get_my_val, 1, datasize, owner, target_disp, 1, datasize, my_window, ierr)
call mpi_win_unlock(owner, my_window, ierr)

end function

!---------------------------------------------------------
!> Get the group number
function get_owner(i, pe)
integer, intent(in) :: i
integer, intent(in) :: pe
integer :: get_owner

integer :: num_groups

num_groups = NX / group_size
get_owner  = (i-1)/num_groups

end function get_owner

!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('simple_test')

call register_module(source,revision,revdate)

call state_vector_io_init()

end subroutine initialize_modules_used

!----------------------------------------------------------------------
!> clean up before exiting

subroutine finalize_modules_used()

! this must be last, and you can't print/write anything
! after this is called.
print*, my_task_id(), ' calling mpi_finialize'
call finalize_mpi_utilities()

end subroutine finalize_modules_used

end program simple_test

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
