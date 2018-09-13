! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program simple_test

use             types_mod, only : r8

use         utilities_mod, only : register_module, E_MSG, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                                  task_count, my_task_id, task_sync, datasize, &
                                  get_dart_mpi_comm, create_groups

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, print_ens_handle, &
                                  all_vars_to_all_copies, all_copies_to_all_vars, &
                                  my_vars_to_group_copies

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename, set_file_metadata

use distributed_state_mod, only : create_mean_window, free_mean_window

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

! model state variables
type(ensemble_type)   :: state_handle
type(ensemble_type)   :: mean_handle
type(ensemble_type)   :: group_mean_handle

! misc. variables
integer  :: i, j
real(r8) ::  u, my_val

integer, allocatable :: group_members(:)
integer :: dart_group
integer :: sub_group
integer :: group_comm
integer :: group_rank

! group window
integer               :: rank
integer               :: my_window
integer(KIND=MPI_ADDRESS_KIND) :: window_size
real(r8), allocatable :: my_array(:) !< local get_my_val info
real(r8)              :: duplicate_array(*)
pointer(aa, duplicate_array)

! index variables
integer :: ii, jj

! timing variables
real(r8) :: t1, t2, max_time, min_time

! MPI variables
integer :: ierr
integer :: sizedouble

!======================================================================
! start of executable code
!======================================================================
call initialize_modules_used()

call find_namelist_in_file("input.nml", "simple_test_nml", iunit)
read(iunit, nml = simple_test_nml, iostat = io)
call check_namelist_read(iunit, io, "simple_test_nml")

!print*, 'SETTING group_size = ', group_size
if (group_size == -1) group_size = task_count()

!----------------------------------------------------------------------
! create groups
!----------------------------------------------------------------------
t1 = MPI_WTIME()
!print*, 'CALLING create_groups'
!call local_create_groups()
call create_groups(group_size, group_comm)
t2 = MPI_WTIME() - t1

call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, get_dart_mpi_comm(), ierr)
call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, get_dart_mpi_comm(), ierr)

!----------------------------------------------------------------------
! create data array on task 0
!----------------------------------------------------------------------
! distribution_type_in :: 1 - round robin, 2 - pair round robin, 3 - block
! layout_type          :: 1 - standard,    2 - round-robbin,     3 - distribute mean in groups
! transpose_type       :: 1 - no vars,     2 - transposable,     3 - transpose and duplicate
!print*, 'INIT_ENSEMBLE_MANAGER state_handle'
call init_ensemble_manager(state_handle,              &
                           num_copies           = 1,  &
                           num_vars             = NX, &
                           distribution_type_in = 1,  &
                           layout_type          = 1,  & 
                           transpose_type_in    = 2,  & 
                           mpi_comm             = get_dart_mpi_comm())  
                           
! Set up the ensemble storage for mean
!print*, 'INIT_ENSEMBLE_MANAGER group_mean_handle'
call init_ensemble_manager(group_mean_handle,            &
                           num_copies           = 1,     &
                           num_vars             = NX,    &
                           distribution_type_in = dtype, & ! 1
                           layout_type          = ltype, & ! 1
                           transpose_type_in    = ttype, & ! 1
                           mpi_comm             = group_comm)  


! first task has the state_handle vars then sends copies to 
! other tasks doing a transpose (i.e all_vars_to_all_copies)
if (my_task_id() == 0) then
   !print*, 'SETTING VALUES FOR state_handle%vars'
   do i = 1, NX
      ! Dimensioned (num_vars, my_num_copies)
      state_handle%vars(i, 1) = i
   enddo
endif

! let all tasks have a subset of the mean copy
call all_vars_to_all_copies(state_handle, label='all_vars_to_all_copies')
call task_sync()

if (debug) &
   call print_ens_handle(state_handle,             &
                         force    = .true.,        &
                         label    = 'state_handle', &
                         contents = .true.)

call task_sync()

! mean_handle is optionally returned when creating mean window.
call create_mean_window(state_handle, mean_copy=1, distribute_mean=.false., return_mean_ens_handle=mean_handle)

call task_sync()

if (debug) &
   call print_ens_handle(mean_handle,             &
                         force    = .true.,        &
                         label    = 'mean_handle', &
                         contents = .true.)

call task_sync()

! task_count = 4, NX = 8
! num_copies_to_receive = 1
! num_vars_to_receive = 4

! vars   is from the mean_handle
! copies is from the group_mean_handle
!>@todo FIXME : this should be in the create_mean_window routine
call my_vars_to_group_copies(mean_handle, group_mean_handle, label='my_vars_to_group_copies')

call task_sync()

if (debug) &
   call print_ens_handle(group_mean_handle,             &
                         force    = .true.,        &
                         label    = 'group_mean_handle', &
                         contents = .true.)

call task_sync()

call free_mean_window()

! want to create the mean window so it uses the group communicator
call mpi_type_size(datasize, sizedouble, ierr)
window_size = (NX/group_size)*sizedouble

!>@todo group_mean_handle needs to be contiguous
call mpi_win_create(group_mean_handle%copies(1,:), window_size, &
                    sizedouble, MPI_INFO_NULL, group_comm, my_window, ierr)


do ii = 1, max_iter
   do rank = 0, task_count()-1
      call random_number(u)
      jj = FLOOR(NX*u)+1
      !if (rank == my_task_id()) then
         t1 = MPI_WTIME()
         my_val = get_my_val(jj)
         t2 = t2 + (MPI_WTIME() - t1)
         if (my_val /= jj) then
            print*, 'jj /= my_val', jj, my_val
         endif
        call task_sync()
     !else
     !  call task_sync()
     !endif
   enddo
enddo

!----------------------------------------------------------------------
! timing test
!----------------------------------------------------------------------
call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, get_dart_mpi_comm(), ierr)
call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, get_dart_mpi_comm(), ierr)
if (my_task_id() == 0) then
   print*, 'get_value     : Max Time = ', max_time
   print*, 'get_value     : Min Time = ', min_time
endif


!print*, 'actav state_handle%copies(:,:)',state_handle%copies(:,:)
!print*, 'actav state_handle%vars  (:,:)',state_handle%vars(:,:)

!call all_vars_to_all_copies(mean_handle)

!print*, 'avtac state_handle%copies(:,:)',state_handle%copies(:,:)
!print*, 'avtac state_handle%vars  (:,:)',state_handle%vars(:,:)


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
owner = get_owner(i, my_task_id(group_comm))
target_disp = (i-1)/(group_size)

if (debug) then
   print*, 'my_task_id() = ', my_task_id(), ' i = ', i, 'owner_group = ', owner, 'owner_world', get_owner(i,my_task_id()),' target_disp', target_disp
endif

! grab the info
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, my_window, ierr)
call mpi_get(get_my_val, 1, datasize, owner, target_disp, 1, datasize, my_window, ierr)
call mpi_win_unlock(owner, my_window, ierr)

end function

!-------------------------------------------------------------
!> create groups for grid
!> make a communicator for the distributed grid
subroutine local_create_groups

allocate(group_members(group_size)) ! this is module global

call mpi_comm_group(  mpi_comm_world, group_comm,                           ierr ) ! get the word group from mpi_comm_world
call build_my_group(  my_task_id(),   group_size, group_members )                 ! create a list of processors in the grid group
call mpi_group_incl(  group_comm,     group_size, group_members, sub_group, ierr )
call mpi_comm_create( mpi_comm_world, sub_group,   group_comm,              ierr )
call mpi_comm_rank(   group_comm,     group_rank,                           ierr ) ! rank within group

if (debug) then
   print*, 'my_task_id(), group_rank, group_size, sub_group', my_task_id(), group_rank, group_size, sub_group
endif
end subroutine local_create_groups

!-----------------------------------------------------------
!> build the group to store the grid
subroutine build_my_group(myrank, group_size, group_members)

implicit none

integer, intent(in)     :: myrank ! why are you passing this in?
integer, intent(inout)  :: group_size ! need to modify this if your #tasks does not divide by group size
integer, intent(out)    :: group_members(group_size)

integer bottom, top !< start and end members of the group
integer i

! integer arithmatic. rouding down to the lowest group size
bottom = (myrank / group_size ) * group_size
top = bottom + group_size - 1
if (top >= task_count()) then
   top = task_count() - 1
   group_size = top - bottom + 1
   print*, 'rank', myrank, 'bottom top', bottom, top, 'group_size', group_size
endif


! fill up group members
group_members(1) = bottom
do i = 2, group_size
   group_members(i) = group_members(i-1) + 1
enddo

end subroutine build_my_group

!---------------------------------------------------------
!> Get the group number
function get_owner(i, pe)
integer, intent(in) :: i
integer, intent(in) :: pe
integer :: get_owner

integer :: my_group
integer :: num_groups
integer :: num_vars

num_groups = NX / group_size
num_vars   = NX / num_groups

my_group = (pe/group_size) * group_size
get_owner  = my_group+mod(i-1,group_size)

end function get_owner

!----------------------------------------------------------------------
!> initialize modules that need it
subroutine initialize_modules_used()

!print*, 'initialize_mpi_utilities'
call initialize_mpi_utilities('simple_test')

!print*, 'register_module'
call register_module(source,revision,revdate)

end subroutine initialize_modules_used

!----------------------------------------------------------------------
!> clean up before exiting
subroutine finalize_modules_used()

! this must be last, and you can't print/write anything
! after this is called.
call finalize_mpi_utilities()

end subroutine finalize_modules_used

end program simple_test

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
