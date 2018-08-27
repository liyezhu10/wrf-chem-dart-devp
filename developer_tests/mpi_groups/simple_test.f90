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
                                  task_count, my_task_id, task_sync, set_group_size,&
                                  get_group_size, set_group_size,    &
                                  group_task_id, get_dart_mpi_comm, get_group_comm, &
                                  get_group_id, datasize

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

namelist /simple_test_nml/ NX, group_size, max_iter, debug

! io variables
integer                   :: iunit, io
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.

! model state variables
type(ensemble_type)   :: ens_handle
type(ensemble_type)   :: mean_handle

type(time_type)       :: model_time
integer(i8)           :: model_size, my_index

! misc. variables
integer :: i, idom, imem, domid, num_domains

! message strings
character(len=512) :: my_base, my_desc
character(len=512) :: string1

character(len=256), allocatable :: file_array_input(:,:)
character(len=256), dimension(1) :: var_names = (/'temp'/)
integer,parameter :: one_domain = 1

real(r8) ::  u, my_val

integer, allocatable :: group_members(:)
integer :: group_all
integer :: subgroup
integer :: mpi_comm_grid
integer :: local_rank

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
call create_groups()
t2 = MPI_WTIME() - t1

call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

if (my_task_id() == 0) then
   print*, 'group_size = ', group_size, ' NX = ', NX
   print*, 'create_groups : Max Time = ', max_time, t1, t2
   print*, 'create_groups : Min Time = ', min_time, t1, t2
endif

!----------------------------------------------------------------------
! create data array
!----------------------------------------------------------------------
allocate(my_array(NX/group_size))

do ii = 1, NX/group_size
   my_array(ii) = local_rank*NX/group_size + ii
end do

!----------------------------------------------------------------------
! create window
!----------------------------------------------------------------------
call create_window()

!----------------------------------------------------------------------
! timing test
!----------------------------------------------------------------------
t2 = 0.0_r8
do ii = 1, max_iter
   call random_number(u)
   jj = FLOOR(NX*u)+1
   t1 = MPI_WTIME()
   my_val = get_my_val(jj)
   t2 = t2 + (MPI_WTIME() - t1)
enddo

call MPI_REDUCE(t2, max_time, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
call MPI_REDUCE(t2, min_time, 1, MPI_REAL8, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
if (my_task_id() == 0) then
   print*, 'get_value     : Max Time = ', max_time, t1, t2
   print*, 'get_value     : Min Time = ', min_time, t1, t2
endif

!----------------------------------------------------------------------
! print array
!----------------------------------------------------------------------
if (debug) then
   do jj = 0, task_count()-1
      if (my_task_id() == jj) then
         print*, my_task_id(), '::', my_array(:)
      else
         call task_sync()
      endif
   enddo
endif

!----------------------------------------------------------------------
! debug get_my_val
!----------------------------------------------------------------------
if (debug) then
   do jj = 0, task_count()-1
      if (my_task_id() == jj) then
         do ii = 1, NX
            print*, 'my_task_id() = ', my_task_id(), 'get_owner(ii)', ii, get_owner(ii, my_task_id()), get_my_val(ii)
         enddo
      else
         call task_sync()
      endif
   enddo
   
   call task_sync()
endif

!----------------------------------------------------------------------
! finalize simple_test
!----------------------------------------------------------------------
call free_window()

call task_sync()

call finalize_modules_used()

!======================================================================
contains
!======================================================================

!-------------------------------------------------------------
!> create groups for grid
!> make a communicator for the distributed grid
subroutine create_groups

allocate(group_members(group_size)) ! this is module global


call mpi_comm_group(  mpi_comm_world, group_all,                           ierr ) ! get the word group from mpi_comm_world
call build_my_group(  my_task_id(),   group_size, group_members )                 ! create a list of processors in the grid group
call mpi_group_incl(  group_all,      group_size, group_members, subgroup, ierr )
call mpi_comm_create( mpi_comm_world, subgroup,   mpi_comm_grid,           ierr )
call mpi_comm_rank(   mpi_comm_grid,  local_rank,                          ierr ) ! rank within group

if (debug) then
   print*, 'my_task_id(), local_rank, group_size, subgroup', my_task_id(), local_rank, group_size, subgroup
endif
end subroutine create_groups

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

!-----------------------------------------------------------
!> create window for get_my_val
!> Need to add a loop for domain.  Just single domain for now.
!> Window is filled with columns contiguous in memory, possibly
!> we will grab a whole column at once.
!> 
!> If you want to change the order in the window, be sure to 
!> change the subroutine who_has_grid_info
subroutine create_window

integer ii, jj, kk, sizedouble, count
integer(KIND=MPI_ADDRESS_KIND) :: window_size

! datasize comes from mpi_utilities_mod
call mpi_type_size(datasize, sizedouble, ierr)

window_size = (NX/group_size)*sizedouble
aa          = malloc(NX)

call MPI_ALLOC_MEM(window_size, mpi_info_null, aa, ierr)

! can't do my_array assignment with a cray pointer, so you need to loop
!#! count = 1
!#! do ii = 1, NX
!#!    duplicate_array(count) = my_array(ii)
!#!    count = count + 1
!#! enddo

call mpi_win_create(my_array,        window_size,   &
                    sizedouble,      MPI_INFO_NULL, &
                    mpi_comm_grid,   my_window,     ierr)

end subroutine create_window

!---------------------------------------------------------
!> Free the mpi windows you created
subroutine free_window
integer :: ierr

call mpi_win_free(my_window, ierr)
call MPI_FREE_MEM(duplicate_array, ierr)

deallocate(group_members, my_array)

end subroutine free_window

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
call finalize_mpi_utilities()

end subroutine finalize_modules_used

end program simple_test

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
