!> Aim: to replace the array phb, with a distributed version
!> Current thinking: Only distribute the array as much as necessary.
!> Ideally, if memory was unliimited, every task would have the whole array.
!>
!> Splitting in x only, because this is simple, and I can't think
!> of a reason to split in more than one dimension yet.

module distributed_static_data_mod

use mpi_utilities_mod,     only : datasize, my_task_id, task_count, send_to, &
                                  receive_from , task_sync

use types_mod,             only : r8

use mpi

implicit none
private

public :: create_window, free_window

!!! only public at the moment for testing
public :: create_groups, build_my_group, initialize_static_data_space, &
          distribute_static_data, get_static_data

! grid window
real(r8), allocatable :: global_static_data(:) !< local static data 

integer  :: sd_window

integer :: group_size = 4 !< should be namelist option
integer :: group_handle   !< mpi_comm_world group
integer :: subgroup       !< subgroup for the grid
integer :: mpi_comm_grid
integer :: local_rank     !< rank within group
integer, allocatable :: group_members(:)

! pntecdf variables
integer(KIND=MPI_OFFSET_KIND) :: x_start, x_length, y_length, num_static_arrays
integer(KIND=MPI_OFFSET_KIND) :: start(4)
integer(KIND=MPI_OFFSET_KIND) :: count(4)
integer(KIND=MPI_OFFSET_KIND) :: stride(4)
integer(KIND=MPI_OFFSET_KIND) :: my_num_x_vals !< splitting up by we

integer :: array_number = 0 ! starting index for array

contains

!-------------------------------------------------------------
!> create groups for grid
!> make a communicator for the distributed grid
subroutine create_groups(grp_size)

integer, intent(in) :: grp_size

integer ierr ! all MPI errors are fatal anyway

group_size = grp_size

allocate(group_members(group_size)) ! this is module global

call mpi_comm_group(mpi_comm_world, group_handle, ierr)  ! get group handle from mpi_comm_world
call build_my_group(group_size, group_members) ! create a list of processors in the grid group
call mpi_group_incl(group_handle, group_size, group_members, subgroup, ierr)
call mpi_comm_create(mpi_comm_world, subgroup, mpi_comm_grid, ierr)
call mpi_comm_rank(mpi_comm_grid, local_rank, ierr) ! rank within group

end subroutine create_groups

!-----------------------------------------------------------
!> build the group to store the grid
subroutine build_my_group(group_size, group_members)

! need to modify this if your #tasks does not divide by group size
integer, intent(inout) :: group_size 
integer, intent(out)   :: group_members(group_size)

integer :: bottom, top !< start and end members of the group
integer :: i 
integer :: myrank

myrank = my_task_id()

! integer arithmatic. rouding down to the lowest group size
bottom = (myrank / group_size ) * group_size
top = bottom + group_size - 1
if (top >= task_count()) then
   top = task_count() - 1
   group_size = top - bottom + 1
   ! print*, 'rank', myrank, 'bottom top', bottom, top, 'group_size', group_size
endif


! fill up group members
group_members(1) = bottom
do i = 2, group_size
   group_members(i) = group_members(i-1) + 1
   !print*, 'group_members(i) ', myrank, i, group_members(i) 
enddo

if (mod(myrank,group_size) == 0) then
  do i = 1, group_size
     !write(*,'(''rank['',I2,''] group_members('',I2,'')'',2I3)') myrank, i, group_members(i)
  enddo
endif

end subroutine build_my_group

!-----------------------------------------------------------
!> create window for phb
!> Need to add a loop for domain.  Just single domain for now.
!> Window is filled with columns contiguous in memory, possibly
!> we will grab a whole column at once.
!> 
!> If you want to change the order in the window, be sure to 
!> change the subroutine who_has_grid_info
subroutine create_window

integer ierr, sizedouble 
integer(KIND=MPI_ADDRESS_KIND) :: window_size

call mpi_type_size(datasize, sizedouble, ierr) ! datasize comes from mpi_utilities_mod
window_size = num_static_arrays*my_num_x_vals*y_length*sizedouble

call mpi_win_create(global_static_data,  window_size, sizedouble, MPI_INFO_NULL, mpi_comm_grid, sd_window, ierr)

end subroutine create_window


!---------------------------------------------------------
!> Free the mpi windows you created
subroutine free_window
integer :: ierr

call mpi_win_free(sd_window, ierr)

deallocate(group_members)

end subroutine free_window

!---------------------------------------------------------
!> Function to get static data
function get_static_data(Static_ID, i, j) result(val)

real(r8)            :: val       !< static data
integer, intent(in) :: i, j      !< 2D location
integer, intent(in) :: Static_ID

integer                          :: owner       !< which task has the part of static data we need
integer(KIND=MPI_ADDRESS_KIND)   :: target_disp !< displacement
integer                          :: ierr

! caluclate who has the info
call who_has_grid_info(Static_ID, i, j, owner, target_disp)

! grab the info
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, sd_window, ierr)
call mpi_get(val, 1, datasize, owner, target_disp, 1, datasize, sd_window, ierr)
call mpi_win_unlock(owner, sd_window, ierr)

end function get_static_data

!---------------------------------------------------------
!> get the owner and location in memory of the grid point
subroutine who_has_grid_info(Static_ID, i, j, owner, target_disp)

integer,                        intent(in)  :: Static_ID
integer,                        intent(in)  :: i
integer,                        intent(in)  :: j
integer,                        intent(out) :: owner
integer(KIND=MPI_ADDRESS_KIND), intent(out) :: target_disp !< displacement

integer ::x_stride, x_local

x_stride = x_length/group_size
owner = (i-1)/x_stride

if(owner == group_size - 1) then
   x_stride = x_length - x_stride*(group_size-1)
endif

x_local = i - x_stride*owner
target_disp = (static_id-1)*x_stride*y_length + (j-1)*x_stride + x_local

end subroutine who_has_grid_info

!---------------------------------------------------------
!> Initalize local static data for each processor
!> @todo FIXME : Need to make this more generic in the for models
!>               with static data that has different shapes.
subroutine initialize_static_data_space(num_arrays, NX, NY)

integer, intent(in) :: num_arrays
integer, intent(in) :: NX
integer, intent(in) :: NY

integer :: x_mod

! set global variables
num_static_arrays = num_arrays
x_length          = NX
y_length          = NY

! split up in x dimension, could also split in y?
my_num_x_vals = x_length / group_size
x_start       = local_rank*my_num_x_vals
x_mod         = mod(x_length,group_size)
if (x_mod /= 0 .and. local_rank == group_size - 1) then
   my_num_x_vals = my_num_x_vals + x_mod
   !write(*,'(''last  rank['',I2,''] my_num_x_vals * y_length = '',I4,'' * '',I4,'' = '',I7)') &
   !       my_task_id(), my_num_x_vals, y_length, my_num_x_vals*y_length
endif

allocate(global_static_data(num_static_arrays*my_num_x_vals*y_length)) ! local static data

end subroutine initialize_static_data_space

!---------------------------------------------------------
!> Distribute the static data to the appropriate processors
function distribute_static_data(local_static_data) result(static_id)

real(r8), intent(inout) :: local_static_data(x_length,y_length)
integer :: static_id

integer :: iy, g_istart, g_iend, l_istart, l_iend

! block data
static_id = array_number + 1
g_istart = (static_id-1)*my_num_x_vals*y_length + 1
g_iend   = g_istart+my_num_x_vals-1
l_istart = x_start+1
l_iend   = l_istart+my_num_x_vals-1

do iy=1,y_length
   !print*, "task_id - ", my_task_id(), ": g_istart ", g_istart
   !print*, "task_id - ", my_task_id(), ": g_iend   ", g_iend
   !print*, "task_id - ", my_task_id(), ": l_istart ", l_istart
   !print*, "task_id - ", my_task_id(), ": l_iend   ", l_iend
   !print*, "task_id - ", my_task_id(), ": iy       ", iy
   !print*, "task_id - ", my_task_id(), ": size     ", size(global_static_data)
   global_static_data(g_istart:g_iend) = local_static_data(l_istart:l_iend,iy) 
   g_istart = g_iend + 1
   g_iend   = g_istart+my_num_x_vals-1
enddo

!print*, "task_id - ", my_task_id(), ": flat-array  ", global_static_data
!do iy=1,my_num_x_vals*y_length
!   write(*,'(A,I3,I3)',advance="no") "task_id - ", my_task_id(), int(global_static_data(iy))
!   write(*,*) ""
!enddo



end function distribute_static_data

end module distributed_static_data_mod

