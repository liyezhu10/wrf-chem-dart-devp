!> Aim: to replace the array phb, with a distributed version
!> Current thinking: Only distribute the array as much as necessary.
!> Ideally, if memory was unliimited, every task would have the whole array.
!>
!> Splitting in x only, because this is simple, and I can't think
!> of a reason to split in more than one dimension yet.

module distributed_static_data_mod

use utilities_mod,         only : register_module,          &
                                  error_handler,            &
                                  E_ERR, E_MSG, E_ALLMSG,   &
                                  find_namelist_in_file,    &
                                  check_namelist_read,      &
                                  do_nml_file, do_nml_term, &
                                  nmlfileunit, do_output

use mpi_utilities_mod,     only : datasize, my_task_id,     &
                                  task_count, send_to,      &
                                  receive_from , task_sync, &
                                  get_dart_mpi_comm

use types_mod,             only : r8

use mpi

implicit none
private

public :: create_window, free_window, collect_static_data_to_array

!!! only public at the moment for testing
public :: create_groups, build_my_group, initialize_static_data_space, &
          distribute_static_data, get_static_data

! grid window
!real(r8), allocatable :: global_static_data(:) !< local static data 
real(r8) :: global_static_data(*)
pointer(aa, global_static_data)

integer  :: sd_window

integer :: group_handle   !< mpi_comm_world group
integer :: subgroup       !< subgroup for the grid
integer :: mpi_comm_grid
integer :: local_rank     !< rank within group
integer, allocatable :: group_members(:)

! pntecdf variables
integer(KIND=MPI_OFFSET_KIND) :: y_start, x_length, y_length, num_static_arrays
integer(KIND=MPI_OFFSET_KIND) :: start(4)
integer(KIND=MPI_OFFSET_KIND) :: count(4)
integer(KIND=MPI_OFFSET_KIND) :: stride(4)
integer(KIND=MPI_OFFSET_KIND) :: my_num_y_vals !< splitting up by we

integer :: array_number = 0 ! starting index for array

!------------------------------------------------------------------
! namelist variables
!------------------------------------------------------------------
integer              :: group_size = 1

namelist /distributed_static_data_nml/  group_size

contains

!-------------------------------------------------------------
!> create groups for grid
!> make a communicator for the distributed grid
subroutine create_groups()

integer ierr ! all MPI errors are fatal anyway

allocate(group_members(group_size)) ! this is module global

call mpi_comm_group(get_dart_mpi_comm(), group_handle, ierr)  
call build_my_group(group_size, group_members) ! create a list of processors in the grid group
call mpi_group_incl(group_handle, group_size, group_members, subgroup, ierr)
call mpi_comm_create(get_dart_mpi_comm(), subgroup, mpi_comm_grid, ierr)
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
   !print*, 'rank', myrank, 'bottom top', bottom, top, 'group_size', group_size
endif
!print*, 'rank', myrank, 'bottom top', bottom, top, 'group_size', group_size
!write(*,*) ""

!call task_sync()

! fill up group members
group_members(1) = bottom
!print*, 'group_members(i) ', myrank, 1, group_members(1) 
!call task_sync()
do i = 2, group_size
   group_members(i) = group_members(i-1) + 1
   !print*, 'group_members(i) ', myrank, i, group_members(i) 
   !call task_sync()
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
window_size = num_static_arrays*my_num_y_vals*y_length*sizedouble

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
!print*, "id =", my_task_id()
!print*, "id =", my_task_id(), "calling who_has_grid"
call who_has_grid_info(Static_ID, i, j, owner, target_disp)
!print*, "rank =", my_task_id(), "owner =", owner, "i =", i, "j =", j

! grab the info
ierr = 0
!print*, "ierr =", ierr
call mpi_win_lock(MPI_LOCK_SHARED, owner, 0, sd_window, ierr)
!print*, " ******* after mip_win_lock: ierr =", ierr
call mpi_get(val, 1, datasize, owner, target_disp, 1, datasize, sd_window, ierr)
!print*, " ******* after mpi_get: ierr =", ierr
!print*, "val =", val
call mpi_win_unlock(owner, sd_window, ierr)
!print*, " ******* after mip_win_unlock "
!print*, "val =", val

end function get_static_data

!---------------------------------------------------------
!> get the owner and location in memory of the grid point
subroutine who_has_grid_info(Static_ID, i, j, owner, target_disp)

integer,                        intent(in)  :: Static_ID
integer,                        intent(in)  :: i
integer,                        intent(in)  :: j
integer,                        intent(out) :: owner
integer(KIND=MPI_ADDRESS_KIND), intent(out) :: target_disp !< displacement

integer ::y_stride, y_local

y_stride = y_length/group_size
owner = (j-1)/y_stride

if(owner >= group_size) owner = group_size-1
   
! The line below must come before the next if statement.
! We need the possibly smaller stride as starting point
y_local = j - y_stride*owner - 1

if(owner == group_size - 1) then
   y_stride = y_length - y_stride*(group_size-1)
endif

target_disp = (static_id-1)*x_length*y_stride + x_length*y_local + i - 1


end subroutine who_has_grid_info

!---------------------------------------------------------
!> Initalize local static data for each processor
!> @todo FIXME : Need to make this more generic in the for models
!>               with static data that has different shapes.
subroutine initialize_static_data_space(num_arrays, NX, NY)

integer, intent(in) :: num_arrays
integer, intent(in) :: NX
integer, intent(in) :: NY

! io file variables
integer io, iunit, ierr, sizedouble

integer(KIND=MPI_ADDRESS_KIND) :: window_size

integer :: y_mod

! Read the namelist entry
call find_namelist_in_file("input.nml", "distributed_static_data_nml", iunit)
read(iunit, nml = distributed_static_data_nml, iostat = io)
call check_namelist_read(iunit, io, "distributed_static_data_nml")

call create_groups()

! set global variables
num_static_arrays = num_arrays
x_length          = NX
y_length          = NY

! split up in x dimension, could also split in y?
my_num_y_vals = y_length / group_size
y_start       = local_rank*my_num_y_vals 
y_mod         = mod(y_length,group_size)
!write(*,'(''last  rank['',I2,'']: x_length , my_num_y_vals = '',I9,'' , '',I9,'': local_rank  = '',I4)') &
!     my_task_id(), x_length, my_num_y_vals, local_rank
!write(*,'(''last  rank['',I2,'']: x_length , my_num_y_vals = '',I9,'' , '',I9,'': y_start  = '',I4)') &
!     my_task_id(), x_length, my_num_y_vals, y_start
!write(*,'(''last  rank['',I2,''] my_num_y_vals * x_length = '',I4,'' * '',I4,'' = '',I7)') &
      !my_task_id(), my_num_y_vals, y_length, my_num_y_vals*x_length
if (y_mod /= 0 .and. local_rank == group_size - 1) then
   my_num_y_vals = my_num_y_vals + y_mod
   !write(*,'(''last  rank['',I2,''] my_num_y_vals * x_length = '',I4,'' * '',I4,'' = '',I7)') &
          !my_task_id(), my_num_y_vals, y_length, my_num_y_vals*x_length
endif
!write(*,'(''last  rank['',I2,''] my_num_y_vals * x_length = '',I4,'' * '',I4,'' = '',I7)') &
!      my_task_id(), my_num_y_vals, x_length, my_num_y_vals*x_length

!allocate(global_static_data(num_static_arrays*x_length*my_num_y_vals)) ! local static data
call mpi_type_size(datasize, sizedouble, ierr) ! datasize comes from mpi_utilities_mod
window_size = num_static_arrays*x_length*my_num_y_vals*sizedouble
aa = malloc(num_static_arrays*x_length*my_num_y_vals)
call MPI_ALLOC_MEM(window_size, mpi_info_null, aa, ierr)

end subroutine initialize_static_data_space

!---------------------------------------------------------
!> Distribute the static data to the appropriate processors
function distribute_static_data(local_static_data) result(static_id)

real(r8), intent(inout) :: local_static_data(x_length,y_length)
integer :: static_id

integer :: iy, istart, iend

! block data
static_id = array_number + 1
array_number = array_number +1

istart = (static_id-1)*x_length*my_num_y_vals + 1
iend   = istart + x_length - 1
!write(*,'(''last  rank['',I2,'']: x_length , my_num_y_vals = '',I9,'' , '',I9,'': y_start  = '',I4)') &
!     my_task_id(), x_length, my_num_y_vals, y_start

do iy=y_start+1,y_start+my_num_y_vals
   !write(*,'(''last  rank['',I2,''] istart , iend = '',I9,'' , '',I9,''iy  = '',I7)') &
   !      my_task_id(), istart, iend, iy
   global_static_data(istart:iend) = local_static_data(:,iy) 
   istart = iend + 1
   iend   = istart + x_length - 1
enddo

!print*, "task_id - ", my_task_id(), ": flat-array  ", global_static_data
!do iy=1,my_num_y_vals*y_length
!   write(*,'(A,I3,I3)',advance="no") "task_id - ", my_task_id(), int(global_static_data(iy))
!   write(*,*) ""
!enddo



end function distribute_static_data

!---------------------------------------------------------
!> collect all the data back into an array -- just a simple hack for now
subroutine collect_static_data_to_array(Static_ID, collected_array)

integer,  intent(in)  :: Static_ID
real(r8), intent(inout) :: collected_array(x_length,y_length)

integer :: i, j

do j=1,y_length
   do i=1,x_length
      collected_array(i,j) = get_static_data(Static_ID, i, j)
   enddo
enddo


end subroutine collect_static_data_to_array

end module distributed_static_data_mod

