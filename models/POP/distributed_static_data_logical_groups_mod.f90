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
public :: create_groups, initialize_static_data_space, &
          distribute_static_data, get_static_data

! grid window
!real(r8), allocatable :: global_static_data(:) !< local static data 
real(r8) :: global_static_data(*)
pointer(aa, global_static_data)

integer  :: sd_window

integer :: group_handle   !< mpi_comm_world group
integer :: subgroup       !< subgroup for the grid
integer :: subgroup_size  !< subgroup for the grid
integer :: mpi_comm_grid
integer :: local_rank     !< rank within group
integer, allocatable :: y_starts(:)
integer, allocatable :: y_ends(:)

! pntecdf variables
integer(KIND=MPI_OFFSET_KIND) :: x_length, y_length, num_static_arrays
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
integer mod_tasks, num_subgroups


subgroup_size = group_size
num_subgroups = task_count()/group_size
subgroup = my_task_id()/group_size
mod_tasks = mod(task_count(),group_size)
local_rank = mod(my_task_id(),group_size)

if(mod_tasks /= 0) then
   if(num_subgroups*group_size-1 < my_task_id()) then
      subgroup_size = mod_tasks
   endif
endif

allocate(y_starts(subgroup_size+1))    ! this is module global
allocate(y_ends(subgroup_size+1))      ! this is module global

end subroutine create_groups

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
window_size = num_static_arrays*x_length*my_num_y_vals*sizedouble

call mpi_win_create(global_static_data,  window_size, sizedouble, MPI_INFO_NULL, mpi_comm_world, sd_window, ierr)

end subroutine create_window


!---------------------------------------------------------
!> Free the mpi windows you created
subroutine free_window
integer :: ierr

call mpi_win_free(sd_window, ierr)

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
ierr = 0
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

integer :: io, local_owner, local_y_num_vals

do io=1,subgroup_size
   if(y_starts(io) <= j .and. j < y_starts(io+1)) local_owner = io-1
enddo
   
owner = local_owner + subgroup*group_size

local_y_num_vals = y_starts(local_owner+2) - y_starts(local_owner+1)

target_disp = x_length*((static_id-1)*local_y_num_vals + j-y_starts(local_owner+1)) + i - 1 


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
integer io, iunit, ierr, sizedouble, my_mod_y_vals

integer(KIND=MPI_ADDRESS_KIND) :: window_size

integer :: i, y_mod

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

my_num_y_vals = y_length / subgroup_size
my_mod_y_vals = mod(y_length,subgroup_size)

y_starts(1) = 1;
do i=1,subgroup_size-1
   y_starts(i+1) = y_starts(i) + my_num_y_vals
   if (my_mod_y_vals /= 0 .and. i-1 < my_mod_y_vals) then
      y_starts(i+1) = y_starts(i+1) + 1
   endif
enddo 
y_starts(subgroup_size+1) = y_length+1

my_num_y_vals = y_starts(local_rank+2) - y_starts(local_rank+1)

call mpi_type_size(datasize, sizedouble, ierr) ! datasize comes from mpi_utilities_mod
window_size = num_static_arrays*x_length*my_num_y_vals*sizedouble
call MPI_ALLOC_MEM(window_size, mpi_info_null, aa, ierr)

end subroutine initialize_static_data_space

!---------------------------------------------------------
!> Distribute the static data to the appropriate processors
function distribute_static_data(local_static_data) result(static_id)

real(r8), intent(inout) :: local_static_data(x_length,y_length)
integer :: static_id

integer :: iy, iy_start, iy_end

! block data
static_id = array_number + 1
array_number = array_number +1

iy_start = x_length*((static_id-1)*my_num_y_vals) + 1
iy_end   = iy_start + x_length - 1

do iy=y_starts(local_rank+1),y_starts(local_rank+2)-1
   global_static_data(iy_start:iy_end) = local_static_data(:,iy) 
   iy_start = iy_end + 1
   iy_end   = iy_start + x_length - 1
enddo

call create_window()

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
