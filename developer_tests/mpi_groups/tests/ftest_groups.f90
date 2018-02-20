! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ftest_mpi

! simple MPI fortran program.  use to test running interactively
! with MPI parallel communication libraries.  warning -- this program
! may compile without obvious errors, but at runtime, unless MPI_Init()
! returns 0 as the error code, there is a good chance the compile and
! link phase did not succeed.

! The following 2 build tips are the 2 places where different installations
! of MPI seem to vary the most.  Some systems have an include file, some
! have a F90 module.  Some require an interface block to use the system()
! function, some give an error if it is here.   You can use this program
! to figure out which combinations work on your system.  Then go into the
! $DART/mpi_utilities and make the same two changes in mpi_utilities_mod.f90,
! and just the system() change (if needed) in null_mpi_utilities_mod.f90.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 1:
! Most fortran MPI implementations provide either a fortran 90 module
! which defines the interfaces to the MPI library routines, or an include
! file which defines constants.  Try to use the module if it is available.

use mpi

implicit none

!include "mpif.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUILD TIP 2:
! Some systems require this interface block in order to use the system()
! function.  However, some other systems complain if this is here... 
! If this is a problem your program will not link and most likely give 
! you an error about an undefined symbol (something like '_system_').  
! Comment this block in or out as needed.

 ! interface block for getting return code back from system() routine
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface
 ! end block

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! integer variables
integer :: ierror, world_rank, world_size, group_rank, group_size
integer :: i, j
integer :: static_id

logical :: do_print

integer, parameter :: nxA = 10
integer, parameter :: nyA = 10
real, allocatable  :: A(:,:)

integer(KIND=MPI_OFFSET_KIND) :: x_length, num_static_arrays
integer(KIND=MPI_OFFSET_KIND) :: my_num_y_vals !< splitting up by we

! real :: B(nxA,nyA) = -1
! real :: C(nxA,nyA) = -1

logical :: not_initialized = .true.

allocate(A(nxA, nyA))
A = -1.0

group_size = 4

ierror = -999
call MPI_Init(ierror)

call get_world_rank(world_rank)
do_print = (world_rank == 0)

if (world_rank == 0) then
   do i = 1,nxA
     do j = 1,nyA
        A(i,j) = i + (j-1)*nxA
     enddo
   enddo
endif

call split_groups()


static_id = distribute_static_data(nxA, nyA, A) 

ierror = -999
call MPI_Finalize(ierror)
if (ierror /= MPI_SUCCESS) then
   print *, "MPI_Finalize() did not succeed, error code = ", ierror
   stop
endif

if (do_print) then
   print *, "All MPI calls succeeded, test passed."
   print *, "program end"
endif

contains

!-------------------------------------------------------------------------------

subroutine initialize_mpi()

if (not_initialized) then
   call MPI_Init(ierror)
   not_initialized = .false.
endif

end subroutine initialize_mpi

!-------------------------------------------------------------------------------

subroutine split_groups()

integer :: row_rank
integer :: row_size
integer :: MPI_ROW_COMM
integer :: color
integer :: world_rank
integer :: world_size

if( .not. not_initialized ) call initialize_mpi()

call get_world_rank(world_rank)
call get_world_size(world_size)

! Testing Split: NOTE can be expensive for large number of processors
color = world_rank / 4
call MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, MPI_ROW_COMM, ierror)
if (ierror /= MPI_SUCCESS) then
   print*, "MPI_Comm_split failed with error code", ierror
   stop
endif

call MPI_Comm_rank(MPI_ROW_COMM, row_rank, ierror)
call MPI_Comm_size(MPI_ROW_COMM, row_size, ierror)

do i = 0, world_size-1
   call MPI_Barrier(MPI_COMM_WORLD, ierror)
   call MPI_Barrier(MPI_ROW_COMM, ierror)
   if(world_rank == i) then
      write(*,'(''WORLD RANK/SIZE:'',I4,''/'',I4,'' ROW RANK/SIZE:'',I4,''/'',I4)') world_rank, world_size, row_rank, row_size 
   endif
enddo

call MPI_Comm_free(MPI_ROW_COMM, ierror)
if (ierror /= MPI_SUCCESS) then
   print*, "MPI_Comm_free failed with error code", ierror
   stop
endif

end subroutine split_groups

!-------------------------------------------------------------------------------

function distribute_static_data(NX, NY, ARRAY) result(static_id)
integer, intent(in) :: NX, NY
real,    intent(in) :: ARRAY(:,:)

integer :: static_id

integer :: MPI_COMM_GRID
integer :: bottom, top
integer :: subgroup
integer, allocatable :: members(:)

! Testing Groups
call MPI_Comm_group(MPI_COMM_WORLD, MPI_COMM_GRID, ierror)
if (ierror /= MPI_SUCCESS) then
   print*, "MPI_Comm_group failed with error code", ierror
   stop
endif

bottom = (world_rank/group_size)*group_size
top    = bottom + group_size - 1

if (top >= world_size) then
   top = world_size - 1
   group_size = top - bottom + 1
endif

allocate(members(group_size))

call MPI_Barrier(MPI_COMM_WORLD, ierror)

do i = 0, world_size-1
   call MPI_Barrier(MPI_COMM_WORLD, ierror)
   if(world_rank == i) then
      write(*,'(''WORLD RANK'',I4,'' bottom - top :'',2I4,'' group_size'', I4)') world_rank, bottom, top, group_size 
   endif
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierror)

members(1) = bottom
do i = 2, group_size
   members(i) = members(i-1) + 1
enddo

do i = 0, world_size-1
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
      if(world_rank == i) then
        do j = 1, group_size
           print*, 'rank', world_rank, 'members(j)', members(j), j 
        enddo
      endif
enddo

call MPI_Barrier(MPI_COMM_WORLD, ierror)

call MPI_Group_incl(MPI_COMM_GRID, group_size, members, subgroup, ierror)
call MPI_Comm_create(MPI_COMM_WORLD, subgroup, MPI_COMM_GRID, ierror)
call MPI_Comm_rank(MPI_COMM_GRID, group_rank, ierror) ! rank within group


call print_array('A1 rank 0', ARRAY, NX, NY, 0)
call MPI_Barrier(MPI_COMM_WORLD, ierror)

call print_array('A1 rank 1', ARRAY, NX, NY, 1)
call MPI_Barrier(MPI_COMM_WORLD, ierror)

call MPI_Group_rank(subgroup, group_size, ierror)
call MPI_Group_size(subgroup, group_rank, ierror)

call MPI_Barrier(MPI_COMM_WORLD, ierror)
do i = 1, (world_size-1)
   call MPI_Barrier(MPI_COMM_WORLD, ierror)
   if(world_rank == i) then
      write(*,'(''GROUP RANK/SIZE:'',I4,''/'',I4,'' GROUP RANK/SIZE:'',I4,''/'',I4)') world_rank, world_size, group_size, group_rank
   endif
enddo

deallocate(members)

static_id = MPI_COMM_GRID

end function distribute_static_data

!-------------------------------------------------------------------------------

subroutine create_window()

integer ierr, sizedouble, datasize, my_window
integer(KIND=MPI_ADDRESS_KIND) :: window_size
!real, allocatable :: global_static_data(:) !< local static data 
real :: global_static_data(*)
pointer(aa, global_static_data)

if( .not. not_initialized ) call initialize_mpi()

datasize = MPI_REAL ! this should come from MPI_UTILITIES

call mpi_type_size(datasize, sizedouble, ierr) ! datasize comes from mpi_utilities_mod
window_size = num_static_arrays*x_length*my_num_y_vals*sizedouble

call mpi_win_create(global_static_data,  window_size, sizedouble, MPI_INFO_NULL, mpi_comm_world, my_window, ierr)

end subroutine create_window

!-------------------------------------------------------------------------------

subroutine get_world_size(world_size)
integer, intent(out) :: world_size

ierror = -999

if( .not. not_initialized ) call initialize_mpi()

world_size = -1
call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierror)
if (ierror /= MPI_SUCCESS) then
   print*, "MPI_Comm_size failed with error code", ierror
   stop
endif

if (do_print) print *, "Total MPI tasks: ", world_size

end subroutine get_world_size
!-------------------------------------------------------------------------------

subroutine get_world_rank(world_rank)
integer, intent(out) :: world_rank

ierror = -999

if( .not. not_initialized ) call initialize_mpi()

world_size = -1
call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierror)
if (ierror /= MPI_SUCCESS) then
   print*, "MPI_Comm_rank failed with error code", ierror
   stop
endif

end subroutine get_world_rank

!-------------------------------------------------------------------------------

subroutine print_array(name_array,A,nx,ny,rank)

character(len=*), intent(in) :: name_array
real,             intent(in) :: A(:,:)
integer,          intent(in) :: nx, ny
integer,          intent(in) :: rank

integer :: i
integer :: world_rank

call get_world_rank(world_rank)

if(rank == world_rank) then
   write(*,'(A,I2,A,I2,A)') trim(name_array)//'(1:',nx, ',1:',ny,') = '
   do i = 1,nx
     do j = 1,ny
        write(*,'(I4)',advance="no") int(A(i,j))
     enddo
     write(*,*) ''
   enddo
endif
end subroutine print_array

!-------------------------------------------------------------------------------

end program ftest_mpi


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
