!****************************************************************************
! FILE: mpi_mm.f
! DESCRIPTION:
!   MPI Matrix Multiply - Fortran Version
!   In this code, the master task distributes a matrix multiply
!   operation to numtasks-1 worker tasks.
!   NOTE1:  C and Fortran versions of this code differ because of the
!   way
!   arrays are stored/passed.  C arrays are row-major order but Fortran
!   arrays are column-major order.
! AUTHOR: Blaise Barney. Adapted from Ros Leibensperger, Cornell Theory
!   Center. Converted to MPI: George L. Gusciora, MHPCC (1/95)
! LAST REVISED: 04/02/05
!****************************************************************************

program test_groups

use mpi

integer :: numtasks, taskid, group_rank, group_size
integer :: ierr
integer :: color
integer :: MPI_GROUP_COMM

call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )

color = taskid / 2

call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, taskid, MPI_GROUP_COMM, ierr)

call MPI_COMM_RANK(MPI_GROUP_COMM, group_rank, ierr)
call MPI_COMM_SIZE(MPI_GROUP_COMM, group_size, ierr)

print *, 'world rank/size:',taskid,'/',numtasks, 'row rank/size', group_rank, group_size
      

call MPI_FINALIZE(ierr)

end
