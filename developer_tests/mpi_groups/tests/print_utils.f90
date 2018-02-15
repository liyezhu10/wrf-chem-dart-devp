module print_utils
   implicit none

   subroutine print_mpi_init_error(ierror)
   integer, intent(in) :: ierror
         print *, "MPI_Init() did not succeed, error code = ", ierror
         print *, "If error code is -999, the most likely problem is that"
         print *, "the right MPI libraries were not found at compile time."
   end subroutine print_mpi_init_error
   
   subroutine print_mpi_comm_rank_error(ierror)
   integer, intent(in) :: ierror
         print *, "MPI_Comm_rank() did not succeed, error code = ", ierror
   subroutine print_mpi_comm_rank_error
   
   subroutine print_mpi_comm_size_error(ierror)
   integer, intent(in) :: ierror
         print *, "MPI_Comm_size() did not succeed, error code = ", ierror
   end subroutine print_mpi_comm_size_error

end module print_utils
