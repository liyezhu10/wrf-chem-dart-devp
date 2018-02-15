! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_distributed_static_data

use types_mod,                   only : r8

use dart_pop_mod,                only : get_horiz_grid_dims,      &
                                        read_horiz_grid,          &
                                        read_topography

use distributed_static_data_mod, only : distribute_static_data,   &
                                        initialize_static_data_space, &
                                        get_static_data,          &
                                        free_window, create_window

use utilities_mod,               only : register_module,          &
                                        error_handler,            &
                                        E_ERR, E_MSG, E_ALLMSG,   &
                                        find_namelist_in_file,    &
                                        check_namelist_read,      &
                                        do_nml_file, do_nml_term, &
                                        nmlfileunit, do_output

use mpi_utilities_mod,           only : initialize_mpi_utilities, &
                                        finalize_mpi_utilities,   &
                                        task_sync, my_task_id,    &
                                        task_count

use mpi

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! message string for error handler
! character(len=128) :: msgstring

! io file variables
integer io, iunit

! indexing variables
integer i, j

! grid dimensions
integer :: NX, NY

 ! These arrays store the longitude and latitude of the lower left corner of
 ! each of the dipole u quadrilaterals and t quadrilaterals.
real(r8), allocatable :: ULAT(:,:), ULON(:,:), TLAT(:,:), TLON(:,:)
 
 ! integer, lowest valid cell number in the vertical
integer, allocatable  :: KMT(:, :), KMU(:, :)

integer  :: ID_ULAT, ID_ULON, ID_TLAT, ID_TLON

real(r8) :: my_val_test
real(r8) :: t1, t2

! initialize the dart libs
call initialize_module()

call get_horiz_grid_dims(NX, NY)

if (do_output()) then
   write(*,*) ''
   write(*,'(''Grid Dimensions NX = '',I4,'', NY = '',I4)') NX, NY
   write(*,*) ''
endif

call initialize_static_data_space(4,NX,NY) 
call create_window()

allocate(ULAT(NX,NY), ULON(NX,NY), TLAT(NX,NY), TLON(NX,NY))
allocate( KMT(NX,NY),  KMU(NX,NY))

! Fill in the grid information 
! horiz grid initializes ULAT/LON, TLAT/LON as well.
call read_horiz_grid(NX, NY, ULAT, ULON, TLAT, TLON)
call read_topography(NX, NY,  KMT,  KMU)

ID_ULAT = distribute_static_data(ULAT)
ID_ULON = distribute_static_data(ULON)
ID_TLAT = distribute_static_data(TLAT)
ID_TLON = distribute_static_data(TLON)


t1 = MPI_Wtime()
!if(my_task_id() == 0) then
   do j = 1, NY
      do i = 1, NX
         !print*, 'getting i,j : ', i,j
         my_val_test = get_static_data(ID_ULAT,i,j)
         !if(ULAT(i,j) /= my_val_test) then
         !   print*, 'error at i,j ', i, j, ULAT(i,j), my_val_test
         !endif
         my_val_test = get_static_data(ID_ULON,i,j)
         !if(ULON(i,j) /= my_val_test) then
         !   print*, 'error at i,j ', i, j, ULON(i,j), my_val_test
         !endif
         my_val_test = get_static_data(ID_TLON,i,j)
         my_val_test = get_static_data(ID_TLAT,i,j)
      enddo
   enddo
!endif
t2 = MPI_Wtime()

! if(my_task_id() == 0) then
   write(*,*) 'time elapsed get_static_data : ', t2-t1, ' PE ', my_task_id()
! endif

call task_sync()

t1 = MPI_Wtime()
! if(my_task_id() == 0) then
   do j = 1, NY
      do i = 1, NX
         my_val_test = ULAT(i,j)
         my_val_test = ULON(i,j)
         my_val_test = TLAT(i,j)
         my_val_test = TLON(i,j)
      enddo
   enddo
! endif
t2 = MPI_Wtime()

! if(my_task_id() == 0) then
   write(*,*) 'time elapsed directly        : ', t2-t1, ' PE ', my_task_id()
! endif

call task_sync()

t1 = MPI_Wtime()
   my_val_test = get_static_data(ID_ULAT,120,140)
t2 = MPI_Wtime()
write(*,*) 'TIME elapsed get_static_data : ', t2-t1, ' PE ', my_task_id()

t1 = MPI_Wtime()
   my_val_test = ULAT(120,140)
t2 = MPI_Wtime()
write(*,*) 'TIME elapsed directly        : ', t2-t1, ' PE ', my_task_id()

deallocate(ULAT, ULON, TLAT, TLON, KMT, KMU)
   
call free_window()
call task_sync()

! finalize test_distributed_static_data
call error_handler(E_MSG,'test_distributed_static_data','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

! end of main code

contains

!----------------------------------------------------------------------

subroutine initialize_module

  call initialize_mpi_utilities('test_distributed_static_data')
  call register_module(source, revision, revdate)

end subroutine initialize_module

!----------------------------------------------------------------------

subroutine print_array(name_array,A,nx,ny)

character(len=*), intent(in) :: name_array
real(r8),         intent(in) :: A(:,:)
integer,          intent(in) :: nx, ny

integer :: i

   write(*,'(A,I2,A,I2,A)') trim(name_array)//'(1:',nx, ',1:',ny,') = '
   do i = 1,nx
     do j = 1,ny
        write(*,'(I4)',advance="no") int(A(i,j))
     enddo
     write(*,*) ''
   enddo

end subroutine print_array

!----------------------------------------------------------------------

end program

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
