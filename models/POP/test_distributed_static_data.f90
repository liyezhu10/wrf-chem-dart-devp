! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: test_distributed_static_data.f90 6256 2013-06-12 16:19:10Z thoar $

program test_distributed_static_data

use types_mod,                   only : r8

use dart_pop_mod,                only : get_horiz_grid_dims,      &
                                        read_horiz_grid,          &
                                        read_topography

use distributed_static_data_mod, only : build_my_group,           &
                                        create_groups,            &
                                        distribute_static_data,   &
                                        initialize_static_data_space, &
                                        get_static_data,          &
                                        create_window

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

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_par_msg/utilities/test_distributed_static_data.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6256 $"
character(len=128), parameter :: revdate  = "$Date: 2013-06-12 10:19:10 -0600 (Wed, 12 Jun 2013) $"

! message string for error handler
! character(len=128) :: msgstring

! io file variables
integer io, iunit

! indexing variables
integer i, j, t, numt

! grid dimensions
integer :: NX, NY

! ! These arrays store the longitude and latitude of the lower left corner of
! ! each of the dipole u quadrilaterals and t quadrilaterals.
! real(r8), allocatable :: ULAT(:,:), ULON(:,:), TLAT(:,:), TLON(:,:)
! 
! ! integer, lowest valid cell number in the vertical
! integer, allocatable  :: KMT(:, :), KMU(:, :)

integer  :: ID
integer, parameter :: nxA = 17
integer, parameter :: nyA = 19
real(r8) :: A(nxA,nyA) = -1

! ! group variables
! integer, allocatable :: group_members(:)

! real(r8) :: my_val_test


!------------------------------------------------------------------
! namelist variables
!------------------------------------------------------------------
integer              :: group_size = 1

namelist /test_distributed_static_data_nml/  group_size

! main code here
 
! initialize the dart libs
call initialize_module()

! write(*,*) ''
! write(*,'(''my_task_id()''I3)') my_task_id()

call get_horiz_grid_dims(NX,NY)


!call initialize_static_data_space(6,NX,NY)
call initialize_static_data_space(1,nxA,nyA)

!if(my_task_id() == 0) then
  do i = 1,nxA
    do j = 1,nyA
       A(i,j) = i + (j-1)*nxA
    enddo
  enddo
  if(my_task_id()==0) call print_array('A0', A, nxA, nyA)
!endif

call task_sync()
print*, ''

!!if(my_task_id() == 1) then
!!   call print_array('A1', A, nxA, nyA)
!!endif
! 
ID = distribute_static_data(A)
call task_sync()

numt = task_count()-1
do t=0,numt
   if(my_task_id() == t) then
     write(*,'(A,I2)') "my_task_id = ", my_task_id()
     do i = 1,nxA
       do j = 1,nyA
         !my_val_test = get_static_data(ID,i,j)
         !write(*,'(''(''i2,'','',i2,'')'',I3)',advance="no") i,j,int(get_static_data(ID,i,j))
         !print*, "my_task_id =", my_task_id()
         write(*,'(I4)',advance="no") int(get_static_data(ID,i,j))
       enddo
       write(*,*) ""
     enddo
   endif
   call task_sync()
enddo

!print*, ''
!
!!if(my_task_id() == 1) then
!!   call print_array('A1', A, nxA, nyA)
!!endif
! 
!ID = distribute_static_data(A)
!print*, "done with distributing"
!
!allocate(ULAT(NX,NY), ULON(NX,NY), TLAT(NX,NY), TLON(NX,NY))
!allocate( KMT(NX,NY),  KMU(NX,NY))
!
!
!! Fill in the grid information 
!! horiz grid initializes ULAT/LON, TLAT/LON as well.
!call read_horiz_grid(NX, NY, ULAT, ULON, TLAT, TLON)
!call read_topography(NX, NY,  KMT,  KMU)
!
!
!if (do_output()) then
!   write(*,*) ''
!   write(*,'(''Grid Dimensions NX = '',I4,'', NY = '',I4)') NX, NY
!   write(*,*) ''
!endif
!
!allocate(group_members(group_size)) ! this is module global
!
!!call build_my_group(group_size, group_members)
!
!deallocate(ULAT, ULON, TLAT, TLON, KMT, KMU)

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
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_par_msg/utilities/test_distributed_static_data.f90 $
! $Id: test_distributed_static_data.f90 6256 2013-06-12 16:19:10Z thoar $
! $Revision: 6256 $
! $Date: 2013-06-12 10:19:10 -0600 (Wed, 12 Jun 2013) $
