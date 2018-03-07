! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_group_mean

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use       assim_model_mod, only : static_init_assim_model

use      time_manager_mod, only : time_type, set_time, print_date, operator(-), &
                                  NO_CALENDAR

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type

use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : get_num_domains, get_model_variable_indices, &
                                  state_structure_info

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename, set_file_metadata, &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use distributed_state_mod, only : create_state_window, free_state_window

use             model_mod, only : static_init_model

use netcdf

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

logical                       :: single_file = .false.
integer                       :: num_ens = 1
character(len=256)            :: input_state_files(MAX_FILES)  = 'null'
character(len=256)            :: output_state_files(MAX_FILES) = 'null'

namelist /test_group_mean_nml/ num_ens, single_file, input_state_files, output_state_files

! io variables
integer                   :: iunit, io
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file = .true.

! model state variables
type(ensemble_type)   :: ens_handle

type(time_type)       :: model_time
integer(i8)           :: model_size

! misc. variables
integer :: idom, imem, num_domains

! message strings
character(len=512) :: my_base, my_desc
character(len=512) :: string1, string2, string3

character(len=256), allocatable :: file_array_input(:,:)

!======================================================================
! start of executable code
!======================================================================

call initialize_modules_used()

call find_namelist_in_file("input.nml", "test_group_mean_nml", iunit)
read(iunit, nml = test_group_mean_nml, iostat = io)
call check_namelist_read(iunit, io, "test_group_mean_nml")

!----------------------------------------------------------------------
! Calling static_init_assim_model() is required for all tests.
! It also calls static_init_model(), so there is no need to explicitly call
! that. Furthermore, the low-order models have no check in them to prevent
! static_init_model() from being called twice, so it BOMBS if you call both.
!----------------------------------------------------------------------

call static_init_assim_model()

!----------------------------------------------------------------------
! read/write restart files
!----------------------------------------------------------------------

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, num_ens+1, model_size)

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)

num_domains = get_num_domains()

allocate(file_array_input( num_ens, num_domains))
file_array_input  = RESHAPE(input_state_files,  (/num_ens,  num_domains/))

! Test the read portion.
call io_filenames_init(file_info_input,              &
                       ncopies       = num_ens,      &
                       cycling       = single_file,  &
                       single_file   = single_file,  &
                       restart_files = file_array_input)

do imem = 1, num_ens
   write(my_base,'(A,I2)') 'inens_',    imem
   write(my_desc,'(A,I2)') 'input ens', imem
   call set_file_metadata(file_info_input,                      &
                          cnum     = imem,                      &
                          fnames   = file_array_input(imem,:),  &
                          basename = my_base,                   &
                          desc     = my_desc)

   call set_io_copy_flag(file_info_input,    &
                         cnum    = imem,     &
                         io_flag = READ_COPY)
enddo

input_restart_files = get_stage_metadata(file_info_input)

do idom = 1, num_domains
   do imem = 1, num_ens
      write(string1, *) '- Reading File : ', trim(get_restart_filename(input_restart_files, imem, domain=idom))
   enddo
enddo

call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

deallocate(file_array_input)

!----------------------------------------------------------------------
! Check window
!----------------------------------------------------------------------

call create_state_window(ens_handle)

call free_state_window(ens_handle)

!----------------------------------------------------------------------
! finalize test_group_mean
!----------------------------------------------------------------------

call finalize_modules_used()

!======================================================================
contains
!======================================================================

!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('test_group_mean')

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

end program test_group_mean

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
