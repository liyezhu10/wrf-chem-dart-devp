! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!----------------------------------------------------------------------
!> purpose: generate initial inflation files so an experiment starting up
!> can always start from a restart file without having to alter the namelist
!> between cycles 1 and 2.
!>
!> an alternative to running this program is to use the nco utilities thus:
!>
!> Here is an example using version 4.4.2 or later of the NCO tools:
!>   ncap2 -s "T=1.0;U=1.0;V=1.0" wrfinput_d01 prior_inflation_mean.nc
!>   ncap2 -s "T=0.6;U=0.6;V=0.6" wrfinput_d01 prior_inflation_sd.nc'
!>
!----------------------------------------------------------------------

program fill_inflation_restart

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  find_namelist_in_file, check_namelist_read,   &
                                  initialize_utilities, finalize_utilities

use       assim_model_mod, only : static_init_assim_model

use      time_manager_mod, only : time_type, set_time, print_time, print_date, operator(-), &
                                  get_calendar_type, NO_CALENDAR

use  ensemble_manager_mod, only : init_ensemble_manager, ensemble_type, get_copy, &
                                  map_pe_to_task, map_task_to_pe
use   state_vector_io_mod, only : state_vector_io_init, read_state, write_state

use   state_structure_mod, only : state_structure_info, get_model_variable_indices, &
                                  get_num_domains, get_num_variables, get_io_num_dims, &
                                  get_variable_name, get_dim_name, get_dim_length

use      io_filenames_mod, only : io_filenames_init, file_info_type,       &
                                  stage_metadata_type, get_stage_metadata, &
                                  get_restart_filename,                    &
                                  set_file_metadata, file_info_dump,       &
                                  set_io_copy_flag, READ_COPY, WRITE_COPY

use direct_netcdf_mod,     only : finalize_single_file_io, write_augmented_state

use netcdf_utilities_mod,  only : nc_open_file_readwrite, nc_add_attribute_to_variable, &
                                  nc_begin_define_mode, nc_end_define_mode,             &
                                  nc_define_real_variable, nc_define_dimension,         &
                                  nc_close_file, nc_put_variable, nc_variable_exists

use             model_mod, only : static_init_model, get_model_size,       &
                                  get_state_meta_data, model_interpolate

use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


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
integer                       :: num_ens = 2 !#! {prior,posterior}_inf_{mean,sd}
character(len=256)            :: input_state_files(MAX_FILES)  = 'null'
character(len=256)            :: output_state_files(MAX_FILES) = 'null'
logical                       :: write_prior_inf, write_post_inf = .false.
logical                       :: verbose = .FALSE.

real(r8) :: prior_inf_mean = MISSING_R8
real(r8) :: prior_inf_sd   = MISSING_R8
real(r8) :: post_inf_mean  = MISSING_R8
real(r8) :: post_inf_sd    = MISSING_R8
integer  :: ss_inflate_index    = 1 
integer  :: ss_inflate_sd_index = 2

!>@todo FIXME
! output_state_files shouldn't be in the namelist, right?  we hardcode this now.
! in the multi-file case we have fixed names.  in the single-file case we can
! update whatever they give us as input.
namelist /fill_inflation_restart_nml/  prior_inf_mean, prior_inf_sd, &
                               post_inf_mean, post_inf_sd, verbose,  &
                               write_prior_inf, write_post_inf, &
                               input_state_files, output_state_files, single_file
! io variables
integer                   :: iunit, io
integer, allocatable      :: ios_out(:)
type(file_info_type)      :: file_info_input, file_info_output
type(stage_metadata_type) :: input_restart_files, output_restart_files
logical :: read_time_from_file     = .true.

! model state variables
type(ensemble_type)   :: ens_handle

type(time_type)       :: model_time
integer(i8)           :: model_size
real(r8), allocatable :: interp_vals(:)

! misc. variables
integer :: idom, imem, num_domains, idomain

! message strings
character(len=512) :: my_base, my_desc, my_stage
character(len=512) :: string1, string2, string3

character(len=256), allocatable :: file_array_input(:,:), file_array_output(:,:)

!======================================================================
! start of executable code
!======================================================================

call initialize_modules_used()

call find_namelist_in_file("input.nml", "fill_inflation_restart_nml", iunit)
read(iunit, nml = fill_inflation_restart_nml, iostat = io)
call check_namelist_read(iunit, io, "fill_inflation_restart_nml")

!----------------------------------------------------------------------
! Calling static_init_assim_model() is required, which also calls
! static_init_model(), so there is no need to explicitly call it.
!----------------------------------------------------------------------

call static_init_assim_model()

!----------------------------------------------------------------------
! initialization code, model size
!----------------------------------------------------------------------

model_size = get_model_size()

!----------------------------------------------------------------------
! read/write restart files
!----------------------------------------------------------------------

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, num_ens, model_size)

! Allocate space for file arrays.  contains a matrix of files (num_ens x num_domains)
! If perturbing from a single instance the number of input files does not have to
! be ens_size but rather a single file (or multiple files if more than one domain)

num_domains = get_num_domains()

file_array_input  = RESHAPE(input_state_files,  (/      1,  num_domains/))
file_array_output = RESHAPE(output_state_files, (/num_ens,  num_domains/))

if(single_file) then
  call error_handler(E_ERR, 'fill_inflation_restart: ', 'single_file not yet supported', &
                     source, revision, revdate)
  file_array_output = file_array_input
endif

! set up the files to be read - the data will not be used.
! the variables and shape of them are a template for the
! inflation values.
call io_filenames_init(file_info_input,             &
                       ncopies      = num_ens,      &
                       cycling      = single_file,  &
                       single_file  = single_file,  &
                       restart_files = file_array_input(:,:))

! Read the template file to get the shape of netCDF file
! and its variables. It is possible to have multiple domains
! but only require one member.
write(my_base,'(A)') 'template'
write(my_desc,'(A)') 'template file'
call set_file_metadata(file_info_input,                  &
                       cnum     = 1,                     &
                       fnames   = file_array_input(1,:), &
                       basename = my_base,               &
                       desc     = my_desc)

call set_io_copy_flag(file_info_input,    &
                      cnum    = 1,     &
                      io_flag = READ_COPY)

input_restart_files = get_stage_metadata(file_info_input)

do idom = 1, num_domains
   imem = 1
   write(*, *) '- Reading File : ', &
      trim(get_restart_filename(input_restart_files, copy=imem, domain=idom))
enddo

call read_state(ens_handle, file_info_input, read_time_from_file, model_time)

if(single_file) then
   call nc_insert_variable('filter_input.nc', varname='prior_inf_mean', copynum=1)
else

   call io_filenames_init(file_info_output,           &
                          ncopies      = num_ens,     &
                          cycling      = single_file, &
                          single_file  = single_file, &
                          restart_files = file_array_output)
   
   ! The inflation file names should match the hardcoded names
   ! that filter requires.  (We no longer allow the user to specify
   ! the inflation names, to simplify the number of namelist items.)
   !
   ! We will create input_priorinf_mean.nc and input_priorinf_sd.nc
   ! and/or input_postinf_mean.nc and input_postinf_sd.nc

   if(write_prior_inf) then  
      call fill_inflation_files(prior_inf_mean, prior_inf_sd, 'prior')
   endif
   
   if(write_post_inf) then  
      call fill_inflation_files(post_inf_mean, post_inf_sd, 'post')
   endif

endif

deallocate(file_array_input, file_array_output)

call exit(0)

!======================================================================
contains
!======================================================================

!----------------------------------------------------------------------
!> fill inflation values in a separate file

subroutine fill_inflation_files(inf_mean, inf_sd, stage)
real(r8),         intent(in) :: inf_mean
real(r8),         intent(in) :: inf_sd
character(len=*), intent(in) :: stage

if (inf_mean == MISSING_R8 .or. inf_sd == MISSING_R8) then
   write(*,*) 'you must specify both inf_mean and inf_sd values'
   write(*,*) 'you have "',trim(stage),'_inf_mean = ', inf_mean,'" and '
   write(*,*) '         "',trim(stage),'_inf_sd   = ', inf_sd,  '"     '
   return
endif
ens_handle%copies(ss_inflate_index   , :) = prior_inf_mean
ens_handle%copies(ss_inflate_sd_index, :) = prior_inf_sd

write(my_stage,'(3A)') 'input_',stage,'inf'
write(my_base, '(A)')  'mean'
write(my_desc, '(2A)') stage, ' inflation mean'
call set_file_metadata(file_info_output,    &
                       cnum     = 1,        &
                       stage    = my_stage, &
                       basename = my_base,  &
                       desc     = my_desc)

call set_io_copy_flag(file_info_output,    &
                      cnum    = 1,         &
                      io_flag = WRITE_COPY)

write(my_base, '(A)') 'sd'
write(my_desc, '(2A)') stage, ' inflation sd'
call set_file_metadata(file_info_output,    &
                       cnum     = 2,        &
                       stage    = my_stage, &
                       basename = my_base,  &
                       desc     = my_desc)

call set_io_copy_flag(file_info_output,    &
                      cnum    = 2,         &
                      io_flag = WRITE_COPY)

output_restart_files = get_stage_metadata(file_info_output)

do idom = 1, num_domains
   do imem = 1,2 ! inflation mean and sd
      write(*, *) '- Writing "',trim(stage), '" Inflation File : ', &
                  trim(get_restart_filename(output_restart_files, imem, domain=idom))
   enddo
enddo

call write_state(ens_handle, file_info_output)

end subroutine fill_inflation_files

!----------------------------------------------------------------------
!> this is code that was lifted from another module and doesn't
!> work correctly.  right now if you have all your ensemble members
!> in single netcdf variables (a "single" or "combined" file)
!> you can't use this tool.

subroutine nc_insert_variable(filename, varname, copynum)
character(len=*), intent(in) :: filename
character(len=*), intent(in) :: varname
integer,          intent(in) :: copynum

character(len=*), parameter :: routine = 'nc_insert_variable'

integer :: ncid, ret, ivar, jdim
integer :: dimlen, domain, ndims, dcount
integer :: timeindex, copy_index
integer, dimension(NF90_MAX_VAR_DIMS) :: dim_lengths
integer, dimension(NF90_MAX_VAR_DIMS) :: start_point
character(len=NF90_MAX_NAME) :: dimname, attname, vname, extraname
character(len=NF90_MAX_NAME), allocatable :: dim_names(:)

integer(i8) :: model_size
real(r8) :: val
real(r8), allocatable :: varvals(:,:,:)
real(r8), allocatable :: temp_ens(:)

model_size = get_model_size()
allocate(temp_ens(model_size))

! NF90_OPEN
ncid = nc_open_file_readwrite(filename, routine)

domain = 1 !>@todo ONLY ONE DOMAIN FOR SINGLE FILE OUTPUT

! NF90_REDEF            ! put it into define mode
call nc_begin_define_mode(ncid, routine, filename)
  ! ...
  ! NF90_DEF_DIM        ! define additional dimensions (if any)
  do ivar = 1, get_num_variables(domain)
     ndims =      get_io_num_dims(  domain, ivar )
     vname = trim(get_variable_name(domain, ivar))
     print*, 'variable name  ', get_variable_name(domain, ivar)
     print*, 'number of dims ', ndims
     allocate(dim_names(ndims))
     ! set the defaults and then change any needed below
     dim_lengths(:) = 1
     start_point(:) = 1
     dcount = 0

     do jdim = 1, ndims
        dimname = get_dim_name(domain, ivar, jdim)
        if (dimname == 'time') then
           dcount = dcount + 1
           start_point(dcount) = jdim
           dim_names(  dcount) = dimname
           print*, 'time id ', jdim
        else if (dimname == 'member') then
           continue   ! extra vars have no member dimension
        else 
           dcount = dcount + 1
           dim_lengths(dcount) = get_dim_length(domain, ivar, jdim)
           dim_names(  dcount) = dimname
        endif
     enddo

     write(extraname,'(a,"_",a)') trim(vname), trim('priorinf_mean')
     write( *       ,'(a,"_",a)') trim(vname), trim('priorinf_mean')
     write( *       , *          ) 'start point ', start_point(1:dcount)
     write( *       , *          ) 'dim lengths ', dim_lengths(1:dcount)
     
     do jdim = 1, dcount
        write( *       , *      ) 'dim names   ', trim(dim_names(jdim))
     enddo
  enddo


  ! NF90_DEF_VAR        ! define additional variables (if any)

  if(.not. nc_variable_exists(ncid, extraname)) then
     write( * , * ) 'writting new variable "', extraname, '"'
     call nc_define_real_variable(ncid, extraname, dim_names, routine, filename)
  endif
  !   ...
  ! NF90_PUT_ATT        ! define other attributes (if any)
  attname = 'prior inflation mean'
  call  nc_add_attribute_to_variable(ncid, extraname, attname, val, routine, filename)
  !   ...


! NF90_ENDDEF           ! check definitions, leave define mode
call nc_end_define_mode(ncid, routine, filename)
  !  ...
  ! NF90_PUT_VAR        ! provide new variable values

temp_ens = prior_inf_mean
print*, 'inquire variable id for '//trim(extraname)
print*, 'put values for '//trim(extraname)
print*, 'varvalues      ', temp_ens
call nc_put_variable(ncid, extraname, temp_ens, routine, filename)

  !   ...
! NF90_CLOSE            ! close netCDF dataset
call nc_close_file(ncid, routine, filename)

deallocate(temp_ens)

end subroutine nc_insert_variable


!----------------------------------------------------------------------
!> initialize modules that need it

subroutine initialize_modules_used()

call initialize_mpi_utilities('fill_inflation_restart')

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

end program fill_inflation_restart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

