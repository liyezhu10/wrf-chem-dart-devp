! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: filter_mod.f90 10343 2016-06-07 21:51:48Z hendric $

module pseudo_orbit_DA_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength

!Du adds operator(+) operator(<)
use time_manager_mod,      only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                  operator(-), operator(+), operator(<), print_time

use utilities_mod,         only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                  logfileunit, nmlfileunit, timestamp,  E_ALLMSG, &
                                  do_output, find_namelist_in_file, check_namelist_read,      &
                                  open_file, close_file, do_nml_file, do_nml_term

use assim_model_mod,       only : static_init_assim_model, get_model_size, &
                                  get_state_meta_data

use obs_model_mod,         only : advance_state

use ensemble_manager_mod,  only : init_ensemble_manager, end_ensemble_manager,     &
                                  ensemble_type, get_my_num_copies,                &
                                  all_vars_to_all_copies, all_copies_to_all_vars,  &
                                  get_copy_owner_index, get_ensemble_time,         &
                                  set_ensemble_time, get_my_num_vars, get_my_vars

use mpi_utilities_mod,     only : initialize_mpi_utilities, finalize_mpi_utilities,           &
                                  my_task_id, task_sync, sum_across_tasks

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state

use io_filenames_mod,      only : io_filenames_init, file_info_type

use location_mod,          only : location_type


!------------------------------------------------------------------------------

implicit none
private

public :: pda_main

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/pda/filter/filter_mod_pda.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 10343 $"
character(len=128), parameter :: revdate  = "$Date: 2016-06-07 15:51:48 -0600 (Tue, 07 Jun 2016) $"

! Some convenient global storage items
character(len=129)      :: msgstring

integer                 :: trace_level, timestamp_level

!Some parameters require to be defined
logical  :: output_restart           = .true.
logical  :: output_restart_mean      = .false.
logical  :: single_restart_file_in   = .false. ! all copies read from 1 file
logical  :: single_restart_file_out  = .false. ! all copies written to 1 file
logical  :: direct_netcdf_read       = .true.  ! default to read from netcdf file
logical  :: direct_netcdf_write      = .true. ! default to write to netcdf file
logical  :: use_restart_list         = .true. ! read the list restart file names from a file

character(len = 129) :: restart_in_file_name  = 'NULL',     &
                        restart_out_file_name = 'pda_output', &
                        adv_ens_command       = './advance_model.csh'

integer :: tasks_per_model_advance = 1

! IO options
logical :: add_domain_extension  = .false. ! add _d0X to output filenames. Note this is always done for X>1
logical :: overwrite_state_input = .false. ! overwrites model netcdf files with output from filter

character(len = 129) :: inf_in_file_name(2)       = 'not_initialized',    &
                        inf_out_file_name(2)      = 'not_initialized'







!----------------------------------------------------------------
! Namelist input with default values, should make a namelist for pda seperately
!
integer  :: async = 0

!restart_list_file(1)='pda_ic_name_list'   !only this one should in the namelist
character(len=512) :: restart_list_file(10)        = 'null' ! name of files containing a list of restart files 

! Number of minimisation should be setup more standardised, or use other criteria 
! to stop the minimazation, should be in the namelist
integer   ::  n_GD=20

! Define assimilation window size
integer   ::  window_size=8

! gd_step_size is the minimisation step size for Gradient Descent, adjusted during 
! the minimisation, double the step size when lower cost function is achieved as 
! long as it is smaller than the gd_max_step_size, shrink step size when it fails.
! gd_initial_step_size=0.001_r8    
! this may use to calculate the number of minimizations, not used now
real(r8)  ::  gd_step_size= 0.0005_r8 !0.05_r8!0.0005_r8
real(r8)  ::  gd_max_step_size= 0.001_r8!0.05_r8!0.001_r8

namelist /pseudo_orbit_DA_nml/ async,             &
                               restart_list_file, &
                               n_GD,              &
                               window_size,       &
                               gd_step_size,      &
                               gd_max_step_size


!----------------------------------------------------------------

contains

!----------------------------------------------------------------
!> The code does pda data assimilation without an adjoint.

subroutine pda_main()

type(time_type)          :: time1
integer(i8)              :: model_size, var_ind
integer                  :: i, j, iunit, io
logical                  :: read_time_from_file
type(file_info_type)     :: file_info
type(ensemble_type)      :: pda_ens_handle, mismatch_ens_handle, ens_update_copy, ens_normalization
integer                  :: n_DA, var_type, my_num_vars, ivar!, seq_len
real(r8)                 :: mis_cost, mis_cost_previous
type(location_type)      :: location
integer(i8), allocatable :: my_vars(:)

call pda_initialize_modules_used() ! static_init_model called in here

! Read the namelist entry
call find_namelist_in_file("input.nml", "pseudo_orbit_DA_nml", iunit)
read(iunit, nml = pseudo_orbit_DA_nml, iostat = io)
call check_namelist_read(iunit, io, "pseudo_orbit_DA_nml")


model_size = get_model_size()


call init_ensemble_manager(pda_ens_handle, window_size, model_size, 1, transpose_type_in=2)

file_info = io_filenames_init(pda_ens_handle, single_restart_file_in, single_restart_file_out,  &
              restart_in_file_name, restart_out_file_name, output_restart, direct_netcdf_read,  &
              direct_netcdf_write, output_restart_mean, add_domain_extension, use_restart_list, &
              restart_list_file, overwrite_state_input, inf_in_file_name, inf_out_file_name)


!read sequence of enkf states
call read_state(pda_ens_handle, file_info, read_time_from_file, time1)

call init_ensemble_manager(mismatch_ens_handle, window_size, model_size, 1, transpose_type_in=2)
call init_ensemble_manager(ens_update_copy,    window_size, model_size, 1, transpose_type_in=2)
call init_ensemble_manager(ens_normalization,            1, model_size, 1, transpose_type_in=2)


!generate scaling vector--------------------------------------
!this is not standardized, should read a scaling vector from a file

my_num_vars = get_my_num_vars(pda_ens_handle)
allocate(my_vars(my_num_vars))
call get_my_vars(pda_ens_handle,my_vars)

do ivar=1,my_num_vars
    var_ind = my_vars(ivar)

    call get_state_meta_data(pda_ens_handle, var_ind, location, var_type)

    if (var_type==1) then
        ens_normalization%copies(1, var_ind)=1!2
    elseif (var_type==2) then
        ens_normalization%copies(1, var_ind)=1!1
    elseif (var_type==3) then
        ens_normalization%copies(1, var_ind)=1!60
    else
        ens_normalization%copies(1, var_ind)=1!2

    endif
end do


!------------------------------------------------------------------



!We might want to do the PDA sequentially in the future, maybe using shell script to do it instead
Sequential_PDA: do n_DA=1,1 !seq_len-window_size+1

    !calculate mismatch cost function---------------------------------------------------
    call cal_mismatch_cost_function(pda_ens_handle, mismatch_ens_handle, async, &
                                    adv_ens_command, tasks_per_model_advance, mis_cost)
    mis_cost_previous=mis_cost

    write(*,*) mis_cost

    !------------------------------------------------------------------------------------

    ens_update_copy%copies = pda_ens_handle%copies
    ens_update_copy%time   = pda_ens_handle%time




    !!GD minimisation update-------------------------------------------------------------

    !The criteria of stoping the minimisation could be based on the mismatch cost function
    GD_runs: do j=1,n_GD

        Gradient_descent: do i=1, window_size

            if (i==1) then
                pda_ens_handle%copies(1,:) = ens_update_copy%copies(1,:) &
                        +gd_step_size*ens_normalization%copies(1,:)*mismatch_ens_handle%copies(1,:)

            else if (i<window_size) then

                pda_ens_handle%copies(i,:)=ens_update_copy%copies(i,:) &
                        -gd_step_size*ens_normalization%copies(1,:)*mismatch_ens_handle%copies(i-1,:) &
                        +gd_step_size*ens_normalization%copies(1,:)*mismatch_ens_handle%copies(i,:)

            else


                pda_ens_handle%copies(window_size,:)=ens_update_copy%copies(window_size,:) &
                        -gd_step_size*ens_normalization%copies(1,:)*mismatch_ens_handle%copies(window_size-1,:)

            endif

        end do Gradient_descent


        !!calculate mismatch cost function for updated sequence state vectors

        call cal_mismatch_cost_function(pda_ens_handle, mismatch_ens_handle, &
                                        async, adv_ens_command, tasks_per_model_advance, mis_cost)

        if ( mis_cost < mis_cost_previous ) then
            !increase the GD minimisation time step by a factor of 2,
            !should have some upper bound for the GD minimisation time step

            ens_update_copy%copies=pda_ens_handle%copies

            if ( gd_step_size < gd_max_step_size ) then
                gd_step_size=gd_step_size*2
            endif
            mis_cost_previous=mis_cost

        else

            !decrease the GD minimisation time step by a factor of 2
            !should have some lower bound for the GD minimisation time step
            pda_ens_handle%copies=ens_update_copy%copies
            gd_step_size     = gd_step_size/2
            gd_max_step_size = gd_step_size

            !may use a vector to store previous forward ensemble copies
            call cal_mismatch_cost_function(pda_ens_handle, mismatch_ens_handle, &
                                            async, adv_ens_command, tasks_per_model_advance, mis_cost)

            write(*,*) 'larger_mis cost', mis_cost


        endif
        write(*,*) mis_cost

    end do GD_runs


    !--------------------------------------------------------------------------------------------------


    !!!output pda final update

    call write_state(pda_ens_handle, file_info)


    !do i=1, window_size
    !    write(*,*) pda_ens_handle%copies(i,1:3)
    !end do

    call end_ensemble_manager(pda_ens_handle)
    call end_ensemble_manager(mismatch_ens_handle)
    call end_ensemble_manager(ens_update_copy)

end do Sequential_PDA


end subroutine pda_main

!-----------------------------------------------------------------------
!calculate mismatch cost function

subroutine cal_mismatch_cost_function(pda_ens_handle, mismatch_ens_handle, async, adv_ens_command, tasks_per_model_advance, mis_cost)

type(ensemble_type),  intent(in)    :: pda_ens_handle
type(ensemble_type),  intent(inout) :: mismatch_ens_handle
integer,              intent(in)    :: async
integer,              intent(in)    :: tasks_per_model_advance
character(len = 129), intent(in)    :: adv_ens_command
real(r8),             intent(out)   :: mis_cost

integer         :: i,j
type(time_type) :: target_time, ens_time
real(r8)        :: result_sum, sum_variable

window_size=pda_ens_handle%num_copies

call get_ensemble_time(pda_ens_handle, 1, ens_time)
call get_ensemble_time(pda_ens_handle, 2, target_time)

mismatch_ens_handle%copies=pda_ens_handle%copies
mismatch_ens_handle%time=pda_ens_handle%time

if(.not. allocated(mismatch_ens_handle%vars)) &
    allocate(mismatch_ens_handle%vars(mismatch_ens_handle%num_vars, mismatch_ens_handle%my_num_copies))

call all_copies_to_all_vars(mismatch_ens_handle)

!!!this is not ideal, shall be able to pass an arrary of target time to advance_state, things needs to be changed in advance_state to adapt this
do i=1, window_size
    mismatch_ens_handle%time(i)=ens_time
end do

call advance_state(mismatch_ens_handle, window_size, target_time, async, adv_ens_command, tasks_per_model_advance)

write(*,*) window_size

call all_vars_to_all_copies(mismatch_ens_handle)

!calculate mismatches
do i=1, window_size-1
    mismatch_ens_handle%copies(i,:)=pda_ens_handle%copies(i+1,:)-mismatch_ens_handle%copies(i,:)
end do

!calculate mismatch cost function
result_sum   = 0
sum_variable = 0
do i=1,window_size-1

    do j=1,mismatch_ens_handle%my_num_vars
        sum_variable=mismatch_ens_handle%copies(i,j)*mismatch_ens_handle%copies(i,j)+sum_variable
    end do

    call sum_across_tasks(sum_variable,result_sum)

end do

mis_cost=result_sum/((window_size-1)*1.0_r8)

end subroutine cal_mismatch_cost_function

!-------------------------------------------------------------------------

subroutine pda_initialize_modules_used()

! Initialize modules used that require it
call initialize_mpi_utilities('pseudo_orbit_DA')

call register_module(source,revision,revdate)

! Initialize the model class data now that obs_sequence is all set up
call trace_message('Before init_model call')
call static_init_assim_model()
call trace_message('After  init_model call')
call state_vector_io_init()
call trace_message('After  init_state_vector_io call')

end subroutine pda_initialize_modules_used

!-------------------------------------------------------------------------

subroutine trace_message(msg, label, threshold)

character(len=*), intent(in)           :: msg
character(len=*), intent(in), optional :: label
integer,          intent(in), optional :: threshold

! Write message to stdout and log file.
integer :: t

t = 0
if (present(threshold)) t = threshold

if (trace_level <= t) return

if (.not. do_output()) return

if (present(label)) then
   call error_handler(E_MSG,trim(label),trim(msg))
else
   call error_handler(E_MSG,'filter trace:',trim(msg))
endif

end subroutine trace_message

!-------------------------------------------------------------------------

subroutine print_ens_time(ens_handle, msg)

type(ensemble_type), intent(in) :: ens_handle
character(len=*), intent(in) :: msg

! Write message to stdout and log file.
type(time_type) :: mtime

if (trace_level <= 0) return

if (do_output()) then
   if (get_my_num_copies(ens_handle) < 1) return
   call get_ensemble_time(ens_handle, 1, mtime)
   call print_time(mtime, ' filter trace: '//msg, logfileunit)
   call print_time(mtime, ' filter trace: '//msg)
endif

end subroutine print_ens_time

!-------------------------------------------------------------------

end module pseudo_orbit_DA_mod

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/pda/filter/filter_mod.f90 $
! $Id: filter_mod.f90 10343 2016-06-07 21:51:48Z hendric $
! $Revision: 10343 $
! $Date: 2016-06-07 15:51:48 -0600 (Tue, 07 Jun 2016) $
