! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: filter_mod.f90 10343 2016-06-07 21:51:48Z hendric $

module filter_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength
use obs_sequence_mod,      only : read_obs_seq, obs_type, obs_sequence_type,                  &
                                  get_obs_from_key, set_copy_meta_data, get_copy_meta_data,   &
                                  get_obs_def, get_time_range_keys, set_obs_values, set_obs,  &
                                  write_obs_seq, get_num_obs, get_obs_values, init_obs,       &
                                  assignment(=), get_num_copies, get_qc, get_num_qc, set_qc,  &
                                  static_init_obs_sequence, destroy_obs, read_obs_seq_header, &
                                  set_qc_meta_data, get_first_obs, get_obs_time_range,        &
                                  delete_obs_from_seq, delete_seq_head,                       &
                                  delete_seq_tail, replace_obs_values, replace_qc,            &
                                  destroy_obs_sequence, get_qc_meta_data, add_qc
                                  
use obs_def_mod,           only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                                  get_obs_kind
use obs_def_utilities_mod, only : set_debug_fwd_op
!Du adds operator(+) operator(<)
use time_manager_mod,      only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                  operator(-), operator(+), operator(<), print_time
use utilities_mod,         only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                  logfileunit, nmlfileunit, timestamp,  &
                                  do_output, find_namelist_in_file, check_namelist_read,      &
                                  open_file, close_file, do_nml_file, do_nml_term
!Du add get_model_time_step
use assim_model_mod,       only : static_init_assim_model, get_model_size,                    &
                                  end_assim_model,  pert_model_copies, get_model_time_step
use assim_tools_mod,       only : filter_assim, set_assim_tools_trace, get_missing_ok_status, &
                                  test_state_copies
use obs_model_mod,         only : move_ahead, advance_state, set_obs_model_trace
use ensemble_manager_mod,  only : init_ensemble_manager, end_ensemble_manager,                &
                                  ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                  all_vars_to_all_copies, all_copies_to_all_vars,             &
                                  compute_copy_mean, compute_copy_mean_sd,                    &
                                  compute_copy_mean_var, duplicate_ens, get_copy_owner_index, &
                                  get_ensemble_time, set_ensemble_time, broadcast_copy,       &
                                  prepare_to_read_from_vars, prepare_to_write_to_vars,        &
                                  prepare_to_read_from_copies,  get_my_num_vars,              &
                                  prepare_to_write_to_copies, get_ensemble_time,              &
                                  map_task_to_pe,  map_pe_to_task, prepare_to_update_copies,  &
                                  copies_in_window, set_num_extra_copies, get_allow_transpose, &
                                  all_copies_to_all_vars, allocate_single_copy,               &
                                  get_single_copy, put_single_copy, deallocate_single_copy
use adaptive_inflate_mod,  only : adaptive_inflate_end, do_varying_ss_inflate,                &
                                  do_single_ss_inflate, inflate_ens, adaptive_inflate_init,   &
                                  do_obs_inflate, adaptive_inflate_type,                      &
                                  output_inflate_diagnostics, log_inflation_info, &
                                  get_minmax_task_zero
use mpi_utilities_mod,     only : initialize_mpi_utilities, finalize_mpi_utilities,           &
                                  my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                  task_count
use smoother_mod,          only : smoother_read_restart, advance_smoother,                    &
                                  smoother_gen_copy_meta_data, smoother_write_restart,        &
                                  init_smoother, do_smoothing, smoother_mean_spread,          &
                                  smoother_assim,            &
                                  smoother_ss_diagnostics, smoother_end, set_smoother_trace

use random_seq_mod,        only : random_seq_type, init_random_seq, random_gaussian

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state

use io_filenames_mod,      only : io_filenames_init, file_info_type

use forward_operator_mod,  only : get_obs_ens_distrib_state
use quality_control_mod,   only : initialize_qc

use state_space_diag_mod,  only : filter_state_space_diagnostics, netcdf_file_type, &
                                  init_diag_output, finalize_diag_output,           &
                                  skip_diag_files

! state copy meta data
use copies_on_off_mod, only : ENS_MEAN_COPY, ENS_SD_COPY, PRIOR_INF_COPY, &
                              PRIOR_INF_SD_COPY,POST_INF_COPY, POST_INF_SD_COPY, &
                              SPARE_PRIOR_MEAN, SPARE_PRIOR_SPREAD, SPARE_PRIOR_INF_MEAN, &
                              SPARE_PRIOR_INF_SPREAD, SPARE_POST_INF_MEAN, &
                              SPARE_POST_INF_SPREAD, query_copy_present


!DuDu adds.....
use model_mod, 		  only  : adv_1step, comp_dt


!------------------------------------------------------------------------------

implicit none
private

public :: filter_sync_keys_time, &
          filter_set_initial_time, &
          pda_main

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/pda/filter/filter_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 10343 $"
character(len=128), parameter :: revdate  = "$Date: 2016-06-07 15:51:48 -0600 (Tue, 07 Jun 2016) $"

! Some convenient global storage items
character(len=129)      :: msgstring
type(obs_type)          :: observation

integer                 :: trace_level, timestamp_level

! Defining whether diagnostics are for prior or posterior
integer, parameter :: PRIOR_DIAG = 0, POSTERIOR_DIAG = 2

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
logical  :: output_restart      = .false.
logical  :: output_restart_mean = .false.
integer  :: tasks_per_model_advance = 1
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days      = -1
integer  :: first_obs_seconds   = -1
integer  :: last_obs_days       = -1
integer  :: last_obs_seconds    = -1
! Assimilation window; defaults to model timestep size.
integer  :: obs_window_days     = -1
integer  :: obs_window_seconds  = -1
! Control diagnostic output for state variables
integer  :: num_output_state_members = 0
integer  :: num_output_obs_members   = 0
integer  :: output_interval     = 1
integer  :: num_groups          = 1
logical  :: output_forward_op_errors = .false.
logical  :: output_timestamps        = .false.
logical  :: trace_execution          = .false.
logical  :: silence                  = .false.
logical  :: distributed_state = .true. ! Default to do state complete forward operators.

! IO options
logical            :: add_domain_extension         = .false. ! add _d0X to output filenames. Note this is always done for X>1
logical            :: use_restart_list             = .false. ! read the list restart file names from a file
character(len=512) :: restart_list_file(10)        = 'null' ! name of files containing a list of restart files (only used if use_restart_list = .true. 1 file per domain
logical            :: overwrite_state_input        = .false. ! overwrites model netcdf files with output from filter
logical            :: perturb_from_single_instance = .false. ! Read in a single file and perturb this to create an ensemble
real(r8)           :: perturbation_amplitude = 0.2_r8
logical            :: single_restart_file_in       = .false. ! all copies read from 1 file
logical            :: single_restart_file_out      = .false. ! all copies written to 1 file
logical            :: direct_netcdf_read = .true.  ! default to read from netcdf file
logical            :: direct_netcdf_write = .true. ! default to write to netcdf file

character(len = 129) :: obs_sequence_in_name  = "obs_seq.out",    &
                        obs_sequence_out_name = "obs_seq.final",  &
                        restart_in_file_name  = 'filter_ics',     &
                        restart_out_file_name = 'filter_restart', &
                        adv_ens_command       = './advance_model.csh'

!                  == './advance_model.csh'    -> advance ensemble using a script

! Inflation namelist entries follow, first entry for prior, second for posterior
! inf_flavor is 0:none, 1:obs space, 2: varying state space, 3: fixed state_space
integer              :: inf_flavor(2)             = 0
logical              :: inf_initial_from_restart(2)    = .false.
logical              :: inf_sd_initial_from_restart(2) = .false.

! old way
logical              :: inf_output_restart(2)     = .false.
! new way
!logical              :: inf_output_prior(2) = .false. ! mean sd
!logical              :: inf_output_post(2)  = .false. ! mean sd

logical              :: inf_deterministic(2)      = .true.
character(len = 129) :: inf_in_file_name(2)       = 'not_initialized',    &
                        inf_out_file_name(2)      = 'not_initialized',    &
                        inf_diag_file_name(2)     = 'not_initialized'
real(r8)             :: inf_initial(2)            = 1.0_r8
real(r8)             :: inf_sd_initial(2)         = 0.0_r8
real(r8)             :: inf_damping(2)            = 1.0_r8
real(r8)             :: inf_lower_bound(2)        = 1.0_r8
real(r8)             :: inf_upper_bound(2)        = 1000000.0_r8
real(r8)             :: inf_sd_lower_bound(2)     = 0.0_r8
logical              :: output_inflation          = .true. ! This is for the diagnostic files, no separate option for prior and posterior

namelist /filter_nml/ async, adv_ens_command, ens_size, tasks_per_model_advance,    &
   output_restart, obs_sequence_in_name, obs_sequence_out_name, &
   restart_in_file_name, restart_out_file_name, init_time_days, init_time_seconds,  &
   first_obs_days, first_obs_seconds, last_obs_days, last_obs_seconds,              &
   obs_window_days, obs_window_seconds, &
   num_output_state_members, num_output_obs_members, output_restart_mean,           &
   output_interval, num_groups, trace_execution,                 &
   output_forward_op_errors, output_timestamps,                 &
   inf_flavor, inf_initial_from_restart, inf_sd_initial_from_restart,               &
   inf_output_restart, inf_deterministic, inf_in_file_name, inf_damping,            &
   inf_out_file_name, inf_diag_file_name, inf_initial, inf_sd_initial,              &
   inf_lower_bound, inf_upper_bound, inf_sd_lower_bound,           &
   silence, direct_netcdf_read, direct_netcdf_write, output_inflation, &
   distributed_state, add_domain_extension, use_restart_list, restart_list_file, &
   overwrite_state_input, single_restart_file_in, single_restart_file_out, &
   perturb_from_single_instance, perturbation_amplitude


!----------------------------------------------------------------

contains

!----------------------------------------------------------------
!> The code does not use %vars arrays except:
!> * Task 0 still writes the obs_sequence file, so there is a transpose (copies to vars) and 
!> sending the obs_fwd_op_ens_handle%vars to task 0. Keys is also size obs%vars.
!> * If you read dart restarts state_ens_handle%vars is allocated.
!> * If you write dart diagnostics state_ens_handle%vars is allocated.
!> * If you are not doing distributed forward operators state_ens_handle%vars is allocated
subroutine pda_main()



type(ensemble_type)         :: state_ens_handle, obs_fwd_op_ens_handle, qc_ens_handle
type(obs_sequence_type)     :: seq
type(netcdf_file_type)      :: PriorStateUnit, PosteriorStateUnit
type(time_type)             :: time1, first_obs_time, last_obs_time
type(time_type)             :: curr_ens_time, next_ens_time, window_time
type(adaptive_inflate_type) :: prior_inflate, post_inflate

integer,    allocatable :: keys(:)
integer(i8)             :: model_size
integer                 :: i, iunit, io, time_step_number, num_obs_in_set
integer                 :: ierr, last_key_used, key_bounds(2)
integer                 :: in_obs_copy, obs_val_index
integer                 :: output_state_mean_index, output_state_spread_index
integer                 :: prior_obs_mean_index, posterior_obs_mean_index
integer                 :: prior_obs_spread_index, posterior_obs_spread_index
! Global indices into ensemble storage - observations
integer                 :: OBS_VAL_COPY, OBS_ERR_VAR_COPY, OBS_KEY_COPY
integer                 :: OBS_GLOBAL_QC_COPY,OBS_EXTRA_QC_COPY
integer                 :: OBS_MEAN_START, OBS_MEAN_END
integer                 :: OBS_VAR_START, OBS_VAR_END, TOTAL_OBS_COPIES
integer                 :: input_qc_index, DART_qc_index
logical                 :: read_time_from_file

integer :: num_extras ! the extra ensemble copies

type(file_info_type) :: file_info

logical                 :: ds, all_gone, allow_missing

! real(r8), allocatable   :: temp_ens(:) ! for smoother
real(r8), allocatable   :: prior_qc_copy(:)


!!DuDu adds....
type(ensemble_type)         :: pda_ens_handle, ens_handle, ens_update, ens_update_copy
logical                 :: duplicate_time, free_descent
integer                 :: j, seq_len, n_GD, window_size, n_DA, i_row, i_col
real(r8)                :: value(1), mis_cost, mis_cost_previous, gd_step_size, gd_max_step_size,gd_initial_step_size
real(r8), allocatable   :: state_vector(:), mismatch_err(:), forward_vector(:)
type(time_type) 	    :: mtime
type(obs_def_type)	    :: obs_def
type(time_type)         :: time_step, target_time, ens_time
real(r8),       pointer :: adjoint(:, :)
real(r8)		:: rtemp

!Du number of minimisation should be setup more standardised, or use other criteria to stop the minimazation
n_GD=4000

!Du define assimilation window size
window_size=984;

call filter_initialize_modules_used() ! static_init_model called in here

!!!Du  should create filterpda namelist in input.nml

! Read the namelist entry
call find_namelist_in_file("input.nml", "filter_nml", iunit)
read(iunit, nml = filter_nml, iostat = io)
call check_namelist_read(iunit, io, "filter_nml")


!DuDu....

model_size = get_model_size()

call init_ensemble_manager(pda_ens_handle, window_size, model_size, 1, transpose_type_in=2)

!shall set pda_ens_handle%num_extras=0

single_restart_file_in=.false.
single_restart_file_out=.false.
restart_in_file_name='Null'
use_restart_list=.true.
restart_list_file(1)='pda_ic_name_list'
direct_netcdf_read=.true.

!write(*,*) restart_in_file_name
!write(*,*) restart_out_file_name
!write(*,*) output_restart
!write(*,*) direct_netcdf_read
!write(*,*) direct_netcdf_write
!write(*,*) output_restart_mean
!write(*,*) add_domain_extension
!write(*,*) use_restart_list
!write(*,*) restart_list_file
!write(*,*) overwrite_state_input
!write(*,*) inf_in_file_name
!write(*,*) inf_out_file_name

file_info = io_filenames_init(pda_ens_handle, single_restart_file_in, single_restart_file_out, &
restart_in_file_name, restart_out_file_name, output_restart, direct_netcdf_read, &
direct_netcdf_write, output_restart_mean, add_domain_extension, use_restart_list, &
restart_list_file, overwrite_state_input, inf_in_file_name, inf_out_file_name)

read_time_from_file=.true.


call read_state(pda_ens_handle, file_info, read_time_from_file, time1)

!do i=1, window_size
!    write(*,*) pda_ens_handle%copies(i,1:3)
!end do

! use free_descent to decide whether use model adjoint or free adjoint (Identity adjoint)
free_descent=.true.

if (.not. free_descent) allocate(adjoint (model_size, model_size))

time_step = get_model_time_step()
allocate(state_vector(model_size))
allocate(mismatch_err(model_size))
allocate(forward_vector(model_size))

!seq_len=200



! do PDA "sequentially"
Sequential_PDA: do n_DA=1,1 !seq_len-window_size+1


    !gd_step_size is the minimisation step size for Gradient Descent, adjusted during the minimisation, double the step size when lower cost function is achieved as long as it is smaller than the gd_max_step_size, shrink step size when it fails.
    gd_initial_step_size=0.0010000000000000    !this may use to calculate the number of minimizations, not used now
    gd_step_size= 0.0010000000000000
    gd_max_step_size=0.0100000000000000




    !!!!!output observations for diagnostic if needed
    !do i=1, window_size
    !    write(*,*) pda_ens_handle%copies(i,1:3)
    !end do



    !calculate mismatch cost function---------------------------------------------------

    call cal_mismatch_cost_function(pda_ens_handle,window_size,mis_cost)
    mis_cost_previous=mis_cost

    write(*,*) mis_cost


    !------------------------------------------------------------------------------------



    !!GD minimisation update-------------------------------------------------------------


    call init_ensemble_manager(ens_update, window_size, model_size, 1)!, transpose_type_in=2)
    call init_ensemble_manager(ens_update_copy, window_size, model_size, 1)!, transpose_type_in=2)

    !duplicate_time=.True.
    !call duplicate_ens(pda_ens_handle, ens_update_copy, duplicate_time)
    !call duplicate_ens(pda_ens_handle, ens_update, duplicate_time)
    ens_update%copies = pda_ens_handle%copies
    ens_update%time = pda_ens_handle%time
    ens_update_copy%copies = pda_ens_handle%copies
    ens_update_copy%time = pda_ens_handle%time


    !The criteria of stoping the minimisation could be based on the mismatch cost function
    GD_runs: do j=1,n_GD

        Gradient_descent: do i=1, window_size

            if (i==1) then
                forward_vector=ens_update_copy%copies(1,:)
                call get_ensemble_time(ens_update_copy, 1, ens_time)
                call get_ensemble_time(ens_update_copy, 2, target_time)

                if (.not. free_descent) then
                    call cal_adjoint(forward_vector,ens_time,target_time,time_step,model_size,adjoint)
                endif


                call cal_model_forward(forward_vector,ens_time,target_time,time_step)

                mismatch_err=ens_update_copy%copies(2,:)-forward_vector
                


                if (free_descent) then
                    ens_update%copies(1,:)=ens_update_copy%copies(1,:)+mismatch_err*gd_step_size
                else
                    do i_row=1, model_size
                        rtemp=0.00_r8
                        do i_col=1, model_size
                            rtemp=rtemp+mismatch_err(i_col)*adjoint(i_col,i_row)
                        end do
                        ens_update%copies(1,i_row)=ens_update_copy%copies(1,i_row)+rtemp*gd_step_size
                    end do
                endif
            else if (i<window_size) then
                call get_ensemble_time(ens_update_copy, i-1, ens_time)
                call get_ensemble_time(ens_update_copy, i, target_time)
                forward_vector=ens_update_copy%copies(i-1,:)

                call cal_model_forward(forward_vector,ens_time,target_time,time_step)

                mismatch_err=ens_update_copy%copies(i,:)-forward_vector
                ens_update%copies(i,:)=ens_update_copy%copies(i,:)-mismatch_err*gd_step_size

                forward_vector=ens_update_copy%copies(i,:)
                call get_ensemble_time(ens_update_copy, i, ens_time)
                call get_ensemble_time(ens_update_copy, i+1, target_time)

                if (.not. free_descent) then
                    call cal_adjoint(forward_vector,ens_time,target_time,time_step,model_size,adjoint)
                endif

                call cal_model_forward(forward_vector,ens_time,target_time,time_step)

                mismatch_err=ens_update_copy%copies(i+1,:)-forward_vector

                if (free_descent) then
                    ens_update%copies(i,:)=ens_update%copies(i,:)+mismatch_err*gd_step_size
                else
                    do i_row=1, model_size
                        rtemp=0.00_r8
                        do i_col=1, model_size
                            rtemp=rtemp+mismatch_err(i_col)*adjoint(i_col,i_row)
                        end do
                            ens_update%copies(i,i_row)=ens_update%copies(i,i_row)+rtemp*gd_step_size
                    end do
                endif

            else
                forward_vector=ens_update_copy%copies(window_size-1,:)
                call get_ensemble_time(ens_update_copy, window_size-1, ens_time)
                call get_ensemble_time(ens_update_copy, window_size, target_time)

                call cal_model_forward(forward_vector,ens_time,target_time,time_step)

                mismatch_err=ens_update_copy%copies(window_size,:)-forward_vector
                ens_update%copies(window_size,:)=ens_update_copy%copies(window_size,:)-mismatch_err*gd_step_size

            endif

        end do Gradient_descent


        !!calculate mismatch cost function for updated sequence state vectors

        call cal_mismatch_cost_function(ens_update,window_size,mis_cost)

        if (mis_cost<=mis_cost_previous) then
            !increase the GD minimisation time step by a factor of 2,
            !should have some upper bound for the GD minimisation time step

            ens_update_copy%copies=ens_update%copies

            if (gd_step_size<gd_max_step_size) then
                gd_step_size=gd_step_size*2
            endif
            mis_cost_previous=mis_cost

        else

            !decrease the GD minimisation time step by a factor of 2
            !should have some lower bound for the GD minimisation time step
            ens_update%copies=ens_update_copy%copies
            gd_step_size=gd_step_size/2
            gd_max_step_size=gd_step_size

        endif
        write(*,*) mis_cost

    end do GD_runs


    !--------------------------------------------------------------------------------------------------

    !!!output pda final update


    do i=1, window_size
        write(*,*) ens_update%copies(i,1:3)
    end do

    !write(*,*) ens_update%vars(1:3,window_size)

    call end_ensemble_manager(pda_ens_handle)
    call end_ensemble_manager(ens_update)
    call end_ensemble_manager(ens_update_copy)

end do Sequential_PDA


end subroutine pda_main

!-----------------------------------------------------------------------


subroutine cal_mismatch_cost_function(ens_handle,window_size,mis_cost)
type(ensemble_type), intent(in)   :: ens_handle
integer, intent(in)   :: window_size
type(time_type)  :: ens_time, target_time, time_step
real(r8), intent(out)  :: mis_cost
real(r8), allocatable   :: state_vector(:), mismatch_err(:), forward_vector(:)
integer :: i

allocate(state_vector(ens_handle%num_vars))
allocate(mismatch_err(ens_handle%num_vars))
allocate(forward_vector(ens_handle%num_vars))

time_step = get_model_time_step()
mismatch_err=0.0_r8
mis_cost=0.0_r8

Mismatch_error: do i=1, window_size-1

call get_ensemble_time(ens_handle, i, ens_time)
call get_ensemble_time(ens_handle, i+1, target_time)

!call print_time(ens_time, ' ens_time: ')
!call print_time(time_step, ' time_step: ')


state_vector=ens_handle%copies(i,:)
forward_vector=state_vector

call cal_model_forward(forward_vector,ens_time,target_time,time_step)

mismatch_err=ens_handle%copies(i+1,:)-forward_vector

mismatch_err=mismatch_err*mismatch_err   !here should have some sort of scalar matrix

mis_cost=sum(mismatch_err)+mis_cost


end do Mismatch_error

mis_cost=mis_cost/((window_size-1)*1.0_r8)

end subroutine cal_mismatch_cost_function

!--------------------------------------------------

subroutine cal_model_forward(model_state,ens_time,target_time,time_step)
real(r8), intent(inout)   :: model_state(:)
type(time_type), intent(inout)         :: ens_time, target_time, time_step

do while(ens_time < target_time)
call adv_1step(model_state, ens_time)
ens_time = ens_time + time_step
end do

end subroutine cal_model_forward

subroutine cal_adjoint(model_state,ens_time,target_time,time_step,model_size,adjoint)

!! this adjoint is only for 2nd order RK-------------------

real(r8), intent(in)   :: model_state(:)
type(time_type), intent(in)         :: ens_time, target_time, time_step
integer(i8), 	 intent(in) 	    :: model_size
real(r8),        intent(inout)      :: adjoint(:,:)

integer :: i,j,nstep,i_row,i_col
real(r8) :: deltat = 0.01_r8
real(r8),        pointer :: jacob(:, :)
real(r8), 	 pointer :: adjoint_1step(:,:), adjoint_temp(:,:)
real(r8) :: x1(3), dx(3) , x(3)
type(time_type)	:: t_start

allocate(jacob (model_size, model_size))
allocate(adjoint_1step (model_size, model_size))
allocate(adjoint_temp (model_size, model_size))


nstep=1
t_start=ens_time
x=model_state


do while(t_start < target_time)

if (nstep>1) call adv_1step(x, ens_time)

do i=1,model_size
do j=1,model_size
adjoint_1step(i,j)=0.000000000
end do
end do

do i=1,model_size
adjoint_1step(i,i)=1.00_r8
end do

call cal_jacob(x,jacob)

do i=1,model_size
do j=1,model_size
adjoint_1step(i,j)=adjoint_1step(i,j)+deltat*jacob(i,j)/2
end do
end do

call comp_dt(x, dx)
x1 = x + deltat * dx

call cal_jacob(x1,jacob)


do i=1,model_size
do j=1,model_size
adjoint_1step(i,j)=adjoint_1step(i,j)+deltat*jacob(i,j)/2
end do
end do


!linear propagate multiplication

if (nstep==1) then
adjoint=adjoint_1step
adjoint_temp=adjoint_1step
else

do i_row=1,model_size
do i_col=1,model_size
adjoint(i_row,i_col)=0.00_r8
do i=1,model_size
adjoint(i_row,i_col)=adjoint(i_row,i_col)+adjoint_temp(i_row,i)*adjoint_1step(i,i_col)
end do
end do
end do
adjoint_temp=adjoint

endif

nstep=nstep+1
t_start = t_start + time_step

end do


end subroutine cal_adjoint


subroutine cal_jacob(model_state,jacob)

real(r8), intent(in)   :: model_state(:)
real(r8),        intent(inout)      :: jacob(:,:)

real(r8) ::  sigma = 10.0_r8
real(r8) ::      r = 28.0_r8
real(r8) ::      b = 8.0_r8 / 3.0_r8

jacob(1,1) = -sigma
jacob(1,2) = sigma
jacob(1,3) = 0.00_r8
jacob(2,1) = r-model_state(3)
jacob(2,2) = -1
jacob(2,3) = -model_state(1)
jacob(3,1) = model_state(2)
jacob(3,2) = model_state(1)
jacob(3,3) = -b

end subroutine cal_jacob


!-----------------------------------------------------------
!> This generates the copy meta data for the diagnostic files.
!> And also creates the state space diagnostic file.
!> Note for the state space diagnostic files the order of copies
!> in the diagnostic file is different from the order of copies
!> in the ensemble handle.
subroutine filter_generate_copy_meta_data(seq, prior_inflate, PriorStateUnit, &
   PosteriorStateUnit, in_obs_copy, output_state_mean_index, &
   output_state_spread_index, prior_obs_mean_index, posterior_obs_mean_index, &
   prior_obs_spread_index, posterior_obs_spread_index)

type(obs_sequence_type),     intent(inout) :: seq
type(adaptive_inflate_type), intent(in)    :: prior_inflate
type(netcdf_file_type),      intent(inout) :: PriorStateUnit, PosteriorStateUnit
integer,                     intent(out)   :: output_state_mean_index, output_state_spread_index
integer,                     intent(in)    :: in_obs_copy
integer,                     intent(out)   :: prior_obs_mean_index, posterior_obs_mean_index
integer,                     intent(out)   :: prior_obs_spread_index, posterior_obs_spread_index

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

character(len=metadatalength) :: prior_meta_data, posterior_meta_data
! The 4 is for ensemble mean and spread plus inflation mean and spread
! The Prior file contains the prior inflation mean and spread only
! Posterior file contains the posterior inflation mean and spread only
character(len=metadatalength) :: state_meta(num_output_state_members + 4)
integer :: i, ensemble_offset, num_state_copies, num_obs_copies

! Section for state variables + other generated data stored with them.

! Ensemble mean goes first 
num_state_copies = num_output_state_members + 2
output_state_mean_index = 1
state_meta(output_state_mean_index) = 'ensemble mean'

! Ensemble spread goes second
output_state_spread_index = 2
state_meta(output_state_spread_index) = 'ensemble spread'

! Check for too many output ensemble members
if(num_output_state_members > 10000) then
   write(msgstring, *)'output metadata in filter needs state ensemble size < 10000, not ', &
                      num_output_state_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Compute starting point for ensemble member output
ensemble_offset = 2

! Set up the metadata for the output state diagnostic files
do i = 1, num_output_state_members
   write(state_meta(i + ensemble_offset), '(a15, 1x, i6)') 'ensemble member', i
end do

! Next two slots are for inflation mean and sd metadata
! To avoid writing out inflation values to the Prior and Posterior netcdf files,
! set output_inflation to false in the filter section of input.nml 
if(output_inflation) then
   num_state_copies = num_state_copies + 2
   state_meta(num_state_copies-1) = 'inflation mean'
   state_meta(num_state_copies)   = 'inflation sd'
endif

! Set up diagnostic output for model state
! All task call init and finalize diag_output.  The choice can then be made 
! in state_space_diag_mod to use a collective call (e.g. pnetcdf) or not.
PriorStateUnit     = init_diag_output('Prior_Diag', &
                     'prior ensemble state', num_state_copies, state_meta)
PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                     'posterior ensemble state', num_state_copies, state_meta)

! Set the metadata for the observations.

! Set up obs ensemble mean
num_obs_copies = in_obs_copy
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_mean_index = num_obs_copies
num_obs_copies = num_obs_copies + 1
posterior_meta_data = 'posterior ensemble mean'
call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
posterior_obs_mean_index = num_obs_copies 

! Set up obs ensemble spread 
num_obs_copies = num_obs_copies + 1
prior_meta_data = 'prior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
prior_obs_spread_index = num_obs_copies
num_obs_copies = num_obs_copies + 1
posterior_meta_data = 'posterior ensemble spread'
call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
posterior_obs_spread_index = num_obs_copies

! Make sure there are not too many copies requested
if(num_output_obs_members > 10000) then
   write(msgstring, *)'output metadata in filter needs obs ensemble size < 10000, not ',&
                      num_output_obs_members
   call error_handler(E_ERR,'filter_generate_copy_meta_data',msgstring,source,revision,revdate)
endif

! Set up obs ensemble members as requested
do i = 1, num_output_obs_members
   write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
   write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   num_obs_copies = num_obs_copies + 1
   call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
end do


end subroutine filter_generate_copy_meta_data

!-------------------------------------------------------------------------

subroutine filter_initialize_modules_used()

! Initialize modules used that require it
call initialize_mpi_utilities('Filter')

call register_module(source,revision,revdate)

! Initialize the obs sequence module
call static_init_obs_sequence()

! Initialize the model class data now that obs_sequence is all set up
call trace_message('Before init_model call')
call static_init_assim_model()
call trace_message('After  init_model call')
call state_vector_io_init()
call trace_message('After  init_state_vector_io call')
call initialize_qc()
call trace_message('After  initialize_qc call')

end subroutine filter_initialize_modules_used

!-------------------------------------------------------------------------

subroutine filter_setup_obs_sequence(seq, in_obs_copy, obs_val_index, &
   input_qc_index, DART_qc_index)

type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(out)   :: in_obs_copy, obs_val_index
integer,                 intent(out)   :: input_qc_index, DART_qc_index

character(len = metadatalength) :: no_qc_meta_data = 'No incoming data QC'
character(len = metadatalength) :: dqc_meta_data   = 'DART quality control'
character(len = 129) :: obs_seq_read_format
integer              :: obs_seq_file_id, num_obs_copies
integer              :: tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, qc_num_inc, num_qc
logical              :: pre_I_format

! Determine the number of output obs space fields
! 4 is for prior/posterior mean and spread, 
! Prior and posterior values for all selected fields (so times 2)
num_obs_copies = 2 * num_output_obs_members + 4

! Input file can have one qc field, none, or more.  note that read_obs_seq_header
! does NOT return the actual metadata values, which would be helpful in trying
! to decide if we need to add copies or qcs.
call read_obs_seq_header(obs_sequence_in_name, tnum_copies, tnum_qc, tnum_obs, tmax_num_obs, &
   obs_seq_file_id, obs_seq_read_format, pre_I_format, close_the_file = .true.)


! if there are less than 2 incoming qc fields, we will need
! to make at least 2 (one for the dummy data qc and one for
! the dart qc).
if (tnum_qc < 2) then
   qc_num_inc = 2 - tnum_qc
else
   qc_num_inc = 0
endif

! Read in with enough space for diagnostic output values and add'l qc field(s)
call read_obs_seq(obs_sequence_in_name, num_obs_copies, qc_num_inc, 0, seq)

! check to be sure that we have an incoming qc field.  if not, look for
! a blank qc field
input_qc_index = get_obs_qc_index(seq)
if (input_qc_index < 0) then
   input_qc_index = get_blank_qc_index(seq)
   if (input_qc_index < 0) then
      ! Need 1 new qc field for dummy incoming qc
      call add_qc(seq, 1)
      input_qc_index = get_blank_qc_index(seq)
      if (input_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', &
            source, revision, revdate)
      endif
   endif
   ! Since we are constructing a dummy QC, label it as such
   call set_qc_meta_data(seq, input_qc_index, no_qc_meta_data)
endif

! check to be sure we either find an existing dart qc field and
! reuse it, or we add a new one.
DART_qc_index = get_obs_dartqc_index(seq)
if (DART_qc_index < 0) then
   DART_qc_index = get_blank_qc_index(seq)
   if (DART_qc_index < 0) then
      ! Need 1 new qc field for the DART quality control
      call add_qc(seq, 1)
      DART_qc_index = get_blank_qc_index(seq)
      if (DART_qc_index < 0) then
         call error_handler(E_ERR,'filter_setup_obs_sequence', &
           'error adding blank qc field to sequence; should not happen', &
            source, revision, revdate)
      endif
   endif
   call set_qc_meta_data(seq, DART_qc_index, dqc_meta_data)
endif

! Get num of obs copies and num_qc
num_qc = get_num_qc(seq)
in_obs_copy = get_num_copies(seq) - num_obs_copies

! Create an observation type temporary for use in filter
call init_obs(observation, get_num_copies(seq), num_qc)

! Set initial DART quality control to 0 for all observations?
! Or leave them uninitialized, since
! obs_space_diagnostics should set them all without reading them

! Determine which copy has actual obs
obs_val_index = get_obs_copy_index(seq)

end subroutine filter_setup_obs_sequence

!-------------------------------------------------------------------------

function get_obs_copy_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_copy_index

integer :: i

! Determine which copy in sequence has actual obs

do i = 1, get_num_copies(seq)
   get_obs_copy_index = i
   ! Need to look for 'observation'
   if(index(get_copy_meta_data(seq, i), 'observation') > 0) return
end do
! Falling of end means 'observations' not found; die
call error_handler(E_ERR,'get_obs_copy_index', &
   'Did not find observation copy with metadata "observation"', &
      source, revision, revdate)

end function get_obs_copy_index

!-------------------------------------------------------------------------

function get_obs_prior_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_prior_index

integer :: i

! Determine which copy in sequence has prior mean, if any.

do i = 1, get_num_copies(seq)
   get_obs_prior_index = i
   ! Need to look for 'prior mean'
   if(index(get_copy_meta_data(seq, i), 'prior ensemble mean') > 0) return
end do
! Falling of end means 'prior mean' not found; not fatal!

get_obs_prior_index = -1

end function get_obs_prior_index

!-------------------------------------------------------------------------

function get_obs_qc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_qc_index

integer :: i

! Determine which qc, if any, has the incoming obs qc
! this is tricky because we have never specified what string
! the metadata has to have.  look for 'qc' or 'QC' and the
! first metadata that matches (much like 'observation' above)
! is the winner.

do i = 1, get_num_qc(seq)
   get_obs_qc_index = i

   ! Need to avoid 'QC metadata not initialized'
   if(index(get_qc_meta_data(seq, i), 'QC metadata not initialized') > 0) cycle
  
   ! Need to look for 'QC' or 'qc'
   if(index(get_qc_meta_data(seq, i), 'QC') > 0) return
   if(index(get_qc_meta_data(seq, i), 'qc') > 0) return
   if(index(get_qc_meta_data(seq, i), 'Quality Control') > 0) return
   if(index(get_qc_meta_data(seq, i), 'QUALITY CONTROL') > 0) return
end do
! Falling off end means 'QC' string not found; not fatal!

get_obs_qc_index = -1

end function get_obs_qc_index

!-------------------------------------------------------------------------

function get_obs_dartqc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_obs_dartqc_index

integer :: i

! Determine which qc, if any, has the DART qc

do i = 1, get_num_qc(seq)
   get_obs_dartqc_index = i
   ! Need to look for 'DART quality control'
   if(index(get_qc_meta_data(seq, i), 'DART quality control') > 0) return
end do
! Falling off end means 'DART quality control' not found; not fatal!

get_obs_dartqc_index = -1

end function get_obs_dartqc_index

!-------------------------------------------------------------------------

function get_blank_qc_index(seq)

type(obs_sequence_type), intent(in) :: seq
integer                             :: get_blank_qc_index

integer :: i

! Determine which qc, if any, is blank

do i = 1, get_num_qc(seq)
   get_blank_qc_index = i
   ! Need to look for 'QC metadata not initialized'
   if(index(get_qc_meta_data(seq, i), 'QC metadata not initialized') > 0) return
end do
! Falling off end means unused slot not found; not fatal!

get_blank_qc_index = -1

end function get_blank_qc_index

!-------------------------------------------------------------------------

subroutine filter_set_initial_time(days, seconds, time, read_time_from_file)

integer,         intent(in)  :: days, seconds
type(time_type), intent(out) :: time
logical,         intent(out) :: read_time_from_file

if(days >= 0) then
   time = set_time(seconds, days)
   read_time_from_file = .false.
else
   time = set_time(0, 0)
   read_time_from_file = .true.
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------

subroutine filter_set_window_time(time)

type(time_type), intent(out) :: time


if(obs_window_days >= 0) then
   time = set_time(obs_window_seconds, obs_window_days)
else
   time = set_time(0, 0)
endif

end subroutine filter_set_window_time

!-------------------------------------------------------------------------

subroutine filter_ensemble_inflate(ens_handle, inflate_copy, inflate, ENS_MEAN_COPY)

type(ensemble_type),         intent(inout) :: ens_handle
integer,                     intent(in)    :: inflate_copy, ENS_MEAN_COPY
type(adaptive_inflate_type), intent(inout) :: inflate

integer :: j, group, grp_bot, grp_top, grp_size

! Assumes that the ensemble is copy complete
call prepare_to_update_copies(ens_handle)

! Inflate each group separately;  Divide ensemble into num_groups groups
grp_size = ens_size / num_groups

do group = 1, num_groups
   grp_bot = (group - 1) * grp_size + 1
   grp_top = grp_bot + grp_size - 1
   ! Compute the mean for this group
   call compute_copy_mean(ens_handle, grp_bot, grp_top, ENS_MEAN_COPY)

   do j = 1, ens_handle%my_num_vars
      call inflate_ens(inflate, ens_handle%copies(grp_bot:grp_top, j), &
         ens_handle%copies(ENS_MEAN_COPY, j), ens_handle%copies(inflate_copy, j))
   end do
end do

end subroutine filter_ensemble_inflate

!-------------------------------------------------------------------------

subroutine obs_space_diagnostics(obs_fwd_op_ens_handle, qc_ens_handle, ens_size, &
   seq, keys, prior_post, num_output_members, members_index, &
   obs_val_index, OBS_KEY_COPY, &
   ens_mean_index, ens_spread_index, num_obs_in_set, &
   OBS_MEAN_START, OBS_VAR_START, OBS_GLOBAL_QC_COPY, OBS_VAL_COPY, &
   OBS_ERR_VAR_COPY, DART_qc_index)

! Do prior observation space diagnostics on the set of obs corresponding to keys

type(ensemble_type),     intent(inout) :: obs_fwd_op_ens_handle, qc_ens_handle
integer,                 intent(in)    :: ens_size
integer,                 intent(in)    :: num_obs_in_set
integer,                 intent(in)    :: keys(num_obs_in_set), prior_post
integer,                 intent(in)    :: num_output_members, members_index
integer,                 intent(in)    :: obs_val_index
integer,                 intent(in)    :: OBS_KEY_COPY
integer,                 intent(in)    :: ens_mean_index, ens_spread_index
type(obs_sequence_type), intent(inout) :: seq
integer,                 intent(in)    :: OBS_MEAN_START, OBS_VAR_START
integer,                 intent(in)    :: OBS_GLOBAL_QC_COPY, OBS_VAL_COPY
integer,                 intent(in)    :: OBS_ERR_VAR_COPY, DART_qc_index

integer               :: j, k, ens_offset
integer               :: ivalue
real(r8), allocatable :: obs_temp(:)
real(r8)              :: rvalue(1)

! Do verbose forward operator output if requested
if(output_forward_op_errors) call verbose_forward_op_output(qc_ens_handle, prior_post, ens_size, keys)

! Make var complete for get_copy() calls below.
! Can you use a gather instead of a transpose and get copy?
call all_copies_to_all_vars(obs_fwd_op_ens_handle)

! allocate temp space for sending data - surely only task 0 needs to allocate this?
allocate(obs_temp(num_obs_in_set))

! Update the ensemble mean
! Get this copy to process 0
call get_copy(map_task_to_pe(obs_fwd_op_ens_handle, 0), obs_fwd_op_ens_handle, OBS_MEAN_START, obs_temp) 
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
     ! Loop through the observations for this time
     do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_obs_values(seq, keys(j), rvalue, ens_mean_index)
     end do
  endif

! Update the ensemble spread
! Get this copy to process 0
call get_copy(map_task_to_pe(obs_fwd_op_ens_handle, 0), obs_fwd_op_ens_handle, OBS_VAR_START, obs_temp)
! Only pe 0 gets to write the sequence
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_fwd_op_ens_handle%num_vars
      ! update the spread in each obs
      if (obs_temp(j) /= missing_r8) then
         rvalue(1) = sqrt(obs_temp(j))
      else
         rvalue(1) = obs_temp(j)
      endif
      call replace_obs_values(seq, keys(j), rvalue, ens_spread_index)
   end do
endif

! May be possible to only do this after the posterior call...
! Update any requested ensemble members
ens_offset = members_index + 4
! Update all of these ensembles that are required to sequence file
do k = 1, num_output_members
   ! Get this copy on pe 0
   call get_copy(map_task_to_pe(obs_fwd_op_ens_handle, 0), obs_fwd_op_ens_handle, k, obs_temp)
   ! Only task 0 gets to write the sequence
   if(my_task_id() == 0) then
      ! Loop through the observations for this time
      do j = 1, obs_fwd_op_ens_handle%num_vars
         ! update the obs values 
         rvalue(1) = obs_temp(j)
         ivalue = ens_offset + 2 * (k - 1)
         call replace_obs_values(seq, keys(j), rvalue, ivalue)
      end do
   endif
end do

! Update the qc global value
call get_copy(map_task_to_pe(obs_fwd_op_ens_handle, 0), obs_fwd_op_ens_handle, OBS_GLOBAL_QC_COPY, obs_temp)
! Only task 0 gets to write the observations for this time
if(my_task_id() == 0) then
   ! Loop through the observations for this time
   do j = 1, obs_fwd_op_ens_handle%num_vars
      rvalue(1) = obs_temp(j)
      call replace_qc(seq, keys(j), rvalue, DART_qc_index)
   end do
endif

! clean up.
deallocate(obs_temp)

end subroutine obs_space_diagnostics

!-------------------------------------------------------------------------

subroutine filter_sync_keys_time(ens_handle, key_bounds, num_obs_in_set, time1, time2)

integer,             intent(inout)  :: key_bounds(2), num_obs_in_set
type(time_type),     intent(inout)  :: time1, time2
type(ensemble_type), intent(inout)     :: ens_handle

! Have owner of copy 1 broadcast these values to all other tasks.
! Only tasks which contain copies have this info; doing it this way
! allows ntasks > nens to work.

real(r8) :: rkey_bounds(2), rnum_obs_in_set(1)
real(r8) :: rtime(4)
integer  :: days, secs
integer  :: copy1_owner, owner_index

call get_copy_owner_index(1, copy1_owner, owner_index)

if( ens_handle%my_pe == copy1_owner) then
   rkey_bounds = key_bounds
   rnum_obs_in_set(1) = num_obs_in_set
   call get_time(time1, secs, days)
   rtime(1) = secs
   rtime(2) = days
   call get_time(time2, secs, days)
   rtime(3) = secs
   rtime(4) = days
   call broadcast_send(map_pe_to_task(ens_handle, copy1_owner), rkey_bounds, rnum_obs_in_set, rtime)
else
   call broadcast_recv(map_pe_to_task(ens_handle, copy1_owner), rkey_bounds, rnum_obs_in_set, rtime)
   key_bounds =     nint(rkey_bounds)
   num_obs_in_set = nint(rnum_obs_in_set(1))
   time1 = set_time(nint(rtime(1)), nint(rtime(2)))
   time2 = set_time(nint(rtime(3)), nint(rtime(4)))
endif

! Every task gets the current time (necessary for the forward operator)
ens_handle%current_time = time1

end subroutine filter_sync_keys_time

!-------------------------------------------------------------------------
! Only copy 1 on task zero has the correct time after reading
! when you read one instance using filter_read_restart. 
! perturb_from_single_instance = .true.
! This routine makes the times consistent across the ensemble.  
! Any task that owns one or more state vectors needs the time for 
! the move ahead call. 
!> @todo This is broadcasting the time to all tasks, not
!> just the tasks that own copies.
subroutine broadcast_time_across_copy_owners(ens_handle, ens_time)

type(ensemble_type), intent(inout) :: ens_handle
type(time_type),     intent(in)    :: ens_time

real(r8) :: rtime(2)
integer  :: days, secs
integer  :: copy1_owner, owner_index
type(time_type) :: time_from_copy1

call get_copy_owner_index(1, copy1_owner, owner_index)

if( ens_handle%my_pe == copy1_owner) then
   call get_time(ens_time, secs, days)
   rtime(1) = secs
   rtime(2) = days
   call broadcast_send(map_pe_to_task(ens_handle, copy1_owner), rtime)
   ens_handle%time(1:ens_handle%my_num_copies) = ens_time
else
   call broadcast_recv(map_pe_to_task(ens_handle, copy1_owner), rtime)
   time_from_copy1 = set_time(nint(rtime(1)), nint(rtime(2)))
   if (ens_handle%my_num_copies > 0) ens_handle%time(1:ens_handle%my_num_copies) = time_from_copy1
endif

end subroutine broadcast_time_across_copy_owners

!-------------------------------------------------------------------------

subroutine set_trace(trace_execution, output_timestamps, silence)

logical, intent(in) :: trace_execution
logical, intent(in) :: output_timestamps
logical, intent(in) :: silence

! Set whether other modules trace execution with messages
! and whether they output timestamps to trace overall performance

! defaults
trace_level     = 0
timestamp_level = 0

! selectively turn stuff back on
if (trace_execution)   trace_level     = 1
if (output_timestamps) timestamp_level = 1

! turn as much off as possible
if (silence) then
   trace_level     = -1
   timestamp_level = -1
endif

call set_smoother_trace(trace_level, timestamp_level)
call set_obs_model_trace(trace_level, timestamp_level)
call set_assim_tools_trace(trace_level, timestamp_level)

end subroutine set_trace

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

subroutine timestamp_message(msg, sync)

character(len=*), intent(in) :: msg
logical, intent(in), optional :: sync

! Write current time and message to stdout and log file. 
! if sync is present and true, sync mpi jobs before printing time.

if (timestamp_level <= 0) return

if (present(sync)) then
  if (sync) call task_sync()
endif

if (do_output()) call timestamp(' '//trim(msg), pos='brief')  ! was debug

end subroutine timestamp_message

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

!-------------------------------------------------------------------------

subroutine print_obs_time(seq, key, msg)

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
character(len=*), intent(in), optional :: msg

! Write time of an observation to stdout and log file.
type(obs_type) :: obs
type(obs_def_type) :: obs_def
type(time_type) :: mtime

if (trace_level <= 0) return

if (do_output()) then
   call init_obs(obs, 0, 0)
   call get_obs_from_key(seq, key, obs)
   call get_obs_def(obs, obs_def)
   mtime = get_obs_def_time(obs_def)
   call print_time(mtime, ' filter trace: '//msg, logfileunit)
   call print_time(mtime, ' filter trace: '//msg)
   call destroy_obs(obs)
endif

end subroutine print_obs_time

!-------------------------------------------------------------------------
!> write out failed forward operators
!> This was part of obs_space_diagnostics
subroutine verbose_forward_op_output(qc_ens_handle, prior_post, ens_size, keys)

type(ensemble_type), intent(inout) :: qc_ens_handle
integer,             intent(in)    :: prior_post
integer,             intent(in)    :: ens_size
integer,             intent(in)    :: keys(:) ! I think this is still var size

character*12 :: task
integer :: j, i
integer :: forward_unit

write(task, '(i6.6)') my_task_id()

! all tasks open file?
if(prior_post == PRIOR_DIAG) then
   forward_unit = open_file('prior_forward_ope_errors' // task, 'formatted', 'append')
else
   forward_unit = open_file('post_forward_ope_errors' // task, 'formatted', 'append')
endif

! qc_ens_handle is a real representing an integer; values /= 0 get written out
do i = 1, ens_size
   do j = 1, qc_ens_handle%my_num_vars
      if(nint(qc_ens_handle%copies(i, j)) /= 0) write(forward_unit, *) i, keys(j), nint(qc_ens_handle%copies(i, j))
   end do
end do

call close_file(forward_unit)

end subroutine verbose_forward_op_output

!------------------------------------------------------------------
!> Produces an ensemble by copying my_vars of the 1st ensemble member
!> and then perturbing the copies array.
!> Mimicks the behaviour of pert_model_state: 
!> pert_model_copies is called:
!>   if no model perturb is provided, perturb_copies_task_bitwise is called.
!> Note: Not enforcing a model_mod to produce a 
!> pert_model_copies that is bitwise across any number of
!> tasks, although there is enough information in the 
!> ens_handle to do this.
!>
!> Some models allow missing_r8 in the state vector.  If missing_r8 is 
!> allowed the locations of missing_r8s are stored before the perturb, 
!> then the missing_r8s are put back in after the perturb.
subroutine create_ensemble_from_single_file(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

integer               :: i ! loop variable
logical               :: interf_provided ! model does the perturbing
logical, allocatable  :: miss_me(:)

! Copy from ensemble member 1 to the other copies
do i = 1, ens_handle%my_num_vars
   ens_handle%copies(2:ens_size, i) = ens_handle%copies(1, i)  ! How slow is this?
enddo

! store missing_r8 locations
if (get_missing_ok_status()) then ! missing_r8 is allowed in the state
   allocate(miss_me(ens_size))
   miss_me = .false.
   where(ens_handle%copies(1, :) == missing_r8) miss_me = .true.
endif

call pert_model_copies(ens_handle, ens_size, perturbation_amplitude, interf_provided)
if (.not. interf_provided) then
   call perturb_copies_task_bitwise(ens_handle)
endif

! Put back in missing_r8
if (get_missing_ok_status()) then
   do i = 1, ens_size
      where(miss_me) ens_handle%copies(i, :) = missing_r8
   enddo
endif

end subroutine create_ensemble_from_single_file


!------------------------------------------------------------------
! Perturb the copies array in a way that is bitwise reproducible 
! no matter how many task you run on.
subroutine perturb_copies_task_bitwise(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

integer               :: i, j ! loop variables
type(random_seq_type) :: r(ens_size)
real(r8)              :: random_number(ens_size) ! array of random numbers
integer               :: local_index

! Need ens_size random number sequences.
do i = 1, ens_size
   call init_random_seq(r(i), i)
enddo

local_index = 1 ! same across the ensemble

! Only one task is going to update per i.  This will not scale at all.
do i = 1, ens_handle%num_vars

   do j = 1, ens_size
     ! Can use %copies here because the random number
     ! is only relevant to the task than owns element i.
     random_number(j)  =  random_gaussian(r(j), ens_handle%copies(j, local_index), perturbation_amplitude)
   enddo

   if (ens_handle%my_vars(local_index) == i) then
      ens_handle%copies(1:ens_size, local_index) = random_number(:)
      local_index = local_index + 1 ! task is ready for the next random number
      local_index = min(local_index, ens_handle%my_num_vars)
   endif

enddo

end subroutine perturb_copies_task_bitwise

!------------------------------------------------------------------
!> Set the time on any extra copies that a pe owns
!> Could we just set the time on all copies?
subroutine set_time_on_extra_copies(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

integer :: copy_num, owner, owners_index
integer :: ens_size

ens_size = ens_handle%num_copies - ens_handle%num_extras

do copy_num = ens_size + 1, ens_handle%num_copies
   ! Set time for a given copy of an ensemble
   call get_copy_owner_index(copy_num, owner, owners_index)
   if(ens_handle%my_pe == owner) then
      call set_ensemble_time(ens_handle, owners_index, ens_handle%current_time)
   endif
enddo

end subroutine  set_time_on_extra_copies

!------------------------------------------------------------------
!> Copy the current mean, sd, inf_mean, inf_sd to spare copies
!> Assuming that if the spare copy is there you should fill it
subroutine store_prior(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

if (query_copy_present(SPARE_PRIOR_MEAN)) &
   ens_handle%copies(SPARE_PRIOR_MEAN, :) = ens_handle%copies(ENS_MEAN_COPY, :)

if (query_copy_present(SPARE_PRIOR_SPREAD)) &
   ens_handle%copies(SPARE_PRIOR_SPREAD, :) = ens_handle%copies(ENS_SD_COPY, :)

if (query_copy_present(SPARE_PRIOR_INF_MEAN)) &
   ens_handle%copies(SPARE_PRIOR_INF_MEAN, :) = ens_handle%copies(PRIOR_INF_COPY, :)

if (query_copy_present(SPARE_PRIOR_INF_SPREAD)) &
   ens_handle%copies(SPARE_PRIOR_INF_SPREAD, :) = ens_handle%copies(PRIOR_INF_SD_COPY, :)

end subroutine store_prior

!------------------------------------------------------------------
!> Copy the current post_inf_mean, post_inf_sd to spare copies
!> Assuming that if the spare copy is there you should fill it
!> No need to store the mean and sd as you would with store_prior because
!> mean and sd are not changed during filter_assim(inflate_only = .true.)
subroutine store_posterior(ens_handle)

type(ensemble_type), intent(inout) :: ens_handle

if (query_copy_present(SPARE_POST_INF_MEAN)) &
   ens_handle%copies(SPARE_POST_INF_MEAN, :) = ens_handle%copies(POST_INF_COPY, :)

if (query_copy_present(SPARE_POST_INF_SPREAD)) &
   ens_handle%copies(SPARE_POST_INF_SPREAD, :) = ens_handle%copies(POST_INF_SD_COPY, :)

end subroutine store_posterior

!------------------------------------------------------------------
!> Calculate how many spare copies are needed for given input options
!> and give a number to the name.
!> The indicies are initailzed to COPY_NOT_PRESENT before they
!> are set here.
!> The copies could be shuffled around depending on inflation options.
!> For now, always have the first 6 extra copies.
!> Note if you remove inflation copies, check read_state() and write_state()
!> for assumptions about which copies are present.
subroutine set_state_copies(ens_size, num_output_state_members, num_extras)

integer, intent(in)  :: num_output_state_members
integer, intent(in)  :: ens_size
integer, intent(out) :: num_extras

! State
ENS_MEAN_COPY        = ens_size + 1
ENS_SD_COPY          = ens_size + 2
PRIOR_INF_COPY       = ens_size + 3
PRIOR_INF_SD_COPY    = ens_size + 4
POST_INF_COPY        = ens_size + 5
POST_INF_SD_COPY     = ens_size + 6

num_extras = 6

! If there are no diagnostic files, we will need to store the
! copies that would have gone in Prior_Diag.nc and Posterior_Diag.nc
! in spare copies in the ensemble.
if (skip_diag_files() .and. num_output_state_members <= 0) then

   ! Not stopping to write prior_members so keep these Prior copies
   ! as extra copies and write them and the end.
   SPARE_PRIOR_MEAN       = ens_size + 7
   SPARE_PRIOR_SPREAD     = ens_size + 8
   SPARE_PRIOR_INF_MEAN   = ens_size + 9
   SPARE_PRIOR_INF_SPREAD = ens_size + 10
   ! need to store posterior inflation mean and inflation spread since
   ! these are overwritten in filter_assim(inflate_only=.true.)
   SPARE_POST_INF_MEAN    = ens_size + 11
   SPARE_POST_INF_SPREAD  = ens_size + 12
   num_extras = num_extras + 6

elseif (skip_diag_files() .and. num_output_state_members > 0) then
   ! Prior members and extra copies are written out prior to filter_assim.
   ! Need to store posterior inflation mean and inflation spread since
   ! these are overwritten in filter_assim(inflate_only=.true.)
   SPARE_POST_INF_MEAN    = ens_size + 7
   SPARE_POST_INF_SPREAD  = ens_size + 8
   num_extras = num_extras + 2
endif

end subroutine set_state_copies

!==================================================================
! TEST FUNCTIONS BELOW THIS POINT
!------------------------------------------------------------------
!> dump out obs_copies to file
subroutine test_obs_copies(obs_fwd_op_ens_handle, information)

type(ensemble_type), intent(in) :: obs_fwd_op_ens_handle
character(len=*),    intent(in) :: information

character*20  :: task_str !< string to hold the task number
character*129 :: file_obscopies !< output file name
integer :: i

write(task_str, '(i10)') obs_fwd_op_ens_handle%my_pe
file_obscopies = TRIM('obscopies_' // TRIM(ADJUSTL(information)) // TRIM(ADJUSTL(task_str)))
open(15, file=file_obscopies, status ='unknown')

do i = 1, obs_fwd_op_ens_handle%num_copies - 4
   write(15, *) obs_fwd_op_ens_handle%copies(i,:)
enddo

close(15)

end subroutine test_obs_copies

!-------------------------------------------------------------------
end module filter_mod

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/pda/filter/filter_mod.f90 $
! $Id: filter_mod.f90 10343 2016-06-07 21:51:48Z hendric $
! $Revision: 10343 $
! $Date: 2016-06-07 15:51:48 -0600 (Tue, 07 Jun 2016) $
