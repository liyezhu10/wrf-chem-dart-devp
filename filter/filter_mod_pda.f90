! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: filter_mod.f90 10343 2016-06-07 21:51:48Z hendric $

module filter_mod

!------------------------------------------------------------------------------
use types_mod,             only : r8, i8, missing_r8, metadatalength
use obs_sequence_mod,      only : static_init_obs_sequence

!Du adds operator(+) operator(<)
use time_manager_mod,      only : time_type, get_time, set_time, operator(/=), operator(>),   &
                                  operator(-), operator(+), operator(<), print_time
use utilities_mod,         only : register_module,  error_handler, E_ERR, E_MSG, E_DBG,       &
                                  logfileunit, nmlfileunit, timestamp,  E_ALLMSG, &
                                  do_output, find_namelist_in_file, check_namelist_read,      &
                                  open_file, close_file, do_nml_file, do_nml_term
!Du add get_model_time_step
use assim_model_mod,       only : static_init_assim_model, get_model_size, get_model_time_step

use obs_model_mod,         only : advance_state
use ensemble_manager_mod,  only : init_ensemble_manager, end_ensemble_manager,                &
                                  ensemble_type, get_copy, get_my_num_copies, put_copy,       &
                                  all_vars_to_all_copies, all_copies_to_all_vars,             &
                                  duplicate_ens, get_copy_owner_index, &
                                  get_ensemble_time, set_ensemble_time, &
                                  prepare_to_read_from_vars, prepare_to_write_to_vars,        &
                                  prepare_to_read_from_copies,  get_my_num_vars,              &
                                  prepare_to_write_to_copies, get_ensemble_time,              &
                                  map_task_to_pe,  map_pe_to_task, prepare_to_update_copies,  &
                                  copies_in_window, set_num_extra_copies, get_allow_transpose, &
                                  all_copies_to_all_vars, allocate_single_copy,               &
                                  get_single_copy, put_single_copy, deallocate_single_copy

use mpi_utilities_mod,     only : initialize_mpi_utilities, finalize_mpi_utilities,           &
                                  my_task_id, task_sync, broadcast_send, broadcast_recv,      &
                                  task_count, sum_across_tasks

use state_vector_io_mod,   only : state_vector_io_init, read_state, write_state

use io_filenames_mod,      only : io_filenames_init, file_info_type

use quality_control_mod,   only : initialize_qc

use state_space_diag_mod,  only : filter_state_space_diagnostics, netcdf_file_type, &
                                  init_diag_output, finalize_diag_output,           &
                                  skip_diag_files

use model_mod, 		  only  : adv_1step, get_state_meta_data
use location_mod,     only  : location_type
use obs_kind_mod,     only  : get_raw_obs_kind_name


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
!> The code does pda data assimilation with free adjoint.


subroutine pda_main()

type(time_type)             :: time1

integer(i8)             :: model_size, var_ind
integer                 :: i, iunit, io
! Global indices into ensemble storage - observations
integer                 :: input_qc_index, DART_qc_index
logical                 :: read_time_from_file

integer :: num_extras ! the extra ensemble copies

type(file_info_type) :: file_info

logical                 :: all_gone, allow_missing



!!DuDu adds....
type(ensemble_type)         :: pda_ens_handle,  forward_ens_handle, ens_update_copy,ens_normalization
logical                 :: duplicate_time, free_descent
integer                 :: j, k, seq_len, n_GD, window_size, n_DA, i_row, i_col, var_type
real(r8)                :: value(1), mis_cost, mis_cost_previous, gd_step_size, gd_max_step_size,gd_initial_step_size,sum_variable,result_sum
real(r8), allocatable   :: state_vector(:), mismatch_err(:), forward_vector(:)
type(time_type) 	    :: mtime
type(time_type)         :: time_step, target_time, ens_time
real(r8),       pointer :: adjoint(:, :)
real(r8)		:: rtemp
type(location_type)     :: location

!Du number of minimisation should be setup more standardised, or use other criteria to stop the minimazation
n_GD=100

!Du define assimilation window size
window_size=120;

call filter_initialize_modules_used() ! static_init_model called in here

!!!Du  should create filterpda namelist in input.nml

! Read the namelist entry
call find_namelist_in_file("input.nml", "filter_nml", iunit)
read(iunit, nml = filter_nml, iostat = io)
call check_namelist_read(iunit, io, "filter_nml")


model_size = get_model_size()

call init_ensemble_manager(pda_ens_handle, window_size, model_size, 1, transpose_type_in=2)

!shall set pda_ens_handle%num_extras=0

!!!!below shall use namelist instead
single_restart_file_in=.false.
single_restart_file_out=.false.
restart_in_file_name='Null'
use_restart_list=.true.
restart_list_file(1)='pda_ic_name_list'
direct_netcdf_read=.true.

file_info = io_filenames_init(pda_ens_handle, single_restart_file_in, single_restart_file_out, &
restart_in_file_name, restart_out_file_name, output_restart, direct_netcdf_read, &
direct_netcdf_write, output_restart_mean, add_domain_extension, use_restart_list, &
restart_list_file, overwrite_state_input, inf_in_file_name, inf_out_file_name)

read_time_from_file=.true.


call read_state(pda_ens_handle, file_info, read_time_from_file, time1)

!! use free_descent to decide whether use model adjoint or free adjoint (Identity adjoint)
!free_descent=.true.


time_step = get_model_time_step()
allocate(mismatch_err(model_size))
call init_ensemble_manager(forward_ens_handle, window_size, model_size, 1, transpose_type_in=2)
call init_ensemble_manager(ens_update_copy, window_size, model_size, 1, transpose_type_in=2)
call init_ensemble_manager(ens_normalization, 1, model_size, 1, transpose_type_in=2)

!generate normalization vector--------------------------------------


!ens_normalization%copies=pda_ens_handle%copies
!ens_normalization%time=pda_ens_handle%time


if(.not. allocated(ens_normalization%vars)) allocate(ens_normalization%vars(ens_normalization%num_vars, ens_normalization%my_num_copies))

call all_copies_to_all_vars(ens_normalization)

if (my_task_id()==0) then
do var_ind=1,model_size




    call get_state_meta_data(pda_ens_handle, var_ind, location, var_type)

    !write(*,*) 'norm_size', shape(ens_normalization%vars(:,:))

    if (var_type==1) then
        ens_normalization%vars(var_ind,1)=2
        elseif (var_type==2) then
            ens_normalization%vars(var_ind,1)=1
        elseif (var_type==3) then
            ens_normalization%vars(var_ind,1)=60
        else
            ens_normalization%vars(var_ind,1)=2

    endif
end do
endif
call all_vars_to_all_copies(ens_normalization)

!write(*,*) ens_normalization%copies(1,:)


!------------------------------------------------------------------






!seq_len=200

! do PDA "sequentially"
Sequential_PDA: do n_DA=1,1 !seq_len-window_size+1


    !gd_step_size is the minimisation step size for Gradient Descent, adjusted during the minimisation, double the step size when lower cost function is achieved as long as it is smaller than the gd_max_step_size, shrink step size when it fails.
    gd_initial_step_size=0.0010000000000000    !this may use to calculate the number of minimizations, not used now
    gd_step_size= 0.00050000000000000
    gd_max_step_size=0.00100000000000000


    !calculate mismatch cost function---------------------------------------------------

    call get_ensemble_time(pda_ens_handle, 1, ens_time)
    call get_ensemble_time(pda_ens_handle, 2, target_time)

    forward_ens_handle%copies=pda_ens_handle%copies
    forward_ens_handle%time=pda_ens_handle%time

    ens_update_copy%copies=pda_ens_handle%copies
    ens_update_copy%time=pda_ens_handle%time

    if(.not. allocated(forward_ens_handle%vars)) allocate(forward_ens_handle%vars(forward_ens_handle%num_vars, forward_ens_handle%my_num_copies))


    call all_copies_to_all_vars(forward_ens_handle)


    !!!this is not ideal, shall be able to pass an arrary of target time to advance_state, things needs to be changed in advance_state to adapt this
    do i=1, window_size
        forward_ens_handle%time(i)=ens_time
    end do

    !call print_time(ens_time, ' ens_time: ')

    call advance_state(forward_ens_handle, window_size, target_time, async, adv_ens_command, tasks_per_model_advance)

    call all_vars_to_all_copies(forward_ens_handle)

    !calculate mismatches
    do i=1, window_size-1
        forward_ens_handle%copies(i,:)=ens_update_copy%copies(i+1,:)-forward_ens_handle%copies(i,:)
    end do

    !calculate mismatch cost function
    result_sum=0
    sum_variable=0
    do i=1,window_size-1

        do j=1,forward_ens_handle%my_num_vars
            sum_variable=forward_ens_handle%copies(i,j)*forward_ens_handle%copies(i,j)+sum_variable
        end do

        call sum_across_tasks(sum_variable,result_sum)

    end do
    mis_cost=result_sum/((window_size-1)*1.0_r8)
    mis_cost_previous=mis_cost


    !write(*,*) mis_cost



    !------------------------------------------------------------------------------------



    !!GD minimisation update-------------------------------------------------------------



    !The criteria of stoping the minimisation could be based on the mismatch cost function
    GD_runs: do j=1,n_GD

        Gradient_descent: do i=1, window_size

            if (i==1) then
                pda_ens_handle%copies(1,:)=ens_update_copy%copies(1,:)+gd_step_size*ens_normalization%copies(1,:)*forward_ens_handle%copies(1,:)

            else if (i<window_size) then

                pda_ens_handle%copies(i,:)=ens_update_copy%copies(i,:) &
                                         -gd_step_size*ens_normalization%copies(1,:)*forward_ens_handle%copies(i-1,:) &
                                         +gd_step_size*ens_normalization%copies(1,:)*forward_ens_handle%copies(i,:)

            else

                pda_ens_handle%copies(window_size,:)=ens_update_copy%copies(window_size,:) &
                            -gd_step_size*ens_normalization%copies(1,:)*forward_ens_handle%copies(window_size-1,:)

                !write(msgstring, *) 'vars=', pda_ens_handle%copies(window_size,:)
                !call error_handler(E_ALLMSG,'filter_main', msgstring, source, revision, revdate)


            endif

        end do Gradient_descent


        !!calculate mismatch cost function for updated sequence state vectors

        call get_ensemble_time(pda_ens_handle, 1, ens_time)
        call get_ensemble_time(pda_ens_handle, 2, target_time)
        forward_ens_handle%copies=pda_ens_handle%copies
        forward_ens_handle%time=pda_ens_handle%time

        if(.not. allocated(forward_ens_handle%vars)) allocate(forward_ens_handle%vars(forward_ens_handle%num_vars, forward_ens_handle%my_num_copies))
        call all_copies_to_all_vars(forward_ens_handle)

        do i=1, window_size
            forward_ens_handle%time(i)=ens_time
        end do

        call advance_state(forward_ens_handle, window_size, target_time, async, &
        adv_ens_command, tasks_per_model_advance)
        call all_vars_to_all_copies(forward_ens_handle)

        do i=1, window_size-1
            forward_ens_handle%copies(i,:)=pda_ens_handle%copies(i+1,:)-forward_ens_handle%copies(i,:)
        end do


        sum_variable=0
        do i=1,window_size-1

            do k=1,forward_ens_handle%my_num_vars
                sum_variable=forward_ens_handle%copies(i,k)*forward_ens_handle%copies(i,k)+sum_variable
            end do

            call sum_across_tasks(sum_variable,result_sum)


        end do
        mis_cost=result_sum/((window_size-1)*1.0_r8)

        write(*,*) mis_cost



        if (mis_cost<mis_cost_previous) then
            !increase the GD minimisation time step by a factor of 2,
            !should have some upper bound for the GD minimisation time step

            !call duplicate_ens(pda_ens_handle, ens_update_copy, duplicate_time)
            ens_update_copy%copies=pda_ens_handle%copies

            if (gd_step_size<gd_max_step_size) then
                gd_step_size=gd_step_size*2
            endif
            mis_cost_previous=mis_cost

        else

            !decrease the GD minimisation time step by a factor of 2
            !should have some lower bound for the GD minimisation time step
            pda_ens_handle%vars=ens_update_copy%vars
            gd_step_size=gd_step_size/2
            gd_max_step_size=gd_step_size

            !may use a vector to store previous forward ensemble copies
            call duplicate_ens(pda_ens_handle, forward_ens_handle, duplicate_time)

            do i=1, window_size
                forward_ens_handle%time(i)=ens_time
            end do

            call advance_state(forward_ens_handle, window_size, target_time, async, &
            adv_ens_command, tasks_per_model_advance)

            write(*,*) 'larger_mis cost', mis_cost


        endif
        write(*,*) mis_cost

    end do GD_runs


    !--------------------------------------------------------------------------------------------------

    !call write_state(pda_ens_handle, file_info)

    !!!output pda final update


    do i=1, window_size
        write(*,*) pda_ens_handle%copies(i,1:4)
    end do

    !write(*,*) ens_update%vars(1:3,window_size)

    call end_ensemble_manager(pda_ens_handle)
    call end_ensemble_manager(forward_ens_handle)
    call end_ensemble_manager(ens_update_copy)

end do Sequential_PDA


end subroutine pda_main

!-----------------------------------------------------------------------

subroutine cal_mismatch_cost_function(ens_handle,forward_ens_handle,window_size,model_size,mis_cost)
type(ensemble_type), intent(in)   :: ens_handle, forward_ens_handle
integer, intent(in)   :: window_size
integer(i8), intent(in)   :: model_size
real(r8), intent(out)  :: mis_cost
real(r8), allocatable   :: state_vector(:), mismatch_err(:), forward_vector(:)
integer :: i


allocate(mismatch_err(model_size))

mismatch_err=0.0_r8
mis_cost=0.0_r8

Mismatch_error: do i=1, window_size-1

mismatch_err=ens_handle%vars(:,i+1)-forward_ens_handle%vars(:,i)

mismatch_err=mismatch_err*mismatch_err   !here should have some sort of scalar matrix

mis_cost=sum(mismatch_err)+mis_cost


end do Mismatch_error


mis_cost=mis_cost/((window_size-1)*1.0_r8)

end subroutine cal_mismatch_cost_function



!-----------------------------------------------------------


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
