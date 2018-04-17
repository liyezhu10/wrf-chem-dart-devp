! Initial graph coloring version - has lots of issues:
!  - very simplified 'filter' routine
!  - needs to handle ob 'skips'
!  - communication is element-wise in chunks, needs optimization
!  - needs a general cleanup... focus right now is just to get it building

!>  A coloring-enabled version of the assim_tools_mod functionality
module assim_graph_tools_mod

use      types_mod,       only : r8, i8, digits12, PI, missing_r8
use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output,    &
                                 find_namelist_in_file, register_module, error_handler,   &
                                 E_ERR, E_MSG, nmlfileunit, do_nml_file, do_nml_term,     &
                                 open_file, close_file, timestamp
use       sort_mod,       only : index_sort 
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq,       &
                                 random_uniform

use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
                                 init_obs, get_obs_from_key, get_obs_def, get_obs_values, &
                                 destroy_obs
   
use          obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time,    &
                                 get_obs_def_error_variance, get_obs_def_type_of_obs

use         obs_kind_mod, only : get_num_types_of_obs, get_index_for_type_of_obs,                   &
                                 get_quantity_for_type_of_obs, assimilate_this_type_of_obs

use       cov_cutoff_mod, only : comp_cov_factor

use       reg_factor_mod, only : comp_reg_factor

use       obs_impact_mod, only : allocate_impact_table, read_impact_table, free_impact_table

use sampling_error_correction_mod, only : get_sampling_error_table_size, &
                                          read_sampling_error_correction

use         location_mod, only : location_type, get_close_type, query_location,           &
                                 operator(==), set_location_missing, write_location,      &
                                 LocationDims, is_vertical, vertical_localization_on,     &
                                 set_vertical, has_vertical_choice, get_close_init,       &
                                 get_vertical_localization_coord, get_close_destroy,      &
                                 set_vertical_localization_coord

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars,             & 
                                 compute_copy_mean_var, get_var_owner_index,              &
                                 prepare_to_update_copies, map_pe_to_task

use mpi_utilities_mod,    only : my_task_id, broadcast_send, broadcast_recv,              & 
                                 sum_across_tasks, task_count, start_mpi_timer,           &
                                 read_mpi_timer, task_sync

use adaptive_inflate_mod, only : do_obs_inflate,  do_single_ss_inflate,                   &
                                 do_varying_ss_inflate,                                   &
                                 update_inflation,                                        &
                                 inflate_ens, adaptive_inflate_type,                      &
                                 deterministic_inflate, solve_quadratic

use time_manager_mod,     only : time_type, get_time

use assim_model_mod,      only : get_state_meta_data,                                     &
                                 get_close_obs,         get_close_state,                  &
                                 convert_vertical_obs,  convert_vertical_state

use distributed_state_mod, only : create_mean_window, free_mean_window

use quality_control_mod, only : good_dart_qc, DARTQC_FAILED_VERT_CONVERT

use assim_tools_mod, only: assim_tools_init, obs_increment, get_my_obs_loc, &
                           log_namelist_selections, update_from_obs_inc

use perf_mod, only: t_startf, t_stopf

implicit none
private

public :: filter_assim_chunks 


! Indicates if module initialization subroutine has been called yet
logical :: module_initialized = .false.
!
integer :: print_timestamps    = 0
integer :: print_trace_details = 0
!
!! True if random sequence needs to be initialized
logical                :: first_inc_ran_call = .true.
type (random_seq_type) :: inc_ran_seq
!
integer                :: num_types = 0
real(r8), allocatable  :: cutoff_list(:)
logical                :: has_special_cutoffs
logical                :: close_obs_caching = .true.
real(r8), parameter    :: small = epsilon(1.0_r8)   ! threshold for avoiding NaNs/Inf
!
!! true if we have multiple vert choices and we're doing vertical localization
!! (make it a local variable so we don't keep making subroutine calls)
logical                :: is_doing_vertical_conversion = .false.
!
character(len = 255)   :: msgstring, msgstring2, msgstring3

! Need to read in table for off-line based sampling correction and store it
integer                :: sec_table_size
real(r8), allocatable  :: exp_true_correl(:), alpha(:)

! if adjust_obs_impact is true, read in triplets from the ascii file
! and fill this 2d impact table. 
real(r8), allocatable  :: obs_impact_table(:,:)

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/assimilation_code/modules/assimilation/assim_tools_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 11799 $"
character(len=128), parameter :: revdate  = "$Date: 2017-07-07 15:08:09 -0600 (Fri, 07 Jul 2017) $"

!

! Graph Coloring type:
type colors_type

   integer :: chunk_size  ! Read from namelist (set to default for now)
   integer :: num_colors  ! Calculated from data

   integer, dimension(:), allocatable :: obs_color ! Color of 1-n observations (read from data)

!   integer, dimension(:), allocatable :: owner       ! Rank that owns this observation (computed from data)
end type colors_type

integer, parameter :: max_chunk_size = 8! for now
type chunk_type
    integer :: num_obs
    integer :: owner

    integer, dimension(max_chunk_size) :: obs_list
end type chunk_type

type chunk_data_type
   integer :: num_obs
   real(r8), dimension(:,:), allocatable :: obs_prior ! (ob in chunk, ens size)
   real(r8), dimension(:,:), allocatable :: obs_inc
   real(r8), dimension(:,:), allocatable :: net_a     ! (ob in chunk, num_groups)
   real(r8), dimension(:), allocatable   :: obs_qc
   real(r8), dimension(:), allocatable   :: vertvalue_obs_in_localization_coord
   real(r8), dimension(:), allocatable   :: whichvert_real

   !type(location_type), dimension(:), allocatable  :: base_obs_loc

   real(r8), dimension(:), allocatable :: bcast_buffer
end type chunk_data_type

!============================================================================

!---- namelist with default values

! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
!      5 = random draw from posterior
!      6 = deterministic draw from posterior with fixed kurtosis
!      8 = Rank Histogram Filter (see Anderson 2011)
!
!  special_localization_obs_types -> Special treatment for the specified observation types
!  special_localization_cutoffs   -> Different cutoff value for each specified obs type
!!
integer  :: filter_kind                     = 1
real(r8) :: cutoff                          = 0.2_r8
logical  :: sort_obs_inc                    = .false.
logical  :: spread_restoration              = .false.
logical  :: sampling_error_correction       = .false.
integer  :: adaptive_localization_threshold = -1
real(r8) :: adaptive_cutoff_floor           = 0.0_r8
integer  :: print_every_nth_obs             = 0
!
!! since this is in the namelist, it has to have a fixed size.
integer, parameter   :: MAX_ITEMS = 300
character(len = 129) :: special_localization_obs_types(MAX_ITEMS)
real(r8)             :: special_localization_cutoffs(MAX_ITEMS)

logical              :: output_localization_diagnostics = .false.
character(len = 129) :: localization_diagnostics_file = "localization_diagnostics"
!
!! Following only relevant for filter_kind = 8
logical  :: rectangular_quadrature          = .true.
logical  :: gaussian_likelihood_tails       = .false.
!
!! Some models are allowed to have MISSING_R8 values in the DART state vector.
!! If they are encountered, it is not necessarily a FATAL error.
!! Most of the time, if a MISSING_R8 is encountered, DART should die.
!! CLM should have allow_missing_in_clm = .true.
!! maybe POP - but in POP the missing values are land and all ensemble members
!! have the same missing values.  CLM is different in that only some ensemble members may
!! have missing values and so we have a deficient ensemble size at those state locations.
logical  :: allow_missing_in_clm = .false.
!
!! False by default; if true, expect to read in an ascii table
!! to adjust the impact of obs on other state vector and obs values.
logical            :: adjust_obs_impact  = .false.
character(len=256) :: obs_impact_filename = ''
logical            :: allow_any_impact_values = .false.
!
!! These next two only affect models with multiple options
!! for vertical localization:
!!
!! "convert_state" is false by default; it depends on the model
!! what is faster - do the entire state up front and possibly
!! do unneeded work, or do the conversion during the assimilation
!! loop. we think this depends heavily on how much of the state
!! is going to be adjusted by the obs.  for a global model
!! we think false may be better; for a regional model with
!! a lot of obs and full coverage true may be better.
!!
!! "convert_obs" is true by default; in general it seems to
!! be better for each task to convert the obs vertical before
!! going into the loop but again this depends on how many
!! obs per task and whether the mean is distributed or 
!! replicated on each task.
logical :: convert_all_state_verticals_first = .false.
logical :: convert_all_obs_verticals_first   = .true.
!
!! Not in the namelist; this var disables the experimental
!! linear and spherical case code in the adaptive localization 
!! sections.  to try out the alternatives, set this to .false.
logical  :: only_area_adapt  = .true.
!
!! Option to distribute the mean.  If 'false' each task will have a full
!! copy of the ensemble mean, which speeds models doing vertical conversion.
!! If 'true' the mean will be spread across all tasks which reduces the
!! memory needed per task but requires communication if the mean is used
!! for vertical conversion.  We have changed the default to be .false.
!! compared to previous versions of this namelist item.
logical  :: distribute_mean  = .false.
!!logical  :: distribute_mean  = .true.  ! this causes hangs?  weird (BPD6)
!
!! Lanai bitwise. This is for unit testing and runs much slower.
!! Only use for when testing against the non-rma trunk.
logical  :: lanai_bitwise = .false.
!
namelist / assim_tools_nml / filter_kind, cutoff, sort_obs_inc, &
   spread_restoration, sampling_error_correction,                          & 
   adaptive_localization_threshold, adaptive_cutoff_floor,                 &
   print_every_nth_obs, rectangular_quadrature, gaussian_likelihood_tails, &
   output_localization_diagnostics, localization_diagnostics_file,         &
   special_localization_obs_types, special_localization_cutoffs,           &
   allow_missing_in_clm, distribute_mean, close_obs_caching,               &
   adjust_obs_impact, obs_impact_filename, allow_any_impact_values,        &
   convert_all_state_verticals_first, convert_all_obs_verticals_first,     &
   lanai_bitwise ! don't document this one -- only used for regression tests

!============================================================================

contains


subroutine filter_assim_chunks(ens_handle, obs_ens_handle, obs_seq, keys,           &
   ens_size, num_groups, obs_val_index, inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,          &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,            &
   OBS_PRIOR_VAR_END, inflate_only)

type(ensemble_type),         intent(inout) :: ens_handle, obs_ens_handle
type(obs_sequence_type),     intent(in)    :: obs_seq
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: ens_size, num_groups, obs_val_index
type(adaptive_inflate_type), intent(inout) :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, ENS_INF_COPY
integer,                     intent(in)    :: ENS_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END
logical,                     intent(in)    :: inflate_only

!>@todo FIXME this routine has a huge amount of local/stack storage.
!>at some point does it need to be allocated instead?  this routine isn't
!>called frequently so doing allocate/deallocate isn't a timing issue.  
!>putting arrays on the stack is fast, but risks running out of stack space 
!>and dying with strange errors.

real(r8) :: obs_prior(ens_size), obs_inc(ens_size), increment(ens_size)
real(r8) :: reg_factor, impact_factor
real(r8) :: net_a(num_groups), reg_coef(num_groups), correl(num_groups)
real(r8) :: cov_factor, obs(1), obs_err_var, my_inflate, my_inflate_sd
real(r8) :: varying_ss_inflate, varying_ss_inflate_sd
real(r8) :: ss_inflate_base, obs_qc, cutoff_rev, cutoff_orig
real(r8) :: gamma, ens_obs_mean, ens_obs_var, ens_var_deflate
real(r8) :: r_mean, r_var
real(r8) :: orig_obs_prior_mean(num_groups), orig_obs_prior_var(num_groups)
real(r8) :: obs_prior_mean(num_groups), obs_prior_var(num_groups)
real(r8) :: close_obs_dist(obs_ens_handle%my_num_vars)
real(r8) :: close_state_dist(ens_handle%my_num_vars)
real(r8) :: last_close_obs_dist(obs_ens_handle%my_num_vars)
real(r8) :: last_close_state_dist(ens_handle%my_num_vars)
real(r8) :: diff_sd, outlier_ratio

integer(i8) :: state_index
integer(i8) :: my_state_indx(ens_handle%my_num_vars)
integer(i8) :: my_obs_indx(obs_ens_handle%my_num_vars)

integer  :: my_num_obs, i, j, owner, owners_index, my_num_state
integer  :: this_obs_key, obs_mean_index, obs_var_index
integer  :: grp_beg(num_groups), grp_end(num_groups), grp_size, grp_bot, grp_top, group
integer  :: close_obs_ind(obs_ens_handle%my_num_vars)
integer  :: close_state_ind(ens_handle%my_num_vars)
integer  :: last_close_obs_ind(obs_ens_handle%my_num_vars)
integer  :: last_close_state_ind(ens_handle%my_num_vars)
integer  :: num_close_obs, obs_index, num_close_states
integer  :: total_num_close_obs, last_num_close_obs, last_num_close_states
integer  :: base_obs_kind, base_obs_type, my_obs_kind(obs_ens_handle%my_num_vars)
integer  :: my_obs_type(obs_ens_handle%my_num_vars)
integer  :: my_state_kind(ens_handle%my_num_vars), nth_obs
integer  :: num_close_obs_cached, num_close_states_cached
integer  :: num_close_obs_calls_made, num_close_states_calls_made
! GSR add new count for only the 'assimilate' type close obs in the tile
integer  :: localization_unit, secs, days, rev_num_close_obs
character(len = 102)  :: base_loc_text   ! longest location formatting possible

type(location_type)  :: my_obs_loc(obs_ens_handle%my_num_vars)
type(location_type)  :: base_obs_loc, last_base_obs_loc, last_base_states_loc
type(location_type)  :: my_state_loc(ens_handle%my_num_vars), dummyloc
type(get_close_type) :: gc_obs, gc_state
type(obs_type)       :: observation

integer :: last_rank

type(obs_def_type)   :: obs_def
type(time_type)      :: obs_time, this_obs_time

logical :: do_adapt_inf_update
logical :: missing_in_state
! for performance, local copies 
logical :: local_single_ss_inflate
logical :: local_varying_ss_inflate
logical :: local_obs_inflate

! HK observation location conversion
real(r8) :: vertvalue_obs_in_localization_coord
integer  :: whichvert_obs_in_localization_coord
real(r8) :: whichvert_real
type(location_type) :: lc(1)
integer             :: kd(1)

! timing - set one or both of the parameters to true
! to get timing info printed out.
real(digits12) :: base, elapsed, base2
logical, parameter :: timing = .false.
logical, parameter :: timing1 = .false.
real(digits12), allocatable :: elapse_array(:)

integer :: istatus 
integer :: vstatus(obs_ens_handle%my_num_vars) !< for vertical conversion status.

! bpd6
type(colors_type) :: colors
integer :: histogram_unit ! bpd6
integer :: list_unit ! bpd6
integer :: obdata_unit ! bpd6
integer :: obdata_unit2 ! bpd6
character(len = 129) :: histogram_file = "histogram_data.txt"
character(len = 129) :: list_file = "list_data.txt"
character(len = 129) :: colors_file = "colors.txt"

! new mods, 2018-01-16:
character(len = 129) :: obdata_file = "obdata.txt"
character(len = 129) :: obdata_file2 = "obdata2.txt"

integer  :: total_close_ranks, rev_close_ranks

integer :: qcd = 0 ! bpd6
integer :: iError
real(r8) :: testval, testval2
integer :: k

integer :: skipped_missing = 0
integer :: skipped_covfactor = 0
real(r8) :: stateupdate_time = 0.0d0
INTEGER(kind=8) :: timer_count, timer_rate, timer_max
INTEGER(kind=8) :: timer_count2, timer_rate2, timer_max2


! coloring - bpd6
!logical :: own_color
!integer :: obs_set_size
!integer, dimension(8000) :: obs_set
!integer, dimension(:), allocatable :: obs_list
type(chunk_type), dimension(:), allocatable :: chunks


type(chunk_data_type) chunk_data

integer (kind=8) :: ob_index

!call task_sync()
call t_startf('ASSIMILATE:Pre.Loop')

! we are going to read/write the copies array
call prepare_to_update_copies(ens_handle)
call prepare_to_update_copies(obs_ens_handle)

! Initialize assim_tools_module if needed
!if (.not. module_initialized) call assim_tools_init()

!HK make window for mpi one-sided communication
! used for vertical conversion in get_close_obs
! Need to give create_mean_window the mean copy
call create_mean_window(ens_handle, ENS_MEAN_COPY, distribute_mean)

! filter kinds 1 and 8 return sorted increments, however non-deterministic
! inflation can scramble these. the sort is expensive, so help users get better 
! performance by rejecting namelist combinations that do unneeded work.
if (sort_obs_inc) then
   if(deterministic_inflate(inflate) .and. ((filter_kind == 1) .or. (filter_kind == 8))) then
      write(msgstring,  *) 'With a deterministic filter [assim_tools_nml:filter_kind = ',filter_kind,']'
      write(msgstring2, *) 'and deterministic inflation [filter_nml:inf_deterministic = .TRUE.]'
      write(msgstring3, *) 'assim_tools_nml:sort_obs_inc = .TRUE. is not needed and is expensive.'
      call error_handler(E_MSG,'', '')  ! whitespace
      call error_handler(E_MSG,'WARNING filter_assim:', msgstring, source, revision, revdate, &
                         text2=msgstring2,text3=msgstring3)
      call error_handler(E_MSG,'', '')  ! whitespace
      sort_obs_inc = .FALSE.
   endif
endif

if (my_task_id() == 0) then
  histogram_unit = open_file(histogram_file) ! bpd6
  list_unit = open_file(list_file) ! bpd6
  obdata_unit = open_file(obdata_file) ! bpd6 ! New, 2018-01-16
  obdata_unit2 = open_file(obdata_file2) ! bpd6 ! New, 2018-01-16
endif

!GSR open the dignostics file
if(output_localization_diagnostics .and. my_task_id() == 0) then
  localization_unit = open_file(localization_diagnostics_file, action = 'append')
endif

! For performance, make local copies of these settings which
! are really in the inflate derived type.
local_single_ss_inflate  = do_single_ss_inflate(inflate)
local_varying_ss_inflate = do_varying_ss_inflate(inflate)
local_obs_inflate        = do_obs_inflate(inflate)

! Default to printing nothing
nth_obs = -1

! Divide ensemble into num_groups groups.
! make sure the number of groups and ensemble size result in 
! at least 2 members in each group (to avoid divide by 0) and 
! that the groups all have the same number of members.
grp_size = ens_size / num_groups
if ((grp_size * num_groups) /= ens_size) then
   write(msgstring,  *) 'The number of ensemble members must divide into the number of groups evenly.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'Change number of groups or ensemble size to avoid remainders.'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, revision, revdate, &
                         text2=msgstring2,text3=msgstring3)
endif
if (grp_size < 2) then
   write(msgstring,  *) 'There must be at least 2 ensemble members in each group.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'results in < 2 members/group.  Decrease number of groups or increase ensemble size'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, revision, revdate, &
                         text2=msgstring2,text3=msgstring3)
endif
do group = 1, num_groups
   grp_beg(group) = (group - 1) * grp_size + 1
   grp_end(group) = grp_beg(group) + grp_size - 1
enddo

! Put initial value of state space inflation in copy normally used for SD
! This is to avoid weird storage footprint in filter
ens_handle%copies(ENS_SD_COPY, :) = ens_handle%copies(ENS_INF_COPY, :)

! For single state or obs space inflation, the inflation is like a token
! Gets passed from the processor with a given obs on to the next
if(local_single_ss_inflate) then
   my_inflate    = ens_handle%copies(ENS_INF_COPY,    1)
   my_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, 1)
end if


! Get info on my number and indices for obs
my_num_obs = get_my_num_vars(obs_ens_handle)
call get_my_vars(obs_ens_handle, my_obs_indx)

! Construct an observation temporary
call init_obs(observation, get_num_copies(obs_seq), get_num_qc(obs_seq))

! Get the locations for all of my observations 
! HK I would like to move this to before the calculation of the forward operator so you could
! overwrite the vertical location with the required localization vertical coordinate when you 
! do the forward operator calculation
call get_my_obs_loc(ens_handle, obs_ens_handle, obs_seq, keys, my_obs_loc, my_obs_kind, my_obs_type, obs_time)

if (convert_all_obs_verticals_first .and. is_doing_vertical_conversion) then
   ! convert the vertical of all my observations to the localization coordinate
   ! this may not be bitwise with Lanai because of a different number of set_location calls
   if (timing) call start_mpi_timer(base)
   call convert_vertical_obs(ens_handle, obs_ens_handle%my_num_vars, my_obs_loc, &
                             my_obs_kind, my_obs_type, get_vertical_localization_coord(), vstatus)
   do i = 1, obs_ens_handle%my_num_vars
      if (good_dart_qc(nint(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)))) then
         !> @todo Can I just use the OBS_GLOBAL_QC_COPY? Is it ok to skip the loop?
         if (vstatus(i) /= 0) obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i) = DARTQC_FAILED_VERT_CONVERT
      endif
   enddo
   if (timing) then
      elapsed = read_mpi_timer(base)
      print*, 'convert_vertical_obs time :', elapsed, 'rank ', my_task_id()
   endif
endif

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state_indx)

! Get the location and kind of all my state variables
if (timing) call start_mpi_timer(base)
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state_indx(i), my_state_loc(i), my_state_kind(i))
end do
if (timing) then
   elapsed = read_mpi_timer(base)
   print*, 'get_state_meta_data time :', elapsed, 'rank ', my_task_id()
endif

!call test_get_state_meta_data(my_state_loc, ens_handle%my_num_vars)

!> optionally convert all state location verticals
if (convert_all_state_verticals_first .and. is_doing_vertical_conversion) then
   if (timing) call start_mpi_timer(base)
   call convert_vertical_state(ens_handle, ens_handle%my_num_vars, my_state_loc, my_state_kind,  &
                                            my_state_indx, get_vertical_localization_coord(), istatus)
   if (timing) then
      elapsed = read_mpi_timer(base)
      print*, 'convert_vertical_state time :', elapsed, 'rank ', my_task_id()
   endif
endif

! PAR: MIGHT BE BETTER TO HAVE ONE PE DEDICATED TO COMPUTING 
! INCREMENTS. OWNING PE WOULD SHIP IT'S PRIOR TO THIS ONE
! BEFORE EACH INCREMENT.

! Get mean and variance of each group's observation priors for adaptive inflation
! Important that these be from before any observations have been used
if(local_varying_ss_inflate .or. local_single_ss_inflate) then
   do group = 1, num_groups
      obs_mean_index = OBS_PRIOR_MEAN_START + group - 1
      obs_var_index  = OBS_PRIOR_VAR_START  + group - 1
         call compute_copy_mean_var(obs_ens_handle, grp_beg(group), grp_end(group), &
           obs_mean_index, obs_var_index) 
   end do
endif

! The computations in the two get_close_maxdist_init are redundant

! Initialize the method for getting state variables close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc)
endif

! Initialize the method for getting obs close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc)
endif

if (close_obs_caching) then
   ! Initialize last obs and state get_close lookups, to take advantage below 
   ! of sequential observations at the same location (e.g. U,V, possibly T,Q)
   ! (this is getting long enough it probably should go into a subroutine. nsc.)
   last_base_obs_loc           = set_location_missing()
   last_base_states_loc        = set_location_missing()
   last_num_close_obs          = -1
   last_num_close_states       = -1
   last_close_obs_ind(:)       = -1
   last_close_state_ind(:)     = -1
   last_close_obs_dist(:)      = 888888.0_r8   ! something big, not small
   last_close_state_dist(:)    = 888888.0_r8   ! ditto
   num_close_obs_cached        = 0
   num_close_states_cached     = 0
   num_close_obs_calls_made    = 0
   num_close_states_calls_made = 0
endif



!bpd6 - get the coloring info
call read_obs_colors(colors_file, obs_ens_handle%num_vars, colors)
call create_chunks(colors, chunks)
call initialize_chunk_data(colors%chunk_size, ens_size, num_groups, chunk_data)

! timing
if (my_task_id() == 0 .and. timing) allocate(elapse_array(obs_ens_handle%num_vars))

call t_stopf('ASSIMILATE:Pre.Loop')
call task_sync()
call t_startf('tmp_loop')

write(*,*) "DEBUG2: About to enter loop 1 -> ", colors%num_colors
write(*,*) "Ensemble Size : ", ens_size


! Loop through all the chunks:
CHUNK_LOOP: do i = 1, size(chunks)

   !write(*,*) "Chunk loop : ", i, size(chunks) 
   !write(*,*) "Loop start/end -> ", chunks(i)%obs_list(1), chunks(i)%obs_list(chunks(i)%num_obs)

   !write(*,'(A5,I4,A3,I4,A4,I4,A3,I7,A3,I7)'), "Rank",my_task_id(),"C#",i,"Own",chunks(i)%owner,"St",chunks(i)%obs_list(1),"En",chunks(i)%obs_list(chunks(i)%num_obs)

   ! This section is only done by one process - the 'owner' of this chunk:
   if (ens_handle%my_pe == chunks(i)%owner) then
      call t_startf('ASSIMILATE:Owned(Compute)')

      write(*,'(A5,I4,A3,I4,A4,I4,A3,I7,A3,I7)'), "Rank",my_task_id(),"C#",i,"Own",chunks(i)%owner,"St",chunks(i)%obs_list(1),"En",chunks(i)%obs_list(chunks(i)%num_obs)
      chunk_data%num_obs = chunks(i)%num_obs

      OBS_LOOP: do j = 1, chunk_data%num_obs
         !write(*,*) "Ob loop : ", j, chunk_data%num_obs
        ! Get the index of this ob into the main sequence:
        ob_index = chunks(i)%obs_list(j)
        !write(*,*) "ob_index = ", ob_index, i, j



        ! Get owners index: (?) - we ignore owner here, so this can be changed,
        ! but need to understand it better first.
        call get_var_owner_index(int(ob_index,i8), owner, owners_index)
        !write(*,*) "DEBUG: ", int(ob_index,i8), owner, owners_index

        ! Get the QC value for this ob:
        !write(*,*) "Checking QC : ", j
        chunk_data%obs_qc(j) = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)

        ! Only value of 0 for DART QC field should be assimilated
        IF_QC_IS_OKAY: if(nint(obs_qc) ==0) then
           chunk_data%obs_prior(j,:) = obs_ens_handle%copies(1:ens_size, owners_index)
           !write(*,*) "DEBUG : obs_prior_sum = ", sum(chunk_data%obs_prior(j,:)), chunks(i)%obs_list(j)

           ! Compute the prior mean and variance for this observation
           orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: OBS_PRIOR_MEAN_END, owners_index) ! unused for now,
           orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:  OBS_PRIOR_VAR_END, owners_index)  ! unused for now
 
           ! Get the value of the observation
           call get_obs_from_key(obs_seq, keys(ob_index), observation)
           call get_obs_def(observation, obs_def)
           call get_obs_values(observation, obs, obs_val_index)
           obs_err_var = get_obs_def_error_variance(obs_def)  ! Add to the chunk data type?


           ! Compute observation space increments for each group
           do group = 1, num_groups
              grp_bot = grp_beg(group)
              grp_top = grp_end(group)
              !!!call obs_increment(obs_prior(grp_bot:grp_top), grp_size, obs(1), obs_err_var, obs_inc(grp_bot:grp_top), inflate, my_inflate, my_inflate_sd, net_a(group))
              call obs_increment(chunk_data%obs_prior(j,grp_bot:grp_top), grp_size, obs(1), obs_err_var, chunk_data%obs_inc(j,grp_bot:grp_top), inflate, my_inflate, my_inflate_sd, chunk_data%net_a(j,group))
           end do

           ! ------- NOTE: Skipping SINGLE_SS_INFLATE section for now -------
         endif IF_QC_IS_OKAY

         ! ----- NOTE: Skipping vertical conversion section for now ------

     enddo OBS_LOOP
     ! Haven't implemented the two other kinds of broadcasts yet, so no 'if' here:
     !!!call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior, obs_inc, net_a, scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord, scalar3=whichvert_real)

    call t_stopf('ASSIMILATE:Owned(Compute)')
    call task_sync()
    call t_startf('ASSIMILATE:Owned(Broadcast)')

     !write(*,*) "Calling broadcast_send_chunk on owner", chunks(i)%owner
     call broadcast_send_chunk(map_pe_to_task(ens_handle, chunks(i)%owner), chunk_data)
    call task_sync()

   else ! (not the owner):
     call task_sync()
     call t_startf('ASSIMILATE:NotOwned(Broadcast)')
     !write(*,*) "Calling broadcast_recv_chunk on ", ens_handle%my_pe, map_pe_to_task(ens_handle, chunks(i)%owner)
     call broadcast_recv_chunk(map_pe_to_task(ens_handle, chunks(i)%owner), chunk_data)
     call t_stopf('ASSIMILATE:NotOwned(Broadcast)')
     call task_sync()
     call t_startf('ASSIMILATE:NotOwned(Compute)')
     call t_stopf('ASSIMILATE:NotOwned(Compute)')
  endif

  call task_sync()

  QC_CHECK: do j = 1, chunk_data%num_obs
     if (nint(chunk_data%obs_qc(j)) /= 0) then
        ! This is a hack - we're going to swap the current value with the last
        ! value, then decrement the count, effectively removing this one.  If
        ! we're at the last one, we cycle to the next chunk:
        if (j == chunk_data%num_obs) then
           cycle CHUNK_LOOP
        else
           chunk_data%obs_prior(j,:) = chunk_data%obs_prior(chunk_data%num_obs,:)
           chunk_data%obs_inc(j,:) = chunk_data%obs_inc(chunk_data%num_obs,:)
           chunk_data%net_a(j,:) = chunk_data%net_a(chunk_data%num_obs,:)
           chunk_data%obs_qc(j) = chunk_data%obs_qc(chunk_data%num_obs)
           chunk_data%vertvalue_obs_in_localization_coord(j) = chunk_data%vertvalue_obs_in_localization_coord(chunk_data%num_obs)
           chunk_data%whichvert_real(j) = chunk_data%whichvert_real(chunk_data%num_obs)
           chunk_data%num_obs = chunk_data%num_obs - 1
           !j = j -1
           write(*,*) "QC_CHECK fail on ", i, j
           cycle QC_CHECK
       endif
    endif
  enddo QC_CHECK


  do j = 1, chunk_data%num_obs
    call t_startf('ASSIMILATE:ComputePriors')
   ! Can compute prior mean and variance of obs for each group just once here
    do group = 1, num_groups
      grp_bot = grp_beg(group)
      grp_top = grp_end(group)
      !write(*,*) "DEBUG : obs_prior_sum(2) = ", sum(chunk_data%obs_prior(j,:)), chunks(i)%obs_list(j)
      obs_prior_mean(group) = sum(chunk_data%obs_prior(j,grp_bot:grp_top)) / grp_size
      obs_prior_var(group) = sum((chunk_data%obs_prior(j,grp_bot:grp_top) - obs_prior_mean(group))**2) / &
         (grp_size - 1)
      if (obs_prior_var(group) < 0.0_r8) obs_prior_var(group) = 0.0_r8
      !write(*,*) "DEBUG : obs_prior_mean = ", obs_prior_mean(group), chunks(i)%obs_list(j)
    end do

   call t_stopf('ASSIMILATE:ComputePriors')
   !call task_sync()
   ! -------------- NOTE: Skipping all adaptive localization stuff for now ---------
   ! ------- NOTE: Turns out, get_close_state was in the            ----
   ! -------       adaptivelocalization section .. we need to do it ----


   ! Do we really need to do get_close_states in the coloring mode?  We already
   ! know the close states!  Something to think about.  Even if storing all
   ! the info is too memory intensive, maybe storing a subgrid of 8 cells would
   ! help?  Measure this to see how long it takes, then decide.


    call t_startf('ASSIMILATE:AdaptiveLocalization')
   ! Every pe has information about the global obs sequence
   if (chunks(i)%obs_list(j) > 192895) then
       write(*,*) "Trace: ", i, j
    endif
   call get_obs_from_key(obs_seq, keys(chunks(i)%obs_list(j)), observation)
   call get_obs_def(observation, obs_def)
   base_obs_loc = get_obs_def_location(obs_def)
   obs_err_var = get_obs_def_error_variance(obs_def)
   base_obs_type = get_obs_def_type_of_obs(obs_def)
   if (base_obs_type > 0) then
     base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
   else
     call get_state_meta_data(-1*int(base_obs_type,i8),dummyloc, base_obs_kind)  ! identity obs
   endif
   
   if (.not. close_obs_caching) then
      if (timing) call start_mpi_timer(base)
      call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                         my_obs_loc, my_obs_kind, my_obs_type, &
                         num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
      if (timing) then
         elapsed = read_mpi_timer(base)
         print*, 'get_close_obs1 time :', elapsed, 'rank ', my_task_id()
      endif

   else
 
      if (base_obs_loc == last_base_obs_loc) then
         num_close_obs     = last_num_close_obs
         close_obs_ind(:)  = last_close_obs_ind(:)
         close_obs_dist(:) = last_close_obs_dist(:)
         num_close_obs_cached = num_close_obs_cached + 1
      else
         if (timing .and. i < 100) call start_mpi_timer(base)
         call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                            my_obs_loc, my_obs_kind, my_obs_type, &
                            num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
         if (timing .and. i < 100) then
            elapsed = read_mpi_timer(base)
            print*, 'get_close_obs2 time :', elapsed, 'rank ', my_task_id()
         endif

         last_base_obs_loc      = base_obs_loc
         last_num_close_obs     = num_close_obs
         last_close_obs_ind(:)  = close_obs_ind(:)
         last_close_obs_dist(:) = close_obs_dist(:)
         num_close_obs_calls_made = num_close_obs_calls_made +1
      endif
   endif

   if (.not. close_obs_caching) then
      call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                           my_state_loc, my_state_kind, my_state_indx, &
                           num_close_states, close_state_ind, close_state_dist, ens_handle)
      if (timing .and. i < 100) then
         elapsed = read_mpi_timer(base)
         print*, 'get_close_state1 time :', elapsed, 'rank ', my_task_id()
      endif
   else
      if (base_obs_loc == last_base_states_loc) then
         num_close_states    = last_num_close_states
         close_state_ind(:)  = last_close_state_ind(:)
         close_state_dist(:) = last_close_state_dist(:)
         num_close_states_cached = num_close_states_cached + 1
     else
         if (timing .and. i < 100) call start_mpi_timer(base)
         call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                              my_state_loc, my_state_kind, my_state_indx, &
                              num_close_states, close_state_ind, close_state_dist, ens_handle)
         if (timing .and. i < 100) then
            elapsed = read_mpi_timer(base)
            print*, 'get_close_state2 time :', elapsed, 'rank ', my_task_id()
         endif

         last_base_states_loc     = base_obs_loc
         last_num_close_states    = num_close_states
         last_close_state_ind(:)  = close_state_ind(:)
         last_close_state_dist(:) = close_state_dist(:)
         num_close_states_calls_made = num_close_states_calls_made + 1
      endif
   endif

   if (base_obs_type > 0) then
      cutoff_orig = cutoff_list(base_obs_type)
   else
      cutoff_orig = cutoff
   endif
  cutoff_rev = cutoff_orig

   call t_stopf('ASSIMILATE:AdaptiveLocalization')

   call t_startf('ASSIMILATE:UpdateState')

  !write(*,*) "DEBUG : num_close_states = ", num_close_states, chunks(i)%obs_list(j)
  !write(*,*) "close_state_dist(1) = ", close_state_dist(1), chunks(i)%obs_list(j)
  !write(*,*) "cutoff_rev = ", cutoff_rev, chunks(i)%obs_list(j)
  !write(*,*) "my_state_loc(1) = ", my_state_loc(close_state_ind(1))%lon, chunks(i)%obs_list(j)

  !testval = 0.0
  !testval2 = 0.0
   STATE_UPDATE: do k = 1, num_close_states
      state_index = close_state_ind(k)

      if ( allow_missing_in_clm ) then
         ! Some models can take evasive action if one or more of the ensembles have
         ! a missing value. Generally means 'do nothing' (as opposed to DIE)
         missing_in_state = any(ens_handle%copies(1:ens_size, state_index) == MISSING_R8)
         if ( missing_in_state ) then
           skipped_missing = skipped_missing + 1
           cycle STATE_UPDATE
         endif
      endif


      ! Get the initial values of inflation for this variable if state varying inflation
      if(local_varying_ss_inflate) then
         varying_ss_inflate    = ens_handle%copies(ENS_INF_COPY,    state_index)
         varying_ss_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, state_index)
      else
         varying_ss_inflate    = 0.0_r8
         varying_ss_inflate_sd = 0.0_r8
      endif
     
      ! Compute the distance and covariance factor 
      cov_factor = comp_cov_factor(close_state_dist(k), cutoff_rev, &
         base_obs_loc, base_obs_type, my_state_loc(state_index), my_state_kind(state_index))
      !testval = testval + cov_factor !debug
      
      ! if external impact factors supplied, factor them in here
      ! FIXME: this would execute faster for 0.0 impact factors if
      ! we check for that before calling comp_cov_factor.  but it makes
      ! the logic more complicated - this is simpler if we do it after.
      if (adjust_obs_impact) then
         impact_factor = obs_impact_table(base_obs_type, my_state_kind(state_index))
         cov_factor = cov_factor * impact_factor
      endif

      ! If no weight is indicated, no more to do with this state variable
      if(cov_factor <= 0.0_r8) then
          skipped_covfactor = skipped_covfactor + 1
          cycle STATE_UPDATE
      endif

      ! Loop through groups to update the state variable ensemble members
      do group = 1, num_groups
         grp_bot = grp_beg(group)
         grp_top = grp_end(group)
         ! Do update of state, correl only needed for varying ss inflate
!!         if(local_varying_ss_inflate .and. varying_ss_inflate > 0.0_r8 .and. &
!!           varying_ss_inflate_sd > 0.0_r8) then
!!           call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
!!               obs_prior_var(group), obs_inc(grp_bot:grp_top), &
!!               ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
!!               increment(grp_bot:grp_top), reg_coef(group), net_a(group), correl(group))
!!         else
            call update_from_obs_inc(chunk_data%obs_prior(j,grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), chunk_data%obs_inc(j,grp_bot:grp_top), &
               ens_handle%copies(grp_bot:grp_top, state_index), grp_size, &
               increment(grp_bot:grp_top), reg_coef(group), net_a(group))
            !write(*,*) "Inc: ", increment(1)
!!         endif
      end do
      !testval2 = testval2 + increment(1)

      ! Compute an information factor for impact of this observation on this state
      if(num_groups == 1) then
          reg_factor = 1.0_r8
      else
         ! Pass the time along with the index for possible diagnostic output
         ! Compute regression factor for this obs-state pair
         reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, my_state_indx(state_index))
      endif

      ! The final factor is the minimum of group regression factor and localization cov_factor
      reg_factor = min(reg_factor, cov_factor)

!PAR NEED TO TURN STUFF OFF MORE EFFICEINTLY
      ! If doing full assimilation, update the state variable ensemble with weighted increments
      if(.not. inflate_only) then
         ens_handle%copies(1:ens_size, state_index) = &
            ens_handle%copies(1:ens_size, state_index) + reg_factor * increment
      endif

      ! Compute spatially-varying state space inflation
      if(local_varying_ss_inflate) then
!!         ! base is the initial inflate value for this state variable
!!         ss_inflate_base = ens_handle%copies(ENS_SD_COPY, state_index)
!!         ! Loop through each group to update inflation estimate
!!         GroupInflate: do group = 1, num_groups
!!            if(varying_ss_inflate > 0.0_r8 .and. varying_ss_inflate_sd > 0.0_r8) then
!!               ! Gamma is less than 1 for varying ss, see adaptive inflate module
!!               gamma = reg_factor * abs(correl(group))
!!               ! Deflate the inflated variance using the INITIAL state inflate
!!               ! value (before these obs started gumming it up).
!!               ens_obs_mean = orig_obs_prior_mean(group)
!!               ens_obs_var =  orig_obs_prior_var(group)
!!
!!               ! Remove the impact of inflation to allow efficient single pass with assim.
!!               if ( abs(gamma) > small ) then
!!                  ens_var_deflate = ens_obs_var / &
!!                     (1.0_r8 + gamma*(sqrt(ss_inflate_base) - 1.0_r8))**2
!!               else
!!                  ens_var_deflate = ens_obs_var
!!               endif
!!                  
!!               ! If this is inflate only (i.e. posterior) remove impact of this obs.
!!               if(inflate_only .and. &
!!                     ens_var_deflate               > small .and. &
!!                     obs_err_var                   > small .and. & 
!!                     obs_err_var - ens_var_deflate > small ) then 
!!                  r_var  = 1.0_r8 / (1.0_r8 / ens_var_deflate - 1.0_r8 / obs_err_var)
!!                  r_mean = r_var *(ens_obs_mean / ens_var_deflate - obs(1) / obs_err_var)
!!               else
!!                  r_var = ens_var_deflate
!!                  r_mean = ens_obs_mean
!!               endif
!!
!!               ! IS A TABLE LOOKUP POSSIBLE TO ACCELERATE THIS?
!!               ! Update the inflation values
!!               call update_inflation(inflate, varying_ss_inflate, varying_ss_inflate_sd, &
!!                  r_mean, r_var, obs(1), obs_err_var, gamma)
!!            else
!!               ! if we don't go into the previous if block, make sure these
!!               ! have good values going out for the block below
!!               r_mean = orig_obs_prior_mean(group)
!!               r_var =  orig_obs_prior_var(group)
!!            endif
!!
!!            ! Update adaptive values if posterior outlier_ratio test doesn't fail.
!!            ! Match code in obs_space_diags() in filter.f90
!!            do_adapt_inf_update = .true.
!!            if (inflate_only) then
!!               diff_sd = sqrt(obs_err_var + r_var) 
!!               if (diff_sd > 0.0_r8) then
!!                  outlier_ratio = abs(obs(1) - r_mean) / diff_sd
!!                  do_adapt_inf_update = (outlier_ratio <= 3.0_r8) 
!!               endif
!!            endif
!!            if (do_adapt_inf_update) then   
!!               ens_handle%copies(ENS_INF_COPY, state_index) = varying_ss_inflate
!!               ens_handle%copies(ENS_INF_SD_COPY, state_index) = varying_ss_inflate_sd
!!            endif
!!         end do GroupInflate
      endif

   end do STATE_UPDATE
   !write(*,*) "DEBUG : cov_factor sum = ", testval, chunks(i)%obs_list(j)
   !write(*,*) "DEBUG : sum_increment = ", testval2, chunks(i)%obs_list(j)
!!   if (timing .and. i < 1000) then
!!      elapsed = read_mpi_timer(base)
!!      print*, 'state_update time :', elapsed, 'rank ', my_task_id()
!!   endif
!!
!! ! stop the timer for this section for this rank:
!!  call system_clock(count=timer_count2)
!!  !write(*,*) "Debug: ", timer_count, timer_count2, timer_count2-timer_count
!!  stateupdate_time = DBLE(timer_count2 - timer_count) / DBLE(timer_rate)
!!!  write(*,*) "Debug: ", stateupdate_time
!!
!!   !bpd6 - 2018-01-16 mod
!!   !call write_obdata(obdata_unit, i, num_close_states, skipped_missing, skipped_covfactor, stateupdate_time)
!!
   call t_stopf('ASSIMILATE:UpdateState')
   !call task_sync()

   !write(*,*) "DEBUG: num_close_obs = ", num_close_obs


!!   call t_startf('ASSIMILATE:UpdateObs')
!!   !call test_state_copies(ens_handle, 'after_state_updates')
!!
!!   !------------------------------------------------------
!!
!!    !bpd6- new obdata update
!!    call write_obdata2(obdata_unit2, i)
!!
!!   ! Now everybody updates their obs priors (only ones after this one)
!!   if (timing .and. i < 1000) call start_mpi_timer(base)
   OBS_UPDATE: do k = 1, num_close_obs
      obs_index = close_obs_ind(k)

      ! Only have to update obs that have not yet been used
      if(my_obs_indx(obs_index) > chunks(i)%obs_list(j)) then

         ! If the forward observation operator failed, no need to 
         ! update the unassimilated observations 
         if (any(obs_ens_handle%copies(1:ens_size, obs_index) == MISSING_R8)) cycle OBS_UPDATE

         ! Compute the distance and the covar_factor
         cov_factor = comp_cov_factor(close_obs_dist(j), cutoff_rev, base_obs_loc, base_obs_type, my_obs_loc(obs_index), my_obs_kind(obs_index))

         ! if external impact factors supplied, factor them in here
         ! FIXME: this would execute faster for 0.0 impact factors if
         ! we check for that before calling comp_cov_factor.  but it makes
         ! the logic more complicated - this is simpler if we do it after.
         if (adjust_obs_impact) then
            impact_factor = obs_impact_table(base_obs_type, my_obs_kind(obs_index))
            cov_factor = cov_factor * impact_factor
         endif

         if(cov_factor <= 0.0_r8) cycle OBS_UPDATE

        ! bpd6 - NOW, finally, we have an ob we're actually updating.. append to
        ! our obdata file:
        !write(*,*) "BRIAN: ", my_obs_indx(obs_index), obs_index, i
!        call append_obdata2(obdata_unit2, obs_index)


         ! Loop through and update ensemble members in each group
         do group = 1, num_groups
            grp_bot = grp_beg(group)
            grp_top = grp_end(group)
            call update_from_obs_inc(chunk_data%obs_prior(j,grp_bot:grp_top), obs_prior_mean(group), &
               obs_prior_var(group), chunk_data%obs_inc(j,grp_bot:grp_top), &
                obs_ens_handle%copies(grp_bot:grp_top, obs_index), grp_size, &
                increment(grp_bot:grp_top), reg_coef(group), net_a(group))
         end do

         ! FIXME: could we move the if test for inflate only to here?

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 1) then
             reg_factor = 1.0_r8
         else
            ! Pass the time along with the index for possible diagnostic output
            ! Compute regression factor for this obs-state pair
            ! Negative indicates that this is an observation index
            reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, i, -1*my_obs_indx(obs_index))
         endif

         ! Final weight is min of group and localization factors
         reg_factor = min(reg_factor, cov_factor)

         ! Only update state if indicated (otherwise just getting inflation)
         if(.not. inflate_only) then
            obs_ens_handle%copies(1:ens_size, obs_index) = &
              obs_ens_handle%copies(1:ens_size, obs_index) + reg_factor * increment
         endif
      endif
   end do OBS_UPDATE
!!   if (timing .and. i < 1000) then
!!      elapsed = read_mpi_timer(base)
!!      print*, 'obs_update time :', elapsed, 'rank ', my_task_id()
!!   endif
!!
!!   !call test_state_copies(ens_handle, 'after_obs_updates')
!!
!!
!!   if (my_task_id() == 0 .and. timing) then
!!      elapse_array(i) = read_mpi_timer(base2)
!!      if (timing1) print*, 'outer sequential obs time :', elapsed, ' obs ', i, ' rank ', my_task_id()
!!   endif
!!
!!   call t_stopf('ASSIMILATE:UpdateObs')
!!   call task_sync()
!!
  enddo ! Obs loop - shouldn't be here, but this is for testing
end do CHUNK_LOOP


! Temporary
!call MPI_Finalize(iError)
!write(*,*) "Exiting.."
!stop

write(*,*) "QC'd obs: ", qcd

call t_stopf('tmp_loop')
call task_sync()
call t_startf('ASSIMILATE:Post.Loop')


! Every pe needs to get the current my_inflate and my_inflate_sd back
if(local_single_ss_inflate) then
   ens_handle%copies(ENS_INF_COPY, :) = my_inflate
   ens_handle%copies(ENS_INF_SD_COPY, :) = my_inflate_sd
end if

! Free up the storage
call destroy_obs(observation)
call get_close_destroy(gc_state)
call get_close_destroy(gc_obs)

! print some stats about the assimilation
if (my_task_id() == 0 .and. timing) then
   write(msgstring, *) 'average assim time: ', sum(elapse_array) / size(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)

   write(msgstring, *) 'minimum assim time: ', minval(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)

   write(msgstring, *) 'maximum assim time: ', maxval(elapse_array)
   call error_handler(E_MSG,'filter_assim:',msgstring)
endif

if (my_task_id() == 0 .and. timing) deallocate(elapse_array)

! Assure user we have done something
write(msgstring, '(A,I8,A)') &
   'Processed', obs_ens_handle%num_vars, ' total observations'
if (print_trace_details >= 0) call error_handler(E_MSG,'filter_assim:',msgstring)

! diagnostics for stats on saving calls by remembering obs at the same location.
! change .true. to .false. in the line below to remove the output completely.
if (close_obs_caching) then
   if (num_close_obs_cached > 0 .and. do_output()) then
      print *, "Total number of calls made    to get_close_obs for obs/states:    ", &
                num_close_obs_calls_made + num_close_states_calls_made
      print *, "Total number of calls avoided to get_close_obs for obs/states:    ", &
                num_close_obs_cached + num_close_states_cached
      if (num_close_obs_cached+num_close_obs_calls_made+ &
          num_close_states_cached+num_close_states_calls_made > 0) then 
         print *, "Percent saved: ", 100.0_r8 * &
                   (real(num_close_obs_cached+num_close_states_cached, r8) /  &
                   (num_close_obs_calls_made+num_close_obs_cached +           &
                    num_close_states_calls_made+num_close_states_cached))
      endif
   endif
endif

!call test_state_copies(ens_handle, 'end')

!GSR close the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) then
  call close_file(localization_unit)
end if

if (my_task_id() == 0) then
  call close_file(histogram_unit) ! bpd6
endif

! get rid of mpi window
call free_mean_window()

call task_sync()
call t_stopf('ASSIMILATE:Post.Loop')

end subroutine filter_assim_chunks

!-------------------------------------------------------------

subroutine read_obs_colors(filename, num_obs, colors)
   character(len=*),  intent(in)  :: filename
   integer(kind=8),           intent(in)  :: num_obs
   type(colors_type), intent(out) :: colors

  integer :: colorfile
  integer :: i

  ! Allocate the array
  allocate(colors%obs_color(num_obs))

  ! open the file 
  colorfile = open_file(filename)

  ! Read in the colors, looping over all the obs:
  do i = 1, num_obs
    read(colorfile,*) colors%obs_color(i)
  enddo

  ! Calculate the number of colors:
  colors%num_colors = maxval(colors%obs_color)

  ! Set the chunk size to be 8 for now - will change this later:
  colors%chunk_size = 3
end subroutine read_obs_colors

!-------------------------------------------------------------

subroutine create_chunks(colors, chunks)
  type(colors_type), intent(in) :: colors
  type(chunk_type), dimension(:), allocatable, intent(out) :: chunks

  integer :: chunk_count
  integer :: i,j
  integer, dimension(:), allocatable :: color_sizes, chunks_per_color
  integer :: chunk_index = 1
  integer :: ob_index = 0
  integer :: remaining_obs, numRanks

  ! Allocate color_sizes (array of sizes for each color)
  allocate(color_sizes(colors%num_colors))
  allocate(chunks_per_color(colors%num_colors))

  write(*,*) "Colors -- number of colors & chunk size = ", colors%num_colors, colors%chunk_size

  ! Get the size and # of chunks of each color - this is a bit hackish, maybe ANY can work better?
  do i = 1, colors%num_colors
    color_sizes(i) = sum(colors%obs_color, colors%obs_color==i) / i
    chunks_per_color(i) = (color_sizes(i) + colors%chunk_size - 1) / colors%chunk_size

    !write(*,*) "Color sizes (",i,") ",color_sizes(i)
    !write(*,*) "Color chunks (",i,") ",chunks_per_color(i)
  enddo 
  
  ! Get the total number of chunks:
  chunk_count = sum(chunks_per_color)

  ! allocate the chunks array:
  allocate(chunks(chunk_count))

  ! Assign chunks
  write(*,*) "Assigning chunks..."
  numRanks = task_count()
  do i = 1, colors%num_colors
    remaining_obs = color_sizes(i)

    do while (remaining_obs > 0)
      if (remaining_obs > colors%chunk_size) then
        chunks(chunk_index)%num_obs = colors%chunk_size
        chunks(chunk_index)%owner = MOD(chunk_index, numRanks)

        do j = 1, colors%chunk_size
          chunks(chunk_index)%obs_list(j) = ob_index + j 
        enddo

        remaining_obs = remaining_obs - colors%chunk_size
        chunk_index = chunk_index + 1
        ob_index = ob_index + colors%chunk_size
      else
        chunks(chunk_index)%num_obs = remaining_obs
        chunks(chunk_index)%owner = MOD(chunk_index, numRanks)

        do j = 1, remaining_obs
          chunks(chunk_index)%obs_list(j) = ob_index + j
        enddo

        remaining_obs = 0
        chunk_index = chunk_index + 1
        ob_index = ob_index + remaining_obs
      endif
    enddo
  enddo 

  ! Debug:
  write(*,*) "Total of colors, chunks : ", sum(color_sizes), chunk_count
  write(*,*) "Colors%chunk_size = ", colors%chunk_size
  
  do i = 1, size(chunks)
    write(*,*) "Chunk Assigment: ",i," -> ",chunks(i)%owner, chunks(i)%num_obs
  enddo

end subroutine create_chunks

!-------------------------------------------------------------

subroutine initialize_chunk_data(chunk_size, ens_size, num_groups, chunk_data)
  integer, intent(in) :: chunk_size
  integer, intent(in) :: ens_size
  integer, intent(in) :: num_groups
  type(chunk_data_type), intent(out) :: chunk_data

  integer :: buffer_size


  allocate(chunk_data%obs_prior(chunk_size, ens_size))
  allocate(chunk_data%obs_inc(chunk_size, ens_size))
  allocate(chunk_data%net_a(chunk_size, num_groups))
  allocate(chunk_data%obs_qc(chunk_size))
  allocate(chunk_data%vertvalue_obs_in_localization_coord(chunk_size))
  allocate(chunk_data%whichvert_real(chunk_size))

 write(*,*) "Size obs_prior : ", size(chunk_data%obs_prior)

  !allocate(chunk_data%base_obs_loc(chunk_size))

  buffer_size = (chunk_size * ens_size * 2) + (chunk_size * num_groups) + (chunk_size * 3) ! Plus base_obs_loc?
  allocate(chunk_data%bcast_buffer(buffer_size))

end subroutine initialize_chunk_data

!-------------------------------------------------------------

subroutine broadcast_send_chunk(from, chunk_data)
    use mpi
   integer, intent(in) :: from
   type(chunk_data_type) :: chunk_data

   integer :: start_offset = 1
   integer :: iError

!   chunk_data%bcast_buffer(start_offset:start_offset+size(chunk_data%obs_prior)) = reshape(chunk_data%obs_prior, [ size(chunk_data%obs_prior) ])
!   start_offset = start_offset + size(chunk_data%obs_prior)

  call MPI_Bcast(chunk_data%obs_prior, size(chunk_data%obs_prior), MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%obs_inc,   size(chunk_data%obs_inc),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%net_a,   size(chunk_data%net_a),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%obs_qc,   size(chunk_data%obs_qc),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%vertvalue_obs_in_localization_coord,   size(chunk_data%vertvalue_obs_in_localization_coord),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%whichvert_real,   size(chunk_data%whichvert_real),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)

end subroutine broadcast_send_chunk

!-------------------------------------------------------------

subroutine broadcast_recv_chunk(from, chunk_data)
    use mpi
   integer, intent(in) :: from
   type(chunk_data_type) :: chunk_data

   integer :: iError

!   write(*,*) "broadcast_recv_chunk - from = ", from
!  chunk_data%obs_prior = reshape(chunk_data%bcast_buffer(1:9), [ 3,3 ])

  call MPI_Bcast(chunk_data%obs_prior, size(chunk_data%obs_prior), MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%obs_inc,   size(chunk_data%obs_inc),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%net_a,   size(chunk_data%net_a),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%obs_qc,   size(chunk_data%obs_qc),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%vertvalue_obs_in_localization_coord,   size(chunk_data%vertvalue_obs_in_localization_coord),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)
  call MPI_Bcast(chunk_data%whichvert_real,   size(chunk_data%whichvert_real),   MPI_DOUBLE_PRECISION, from, MPI_COMM_WORLD, iError)

end subroutine broadcast_recv_chunk

!-------------------------------------------------------------

subroutine get_obs_from_color(colors, i, obs_list, last_rank)
  type(colors_type), intent(in) :: colors
  integer,           intent(in) :: i
  integer, dimension(:), allocatable, intent(inout) :: obs_list
  integer, intent(out) :: last_rank ! hack for now

  integer :: color_size
  integer :: chunks, chunks_per_rank, my_num_obs

  integer :: numRanks, mpiRank, iError
  integer :: k, remaining_obs
  integer, save :: current_rank = 0  ! Save so we resume where we left off at the next round

  integer, dimension(:), allocatable :: obs_per_rank

  integer :: start_ob, end_ob

  ! Deallocate if we'd already allocated:
  if (allocated(obs_list)) then
     deallocate(obs_list)
  endif

  ! Get the size of this color -- this is a bit hackish, would ANY work better?
  color_size = sum(colors%obs_color, colors%obs_color==i) / i

  ! Get the number of chunks
  chunks = (color_size + colors%chunk_size - 1) / colors%chunk_size

  ! Get the number of chunks per rank
  chunks_per_rank = (chunks + task_count() - 1) / task_count()
  if (chunks_per_rank > 1) then
       write(*,*) "Error - don't yet support >1 chunk per rank.. quitting."
       call MPI_Finalize(iError)
       stop
  endif


  ! Now allocate based on chunks_per_rank * chunk_size
  allocate(obs_list(chunks_per_rank * colors%chunk_size))

  ! set it to 0 (invalid ob): 
  obs_list = 0

  ! Now assign... for now, with no offset - always start at rank 0:
  ! Note: We're doing this the dumb way with a loop for now; there's a better
  ! way, but I don't want to waste time figuring it out until I have something
  ! working!

  remaining_obs = color_size
  numRanks = task_count()
  mpiRank = my_task_id()

  allocate(obs_per_rank(numRanks)) ! A value per rank so we know what each rank has
  obs_per_rank = 0

  ! Cheat - setting current_rank to 0 ensures we always start at rank 0... this
  ! is easier for now since we only will need to loop from 0->x, instead of
  ! (potentially) y->(some value less than y due to wrap-around in assignment).
  current_rank = 0

  do while (remaining_obs > 0) 
   if (remaining_obs > colors%chunk_size) then
      obs_per_rank(current_rank+1) = obs_per_rank(current_rank+1) + colors%chunk_size
      remaining_obs = remaining_obs - colors%chunk_size
      current_rank = MOD(current_rank + 1, numRanks)
    else
      obs_per_rank(current_rank+1) = obs_per_rank(current_rank+1) + remaining_obs
      remaining_obs = 0
      current_rank = MOD(current_rank + 1, numRanks)
    endif
  enddo

  last_rank = current_rank 

  ! debug printing:
  do k = 0, numRanks-1
     if (k == mpiRank) then
       !write(*,*) "Rank: ",mpiRank," has ",obs_per_rank(k+1)
     endif
  enddo

  ! Actual assignment:
  if (mpiRank < current_rank) then
    start_ob = sum(obs_per_rank(1:mpiRank)) + 1
    end_ob = sum(obs_per_rank(1:mpiRank+1))
    do k = 1, obs_per_rank(mpiRank+1)
       obs_list(k) = start_ob + k - 1
    enddo
  endif



  write(*,*) "C(",i,") Rank# ",mpiRank," Start/End => ",start_ob,end_ob,last_rank
  write(*,*) "C(",i,") Rank# ",mpiRank," LIST => ",obs_list(:)

end subroutine get_obs_from_color

!----------------------------------------------------------------------

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
!-----------------------------------------------------------
!> test get_state_meta_data
!> Write out the resutls of get_state_meta_data for each task
!> They should be the same as the Trunk version
subroutine test_get_state_meta_data(locations, num_vars)

type(location_type), intent(in) :: locations(:)
integer,             intent(in) :: num_vars

character*20  :: task_str !< string to hold the task number
character*129 :: file_meta !< output file name
character(len=128) :: locinfo
integer :: i

write(task_str, '(i10)') my_task_id()
file_meta = TRIM('test_get_state_meta_data' // TRIM(ADJUSTL(task_str)))

open(15, file=file_meta, status = 'unknown')

do i = 1, num_vars
   call write_location(-1, locations(i), charstring=locinfo)
   write(15,*) trim(locinfo)
enddo

close(15)


end subroutine test_get_state_meta_data

!--------------------------------------------------------
!> dump out the copies array for the state ens handle
subroutine test_state_copies(state_ens_handle, information)

type(ensemble_type), intent(in) :: state_ens_handle
character(len=*),        intent(in) :: information

character*20  :: task_str !< string to hold the task number
character*129 :: file_copies !< output file name
integer :: i

write(task_str, '(i10)') state_ens_handle%my_pe
file_copies = TRIM('statecopies_'  // TRIM(ADJUSTL(information)) // '.' // TRIM(ADJUSTL(task_str)))
open(15, file=file_copies, status ='unknown')

do i = 1, state_ens_handle%num_copies - state_ens_handle%num_extras
   write(15, *) state_ens_handle%copies(i,:)
enddo

close(15)

end subroutine test_state_copies

!--------------------------------------------------------
!> dump out the distances calculated in get_close_obs
subroutine test_close_obs_dist(distances, num_close, ob)

real(r8), intent(in) :: distances(:) !< array of distances calculated in get_close
integer,  intent(in) :: num_close !< number of close obs
integer,  intent(in) :: ob

character*20  :: task_str !< string to hold the task number
character*20  :: ob_str !< string to hold ob number
character*129 :: file_dist !< output file name
integer :: i

write(task_str, '(i10)') my_task_id()
write(ob_str, '(i20)') ob
file_dist = TRIM('distances'   // TRIM(ADJUSTL(task_str)) // '.' // TRIM(ADJUSTL(ob_str)))
open(15, file=file_dist, status ='unknown')

write(15, *) num_close

do i = 1, num_close
   write(15, *) distances(i)
enddo

close(15)

end subroutine test_close_obs_dist

!> @}

!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_graph_tools_mod

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/assimilation_code/modules/assimilation/assim_tools_mod.f90 $
! $Id: assim_tools_mod.f90 11799 2017-07-07 21:08:09Z nancy@ucar.edu $
! $Revision: 11799 $
! $Date: 2017-07-07 15:08:09 -0600 (Fri, 07 Jul 2017) $
