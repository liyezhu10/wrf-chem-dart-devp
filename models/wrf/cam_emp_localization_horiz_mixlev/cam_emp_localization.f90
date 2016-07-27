
PROGRAM emp_localization

use               types_mod, only : r8, missing_r8, DEG2RAD, PI, gravity, ps0,              &
                                    ps0, gas_constant, gas_constant_v, t_kelvin
use           utilities_mod, only : find_namelist_in_file, check_namelist_read, nc_check
use               model_mod, only : static_init_model, compute_geometric_height, toGrid,    &
                                    pres_to_zk
use               map_utils, only : proj_info, map_init, map_set, latlon_to_ij,             &
                                    ij_to_latlon, gridwind_to_truewind
use       mpi_utilities_mod, only : sleep_seconds
use   obs_def_altimeter_mod, only : compute_altimeter
use              meteor_mod, only : temp_and_dewpoint_to_rh, sat_vapor_pressure,            &
                                    specific_humidity
use             obs_err_mod, only : land_pres_error, land_temp_error, land_wind_error
use    dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, rh_error_from_dewpt_and_temp
use misc_definitions_module, only : PROJ_LC

use           netcdf

implicit none

! cam_emp_localization.f90_ELF is the code ofr 1 group with computation of ELF (including obs type 'TPW')
! 2012-08-03: update the code for global group filter 
!             for real obs experiment, there is no truth, thus only do computation of GGF, no ELF

! 2012-11-01: 
! 1. NOTE: subroutine 'find_index_in_locind' is updated, including 0distance, and w/o searching
!          also localization%locind is changed to include 0dist for index 1
! 2. read posterior files
!    use posterior of another group as an estimate of "true" for this group, when computing ELF

integer           :: emp_loc_kind     = 1
integer           :: locnum           = 160
real(r8)          :: locrad           = 0.02  ! 200 meter instead of radiance of 0.0002 (~1km) 
integer           :: group_size       = 2
integer           :: nbdy_nouse       = 5
integer           :: num_state_vars   = 7
integer           :: ens_size         = 80     ! change ens_copy when changing ens_size
integer           :: ens_copy         = 84     ! ens mean + ens spd + ens_size + inf mean + inf spd
integer           :: lat_s            = 96+33
integer           :: lat_e            = 96+75
! integer           :: lev_s            = 7
! integer           :: lev_e            = 7
integer           :: obslev           = 1      ! change obs_ind below
! integer           :: obslev_s         = 1
! integer           :: obslev_e         = 1
integer           :: varlev           = 1      ! change var_ind below
! integer           :: varlev_s         = 7
! integer           :: varlev_e         = 7
! dimlat, dimlon are used for subregion. (nlon/nslon=288)
integer           :: dimlat           = 1
integer           :: dimlon           = 8
integer           :: indobs_s         = 2
integer           :: indobs_e         = 4
integer           :: indvar_s         = 1
integer           :: indvar_e         = 5

real(r8)          :: R_P              = 40000.0_r8 ! obs error variance of Psfc
real(r8)          :: R_T              = 1.0_r8   ! obs error variance of Temperature
real(r8)          :: R_U              = 4.0_r8   ! obs error variance of Wind
real(r8)          :: R_TPW            = 0.015    ! 0.0146 computed from RMSE fit

real(r8)          :: pressure_top     = 20000.0  ! top pressure for 'TPW'

character(len=80)  :: truth_name_input  = 'truth.nc',        &
                      prior_name_input  = 'prior',           &
                      posterior_name_input = 'posterior',    &
                      obs_name_input    = 'obs.nc'
character(len=80), allocatable  :: prior_name_tmp(:)
character(len=80), allocatable  :: posterior_name_tmp(:)

type cam_data 
   integer  :: nlat, nlon, nslat, nslon, nlev, nilev         ! dimensions
   integer  :: ncopy                                         ! copies of truth and prior; nmhgt
   real(r8), dimension(:),  pointer        :: lat, lon, slat, slon
   real(r8), dimension(:),  pointer        :: lev, ilev, hyam, hybm, hyai, hybi, P0
   ! lat lon for each variable (US staggered on lat, VS staggered on lon), NOTE no variable staggered on ilev
   real(r8), dimension(:,:),  pointer      :: var_lat, var_lon 
   integer                                 :: number_of_state_variables
   integer, dimension(:,:), pointer        :: var_size
   character(len=129),dimension(:),pointer :: description
   real(r8), dimension(:,:,:,:,:), pointer :: variables 
   real(r8), dimension(:,:,:,:,:), pointer :: variables_post
end type cam_data

type(cam_data)    :: cam 


type obs_data
    integer    :: obsindex, obscopysize, obslocsize, obsqcsize, obstypesize, type_size
    integer    :: numobstype
    integer,  dimension(:), pointer   ::  obstype
    character(len=129), dimension(:), pointer :: obstype_metadata
    integer,  dimension(:,:), pointer :: obsqc
    integer,  dimension(:), pointer   :: which_vert
    real(r8), dimension(:,:), pointer :: location
    real(r8), dimension(:,:), pointer :: observations
    integer,  dimension(:), pointer   :: obstypeindex

    integer,  dimension(:,:), pointer :: obs_size
    real(r8), dimension(:,:,:,:,:), pointer :: obs_val
    real(r8), dimension(:),   pointer :: obsvar
end type obs_data

type(obs_data)   :: obs


type localization_data
    real(r8), dimension(:),         pointer  :: locind
    integer,  dimension(:,:,:,:),   pointer  :: num
    real(r8), dimension(:,:,:,:,:), pointer  :: correlation
    real(r8), dimension(:,:,:,:,:), pointer  :: sumcorr
    real(r8), dimension(:,:,:,:,:), pointer  :: sumcorrsq

    ! for GGF_reg and DGF_reg
    real(r8), dimension(:,:,:,:),   pointer   :: beta1           ! (sum(beta_reg))**2
    real(r8), dimension(:,:,:,:),   pointer   :: beta2           ! sum(beta_reg**2)
    real(r8), dimension(:,:,:,:),   pointer   :: meangf

    ! for GGF_inc, DGF_inc, ELF with cross observation (unknown "true")
    integer,  dimension(:,:,:,:),   pointer   :: num_cobs
    real(r8), dimension(:,:,:,:),   pointer   :: incbeta1_cobs           ! (sum(beta_increment))**2
    real(r8), dimension(:,:,:,:),   pointer   :: incbeta2_cobs           ! sum(beta_increment**2)
    real(r8), dimension(:,:,:,:),   pointer   :: incmeangf_cobs
    real(r8), dimension(:,:,:,:,:), pointer   :: sumx_cobs               ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:,:,:), pointer   :: sumy_cobs               ! used to get mean(y)
    real(r8), dimension(:,:,:,:,:), pointer   :: sumx2_cobs              ! used to get mean(x^2)
    real(r8), dimension(:,:,:,:,:), pointer   :: sumy2_cobs              ! used to get mean(y^2)
    real(r8), dimension(:,:,:,:,:), pointer   :: numerator_cobs          ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:,:,:), pointer   :: denominator_cobs        ! i.e., used to get mean(y^2) (note: it contains Robs)

    ! for GGF_inc, DGF_inc, ELF with fake perfect observation (unknown "true", Robs = 0)
    integer,  dimension(:,:,:,:),   pointer   :: num_pobs
    real(r8), dimension(:,:,:,:),   pointer   :: incbeta1_pobs           ! (sum(beta_increment))**2
    real(r8), dimension(:,:,:,:),   pointer   :: incbeta2_pobs           ! sum(beta_increment**2)
    real(r8), dimension(:,:,:,:),   pointer   :: incmeangf_pobs
    real(r8), dimension(:,:,:,:,:), pointer   :: sumx_pobs               ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:,:,:), pointer   :: sumy_pobs               ! used to get mean(y)
    real(r8), dimension(:,:,:,:,:), pointer   :: sumx2_pobs              ! used to get mean(x^2)
    real(r8), dimension(:,:,:,:,:), pointer   :: sumy2_pobs              ! used to get mean(y^2)
    real(r8), dimension(:,:,:,:,:), pointer   :: numerator_pobs          ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:,:,:), pointer   :: denominator_pobs        ! i.e., used to get mean(y^2)

end type localization_data

type(localization_data) localization

type localization_omp_data
  type(localization_data), pointer  :: obs(:)
  integer                           :: nobs
end type localization_omp_data

type(localization_omp_data) loc_for_omp


integer          :: i, j, k, id, idobs, ind, inddim, indens,               &
                    ii, jj, kk, indk1, indk2, var_id, ndims, dimids(10),   &
                    obslevk, varlevk, ig, igs, ige, igg
integer          :: obs_lat_s, obs_lat_e, obs_lon_s, obs_lon_e,            &
                    var_lat_s, var_lat_e, var_lon_s, var_lon_e
integer          :: dim_region, dim_ind, dim_lat_interval, dim_lon_interval, ireg
integer          :: ncid, ncidtruth, ncidprior, ncidobs, iunit
integer          :: proj_code
integer          :: ilocind, sublocnum, writelenloc
integer          :: obs_ind(4)    ! for mixlev us ( lev1, 7, 12, 25)
integer          :: var_ind(4,4)  ! for mixlev use (lev1, 7, 12, 25)

real(r8)         :: stdlon,truelat1,truelat2,latinc,loninc
real(r8)         :: obslat, obslon, obspres, obsi, obsj, obsk, obspalt, obshgt
real(r8)         :: varlat, varlon
!real(r8)         :: yobs, ytrue, yfmean, oerr, Robs, yprior_var, ypost_var, yvar, yinno_norm,              &
!                    xfmean, xtrue, xprior_var,xt_xfmean, covxy, deltaxmean, ccc, ddd, dist
real(r8)         :: yobs, ytrue, oerr, Robs, yvar, yinno_norm, xtrue, xt_xfmean, deltaxmean, dist
real(r8)         :: obsdiag1, obsdiag2, obsdiag3, obsdiag4
real(r8)         :: sumbeta1, sumbeta2, incsumbeta1, incsumbeta2
real(r8)         :: ref_P, surf_P

logical          :: is_lev0
character(len=2) :: idchar1, idchar2
character(len=1) :: igchar
character(len=NF90_MAX_NAME) :: var_name_tmp, name
character(len=100) :: obsnum_format, locint_format, locreal_format, mhgt_format

integer,  allocatable :: ncidpriors(:), ncidposteriors(:)
integer,  allocatable :: dim_latlon(:,:)
real(r8), allocatable :: cam_var_1d(:), cam_var_2d(:,:), cam_var_3d(:,:,:), cam_var_4d(:,:,:,:),           &
                                                         cam_var_3d_2(:,:,:)
real(r8), allocatable :: cam_hyam(:), cam_hybm(:), qv_1d(:)
real(r8), allocatable :: obs_var_1d(:), obs_var_2d(:,:), obs_var_3d(:,:,:), obs_var_4d(:,:,:,:)
real(r8), allocatable :: wrf_mu_d01(:,:), wrf_psfc_d01(:,:), wrf_ph_d01(:,:,:),       &
                         wrf_t_d01(:,:,:), wrf_qvapor_d01(:,:,:)
real(r8), allocatable :: xprior(:), yprior(:)
real(r8), allocatable :: xtrueprf(:), xfmeanprf(:), xfspdprf(:), xpriorprf(:,:)
real(r8), allocatable :: xtrue3d(:,:,:), xfmean3d(:,:,:), xfspd3d(:,:,:), xprior3d(:,:,:,:) 
real(r8), allocatable :: numtmp(:)
real(r8), allocatable :: yfmean(:), yamean(:), yprior_var(:), ypost_var(:)
real(r8), allocatable :: xfmean(:), xamean(:), xprior_var(:), xpost_var(:)
real(r8), allocatable :: covxy(:), corr(:), ccc(:), ddd(:), reg_coef(:)
real(r8), allocatable :: incnumer(:), incdenom(:)
real(r8), allocatable :: sublocx(:), sublocy(:), sublocx2(:), sublocy2(:),           &
                         sublocnumer(:), sublocdenom(:), sublocalpha(:)


! define format to write "localization"
locint_format = ' ( 2000i10 ) '
locreal_format = ' ( 2000e15.7 ) '
! define format to write obs number
!obsnum_format = ' ( i10,4e15.7 )'
obsnum_format = ' ( i10 ) '
mhgt_format = ' ( 49e15.7 ) '

obs_ind(1) = 1

!var_ind(1,1) = 1
var_ind(1,1) = 1   ! maxi veritcal weight for TPW at model level 8

!!-------------------------------------------------------
!! open obs file (obs_epoch.nc) 
!!-------------------------------------------------------
!   print *, 'starting reading obs file'
!   call nc_check( nf90_open(obs_name_input, NF90_NOWRITE, ncidobs),       &
!                  'wrf_emp_localization','obs file')
!
!!-------------------------------------------------------
!! read obs file dimensions
!!-------------------------------------------------------
!   call nc_check( nf90_inq_dimid(ncidobs, "ObsIndex", var_id),       &
!                  'read_obs_dimensions','inq_dimid ObsIndex')
!   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obsindex),  &
!                  'read_obs_dimensions','inquire_dimension'//trim(name))
!
!   call nc_check( nf90_inq_dimid(ncidobs, "copy", var_id),       &
!                  'read_obs_dimensions','inq_dimid copy')
!   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obscopysize),  &
!                  'read_obs_dimensions','inquire_dimension'//trim(name))
!
!   call nc_check( nf90_inq_dimid(ncidobs, "qc_copy", var_id),       &
!                  'read_obs_dimensions','inq_dimid qc_copy')
!   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obsqcsize),  &
!                  'read_obs_dimensions','inquire_dimension'//trim(name))
!
!   call nc_check( nf90_inq_dimid(ncidobs, "location", var_id),       &
!                  'read_obs_dimensions','inq_dimid location')
!   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obslocsize),  &
!                  'read_obs_dimensions','inquire_dimension'//trim(name))
!
!   call nc_check( nf90_inq_dimid(ncidobs, "ObsTypes", var_id),       &
!                  'read_obs_dimensions','inq_dimid ObsTypes')
!   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obstypesize),  &
!                  'read_obs_dimensions','inquire_dimension'//trim(name))
!
!!-------------------------------------------------------
!! read obs file variable
!!-------------------------------------------------------
!   allocate(obs%obstype_metadata(obs%obstypesize))
!   call nc_check( nf90_inq_varid(ncidobs, "ObsTypesMetaData", var_id), &
!                     'read_obs_variable','inq_varid ObsTypesMetaData')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%obstype_metadata), &
!                     'read_obs_variable','get_var ObsTypesMetaData')
!
!   allocate(obs%obstype(obs%obsindex))
!   call nc_check( nf90_inq_varid(ncidobs, "obs_type", var_id), &
!                     'read_obs_variable','inq_varid obs_type')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%obstype), &
!                     'read_obs_variable','get_var obs_type')
!
!   ! NOTE: nf90_get_var has the dimension inversely after read in than in the original .nc file
!   allocate(obs%obsqc(obs%obsqcsize,obs%obsindex))
!   call nc_check( nf90_inq_varid(ncidobs, "qc", var_id), &
!                     'read_obs_variable','inq_varid qc')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%obsqc), &
!                     'read_obs_variable','get_var qc')
!
!   allocate(obs%which_vert(obs%obsindex))
!   call nc_check( nf90_inq_varid(ncidobs, "which_vert", var_id), &
!                     'read_obs_variable','inq_varid which_vert')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%which_vert), &
!                     'read_obs_variable','get_var which_vert')
!
!   allocate(obs%location(obs%obslocsize,obs%obsindex))
!   call nc_check( nf90_inq_varid(ncidobs, "location", var_id), &
!                     'read_obs_variable','inq_varid location')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%location), &
!                     'read_obs_variable','get_var location')
!
!   allocate(obs%observations(obs%obscopysize,obs%obsindex))
!   call nc_check( nf90_inq_varid(ncidobs, "observations", var_id), &
!                     'read_obs_variable','inq_varid observations')  ! reuse var_id, no harm
!   call nc_check( nf90_get_var(ncidobs, var_id, obs%observations), &
!                     'read_obs_variable','get_var observations')


!-------------------------------------------------------
! open truth file (True_State.nc) 
!-------------------------------------------------------
! Lili test
!print *, 'starting reading truth and prior file'
!call nc_check( nf90_open(truth_name_input, NF90_NOWRITE, ncidtruth),   &
!               'cam_emp_localization','open truth file')

!-------------------------------------------------------
! open prior file (Prior_Diag.nc) 
!-------------------------------------------------------
allocate(prior_name_tmp(group_size))
allocate(ncidpriors(group_size))
do i = 1, group_size
   write(igchar,'(i1.1)') i
   prior_name_tmp(i) = trim(prior_name_input)//igchar//'.nc'
   call nc_check( nf90_open(trim(prior_name_tmp(i)), NF90_NOWRITE, ncidpriors(i)),   &
                  'cam_emp_localization','open prior file')
enddo

!-------------------------------------------------------
! open posterior file (Posterior_Diag.nc) 
!-------------------------------------------------------
allocate(posterior_name_tmp(group_size))
allocate(ncidposteriors(group_size))
do i = 1, group_size
   write(igchar,'(i1.1)') i
   posterior_name_tmp(i) = trim(posterior_name_input)//igchar//'.nc'
   call nc_check( nf90_open(trim(posterior_name_tmp(i)), NF90_NOWRITE, ncidposteriors(i)),   &
                  'cam_emp_localization','open posterior file')
enddo

!-------------------------------------------------------
! fill in state variable information 
!-------------------------------------------------------
   cam%number_of_state_variables = num_state_vars
   allocate(cam%description(cam%number_of_state_variables))
   cam%description(1)  = 'PS'
   cam%description(2)  = 'T'
   cam%description(3)  = 'US'
   cam%description(4)  = 'VS'
   cam%description(5)  = 'Q'
   cam%description(6)  = 'CLDLIQ'
   cam%description(7)  = 'CLDICE'


!-------------------------------------------------------
! read CAM dimensions
!-------------------------------------------------------
   ! dimension of copies (2(ens mean, ens spread)+ens_size)
   ! NOTE: (inflation is turned off in the OSSE)
   cam%ncopy = ens_copy * group_size

   var_name_tmp = 'lat'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid lat')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nlat),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'slat'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nslat')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nslat),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'lon'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid lon')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nlon),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'slon'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nslon')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nslon),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'lev'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nlev')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nlev),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'ilev'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nilev')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name, cam%nilev),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

!-------------------------------------------------------
! read CAM dimension data 
!-------------------------------------------------------
! 1. 1D array

   allocate(cam%lat(cam%nlat))
   var_name_tmp = 'lat'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lat')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%lat),  &
                     'read_cam_dimension_data','get_var lat')

   allocate(cam%slat(cam%nslat))
   var_name_tmp = 'slat'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid slat')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%slat),  &
                     'read_cam_dimension_data','get_var slat')

   allocate(cam%lon(cam%nlon))
   var_name_tmp = 'lon'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lon')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%lon),  &
                     'read_cam_dimension_data','get_var lon')

   allocate(cam%slon(cam%nslon))
   var_name_tmp = 'slon'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid slon')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%slon),  &
                     'read_cam_dimension_data','get_var slon')

   allocate(cam%lev(cam%nlev))
   var_name_tmp = 'lev'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lev')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%lev),  &
                     'read_cam_dimension_data','get_var lev')

   allocate(cam%ilev(cam%nilev))
   var_name_tmp = 'ilev'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid ilev')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%ilev),  &
                     'read_cam_dimension_data','get_var ilev')

   allocate(cam%hyam(cam%nlev))
   var_name_tmp = 'hyam'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid hyam')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%hyam),  &
                     'read_cam_dimension_data','get_var hyam')
   ! reshape data from top-bottom to bottom-top
   allocate(cam_var_1d(cam%nlev))
   cam_var_1d = cam%hyam
   do i = 1, cam%nlev
      cam%hyam(i) = cam_var_1d(cam%nlev-i+1)
   enddo
   deallocate(cam_var_1d)

   allocate(cam%hybm(cam%nlev))
   var_name_tmp = 'hybm'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid hybm')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%hybm),  &
                     'read_cam_dimension_data','get_var hybm')
   ! reshape data from top-bottom to bottom-top
   allocate(cam_var_1d(cam%nlev))
   cam_var_1d = cam%hybm
   do i = 1, cam%nlev
      cam%hybm(i) = cam_var_1d(cam%nlev-i+1)
   enddo
   deallocate(cam_var_1d)

   allocate(cam%hyai(cam%nilev))
   var_name_tmp = 'hyai'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid hyai')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%hyai),  &
                     'read_cam_dimension_data','get_var hyai')
   ! reshape data from top-bottom to bottom-top
   allocate(cam_var_1d(cam%nilev))
   cam_var_1d = cam%hyai
   do i = 1, cam%nilev
      cam%hyai(i) = cam_var_1d(cam%nilev-i+1)
   enddo
   deallocate(cam_var_1d)

   allocate(cam%hybi(cam%nilev))
   var_name_tmp = 'hybi'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid hybi')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%hybi),  &
                     'read_cam_dimension_data','get_var hybi')
   ! reshape data from top-bottom to bottom-top
   allocate(cam_var_1d(cam%nilev))
   cam_var_1d = cam%hybi
   do i = 1, cam%nilev
      cam%hybi(i) = cam_var_1d(cam%nilev-i+1)
   enddo
   deallocate(cam_var_1d)

   allocate(cam%P0(cam%nlev))
   var_name_tmp = 'P0'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid P0')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, cam%P0),  &
                     'read_cam_dimension_data','get_var P0')

!-------------------------------------------------------
! read state variables of True_State.nc (copy 1) 
!-------------------------------------------------------
   ! var_size(1,:) - dimension of lon, (i) 
   ! var_size(2,:) - dimension of lat, (j)
   ! var_size(3,:) - dimension of lev  (k) ( = 1, means 2D variable, others 3D variable)

   ! variables dimensions: 1 - lat, 2 - lon, 3 - lev, 
   !                       4 - copies, 5 - state variables

   allocate(cam%var_size(3,cam%number_of_state_variables))
   allocate(cam%variables(cam%nlon,cam%nlat,cam%nilev,cam%ncopy,    &
                             cam%number_of_state_variables))
   allocate(cam%variables_post(cam%nlon,cam%nlat,cam%nilev,cam%ncopy,    &
                             cam%number_of_state_variables))
   allocate(cam%var_lat(cam%nlat,cam%number_of_state_variables))
   allocate(cam%var_lon(cam%nlon,cam%number_of_state_variables))
   
   do ind = 1,cam%number_of_state_variables

      ! Get the dimension size 
      ! Once this is done for variable of True_State.nc, there is no need to do this for
      ! variables in Prior_Diag.nc, since they should be consistent
      var_name_tmp = trim(cam%description(ind))
      call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'get_variable_size_from_file',                    &
                     'inq_varid '//var_name_tmp)
      call nc_check( nf90_inquire_variable(ncidpriors(1), var_id,          &
                     ndims=ndims, dimids=dimids),                      &
                     'get_variable_size_from_file',                    &
                     'inquire_variable '//var_name_tmp)
      do inddim = 1,ndims-2
         call nc_check( nf90_inquire_dimension(ncidpriors(1), dimids(inddim),  &
                        len=cam%var_size(inddim,ind)),                     &
                        'get_variable_size_from_file',                     &
                        'inquire_dimension '//var_name_tmp)
      enddo

      if ( ndims < 5 ) then

            ! 2D variable, its vertical dimension is 1
            cam%var_size(ndims-1,ind) = 1

      else

            ! reshape var_size
            j = cam%var_size(1,ind)
            cam%var_size(1,ind) = cam%var_size(2,ind)
            cam%var_size(2,ind) = cam%var_size(3,ind)
            cam%var_size(3,ind) = j

      endif

      ! get the lat, lon ready for each variable (US, VS are staggered)
      if ( var_name_tmp == 'US' ) then
         cam%var_lat(1:cam%nslat,ind) = cam%slat(1:cam%nslat)
         cam%var_lon(1:cam%nlon,ind)  = cam%lon(1:cam%nlon)

      elseif ( var_name_tmp == 'VS' ) then
         cam%var_lat(1:cam%nlat,ind)  = cam%lat(1:cam%nlat)
         cam%var_lon(1:cam%nslon,ind) = cam%slon(1:cam%nslon)
      else
         cam%var_lat(1:cam%nlat,ind)  = cam%lat(1:cam%nlat)
         cam%var_lon(1:cam%nlon,ind)  = cam%lon(1:cam%nlon)
      endif

   enddo   ! end ind for number_of_state_variables


!-------------------------------------------------------
! read state variables of Prior_Diag.nc (copy ens_mean, ens_spd, ens_member) 
!-------------------------------------------------------

   do ind = 1,cam%number_of_state_variables

      do i = 1, group_size

         igs = (i-1)*ens_copy + 1
         ige = i*ens_copy

         if ( i == group_size .and. ige /= cam%ncopy ) then
            print *, 'ERROR in read pirors, index WRONG'
         endif
   
         var_name_tmp = trim(cam%description(ind))

         call nc_check( nf90_inq_varid(ncidpriors(i), var_name_tmp, var_id),  &
                        'get_variable_name_from_file',                        &
                        'inq_varid '//var_name_tmp)

         if ( cam%var_size(3,ind) == 1 ) then

            ! Get 2D variable 
            ! NOTE: nf90_get_var automatically ignor the dimension with length 1
            allocate(cam_var_3d(cam%var_size(1,ind),                          &
                                cam%var_size(2,ind),                          &
                                ens_copy))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, cam_var_3d),   &
                        'read_in_cam_var','get_var '//var_name_tmp )
            cam%variables(1:cam%var_size(1,ind),                              &
                          1:cam%var_size(2,ind),                              &
                          1,                                                  &
                          igs:ige,                                            &
                          ind) =                                              &
                       cam_var_3d(1:cam%var_size(1,ind),                      &
                                  1:cam%var_size(2,ind),                      &
                                  1:ens_copy)
            deallocate(cam_var_3d)

         else

            ! Get 3D variable
            ! NOTE: var_size was reshaped when read in True_State.nc,
            !       but nc_varget starts from lev. this is how to set cam_var_4d
            allocate(cam_var_4d(cam%var_size(3,ind),                          &
                                cam%var_size(1,ind),                          &
                                cam%var_size(2,ind),                          &
                                ens_copy))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, cam_var_4d),   &
                        'read_in_cam_var','get_var '//var_name_tmp )
            do k = 1, cam%var_size(3,ind)
               cam%variables(1:cam%var_size(1,ind),                           &
                             1:cam%var_size(2,ind),                           &
                             cam%var_size(3,ind)-k+1,                         &
                             igs:ige,                                         &
                             ind) =                                           &
                          cam_var_4d(k,                                       &
                                     1:cam%var_size(1,ind),                   &
                                     1:cam%var_size(2,ind),                   &
                                     1:ens_copy)
            enddo
            deallocate(cam_var_4d)

         endif

      enddo   ! i for group_size

   enddo   ! ind for number_of_state_variables

!-------------------------------------------------------
! read state variables of Posterior_Diag.nc (copy ens_mean, ens_spd, ens_member) 
!-------------------------------------------------------

   do ind = 1,cam%number_of_state_variables

      do i = 1, group_size

         igs = (i-1)*ens_copy + 1
         ige = i*ens_copy

         if ( i == group_size .and. ige /= cam%ncopy ) then
            print *, 'ERROR in read pirors, index WRONG'
         endif

         var_name_tmp = trim(cam%description(ind))

         call nc_check( nf90_inq_varid(ncidposteriors(i), var_name_tmp, var_id),  &
                        'get_variable_name_from_file',                        &
                        'inq_varid '//var_name_tmp)

         if ( cam%var_size(3,ind) == 1 ) then

            ! Get 2D variable 
            ! NOTE: nf90_get_var automatically ignor the dimension with length 1
            allocate(cam_var_3d(cam%var_size(1,ind),                          &
                                cam%var_size(2,ind),                          &
                                ens_copy))
            call nc_check( nf90_get_var(ncidposteriors(i), var_id, cam_var_3d),   &
                        'read_in_cam_var','get_var '//var_name_tmp )
            cam%variables_post(1:cam%var_size(1,ind),                             &
                          1:cam%var_size(2,ind),                              &
                          1,                                                  &
                          igs:ige,                                            &
                          ind) =                                              &
                       cam_var_3d(1:cam%var_size(1,ind),                      &
                                  1:cam%var_size(2,ind),                      &
                                  1:ens_copy)
            deallocate(cam_var_3d)

         else

            ! Get 3D variable
            ! NOTE: var_size was reshaped when read in True_State.nc,
            !       but nc_varget starts from lev. this is how to set cam_var_4d
            allocate(cam_var_4d(cam%var_size(3,ind),                          &
                                cam%var_size(1,ind),                          &
                                cam%var_size(2,ind),                          &
                                ens_copy))
            call nc_check( nf90_get_var(ncidposteriors(i), var_id, cam_var_4d),   &
                        'read_in_cam_var','get_var '//var_name_tmp )
            do k = 1, cam%var_size(3,ind)
               cam%variables_post(1:cam%var_size(1,ind),                          &
                             1:cam%var_size(2,ind),                           &
                             cam%var_size(3,ind)-k+1,                         &
                             igs:ige,                                         &
                             ind) =                                           &
                          cam_var_4d(k,                                       &
                                     1:cam%var_size(1,ind),                   &
                                     1:cam%var_size(2,ind),                   &
                                     1:ens_copy)
            enddo
            deallocate(cam_var_4d)

         endif

      enddo   ! i for group_size

   enddo   ! ind for number_of_state_variables


! ----------------------------------------------------
! close netcdf files
! ----------------------------------------------------
   do i = 1, group_size
      call nc_check(nf90_close(ncidpriors(i)),'cam_emp_localization','close prior file')
   enddo
   deallocate(prior_name_tmp, ncidpriors)

  do i = 1, group_size
      call nc_check(nf90_close(ncidposteriors(i)),'cam_emp_localization','close posterior file')
   enddo
   deallocate(posterior_name_tmp, ncidposteriors)

!-------------------------------------------------------
! define obs type and obs error variance 
!-------------------------------------------------------

   obs%type_size = 5
   allocate(obs%obstype_metadata(obs%type_size))
   obs%obstype_metadata(1)  = 'PRESSURE'
   obs%obstype_metadata(2)  = 'TEMPERATURE'
   obs%obstype_metadata(3)  = 'U_WIND_COMPONENT'
   obs%obstype_metadata(4)  = 'V_WIND_COMPONENT'
   obs%obstype_metadata(5)  = 'TPW'

   allocate(obs%obsvar(obs%type_size))
   obs%obsvar(1) = R_P
   obs%obsvar(2) = R_T
   obs%obsvar(3) = R_U
   obs%obsvar(4) = R_U
   obs%obsvar(5) = R_TPW


allocate(cam_hyam(cam%nlev))
allocate(cam_hybm(cam%nlev))
cam_hyam = cam%hyam
cam_hybm = cam%hybm

dim_region = dimlat*dimlon
allocate(dim_latlon(4,dim_region))

! Initialize locind which contains the radiance for each subset
!   print *, 'empirical localization kind: horizontal'
   allocate(localization%locind(locnum))
   do ind = 1, locnum
      if ( ind == 1 ) then
         localization%locind(ind) = 0.0_r8
      else
         localization%locind(ind) = locrad * real(ind-1)
      endif
   enddo

LoopObsTypes: do idobs = indobs_s, indobs_e
   print *, 'Now processing obs type ', trim(obs%obstype_metadata(idobs))

   if ( idobs == 3 ) then
      ! U stagger on lat cam%nslat
      obs_lat_s = lat_s
      obs_lat_e = lat_e
      obs_lon_s = 1
      obs_lon_e = cam%nlon
   elseif ( idobs == 4 ) then
      ! V stagger on lon cam%nslon
      obs_lat_s = lat_s
      obs_lat_e = lat_e
      obs_lon_s = 1
      obs_lon_e = cam%nslon
   else
      obs_lat_s = lat_s
      obs_lat_e = lat_e
      obs_lon_s = 1
      obs_lon_e = cam%nlon
   endif

   ! assign (lat lon) index for subregions
   ! dim_latlon(4,dim_region), 4 for lat_beg, lat_end, lon_beg, lon_end
   dim_lat_interval = (obs_lat_e-obs_lat_s+1)/dimlat
   dim_lon_interval = (obs_lon_e-obs_lon_s+1)/dimlon
   do i = 1, dimlat
      do j = 1, dimlon
         dim_ind = (i-1)*dimlon + j
         dim_latlon(1,dim_ind) = (i-1)*dim_lat_interval+obs_lat_s
         dim_latlon(2,dim_ind) = i*dim_lat_interval+obs_lat_s-1
         dim_latlon(3,dim_ind) = (j-1)*dim_lon_interval+1
         dim_latlon(4,dim_ind) = j*dim_lon_interval
         if ( i == dimlat ) then
            dim_latlon(2,dim_ind) = max(dim_latlon(2,dim_ind),obs_lat_e)
         endif
         if ( j == dimlon ) then
            dim_latlon(4,dim_ind) = max(dim_latlon(4,dim_ind),obs_lon_e)
         endif
      enddo
   enddo
!   do i = 1, dim_region
!      print *, dim_latlon(1:4,i)
!   enddo
!   pause


! allocate and initialize localization variables
   ndims = obslev * varlev
   allocate(localization%num(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%correlation(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumcorr(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumcorrsq(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))

   allocate(localization%beta1(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%beta2(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%meangf(locnum,dim_region,ndims,cam%number_of_state_variables))

   allocate(localization%num_cobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incbeta1_cobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incbeta2_cobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incmeangf_cobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumx_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumy_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumx2_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumy2_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%numerator_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%denominator_cobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))

   allocate(localization%num_pobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incbeta1_pobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incbeta2_pobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%incmeangf_pobs(locnum,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumx_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumy_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumx2_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%sumy2_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%numerator_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))
   allocate(localization%denominator_pobs(locnum,group_size,dim_region,ndims,cam%number_of_state_variables))

   localization%num                     = 0
   localization%correlation             = 0.0_r8
   localization%sumcorr                 = 0.0_r8
   localization%sumcorrsq               = 0.0_r8

   localization%beta1                   = 0.0_r8
   localization%beta2                   = 0.0_r8
   localization%meangf                  = 0.0_r8

   localization%num_cobs                = 0
   localization%incbeta1_cobs           = 0.0_r8
   localization%incbeta2_cobs           = 0.0_r8
   localization%incmeangf_cobs          = 0.0_r8
   localization%sumx_cobs               = 0.0_r8
   localization%sumy_cobs               = 0.0_r8
   localization%sumx2_cobs              = 0.0_r8
   localization%sumy2_cobs              = 0.0_r8
   localization%numerator_cobs          = 0.0_r8
   localization%denominator_cobs        = 0.0_r8

   localization%num_pobs                = 0
   localization%incbeta1_pobs           = 0.0_r8
   localization%incbeta2_pobs           = 0.0_r8
   localization%incmeangf_pobs          = 0.0_r8
   localization%sumx_pobs               = 0.0_r8
   localization%sumy_pobs               = 0.0_r8
   localization%sumx2_pobs              = 0.0_r8
   localization%sumy2_pobs              = 0.0_r8
   localization%numerator_pobs          = 0.0_r8
   localization%denominator_pobs        = 0.0_r8

!$OMP PARALLEL DO PRIVATE(ind,i,j,k,ii,jj,kk,obslevk,varlevk,inddim,ndims,indk1,indk2,ig,igs,ige,igg,ireg,iunit,writelenloc,ref_P,surf_P,var_name_tmp,obslon,obslat,varlon,varlat,obsi,obsj,obsk,obspres,obshgt,is_lev0,yobs,ytrue,yfmean,yamean,Robs,yvar,yprior_var,ypost_var,corr,reg_coef,ccc,ddd,dist,yinno_norm,xtrue,xfmean,xamean,xprior_var,xpost_var,xt_xfmean,covxy,deltaxmean,ilocind,sublocx,sublocy,sublocx2,sublocy2,sublocnum,sublocnumer,sublocalpha,sublocdenom,incnumer,incdenom,cam_var_1d,obs_var_1d,qv_1d,xprior,yprior,xtrueprf,xfmeanprf,xfspdprf,xpriorprf,xtrue3d,xfmean3d,xfspd3d,xprior3d,sumbeta1,sumbeta2)

!   LoopVarTypes: do inddim = 1, cam%number_of_state_variables
   LoopVarTypes: do inddim = indvar_s, indvar_e
      print *, 'processing variable ', trim(cam%description(inddim))

      if ( inddim == 3 ) then
         ! U stagger on lat cam%nslat
         var_lat_s = lat_s
         var_lat_e = lat_e
         var_lon_s = 1
         var_lon_e = cam%nlon
      elseif ( inddim == 4 ) then
         ! V stagger on lon cam%nslon
         var_lat_s = lat_s
         var_lat_e = lat_e
         var_lon_s = 1
         var_lon_e = cam%nslon
      else
         var_lat_s = lat_s
         var_lat_e = lat_e
         var_lon_s = 1
         var_lon_e = cam%nlon
      endif

      allocate(yprior(ens_size))
      allocate(xprior(ens_size))

      allocate(xfmean(group_size),xamean(group_size),xprior_var(group_size),xpost_var(group_size))
      allocate(yfmean(group_size),yamean(group_size),yprior_var(group_size),ypost_var(group_size))
      allocate(covxy(group_size),corr(group_size),reg_coef(group_size))
      allocate(ccc(group_size),ddd(group_size))
      allocate(sublocx(group_size),sublocy(group_size))
      allocate(sublocx2(group_size),sublocy2(group_size))
      allocate(sublocnumer(group_size),sublocdenom(group_size),sublocalpha(group_size))
      allocate(incnumer(group_size),incdenom(group_size))


      ! loop obs grid points
!      do k = lev_s, lev_e
      do obslevk = 1, obslev
         k = obs_ind(obslevk)
! 2012-08-06: add PS even PS is not at variable level
!             treat PS as integral variable
         if ( idobs == 1 ) then 
              k = 1
         endif

         ! loop target variable grid points
!         do kk = lev_s, lev_e
         do varlevk = 1, varlev
            kk = var_ind(obslevk,varlevk)
! 2012-08-06: add PS even PS is not at variable level
!             treat PS as integral variable
            if ( inddim == 1 ) then
                 kk = 1
            endif

            ndims = (obslevk-1)*varlev + varlevk

            do ireg = 1, dim_region

!            do j = lat_s, lat_e
!               do i = 1, cam%nlon
            do j = dim_latlon(1,ireg), dim_latlon(2,ireg)
               do i = dim_latlon(3,ireg), dim_latlon(4,ireg)

                  do jj = lat_s, lat_e
                     do ii = 1, cam%nlon

                        ! computation for each group
                        do ig = 1, group_size

                           igs = (ig-1)*ens_copy + 1
                           ige = ig*ens_copy

                           ! get obs information ready
                           if ( trim(obs%obstype_metadata(idobs)) == 'TPW' ) then
!                              ! compute 'TPW' truth
!                              allocate(cam_var_1d(cam%nlev))
!                              allocate(qv_1d(cam%nlev))
!                              ref_P  = cam%P0(1)
!                              surf_P = cam%variables(i,j,1,1,1)
!                              call get_pressure(ref_P, surf_P, cam%nlev,                     &
!                                                cam_hyam, cam_hybm, cam_var_1d)
!                              qv_1d = cam%variables(i,j,:,1,5)
!                              call get_expected_tpw(cam_var_1d, qv_1d, cam%nlev,             &
!                                                    ytrue, pressure_top)
!                              deallocate(cam_var_1d, qv_1d)

                              ! compute 'TPW' prior   
                              do indk1 = 1, ens_size
                                 allocate(cam_var_1d(cam%nlev))
                                 allocate(qv_1d(cam%nlev))
                                 ref_P  = cam%P0(1)
                                 surf_P = cam%variables(i,j,1,igs+1+indk1,1)
                                 call get_pressure(ref_P, surf_P, cam%nlev,                  &
                                                   cam_hyam, cam_hybm, cam_var_1d)
                                 qv_1d =  cam%variables(i,j,:,igs+1+indk1,5)
                                 call get_expected_tpw(cam_var_1d, qv_1d, cam%nlev,          &
                                                       yprior(indk1), pressure_top)
                                 deallocate(cam_var_1d, qv_1d)
                              enddo
                              yfmean(ig) = sum(yprior)/ens_size
                              yprior_var(ig) = 0.0_r8
                              do indk1 = 1, ens_size
                                 yprior_var(ig) = yprior_var(ig) + (yprior(indk1)-yfmean(ig))**2
                              enddo
                              yprior_var(ig) = yprior_var(ig)/(ens_size-1) 
                              Robs = obs%obsvar(idobs)
!                              ypost_var    = 1.0_r8/(1.0_r8/yprior_var + 1.0_r8/Robs)

                              ! compute 'TPW' posterior
                              allocate(obs_var_1d(ens_size))
                              do indk1 = 1, ens_size
                                 allocate(cam_var_1d(cam%nlev))
                                 allocate(qv_1d(cam%nlev))
                                 ref_P  = cam%P0(1)
                                 surf_P = cam%variables_post(i,j,1,igs+1+indk1,1)
                                 call get_pressure(ref_P, surf_P, cam%nlev,                  &
                                                   cam_hyam, cam_hybm, cam_var_1d)
                                 qv_1d =  cam%variables_post(i,j,:,igs+1+indk1,5)
                                 call get_expected_tpw(cam_var_1d, qv_1d, cam%nlev,          &
                                                       obs_var_1d(indk1), pressure_top)
                                 deallocate(cam_var_1d, qv_1d)
                              enddo
                              yamean(ig) = sum(obs_var_1d)/ens_size
                              ypost_var(ig) = 0.0_r8
                              do indk1 = 1, ens_size
                                 ypost_var(ig) = ypost_var(ig) + (obs_var_1d(indk1)-yamean(ig))**2
                              enddo
                              ypost_var(ig) = ypost_var(ig)/(ens_size-1)
                              deallocate(obs_var_1d)

                              ! get [lon, lat] and convert to radiance
                              obslon = cam%var_lon(i,2) * DEG2RAD
                              obslat = cam%var_lat(j,2) * DEG2RAD
   
!print *, 'y info = ', ytrue, yfmean, yprior_var, Robs, ypost_var
!print *, 'yprior = ', yprior(1:ens_size:10)
!
!print *, 'y loc = ', cam%var_lon(i,idobs), cam%var_lat(j,idobs), obslon, obslat
!pause

                           else
!                              ytrue      = cam%variables(i,j,k,1,idobs)
!                              ypost_var    = 1.0_r8/(1.0_r8/yprior_var + 1.0_r8/Robs)
                              yfmean(ig)     = cam%variables(i,j,k,igs,idobs)
                              yprior_var(ig) = cam%variables(i,j,k,igs+1,idobs)**2
                              yprior(1:ens_size) = cam%variables(i,j,k,igs+1+1:igs+1+ens_size,idobs)
                              Robs           = obs%obsvar(idobs)
                              yamean(ig)     = cam%variables_post(i,j,k,igs,idobs)
                              ypost_var(ig)  = cam%variables_post(i,j,k,igs+1,idobs)**2

!print *, 'y info = ', ytrue, yfmean, yprior_var, Robs, ypost_var
!print *, 'yprior = ', yprior(1:ens_size:10)

                  ! get [lon, lat] and convert to radiance
                              obslon = cam%var_lon(i,idobs) * DEG2RAD
                              obslat = cam%var_lat(j,idobs) * DEG2RAD
                           endif

!print *, 'y loc = ', cam%var_lon(i,idobs), cam%var_lat(j,idobs), obslon, obslat
!pause


                           ! get target variable information ready
!                           xtrue      = cam%variables(ii,jj,kk,1,inddim)
                           xfmean(ig)     = cam%variables(ii,jj,kk,igs,inddim)
                           xprior_var(ig) = cam%variables(ii,jj,kk,igs+1,inddim)**2
                           xprior(1:ens_size) = cam%variables(ii,jj,kk,igs+1+1:igs+1+ens_size,inddim)
                           xamean(ig)     = cam%variables_post(ii,jj,kk,igs,inddim)
                           xpost_var(ig)  = cam%variables_post(ii,jj,kk,igs+1,inddim)**2

!print *, 'var info = ', xtrue, xfmean, xprior_var
!print *, 'var prior = ', xprior(1:ens_size:10)

                           ! get [lon, lat] and convert to radiance
                           varlon = cam%var_lon(ii,inddim) * DEG2RAD
                           varlat = cam%var_lat(jj,inddim) * DEG2RAD

!print *, 'var loc = ', cam%var_lon(ii,inddim), cam%var_lat(jj,inddim), varlon, varlat
!pause

                           ! get the horizontal distance
                           dist = get_horiz_dist(varlon, varlat, obslon, obslat)

!print *, 'dist = ', dist

                           call find_index_in_locind(ilocind, dist, locnum, localization%locind)

!print *, 'ilocind = ', ilocind
!pause

                           ! localization computation
!                           xt_xfmean = xtrue - xfmean
                           covxy(ig) = comp_cov(ens_size,xfmean(ig),xprior,yfmean(ig),yprior)
                           corr(ig) = covxy(ig) / (sqrt(yprior_var(ig)) * sqrt(xprior_var(ig)))
                           reg_coef(ig) = covxy(ig) / yprior_var(ig)

                        enddo   ! ig

                        ! record in distance
                        localization%num(ilocind,ireg,ndims,inddim) =                   &
                            localization%num(ilocind,ireg,ndims,inddim) + 1

                        ! correlation (this has no impact from inflation_correction, different ELF algorithm)
                        do ig = 1, group_size
                           localization%correlation(ilocind,ig,ireg,ndims,inddim) =                     &
                               localization%correlation(ilocind,ig,ireg,ndims,inddim) + abs(corr(ig))
                           localization%sumcorr(ilocind,ig,ireg,ndims,inddim) =                         &
                               localization%sumcorr(ilocind,ig,ireg,ndims,inddim) + corr(ig)
                           localization%sumcorrsq(ilocind,ig,ireg,ndims,inddim) =                       &
                               localization%sumcorrsq(ilocind,ig,ireg,ndims,inddim) + corr(ig)**2
                        enddo


                        !-----------------------------------------------!
                        ! GGF_reg and DGF_reg
                        !-----------------------------------------------!
                        sumbeta1 = (sum(reg_coef))**2
                        sumbeta2 = sum(reg_coef*reg_coef)
                        localization%beta1(ilocind,ireg,ndims,inddim) =                  &
                            localization%beta1(ilocind,ireg,ndims,inddim) + sumbeta1
                        localization%beta2(ilocind,ireg,ndims,inddim) =                  &
                            localization%beta2(ilocind,ireg,ndims,inddim) + sumbeta2
                        localization%meangf(ilocind,ireg,ndims,inddim) =                 &
                            localization%meangf(ilocind,ireg,ndims,inddim) + (sumbeta1/sumbeta2-1.0_r8)/(group_size-1.0_r8)
                        !-----------------------------------------------!


                        !-----------------------------------------------!
                        ! GGF_inc, DGF_inc, and ELF with cross observation (_cobs)
                        ! unknown the "true", use the xamean/yamean of another group as an estimate of "true",
                        ! also integrate over all possibilities, that is having the term with Robs
                        !-----------------------------------------------!
                        indk2 = 1
                        do igg = 1, group_size
                           ! indk1 is the index of the estimate "true"
                           indk1 = igg + 1
                           if ( indk1 == group_size + 1 ) then
                              indk1 = 1
                           endif
                           if ( indk1 < 1 .or. indk1 > group_size ) then
                              print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                           endif
                           ddd(igg) = covxy(igg)/(yprior_var(igg)+Robs)
                           sublocdenom(igg) = (ddd(igg)*(yamean(indk1)-yfmean(igg)))**2 + ddd(igg)**2 * Robs
                           if ( abs(sublocdenom(igg)) < 0.0000000001_r8 ) then
                               indk2 = -1
                           endif
                        enddo
                        if ( indk2 > 0 ) then
                           localization%num_cobs(ilocind,ireg,ndims,inddim) =                   &
                               localization%num_cobs(ilocind,ireg,ndims,inddim) + 1

                           do igg = 1, group_size
                              ! indk1 is the index of the estimate "true"
                              indk1 = igg + 1
                              if ( indk1 == group_size + 1 ) then
                                 indk1 = 1
                              endif
                              if ( indk1 < 1 .or. indk1 > group_size ) then
                                 print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                              endif
                              ddd(igg) = covxy(igg)/(yprior_var(igg)+Robs)
                              incnumer(igg) = ddd(igg)*(yamean(indk1)-yfmean(igg))
                           enddo
                           incsumbeta1 = (sum(incnumer))**2 + Robs * (sum(ddd))**2
                           incsumbeta2 = sum(incnumer*incnumer) + Robs * sum(ddd*ddd)
                           localization%incbeta1_cobs(ilocind,ireg,ndims,inddim) =             &
                               localization%incbeta1_cobs(ilocind,ireg,ndims,inddim) + incsumbeta1
                           localization%incbeta2_cobs(ilocind,ireg,ndims,inddim) =             &
                               localization%incbeta2_cobs(ilocind,ireg,ndims,inddim) + incsumbeta2
                           localization%incmeangf_cobs(ilocind,ireg,ndims,inddim) =            &
                               localization%incmeangf_cobs(ilocind,ireg,ndims,inddim) +        &
                               (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                              ! indk1 is the index of the estimate "true"
                              indk1 = ig + 1
                              if ( indk1 == group_size + 1 ) then
                                 indk1 = 1
                              endif
                              if ( indk1 < 1 .or. indk1 > group_size ) then
                                 print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                              endif
                              sublocx(ig) = xamean(indk1) - xfmean(ig)
                              ddd(ig)     = covxy(ig)/(yprior_var(ig)+Robs)
                              sublocy(ig) = ddd(ig)*(yamean(indk1)-yfmean(ig))
                              localization%sumx_cobs(ilocind,ig,ireg,ndims,inddim) =           &
                                  localization%sumx_cobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)
                              localization%sumy_cobs(ilocind,ig,ireg,ndims,inddim) =           &
                                  localization%sumy_cobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)
                              localization%sumx2_cobs(ilocind,ig,ireg,ndims,inddim) =          &
                                  localization%sumx2_cobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)**2
                              localization%sumy2_cobs(ilocind,ig,ireg,ndims,inddim) =          &
                                  localization%sumy2_cobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)**2
                              localization%numerator_cobs(ilocind,ig,ireg,ndims,inddim) =      &
                                  localization%numerator_cobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)*sublocy(ig)
                              localization%denominator_cobs(ilocind,ig,ireg,ndims,inddim) =    &
                                  localization%denominator_cobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)**2 + ddd(ig)**2 * Robs
                           enddo   ! ig

                        endif   ! indk2 > 0
                        !-----------------------------------------------!


                        !-----------------------------------------------!
                        ! GGF_inc, DGF_inc, and ELF with fake perfect observation (_pobs)
                        ! unknown the "true", use the xamean/yamean of another group as an estimate of "true",
                        ! Robs = 0
                        !-----------------------------------------------!
                        indk2 = 1
                        do igg = 1, group_size
                           ! indk1 is the index of the estimate "true"
                           indk1 = igg + 1
                           if ( indk1 == group_size + 1 ) then
                              indk1 = 1
                           endif
                           if ( indk1 < 1 .or. indk1 > group_size ) then
                              print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                           endif
                           ddd(igg) = covxy(igg)/yprior_var(igg)
                           sublocdenom(igg) = (ddd(igg)*(yamean(indk1)-yfmean(igg)))**2
                           if ( abs(sublocdenom(igg)) < 0.0000000001_r8 ) then
                               indk2 = -1
                           endif
                        enddo
                        if ( indk2 > 0 ) then
                           localization%num_pobs(ilocind,ireg,ndims,inddim) =                   &
                               localization%num_pobs(ilocind,ireg,ndims,inddim) + 1

                           do igg = 1, group_size
                              ! indk1 is the index of the estimate "true"
                              indk1 = igg + 1
                              if ( indk1 == group_size + 1 ) then
                                 indk1 = 1
                              endif
                              if ( indk1 < 1 .or. indk1 > group_size ) then
                                 print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                              endif
                              ddd(igg) = covxy(igg)/yprior_var(igg)
                              incnumer(igg) = ddd(igg)*(yamean(indk1)-yfmean(igg))
                           enddo
                           incsumbeta1 = (sum(incnumer))**2
                           incsumbeta2 = sum(incnumer*incnumer)
                           localization%incbeta1_pobs(ilocind,ireg,ndims,inddim) =             &
                               localization%incbeta1_pobs(ilocind,ireg,ndims,inddim) + incsumbeta1
                           localization%incbeta2_pobs(ilocind,ireg,ndims,inddim) =             &
                               localization%incbeta2_pobs(ilocind,ireg,ndims,inddim) + incsumbeta2
                           localization%incmeangf_pobs(ilocind,ireg,ndims,inddim) =            &
                               localization%incmeangf_pobs(ilocind,ireg,ndims,inddim) +        &
                               (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                              ! indk1 is the index of the estimate "true"
                              indk1 = ig + 1
                              if ( indk1 == group_size + 1 ) then
                                 indk1 = 1
                              endif
                              if ( indk1 < 1 .or. indk1 > group_size ) then
                                 print *, 'WRONG index for estimate of "true": ', igg, group_size, indk1
                              endif
                              sublocx(ig) = xamean(indk1) - xfmean(ig)
                              ddd(ig)     = covxy(ig)/yprior_var(ig)
                              sublocy(ig) = ddd(ig)*(yamean(indk1)-yfmean(ig))
                              localization%sumx_pobs(ilocind,ig,ireg,ndims,inddim) =           &
                                  localization%sumx_pobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)
                              localization%sumy_pobs(ilocind,ig,ireg,ndims,inddim) =           &
                                  localization%sumy_pobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)
                              localization%sumx2_pobs(ilocind,ig,ireg,ndims,inddim) =          &
                                  localization%sumx2_pobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)**2
                              localization%sumy2_pobs(ilocind,ig,ireg,ndims,inddim) =          &
                                  localization%sumy2_pobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)**2
                              localization%numerator_pobs(ilocind,ig,ireg,ndims,inddim) =      &
                                  localization%numerator_pobs(ilocind,ig,ireg,ndims,inddim) + sublocx(ig)*sublocy(ig)
                              localization%denominator_pobs(ilocind,ig,ireg,ndims,inddim) =    &
                                  localization%denominator_pobs(ilocind,ig,ireg,ndims,inddim) + sublocy(ig)**2
                           enddo   ! ig

                        endif   ! indk2 > 0
                        !-----------------------------------------------!

                     enddo  ! ii
                  enddo     ! jj
               enddo        ! i
            enddo           ! j

            enddo  ! ireg (subregion)

         enddo     ! varlevk
      enddo        ! obslevk

      deallocate(yprior,xprior)
      deallocate(xfmean,xamean,xprior_var,xpost_var)
      deallocate(yfmean,yamean,yprior_var,ypost_var)
      deallocate(covxy,corr,reg_coef)
      deallocate(ccc,ddd)
      deallocate(sublocx,sublocy,sublocx2,sublocy2)
      deallocate(sublocnumer,sublocdenom,sublocalpha)
      deallocate(incnumer,incdenom)

   enddo LoopVarTypes
!$OMP END PARALLEL DO
  
! NOTE: the writting in OMP does not work...
!       thus put the writting out of OMP loop
!   do inddim = 1, cam%number_of_state_variables
   do inddim = indvar_s, indvar_e
      writelenloc = locnum

      do obslevk = 1, obslev
         do varlevk = 1, varlev
            ndims = (obslevk-1)*varlev + varlevk

! 2012-08-06: add PS even PS is not at variable level
!             treat PS as integral variable
            if ( idobs == 1 ) then
               write(idchar1,'(i1)') 1
            else
               if ( obs_ind(obslevk) < 10 ) then
                  write(idchar1,'(i1)') obs_ind(obslevk)
               else
                  write(idchar1,'(i2)') obs_ind(obslevk)
               endif
            endif
! 2012-08-06: add PS even PS is not at variable level
!             treat PS as integral variable
            if ( inddim == 1 ) then
               write(idchar2,'(i1)') 1
            else
               if ( var_ind(obslevk,varlevk) < 10 ) then
                  write(idchar2,'(i1)') var_ind(obslevk,varlevk)
               else
                  write(idchar2,'(i2)') var_ind(obslevk,varlevk)
               endif
            endif
            var_name_tmp = trim(obs%obstype_metadata(idobs))//'_'//         &
                           trim(cam%description(inddim))//'_obslev'//trim(idchar1)//'_varlev'//trim(idchar2)
            iunit = 10 + inddim
            open(iunit, FILE=var_name_tmp, FORM='FORMATTED')

            do ireg = 1, dim_region

               allocate(numtmp(writelenloc))
               numtmp(1:writelenloc) = localization%num(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) numtmp

               do ig = 1, group_size
                  write(iunit, FMT=locreal_format) localization%correlation(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumcorr(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumcorrsq(1:writelenloc,ig,ireg,ndims,inddim)
               enddo

               write(iunit, FMT=locreal_format) localization%beta1(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%beta2(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%meangf(1:writelenloc,ireg,ndims,inddim)

               ! ELF of cross observation with Robs
               numtmp(1:writelenloc) = localization%num_cobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) numtmp

               write(iunit, FMT=locreal_format) localization%incbeta1_cobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%incbeta2_cobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%incmeangf_cobs(1:writelenloc,ireg,ndims,inddim)

               do ig = 1, group_size
                  write(iunit, FMT=locreal_format) localization%sumx_cobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy_cobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumx2_cobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy2_cobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%numerator_cobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%denominator_cobs(1:writelenloc,ig,ireg,ndims,inddim)
               enddo

               ! ELF of fake perfect observation (Robs = 0)
               numtmp(1:writelenloc) = localization%num_pobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) numtmp

               write(iunit, FMT=locreal_format) localization%incbeta1_pobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%incbeta2_pobs(1:writelenloc,ireg,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%incmeangf_pobs(1:writelenloc,ireg,ndims,inddim)

               do ig = 1, group_size
                  write(iunit, FMT=locreal_format) localization%sumx_pobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy_pobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumx2_pobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy2_pobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%numerator_pobs(1:writelenloc,ig,ireg,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%denominator_pobs(1:writelenloc,ig,ireg,ndims,inddim)
               enddo

               deallocate(numtmp)

            enddo   ! ireg

            close(iunit)

         enddo   ! varlevk
      enddo      ! obslevk
   enddo         ! inddim

   deallocate(localization%num)
   deallocate(localization%correlation,localization%sumcorr,localization%sumcorrsq)
   deallocate(localization%beta1,localization%beta2,localization%meangf)

   deallocate(localization%num_cobs)
   deallocate(localization%incbeta1_cobs,localization%incbeta2_cobs,localization%incmeangf_cobs)
   deallocate(localization%sumx_cobs,localization%sumy_cobs)
   deallocate(localization%sumx2_cobs,localization%sumy2_cobs)
   deallocate(localization%numerator_cobs,localization%denominator_cobs)

   deallocate(localization%num_pobs)
   deallocate(localization%incbeta1_pobs,localization%incbeta2_pobs,localization%incmeangf_pobs)
   deallocate(localization%sumx_pobs,localization%sumy_pobs)
   deallocate(localization%sumx2_pobs,localization%sumy2_pobs)
   deallocate(localization%numerator_pobs,localization%denominator_pobs)
 
enddo LoopObstypes

deallocate(dim_latlon)
deallocate(cam_hyam, cam_hybm)




contains








subroutine draw_perfect_obs_2d(we,sn,nbdyoff,ens_size,description,wrf_var_3d,obs_var_3d)

integer,            intent(in)  :: we, sn, nbdyoff, ens_size
character(len=129), intent(in)  :: description
real(r8),           intent(in)  :: wrf_var_3d(:,:,:)
real(r8),           intent(out) :: obs_var_3d(:,:,:) 
! real(r8),           intent(in)  :: Robs 

integer    :: i, j, k, ind, obscopysize
real(r8)   :: pres, oerr

obscopysize = ens_size + 4

obs_var_3d(1:we,1:sn,1:obscopysize) = missing_r8
   
if ( (trim(description) == 'SFC_U_WIND_COMPONENT') .or.            &
     (trim(description) == 'SFC_V_WIND_COMPONENT') ) then 

    obs_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,1:obscopysize-1) =    &
        wrf_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,1:obscopysize-1)

    pres = 1000.00_r8              ! this pressure is not really used in land_wind_error
    oerr = land_wind_error(pres)
    obs_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,obscopysize) = oerr * oerr 

elseif ( (trim(description) == 'SFC_TEMPERATURE') ) then

    obs_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,1:obscopysize-1) =    &
        wrf_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,1:obscopysize-1)

    pres = 1000.00_r8              ! this pressure is not really used in land_temp_error
    oerr = land_temp_error(pres)
    obs_var_3d(nbdyoff+1:we-nbdyoff,nbdyoff+1:sn-nbdyoff,obscopysize) = oerr * oerr

else
    print *, 'ERROR: draw_perfect_obs_2d'
    print *, 'ERROR: wrong obs type', trim(description)
endif

return
end subroutine draw_perfect_obs_2d


subroutine draw_perfect_obs_altimeter(we,sn,nbdyoff,ens_size,description,wrf_var_2d,wrf_var_3d,obs_var_3d)

integer,            intent(in)  :: we, sn, nbdyoff, ens_size
character(len=129), intent(in)  :: description
real(r8),           intent(in)  :: wrf_var_2d(:,:)      ! hgt
real(r8),           intent(in)  :: wrf_var_3d(:,:,:)    ! psfc
real(r8),           intent(out) :: obs_var_3d(:,:,:)
! real(r8),           intent(in)  :: Robs

integer    :: i, j, k, ind, obscopysize
real(r8)   :: pres, hsfc, altr, var1, oerr
real(r8), allocatable :: temp_var_1d(:)

allocate(temp_var_1d(ens_size))

obscopysize = ens_size + 4

obs_var_3d(1:we,1:sn,1:obscopysize) = missing_r8

if ( (trim(description) == "SFC_ALTIMETER") ) then 

    do j = nbdyoff+1, sn-nbdyoff
       do i = nbdyoff+1, we-nbdyoff
          ! compute the truth
          pres = wrf_var_3d(i,j,1)
          hsfc = wrf_var_2d(i,j)
          altr = compute_altimeter(pres*0.01_r8,hsfc)          
          if ( (altr >= 880.0_r8) .and. (altr <1100.0_r8) )  then
              obs_var_3d(i,j,1) = altr
          else
              return
          endif

          ! assign obs error variance
          oerr = land_pres_error(pres)
          obs_var_3d(i,j,obscopysize) = oerr * oerr

          ! compute ensemble members
          do ind = 1, ens_size
             pres = wrf_var_3d(i,j,ind+3)
             altr = compute_altimeter(pres*0.01_r8,hsfc)
             if ( (altr >= 880.0_r8) .and. (altr <1100.0_r8) )  then
                 obs_var_3d(i,j,ind+3) = altr
             endif
          enddo

          ! compute ensemble mean
          temp_var_1d(1:ens_size) = obs_var_3d(i,j,4:ens_size+3)
          obs_var_3d(i,j,2) = sum(temp_var_1d) / ens_size

          ! compute ensemble spread
          var1 = 0.0_r8
          do k = 1, ens_size
             var1 = var1 + (temp_var_1d(k)-obs_var_3d(i,j,2))**2
          enddo
          var1 = var1 / (ens_size - 1)
          obs_var_3d(i,j,3) = sqrt(var1)

       enddo
    enddo

else
    print *, 'ERROR: draw_perfect_obs_altimeter'
    print *, 'ERROR: wrong obs type', trim(description)
endif

deallocate(temp_var_1d)

return
end subroutine draw_perfect_obs_altimeter


subroutine draw_perfect_obs_2d_dewpoint(we,sn,nbdyoff,ens_size,description,pres_var_3d,qv_var_3d,wrf_var_2d,obs_var_3d)

integer,            intent(in)  :: we, sn, nbdyoff, ens_size
character(len=129), intent(in)  :: description
real(r8),           intent(in)  :: pres_var_3d(:,:,:)  ! Psfc     
real(r8),           intent(in)  :: qv_var_3d(:,:,:)    ! Q2
real(r8),           intent(in)  :: wrf_var_2d(:,:)     ! T2
real(r8),           intent(out) :: obs_var_3d(:,:,:)
! real(r8),           intent(in)  :: Robs

integer    :: i, j, k, ind, obscopysize
real(r8)   :: pres, qv, td, tair, rh, var1, qerr, oerr
real(r8), allocatable :: temp_var_1d(:)

allocate(temp_var_1d(ens_size))

obscopysize = ens_size + 4

obs_var_3d(1:we,1:sn,1:obscopysize) = missing_r8

if ( (trim(description) == "SFC_DEWPOINT") ) then

    do j = nbdyoff+1, sn-nbdyoff
       do i = nbdyoff+1, we-nbdyoff
          ! compute the truth
          pres = pres_var_3d(i,j,1)
          qv   = qv_var_3d(i,j,1)
          if ( qv >= 0.0_r8 .and. qv < 1.0_r8) then
             td = compute_dewpoint(pres, qv) 
             obs_var_3d(i,j,1) = td
          endif

          ! assign obs error variance
          tair = wrf_var_2d(i,j)
          rh = temp_and_dewpoint_to_rh(tair, td)
          oerr = dewpt_error_from_rh_and_temp(tair,rh)
          obs_var_3d(i,j,obscopysize) = oerr * oerr

          ! compute ensemble members
          do ind = 1, ens_size
             pres = pres_var_3d(i,j,ind+3)
             qv   = qv_var_3d(i,j,ind+3)
             if ( qv >= 0.0_r8 .and. qv < 1.0_r8) then
                td = compute_dewpoint(pres, qv)
                obs_var_3d(i,j,ind+3) = td
             endif
          enddo

          ! compute ensemble mean
          temp_var_1d(1:ens_size) = obs_var_3d(i,j,4:ens_size+3)
          obs_var_3d(i,j,2) = sum(temp_var_1d) / ens_size

          ! compute ensemble spread
          var1 = 0.0_r8
          do k = 1, ens_size
             var1 = var1 + (temp_var_1d(k)-obs_var_3d(i,j,2))**2
          enddo
          var1 = var1 / (ens_size - 1)
          obs_var_3d(i,j,3) = sqrt(var1)

       enddo
    enddo

else
    print *, 'ERROR: draw_perfect_obs_2d_dewpoint'
    print *, 'ERROR: wrong obs type', trim(description)
endif

deallocate(temp_var_1d)

return
end subroutine draw_perfect_obs_2d_dewpoint


subroutine draw_perfect_obs_2d_q(we,sn,nbdyoff,ens_size,description,pres_var_3d,qv_var_3d,wrf_var_2d,obs_var_3d)

integer,            intent(in)  :: we, sn, nbdyoff, ens_size
character(len=129), intent(in)  :: description
real(r8),           intent(in)  :: pres_var_3d(:,:,:)  ! Psfc     
real(r8),           intent(in)  :: qv_var_3d(:,:,:)    ! Q2
real(r8),           intent(in)  :: wrf_var_2d(:,:)     ! T2
real(r8),           intent(out) :: obs_var_3d(:,:,:)
! real(r8),           intent(in)  :: Robs

integer    :: i, j, k, ind, obscopysize
real(r8)   :: pres, qv, qsat, td, tair, rh, var1, qerr, oerr
real(r8), allocatable :: temp_var_1d(:)

allocate(temp_var_1d(ens_size))

obscopysize = ens_size + 4

obs_var_3d(1:we,1:sn,1:obscopysize) = missing_r8

if ( (trim(description) == "SFC_SPECIFIC_HUMIDITY") ) then

    do j = nbdyoff+1, sn-nbdyoff
       do i = nbdyoff+1, we-nbdyoff
          ! compute the truth
          qv   = qv_var_3d(i,j,1)
          if ( qv >= 0.0_r8 .and. qv < 1.0_r8) then
             obs_var_3d(i,j,1) = qv
          endif

          ! assign obs error variance
          tair = wrf_var_2d(i,j)
          pres = pres_var_3d(i,j,1)
          qv   = qv_var_3d(i,j,1)
          if ( qv >= 0.0_r8 .and. qv < 1.0_r8) then
             td = compute_dewpoint(pres, qv)
             qsat = specific_humidity(sat_vapor_pressure(tair), pres * 100.0_r8)
             if ( tair >= td ) then
                qerr = rh_error_from_dewpt_and_temp(tair, td)
                oerr = max(qerr * qsat, 0.0001_r8)
                obs_var_3d(i,j,obscopysize) = oerr * oerr
             else
                obs_var_3d(i,j,1) = missing_r8
             endif
          endif

          ! compute ensemble members
          do ind = 1, ens_size
             pres = pres_var_3d(i,j,ind+3)
             qv   = qv_var_3d(i,j,ind+3)
             if ( qv >= 0.0_r8 .and. qv < 1.0_r8) then
                obs_var_3d(i,j,ind+3) = qv 
             endif
          enddo

          ! compute ensemble mean
          temp_var_1d(1:ens_size) = obs_var_3d(i,j,4:ens_size+3)
          obs_var_3d(i,j,2) = sum(temp_var_1d) / ens_size

          ! compute ensemble spread
          var1 = 0.0_r8
          do k = 1, ens_size
             var1 = var1 + (temp_var_1d(k)-obs_var_3d(i,j,2))**2
          enddo
          var1 = var1 / (ens_size - 1)
          obs_var_3d(i,j,3) = sqrt(var1)

       enddo
    enddo
    
else
    print *, 'ERROR: draw_perfect_obs_2d_q'
    print *, 'ERROR: wrong obs type', trim(description)
endif

deallocate(temp_var_1d)

return
end subroutine draw_perfect_obs_2d_q








function compute_dewpoint(pres,qv)
! compute dewpoint
! adapted from obs_def_mod.f90/get_expected_dew_point 

real(r8),   intent(in)   :: pres, qv 
real(r8),   PARAMETER    :: e_min = 0.001_r8    ! threshold for minimum vapor pressure (mb),
                                                !   to avoid problems near zero in Bolton's equation
real(r8)                 :: compute_dewpoint 

real(r8)      :: p_mb, e_mb
real(r8)      :: ts0, cpovcv, rd_over_rv, qvf1, ph_e, rho

compute_dewpoint = missing_r8

p_mb = pres * 0.01_r8

e_mb = qv * p_mb / (0.622_r8 + qv)
e_mb = max(e_mb, e_min)

compute_dewpoint = t_kelvin + (243.5_r8 / ((17.67_r8 / log(e_mb/6.112_r8)) - 1.0_r8) )

end function compute_dewpoint 












subroutine count_numobstype(obssize, obstypein, obsqcin, nobstype, obstypeindexout)
! read in one obs_epoch file, count the number of obs type

integer, intent(in)         :: obssize
integer, intent(in)         :: obstypein(:)
integer, intent(in)         :: obsqcin(:,:)
integer, intent(out)        :: nobstype
integer, intent(out)        :: obstypeindexout(:)

integer  :: ind, indtmp, numobs, obsindextmp(100)

! initialize obstypeindex
obstypeindexout(1:100) = 0

nobstype = 0
do numobs = 1, obssize
   if ( (nobstype == 0) .and. (obsqcin(1,numobs) < 1.0) .and. (obsqcin(2,numobs) < 1.0) ) then
       nobstype = 1
       ind = numobs
       obsindextmp(1) = obstypein(numobs)
       goto 101
   endif
enddo

101 continue

if ( ind < obssize ) then
   do  numobs = ind+1, obssize
       do  indtmp = 1, nobstype
           if ( obstypein(numobs) == obsindextmp(indtmp) ) then
               goto 102
           endif
       enddo

       if ( (obsqcin(1,numobs) < 1.0) .and. (obsqcin(2,numobs) < 1.0) ) then
          nobstype = nobstype + 1
          obsindextmp(nobstype) = obstypein(numobs)
       endif

102 continue

   enddo
endif

obstypeindexout(1:nobstype) = obsindextmp(1:nobstype)

return
end subroutine count_numobstype


subroutine compute_model_height(we, wes, sn, sns, bt, bts, hgt, phb, ph, lon, lat, nmhgt, mhgt)
! compute the model height for each model height type
  
integer, intent(in)          :: we, wes, sn, sns, bt, bts
real(r8), intent(in)         :: hgt(:, :)
real(r8), intent(in)         :: phb(:, :, :), ph(:, :, :)
real(r8), intent(in)         :: lon(:, :), lat(:, :)
integer, intent(in)          :: nmhgt
real(r8), intent(out)        :: mhgt(:, :, :, :)

integer      :: i,j,k,ind
real(r8)     :: geop, lattmp

! mhgtindex: 1 - surface mass
!            2 - surface U
!            3 - surface T
!!            4 - 3d U
!!            5 - 3d V
!!            6 - 3d W (vertical half level)
!!            7 - 3d T (vertical full level, mass)
!            4 - 3d T

do ind = 1, nmhgt
   if ( ind == 1 ) then
        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn)
   elseif ( ind == 2 ) then
        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn) + 10.0
   elseif ( ind == 3 ) then
        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn) + 2.0
   elseif ( ind == 4 ) then
        do i = 1, we
           do j = 1, sn
              do k = 1, bt
                 geop = ( (phb(i,j,k)+ph(i,j,k))  + (phb(i,j,k+1)+ph(i,j,k+1)) ) / (2.0*gravity)
                 mhgt(i,j,k,ind) = compute_geometric_height(geop, lat(i,j))
              enddo
           enddo
        enddo 
   else
        print *, 'ERROR in compute_model_height'
        print *, 'ERROR: do not support model height type', ind
   endif
enddo

return
end subroutine compute_model_height

subroutine get_obs_model_pressure_profile(bt, vp, obsi, obsj, dnw, mub, mu, phb, ph, t, qvapor, psfc)
! compute the model_pressure_profile at obs point

integer,     intent(in)     :: bt 
real(r8),    intent(out)    :: vp(0:bt)
real(r8),    intent(in)     :: obsi, obsj
real(r8),    intent(in)     :: dnw(:)
real(r8),    intent(in)     :: mub(:,:), mu(:,:)
real(r8),    intent(in)     :: phb(:,:,:), ph(:,:,:)
real(r8),    intent(in)     :: t(:,:,:), qvapor(:,:,:)
real(r8),    intent(in)     :: psfc(:,:)

integer      :: ind, indi, indj
real(r8)     :: disti, disti2, distj, distj2
real(r8)     :: pres1, pres2, pres3, pres4

vp(0:bt) = missing_r8

call toGrid(obsi, indi, disti, disti2)
call toGrid(obsj, indj, distj, distj2)

! put PSFC to vp(0)
pres1 = psfc(indi,   indj)
pres2 = psfc(indi+1, indj)
pres3 = psfc(indi,   indj+1)
pres4 = psfc(indi+1, indj+1)
vp(0) = distj2*(disti2*pres1 + disti*pres2) + distj*(disti2*pres3 + disti*pres4)

! compute vp(1:bt)
do ind = 1, bt
   pres1 = model_pressure_t(indi,   indj,   ind, dnw, mub, mu, phb, ph, t, qvapor) 
   pres2 = model_pressure_t(indi+1, indj,   ind, dnw, mub, mu, phb, ph, t, qvapor)
   pres3 = model_pressure_t(indi,   indj+1, ind, dnw, mub, mu, phb, ph, t, qvapor)
   pres4 = model_pressure_t(indi+1, indj+1, ind, dnw, mub, mu, phb, ph, t, qvapor)

   vp(ind) = distj2*(disti2*pres1 + disti*pres2) + distj*(disti2*pres3 + disti*pres4)

enddo

if ( vp(0) < vp(1) ) then
   print *, 'ERROR in get_obs_model_pressure_profile'
   print *, 'ERROR: vp(0) < vp(1)', vp(0), vp(1)
endif

return
end subroutine get_obs_model_pressure_profile


function model_pressure_t(i,j,k,dnw,mub,mu,phb,ph,t,qvapor)
! compute model_pressure give (i,j,k) on mass point full level (mass level)

integer,    intent(in)   :: i, j, k
real(r8),   intent(in)   :: dnw(:)
real(r8),   intent(in)   :: mub(:,:), mu(:,:)
real(r8),   intent(in)   :: phb(:,:,:), ph(:,:,:) 
real(r8),   intent(in)   :: t(:,:,:), qvapor(:,:,:)
real(r8)                 :: model_pressure_t

real(r8)      :: ts0, cpovcv, rd_over_rv, qvf1, ph_e, rho

model_pressure_t = missing_r8

ts0 = 300.0
cpovcv = 1.4
rd_over_rv = gas_constant / gas_constant_v

qvf1 = 1.0 + qvapor(i,j,k) / rd_over_rv

ph_e = ( (phb(i,j,k+1) + ph(i,j,k+1)) - (phb(i,j,k)+ph(i,j,k)) ) / dnw(k) 

rho  = -( mub(i,j) + mu(i,j) ) / ph_e

model_pressure_t = ps0 * ( (gas_constant * (ts0 + t(i,j,k)) * qvf1) / (ps0/rho) )**cpovcv 

end function model_pressure_t


function grid_mhgt_from_vp(obsi, obsj, obsk, mhgt)
! get the model height of the observation

real(r8),    intent(in)   :: obsi, obsj, obsk
real(r8),    intent(in)   :: mhgt(:,:,:)
real(r8)                  :: grid_mhgt_from_vp

integer      :: ind, indi, indj, indk
real(r8)     :: disti, disti2, distj, distj2, distk, distk2
real(r8)     :: hgt1, hgt2, hgt3, hgt4, hgt_lev1, hgt_lev2

call toGrid(obsi, indi, disti, disti2)
call toGrid(obsj, indj, distj, distj2)
call toGrid(obsk, indk, distk, distk2)

hgt1 = mhgt(indi,   indj,   indk)
hgt2 = mhgt(indi+1, indj,   indk)
hgt3 = mhgt(indi,   indj+1, indk)
hgt4 = mhgt(indi+1, indj+1, indk)

hgt_lev1 = distj2*(disti2*hgt1 + disti*hgt2) + distj*(disti2*hgt3 + disti*hgt4)

hgt1 = mhgt(indi,   indj,   indk+1)
hgt2 = mhgt(indi+1, indj,   indk+1)
hgt3 = mhgt(indi,   indj+1, indk+1)
hgt4 = mhgt(indi+1, indj+1, indk+1)

hgt_lev2 = distj2*(disti2*hgt1 + disti*hgt2) + distj*(disti2*hgt3 + disti*hgt4)

grid_mhgt_from_vp = distk2*hgt_lev1 + distk*hgt_lev1


end function grid_mhgt_from_vp


subroutine get_model_variable_at_ijk(i,j,k,variable,description,nens,xtrue,xfmean,xprior)
! some state variables are at mass points
! some are not, these need interpolation

integer,            intent(in)  :: i, j, k
real(r8),           intent(in)  :: variable(:,:,:,:)
character(len=129), intent(in)  :: description
integer,            intent(in)  :: nens
real(r8),           intent(out) :: xtrue, xfmean
real(r8),           intent(out) :: xprior(:)

integer    :: ind

if ( (trim(description) == "U") ) then
   ! U is west-east stag
   xtrue = (variable(i,j,k,1) + variable(i+1,j,k,1)) / 2.0
   xfmean = (variable(i,j,k,2) + variable(i+1,j,k,2)) / 2.0
   do ind = 1, nens
      xprior(ind) = (variable(i,j,k,3+ind) + variable(i+1,j,k,3+ind)) / 2.0
   enddo
elseif ( (trim(description) == "V") ) then
   ! V is south-north stag
   xtrue = (variable(i,j,k,1) + variable(i,j+1,k,1)) / 2.0
   xfmean = (variable(i,j,k,2) + variable(i,j+1,k,2)) / 2.0
   do ind = 1, nens
      xprior(ind) = (variable(i,j,k,3+ind) + variable(i,j+1,k,3+ind)) / 2.0
   enddo
elseif ( (trim(description) == "W")     .or.       &
         (trim(description) == "PH") ) then
   ! W and PH are bottom-top stag
   xtrue = (variable(i,j,k,1) + variable(i,j,k+1,1)) / 2.0
   xfmean = (variable(i,j,k,2) + variable(i,j,k+1,2)) / 2.0
   do ind = 1, nens
      xprior(ind) = (variable(i,j,k,3+ind) + variable(i,j,k+1,3+ind)) / 2.0
   enddo
else
   ! other variables are on mass point
   xtrue = variable(i,j,k,1)
   xfmean = variable(i,j,k,2)
   xprior(1:nens) = variable(i,j,k,4:nens+3)
endif

return
end subroutine get_model_variable_at_ijk


subroutine get_model_variable(k1,k2,we,sn,bt,indk,ens_size,description,variable,xtrue3d,xfmean3d,xfspd3d,xprior3d)
! some state variables are at mass points
! some are not, these need interpolation
   
integer,            intent(in)  :: k1, k2, we, sn, bt, indk, ens_size
character(len=129), intent(in)  :: description
real(r8),           intent(in)  :: variable(:,:,:,:)
real(r8),           intent(out) :: xtrue3d(:,:,:), xfmean3d(:,:,:), xfspd3d(:,:,:)
real(r8),           intent(out) :: xprior3d(:,:,:,:)

integer    :: i, j, k, ind
real(r8)   :: var1, var2, var3, var4, varspd

if ( (trim(description) == "U10")      .or.      &
     (trim(description) == "V10")      .or.      &
     (trim(description) == "T2")       .or.      &
     (trim(description) == "Q2")       .or.      &
     (trim(description) == "TH2")      .or.      &
     (trim(description) == "MU")       .or.      &
     (trim(description) == "PSFC") ) then

   if ( (k1 /= 1) .or. (k2 /=1 ) .or. (indk > 1) ) then
      print *, 'ERROR: get_model_variable'
      print *, 'ERROR: layer higher than sfc', k2, indk
      return
   else

      xtrue3d(1:we,1:sn,1)  = variable(1:we,1:sn,1,1)
      xfmean3d(1:we,1:sn,1) = variable(1:we,1:sn,1,2)
      xfspd3d(1:we,1:sn,1)  = variable(1:we,1:sn,1,3)
      xprior3d(1:we,1:sn,1,1:ens_size) = variable(1:we,1:sn,1,4:ens_size+3)

   endif

elseif ( trim(description) == "U" )  then

   do k = k1, k2
      do j = 1, sn
         do i = 1, we

            xtrue3d(i,j,k)  = (variable(i,j,k,1) + variable(i+1,j,k,1)) / 2.0_r8
            xfmean3d(i,j,k) = (variable(i,j,k,2) + variable(i+1,j,k,2)) / 2.0_r8
            do ind = 1, ens_size            
               xprior3d(i,j,k,ind) = (variable(i,j,k,3+ind) + variable(i+1,j,k,3+ind)) / 2.0_r8
            enddo
            ! compute ensemble spread
            varspd = 0.0_r8
            do ind = 1, ens_size
              varspd = varspd + (xprior3d(i,j,k,ind)-xfmean3d(i,j,k))**2
            enddo
            varspd = varspd / (ens_size - 1)
            xfspd3d(i,j,k) = sqrt(varspd)

         enddo 
      enddo
   enddo

elseif ( trim(description) == "V" )  then

   do k = k1, k2
      do j = 1, sn
         do i = 1, we

            xtrue3d(i,j,k)  = (variable(i,j,k,1) + variable(i,j+1,k,1)) / 2.0_r8
            xfmean3d(i,j,k) = (variable(i,j,k,2) + variable(i,j+1,k,2)) / 2.0_r8
            do ind = 1, ens_size
               xprior3d(i,j,k,ind) = (variable(i,j,k,3+ind) + variable(i,j+1,k,3+ind)) / 2.0_r8
            enddo
            ! compute ensemble spread
            varspd = 0.0_r8
            do ind = 1, ens_size
              varspd = varspd + (xprior3d(i,j,k,ind)-xfmean3d(i,j,k))**2
            enddo
            varspd = varspd / (ens_size - 1)
            xfspd3d(i,j,k) = sqrt(varspd)

         enddo
      enddo
   enddo

elseif ( (trim(description) == "W")       .or.        &
         (trim(description) == "PH") ) then

   do k = k1, k2
      do j = 1, sn
         do i = 1, we

            xtrue3d(i,j,k)  = (variable(i,j,k,1) + variable(i,j,k+1,1)) / 2.0_r8
            xfmean3d(i,j,k) = (variable(i,j,k,2) + variable(i,j,k+1,2)) / 2.0_r8
            do ind = 1, ens_size
               xprior3d(i,j,k,ind) = (variable(i,j,k,3+ind) + variable(i,j,k+1,3+ind)) / 2.0_r8
            enddo
            ! compute ensemble spread
            varspd = 0.0_r8
            do ind = 1, ens_size
              varspd = varspd + (xprior3d(i,j,k,ind)-xfmean3d(i,j,k))**2
            enddo
            varspd = varspd / (ens_size - 1)
            xfspd3d(i,j,k) = sqrt(varspd)

         enddo
      enddo
   enddo

else

   xtrue3d(1:we,1:sn,k1:k2)  = variable(1:we,1:sn,k1:k2,1)
   xfmean3d(1:we,1:sn,k1:k2) = variable(1:we,1:sn,k1:k2,2)
   xfspd3d(1:we,1:sn,k1:k2)  = variable(1:we,1:sn,k1:k2,3)
   xprior3d(1:we,1:sn,k1:k2,1:ens_size) = variable(1:we,1:sn,k1:k2,4:ens_size+3)

endif

return
end subroutine get_model_variable


subroutine get_model_variable_profile(obsi,obsj,bt,variable,description,nens,xtrueprf,xfmeanprf,xfspdprf,xpriorprf)
! some state variables are at mass points
! some are not, these need interpolation

! given observation (i, j), find the model variable profile at (obsi, obsj)

real(r8),           intent(in)  :: obsi, obsj
integer,            intent(in)  :: bt
real(r8),           intent(in)  :: variable(:,:,:,:)
character(len=129), intent(in)  :: description
integer,            intent(in)  :: nens
real(r8),           intent(out) :: xtrueprf(:), xfmeanprf(:), xfspdprf(:)
real(r8),           intent(out) :: xpriorprf(:,:)

integer    :: indi, indj, k, ind
real(r8)   :: obsinew, obsjnew, disti, disti2, distj, distj2
real(r8)   :: var1, var2, var3, var4, var_k1, var_k2, varspd

obsinew = obsi
obsjnew = obsj

if ( (trim(description) == "U") ) then
   ! U is west-east stag
   obsinew = obsi + 0.5
   call toGrid(obsinew, indi, disti, disti2)
   call toGrid(obsjnew, indj, distj, distj2)
   do k = 1, bt
      ! copy truth
      var1 = variable(indi,  indj,  k,1)
      var2 = variable(indi+1,indj,  k,1)
      var3 = variable(indi,  indj+1,k,1)
      var4 = variable(indi+1,indj+1,k,1)
      xtrueprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy prior mean
      var1 = variable(indi,  indj,  k,2)
      var2 = variable(indi+1,indj,  k,2)
      var3 = variable(indi,  indj+1,k,2)
      var4 = variable(indi+1,indj+1,k,2)
      xfmeanprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy ensemble members
      do ind = 1, nens
         var1 = variable(indi,  indj,  k,3+ind)
         var2 = variable(indi+1,indj,  k,3+ind)
         var3 = variable(indi,  indj+1,k,3+ind)
         var4 = variable(indi+1,indj+1,k,3+ind)
         xpriorprf(k,ind) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      enddo
      ! compute ensemble spread
      varspd = 0.0_r8
      do ind = 1, nens
         varspd = varspd + (xpriorprf(k,ind)-xfmeanprf(k))**2
      enddo
      varspd = varspd / (nens - 1)
      xfspdprf(k) = sqrt(varspd)
   enddo
elseif ( (trim(description) == "V") ) then
   ! V is south-north stag
   obsjnew = obsj + 0.5
   call toGrid(obsinew, indi, disti, disti2)
   call toGrid(obsjnew, indj, distj, distj2)
   do k = 1, bt
      ! copy truth
      var1 = variable(indi,  indj,  k,1)
      var2 = variable(indi+1,indj,  k,1)
      var3 = variable(indi,  indj+1,k,1)
      var4 = variable(indi+1,indj+1,k,1)
      xtrueprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy prior mean
      var1 = variable(indi,  indj,  k,2)
      var2 = variable(indi+1,indj,  k,2)
      var3 = variable(indi,  indj+1,k,2)
      var4 = variable(indi+1,indj+1,k,2)
      xfmeanprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy ensemble members
      do ind = 1, nens
         var1 = variable(indi,  indj,  k,3+ind)
         var2 = variable(indi+1,indj,  k,3+ind)
         var3 = variable(indi,  indj+1,k,3+ind)
         var4 = variable(indi+1,indj+1,k,3+ind)
         xpriorprf(k,ind) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      enddo
      ! compute ensemble spread
      varspd = 0.0_r8
      do ind = 1, nens
         varspd = varspd + (xpriorprf(k,ind)-xfmeanprf(k))**2
      enddo
      varspd = varspd / (nens - 1)
      xfspdprf(k) = sqrt(varspd)
   enddo
elseif ( (trim(description) == "W")     .or.       &
         (trim(description) == "PH") ) then
   ! W and PH are bottom-top stag
   call toGrid(obsinew, indi, disti, disti2)
   call toGrid(obsjnew, indj, distj, distj2)
   do k = 1, bt
      ! copy truth
      var1 = variable(indi,  indj,  k,1)
      var2 = variable(indi+1,indj,  k,1)
      var3 = variable(indi,  indj+1,k,1)
      var4 = variable(indi+1,indj+1,k,1)
      var_k1 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      var1 = variable(indi,  indj,  k+1,1)
      var2 = variable(indi+1,indj,  k+1,1)
      var3 = variable(indi,  indj+1,k+1,1)
      var4 = variable(indi+1,indj+1,k+1,1)
      var_k2 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      xtrueprf(k) = (var_k1 + var_k2) / 2.0 
      ! copy prior mean
      var1 = variable(indi,  indj,  k,2)
      var2 = variable(indi+1,indj,  k,2)
      var3 = variable(indi,  indj+1,k,2)
      var4 = variable(indi+1,indj+1,k,2)
      var_k1 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      var1 = variable(indi,  indj,  k+1,2)
      var2 = variable(indi+1,indj,  k+1,2)
      var3 = variable(indi,  indj+1,k+1,2)
      var4 = variable(indi+1,indj+1,k+1,2)
      var_k2 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      xfmeanprf(k) = (var_k1 + var_k2) / 2.0
      ! copy ensemble members
      do ind = 1, nens
         var1 = variable(indi,  indj,  k,3+ind)
         var2 = variable(indi+1,indj,  k,3+ind)
         var3 = variable(indi,  indj+1,k,3+ind)
         var4 = variable(indi+1,indj+1,k,3+ind)
         var_k1 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
         var1 = variable(indi,  indj,  k+1,3+ind)
         var2 = variable(indi+1,indj,  k+1,3+ind)
         var3 = variable(indi,  indj+1,k+1,3+ind)
         var4 = variable(indi+1,indj+1,k+1,3+ind)
         var_k2 = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
         xpriorprf(k,ind) = (var_k1 + var_k2) / 2.0 
      enddo
      ! compute ensemble spread
      varspd = 0.0_r8
      do ind = 1, nens
         varspd = varspd + (xpriorprf(k,ind)-xfmeanprf(k))**2
      enddo
      varspd = varspd / (nens - 1)
      xfspdprf(k) = sqrt(varspd)
   enddo
else
   ! other variables are on mass point
   call toGrid(obsinew, indi, disti, disti2)
   call toGrid(obsjnew, indj, distj, distj2)
   do k = 1, bt
   ! copy truth
      var1 = variable(indi,  indj,  k,1)
      var2 = variable(indi+1,indj,  k,1)
      var3 = variable(indi,  indj+1,k,1)
      var4 = variable(indi+1,indj+1,k,1)
      xtrueprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy prior mean
      var1 = variable(indi,  indj,  k,2)
      var2 = variable(indi+1,indj,  k,2)
      var3 = variable(indi,  indj+1,k,2)
      var4 = variable(indi+1,indj+1,k,2)
      xfmeanprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      ! copy ensemble members
      do ind = 1, nens
         var1 = variable(indi,  indj,  k,3+ind)
         var2 = variable(indi+1,indj,  k,3+ind)
         var3 = variable(indi,  indj+1,k,3+ind)
         var4 = variable(indi+1,indj+1,k,3+ind)
         xpriorprf(k,ind) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
      enddo
      ! compute ensemble spread
      varspd = 0.0_r8
      do ind = 1, nens
         varspd = varspd + (xpriorprf(k,ind)-xfmeanprf(k))**2
      enddo
      varspd = varspd / (nens - 1)
      xfspdprf(k) = sqrt(varspd)
   enddo
endif

return
end subroutine get_model_variable_profile


subroutine get_model_height_profile(obsi,obsj,bt,mhgt,hgtprf)
! given observation (i, j), find the model height profile at (obsi, obsj)

real(r8),           intent(in)  :: obsi, obsj
integer,            intent(in)  :: bt
real(r8),           intent(in)  :: mhgt(:,:,:)
real(r8),           intent(out) :: hgtprf(:)

integer    :: indi, indj, k, ind
real(r8)   :: obsinew, obsjnew, disti, disti2, distj, distj2
real(r8)   :: var1, var2, var3, var4, var_k1, var_k2

obsinew = obsi
obsjnew = obsj

! 2012-03-07, 3D model height only have one type --- on T (mass, half level)
call toGrid(obsinew, indi, disti, disti2)
call toGrid(obsjnew, indj, distj, distj2)
do k = 1, bt
   var1 = mhgt(indi,  indj,  k)
   var2 = mhgt(indi+1,indj,  k)
   var3 = mhgt(indi,  indj+1,k)
   var4 = mhgt(indi+1,indj+1,k)
   hgtprf(k) = distj2*(disti2*var1+disti*var2) + distj*(disti2*var3+disti*var4)
enddo

return
end subroutine get_model_height_profile


subroutine comp_emp_loc_horiz(ilocind, sublocnum, sublocnomi, sublocdenomi,      &
                              locnum, locind,                                    &
                              obslongitude, obslatitude,                         &
                              varlongitude, varlatitude,                         &
                              xt_xfmean, deltaxmean)

!                                truth, state, hgt, mhgt, cutoff)
! horizontal empirical localization for 2D variable

integer,    intent(out)  :: ilocind, sublocnum                 
real(r8),   intent(out)  :: sublocnomi, sublocdenomi
integer,    intent(in)   :: locnum
real(r8),   intent(in)   :: locind(:)
real(r8),   intent(in)   :: obslongitude, obslatitude
real(r8),   intent(in)   :: varlongitude, varlatitude
real(r8),   intent(in)   :: xt_xfmean, deltaxmean

real(r8)     :: obslon, obslat, varlon, varlat, dist

! initialize output variables
ilocind      = -1 
sublocnum    = 0.0
sublocnomi   = 0.0
sublocdenomi = 0.0

! convert [lon, lat] to radiance
obslon = obslongitude * DEG2RAD
obslat = obslatitude  * DEG2RAD
varlon = varlongitude * DEG2RAD
varlat = varlatitude  * DEG2RAD
! get the horizontal distance
dist = get_horiz_dist(varlon, varlat, obslon, obslat)

!! Lili test
!print *, 'find dadada ', dist, locnum, locind(1:2000:200)

! find the index in localization%locind
call find_index_in_locind(ilocind, dist, locnum, locind)
if ( ilocind == -1 ) then
   print *, 'ERROR in subroutine comp_emp_loc_horiz'
   print *, 'ERROR: no index found in localization%obsloc', dist
   return
else
   sublocnum = 1
   sublocnomi = xt_xfmean * deltaxmean
   sublocdenomi = deltaxmean**2
endif

return
end subroutine comp_emp_loc_horiz


subroutine comp_emp_loc_vert(ilocind, sublocnum, sublocnomi, sublocdenomi,      &
                             locnumvert, locind,                                &
                             obslatitude, obshgt, varhgt,                       &
                             xt_xfmean, deltaxmean)

!                                truth, state, hgt, mhgt, cutoff)
! horizontal empirical localization for 2D variable

integer,    intent(out)  :: ilocind, sublocnum
real(r8),   intent(out)  :: sublocnomi, sublocdenomi
integer,    intent(in)   :: locnumvert
real(r8),   intent(in)   :: locind(:)
real(r8),   intent(in)   :: obslatitude, obshgt, varhgt
real(r8),   intent(in)   :: xt_xfmean, deltaxmean

real(r8)     :: earth_radius = 6370000.0       ! unit: meter 
real(r8)     :: obslon, obslat, varlon, varlat, dist

! initialize output variables
ilocind      = -1
sublocnum    = 0.0
sublocnomi   = 0.0
sublocdenomi = 0.0

! distance in radiance
! drad/2pi = dx/(2pi*R) => drad = dx/R
! R changes with latitude
!dist = (varhgt - obshgt) / (earth_radius * cos(obslatitude*DEG2RAD))     

! distance in meter
dist = varhgt - obshgt

! find the index in localization%locind
call find_index_in_locind_2d(ilocind, dist, locnumvert, locind)
if ( ilocind == -1 ) then
   print *, 'ERROR in subroutine comp_emp_loc_vert'
   print *, 'ERROR: no index found in localization%obsloc', dist
   return
else
   sublocnum = 1
   sublocnomi = xt_xfmean * deltaxmean
   sublocdenomi = deltaxmean**2
endif

return
end subroutine comp_emp_loc_vert


subroutine find_index_in_vertical(indk1, indk2, bt, varsize3, mhgt, obshgt,     &
                                  cutoffsfc, cutoffuppr)
! find the start/end vertical index (k) for one state variable
integer,     intent(out)  :: indk1, indk2
integer,     intent(in)   :: bt, varsize3
real(r8),    intent(in)   :: mhgt(:)
real(r8),    intent(in)   :: obshgt
real(r8),    intent(in)   :: cutoffsfc, cutoffuppr

integer   :: k

indk1 = -1 
indk2 = -1 

if ( varsize3 == 1 ) then
   ! 2D variable
   if ( abs(obshgt - mhgt(1)) <= cutoffsfc ) then
      indk1 = 1
      indk2 = 1
   endif
   return
else
   ! 3D variable
   if ( obshgt - mhgt(1) <= cutoffuppr ) then
      indk1 = 1
      indk2 = indk1 
  else
      do k = 1, bt-1
         if ( (obshgt - mhgt(k) > cutoffuppr) .and. (obshgt - mhgt(k+1) <= cutoffuppr) ) then
            indk1 = k+1
            indk2 = indk1
            goto 201
         endif
      enddo
   endif
201 continue
   
   if ( mhgt(bt) - obshgt < cutoffuppr ) then
      indk2 = bt
   else
      do k = bt, indk1+1, -1 
         if ( (mhgt(k) - obshgt > cutoffuppr) .and. (mhgt(k-1) -obshgt <= cutoffuppr ) ) then
            indk2 = k-1
            goto 202
         endif
      enddo
   endif
202 continue

   return 

endif

end subroutine find_index_in_vertical


function get_horiz_dist(lon1, lat1, lon2, lat2)
real(r8),     intent(in)  :: lon1, lat1, lon2, lat2
real(r8)                  :: get_horiz_dist

real(r8)   :: londiff, rtemp

londiff = lon1 - lon2

if ( londiff == 0.0 ) then
   get_horiz_dist = abs(lat2 - lat1)
else
   rtemp =  sin(lat2) * sin(lat1) +                 &
            cos(lat2) * cos(lat1) * cos(londiff) 
   if (rtemp < -1.0) then
      get_horiz_dist = PI
   else if (rtemp > 1.0) then
      get_horiz_dist = 0.0
   else
      get_horiz_dist = acos(rtemp)
   endif
endif

end function get_horiz_dist


subroutine find_index_in_locind(indfind, dist, locnum, locind)
! deal with locind has only positive subset

integer,      intent(out)  :: indfind
real(r8),     intent(in)   :: dist
integer,      intent(in)   :: locnum
real(r8),     intent(in)   :: locind(:)

integer    :: ind

indfind = -1

! 2012-09-06, find index without searching, because the distance array
!             is equally distributed

!!if ( abs(dist) > max(abs(locind(1)),abs(locind(locnum))) ) then
!!   return
!!endif 
!
!if ( dist <= locind(1) ) then
!   indfind = 1
!   return
!else
!   do ind = 2, locnum
!      if ( (dist > locind(ind-1)) .and. (dist <= locind(ind)) ) then
!         indfind = ind
!         return
!      endif
!   enddo
!endif
if ( dist < 0.0_r8 ) then
   print *, 'ERROR: distance is negative'
   return
endif

if ( dist <= 0.000000000001_r8 ) then
   indfind = 1
   return
else
   indfind = CEILING(dist/locrad) + 1
   if ( dist <= locind(indfind-1) .or. dist > locind(indfind) ) then
      print *, 'ERROR: wrong distance index: ', dist, locind(indfind-1), locind(indfind)
   elseif ( indfind > locnum ) then
      print *, 'ERROR: larger locnum needed: ', locnum, locrad, dist, locind(locnum)
   endif
endif


return
end subroutine find_index_in_locind


subroutine find_index_in_locind_2d(indfind, dist, locnumvert, locind)
! deal with locind has both positive and negative subset

! NOTE: locnumvert is the length of half index in locind
integer,      intent(out)  :: indfind
real(r8),     intent(in)   :: dist
integer,      intent(in)   :: locnumvert
real(r8),     intent(in)   :: locind(:)

integer    :: ind

indfind = -1

if ( abs(dist) > max(abs(locind(1)),abs(locind(locnum))) ) then
   return
endif

!do ind = 1, locnumvert
!   if ( (dist >= locind(ind)) .and. (dist < locind(ind+1)) ) then
!         indfind = ind
!   endif
!   return
!enddo
!
!if ( abs(dist) <= 0.00000002 ) then
!   indfind = locnumvert+1
!   return
!endif
!
!do ind = 1, locnumvert
!   if ( (dist > locind(locnumvert+1+ind-1)) .and. (dist <= locind(locnumvert+1+ind)) ) then
!      indfind = locnumvert+1+ind 
!      return
!   endif
!enddo

if ( abs(dist) <= 0.00000002 ) then
   indfind = locnumvert + 1
   return
else
   if ( dist < 0.0 ) then
      do ind = 1, locnumvert
         if ( (dist >= locind(ind)) .and. (dist < locind(ind+1)) ) then
            indfind = ind
            return
         endif
      enddo
   else
      do ind = 1, locnumvert
         if ( (dist > locind(locnumvert+1+ind-1)) .and. (dist <= locind(locnumvert+1+ind)) ) then
            indfind = locnumvert + 1 + ind
            return
         endif
      enddo
   endif
endif

return
end subroutine find_index_in_locind_2d


function comp_cov(numens, xmean, x, ymean, y)
integer,   intent(in)  :: numens
real(r8),  intent(in)  :: xmean
real(r8),  intent(in)  :: x(:)
real(r8),  intent(in)  :: ymean
real(r8),  intent(in)  :: y(:)
real(r8)               :: comp_cov

integer   :: i

comp_cov = 0.0
do i = 1, numens
   comp_cov = comp_cov + ( x(i) - xmean ) * ( y(i) - ymean ) 
enddo
comp_cov = comp_cov / (real(numens - 1))

end function comp_cov


subroutine get_pressure(P0,PS,nlev,hyam,hybm,pressure)
! compure the pressure in model levels
real(r8),    intent(in)  :: P0, PS
integer,     intent(in)  :: nlev
real(r8),    intent(in)  :: hyam(:), hybm(:)
real(r8),    intent(out) :: pressure(:)

integer   :: i

pressure(1:nlev) = missing_r8
do i = 1, nlev
   pressure(i) = hyam(i) * P0 + hybm(i) * PS
enddo

end subroutine get_pressure


subroutine get_expected_tpw(pressure,qvapor,nlev,tpw,pressure_top)
! adapted from obs_def_tpw_mod.f90
! compute total precipitation water by model levels
real(r8),      intent(in)  :: pressure(:)
real(r8),      intent(in)  :: qvapor(:)
integer,       intent(in)  :: nlev
real(r8),      intent(out) :: tpw
real(r8),      intent(in)  :: pressure_top

integer    :: i,k,lastk
real(r8)   :: density

density = 1000.0_r8   ! water density in kg/m^3

tpw = missing_r8

lastk = 0

LEVELS: do k = 1, nlev

!print *, k, lastk, pressure(k), qvapor(k) 
   if ( pressure(k) == missing_r8 .or. pressure(k) < pressure_top) exit LEVELS

   if ( qvapor(k) == missing_r8 ) return

   lastk = lastk + 1

enddo LEVELS

if ( lastk == 1 ) then
   return
endif

!print *, 'lastk = ', lastk

! whichever way the column was made (pressure levels or model levels),
! sum the values in the column, computing the area under the curve.
! pressure is in pascals (not hPa or mb), and moisture is in kg/kg.
tpw = 0.0
do k = 1, lastk - 1
   tpw = tpw + 0.5 * (qvapor(k) + qvapor(k+1) ) * (pressure(k) - pressure(k+1) )
enddo

! convert to centimeters of water and return
tpw = 100.0 * tpw /(density*gravity)   ! -> cm

end subroutine get_expected_tpw






end PROGRAM emp_localization



