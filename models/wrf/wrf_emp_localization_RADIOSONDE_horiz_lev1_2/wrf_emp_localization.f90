
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

! wrf_emp_localization.f90_ELF is the code for 1 group with computation of ELF
! 2012-07-10: update the code for global group filter


integer           :: emp_loc_kind     = 1
integer           :: locnum           = 2000
real(r8)          :: locrad           = 0.0002  ! ~1km 
integer           :: group_size       = 2       ! number of groups 
integer           :: wrf_ndom         = 1
integer           :: dom_start        = 1
integer           :: nbdy_nouse       = 10
integer           :: num_state_vars   = 20
integer           :: ens_size         = 60   ! change ens_copy when you change ens_size
integer           :: ens_copy         = 64   ! ens mean + ens spd + ens_size + inf mean + inf spd
integer           :: obslev_num       = 1     ! change obslev_ind below
integer           :: varlev_num       = 1     ! change varlev_ind below
!integer           :: obskind_num      = 4
!integer           :: varkind_num      = 7
! integer           :: lev_s            = 1
! integer           :: lev_e            = 1
 integer           :: indobs_s         = 2
 integer           :: indobs_e         = 4
 integer           :: indvar_s         = 6   !6
 integer           :: indvar_e         = 13   !13

real              :: R_P              = 40000.0_r8 ! obs error variance of Psfc
real              :: R_T              = 1.0_r8     ! obs error variance of Temperature
real              :: R_U              = 4.0_r8     ! obs error variance of Wind


!namelist /wrf_emp_localization_nml/emp_loc_kind, wrf_ndom,   &
!         cutoff_sfcobs, cutoff_upprobs 


character(len=80)  :: truth_name_input  = 'truth.nc',        &
                      prior_name_input  = 'prior',        &
                      obs_name_input    = 'obs.nc'
character(len=80),allocatable  :: prior_name_tmp(:)


type wrf_data 
   integer  :: bt, bts, sn, sns, we, wes, sls
   integer  :: ncopy, nmhgt           ! ncopy, copies of truth and prior; nmhgt, copies of model height
   real(r8) :: dx, dy                           ! dt, p_top
   integer  :: map_proj
   real(r8) :: cen_lat,cen_lon
!   truelat1,        &
!               truelat2,stdlon,re_m
   type(proj_info) :: proj
   real(r8), dimension(:),     pointer :: znu, znw, dn, dnw      ! zs
   real(r8), dimension(:,:),   pointer :: mub, hgt
   real(r8), dimension(:,:),   pointer :: latitude, latitude_u, latitude_v
   real(r8), dimension(:,:),   pointer :: longitude, longitude_u, longitude_v
   real(r8), dimension(:,:,:), pointer :: phb

!   integer :: type_u, type_v, type_w, type_t, type_qv, type_qr, type_hdiab, &
!              type_qndrp, type_qnsnow, type_qnrain, type_qngraupel, type_qnice, &
!              type_qc, type_qg, type_qi, type_qs, type_gz, type_refl, type_fall_spd
!   integer :: type_u10, type_v10, type_t2, type_th2, type_q2, &
!              type_ps, type_mu, type_tsk, type_tslb, type_sh2o, type_smois

   integer :: number_of_state_variables
!   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size
!   integer, dimension(:),   pointer :: var_type
!   integer, dimension(:),   pointer :: var_index_list
!   integer, dimension(:),   pointer :: dart_kind
!   integer, dimension(:,:), pointer :: land
!   real(r8), dimension(:), pointer  :: lower_bound,upper_bound
!   character(len=10), dimension(:),pointer :: clamp_or_fail
   character(len=129),dimension(:),pointer :: description
!   character(len=129),dimension(:),pointer :: units, stagger, coordinates

   integer, dimension(:),   pointer :: mhgtind
   real(r8), dimension(:,:,:,:), pointer   :: mhgt

   real(r8), dimension(:,:,:,:,:), pointer :: variables 

end type wrf_data

type wrf_dom
   type(wrf_data), pointer  :: dom(:)
   integer                  :: model_size
end type wrf_dom

type(wrf_dom)    :: wrf 


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
    real(r8), dimension(:,:,:,:,:), pointer :: obs_val     ! obs values
    real(r8), dimension(:),   pointer :: obsvar            ! obs error variance
    integer,  dimension(:), pointer   :: obs_to_statevar   ! index of state variable for the obs
end type obs_data

type obs_dom
   type(obs_data), pointer  :: dom(:)
   integer                  :: model_size
end type obs_dom

type(obs_dom)   :: obs


type subobs_data
    integer    :: num
    integer,  dimension(:), pointer     :: obsind
    real(r8), dimension(:,:), pointer   :: obsloc
    integer,  dimension(:,:), pointer   :: obsijk
    real(r8), dimension(:,:), pointer   :: obsvp

end type subobs_data

type(subobs_data) :: subobs

type localization_data
    real(r8), dimension(:), pointer         :: locind
    integer,  dimension(:,:,:), pointer     :: num
! 2012-07-10, add two more variables to store globalgroup filter information
    real(r8), dimension(:,:,:), pointer     :: beta1       ! sum beta
    real(r8), dimension(:,:,:), pointer     :: beta2       ! sum beta^2
! 2012-07-10, add one more dimension, save for each group
    real(r8), dimension(:,:,:,:), pointer   :: sumx        ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:,:), pointer   :: sumy        ! used to get mean(y)
    real(r8), dimension(:,:,:,:), pointer   :: sumx2       ! used to get mean(x^2)
    real(r8), dimension(:,:,:,:), pointer   :: sumy2       ! used to get mean(y^2)
    real(r8), dimension(:,:,:,:), pointer   :: numerator   ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:,:), pointer   :: denominator ! i.e., used to get mean(y^2) (note: it contains Robs)
    real(r8), dimension(:,:,:,:), pointer   :: alpha1
    real(r8), dimension(:,:,:,:), pointer   :: alpha2
    real(r8), dimension(:,:,:,:), pointer   :: correlation
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
integer          :: ncid, ncidtruth, ncidprior, ncidobs, iunit
integer          :: proj_code
integer          :: ilocind, sublocnum, writelenloc
integer          :: obslev_ind(50), varlev_ind(50,50)
integer          :: obskind_ind(4), varkind_ind(4,20)

real(r8)         :: stdlon,truelat1,truelat2,latinc,loninc
real(r8)         :: obslat, obslon, obspres, obsi, obsj, obsk, obspalt, obshgt
real(r8)         :: varlat, varlon, varhgt
real(r8)         :: yobs, ytrue, yfmean, oerr, Robs, yprior_var, ypost_var, yvar, yinno_norm,              &
                    xfmean, xtrue, xprior_var, xt_xfmean, covxy, deltaxmean, ccc, ddd, dist
real(r8)         :: obsdiag1, obsdiag2, obsdiag3, obsdiag4
real(r8)         :: sumbeta1, sumbeta2
logical          :: is_lev0
character(len=1) :: idchar, igchar
character(len=2) :: levchar1, levchar2
character(len=NF90_MAX_NAME) :: var_name_tmp, name
character(len=100) :: obsnum_format, locint_format, locreal_format, mhgt_format

integer,  allocatable :: ncidpriors(:)
real(r8), allocatable :: wrf_var_1d(:), wrf_var_2d(:,:), wrf_var_3d(:,:,:), wrf_var_4d(:,:,:,:),           &
                                                         wrf_var_3d_2(:,:,:)
real(r8), allocatable :: obs_var_1d(:), obs_var_2d(:,:), obs_var_3d(:,:,:), obs_var_4d(:,:,:,:)
real(r8), allocatable :: wrf_mu_d01(:,:), wrf_psfc_d01(:,:), wrf_ph_d01(:,:,:),       &
                         wrf_t_d01(:,:,:), wrf_qvapor_d01(:,:,:)
real(r8), allocatable :: xprior(:), yprior(:)
real(r8), allocatable :: xtrueprf(:), xfmeanprf(:), xfspdprf(:), xpriorprf(:,:)
real(r8), allocatable :: xtrue3d(:,:,:), xfmean3d(:,:,:), xfspd3d(:,:,:), xprior3d(:,:,:,:) 
real(r8), allocatable :: numtmp(:)
real(r8), allocatable :: reg_coef(:), corr(:), sublocx(:), sublocy(:), sublocx2(:), sublocy2(:),           &
                         sublocnumer(:), sublocdenom(:), sublocalpha(:)

! define format to write "localization"
locint_format = ' ( 2000i10 ) '
locreal_format = ' ( 2000e15.7 ) '
! define format to write obs number
!obsnum_format = ' ( i10,4e15.7 )'
obsnum_format = ' ( i10 ) '
mhgt_format = ' ( 49e15.7 ) '


! ! define obskind and varkind to be analyzed
! obskind_ind(1)  = 1   ! RADIOSONDE_PRESSURE
! obskind_ind(2)  = 2   ! RADIOSONDE_T
! obskind_ind(3)  = 3   ! RADIOSONDE_U
! obskind_ind(4)  = 4   ! RADIOSONDE_V
! varkind_ind(1)  = 1   
! varkind_ind(2)  = 2
! varkind_ind(3)  = 3
! varkind_ind(4)  = 4
! varkind_ind(5)  = 5
! varkind_ind(6)  = 6
! varkind_ind(7)  = 7   
! varkind_ind(8)  = 8
! varkind_ind(9)  = 9
! varkind_ind(10) = 10
! varkind_ind(11) = 11
! varkind_ind(12) = 12
! varkind_ind(13) = 13
! ! varkind_ind(14) = 14
! varkind_ind(15) = 16
! varkind_ind(16) = 18

! define obs & state variable levels
obslev_ind(1) = 26
varlev_ind(1,1) = 26


wrf%model_size = wrf_ndom

allocate(wrf%dom(wrf_ndom))

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
call nc_check( nf90_open(truth_name_input, NF90_NOWRITE, ncidtruth),   &
               'wrf_emp_localization','open truth file')

!-------------------------------------------------------
! open prior files (Prior_Diag.nc) 
!-------------------------------------------------------
allocate(prior_name_tmp(group_size))
allocate(ncidpriors(group_size))
do i = 1, group_size
   write(igchar,'(i1.1)') i
   prior_name_tmp(i) = trim(prior_name_input)//igchar//'.nc'
   call nc_check( nf90_open(trim(prior_name_tmp(i)), NF90_NOWRITE, ncidpriors(i)),   &
                  'wrf_emp_localization','open prior file')
enddo

WRFDomains : do id = 1, wrf_ndom 

!-------------------------------------------------------
! fill in state variable information 
!-------------------------------------------------------
   wrf%dom(id)%number_of_state_variables = num_state_vars
   allocate(wrf%dom(id)%description(wrf%dom(id)%number_of_state_variables))
   wrf%dom(id)%description(1)  = 'U10'
   wrf%dom(id)%description(2)  = 'V10'
   wrf%dom(id)%description(3)  = 'T2'
   wrf%dom(id)%description(4)  = 'Q2'
   wrf%dom(id)%description(5)  = 'TH2'
   wrf%dom(id)%description(6)  = 'MU'
   wrf%dom(id)%description(7)  = 'PSFC'
   wrf%dom(id)%description(8)  = 'U'
   wrf%dom(id)%description(9)  = 'V'
   wrf%dom(id)%description(10) = 'W'
   wrf%dom(id)%description(11) = 'T'
   wrf%dom(id)%description(12) = 'PH'
   wrf%dom(id)%description(13) = 'QVAPOR'
   wrf%dom(id)%description(14) = 'QCLOUD'
   wrf%dom(id)%description(15) = 'QNRAIN'
   wrf%dom(id)%description(16) = 'QNICE'
   wrf%dom(id)%description(17) = 'QRAIN'
   wrf%dom(id)%description(18) = 'QICE'
   wrf%dom(id)%description(19) = 'QSNOW'
   wrf%dom(id)%description(20) = 'QGRAUP'

   write(idchar,'(i1.1)') id

!-------------------------------------------------------
! read WRF dimensions
!-------------------------------------------------------
   ! dimension of copies (1(truth)+ [ 2(ens mean, ens spread)+ens_size+2(inf mean, inf spread)) ] for each group
!   wrf%dom(id)%ncopy = 1 + (ens_size + 4) * group_size
   wrf%dom(id)%ncopy = 1 + ens_copy * group_size

   var_name_tmp = 'bottom_top_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_wrf_dimensions','inq_dimid bottom_top')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%bt),  &
                  'read_wrf_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'bottom_top_stag_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid bottom_top_stag') ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%bts), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   var_name_tmp = 'south_north_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid south_north')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%sn), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   var_name_tmp = 'south_north_stag_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid south_north_stag') ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%sns), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   var_name_tmp = 'west_east_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid west_east')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%we), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   var_name_tmp = 'west_east_stag_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid west_east_stag')  ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%wes), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   var_name_tmp = 'soil_layers_stag_d0'//idchar
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id), &
                     'read_wrf_dimensions','inq_dimid soil_layers_stag')  ! reuse var_id, no harm
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, wrf%dom(id)%sls), &
                     'read_wrf_dimensions','inquire_dimension '//trim(name))

   call nc_check( nf90_inq_varid(ncidtruth, "DX", var_id), &
                     'read_wrf_dimensions','inq_varid DX')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%dx), &
                     'read_wrf_dimensions','get_var DX')

   call nc_check( nf90_inq_varid(ncidtruth, "DY", var_id), &
                     'read_wrf_dimensions','inq_varid DY')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%dy), &
                     'read_wrf_dimensions','get_var DY')

!-------------------------------------------------------
! read WRF map parameters 
!
! NOTE: for nested domains, CEN_LAT and CEN_LON are different,
!       MAP_PROJ, TRUELAT1, TRUELAT2, STAND_LON are the same
!-------------------------------------------------------
!   wrf%dom(id)%re_m = earth_radius

   call nc_check( nf90_inq_varid(ncidtruth, "MAP_PROJ", var_id), &
                     'read_wrf_map_para','inq_varid MAP_PROJ')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%map_proj), &
                     'read_wrf_map_para','get_var MAP_PROJ')

   call nc_check( nf90_inq_varid(ncidtruth, "CEN_LAT", var_id), &
                     'read_wrf_map_para','inq_varid CEN_LAT')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%cen_lat), &
                     'read_wrf_map_para','get_var CEN_LAT')

   call nc_check( nf90_inq_varid(ncidtruth, "CEN_LON", var_id), &
                     'read_wrf_map_para','inq_varid CEN_LON')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%cen_lon), &
                     'read_wrf_map_para','get_var CEN_LON')
!   if ( wrf%dom(id)%cen_lon < 0.0 )         & 
!        wrf%dom(id)%cen_lon = wrf%dom(id)%cen_lon + 360.0

   call nc_check( nf90_inq_varid(ncidtruth, "TRUELAT1", var_id), &
                     'read_wrf_map_para','inq_varid TRUELAT1')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, truelat1), &
                     'read_wrf_map_para','get_var TRUELAT1')

   call nc_check( nf90_inq_varid(ncidtruth, "TRUELAT2", var_id), &
                     'read_wrf_map_para','inq_varid TRUELAT2')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, truelat2), &
                     'read_wrf_map_para','get_var TRUELAT2')

   call nc_check( nf90_inq_varid(ncidtruth, "STAND_LON", var_id), &
                     'read_wrf_map_para','inq_varid STAND_LON')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidtruth, var_id, stdlon), &
                     'read_wrf_map_para','get_var STAND_LON')

!-------------------------------------------------------
! read WRF static data 
!-------------------------------------------------------
! 1. 1D array

   allocate(wrf%dom(id)%dn(1:wrf%dom(id)%bt))
   var_name_tmp = 'DN_d0'//idchar
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid DN')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%dn),  &
                     'read_wrf_static_data','get_var DN')

   var_name_tmp = 'DNW_d0'//idchar
   allocate(wrf%dom(id)%dnw(1:wrf%dom(id)%bt))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid DNW')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%dnw), &
                     'read_wrf_static_data','get_var DNW')

   var_name_tmp = 'ZNU_d0'//idchar
   allocate(wrf%dom(id)%znu(1:wrf%dom(id)%bt))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid ZNU')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%znu), &
                     'read_wrf_static_data','get_var ZNU')

   var_name_tmp = 'ZNW_d0'//idchar
   allocate(wrf%dom(id)%znw(1:wrf%dom(id)%bts))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid ZNW')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%znw), &
                     'read_wrf_static_data','get_var ZNW')

! 2. 2D array

   var_name_tmp = 'MUB_d0'//idchar
   allocate(wrf%dom(id)%mub(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid MUB')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%mub), &
                     'read_wrf_static_data','get_var MUB')

   var_name_tmp = 'HGT_d0'//idchar
   allocate(wrf%dom(id)%hgt(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_wrf_static_data','inq_varid HGT')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%hgt), &
                     'read_wrf_static_data','get_var HGT')

   var_name_tmp = 'XLAT_d0'//idchar
   allocate(wrf%dom(id)%latitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLAT')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%latitude), &
                     'read_wrf_static_data','get_var XLAT')

   var_name_tmp = 'XLAT_U_d0'//idchar
   allocate(wrf%dom(id)%latitude_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLAT_U')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%latitude_u), &
                     'read_wrf_static_data','get_var XLAT_U')

   var_name_tmp = 'XLAT_V_d0'//idchar
   allocate(wrf%dom(id)%latitude_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLAT_V')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%latitude_v), &
                     'read_wrf_static_data','get_var XLAT_V')

   var_name_tmp = 'XLONG_d0'//idchar
   allocate(wrf%dom(id)%longitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLONG')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%longitude), &
                     'read_wrf_static_data','get_var XLONG')
!   do i = 1, wrf%dom(id)%we
!      do j = 1, wrf%dom(id)%sn
!         if ( wrf%dom(id)%longitude(i,j) < 0.0 ) then
!              wrf%dom(id)%longitude(i,j) = wrf%dom(id)%longitude(i,j) + 360.0
!         endif
!      enddo
!   enddo

   var_name_tmp = 'XLONG_U_d0'//idchar
   allocate(wrf%dom(id)%longitude_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLONG_U')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%longitude_u), &
                     'read_wrf_static_data','get_var XLONG_U')
!   do i = 1, wrf%dom(id)%wes
!      do j = 1, wrf%dom(id)%sn
!         if ( wrf%dom(id)%longitude_u(i,j) < 0.0 ) then
!              wrf%dom(id)%longitude_u(i,j) = wrf%dom(id)%longitude_u(i,j) + 360.0
!         endif
!      enddo
!   enddo

   var_name_tmp = 'XLONG_V_d0'//idchar
   allocate(wrf%dom(id)%longitude_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLONG_V')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%longitude_v), &
                     'read_wrf_static_data','get_var XLONG_V')
!   do i = 1, wrf%dom(id)%we
!      do j = 1, wrf%dom(id)%sns
!         if ( wrf%dom(id)%longitude_v(i,j) < 0.0 ) then
!              wrf%dom(id)%longitude_v(i,j) = wrf%dom(id)%longitude_v(i,j) + 360.0
!         endif
!      enddo
!   enddo

! 3. 3D array

   var_name_tmp = 'PHB_d0'//idchar
   allocate(wrf%dom(id)%phb(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%bts))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid PHB')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%phb), &
                     'read_wrf_static_data','get_var PHB')

!-------------------------------------------------------
! read WRF variables of True_State.nc (copy 1) 
!-------------------------------------------------------
   ! var_size(1,:) - dimension of we, 
   ! var_size(2,:) - dimension of sn,
   ! var_size(3,:) - dimension of bt ( = 1, means 2D variable, others 3D variable)
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_state_variables))
   allocate(wrf%dom(id)%variables(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,    &
                                  wrf%dom(id)%ncopy,wrf%dom(id)%number_of_state_variables))

   do ind = 1,wrf%dom(id)%number_of_state_variables

      ! Get the dimension size (ignore the dimension of time) 
      ! Once this is done for variable of True_State.nc, there is no need to do this for
      ! variables in Prior_Diag.nc, since they should be consistent
      var_name_tmp = trim(wrf%dom(id)%description(ind))//'_d0'//idchar
      call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'get_variable_size_from_file',                    &
                     'inq_varid '//var_name_tmp)
      call nc_check( nf90_inquire_variable(ncidtruth, var_id,          &
                     ndims=ndims, dimids=dimids),                      &
                     'get_variable_size_from_file',                    &
                     'inquire_variable '//var_name_tmp)
      do inddim = 1,ndims-2
         call nc_check( nf90_inquire_dimension(ncidtruth, dimids(inddim),  &
                        len=wrf%dom(id)%var_size(inddim,ind)),             &
                        'get_variable_size_from_file',                     &
                        'inquire_dimension '//var_name_tmp)
      enddo

      if ( ndims < 5 ) then
         ! 2D variable, its vertical dimension is 1
         wrf%dom(id)%var_size(ndims-1,ind) = 1

         ! Get 2D variable
         allocate(wrf_var_2d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind)))
         call nc_check( nf90_get_var(ncidtruth, var_id, wrf_var_2d),       &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1,                                          &
                               1,                                          &
                               ind) =                                      &
                    wrf_var_2d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind))
         deallocate(wrf_var_2d)

      else

         ! Get 3D variable
         allocate(wrf_var_3d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind),                  &
                             wrf%dom(id)%var_size(3,ind)))
         call nc_check( nf90_get_var(ncidtruth, var_id, wrf_var_3d),       &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind),              &
                               1,                                          &
                               ind) =                                      &
                    wrf_var_3d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind))

         deallocate(wrf_var_3d)

      endif

   enddo

!-------------------------------------------------------
! read WRF variables of Prio_Diag.nc (copy 2-65) 
! read each group
!-------------------------------------------------------
   do ind = 1,wrf%dom(id)%number_of_state_variables

      do i = 1, group_size

         igs = 1 + (i-1)*ens_copy + 1
         ige = 1 + i*ens_copy

         if ( i==group_size .and. ige/=wrf%dom(id)%ncopy ) then
            print *, 'ERROR in read priors, index WRONG'
         endif

         var_name_tmp = trim(wrf%dom(id)%description(ind))//'_d0'//idchar

         call nc_check( nf90_inq_varid(ncidpriors(i), var_name_tmp, var_id),  &
                        'get_variable_name_from_file',                        &
                        'inq_varid '//var_name_tmp)

         if ( wrf%dom(id)%var_size(3,ind) == 1 ) then

            ! Get 2D variable 
            ! NOTE: nf90_get_var automatically ignor the dimension with length 1
            allocate(wrf_var_3d(wrf%dom(id)%var_size(1,ind),                  &
                                wrf%dom(id)%var_size(2,ind),                  &
                                ens_copy))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, wrf_var_3d),   &
                        'read_in_wrf_var','get_var '//var_name_tmp )
            wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                                  1:wrf%dom(id)%var_size(2,ind),              &
                                  1,                                          &
                                  igs:ige,                                    &
                                  ind) =                                      &
                       wrf_var_3d(1:wrf%dom(id)%var_size(1,ind),              &
                                  1:wrf%dom(id)%var_size(2,ind),              &
                                  1:ens_copy)
            deallocate(wrf_var_3d)

         else
   
            ! Get 3D variable
            allocate(wrf_var_4d(wrf%dom(id)%var_size(1,ind),                  &
                                wrf%dom(id)%var_size(2,ind),                  &
                                wrf%dom(id)%var_size(3,ind),                  &
                                ens_copy))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, wrf_var_4d),   &
                        'read_in_wrf_var','get_var '//var_name_tmp )
            wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                                  1:wrf%dom(id)%var_size(2,ind),              &
                                  1:wrf%dom(id)%var_size(3,ind),              &
                                  igs:ige,                                    &
                                  ind) =                                      &
                       wrf_var_4d(1:wrf%dom(id)%var_size(1,ind),              &
                                  1:wrf%dom(id)%var_size(2,ind),              &
                                  1:wrf%dom(id)%var_size(3,ind),              &
                                  1:ens_copy)

            deallocate(wrf_var_4d)
   
         endif

      enddo

   enddo


!-------------------------------------------------------
! setup map projection 
!-------------------------------------------------------
   call map_init(wrf%dom(id)%proj)
   
   latinc = 180.0_r8/wrf%dom(id)%sn
   loninc = 180.0_r8/wrf%dom(id)%we
   
   call map_set( proj_code=PROJ_LC,               &
                 proj=wrf%dom(id)%proj,           &
                 lat1=wrf%dom(id)%latitude(1,1),  &
                 lon1=wrf%dom(id)%longitude(1,1), &
                 lat0=90.0_r8,                    &
                 lon0=0.0_r8,                     &
                 knowni=1.0_r8,                   &
                 knownj=1.0_r8,                   &
                 dx=wrf%dom(id)%dx,               &
                 latinc=latinc,                   &
                 loninc=loninc,                   &
                 stdlon=stdlon,                   &
                 truelat1=truelat1,               &
                 truelat2=truelat2  )

!-------------------------------------------------------
! compute height on model grid points for later use
!-------------------------------------------------------
! 1. sigh the model height index to the number_of_state_varaibles
  wrf%dom(id)%nmhgt = 4
  allocate(wrf%dom(id)%mhgtind(wrf%dom(id)%number_of_state_variables))

  do ind = 1,wrf%dom(id)%number_of_state_variables
     if ( (trim(wrf%dom(id)%description(ind)) == "MU") .or.        &
          (trim(wrf%dom(id)%description(ind)) == "PSFC")  ) then
         wrf%dom(id)%mhgtind(ind) = 1
     elseif ( (trim(wrf%dom(id)%description(ind)) == "U10") .or.   &
              (trim(wrf%dom(id)%description(ind)) == "V10") ) then
         wrf%dom(id)%mhgtind(ind) = 2
     elseif ( (trim(wrf%dom(id)%description(ind)) == "T2") .or.    & 
              (trim(wrf%dom(id)%description(ind)) == "Q2") .or.    &
              (trim(wrf%dom(id)%description(ind)) == "TH2") ) then
         wrf%dom(id)%mhgtind(ind) = 3
!     elseif ( trim(wrf%dom(id)%description(ind)) == "U" ) then
!         wrf%dom(id)%mhgtind(ind) = 4
!     elseif ( trim(wrf%dom(id)%description(ind)) == "V" ) then
!         wrf%dom(id)%mhgtind(ind) = 5
!     elseif ( (trim(wrf%dom(id)%description(ind)) == "W") .or.     &
!              (trim(wrf%dom(id)%description(ind)) == "PH") ) then
!         wrf%dom(id)%mhgtind(ind) = 6
!     else
!         wrf%dom(id)%mhgtind(ind) = 7
     else
         wrf%dom(id)%mhgtind(ind) = 4
     endif
  enddo

! 2. compute the model height given each model height type
  allocate(wrf%dom(id)%mhgt(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,    &
                            wrf%dom(id)%nmhgt))
  allocate(wrf_var_3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts))

! NOTE: give truth PH to compute the model height (truth copy - 1, PH index - 12)
  wrf_var_3d(1:wrf%dom(id)%we, 1:wrf%dom(id)%sn, 1:wrf%dom(id)%bts) =           &
     wrf%dom(id)%variables(1:wrf%dom(id)%we, 1:wrf%dom(id)%sn, 1:wrf%dom(id)%bts, 1, 12)
 
  call compute_model_height(wrf%dom(id)%we, wrf%dom(id)%wes,                    &
                            wrf%dom(id)%sn, wrf%dom(id)%sns,                    &
                            wrf%dom(id)%bt, wrf%dom(id)%bts,                    &
                            wrf%dom(id)%hgt, wrf%dom(id)%phb, wrf_var_3d,       &
                            wrf%dom(id)%longitude, wrf%dom(id)%latitude,        &
                            wrf%dom(id)%nmhgt, wrf%dom(id)%mhgt) 

!  wrf%dom(id)%mhgt(1:wrf%dom(id)%wes, 1:wrf%dom(id)%sns,                        &
!                   1: wrf%dom(id)%bts, 1:wrf%dom(id)%nmhgt) =                   &
!        wrf_var_4d(1:wrf%dom(id)%wes, 1:wrf%dom(id)%sns,                        &
!                   1: wrf%dom(id)%bts, 1:wrf%dom(id)%nmhgt)           
  deallocate(wrf_var_3d)

enddo WRFDomains


!-------------------------------------------------------
! close netcdf files 
!-------------------------------------------------------
call nc_check(nf90_close(ncidtruth),'wrf_emp_localization','close truth file')
do i = 1, group_size
   call nc_check(nf90_close(ncidpriors(i)),'wrf_emp_localization','close prior file')
enddo
deallocate(prior_name_tmp, ncidpriors)


!!-------------------------------------------------------
!! find how many obs types in obs_epoch file 
!!-------------------------------------------------------
!print *, 'finding the number of obs types'
!allocate(obs%obstypeindex(100))
!
!call count_numobstype(obs%obsindex, obs%obstype, obs%obsqc, obs%numobstype, obs%obstypeindex)


!-------------------------------------------------------
! define identity obs information for each domain 
!-------------------------------------------------------
obs%model_size = wrf_ndom
allocate(obs%dom(wrf_ndom))

do id = 1, wrf_ndom 

!   ! 1. define obs copy size: truth, prior mean, prior spread, mem1, mem2, ..., error variance
!   obs%dom(id)%obscopysize = ens_size + 4

   ! 2. define obs type
   obs%dom(id)%type_size = 4 
   allocate(obs%dom(id)%obstype_metadata(obs%dom(id)%type_size))
!   obs%dom(id)%obstype_metadata(1)   = 'SFC_ALTIMETER'
!   obs%dom(id)%obstype_metadata(2)   = 'SFC_U_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(3)   = 'SFC_V_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(4)   = 'SFC_TEMPERATURE'
!   obs%dom(id)%obstype_metadata(5)   = 'SFC_DEWPOINT'
!   obs%dom(id)%obstype_metadata(6)   = 'SFC_SPECIFIC_HUMIDITY'
!   obs%dom(id)%obstype_metadata(7)   = 'RADIOSONDE_SURFACE_ALTIMETER'
   obs%dom(id)%obstype_metadata(1)   = 'RADIOSONDE_PRESSURE' 
   obs%dom(id)%obstype_metadata(2)   = 'RADIOSONDE_TEMPERATURE'
   obs%dom(id)%obstype_metadata(3)   = 'RADIOSONDE_U_WIND_COMPONENT'
   obs%dom(id)%obstype_metadata(4)   = 'RADIOSONDE_V_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(11)  = 'RADIOSONDE_DEWPOINT'
!   obs%dom(id)%obstype_metadata(12)  = 'SAT_U_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(13)  = 'SAT_V_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(14)  = 'ACARS_U_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(15)  = 'ACARS_V_WIND_COMPONENT'
!   obs%dom(id)%obstype_metadata(16)  = 'ACARS_TEMPERATURE'
!   obs%dom(id)%obstype_metadata(17)  = 'ACARS_DEWPOINT'

   allocate(obs%dom(id)%obsvar(obs%dom(id)%type_size))
   obs%dom(id)%obsvar(1) = R_P
   obs%dom(id)%obsvar(2) = R_T
   obs%dom(id)%obsvar(3) = R_U
   obs%dom(id)%obsvar(4) = R_U


   allocate(obs%dom(id)%obs_to_statevar(obs%dom(id)%type_size))
   do i = 1, obs%dom(id)%type_size
      if ( trim(obs%dom(id)%obstype_metadata(i)) == "RADIOSONDE_PRESSURE" ) then
         obs%dom(id)%obs_to_statevar(i) = 7
      elseif ( trim(obs%dom(id)%obstype_metadata(i)) == "RADIOSONDE_TEMPERATURE" ) then
         obs%dom(id)%obs_to_statevar(i) = 11
      elseif ( trim(obs%dom(id)%obstype_metadata(i)) == "RADIOSONDE_U_WIND_COMPONENT" ) then
         obs%dom(id)%obs_to_statevar(i) = 8
      elseif ( trim(obs%dom(id)%obstype_metadata(i)) == "RADIOSONDE_V_WIND_COMPONENT" ) then
         obs%dom(id)%obs_to_statevar(i) = 9
      else
         print *, 'ERROR: obs_to_statevar'
      endif
   enddo

!   ! 3. define obs size (obs_size is similar to var_size)
!   ! obs_size(1,:) - dimension of we, 
!   ! obs_size(2,:) - dimension of sn,
!   ! obs_size(3,:) - dimension of bt ( = 1, means 2D variable, others 3D variable) 
!   allocate(obs%dom(id)%obs_size(3,obs%dom(id)%type_size))
!!   allocate(obs%dom(id)%obs_err_var(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt,     &
!!                                    obs%dom(id)%type_size)
!   allocate(obs%dom(id)%obs_val(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt, &
!                                obs%dom(id)%obscopysize,obs%dom(id)%type_size))
!   allocate(wrf_var_2d(wrf%dom(id)%we,wrf%dom(id)%sn))
!   allocate(wrf_var_3d(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%ncopy))
!   allocate(wrf_var_3d_2(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%ncopy))
!   allocate(obs_var_3d(wrf%dom(id)%we,wrf%dom(id)%sn,obs%dom(id)%obscopysize))
!
!! Lili test
!!   do ind = 1, obs%dom(id)%type_size 
!   do ind = 1, 1
!! print *, 'ind = ', ind
!      obs%dom(id)%obs_size(1,ind) = wrf%dom(id)%we
!      obs%dom(id)%obs_size(2,ind) = wrf%dom(id)%sn
!
!      if ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_ALTIMETER' )  then  
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! PSFC index --- 7
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,7)
!         ! HGT
!         wrf_var_2d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn) =                               &
!              wrf%dom(id)%hgt(1:wrf%dom(id)%we,1:wrf%dom(id)%sn)
!         call draw_perfect_obs_altimeter(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,          &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_2d,wrf_var_3d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_U_WIND_COMPONENT' ) then
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! U10 index --- 1
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,1)
!         call draw_perfect_obs_2d(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,  &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_3d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_V_WIND_COMPONENT' ) then
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! V10 index --- 2
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,2)
!         call draw_perfect_obs_2d(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,  &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_3d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_TEMPERATURE' ) then
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! T2 index --- 3
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,3)
!         call draw_perfect_obs_2d(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,  &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_3d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_DEWPOINT' ) then
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! PSFC index --- 7
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,7)
!         ! T2 index --- 3, only use 1st copy 
!         wrf_var_2d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn) =                               &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1,3)
!         ! Q2 index --- 4
!         wrf_var_3d_2(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =         &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,4)
!         call draw_perfect_obs_2d_dewpoint(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,        &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_3d,wrf_var_3d_2,wrf_var_2d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SFC_SPECIFIC_HUMIDITY' ) then
!         obs%dom(id)%obs_size(3,ind) = 1
!         ! PSFC index --- 7
!         wrf_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =           &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,7)
!         ! T2 index --- 3, only use 1st copy 
!         wrf_var_2d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn) =                               &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1,3)
!         ! Q2 index --- 4
!         wrf_var_3d_2(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%ncopy) =         &
!              wrf%dom(id)%variables(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:wrf%dom(id)%ncopy,4)
!         call draw_perfect_obs_2d_q(wrf%dom(id)%we,wrf%dom(id)%sn,nbdy_nouse, ens_size,        &
!                    obs%dom(id)%obstype_metadata(ind),wrf_var_3d,wrf_var_3d_2,wrf_var_2d,obs_var_3d)
!         obs%dom(id)%obs_val(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1,1:obs%dom(id)%obscopysize,ind) =     &
!              obs_var_3d(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:obs%dom(id)%obscopysize)
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'RADIOSONDE_SURFACE_ALTIMETER' ) then
!
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'RADIOSONDE_U_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'RADIOSONDE_V_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'RADIOSONDE_TEMPERATURE' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'RADIOSONDE_DEWPOINT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SAT_U_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'SAT_V_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'ACARS_U_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'ACARS_V_WIND_COMPONENT' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'ACARS_TEMPERATURE' ) then
!
!      elseif ( trim(obs%dom(id)%obstype_metadata(ind)) == 'ACARS_DEWPOINT' ) then
!
!
!         obs%dom(id)%obs_size(3,ind) = wrf%dom(id)%bt
!      endif
!
!   enddo
!
!   deallocate(wrf_var_2d, wrf_var_3d, wrf_var_3d_2, obs_var_3d)

enddo


! Initialize locind which contains the radiance for each subset
   print *, 'empirical localization kind: horizontal'
   allocate(localization%locind(locnum))
   do ind = 1, locnum
      localization%locind(ind) = locrad * real(ind)
   enddo

!-------------------------------------------------------
! Big loop: for each domain  
!-------------------------------------------------------
LoopDomain: do id = dom_start, wrf_ndom 

   LoopObsTypes: do idobs = indobs_s, indobs_e
   print *, 'Now processing obs type ', trim(obs%dom(id)%obstype_metadata(idobs))

      ndims = 2 * (obslev_num * varlev_num)    ! 2 for whole domain and ocean 
      allocate(localization%num(locnum,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%beta1(locnum,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%beta2(locnum,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%sumx(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%sumy(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%sumx2(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%sumy2(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%numerator(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%denominator(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%alpha1(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%alpha2(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))
      allocate(localization%correlation(locnum,group_size,ndims,wrf%dom(id)%number_of_state_variables))

      localization%num(1:locnum,1:ndims,1:wrf%dom(id)%number_of_state_variables)                      = 0
      localization%beta1(1:locnum,1:ndims,1:wrf%dom(id)%number_of_state_variables)                    = 0.0
      localization%beta2(1:locnum,1:ndims,1:wrf%dom(id)%number_of_state_variables)                    = 0.0
      localization%sumx(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)        = 0.0
      localization%sumy(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)        = 0.0
      localization%sumx2(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)       = 0.0
      localization%sumy2(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)       = 0.0
      localization%numerator(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)   = 0.0
      localization%denominator(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables) = 0.0
      localization%alpha1(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)      = 0.0
      localization%alpha2(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables)      = 0.0
      localization%correlation(1:locnum,1:group_size,1:ndims,1:wrf%dom(id)%number_of_state_variables) = 0.0



!$OMP PARALLEL DO PRIVATE(ind,i,j,k,ii,jj,kk,obslevk,varlevk,inddim,ndims,indk1,indk2,ig,igs,ige,igg,iunit,writelenloc,var_name_tmp,obslon,obslat,varlon,varlat,varhgt,obsi,obsj,obsk,obspres,obshgt,is_lev0,yobs,ytrue,yfmean,Robs,yvar,yprior_var,ypost_var,corr,reg_coef,ccc,ddd,dist,yinno_norm,xtrue,xfmean,xprior_var,xt_xfmean,covxy,deltaxmean,ilocind,sublocx,sublocy,sublocx2,sublocy2,sublocnum,sublocnumer,sublocalpha,sublocdenom,wrf_var_1d,xprior,yprior,xtrueprf,xfmeanprf,xfspdprf,xpriorprf,xtrue3d,xfmean3d,xfspd3d,xprior3d,sumbeta1,sumbeta2)

      LoopVarTypes: do inddim = indvar_s, indvar_e
!      LoopVarTypes: do inddim = 6, wrf%dom(id)%number_of_state_variables 
         print *, 'processing variable ', trim(wrf%dom(id)%description(inddim))

      allocate(yprior(ens_size))
      allocate(xprior(ens_size))

      allocate(corr(group_size),reg_coef(group_size))
      allocate(sublocx(group_size),sublocy(group_size))
      allocate(sublocx2(group_size),sublocy2(group_size))
      allocate(sublocnumer(group_size),sublocdenom(group_size),sublocalpha(group_size))

      ! loop obs grid points
!      do k = lev_s, lev_e
      do obslevk = 1, obslev_num
         k = obslev_ind(obslevk)
!         print *, 'obs lev = ', obslevk, k
         if ( idobs == 1 ) then
              k = 1
              print *, 'OBS is ', trim(obs%dom(id)%obstype_metadata(idobs))
              print *, 'NOTE: change obs lev index ', obslev_ind(obslevk), k
         endif

         do j = nbdy_nouse+1, wrf%dom(id)%sn-nbdy_nouse
            do i = nbdy_nouse+1, wrf%dom(id)%we-nbdy_nouse

               ind = obs%dom(id)%obs_to_statevar(idobs)

               ! loop target variable grid points
!               do kk = lev_s, lev_e
               do varlevk = 1, varlev_num
                  kk = varlev_ind(obslevk,varlevk)
!                  print *, 'var lev = ', varlevk, kk, ndims
                  if ( inddim == 12 .and. kk == 1 ) then
                       kk = 2
                       if ( j == nbdy_nouse+1 .and. i == nbdy_nouse+1 ) then
                       print *, 'VAR is ', trim(wrf%dom(id)%description(inddim))
                       print *, 'NOTE: change var lev index ', varlev_ind(obslevk,varlevk), kk 
                       endif
                  elseif ( inddim < 8 ) then
                       kk = 1
                       if ( j == nbdy_nouse+1 .and. i == nbdy_nouse+1 ) then
                       print *, 'VAR is ', trim(wrf%dom(id)%description(inddim))
                       print *, 'NOTE: change var lev index ', varlev_ind(obslevk,varlevk), kk
                       endif
                  endif

                  ndims = (obslevk - 1) * varlev_num + varlevk
                  do jj = nbdy_nouse+1, wrf%dom(id)%sn-nbdy_nouse 
                     do ii = nbdy_nouse+1, wrf%dom(id)%we-nbdy_nouse

                        ! computation for each group
                        do ig = 1, group_size

                           igs = 1 + (ig-1)*ens_copy + 1
                           ige = 1 + ig*ens_copy
!print *, 'ig index ', igs, ige
!pause

                           ! get obs information ready
                           ytrue      = wrf%dom(id)%variables(i,j,k,1,ind)
                           yfmean     = wrf%dom(id)%variables(i,j,k,igs,ind)
                           yprior_var = wrf%dom(id)%variables(i,j,k,igs+1,ind)**2
                           Robs       = obs%dom(id)%obsvar(idobs)
                           ypost_var    = 1.0_r8/(1.0_r8/yprior_var + 1.0_r8/Robs)
                           yprior(1:ens_size) = wrf%dom(id)%variables(i,j,k,igs+1+1:igs+1+ens_size,ind)

                           ! get [lon, lat] and convert to radian
                           if ( (trim(obs%dom(id)%obstype_metadata(idobs)) == "RADIOSONDE_U_WIND_COMPONENT") ) then
                              obslon = wrf%dom(id)%longitude_u(i,j) * DEG2RAD
                              obslat = wrf%dom(id)%latitude_u(i,j) * DEG2RAD
                           elseif ( (trim(obs%dom(id)%obstype_metadata(idobs)) == "RADIOSONDE_V_WIND_COMPONENT") ) then
                              obslon = wrf%dom(id)%longitude_v(i,j) * DEG2RAD
                              obslat = wrf%dom(id)%latitude_v(i,j) * DEG2RAD
                           else
                              obslon = wrf%dom(id)%longitude(i,j) * DEG2RAD
                              obslat = wrf%dom(id)%latitude(i,j) * DEG2RAD
                           endif

                           obshgt = wrf%dom(id)%mhgt(i,j,1,1)

!print *, 'y info ', ytrue, yfmean, yprior_var
!print *, 'y info ', yprior(1:ens_size:10)
!print *, 'y loc  ', wrf%dom(id)%longitude(i,j), wrf%dom(id)%latitude(i,j), obshgt
!pause

                           ! get target variable information ready
                           xtrue      = wrf%dom(id)%variables(ii,jj,kk,1,inddim)
                           xfmean     = wrf%dom(id)%variables(ii,jj,kk,igs,inddim)
                           xprior_var = wrf%dom(id)%variables(ii,jj,kk,igs+1,inddim)**2
                           xprior(1:ens_size) = wrf%dom(id)%variables(ii,jj,kk,igs+1+1:igs+1+ens_size,inddim)

                           ! get [lon, lat] and convert to radian
                           ! get [lon, lat] and convert to radian
                           if ( (trim(wrf%dom(id)%description(inddim)) == "U") ) then
                              varlon = wrf%dom(id)%longitude_u(ii,jj) * DEG2RAD
                              varlat = wrf%dom(id)%latitude_u(ii,jj) * DEG2RAD
                           elseif ( (trim(wrf%dom(id)%description(inddim)) == "V") ) then
                              varlon = wrf%dom(id)%longitude_v(ii,jj) * DEG2RAD
                              varlat = wrf%dom(id)%latitude_v(ii,jj) * DEG2RAD
                           else
                              varlon = wrf%dom(id)%longitude(ii,jj) * DEG2RAD
                              varlat = wrf%dom(id)%latitude(ii,jj) * DEG2RAD
                           endif

                           varhgt = wrf%dom(id)%mhgt(ii,jj,1,1)
!print *, 'x info ', xtrue, xfmean, xprior_var
!print *, 'x info ', xprior(1:ens_size:10)
!print *, 'x loc  ', wrf%dom(id)%longitude(ii,jj), wrf%dom(id)%latitude(ii,jj), varhgt
!pause


                           ! get the horizontal distance
                           dist = get_horiz_dist(varlon, varlat, obslon, obslat)
 
                           call find_index_in_locind(ilocind, dist, locnum, localization%locind)
!print *, 'dist diag ', dist
!print *, 'index diag ', ilocind
!pause

                           ! localization computation
                           xt_xfmean = xtrue - xfmean
                           covxy = comp_cov(ens_size,xfmean,xprior,yfmean,yprior)
                           corr(ig) = covxy / (sqrt(yprior_var) * sqrt(xprior_var))

                           reg_coef(ig) = covxy / yprior_var
                           ccc = reg_coef(ig) * (ypost_var / yprior_var - 1.0_r8) * yfmean
                           ddd = reg_coef(ig) * ypost_var / Robs

                           sublocx(ig)  = xt_xfmean
                           sublocy(ig)  = ccc + ddd * ytrue
                           sublocx2(ig) = xt_xfmean**2
                           sublocy2(ig) = (ccc + ddd*ytrue)**2
 
                           sublocnumer(ig) = xt_xfmean * (ccc + ddd * ytrue)
                           sublocdenom(ig) = (ccc + ddd*ytrue)**2 + ddd**2 * Robs

!print *, ' xy   info', xtrue, xfmean, xt_xfmean

                        enddo   ! ig

!print *, 'stat info ', corr(:), reg_coef(:)
!print *, 'stat info ', sublocx(:), sublocy(:)
!print *, 'stat info ', sublocx2(:), sublocy2(:)
!print *, 'stat info ', sublocnumer(:), sublocdenom(:)
!pause

                        ! storage the computations
                        ! discard the pair that has very small denominator for at least one group
                        indk1 = 1
                        do igg = 1, group_size
                           if ( abs(sublocdenom(igg)) < 0.0000000001_r8 ) then
                               indk1 = -1
                           endif
                        enddo

                        if ( indk1 > 0 ) then
                           localization%num(ilocind,ndims,inddim) =                                &
                               localization%num(ilocind,ndims,inddim) + 1

                           sumbeta1 = (sum(reg_coef))**2
                           sumbeta2 = 0.0_r8
                           do igg = 1, group_size
                              sumbeta2 = sumbeta2 + reg_coef(igg)**2
                           enddo

                           do igg = 1, group_size
                              sublocalpha(igg) = sublocnumer(igg)/sublocdenom(igg)
                           enddo

                           localization%beta1(ilocind,ndims,inddim) =                              &
                               localization%beta1(ilocind,ndims,inddim) + sumbeta1
                           localization%beta2(ilocind,ndims,inddim) =                              &
                               localization%beta2(ilocind,ndims,inddim) + sumbeta2

                           do ig = 1, group_size
                           localization%sumx(ilocind,ig,ndims,inddim) =                            &
                               localization%sumx(ilocind,ig,ndims,inddim) + sublocx(ig)
                           localization%sumy(ilocind,ig,ndims,inddim) =                            &
                               localization%sumy(ilocind,ig,ndims,inddim) + sublocy(ig)

                           localization%sumx2(ilocind,ig,ndims,inddim) =                           &
                               localization%sumx2(ilocind,ig,ndims,inddim) + sublocx2(ig)
                           localization%sumy2(ilocind,ig,ndims,inddim) =                           &
                               localization%sumy2(ilocind,ig,ndims,inddim) + sublocy2(ig)

                           localization%numerator(ilocind,ig,ndims,inddim) =                       &
                               localization%numerator(ilocind,ig,ndims,inddim) + sublocnumer(ig)
                           localization%denominator(ilocind,ig,ndims,inddim) =                     &
                               localization%denominator(ilocind,ig,ndims,inddim) + sublocdenom(ig)

                           localization%alpha1(ilocind,ig,ndims,inddim) =                          &
                               localization%alpha1(ilocind,ig,ndims,inddim) + sublocalpha(ig)
                           localization%alpha2(ilocind,ig,ndims,inddim) =                          &
                               localization%alpha2(ilocind,ig,ndims,inddim) + sublocalpha(ig)**2

                           localization%correlation(ilocind,ig,ndims,inddim) =                     &
                               localization%correlation(ilocind,ig,ndims,inddim) + abs(corr(ig))
                           enddo

!print *, 'loc info ', localization%sumx(ilocind,:,ndims,inddim)
!pause

                           ! record terrain less than 100m for both obs and var points
                           if ( obshgt <= 100.0 .and. varhgt <= 100.0 ) then
                              indk2 = obslev_num*varlev_num+ndims

                              localization%num(ilocind,indk2,inddim) =                                &
                                  localization%num(ilocind,indk2,inddim) + 1

                              localization%beta1(ilocind,indk2,inddim) =                              &
                                  localization%beta1(ilocind,indk2,inddim) + sumbeta1
                              localization%beta2(ilocind,indk2,inddim) =                              &
                                  localization%beta2(ilocind,indk2,inddim) + sumbeta2

                              do ig = 1, group_size   
                              localization%sumx(ilocind,ig,indk2,inddim) =                            &
                                  localization%sumx(ilocind,ig,indk2,inddim) + sublocx(ig)
                              localization%sumy(ilocind,ig,indk2,inddim) =                            &
                                  localization%sumy(ilocind,ig,indk2,inddim) + sublocy(ig)
   
                              localization%sumx2(ilocind,ig,indk2,inddim) =                           &
                                  localization%sumx2(ilocind,ig,indk2,inddim) + sublocx2(ig)
                              localization%sumy2(ilocind,ig,indk2,inddim) =                           &
                                  localization%sumy2(ilocind,ig,indk2,inddim) + sublocy2(ig)
   
                              localization%numerator(ilocind,ig,indk2,inddim) =                       &
                                  localization%numerator(ilocind,ig,indk2,inddim) + sublocnumer(ig)
                              localization%denominator(ilocind,ig,indk2,inddim) =                     &
                                  localization%denominator(ilocind,ig,indk2,inddim) + sublocdenom(ig)
   
                              localization%alpha1(ilocind,ig,indk2,inddim) =                          &
                                  localization%alpha1(ilocind,ig,indk2,inddim) + sublocalpha(ig)
                              localization%alpha2(ilocind,ig,indk2,inddim) =                          &
                                  localization%alpha2(ilocind,ig,indk2,inddim) + sublocalpha(ig)**2
   
                              localization%correlation(ilocind,ig,indk2,inddim) =                     &
                                  localization%correlation(ilocind,ig,indk2,inddim) + abs(corr(ig))
                              enddo
                           endif
                        endif

                     enddo   ! ii   
                  enddo      ! jj
               enddo         ! kk

            enddo   !i
         enddo      !j
      enddo         !k

      deallocate(yprior,xprior)

      deallocate(corr,reg_coef,sublocx,sublocy,sublocx2,sublocy2,sublocnumer,sublocdenom,sublocalpha)

      enddo LoopVarTypes
!$OMP END PARALLEL DO 

! Lili test
print *, 'computation completed'

! NOTE: the writting in OMP does not work...
!       thus put the writting out of OMP loop
      do inddim = indvar_s, indvar_e
         writelenloc = locnum
         write(idchar, '(i1.1)') id

         do obslevk = 1, obslev_num
            k = obslev_ind(obslevk)
            if ( idobs == 1 ) then
                 k = 1
                 print *, 'OBS is ', obs%dom(id)%obstype_metadata(idobs)
                 print *, 'NOTE: change obs lev index ', obslev_ind(obslevk), k
            endif

            do varlevk = 1, varlev_num
               kk = varlev_ind(obslevk,varlevk)
               if ( inddim == 12 .and. kk == 1 ) then
                    kk = 2
                    print *, 'VAR is ', trim(wrf%dom(id)%description(inddim))
                    print *, 'NOTE: change var lev index ', varlev_ind(obslevk,varlevk), kk
               elseif ( inddim < 8 ) then
                    kk = 1
                    print *, 'VAR is ', trim(wrf%dom(id)%description(inddim))
                    print *, 'NOTE: change var lev index ', varlev_ind(obslevk,varlevk), kk
               endif

               ndims = (obslevk-1)*varlev_num + varlevk

               if ( k < 10 ) then
                  write(levchar1,'(i1)') k
               else
                  write(levchar1,'(i2)') k
               endif
               if ( kk < 10 ) then
                  write(levchar2,'(i1)') kk
               else
                  write(levchar2,'(i2)') kk
               endif

               ! write pairs in whole domain
               var_name_tmp = trim(obs%dom(id)%obstype_metadata(idobs))//'_'//         &
                              trim(wrf%dom(id)%description(inddim))//'_d0'//idchar     &
                              //'_obslev'//trim(levchar1)//'_varlev'//trim(levchar2)

               iunit = 10 + inddim
               open(iunit, FILE=var_name_tmp, FORM='FORMATTED')

               allocate(numtmp(writelenloc))
               numtmp(1:writelenloc) = localization%num(1:writelenloc,ndims,inddim)
               write(iunit, FMT=locreal_format) numtmp 
               write(iunit, FMT=locreal_format) localization%beta1(1:writelenloc,ndims,inddim)
               write(iunit, FMT=locreal_format) localization%beta2(1:writelenloc,ndims,inddim)

               do ig = 1, group_size
                  write(iunit, FMT=locreal_format) localization%sumx(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumx2(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy2(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%numerator(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%denominator(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%alpha1(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%alpha2(1:writelenloc,ig,ndims,inddim)
                  write(iunit, FMT=locreal_format) localization%correlation(1:writelenloc,ig,ndims,inddim)
               enddo
               deallocate(numtmp)
               close(iunit)

               ! write pairs in ocean
               var_name_tmp = trim(obs%dom(id)%obstype_metadata(idobs))//'_'//               &
                              trim(wrf%dom(id)%description(inddim))//'_ocean_d0'//idchar     &
                              //'_obslev'//trim(levchar1)//'_varlev'//trim(levchar2)

               iunit = 10 + inddim
               open(iunit, FILE=var_name_tmp, FORM='FORMATTED')
               indk2 = obslev_num*varlev_num+ndims

               allocate(numtmp(writelenloc))
               numtmp(1:writelenloc) = localization%num(1:writelenloc,indk2,inddim)
               write(iunit, FMT=locreal_format) numtmp
               write(iunit, FMT=locreal_format) localization%beta1(1:writelenloc,indk2,inddim)
               write(iunit, FMT=locreal_format) localization%beta2(1:writelenloc,indk2,inddim)

               do ig = 1, group_size
                  write(iunit, FMT=locreal_format) localization%sumx(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%sumx2(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%sumy2(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%numerator(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%denominator(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%alpha1(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%alpha2(1:writelenloc,ig,indk2,inddim)
                  write(iunit, FMT=locreal_format) localization%correlation(1:writelenloc,ig,indk2,inddim)
               enddo
               deallocate(numtmp)
               close(iunit)

            enddo
         enddo
      enddo

      deallocate(localization%num, localization%numerator, localization%denominator)
      deallocate(localization%sumx, localization%sumy, localization%sumx2, localization%sumy2)
      deallocate(localization%alpha1, localization%alpha2, localization%correlation)

   enddo LoopObstypes

enddo LoopDomain


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


!subroutine draw_perfect_obs_altimeter(we,sn,nbdyoff,ens_size,description,wrf_var_2d,wrf_var_3d,obs_var_3d)
!
!integer,            intent(in)  :: we, sn, nbdyoff, ens_size
!character(len=129), intent(in)  :: description
!real(r8),           intent(in)  :: wrf_var_2d(:,:)      ! hgt
!real(r8),           intent(in)  :: wrf_var_3d(:,:,:)    ! psfc
!real(r8),           intent(out) :: obs_var_3d(:,:,:)
!! real(r8),           intent(in)  :: Robs
!
!integer    :: i, j, k, ind, obscopysize
!real(r8)   :: pres, hsfc, altr, var1, oerr
!real(r8), allocatable :: temp_var_1d(:)
!
!allocate(temp_var_1d(ens_size))
!
!obscopysize = ens_size + 4
!
!obs_var_3d(1:we,1:sn,1:obscopysize) = missing_r8
!
!if ( (trim(description) == "SFC_ALTIMETER") ) then 
!
!    do j = nbdyoff+1, sn-nbdyoff
!       do i = nbdyoff+1, we-nbdyoff
!          ! compute the truth
!          pres = wrf_var_3d(i,j,1)
!          hsfc = wrf_var_2d(i,j)
!          altr = compute_altimeter(pres*0.01_r8,hsfc)          
!          if ( (altr >= 880.0_r8) .and. (altr <1100.0_r8) )  then
!              obs_var_3d(i,j,1) = altr
!          else
!              return
!          endif
!
!          ! assign obs error variance
!          oerr = land_pres_error(pres)
!          obs_var_3d(i,j,obscopysize) = oerr * oerr
!
!          ! compute ensemble members
!          do ind = 1, ens_size
!             pres = wrf_var_3d(i,j,ind+3)
!             altr = compute_altimeter(pres*0.01_r8,hsfc)
!             if ( (altr >= 880.0_r8) .and. (altr <1100.0_r8) )  then
!                 obs_var_3d(i,j,ind+3) = altr
!             endif
!          enddo
!
!          ! compute ensemble mean
!          temp_var_1d(1:ens_size) = obs_var_3d(i,j,4:ens_size+3)
!          obs_var_3d(i,j,2) = sum(temp_var_1d) / ens_size
!
!          ! compute ensemble spread
!          var1 = 0.0_r8
!          do k = 1, ens_size
!             var1 = var1 + (temp_var_1d(k)-obs_var_3d(i,j,2))**2
!          enddo
!          var1 = var1 / (ens_size - 1)
!          obs_var_3d(i,j,3) = sqrt(var1)

!       enddo
!    enddo
!
!else
!    print *, 'ERROR: draw_perfect_obs_altimeter'
!    print *, 'ERROR: wrong obs type', trim(description)
!endif
!
!deallocate(temp_var_1d)
!
!return
!end subroutine draw_perfect_obs_altimeter


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

!if ( abs(dist) > max(abs(locind(1)),abs(locind(locnum))) ) then
!   return
!endif 

if ( dist <= locind(1) ) then
   indfind = 1
   return
else
   do ind = 2, locnum
      if ( (dist > locind(ind-1)) .and. (dist <= locind(ind)) ) then
         indfind = ind
         return
      endif
   enddo
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


end PROGRAM emp_localization



