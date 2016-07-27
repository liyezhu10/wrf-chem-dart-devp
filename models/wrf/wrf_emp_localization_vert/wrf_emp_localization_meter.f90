
PROGRAM emp_localization

use        types_mod, only : r8, missing_r8, DEG2RAD, PI, gravity, ps0,              &
                             ps0, gas_constant, gas_constant_v
use    utilities_mod, only : find_namelist_in_file, check_namelist_read, nc_check
use        model_mod, only : static_init_model, compute_geometric_height, toGrid,    &
                             pres_to_zk
use         map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                              ij_to_latlon, gridwind_to_truewind
use  mpi_utilities_mod, only: sleep_seconds

use misc_definitions_module, only : PROJ_LC

use           netcdf

implicit none



integer           :: emp_loc_kind     = 2
integer           :: locnum           = 2000
integer           :: locnumvert       = 40     ! this is one way number of index points (total: 2*locnumvert+1)
real(r8)          :: locrad           = 1000     ! 200 meter instead of radiance of 0.0002 (~1km) 
integer           :: wrf_ndom         = 1
real(r8)          :: cutoff_sfcobs    = 1000.0 
real(r8)          :: cutoff_upprobs   = 2000.0
integer           :: nbdy_nouse       = 5
integer           :: num_state_vars   = 20
integer           :: ens_size         = 60

namelist /wrf_emp_localization_nml/emp_loc_kind, wrf_ndom,   &
         cutoff_sfcobs, cutoff_upprobs 




character(len=80)  :: truth_name_input  = 'truth.nc',        &
                      prior_name_input  = 'prior.nc',        &
                      posterior_name_input = 'posterior.nc', &
                      obs_name_input    = 'obs.nc'
!character(len=80)  :: staticvarname(100), statevarname(100), modelhgtname(100)



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

   integer :: number_of_wrf_variables
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
    integer    :: obsindex, obscopysize, obslocsize, obsqcsize, obstypesize
    integer    :: numobstype
    integer,  dimension(:), pointer   ::  obstype
    character(len=32), dimension(:), pointer :: obstype_metadata
    integer,  dimension(:,:), pointer :: obsqc
    integer,  dimension(:), pointer   :: which_vert
    real(r8), dimension(:,:), pointer :: location
    real(r8), dimension(:,:), pointer :: observations
    integer,  dimension(:), pointer   :: obstypeindex
end type obs_data

type(obs_data)   :: obs

type subobs_data
    integer    :: num
    integer,  dimension(:), pointer     :: obsind
    real(r8), dimension(:,:), pointer   :: obsloc
    integer,  dimension(:,:), pointer   :: obsijk
    real(r8), dimension(:,:), pointer   :: obsvp

end type subobs_data

type(subobs_data) :: subobs

type localization_data
    real(r8), dimension(:), pointer     :: locind
    integer,  dimension(:,:), pointer   :: num
    real(r8), dimension(:,:), pointer   :: nominator
    real(r8), dimension(:,:), pointer   :: denominator
end type localization_data

type(localization_data) localization

type localization_omp_data
  type(localization_data), pointer  :: obs(:)
  integer                           :: nobs
end type localization_omp_data

type(localization_omp_data) loc_for_omp


integer          :: i, j, k, id, idobs, ind, inddim, indens,           &
                    indk1, indk2, var_id, ndims, dimids(10)
integer          :: ncid, ncidtruth, ncidprior, ncidposterior, ncidobs
integer          :: proj_code
integer          :: ilocind, sublocnum, writelenloc
real(r8)         :: stdlon,truelat1,truelat2,latinc,loninc
real(r8)         :: obslat, obslon, obspres, obsi, obsj, obsk
real(r8)         :: yobs, yfmean, Robs, yvar, yinno_norm,              &
                    xfmean, xtrue, xt_xfmean, covxy, deltaxmean
real(r8)         :: sublocnomi, sublocdenomi
logical          :: is_lev0
character(len=1) :: idchar
character(len=NF90_MAX_NAME) :: var_name_tmp, name
character(len=100) :: obsnum_format, locint_format, locreal_format

real(r8), allocatable :: wrf_var_1d(:), wrf_var_2d(:,:), wrf_var_3d(:,:,:), wrf_var_4d(:,:,:,:)
real(r8), allocatable :: wrf_mu_d01(:,:), wrf_psfc_d01(:,:), wrf_ph_d01(:,:,:),       &
                         wrf_t_d01(:,:,:), wrf_qvapor_d01(:,:,:)
real(r8), allocatable :: xprior(:), yprior(:)
real(r8), allocatable :: xtrueprf(:), xfmeanprf(:), xpriorprf(:,:)

! define format to write "localization"
locint_format = ' ( 2000i10 ) '
locreal_format = ' ( 2000e15.7 ) '
! define format to write obs number
obsnum_format = ' ( i10 )'

wrf%model_size = wrf_ndom

allocate(wrf%dom(wrf_ndom))

!-------------------------------------------------------
! open obs file (obs_epoch.nc) 
!-------------------------------------------------------
   print *, 'starting reading obs file'
   call nc_check( nf90_open(obs_name_input, NF90_NOWRITE, ncidobs),       &
                  'wrf_emp_localization','obs file')

!-------------------------------------------------------
! read obs file dimensions
!-------------------------------------------------------
   call nc_check( nf90_inq_dimid(ncidobs, "ObsIndex", var_id),       &
                  'read_obs_dimensions','inq_dimid ObsIndex')
   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obsindex),  &
                  'read_obs_dimensions','inquire_dimension'//trim(name))

   call nc_check( nf90_inq_dimid(ncidobs, "copy", var_id),       &
                  'read_obs_dimensions','inq_dimid copy')
   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obscopysize),  &
                  'read_obs_dimensions','inquire_dimension'//trim(name))

   call nc_check( nf90_inq_dimid(ncidobs, "qc_copy", var_id),       &
                  'read_obs_dimensions','inq_dimid qc_copy')
   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obsqcsize),  &
                  'read_obs_dimensions','inquire_dimension'//trim(name))

   call nc_check( nf90_inq_dimid(ncidobs, "location", var_id),       &
                  'read_obs_dimensions','inq_dimid location')
   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obslocsize),  &
                  'read_obs_dimensions','inquire_dimension'//trim(name))

   call nc_check( nf90_inq_dimid(ncidobs, "ObsTypes", var_id),       &
                  'read_obs_dimensions','inq_dimid ObsTypes')
   call nc_check( nf90_inquire_dimension(ncidobs, var_id, name, obs%obstypesize),  &
                  'read_obs_dimensions','inquire_dimension'//trim(name))

!-------------------------------------------------------
! read obs file variable
!-------------------------------------------------------
   allocate(obs%obstype_metadata(obs%obstypesize))
   call nc_check( nf90_inq_varid(ncidobs, "ObsTypesMetaData", var_id), &
                     'read_obs_variable','inq_varid ObsTypesMetaData')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%obstype_metadata), &
                     'read_obs_variable','get_var ObsTypesMetaData')

   allocate(obs%obstype(obs%obsindex))
   call nc_check( nf90_inq_varid(ncidobs, "obs_type", var_id), &
                     'read_obs_variable','inq_varid obs_type')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%obstype), &
                     'read_obs_variable','get_var obs_type')

   ! NOTE: nf90_get_var has the dimension inversely after read in than in the original .nc file
   allocate(obs%obsqc(obs%obsqcsize,obs%obsindex))
   call nc_check( nf90_inq_varid(ncidobs, "qc", var_id), &
                     'read_obs_variable','inq_varid qc')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%obsqc), &
                     'read_obs_variable','get_var qc')

   allocate(obs%which_vert(obs%obsindex))
   call nc_check( nf90_inq_varid(ncidobs, "which_vert", var_id), &
                     'read_obs_variable','inq_varid which_vert')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%which_vert), &
                     'read_obs_variable','get_var which_vert')

   allocate(obs%location(obs%obslocsize,obs%obsindex))
   call nc_check( nf90_inq_varid(ncidobs, "location", var_id), &
                     'read_obs_variable','inq_varid location')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%location), &
                     'read_obs_variable','get_var location')

   allocate(obs%observations(obs%obscopysize,obs%obsindex))
   call nc_check( nf90_inq_varid(ncidobs, "observations", var_id), &
                     'read_obs_variable','inq_varid observations')  ! reuse var_id, no harm
   call nc_check( nf90_get_var(ncidobs, var_id, obs%observations), &
                     'read_obs_variable','get_var observations')


!-------------------------------------------------------
! open truth file (True_State.nc) 
!-------------------------------------------------------
print *, 'starting reading truth and prior file'
call nc_check( nf90_open(truth_name_input, NF90_NOWRITE, ncidtruth),   &
               'wrf_emp_localization','open truth file')

!-------------------------------------------------------
! open prior file (Prior_Diag.nc) 
!-------------------------------------------------------
call nc_check( nf90_open(prior_name_input, NF90_NOWRITE, ncidprior),   &
               'wrf_emp_localization','open prior file')

!-------------------------------------------------------
! open prior file (Prior_Diag.nc) 
!-------------------------------------------------------
call nc_check( nf90_open(posterior_name_input, NF90_NOWRITE, ncidposterior),   &
               'wrf_emp_localization','open prior file')


WRFDomains : do id = 1, wrf_ndom 

!-------------------------------------------------------
! fill in state variable information 
!-------------------------------------------------------
   wrf%dom(id)%number_of_wrf_variables = num_state_vars
   allocate(wrf%dom(id)%description(wrf%dom(id)%number_of_wrf_variables))
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
   wrf%dom(id)%description(15) = 'QRAIN'
   wrf%dom(id)%description(16) = 'QNRAIN'
   wrf%dom(id)%description(17) = 'QICE'
   wrf%dom(id)%description(18) = 'QNICE'
   wrf%dom(id)%description(19) = 'QSNOW'
   wrf%dom(id)%description(20) = 'QGRAUP'

!   print *, 'state_var_discription', trim(wrf%dom(id)%description(1)), trim(wrf%dom(id)%description(20))

   write(idchar,'(i1.1)') id

!-------------------------------------------------------
! read WRF dimensions
!-------------------------------------------------------
   ! dimension of copies (1(truth)+2(ens mean, ens spread)+ens_size+2(inf mean, inf spread))
   wrf%dom(id)%ncopy = 1 +  ens_size + 4

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

   var_name_tmp = 'XLONG_U_d0'//idchar
   allocate(wrf%dom(id)%longitude_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLONG_U')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%longitude_u), &
                     'read_wrf_static_data','get_var XLONG_U')

   var_name_tmp = 'XLONG_V_d0'//idchar
   allocate(wrf%dom(id)%longitude_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id), &
                     'read_wrf_static_data','inq_varid XLONG_V')
   call nc_check( nf90_get_var(ncidtruth, var_id, wrf%dom(id)%longitude_v), &
                     'read_wrf_static_data','get_var XLONG_V')

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
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))
   allocate(wrf%dom(id)%variables(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,    &
                                  wrf%dom(id)%ncopy,wrf%dom(id)%number_of_wrf_variables))

   do ind = 1,wrf%dom(id)%number_of_wrf_variables

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
!-------------------------------------------------------
   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      var_name_tmp = trim(wrf%dom(id)%description(ind))//'_d0'//idchar

      call nc_check( nf90_inq_varid(ncidprior, var_name_tmp, var_id),  &
                     'get_variable_name_from_file',                    &
                     'inq_varid '//var_name_tmp)

      if ( wrf%dom(id)%var_size(3,ind) == 1 ) then

         ! Get 2D variable 
         ! NOTE: nf90_get_var automatically ignor the dimension with length 1
         allocate(wrf_var_3d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind),                  &
                             wrf%dom(id)%ncopy-1))
         call nc_check( nf90_get_var(ncidprior, var_id, wrf_var_3d),       &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1,                                          &
                               2:wrf%dom(id)%ncopy,                        &
                               ind) =                                      &
                    wrf_var_3d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%ncopy-1)
         deallocate(wrf_var_3d)

      else

         ! Get 3D variable
         allocate(wrf_var_4d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind),                  &
                             wrf%dom(id)%var_size(3,ind),                  &
                             wrf%dom(id)%ncopy-1))
         call nc_check( nf90_get_var(ncidprior, var_id, wrf_var_4d),       &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind),              &
                               2:wrf%dom(id)%ncopy,                        &
                               ind) =                                      &
                    wrf_var_4d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind),              &
                               1:wrf%dom(id)%ncopy-1)

         deallocate(wrf_var_4d)

      endif

   enddo


!-------------------------------------------------------
! read WRF variables of Posterior_Diag.nc  
! using posterior mean to replace truth copy (1)
!-------------------------------------------------------
   do ind = 1,wrf%dom(id)%number_of_wrf_variables

      var_name_tmp = trim(wrf%dom(id)%description(ind))//'_d0'//idchar

      call nc_check( nf90_inq_varid(ncidposterior, var_name_tmp, var_id),  &
                     'get_variable_name_from_file',                    &
                     'inq_varid '//var_name_tmp)

      if ( wrf%dom(id)%var_size(3,ind) == 1 ) then

         ! Get 2D variable 
         ! NOTE: nf90_get_var automatically ignor the dimension with length 1
         allocate(wrf_var_3d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind),                  &
                             wrf%dom(id)%ncopy-1))
         call nc_check( nf90_get_var(ncidposterior, var_id, wrf_var_3d),   &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1,                                          &
                               1,                                          &
                               ind) =                                      &
                    wrf_var_3d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1)
         deallocate(wrf_var_3d)

      else

         ! Get 3D variable
         allocate(wrf_var_4d(wrf%dom(id)%var_size(1,ind),                  &
                             wrf%dom(id)%var_size(2,ind),                  &
                             wrf%dom(id)%var_size(3,ind),                  &
                             wrf%dom(id)%ncopy-1))
         call nc_check( nf90_get_var(ncidposterior, var_id, wrf_var_4d),   &
                     'read_in_wrf_var','get_var '//var_name_tmp )
         wrf%dom(id)%variables(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind),              &
                               1,                                          &
                               ind) =                                      &
                    wrf_var_4d(1:wrf%dom(id)%var_size(1,ind),              &
                               1:wrf%dom(id)%var_size(2,ind),              &
                               1:wrf%dom(id)%var_size(3,ind),              &
                               1)

         deallocate(wrf_var_4d)

      endif

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
! 1. sigh the model height index to the number_of_wrf_varaibles
  wrf%dom(id)%nmhgt = 4
  allocate(wrf%dom(id)%mhgtind(wrf%dom(id)%number_of_wrf_variables))

  do ind = 1,wrf%dom(id)%number_of_wrf_variables
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

  deallocate(wrf_var_3d)

enddo WRFDomains

!-------------------------------------------------------
! find how many obs types in obs_epoch file 
!-------------------------------------------------------
print *, 'finding the number of obs types'
allocate(obs%obstypeindex(100))

call count_numobstype(obs%obsindex, obs%obstype, obs%obsqc, obs%numobstype, obs%obstypeindex)

!print *, 'num obs type is ', obs%numobstype
!print *, 'obs type index', obs%obstypeindex(1:obs%numobstype)
!print *, 'obs type meta', obs%obstype_metadata(obs%obstypeindex(1))

!-------------------------------------------------------
! Big loop: for each obs type  
!-------------------------------------------------------
! allocate the variables needed when computing model_pressure_profile
   allocate(wrf_mu_d01(wrf%dom(1)%we, wrf%dom(1)%sn))
   allocate(wrf_psfc_d01(wrf%dom(1)%we, wrf%dom(1)%sn))
   allocate(wrf_ph_d01(wrf%dom(1)%we, wrf%dom(1)%sn, wrf%dom(1)%bts))
   allocate(wrf_t_d01(wrf%dom(1)%we, wrf%dom(1)%sn, wrf%dom(1)%bt))
   allocate(wrf_qvapor_d01(wrf%dom(1)%we, wrf%dom(1)%sn, wrf%dom(1)%bt))

! NOTE: give truth MU, PH, T, QVAPOR, PSFC to compute model_pressure_file 
!   (truth copy - 1; MU index - 6, PH index - 12, T index - 11, QVAPOR index - 13, PSFC index - 7)  
   wrf_mu_d01(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn) =                          &
      wrf%dom(1)%variables(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1, 1, 6)  
   wrf_psfc_d01(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn) =                        &
      wrf%dom(1)%variables(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1, 1, 7) 
   wrf_ph_d01(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bts) =        &
      wrf%dom(1)%variables(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bts, 1, 12) 
   wrf_t_d01(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bt) =          &
      wrf%dom(1)%variables(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bt, 1, 11)
   wrf_qvapor_d01(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bt) =     &
      wrf%dom(1)%variables(1:wrf%dom(1)%we, 1:wrf%dom(1)%sn, 1:wrf%dom(1)%bt, 1, 13)

! before big loop for each obs type, initialize locind which contains the radiance for each subset
   if ( emp_loc_kind == 1 ) then
      print *, 'empirical localization kind: horizontal'
      allocate(localization%locind(locnum))
      do ind = 1, locnum
         localization%locind(ind) = locrad * real(ind)
      enddo
      !print *, 'loc index ', localization%locind(1:locnum)
   elseif ( emp_loc_kind == 2 ) then
      print *, 'empirical localization kind: vertical'
      allocate(localization%locind(2*locnumvert+1))
      do ind = 1, locnumvert 
         localization%locind(ind) = -1.0 * locrad * real(locnumvert - ind + 1.0)
      enddo
      localization%locind(locnumvert+1) = 0.0
      do ind = 1, locnumvert
         localization%locind(locnumvert+1+ind) = locrad * real(ind)
      enddo
      print *, 'loc index ', localization%locind(1:2*locnumvert+1)

!pause

   elseif ( emp_loc_kind == 3 ) then
      print *, 'empirical localization kind: 3D'
   else
      print *, 'ERROR: not support empirical localization kind', emp_loc_kind
   endif

! Lili test
LoopObsTypes: do idobs = 1, obs%numobstype
!LoopObsTypes: do idobs = 16, 16

! Lili test
if ( trim(obs%obstype_metadata(obs%obstypeindex(idobs))) == "SFC_ALTIMETER" ) then

print *, 'Now processing obs type ', trim(obs%obstype_metadata(obs%obstypeindex(idobs)))

! 1. find the number of observations given obs type
!    set the index to subobs
   allocate(subobs%obsind(obs%obsindex))
   subobs%num = 0
   do ind = 1, obs%obsindex
      if ( obs%obstype(ind) == obs%obstypeindex(idobs) ) then
         subobs%num = subobs%num + 1
         subobs%obsind(subobs%num) = ind
      endif
   enddo    

   print *, 'subobs ', subobs%num
!   print *, 'subobs ', subobs%obsind(1), subobs%obsind(subobs%num)

! 2. loop through subobs 
   ! allocate 'subobs' variable
   allocate(subobs%obsloc(3,subobs%num))
   allocate(subobs%obsijk(3,subobs%num))
   allocate(subobs%obsvp(wrf%dom(1)%bt,subobs%num))

   ! allocate and initialize 'localization' variable
   ndims = 0
   do id = 1, wrf_ndom
      ndims = ndims + wrf%dom(id)%number_of_wrf_variables
   enddo
   if ( emp_loc_kind == 1 ) then

   elseif ( emp_loc_kind == 2 ) then
      allocate(localization%num(2*locnumvert+1,ndims))
      allocate(localization%nominator(2*locnumvert+1,ndims))
      allocate(localization%denominator(2*locnumvert+1,ndims))
      localization%num(1:2*locnumvert+1, 1:ndims)         = 0
      localization%nominator(1:2*locnumvert+1, 1:ndims)   = 0.0
      localization%denominator(1:2*locnumvert+1, 1:ndims) = 0.0

! Lili test: for omp
      loc_for_omp%nobs = subobs%num
      allocate(loc_for_omp%obs(subobs%num))
      do id = 1, loc_for_omp%nobs
         allocate(loc_for_omp%obs(id)%num(2*locnumvert+1,ndims))
         allocate(loc_for_omp%obs(id)%nominator(2*locnumvert+1,ndims))
         allocate(loc_for_omp%obs(id)%denominator(2*locnumvert+1,ndims))
         loc_for_omp%obs(id)%num(1:2*locnumvert+1,1:ndims) = 0
         loc_for_omp%obs(id)%nominator(1:2*locnumvert+1,1:ndims) = 0.0
         loc_for_omp%obs(id)%denominator(1:2*locnumvert+1,1:ndims) = 0.0
      enddo

   elseif ( emp_loc_kind == 3 ) then

   else
      print *, 'ERROR: not supprot empirical localization kind ', emp_loc_kind
   endif


! Lili test   LoopSubObs: do ind = 1, subobs%num
!$OMP PARALLEL DO PRIVATE(ind,id,i,j,k,inddim,ndims,indk1,indk2,obsi,obsj,obsk,obspres,is_lev0,yobs,yfmean,Robs,yvar,yinno_norm,xtrue,xfmean,xt_xfmean,covxy,deltaxmean,ilocind,sublocnum,sublocnomi,sublocdenomi,wrf_var_1d,xprior,yprior,xtrueprf,xfmeanprf,xpriorprf)
   LoopSubObs: do ind = 1, subobs%num 

   print *, 'processing obs number ', ind

   ! 2.1. given one obs, compute its height, model zk
   !      NOTE: this is done on D01
      if ( obs%which_vert(subobs%obsind(ind)) == -1 ) then
!        ! NOTE: locations from obs_epoch has index of (lon, lat, vert)
         subobs%obsloc(1,ind) = obs%location(1,subobs%obsind(ind))
!         ! NOTE: wrf has lon [-180, 180], to be consistent with wrf when compute the distance,
!         !       convert lon from obs in [0, 360] to [-180,180]
!         !       latlon_to_ij actually can handle lon both in [0, 360] and [-180, 180]
!         ! (NOTE: DART has lon [0, 360])
         if ( subobs%obsloc(1,ind) > 180.0 ) then
            subobs%obsloc(1,ind) = subobs%obsloc(1,ind) - 360.0
         endif
         subobs%obsloc(2,ind) = obs%location(2,subobs%obsind(ind))
         subobs%obsloc(3,ind) = obs%location(3,subobs%obsind(ind))
         call latlon_to_ij(wrf%dom(1)%proj, subobs%obsloc(2,ind), subobs%obsloc(1,ind), obsi, obsj)
         subobs%obsijk(1,ind) = obsi
         subobs%obsijk(2,ind) = obsj

!! Lili test: for omp
!print *, 'obs loc ', obs%location(1:3,subobs%obsind(ind))
!print *, 'subobs loc ', subobs%obsloc(1:3,ind)
!print *, 'subobs ijk ', subobs%obsijk(1:3,ind)

! NOTE: following uppr obs code has not been tested!
      elseif ( obs%which_vert(subobs%obsind(ind)) == 2 ) then
!        ! NOTE: locations from obs_epoch has index of (lon, lat, vert)
         subobs%obsloc(1,ind) = obs%location(1,subobs%obsind(ind)) 
!         ! NOTE: wrf has lon [-180, 180], to be consistent with wrf when compute the distance,
!         !       convert lon from obs in [0, 360] to [-180,180]
!         !       latlon_to_ij actually can handle lon both in [0, 360] and [-180, 180]
!         ! (NOTE: DART has lon [0, 360])
         if ( subobs%obsloc(1,ind) > 180.0 ) then
            subobs%obsloc(1,ind) = subobs%obsloc(1,ind) - 360.0
         endif
         subobs%obsloc(2,ind) = obs%location(2,subobs%obsind(ind))
         obspres = obs%location(3,subobs%obsind(ind))

         call latlon_to_ij(wrf%dom(1)%proj, subobs%obsloc(2,ind), subobs%obsloc(1,ind), obsi, obsj)          
         subobs%obsijk(1,ind) = obsi
         subobs%obsijk(2,ind) = obsj
         !print *, 'lat lon i j ', subobs%obsloc(2,ind), subobs%obsloc(1,ind), obsi, obsj

         ! NOTE: give truth MU, PH, T, QVAPOR, PSFC to compute the model pressure file 
         ! NOTE: vp(0) is computed from PSFC, vp(1:bt) is on model full level (mass)
         allocate(wrf_var_1d(0:wrf%dom(1)%bt))
         call get_obs_model_pressure_profile(wrf%dom(1)%bt, wrf_var_1d, obsi, obsj,              &
                 wrf%dom(1)%dnw, wrf%dom(1)%mub, wrf_mu_d01, wrf%dom(1)%phb, wrf_ph_d01,         &
                 wrf_t_d01, wrf_qvapor_d01, wrf_psfc_d01) 
         !print *, wrf_var_1d(0), wrf_var_1d(1), wrf_var_1d(wrf%dom(1)%bt)

         call pres_to_zk(obspres, wrf_var_1d, wrf%dom(1)%bt, obsk, is_lev0)
         if ( is_lev0 ) then
             ! NOTE: is_lev0 = .true. means obs above surface but below lowest sigma level
             !   since 'allow_obs_below' is .false. in input.nml, set obs vertical location as missing
             subobs%obsloc(3,ind) = missing_r8
         else
             ! NOTE: find the model height of this observation by using 3d model height on mass point (full level)
             subobs%obsloc(3,ind) = grid_mhgt_from_vp(obsi, obsj, obsk, wrf%dom(1)%mhgt(:,:,:,4))
         endif
         !print *, 'obsk, is_lev0 ', obsk, is_lev0
         !print *, 'obs height', wrf_psfc_d01(obsi,obsj), obspres, obsk, subobs%obsloc(3,ind)

         deallocate(wrf_var_1d)
      else
         print *, 'ERROR: do not support obs vertical type ', obs%which_vert(subobs%obsind(ind))
      endif

   ! 2.2. loop through domains 
      if ( subobs%obsloc(3,ind) /= missing_r8 ) then 

         ! get part of deltaxbar: yinno_norm = ( (yobs - yfmean) / (var(y) + Robs) )
         yobs = obs%observations(1,subobs%obsind(ind))      ! observation is copy 1
         yfmean = obs%observations(3,subobs%obsind(ind))    ! prior mean is copy 3
         Robs = obs%observations(obs%obscopysize,subobs%obsind(ind))  ! obs error variance is the last copy
         yvar = obs%observations(5,subobs%obsind(ind))**2   ! obs prior spread is copy 5
         yinno_norm = (yobs - yfmean) / (yvar + Robs)
!         print *, 'obs test ', yobs, yfmean, Robs, yvar, yinno_norm

         ! get priors of y to compute the covariance with state variable x later
         allocate(yprior(ens_size))
         yprior(1:ens_size) = obs%observations(7:obs%obsindex-1:2,subobs%obsind(ind)) ! prior copies
!         print *, 'yprior ', yprior(1), yprior(20), yprior(ens_size)

         ndims = 0       ! ndims now is the index for state-variable in 'localization'
         LoopWRFDomains: do id = 1, wrf_ndom

            if ( emp_loc_kind == 1 ) then        ! horizontal empirical localization


            elseif ( emp_loc_kind == 2 ) then           ! vertical empirical localization

               allocate(xtrueprf(wrf%dom(id)%bt))
               allocate(xfmeanprf(wrf%dom(id)%bt))
               allocate(xpriorprf(wrf%dom(id)%bt, ens_size))
               allocate(wrf_var_1d(wrf%dom(id)%bt))

! Lili test
!               do inddim = 1, wrf%dom(id)%number_of_wrf_variables
               do inddim = 12, 12

                  ! verticala empirical localization is only computed for 3D state variables
                  if ( wrf%dom(id)%var_size(3,inddim) > 1 ) then 
!print *, 'wrf var ', trim(wrf%dom(id)%description(inddim))
                     call get_model_variable_profile(obsi,obsj,wrf%dom(id)%bt,                      &
                                                     wrf%dom(id)%variables(:,:,:,:,inddim),         &
                                                     wrf%dom(id)%description(inddim), ens_size,     &
                                                     xtrueprf, xfmeanprf, xpriorprf)

                     call get_model_height_profile(obsi,obsj,wrf%dom(id)%bt,                             &
                                                   wrf%dom(id)%mhgt(:,:,:,wrf%dom(id)%mhgtind(inddim)),  &
                                                   wrf_var_1d)

                     do k = 1, wrf%dom(id)%bt
                        ! get xt_xfmean = xtrue - xfmean
                        xt_xfmean = xtrueprf(k) - xfmeanprf(k)
                        covxy = comp_cov(ens_size, xfmeanprf(k), xpriorprf(k,:), yfmean, yprior)
                        ! compute deltaxmean = cov(x,y) * (yobs - yfmean) / (var(y) + Robs)
                        !                    = covxy * yinno_norm
                        deltaxmean = covxy * yinno_norm
                        call comp_emp_loc_vert(ilocind,sublocnum,sublocnomi,sublocdenomi,locnumvert,localization%locind,  &
                                               subobs%obsloc(2,ind),subobs%obsloc(3,ind),wrf_var_1d(k),                   &
                                               xt_xfmean, deltaxmean)

                        !print *, 'loc ', ilocind, sublocnomi, sublocdenomi

                        ! record in model level
                           loc_for_omp%obs(ind)%num(k,ndims+inddim)         =    &
                               loc_for_omp%obs(ind)%num(k,ndims+inddim) + sublocnum
                           loc_for_omp%obs(ind)%nominator(k,ndims+inddim)   =    &
                               loc_for_omp%obs(ind)%nominator(k,ndims+inddim) + sublocnomi
                           loc_for_omp%obs(ind)%denominator(k,ndims+inddim) =    &
                               loc_for_omp%obs(ind)%denominator(k,ndims+inddim) + sublocdenomi


!                        if ( ilocind /= -1 ) then
!
!                           loc_for_omp%obs(ind)%num(ilocind,ndims+inddim)         =    &
!                               loc_for_omp%obs(ind)%num(ilocind,ndims+inddim) + sublocnum
!                           loc_for_omp%obs(ind)%nominator(ilocind,ndims+inddim)   =    &
!                               loc_for_omp%obs(ind)%nominator(ilocind,ndims+inddim) + sublocnomi
!                           loc_for_omp%obs(ind)%denominator(ilocind,ndims+inddim) =    &
!                               loc_for_omp%obs(ind)%denominator(ilocind,ndims+inddim) + sublocdenomi
!
!                        endif

                     enddo

                  endif
               enddo

               deallocate(xtrueprf, xfmeanprf, xpriorprf)
               deallocate(wrf_var_1d)

            elseif ( emp_loc_kind == 3 ) then           ! 3D empirical localization

            else

               print *, 'ERROR in wrf_emp_localization'
               print *, 'ERROR: emp_loc_kind does not supported' 

            endif          ! endif choice of emp_loc_kind 

            ndims = ndims + wrf%dom(id)%number_of_wrf_variables

         enddo LoopWRFDomains

         deallocate(yprior)

      endif      ! endif subobs%obsloc(3,ind) /= missing_r8

   enddo LoopSubObs
!$OMP END PARALLEL DO 

! Lili test: for omp
   do id = 1, loc_for_omp%nobs
      localization%num = localization%num + loc_for_omp%obs(id)%num
      localization%nominator = localization%nominator + loc_for_omp%obs(id)%nominator
      localization%denominator = localization%denominator + loc_for_omp%obs(id)%denominator
   enddo

   ! 2.4. write "localizaton" out for each obs type and each state variable for each domain
   if ( emp_loc_kind == 1 ) then
      writelenloc = locnum
   elseif ( emp_loc_kind == 2 ) then
      writelenloc = 2*locnumvert+1
   elseif ( emp_loc_kind == 3 ) then


   else
      print *, 'ERROR: not support emp_loc_kind', emp_loc_kind
   endif


   ndims = 0
   do id = 1, wrf_ndom 
      write(idchar,'(i1.1)') id
! Lili test
!      do inddim = 1, wrf%dom(id)%number_of_wrf_variables
      do inddim = 12, 12
      var_name_tmp = trim(obs%obstype_metadata(obs%obstypeindex(idobs)))//'_'//      &
                     trim(wrf%dom(id)%description(inddim))//'_d0'//idchar//'.dat'
      ! print *, 'var_name_tmp ', var_name_tmp
      open(20, FILE=var_name_tmp, FORM='FORMATTED')
      write(20, FMT=locint_format) localization%num(1:writelenloc,ndims+inddim)
      write(20, FMT=locreal_format) localization%nominator(1:writelenloc,ndims+inddim)
      write(20, FMT=locreal_format) localization%denominator(1:writelenloc,ndims+inddim)
      close(20)
      enddo
      ndims = ndims + wrf%dom(id)%number_of_wrf_variables
   enddo

   ! 2.5. write obs number out for each obs type
   var_name_tmp = trim(obs%obstype_metadata(obs%obstypeindex(idobs)))//'_numbers.dat'
   open(20, FILE=var_name_tmp, FORM='FORMATTED')
   write(20, FMT=obsnum_format) subobs%num
   close(20)

!   print *, 'obsloc ', subobs%obsloc(1,1), subobs%obsloc(2,1), subobs%obsloc(3,1)
!   print *, 'obsloc ', subobs%obsloc(1,subobs%num), subobs%obsloc(2,subobs%num), subobs%obsloc(3,subobs%num)

   deallocate(subobs%obsloc)
   deallocate(subobs%obsijk)
   deallocate(subobs%obsvp)
   deallocate(subobs%obsind)
   deallocate(localization%num, localization%nominator, localization%denominator)

! Lili test: for omp
   deallocate(loc_for_omp%obs)

endif

enddo LoopObsTypes

   deallocate(wrf_mu_d01, wrf_psfc_d01, wrf_ph_d01, wrf_t_d01, wrf_qvapor_d01)



contains


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


subroutine get_model_variable_profile(obsi,obsj,bt,variable,description,nens,xtrueprf,xfmeanprf,xpriorprf)
! some state variables are at mass points
! some are not, these need interpolation

! given observation (i, j), find the model variable profile at (obsi, obsj)

real(r8),           intent(in)  :: obsi, obsj
integer,            intent(in)  :: bt
real(r8),           intent(in)  :: variable(:,:,:,:)
character(len=129), intent(in)  :: description
integer,            intent(in)  :: nens
real(r8),           intent(out) :: xtrueprf(:), xfmeanprf(:)
real(r8),           intent(out) :: xpriorprf(:,:)

integer    :: indi, indj, k, ind
real(r8)   :: obsinew, obsjnew, disti, disti2, distj, distj2
real(r8)   :: var1, var2, var3, var4, var_k1, var_k2

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



