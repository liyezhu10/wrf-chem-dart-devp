
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



integer           :: emp_loc_kind     = 1
integer           :: locnum           = 1600
real(r8)          :: locrad           = 0.002  ! 200 meter instead of radiance of 0.0002 (~1km) 
integer           :: wrf_ndom         = 1
integer           :: nbdy_nouse       = 5
integer           :: num_state_vars   = 7
integer           :: ens_size         = 80
integer           :: lat_s            = 96+33
integer           :: lat_e            = 96+75
integer           :: lev_s            = 25
integer           :: lev_e            = 25
integer           :: indobs_s         = 2
integer           :: indobs_e         = 4
integer           :: indvar_s         = 2
integer           :: indvar_e         = 5

real              :: R_P              = 40000.0_r8 ! obs error variance of Psfc
real              :: R_T              = 1.0_r8   ! obs error variance of Temperature
real              :: R_U              = 4.0_r8   ! obs error variance of Wind



character(len=80)  :: truth_name_input  = 'truth.nc',        &
                      prior_name_input  = 'prior.nc',        &
                      obs_name_input    = 'obs.nc'


type cam_data 
   integer  :: nlat, nlon, nslat, nslon, nlev, nilev         ! dimensions
   integer  :: ncopy                                         ! copies of truth and prior; nmhgt
   real(r8), dimension(:),  pointer        :: lat, lon, lev, slat, slon, ilev
   ! lat lon for each variable (US staggered on lat, VS staggered on lon), NOTE no variable staggered on ilev
   real(r8), dimension(:,:),  pointer      :: var_lat, var_lon 
   integer                                 :: number_of_state_variables
   integer, dimension(:,:), pointer        :: var_size
   character(len=129),dimension(:),pointer :: description
   real(r8), dimension(:,:,:,:,:), pointer :: variables 
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
    real(r8), dimension(:), pointer       :: locind
    integer,  dimension(:,:,:), pointer   :: num
    real(r8), dimension(:,:,:), pointer   :: sumx        ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:), pointer   :: sumy        ! used to get mean(y)
    real(r8), dimension(:,:,:), pointer   :: sumx2       ! used to get mean(x^2)
    real(r8), dimension(:,:,:), pointer   :: sumy2       ! used to get mean(y^2)
    real(r8), dimension(:,:,:), pointer   :: numerator   ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:), pointer   :: denominator ! i.e., used to get mean(y^2) (note: it contains Robs)
    real(r8), dimension(:,:,:), pointer   :: alpha1
    real(r8), dimension(:,:,:), pointer   :: alpha2
    real(r8), dimension(:,:,:), pointer   :: correlation
end type localization_data

type(localization_data) localization

type localization_omp_data
  type(localization_data), pointer  :: obs(:)
  integer                           :: nobs
end type localization_omp_data

type(localization_omp_data) loc_for_omp


integer          :: i, j, k, id, idobs, ind, inddim, indens,           &
                    ii, jj, kk, indk1, indk2, var_id, ndims, dimids(10)
integer          :: ncid, ncidtruth, ncidprior, ncidobs, iunit
integer          :: proj_code
integer          :: ilocind, sublocnum, writelenloc
real(r8)         :: stdlon,truelat1,truelat2,latinc,loninc
real(r8)         :: obslat, obslon, obspres, obsi, obsj, obsk, obspalt, obshgt
real(r8)         :: varlat, varlon
real(r8)         :: yobs, ytrue, yfmean, oerr, Robs, yprior_var, ypost_var, yvar, yinno_norm,              &
                    xfmean, xtrue, xprior_var,xt_xfmean, covxy, deltaxmean, corr, reg_coef, ccc, ddd, dist
real(r8)         :: obsdiag1, obsdiag2, obsdiag3, obsdiag4
real(r8)         :: sublocx, sublocy, sublocx2, sublocy2, sublocnumer, sublocdenom, sublocalpha
logical          :: is_lev0
character(len=1) :: idchar
character(len=NF90_MAX_NAME) :: var_name_tmp, name
character(len=100) :: obsnum_format, locint_format, locreal_format, mhgt_format

real(r8), allocatable :: cam_var_1d(:), cam_var_2d(:,:), cam_var_3d(:,:,:), cam_var_4d(:,:,:,:),           &
                                                         cam_var_3d_2(:,:,:)
real(r8), allocatable :: obs_var_1d(:), obs_var_2d(:,:), obs_var_3d(:,:,:), obs_var_4d(:,:,:,:)
real(r8), allocatable :: wrf_mu_d01(:,:), wrf_psfc_d01(:,:), wrf_ph_d01(:,:,:),       &
                         wrf_t_d01(:,:,:), wrf_qvapor_d01(:,:,:)
real(r8), allocatable :: xprior(:), yprior(:)
real(r8), allocatable :: xtrueprf(:), xfmeanprf(:), xfspdprf(:), xpriorprf(:,:)
real(r8), allocatable :: xtrue3d(:,:,:), xfmean3d(:,:,:), xfspd3d(:,:,:), xprior3d(:,:,:,:) 
real(r8), allocatable :: numtmp(:)

! define format to write "localization"
locint_format = ' ( 2000i10 ) '
locreal_format = ' ( 2000e15.7 ) '
! define format to write obs number
!obsnum_format = ' ( i10,4e15.7 )'
obsnum_format = ' ( i10 ) '
mhgt_format = ' ( 49e15.7 ) '


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
               'cam_emp_localization','open truth file')

!-------------------------------------------------------
! open prior file (Prior_Diag.nc) 
!-------------------------------------------------------
call nc_check( nf90_open(prior_name_input, NF90_NOWRITE, ncidprior),   &
               'cam_emp_localization','open prior file')


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
   ! dimension of copies (1(truth)+2(ens mean, ens spread)+ens_size)
   ! NOTE: (inflation is turned off in the OSSE)
   cam%ncopy = 1 + ens_size + 2

   var_name_tmp = 'lat'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid lat')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nlat),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'slat'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nslat')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nslat),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'lon'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid lon')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nlon),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'slon'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nslon')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nslon),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'lev'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nlev')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nlev),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

   var_name_tmp = 'ilev'
   call nc_check( nf90_inq_dimid(ncidtruth, trim(var_name_tmp), var_id),       &
                  'read_cam_dimensions','inq_dimid nilev')
   call nc_check( nf90_inquire_dimension(ncidtruth, var_id, name, cam%nilev),  &
                  'read_cam_dimensions','inquire_dimension'//trim(name))

!-------------------------------------------------------
! read CAM dimension data 
!-------------------------------------------------------
! 1. 1D array

   allocate(cam%lat(1:cam%nlat))
   var_name_tmp = 'lat'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lat')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%lat),  &
                     'read_cam_dimension_data','get_var lat')

   allocate(cam%slat(1:cam%nslat))
   var_name_tmp = 'slat'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid slat')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%slat),  &
                     'read_cam_dimension_data','get_var slat')

   allocate(cam%lon(1:cam%nlon))
   var_name_tmp = 'lon'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lon')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%lon),  &
                     'read_cam_dimension_data','get_var lon')

   allocate(cam%slon(1:cam%nslon))
   var_name_tmp = 'slon'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid slon')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%slon),  &
                     'read_cam_dimension_data','get_var slon')

   allocate(cam%lev(1:cam%nlev))
   var_name_tmp = 'lev'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid lev')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%lev),  &
                     'read_cam_dimension_data','get_var lev')

   allocate(cam%ilev(1:cam%nilev))
   var_name_tmp = 'ilev'
   call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'read_cam_dimension_data','inq_varid ilev')
   call nc_check( nf90_get_var(ncidtruth, var_id, cam%ilev),  &
                     'read_cam_dimension_data','get_var ilev')

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
   allocate(cam%var_lat(1:cam%nlat,1:cam%number_of_state_variables))
   allocate(cam%var_lon(1:cam%nlon,1:cam%number_of_state_variables))
   
   do ind = 1,cam%number_of_state_variables
!   do ind = 1, 1

      ! Get the dimension size 
      ! Once this is done for variable of True_State.nc, there is no need to do this for
      ! variables in Prior_Diag.nc, since they should be consistent
      var_name_tmp = trim(cam%description(ind))
      call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                     'get_variable_size_from_file',                    &
                     'inq_varid '//var_name_tmp)
      call nc_check( nf90_inquire_variable(ncidtruth, var_id,          &
                     ndims=ndims, dimids=dimids),                      &
                     'get_variable_size_from_file',                    &
                     'inquire_variable '//var_name_tmp)
      do inddim = 1,ndims-2
         call nc_check( nf90_inquire_dimension(ncidtruth, dimids(inddim),  &
                        len=cam%var_size(inddim,ind)),                     &
                        'get_variable_size_from_file',                     &
                        'inquire_dimension '//var_name_tmp)
      enddo

      ! Get the data
      if ( ndims < 5 ) then

         ! 2D variable, its vertical dimension is 1
         cam%var_size(ndims-1,ind) = 1

         ! Get 2D variable
         allocate(cam_var_2d(cam%var_size(1,ind),                       &
                             cam%var_size(2,ind)))                      
         call nc_check( nf90_get_var(ncidtruth, var_id, cam_var_2d),    &
                        'read_in_cam_var','get_var '//var_name_tmp )
         cam%variables(1:cam%var_size(1,ind),                           &
                       1:cam%var_size(2,ind),                           &
                       1,                                               &
                       1,                                               &
                       ind) =                                           &
                    cam_var_2d(1:cam%var_size(1,ind),                   &
                               1:cam%var_size(2,ind))
         deallocate(cam_var_2d)

      else

         ! Get 3D variable
         allocate(cam_var_3d(cam%var_size(1,ind),                       &
                             cam%var_size(2,ind),                       &
                             cam%var_size(3,ind)))
         call nc_check( nf90_get_var(ncidtruth, var_id, cam_var_3d),    &
                     'read_in_cam_var','get_var '//var_name_tmp )
         ! reshape the data:
         ! cam_var_3d(lev,lon,lat) -> cam%variables(lon,lat,lev)
         ! also reshape the lev: 
         ! cam has lev(1->30) top->bottom, here shape to bottom->top
         do i = 1, cam%var_size(1,ind)
            cam%variables(1:cam%var_size(2,ind),                        &
                          1:cam%var_size(3,ind),                        &
                          cam%var_size(1,ind)-i+1,                    &
                          1,                                            &
                          ind) =                                        &
                       cam_var_3d(i,                                    &
                                  1:cam%var_size(2,ind),                &
                                  1:cam%var_size(3,ind))
         enddo
         ! reshape var_size 
         j = cam%var_size(1,ind)
         cam%var_size(1,ind) = cam%var_size(2,ind)
         cam%var_size(2,ind) = cam%var_size(3,ind)
         cam%var_size(3,ind) = j
         deallocate(cam_var_3d)

!print *, 'var ', cam%variables(10,10,:,1,ind)
!pause

      endif

      ! get the lat, lon ready for each variable (US, VS are staggered)
      if ( var_name_tmp == 'US' ) then
         cam%var_lat(1:cam%nslat,ind) = cam%slat(1:cam%nslat)
         cam%var_lon(1:cam%nlon,ind)  = cam%lon(1:cam%nlon)
!print *, 'lat US = ', cam%var_lat(:,ind)
!pause
!print *, 'lon US = ', cam%var_lon(:,ind)
!pause

      elseif ( var_name_tmp == 'VS' ) then
         cam%var_lat(1:cam%nlat,ind)  = cam%lat(1:cam%nlat)
         cam%var_lon(1:cam%nslon,ind) = cam%slon(1:cam%nslon)
      else
         cam%var_lat(1:cam%nlat,ind)  = cam%lat(1:cam%nlat)
         cam%var_lon(1:cam%nlon,ind)  = cam%lon(1:cam%nlon)
      endif

   enddo


!-------------------------------------------------------
! read state variables of Prio_Diag.nc (copy ens_mean, ens_spd, ens_member) 
!-------------------------------------------------------
   do ind = 1,cam%number_of_state_variables
   
      var_name_tmp = trim(cam%description(ind))

      call nc_check( nf90_inq_varid(ncidprior, var_name_tmp, var_id),  &
                     'get_variable_name_from_file',                    &
                     'inq_varid '//var_name_tmp)

      if ( cam%var_size(3,ind) == 1 ) then

         ! Get 2D variable 
         ! NOTE: nf90_get_var automatically ignor the dimension with length 1
         allocate(cam_var_3d(cam%var_size(1,ind),                      &
                             cam%var_size(2,ind),                      &
                             cam%ncopy-1))
         call nc_check( nf90_get_var(ncidprior, var_id, cam_var_3d),   &
                     'read_in_cam_var','get_var '//var_name_tmp )
         cam%variables(1:cam%var_size(1,ind),                          &
                       1:cam%var_size(2,ind),                          &
                       1,                                              &
                       2:cam%ncopy,                                    &
                       ind) =                                          &
                    cam_var_3d(1:cam%var_size(1,ind),                  &
                               1:cam%var_size(2,ind),                  &
                               1:cam%ncopy-1)
         deallocate(cam_var_3d)
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,1,2,ind)
!pause
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,1,83,ind)
!pause

      else

         ! Get 3D variable
         ! NOTE: var_size was reshaped when read in True_State.nc,
         !       but nc_varget starts from lev. this is how to set cam_var_4d
         allocate(cam_var_4d(cam%var_size(3,ind),                      &
                             cam%var_size(1,ind),                      &
                             cam%var_size(2,ind),                      &
                             cam%ncopy-1))
         call nc_check( nf90_get_var(ncidprior, var_id, cam_var_4d),   &
                     'read_in_cam_var','get_var '//var_name_tmp )
         do i = 1, cam%var_size(3,ind)
            cam%variables(1:cam%var_size(1,ind),                       &
                          1:cam%var_size(2,ind),                       &
                          cam%var_size(3,ind)-i+1,                     &
                          2:cam%ncopy,                                 &
                          ind) =                                       &
                       cam_var_4d(i,                                   &
                                  1:cam%var_size(1,ind),               &
                                  1:cam%var_size(2,ind),               &
                                  1:cam%ncopy-1)
         enddo
         deallocate(cam_var_4d)

!print *, 'var ', cam%variables(10,10,:,2,ind)
!pause
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,1,2,ind)
!pause
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,1,83,ind)
!pause
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,10,2,ind)
!pause
!print *, 'var = ', cam%variables(1:cam%var_size(1,ind),1,10,83,ind)
!pause

      endif

   enddo

!print *, 'dimension ', cam%nlat, cam%nslat, cam%nlon, cam%nslon, cam%nlev, cam%nilev
!pause
!print *, 'lat = ', cam%lat
!pause
!print *, 'slat = ', cam%slat
!pause
!print *, 'lon = ', cam%lon
!pause
!print *, 'slon = ', cam%slon
!pause
!print *, 'lev = ', cam%lev
!pause
!print *, 'ilev = ', cam%ilev


!-------------------------------------------------------
! define obs type and obs error variance 
!-------------------------------------------------------

   obs%type_size = 4
   allocate(obs%obstype_metadata(obs%type_size))
   obs%obstype_metadata(1)  = 'RADIOSONDE_PRESSURE'
   obs%obstype_metadata(2)  = 'RADIOSONDE_TEMPERATURE'
   obs%obstype_metadata(3)  = 'RADIOSONDE_U_WIND_COMPONENT'
   obs%obstype_metadata(4)  = 'RADIOSONDE_V_WIND_COMPONENT'

   allocate(obs%obsvar(obs%type_size))
   obs%obsvar(1) = R_P
   obs%obsvar(2) = R_T
   obs%obsvar(3) = R_U
   obs%obsvar(4) = R_U


LoopObsTypes: do idobs = indobs_s, indobs_e
   print *, 'Now processing obs type ', trim(obs%obstype_metadata(idobs))

! Initialize locind which contains the radiance for each subset
   print *, 'empirical localization kind: horizontal'
   allocate(localization%locind(locnum))
   do ind = 1, locnum
      localization%locind(ind) = locrad * real(ind)
   enddo

! allocate and initialize localization variables
   ndims = 1
   allocate(localization%num(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%sumx(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%sumy(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%sumx2(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%sumy2(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%numerator(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%denominator(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%alpha1(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%alpha2(locnum,ndims,cam%number_of_state_variables))
   allocate(localization%correlation(locnum,ndims,cam%number_of_state_variables))
   localization%num(1:locnum,1:ndims,1:cam%number_of_state_variables)          = 0
   localization%sumx(1:locnum,1:ndims,1:cam%number_of_state_variables)         = 0.0
   localization%sumy(1:locnum,1:ndims,1:cam%number_of_state_variables)         = 0.0
   localization%sumx2(1:locnum,1:ndims,1:cam%number_of_state_variables)        = 0.0
   localization%sumy2(1:locnum,1:ndims,1:cam%number_of_state_variables)        = 0.0
   localization%numerator(1:locnum,1:ndims,1:cam%number_of_state_variables)    = 0.0
   localization%denominator(1:locnum,1:ndims,1:cam%number_of_state_variables)  = 0.0
   localization%alpha1(1:locnum,1:ndims,1:cam%number_of_state_variables)       = 0.0
   localization%alpha2(1:locnum,1:ndims,1:cam%number_of_state_variables)       = 0.0
   localization%correlation(1:locnum,1:ndims,1:cam%number_of_state_variables)  = 0.0


!$OMP PARALLEL DO PRIVATE(ind,i,j,k,ii,jj,kk,inddim,ndims,indk1,indk2,iunit,writelenloc,var_name_tmp,obslon,obslat,varlon,varlat,obsi,obsj,obsk,obspres,obshgt,is_lev0,yobs,ytrue,yfmean,Robs,yvar,yprior_var,ypost_var,corr,reg_coef,ccc,ddd,dist,yinno_norm,xtrue,xfmean,xprior_var,xt_xfmean,covxy,deltaxmean,ilocind,sublocx,sublocy,sublocx2,sublocy2,sublocnum,sublocnumer,sublocalpha,sublocdenom,cam_var_1d,xprior,yprior,xtrueprf,xfmeanprf,xfspdprf,xpriorprf,xtrue3d,xfmean3d,xfspd3d,xprior3d)

!   LoopVarTypes: do inddim = 1, cam%number_of_state_variables
   LoopVarTypes: do inddim = indvar_s, indvar_e
      print *, 'processing variable ', trim(cam%description(inddim))

      ! loop obs grid points
      do k = lev_s, lev_e
         do j = lat_s, lat_e 
            do i = 1, cam%nlon

               ! get obs information ready
               ytrue      = cam%variables(i,j,k,1,idobs)
               yfmean     = cam%variables(i,j,k,2,idobs)
               yprior_var = cam%variables(i,j,k,3,idobs)**2
               Robs       = obs%obsvar(idobs)
               ypost_var    = 1.0_r8/(1.0_r8/yprior_var + 1.0_r8/Robs)
               allocate(yprior(ens_size))
               yprior(1:ens_size) = cam%variables(i,j,k,3+1:3+ens_size,idobs)

!print *, 'y info = ', ytrue, yfmean, yprior_var, Robs, ypost_var
!print *, 'yprior = ', yprior(1:ens_size:10)

               ! get [lon, lat] and convert to radiance
               obslon = cam%var_lon(i,idobs) * DEG2RAD
               obslat = cam%var_lat(j,idobs) * DEG2RAD

!print *, 'y loc = ', cam%var_lon(i,idobs), cam%var_lat(j,idobs), obslon, obslat
!pause

               ! loop target variable grid points
               do kk = lev_s, lev_e
                  do jj = lat_s, lat_e 
                     do ii = 1, cam%nlon

                     ! get target variable information ready
                     xtrue      = cam%variables(ii,jj,kk,1,inddim)
                     xfmean     = cam%variables(ii,jj,kk,2,inddim)
                     xprior_var = cam%variables(ii,jj,kk,3,inddim)**2
                     allocate(xprior(ens_size))
                     xprior(1:ens_size) = cam%variables(ii,jj,kk,3+1:3+ens_size,inddim)

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
                     xt_xfmean = xtrue - xfmean
                     covxy = comp_cov(ens_size,xfmean,xprior,yfmean,yprior)
                     corr = covxy / (sqrt(yprior_var) * sqrt(xprior_var))

                     reg_coef = covxy / yprior_var
                     ccc = reg_coef * (ypost_var / yprior_var - 1.0_r8) * yfmean
                     ddd = reg_coef * ypost_var / Robs

!print *, 'stat = ', covxy, corr, reg_coef, ccc, ddd

                     sublocx  = xt_xfmean
                     sublocy  = ccc + ddd * ytrue
                     sublocx2 = xt_xfmean**2
                     sublocy2 = (ccc + ddd*ytrue)**2

                     sublocnumer = xt_xfmean * (ccc + ddd * ytrue)
!                     sublocdenom = ccc**2 + 2.0_r8 * ccc * ddd * ytrue + ddd**2 * (Robs + ytrue**2)
                     sublocdenom = (ccc + ddd*ytrue)**2 + ddd**2 * Robs
                     
                     if ( sublocdenom >= 0.0000000001_r8 ) then
                     sublocalpha = sublocnumer/sublocdenom
!print *, 'tutu = ', xt_xfmean*covxy*(ytrue-yfmean)/(yprior_var+Robs), (covxy/(yprior_var+Robs))**2 * ((ytrue-yfmean)**2 + Robs)
!print *, 'lili diag', ccc + ddd*ytrue, covxy/(yprior_var+Robs)*(ytrue-yfmean)
!print *, 'lili diag', ddd, covxy/(yprior_var+Robs)
!print *, 'lili diag', ddd**2 * Robs, (covxy/(yprior_var+Robs))**2 * Robs
!print *, 'lili diag', ccc**2 + 2.0_r8*ccc*ddd*ytrue + (ddd**2) * (ytrue**2), (ccc+ddd*ytrue)**2, (covxy/(yprior_var+Robs))**2 * ((ytrue-yfmean)**2)
!print *, 'denom diag', ccc**2, 2*ccc*ddd*ytrue, ddd**2 * ytrue**2

!print *, 'loc = ', corr, sublocnumer, sublocdenom, sublocalpha
!pause

                     ! record in distance
                     localization%num(ilocind,1,inddim) = localization%num(ilocind,1,inddim) + 1

                     localization%sumx(ilocind,1,inddim) = localization%sumx(ilocind,1,inddim) + sublocx
                     localization%sumy(ilocind,1,inddim) = localization%sumy(ilocind,1,inddim) + sublocy

                     localization%sumx2(ilocind,1,inddim) = localization%sumx2(ilocind,1,inddim) + sublocx2
                     localization%sumy2(ilocind,1,inddim) = localization%sumy2(ilocind,1,inddim) + sublocy2

                     localization%numerator(ilocind,1,inddim) = localization%numerator(ilocind,1,inddim)     &
                                                                + sublocnumer
                     localization%denominator(ilocind,1,inddim) = localization%denominator(ilocind,1,inddim) &
                                                                  + sublocdenom
                     localization%alpha1(ilocind,1,inddim) = localization%alpha1(ilocind,1,inddim) + sublocalpha
                     localization%alpha2(ilocind,1,inddim) = localization%alpha2(ilocind,1,inddim) + sublocalpha**2
                     localization%correlation(ilocind,1,inddim) = localization%correlation(ilocind,1,inddim) + abs(corr)
                     endif

                     deallocate(xprior)

                     enddo  ! ii
                  enddo     ! jj
               enddo        ! kk

               deallocate(yprior)

            enddo  ! i
         enddo     ! j
      enddo        ! k


   enddo LoopVarTypes
!$OMP END PARALLEL DO
  
! NOTE: the writting in OMP does not work...
!       thus put the writting out of OMP loop
!   do inddim = 1, cam%number_of_state_variables
   do inddim = indvar_s, indvar_e
      writelenloc = locnum
      var_name_tmp = trim(obs%obstype_metadata(idobs))//'_'//         &
                     trim(cam%description(inddim))

      iunit = 10 + inddim
      open(iunit, FILE=var_name_tmp, FORM='FORMATTED')

      do ind = 1, ndims
         allocate(numtmp(writelenloc))
         numtmp(1:writelenloc) = localization%num(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) numtmp
         write(iunit, FMT=locreal_format) localization%sumx(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%sumy(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%sumx2(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%sumy2(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%numerator(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%denominator(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%alpha1(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%alpha2(1:writelenloc,ind,inddim)
         write(iunit, FMT=locreal_format) localization%correlation(1:writelenloc,ind,inddim)
         deallocate(numtmp)
      enddo

      close(iunit)
   enddo

   deallocate(localization%num, localization%numerator, localization%denominator)
   deallocate(localization%sumx, localization%sumy, localization%sumx2, localization%sumy2)
   deallocate(localization%alpha1, localization%alpha2, localization%correlation)

 
enddo LoopObstypes





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



