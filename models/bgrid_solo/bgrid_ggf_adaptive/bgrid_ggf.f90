
PROGRAM emp_localization

use               types_mod, only : r8, missing_r8, DEG2RAD, PI, gravity, ps0,              &
                                    ps0, gas_constant, gas_constant_v, t_kelvin
use           utilities_mod, only : find_namelist_in_file, check_namelist_read, nc_check
use       mpi_utilities_mod, only : sleep_seconds
use              meteor_mod, only : temp_and_dewpoint_to_rh, sat_vapor_pressure,            &
                                    specific_humidity
use             obs_err_mod, only : land_pres_error, land_temp_error, land_wind_error
use    dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, rh_error_from_dewpt_and_temp

use           netcdf

implicit none



integer           :: emp_loc_kind     = 1
integer           :: locnum           = 64
integer           :: ifELF            = 1
integer           :: ifINFCOR         = 1
! Bgrid model, lat (30) * lon (60), minimal distance between two grid points
! is the distance on the highest latitude within nearby two grid points, that is 0.0055;
! the 2nd minimal distance is the distance on the 2nd high latitude within nearby two grid 
! points, that is 0.0164;
! the maximal distance between two grid points is the distance on the closest to equator
! withn nearby two grid points, that is 0.1046;
! the distance betweeen two nearby grid poins with the same longitude, that is 0.1047;
real(r8)          :: locrad           = 0.05 
real(r8)          :: drho             = 0.01
integer           :: num_state_vars   = 4
integer           :: group_size       = 2
integer           :: ens_size         = 32
integer           :: ens_copy         = 36
integer           :: dimtime_start    = 348   ! start from 1 instead of 0 as ncview
integer           :: dimtime_length   = 20
!integer           :: lat_s            = 96+33
!integer           :: lat_e            = 96+75
! integer           :: lev_s            = 7
! integer           :: lev_e            = 7
!integer           :: obslev           = 5      ! change obs_ind below
!! integer           :: obslev_s         = 1
!! integer           :: obslev_e         = 1
!integer           :: varlev           = 5      ! change var_ind below
!! integer           :: varlev_s         = 7
!! integer           :: varlev_e         = 7
integer           :: indobs_s         = 1
integer           :: indobs_e         = 1
integer           :: indvar_s         = 1
integer           :: indvar_e         = 4

real              :: R_P              = 20000.0_r8 ! obs error variance of Psfc
real              :: R_T              = 2.0_r8   ! obs error variance of Temperature
real              :: R_U              = 8.0_r8   ! obs error variance of Wind
real              :: R_TPW            = 0.015    ! 0.0146 computed from RMSE fit

real              :: pressure_top     = 20000.0  ! top pressure for 'TPW'

character(len=80)  :: truth_name_input  = 'truth.nc',        &
                      prior_name_input  = 'prior',        &
                      obs_name_input    = 'obs.nc'
character(len=80), allocatable   :: prior_name_tmp(:)


type bgrid_data 
   integer  :: ntime, nlat, nlon, nslat, nslon, nlev         ! dimensions
   integer  :: filter_copy, ncopy                            ! copies of truth and prior; nmhgt
   real(r8), dimension(:),  pointer        :: time, lat, lon, slat, slon, lev
   ! lat lon for ps and t, slat and slon for u and v.
   real(r8), dimension(:,:),  pointer      :: var_lat, var_lon 
   integer                                 :: number_of_state_variables
   integer, dimension(:,:), pointer        :: var_size
   character(len=129),dimension(:),pointer :: description
   real(r8), dimension(:,:,:,:,:,:), pointer :: variables 
   real(r8), dimension(:,:,:,:,:,:), pointer :: var_inf_cor
end type bgrid_data

type(bgrid_data)    :: bgrid 


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
    real(r8), dimension(:), pointer     :: locind
    integer,  dimension(:,:), pointer   :: num
    real(r8), dimension(:,:), pointer   :: beta1     ! (sum(beta))**2
    real(r8), dimension(:,:), pointer   :: beta2     ! sum(beta**2)
    real(r8), dimension(:,:), pointer   :: meangf

    real(r8), dimension(:,:), pointer   :: incbeta1     ! (sum(beta))**2
    real(r8), dimension(:,:), pointer   :: incbeta2     ! sum(beta**2)
    real(r8), dimension(:,:), pointer   :: incmeangf

    real(r8), dimension(:,:), pointer   :: elfbeta1
    real(r8), dimension(:,:), pointer   :: elfbeta2
    real(r8), dimension(:,:), pointer   :: elfmeangf

    integer,  dimension(:,:,:), pointer   :: num2        ! num for each group
    real(r8), dimension(:,:,:), pointer   :: sumx        ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:), pointer   :: sumy        ! used to get mean(y)
    real(r8), dimension(:,:,:), pointer   :: sumx2       ! used to get mean(x^2)
    real(r8), dimension(:,:,:), pointer   :: sumy2       ! used to get mean(y^2)
    real(r8), dimension(:,:,:), pointer   :: numerator   ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:), pointer   :: denominator ! i.e., used to get mean(y^2) (note: it contains Robs)
    real(r8), dimension(:,:,:), pointer   :: alpha1
    real(r8), dimension(:,:,:), pointer   :: alpha2
    real(r8), dimension(:,:,:), pointer   :: correlation
    real(r8), dimension(:,:,:), pointer   :: sumcorr
    real(r8), dimension(:,:,:), pointer   :: sumcorrsq

    integer,  dimension(:,:,:,:), pointer :: num_rf
    real(r8), dimension(:,:,:,:), pointer :: gbetasq
    real(r8), dimension(:,:,:,:), pointer :: tgbeta
    real(r8), dimension(:,:,:,:), pointer :: gbetasq_infcor
    real(r8), dimension(:,:,:,:), pointer :: tgbeta_infcor

    integer,  dimension(:,:), pointer   :: num_infcor
    real(r8), dimension(:,:), pointer   :: beta1_infcor     ! (sum(beta))**2
    real(r8), dimension(:,:), pointer   :: beta2_infcor     ! sum(beta**2)
    real(r8), dimension(:,:), pointer   :: meangf_infcor

    real(r8), dimension(:,:), pointer   :: incbeta1_infcor
    real(r8), dimension(:,:), pointer   :: incbeta2_infcor
    real(r8), dimension(:,:), pointer   :: incmeangf_infcor

    real(r8), dimension(:,:), pointer   :: elfbeta1_infcor
    real(r8), dimension(:,:), pointer   :: elfbeta2_infcor
    real(r8), dimension(:,:), pointer   :: elfmeangf_infcor

    integer,  dimension(:,:,:), pointer   :: num2_infcor        ! num for each group
    real(r8), dimension(:,:,:), pointer   :: sumx_infcor        ! used to get mean(x) in least squares fit
    real(r8), dimension(:,:,:), pointer   :: sumy_infcor        ! used to get mean(y)
    real(r8), dimension(:,:,:), pointer   :: sumx2_infcor       ! used to get mean(x^2)
    real(r8), dimension(:,:,:), pointer   :: sumy2_infcor       ! used to get mean(y^2)
    real(r8), dimension(:,:,:), pointer   :: numerator_infcor   ! i.e., used to get mean(xy)
    real(r8), dimension(:,:,:), pointer   :: denominator_infcor ! i.e., used to get mean(y^2) (note: it contains Robs)
    real(r8), dimension(:,:,:), pointer   :: alpha1_infcor
    real(r8), dimension(:,:,:), pointer   :: alpha2_infcor
    real(r8), dimension(:,:,:), pointer   :: correlation_infcor
    real(r8), dimension(:,:,:), pointer   :: sumcorr_infcor
    real(r8), dimension(:,:,:), pointer   :: sumcorrsq_infcor

end type localization_data

type(localization_data) localization

type localization_omp_data
  type(localization_data), pointer  :: obs(:)
  integer                           :: nobs
end type localization_omp_data

type(localization_omp_data) loc_for_omp


integer          :: i, j, k, id, idobs, ind, inddim, indvars, indens,      &
                    ii, jj, kk, itime, indk1, indk2, var_id, ndims,        &
                    obslevk, varlevk, numobslev, numvarlev
integer          :: nrho, indcorr
integer          :: ig, igs, ige, igg
integer          :: obs_lat_s, obs_lat_e, obs_lon_s, obs_lon_e,            &
                    var_lat_s, var_lat_e, var_lon_s, var_lon_e
integer          :: ncid, ncidtruth, ncidprior, ncidobs, iunit
integer          :: proj_code
integer          :: ilocind, sublocnum, writelenloc
integer          :: obs_ind(4)    ! for mixlev us ( lev1, 7, 12, 25)
integer          :: var_ind(4,4)  ! for mixlev use (lev1, 7, 12, 25)

real(r8)         :: stdlon,truelat1,truelat2,latinc,loninc
real(r8)         :: obslat, obslon, obspres, obsi, obsj, obsk, obspalt, obshgt
real(r8)         :: varlat, varlon
real(r8)         :: yobs, ytrue, yfmean, oerr, Robs, yprior_var, ypost_var, yvar, yinno_norm,              &
                    xfmean, xtrue, xprior_var,xt_xfmean, covxy, deltaxmean, ccc1, ddd1, dist
real(r8)         :: yprior_var_infcor, ypost_var_infcor, xprior_var_infcor, covxy_infcor
real(r8)         :: xvar_ratio
real(r8)         :: obsdiag1, obsdiag2, obsdiag3, obsdiag4
real(r8)         :: sumbeta1, sumbeta2, elfsumbeta1, elfsumbeta2, incsumbeta1, incsumbeta2
real(r8)         :: ref_P, surf_P

logical          :: is_lev0
character(len=2) :: idchar1, idchar2
character(len=1) :: igchar
character(len=100) :: obsnum_format, locint_format, locreal_format, mhgt_format

integer          :: TimeDimID, CopyDimID
integer          :: dimids(NF90_MAX_DIMS)
character(len=NF90_MAX_NAME) :: var_name_tmp, dimname
!integer, dimension(NF90_MAX_DIMS) :: ncstart, nccount


integer,  allocatable :: ncidpriors(:)
integer,  allocatable :: ncdims(:),ncstart(:,:), nccount(:,:), ncdims_true(:), ncstart_true(:,:), nccount_true(:,:)
real(r8), allocatable :: bgrid_var_1d(:), bgrid_var_2d(:,:), bgrid_var_3d(:,:,:), bgrid_var_4d(:,:,:,:),    &
                                                             bgrid_var_3d_2(:,:,:), bgrid_var_5d(:,:,:,:,:)
real(r8), allocatable :: bgrid_hyam(:), bgrid_hybm(:), qv_1d(:)
real(r8), allocatable :: obs_var_1d(:), obs_var_2d(:,:), obs_var_3d(:,:,:), obs_var_4d(:,:,:,:)
real(r8), allocatable :: wrf_mu_d01(:,:), wrf_psfc_d01(:,:), wrf_ph_d01(:,:,:),       &
                         wrf_t_d01(:,:,:), wrf_qvapor_d01(:,:,:)
real(r8), allocatable :: xprior(:), yprior(:)
real(r8), allocatable :: xtrueprf(:), xfmeanprf(:), xfspdprf(:), xpriorprf(:,:)
real(r8), allocatable :: xtrue3d(:,:,:), xfmean3d(:,:,:), xfspd3d(:,:,:), xprior3d(:,:,:,:) 
real(r8), allocatable :: numtmp(:)
real(r8), allocatable :: ccc(:), ddd(:), ccc_infcor(:), ddd_infcor(:)
real(r8), allocatable :: reg_coef(:), corr(:), sublocx(:), sublocy(:), sublocx2(:), sublocy2(:),            &
                         sublocnumer(:), sublocdenom(:), sublocalpha(:),                                    &
                         incnumer(:), incdenom(:), incalpha(:)
real(r8), allocatable :: reg_coef_infcor(:),     &
                         sublocnumer_infcor(:), sublocdenom_infcor(:), sublocalpha_infcor(:),               &
                         incnumer_infcor(:), incdenom_infcor(:), incalpha_infcor(:)
real(r8), allocatable :: subset_regcoef(:)


! define format to write "localization"
locint_format = ' ( 2000i10 ) '
locreal_format = ' ( 2000e15.7 ) '
! define format to write obs number
!obsnum_format = ' ( i10,4e15.7 )'
obsnum_format = ' ( i10 ) '
mhgt_format = ' ( 49e15.7 ) '


!-------------------------------------------------------
! open truth file (True_State.nc) 
!-------------------------------------------------------
if ( ifELF == 1 ) then
call nc_check( nf90_open(truth_name_input, NF90_NOWRITE, ncidtruth),   &
               'bgrid_emp_localization','open truth file')
endif

!-------------------------------------------------------
! open prior file (Prior_Diag.nc) 
!-------------------------------------------------------
allocate(prior_name_tmp(group_size))
allocate(ncidpriors(group_size))
do i = 1, group_size
   write(igchar,'(i1.1)') i
   prior_name_tmp(i) = trim(prior_name_input)//igchar//'.nc'
   call nc_check( nf90_open(trim(prior_name_tmp(i)), NF90_NOWRITE, ncidpriors(i)),   &
                  'bgrid_emp_localization','open prior file')
enddo

!-------------------------------------------------------
! fill in state variable information 
!-------------------------------------------------------
   bgrid%number_of_state_variables = num_state_vars
   allocate(bgrid%description(bgrid%number_of_state_variables))
   bgrid%description(1)  = 'ps'
   bgrid%description(2)  = 't'
   bgrid%description(3)  = 'u'
   bgrid%description(4)  = 'v'


!-------------------------------------------------------
! read CAM dimensions
!-------------------------------------------------------
   ! dimension of copies (1(truth)+2(ens mean, ens spread)+ens_size)
   ! NOTE: (inflation is turned off in the OSSE)
   bgrid%ncopy = ens_copy * group_size + 1

   var_name_tmp = 'time'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), TimeDimID),       &
                  'read_bgrid_dimensions','inq_dimid time')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), TimeDimID, name=dimname, len=bgrid%ntime),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'copy'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid TmpJ')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%filter_copy),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'TmpJ'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid TmpJ')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%nlat),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'VelJ'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid VelJ')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%nslat),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'TmpI'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid TmpI')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%nlon),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'VelI'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid VelI')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%nslon),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

   var_name_tmp = 'lev'
   call nc_check( nf90_inq_dimid(ncidpriors(1), trim(var_name_tmp), var_id),       &
                  'read_bgrid_dimensions','inq_dimid nlev')
   call nc_check( nf90_inquire_dimension(ncidpriors(1), var_id, name=dimname, len=bgrid%nlev),  &
                  'read_bgrid_dimensions','inquire_dimension'//trim(dimname))

print *, 'file dimenisons of time = ', bgrid%ntime
print *, 'file dimensions of lat and lon = ', bgrid%nlat, bgrid%nlon
print *, 'file dimensions of staggered lat and lon = ', bgrid%nslat, bgrid%nslon
print *, 'file dimensions of lev = ', bgrid%nlev
print *, 'timedimid = ', TimeDimID
!pause

!-------------------------------------------------------
! read CAM dimension data 
!-------------------------------------------------------
! 1. 1D array

   allocate(bgrid%lat(bgrid%nlat))
   var_name_tmp = 'TmpJ'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_bgrid_dimension_data','inq_varid TmpJ')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, bgrid%lat),  &
                     'read_bgrid_dimension_data','get_var TmpJ')

   allocate(bgrid%slat(bgrid%nslat))
   var_name_tmp = 'VelJ'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_bgrid_dimension_data','inq_varid VelJ')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, bgrid%slat),  &
                     'read_bgrid_dimension_data','get_var VelJ')

   allocate(bgrid%lon(bgrid%nlon))
   var_name_tmp = 'TmpI'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_bgrid_dimension_data','inq_varid TmpI')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, bgrid%lon),  &
                     'read_bgrid_dimension_data','get_var TmpI')

   allocate(bgrid%slon(bgrid%nslon))
   var_name_tmp = 'VelI'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_bgrid_dimension_data','inq_varid VelI')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, bgrid%slon),  &
                     'read_bgrid_dimension_data','get_var VelI')

   allocate(bgrid%lev(bgrid%nlev))
   var_name_tmp = 'lev'
   call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'read_bgrid_dimension_data','inq_varid lev')
   call nc_check( nf90_get_var(ncidpriors(1), var_id, bgrid%lev),  &
                     'read_bgrid_dimension_data','get_var lev')


!-------------------------------------------------------
! read state variables of True_State.nc (copy 1) 
!-------------------------------------------------------
   ! var_size(1,:) - dimension of lon, (i) 
   ! var_size(2,:) - dimension of lat, (j)
   ! var_size(3,:) - dimension of lev  (k) ( = 1, means 2D variable, others 3D variable)

   ! variables dimensions: 1 - lat, 2 - lon, 3 - lev, 
   !                       4 - copies, 5 - state variables

   allocate(ncdims(bgrid%number_of_state_variables))
   allocate(ncstart(NF90_MAX_DIMS,bgrid%number_of_state_variables))
   allocate(nccount(NF90_MAX_DIMS,bgrid%number_of_state_variables))
   allocate(bgrid%var_size(NF90_MAX_DIMS,bgrid%number_of_state_variables))
   allocate(bgrid%variables(bgrid%nlon,bgrid%nlat,bgrid%nlev,bgrid%ncopy,    &
                            dimtime_length,bgrid%number_of_state_variables))
   allocate(bgrid%var_lat(bgrid%nlat,bgrid%number_of_state_variables))
   allocate(bgrid%var_lon(bgrid%nlon,bgrid%number_of_state_variables))
   
   do ind = 1,bgrid%number_of_state_variables
!   do ind = 1, 1

      ! Get the dimension size 
      ! Once this is done for variable of True_State.nc, there is no need to do this for
      ! variables in Prior_Diag.nc, since they should be consistent
      var_name_tmp = trim(bgrid%description(ind))
      call nc_check( nf90_inq_varid(ncidpriors(1), var_name_tmp, var_id),  &
                     'get_variable_size_from_file',                    &
                     'inq_varid '//var_name_tmp)
      call nc_check( nf90_inquire_variable(ncidpriors(1), var_id,          &
                     ndims=ndims, dimids=dimids),                      &
                     'get_variable_size_from_file',                    &
                     'inquire_variable '//var_name_tmp)

      ncdims(ind) = ndims
      ncstart(:,ind) = 1
      do inddim = 1,ndims
         call nc_check( nf90_inquire_dimension(ncidpriors(1), dimids(inddim),  &
                        len=bgrid%var_size(inddim,ind)),                     &
                        'get_variable_size_from_file',                     &
                        'inquire_dimension '//var_name_tmp)
         nccount(inddim,ind) = bgrid%var_size(inddim,ind)
      enddo

      where(dimids == TimeDimID) ncstart(:,ind) = dimtime_start
      where(dimids == TimeDimID) nccount(:,ind) = dimtime_length

      ! Get the data
      if ( ndims < 5 ) then

         ! 2D variable, its vertical dimension is 1
         bgrid%var_size(ndims+1,ind) = bgrid%var_size(ndims,ind)
         bgrid%var_size(ndims,ind)   = bgrid%var_size(ndims-1,ind)
         bgrid%var_size(ndims-1,ind) = 1

      elseif ( ndims > 5 ) then

         print *, 'ERROR: variable ', trim(var_name_tmp),' has ',ndims,' dimensions.'
         print *, 'DO NOT know what to do... Stopping.'
         stop

      endif

      ! get the lat, lon ready for each variable (US, VS are staggered)
      if ( var_name_tmp == 'u' .or. var_name_tmp == 'v' ) then
         bgrid%var_lat(1:bgrid%nslat,ind) = bgrid%slat(1:bgrid%nslat)
         bgrid%var_lon(1:bgrid%nslon,ind) = bgrid%slon(1:bgrid%nslon)
      else
         bgrid%var_lat(1:bgrid%nlat,ind)  = bgrid%lat(1:bgrid%nlat)
         bgrid%var_lon(1:bgrid%nlon,ind)  = bgrid%lon(1:bgrid%nlon)
      endif

print *, 'dimension of variable ', trim(var_name_tmp)
print *, 'ncdims =  ', ncdims(ind)
print *, 'ncstart = ', ncstart(1:5,ind)
print *, 'nccount = ', nccount(1:5,ind)
print *, 'var size = ', bgrid%var_size(1:5,ind)

   enddo

print *, 'end of reading dimension'



!-------------------------------------------------------
! read state variables of Prio_Diag.nc (copy ens_mean, ens_spd, ens_member) 
!-------------------------------------------------------
   do ind = 1,bgrid%number_of_state_variables

      do i = 1, group_size

         igs = (i-1)*ens_copy + 1
         ige = i*ens_copy
   
         var_name_tmp = trim(bgrid%description(ind))

         call nc_check( nf90_inq_varid(ncidpriors(i), var_name_tmp, var_id),  &
                        'get_variable_name_from_file',                        &
                        'inq_varid '//var_name_tmp)

         if ( bgrid%var_size(3,ind) == 1 ) then

            ! Get 2D variable 
            ! NOTE: nf90_get_var automatically ignor the dimension with length 1
            allocate(bgrid_var_4d(bgrid%var_size(1,ind),                    &
                                  bgrid%var_size(2,ind),                    &
                                  ens_copy,                                 &
                                  dimtime_length))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, bgrid_var_4d, &
                        start=ncstart(1:ncdims(ind),ind),                    &
                        count=nccount(1:ncdims(ind),ind)),   &
                        'read_in_bgrid_var','get_var '//var_name_tmp )
            bgrid%variables(1:bgrid%var_size(1,ind),                        &
                            1:bgrid%var_size(2,ind),                        &
                            1,                                              &
                            igs:ige,                                        &
                            1:dimtime_length,                               &
                            ind) =                                          &
                       bgrid_var_4d(1:bgrid%var_size(1,ind),                &
                                    1:bgrid%var_size(2,ind),                &
                                    1:ens_copy,                             &
                                    1:dimtime_length)
            deallocate(bgrid_var_4d)

         else

            ! Get 3D variable
            ! NOTE: var_size was reshaped when read in True_State.nc,
            !       but nc_varget starts from lev. this is how to set bgrid_var_4d
            allocate(bgrid_var_5d(bgrid%var_size(1,ind),                    &
                                  bgrid%var_size(2,ind),                    &
                                  bgrid%var_size(3,ind),                    &
                                  ens_copy,                                 &
                                  dimtime_length))
            call nc_check( nf90_get_var(ncidpriors(i), var_id, bgrid_var_5d, &
                        start=ncstart(1:ncdims(ind),ind),                    &
                        count=nccount(1:ncdims(ind),ind)),                   &
                           'read_in_bgrid_var','get_var '//var_name_tmp )
            bgrid%variables(1:bgrid%var_size(1,ind),                       &
                            1:bgrid%var_size(2,ind),                       &
                            1:bgrid%var_size(3,ind),                       &
                            igs:ige,                                       &
                            1:dimtime_length,                              &
                            ind) =                                         &
                       bgrid_var_5d(1:bgrid%var_size(1,ind),               &
                                    1:bgrid%var_size(2,ind),               &
                                    1:bgrid%var_size(3,ind),               &
                                    1:ens_copy,                            &
                                    1:dimtime_length)
            deallocate(bgrid_var_5d)

         endif

      enddo  ! i fro group_size

   enddo   ! ind for number_of_state_variables

print *, 'end of reading prior.nc'

!-------------------------------------------------------
! close netcdf files
!-------------------------------------------------------
  do i = 1, group_size
     call nc_check(nf90_close(ncidpriors(i)),'bgrid_emp_localization','close prior file')
  enddo
  deallocate(prior_name_tmp,ncidpriors)


!-------------------------------------------------------
! ifELF = 1, read in True_State.nc 
!-------------------------------------------------------
if ( ifELF == 1 ) then

   allocate(ncdims_true(bgrid%number_of_state_variables))
   allocate(ncstart_true(NF90_MAX_DIMS,bgrid%number_of_state_variables))
   allocate(nccount_true(NF90_MAX_DIMS,bgrid%number_of_state_variables))

   ncdims_true = ncdims
   ncstart_true = ncstart
   nccount_true = nccount
   do i = 1, bgrid%number_of_state_variables
      if ( bgrid%var_size(3,i) == 1 ) then
         nccount_true(3,i) = 1
      else
         nccount_true(4,i) = 1
      endif
   enddo


print *, 'ncdims_true =  ', ncdims_true
print *, 'ncstart_true = ', ncstart_true(1:5,:)
print *, 'nccount_true = ', nccount_true(1:5,:)
!pause

   do ind = 1,bgrid%number_of_state_variables

         var_name_tmp = trim(bgrid%description(ind))

         call nc_check( nf90_inq_varid(ncidtruth, var_name_tmp, var_id),  &
                        'get_variable_name_from_file',                    &
                        'inq_varid '//var_name_tmp)

         if ( bgrid%var_size(3,ind) == 1 ) then

            ! Get 2D variable 
            ! NOTE: nf90_get_var automatically ignor the dimension with length 1
            allocate(bgrid_var_3d(bgrid%var_size(1,ind),                    &
                                  bgrid%var_size(2,ind),                    &
                                  dimtime_length))
            call nc_check( nf90_get_var(ncidtruth, var_id, bgrid_var_3d,    &
                        start=ncstart_true(1:ncdims_true(ind),ind),         &
                        count=nccount_true(1:ncdims_true(ind),ind)),        &
                        'read_in_bgrid_var','get_var '//var_name_tmp )
            bgrid%variables(1:bgrid%var_size(1,ind),                        &
                            1:bgrid%var_size(2,ind),                        &
                            1,                                              &
                            bgrid%ncopy,                                    &
                            1:dimtime_length,                               &
                            ind) =                                          &
                       bgrid_var_3d(1:bgrid%var_size(1,ind),                &
                                    1:bgrid%var_size(2,ind),                &
                                    1:dimtime_length)
            deallocate(bgrid_var_3d)

         else

            ! Get 3D variable
            ! NOTE: var_size was reshaped when read in True_State.nc,
            !       but nc_varget starts from lev. this is how to set bgrid_var_4d
            allocate(bgrid_var_4d(bgrid%var_size(1,ind),                    &
                                  bgrid%var_size(2,ind),                    &
                                  bgrid%var_size(3,ind),                    &
                                  dimtime_length))
            call nc_check( nf90_get_var(ncidtruth, var_id, bgrid_var_4d,    &
                        start=ncstart_true(1:ncdims_true(ind),ind),         &
                        count=nccount_true(1:ncdims_true(ind),ind)),        &
                           'read_in_bgrid_var','get_var '//var_name_tmp )
            bgrid%variables(1:bgrid%var_size(1,ind),                       &
                            1:bgrid%var_size(2,ind),                       &
                            1:bgrid%var_size(3,ind),                       &
                            bgrid%ncopy,                                   &
                            1:dimtime_length,                              &
                            ind) =                                         &
                       bgrid_var_4d(1:bgrid%var_size(1,ind),               &
                                    1:bgrid%var_size(2,ind),               &
                                    1:bgrid%var_size(3,ind),               &
                                    1:dimtime_length)
            deallocate(bgrid_var_4d)

       endif

   enddo

   print *, 'end of reading truth.nc'

endif


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


!-------------------------------------------------------
! compute the inflation correction (correct spread)
!-------------------------------------------------------
if ( ifELF == 1 .and. ifINFCOR == 1 ) then
   allocate(bgrid%var_inf_cor(bgrid%nlon,bgrid%nlat,bgrid%nlev,group_size,        &
                              dimtime_length,bgrid%number_of_state_variables))
   allocate(sublocnumer(group_size),sublocdenom(group_size),sublocalpha(group_size))
   allocate(reg_coef(group_size))

   do inddim = 1, bgrid%number_of_state_variables

      if ( inddim == 1 .or. inddim == 2 ) then
           var_lat_s = 1
           var_lat_e = bgrid%nlat
           var_lon_s = 1
           var_lon_e = bgrid%nlon
      else
           var_lat_s = 1
           var_lat_e = bgrid%nslat
           var_lon_s = 1
           var_lon_e = bgrid%nslon
      endif

      if ( inddim == 1 ) then
           numvarlev = 1
      else
           numvarlev = 5
      endif

      do kk = 1, numvarlev

         sublocnumer = 0.0_r8
         sublocdenom = 0.0_r8

         do jj = var_lat_s, var_lat_e
         do ii = var_lon_s, var_lon_e
            do itime = 1, dimtime_length
            do ig = 1, group_size

               igs = (ig-1)*ens_copy + 1
               ige = ig*ens_copy

! NOTE: Robs should be changed when state variables is not consistent with observations
               ! get obs information ready
               Robs       = obs%obsvar(inddim)

               ! get target variable information ready
               xfmean     = bgrid%variables(ii,jj,kk,igs,itime,inddim)
               xprior_var = bgrid%variables(ii,jj,kk,igs+1,itime,inddim)**2
               xtrue      = bgrid%variables(ii,jj,kk,bgrid%ncopy,itime,inddim)

               xt_xfmean    = xtrue-xfmean
               reg_coef(ig) = xprior_var/(xprior_var + Robs)
               ccc1         = xt_xfmean * reg_coef(ig) * xt_xfmean
               ddd1         = ( reg_coef(ig)**2 ) * ( xt_xfmean**2 + Robs)

               sublocnumer(ig)  = sublocnumer(ig) + ccc1
               sublocdenom(ig)  = sublocdenom(ig) + ddd1

            enddo   ! ig
            enddo   ! itime
         enddo   ! ii
         enddo   ! jj

         sublocalpha = sublocnumer / sublocdenom
print *, 'sublocalpha = ', sublocalpha


         ! compute the inflation correction for state variable at each grid point at a given time
         do jj = var_lat_s, var_lat_e
         do ii = var_lon_s, var_lon_e
            do itime = 1, dimtime_length
            do ig = 1, group_size

               igs = (ig-1)*ens_copy + 1
               ige = ig*ens_copy

! NOTE: Robs should be changed when state variables is not consistent with observations
               Robs       = obs%obsvar(inddim)
               xprior_var = bgrid%variables(ii,jj,kk,igs+1,itime,inddim)**2

               xvar_ratio = Robs/xprior_var

               if ( xvar_ratio > (sublocalpha(ig) - 1.0_r8) ) then               
                  bgrid%var_inf_cor(ii,jj,kk,ig,itime,inddim) =                 &
                     sublocalpha(ig) / ((1-sublocalpha(ig))*xprior_var/Robs + 1)
               else
                  bgrid%var_inf_cor(ii,jj,kk,ig,itime,inddim) = 1.0_r8
               endif

            enddo   ! ig
            enddo   ! itime
         enddo   ! ii
         enddo   ! jj

      enddo   ! kk

   enddo   ! inddim

   deallocate(sublocnumer, sublocdenom, sublocalpha)
   deallocate(reg_coef)

endif   ! ifELF, inINFCOR



!-------------------------------------------------------
! Section for computing the localization
!-------------------------------------------------------

! Initialize locind which contains the radiance for each subset
   allocate(localization%locind(locnum))
   do ind = 1, locnum
      if ( ind == 1 ) then
         localization%locind(ind) = 0.0_r8
      else
         localization%locind(ind) = locrad * real(ind-1)
      endif
   enddo

!  the number of correlation bin values
!  NOTE: for bin j, correlation dropped in with binvalue(j) <= corr < binvalue(j+1)
!        for bin nrho, correlation dropped in with corr == bin_value(nrho)
!        corr less than binvalue(1) and corr larger than binvalue(nrho) will not be counted
   nrho = int((1.0_r8+drho/2.0_r8)*2.0_r8/drho) + 1


LoopObsTypes: do idobs = indobs_s, indobs_e
   print *, 'Now processing obs type ', trim(obs%obstype_metadata(idobs))

   if ( idobs == 1 .or. idobs == 2 ) then
        obs_lat_s = 1
        obs_lat_e = bgrid%nlat
        obs_lon_s = 1
        obs_lon_e = bgrid%nlon
   else
        obs_lat_s = 1
        obs_lat_e = bgrid%nslat
        obs_lon_s = 1
        obs_lon_e = bgrid%nslon
   endif

! allocate and initialize localization variables
!  ndims for each level give a obs and a variable, thus 16 = 1 (ps) + 3 (t,u,v) * 5
   if ( idobs == 1 ) then
        numobslev = 1 
   else
        numobslev = 5
   endif

   ! loop through observation vertical levels
   LoopObsLevels: do obslevk = 1, numobslev
      k = obslevk

      ndims = 16
      allocate(localization%num(locnum,ndims))
      allocate(localization%beta1(locnum,ndims))
      allocate(localization%beta2(locnum,ndims))
      allocate(localization%meangf(locnum,ndims))

      allocate(localization%incbeta1(locnum,ndims))
      allocate(localization%incbeta2(locnum,ndims))
      allocate(localization%incmeangf(locnum,ndims))
   
      allocate(localization%elfbeta1(locnum,ndims))
      allocate(localization%elfbeta2(locnum,ndims))
      allocate(localization%elfmeangf(locnum,ndims))
   
      allocate(localization%num2(locnum,group_size,ndims))
      allocate(localization%sumx(locnum,group_size,ndims))
      allocate(localization%sumy(locnum,group_size,ndims))
      allocate(localization%sumx2(locnum,group_size,ndims))
      allocate(localization%sumy2(locnum,group_size,ndims))
      allocate(localization%numerator(locnum,group_size,ndims))
      allocate(localization%denominator(locnum,group_size,ndims))
      allocate(localization%alpha1(locnum,group_size,ndims))
      allocate(localization%alpha2(locnum,group_size,ndims))
      allocate(localization%correlation(locnum,group_size,ndims))
      allocate(localization%sumcorr(locnum,group_size,ndims))
      allocate(localization%sumcorrsq(locnum,group_size,ndims))

      allocate(localization%num_rf(locnum,nrho,group_size,ndims))
      allocate(localization%gbetasq(locnum,nrho,group_size,ndims))
      allocate(localization%tgbeta(locnum,nrho,group_size,ndims))
      allocate(localization%gbetasq_infcor(locnum,nrho,group_size,ndims))
      allocate(localization%tgbeta_infcor(locnum,nrho,group_size,ndims))

      allocate(localization%num_infcor(locnum,ndims))
      allocate(localization%beta1_infcor(locnum,ndims))
      allocate(localization%beta2_infcor(locnum,ndims))
      allocate(localization%meangf_infcor(locnum,ndims))

      allocate(localization%incbeta1_infcor(locnum,ndims))
      allocate(localization%incbeta2_infcor(locnum,ndims))
      allocate(localization%incmeangf_infcor(locnum,ndims))

      allocate(localization%elfbeta1_infcor(locnum,ndims))
      allocate(localization%elfbeta2_infcor(locnum,ndims))
      allocate(localization%elfmeangf_infcor(locnum,ndims))

      allocate(localization%num2_infcor(locnum,group_size,ndims))
      allocate(localization%sumx_infcor(locnum,group_size,ndims))
      allocate(localization%sumy_infcor(locnum,group_size,ndims))
      allocate(localization%sumx2_infcor(locnum,group_size,ndims))
      allocate(localization%sumy2_infcor(locnum,group_size,ndims))
      allocate(localization%numerator_infcor(locnum,group_size,ndims))
      allocate(localization%denominator_infcor(locnum,group_size,ndims))
      allocate(localization%alpha1_infcor(locnum,group_size,ndims))
      allocate(localization%alpha2_infcor(locnum,group_size,ndims))
      allocate(localization%correlation_infcor(locnum,group_size,ndims))
      allocate(localization%sumcorr_infcor(locnum,group_size,ndims))
      allocate(localization%sumcorrsq_infcor(locnum,group_size,ndims))  
 
      localization%num          = 0
      localization%beta1        = 0.0_r8
      localization%beta2        = 0.0_r8
      localization%meangf       = 0.0_r8

      localization%incbeta1        = 0.0_r8
      localization%incbeta2        = 0.0_r8
      localization%incmeangf       = 0.0_r8

      localization%elfbeta1     = 0.0_r8
      localization%elfbeta2     = 0.0_r8
      localization%elfmeangf    = 0.0_r8
   
      localization%num2         = 0
      localization%sumx         = 0.0_r8
      localization%sumy         = 0.0_r8
      localization%sumx2        = 0.0_r8
      localization%sumy2        = 0.0_r8
      localization%numerator    = 0.0_r8
      localization%denominator  = 0.0_r8
      localization%alpha1       = 0.0_r8
      localization%alpha2       = 0.0_r8
      localization%correlation  = 0.0_r8
      localization%sumcorr      = 0.0_r8
      localization%sumcorrsq    = 0.0_r8

      localization%num_rf         = 0
      localization%gbetasq        = 0.0_r8
      localization%tgbeta         = 0.0_r8
      localization%gbetasq_infcor = 0.0_r8
      localization%tgbeta_infcor  = 0.0_r8

      localization%num_infcor          = 0
      localization%beta1_infcor        = 0.0_r8
      localization%beta2_infcor        = 0.0_r8
      localization%meangf_infcor       = 0.0_r8

      localization%incbeta1_infcor        = 0.0_r8
      localization%incbeta2_infcor        = 0.0_r8
      localization%incmeangf_infcor       = 0.0_r8

      localization%elfbeta1_infcor     = 0.0_r8
      localization%elfbeta2_infcor     = 0.0_r8
      localization%elfmeangf_infcor    = 0.0_r8

      localization%num2_infcor         = 0
      localization%sumx_infcor         = 0.0_r8
      localization%sumy_infcor         = 0.0_r8
      localization%sumx2_infcor        = 0.0_r8
      localization%sumy2_infcor        = 0.0_r8
      localization%numerator_infcor    = 0.0_r8
      localization%denominator_infcor  = 0.0_r8
      localization%alpha1_infcor       = 0.0_r8
      localization%alpha2_infcor       = 0.0_r8
      localization%correlation_infcor  = 0.0_r8
      localization%sumcorr_infcor      = 0.0_r8
      localization%sumcorrsq_infcor    = 0.0_r8

! deleted {k, inddim} from PRIVATE below

!!!$OMP PARALLEL DO PRIVATE(ind,i,j,ii,jj,kk,obslevk,varlevk,indvars,inddim,indk1,indk2,ig,igs,ige,igg,itime,iunit,writelenloc,ref_P,surf_P,var_name_tmp,obslon,obslat,varlon,varlat,var_lat_s,var_lat_e,var_lon_s,var_lon_e,numvarlev,obsi,obsj,obsk,obspres,obshgt,is_lev0,yobs,ytrue,yfmean,Robs,yvar,yprior_var,ypost_var,corr,reg_coef,subset_regcoef,ccc,ddd,ccc_infcor,ddd_infcor,dist,yinno_norm,xtrue,xfmean,xprior_var,xt_xfmean,covxy,deltaxmean,indcorr,ilocind,sublocx,sublocy,sublocx2,sublocy2,sublocnum,sublocnumer,sublocalpha,sublocdenom,bgrid_var_1d,qv_1d,xprior,yprior,xtrueprf,xfmeanprf,xfspdprf,xpriorprf,xtrue3d,xfmean3d,xfspd3d,xprior3d,sumbeta1,sumbeta2,elfsumbeta1,elfsumbeta2,incsumbeta1,incsumbeta2,yprior_var_infcor,ypost_var_infcor,xprior_var_infcor,covxy_infcor,reg_coef_infcor,sublocnumer_infcor,sublocdenom_infcor,sublocalpha_infcor,incnumer,incdenom,incalpha,incnumer_infcor,incdenom_infcor,incalpha_infcor)

!      LoopVarIndex: do indvars = 1, ndims
      LoopVarIndex: do indvars = 1, 1

         ! find the state variable index and its vertical index
         if ( indvars == 1 ) then
              inddim = 1
              kk = 1
         else
              if ( indvars <= 6 ) then
                   inddim = 2
              elseif ( indvars <= 11 ) then
                   inddim = 3
              else
                   inddim = 4
              endif
              indk1 = mod(indvars-1,5)
              if ( indk1 == 0 ) then
                   kk = 5
              else
                   kk = indk1
              endif
         endif

print *, 'ind test = ', indvars, inddim, kk, k, ndims
!pause

         if ( indvars <=6 ) then
              var_lat_s = 1
              var_lat_e = bgrid%nlat
              var_lon_s = 1
              var_lon_e = bgrid%nlon
         else
              var_lat_s = 1
              var_lat_e = bgrid%nslat
              var_lon_s = 1
              var_lon_e = bgrid%nslon
         endif

         allocate(yprior(ens_size))
         allocate(xprior(ens_size))

         allocate(ccc(group_size),ddd(group_size))
         allocate(corr(group_size),reg_coef(group_size))
         allocate(sublocx(group_size),sublocy(group_size))
         allocate(sublocx2(group_size),sublocy2(group_size))
         allocate(sublocnumer(group_size),sublocdenom(group_size),sublocalpha(group_size))
         allocate(incnumer(group_size),incdenom(group_size),incalpha(group_size))

         allocate(ccc_infcor(group_size),ddd_infcor(group_size))
         allocate(reg_coef_infcor(group_size))
         allocate(sublocnumer_infcor(group_size),sublocdenom_infcor(group_size),sublocalpha_infcor(group_size))
         allocate(incnumer_infcor(group_size),incdenom_infcor(group_size),incalpha_infcor(group_size))

         allocate(subset_regcoef(group_size-1))

!         ! loop obs grid points
!         do j = obs_lat_s, obs_lat_e 
!            do i = obs_lon_s, obs_lon_e
!
!            ! loop target variable grid points
!               do jj = var_lat_s, var_lat_e 
!                  do ii = var_lon_s, var_lon_e

!! obs along a latitude
!         do j = 17, 17 
!            do i = obs_lon_s, obs_lon_e
!               do jj = 17, 17 
!                  do ii = var_lon_s, var_lon_e

!! obs along a longitude
!         do j = obs_lat_s, obs_lat_e 
!            do i = 1, 1
!               do jj = var_lat_s, var_lat_e 
!                  do ii = 1, 1

!! obs along latitudes and var on the same latitude as obs
!         do j = obs_lat_s, obs_lat_e 
!            do i = obs_lon_s, obs_lon_e
!               do jj = j, j 
!                  do ii = var_lon_s, var_lon_e

! obs along longitudes and var on the same longitude as obs
         do j = obs_lat_s, obs_lat_e 
            do i = obs_lon_s, obs_lon_e 
               do jj = var_lat_s, var_lat_e 
                  do ii = i, i


                     do itime = 1, dimtime_length
                        ! computation for each group
                        do ig = 1, group_size

                           igs = (ig-1)*ens_copy + 1
                           ige = ig*ens_copy

                           ! get obs information ready
                           yfmean     = bgrid%variables(i,j,k,igs,itime,idobs)
                           yprior_var = bgrid%variables(i,j,k,igs+1,itime,idobs)**2
                           yprior(1:ens_size) = bgrid%variables(i,j,k,igs+1+1:igs+1+ens_size,itime,idobs)
                           if ( ifELF == 1 ) then
                              ytrue      = bgrid%variables(i,j,k,bgrid%ncopy,itime,idobs)
                              Robs       = obs%obsvar(idobs)
                              ypost_var  = 1.0_r8/(1.0_r8/yprior_var + 1.0_r8/Robs)

                              if ( ifINFCOR == 1 ) then
                                 yprior_var_infcor = yprior_var * bgrid%var_inf_cor(i,j,k,ig,itime,idobs)
                                 ypost_var_infcor  = 1.0_r8/(1.0_r8/yprior_var_infcor + 1.0_r8/Robs)
                              endif
                           endif

                           ! get [lon, lat] and convert to radiance
                           obslon = bgrid%var_lon(i,idobs) * DEG2RAD
                           obslat = bgrid%var_lat(j,idobs) * DEG2RAD

                           ! get target variable information ready
                           xfmean     = bgrid%variables(ii,jj,kk,igs,itime,inddim)
                           xprior_var = bgrid%variables(ii,jj,kk,igs+1,itime,inddim)**2
                           xprior(1:ens_size) = bgrid%variables(ii,jj,kk,igs+1+1:igs+1+ens_size,itime,inddim)
                           if ( ifELF == 1 ) then
                              xtrue    = bgrid%variables(ii,jj,kk,bgrid%ncopy,itime,inddim)

                              if ( ifINFCOR == 1 ) then
                                 xprior_var_infcor = xprior_var * bgrid%var_inf_cor(ii,jj,kk,ig,itime,inddim)
                              endif
                           endif

                           ! get [lon, lat] and convert to radiance
                           varlon = bgrid%var_lon(ii,inddim) * DEG2RAD
                           varlat = bgrid%var_lat(jj,inddim) * DEG2RAD

                           ! get the horizontal distance
                           dist = get_horiz_dist(varlon, varlat, obslon, obslat)

                           call find_index_in_locind(ilocind, dist, locnum, localization%locind)

                           ! localization computation
                           covxy = comp_cov(ens_size,xfmean,xprior,yfmean,yprior)
                           corr(ig) = covxy / (sqrt(yprior_var) * sqrt(xprior_var))
                           reg_coef(ig) = covxy / yprior_var

                           if ( ifINFCOR == 1 ) then
                              covxy_infcor = covxy * sqrt(bgrid%var_inf_cor(i,j,k,ig,itime,idobs)) *      &
                                             sqrt(bgrid%var_inf_cor(ii,jj,kk,ig,itime,inddim))
                              reg_coef_infcor(ig) = covxy_infcor / yprior_var_infcor
                           endif

                           if ( ifELF == 1 ) then
                              xt_xfmean = xtrue - xfmean
                              ccc(ig) = reg_coef(ig) * (ypost_var / yprior_var - 1.0_r8) * yfmean
                              ddd(ig) = reg_coef(ig) * ypost_var / Robs
!                              sublocx(ig)  = xt_xfmean
!                              sublocy(ig)  = ccc + ddd * ytrue
!                              sublocx2(ig) = xt_xfmean**2
!                              sublocy2(ig) = (ccc + ddd*ytrue)**2
                              sublocnumer(ig) = xt_xfmean * (ccc(ig) + ddd(ig) * ytrue)
                              sublocdenom(ig) = (ccc(ig) + ddd(ig)*ytrue)**2 + ddd(ig)**2 * Robs

                              incnumer(ig) = ccc(ig) + ddd(ig)*ytrue
                              incdenom(ig) = (ccc(ig) + ddd(ig)*ytrue)**2

                              if ( ifINFCOR == 1 ) then
                                 ccc_infcor(ig) = reg_coef_infcor(ig) * (ypost_var_infcor / yprior_var_infcor - 1.0_r8) * yfmean
                                 ddd_infcor(ig) = reg_coef_infcor(ig) * ypost_var_infcor / Robs
                                 sublocnumer_infcor(ig) = xt_xfmean * (ccc_infcor(ig) + ddd_infcor(ig) * ytrue)
                                 sublocdenom_infcor(ig) = (ccc_infcor(ig) + ddd_infcor(ig)*ytrue)**2 + ddd_infcor(ig)**2 * Robs 

                                 incnumer_infcor(ig) = ccc_infcor(ig) + ddd_infcor(ig) * ytrue
!                                 incdenom_infcor(ig) = (ccc_infcor(ig) + ddd_infcor(ig)*ytrue)**2
                              endif
                           endif

                        enddo   !ig 

                        ! record in distance
                        localization%num(ilocind,indvars) =                    &
                            localization%num(ilocind,indvars) + 1

                        sumbeta1 = (sum(reg_coef))**2
                        sumbeta2 = 0.0_r8
                        do igg = 1, group_size
                           sumbeta2 = sumbeta2 + reg_coef(igg)**2
                        enddo
                        localization%beta1(ilocind,indvars) =                  &
                            localization%beta1(ilocind,indvars) + sumbeta1
                        localization%beta2(ilocind,indvars) =                  &
                            localization%beta2(ilocind,indvars) + sumbeta2

                        localization%meangf(ilocind,indvars) =                 &
                            localization%meangf(ilocind,indvars) + (sumbeta1/sumbeta2-1.0_r8)/(group_size-1.0_r8)

!print *, 'reg ', reg_coef(:)
!print *, 'regi', reg_coef_infcor(:)
!                       correlation space given one subset
                        do ig = 1, group_size
                           indcorr = floor((corr(ig)+1.0_r8+drho/2.0_r8)/drho + 1.0_r8);
                           if ( indcorr < 1 .or. indcorr > nrho ) then
                              print *, 'ERROR: wrong index for correlation ', ig, indcorr, corr(ig)
                           endif

                           localization%num_rf(ilocind,indcorr,ig,indvars) =            &
                               localization%num_rf(ilocind,indcorr,ig,indvars) + 1

                           sumbeta1 = reg_coef(ig)
                           indk2 = 0
                           do indk1 = 1, group_size
                              if ( indk1 /= ig ) then
                                 indk2 = indk2 + 1
                                 subset_regcoef(indk2) = reg_coef(indk1)
                              endif
                           enddo
 
                           localization%gbetasq(ilocind,indcorr,ig,indvars) =           &
                               localization%gbetasq(ilocind,indcorr,ig,indvars) + sum(subset_regcoef**2)
                           localization%tgbeta(ilocind,indcorr,ig,indvars) =            &
                               localization%tgbeta(ilocind,indcorr,ig,indvars) + sum(sumbeta1*subset_regcoef)

!print *, 'ig = ', ig
!print *, sumbeta1, subset_regcoef
!print *, sum(subset_regcoef**2), sum(sumbeta1*subset_regcoef)

                           if ( ifELF == 1 .and. ifINFCOR == 1 ) then

                              sumbeta1 = reg_coef_infcor(ig)
                              indk2 = 0
                              do indk1 = 1, group_size
                                 if ( indk1 /= ig ) then
                                    indk2 = indk2 + 1
                                    subset_regcoef(indk2) = reg_coef_infcor(indk1)
                                 endif
                              enddo

                              localization%gbetasq_infcor(ilocind,indcorr,ig,indvars) =           &
                                  localization%gbetasq_infcor(ilocind,indcorr,ig,indvars) + sum(subset_regcoef**2)
                              localization%tgbeta_infcor(ilocind,indcorr,ig,indvars) =            &
                                  localization%tgbeta_infcor(ilocind,indcorr,ig,indvars) + sum(sumbeta1*subset_regcoef)

                           endif

!print *, sumbeta1, subset_regcoef
!print *, sum(subset_regcoef**2), sum(sumbeta1*subset_regcoef)

                        enddo
!pause

!                        incsumbeta1 = (sum(incnumer))**2 + Robs * (sum(ddd))**2
!                        incsumbeta2 = 0.0_r8
!                        do igg = 1, group_size
!                           incsumbeta2 = incsumbeta2 + incdenom(igg) + Robs * ddd(igg)**2
!                        enddo
!                        localization%incbeta1(ilocind,indvars) =               &
!                            localization%incbeta1(ilocind,indvars) + incsumbeta1
!                        localization%incbeta2(ilocind,indvars) =               &
!                            localization%incbeta2(ilocind,indvars) + incsumbeta2
!
!                        localization%incmeangf(ilocind,indvars) =                 &
!                            localization%incmeangf(ilocind,indvars) + (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                        if ( ifELF == 1 ) then
                        indk2 = 1
                        do igg = 1, group_size
                           if ( abs(sublocdenom(igg)) < 0.0000000001_r8 ) then
                               indk2 = -1
                           endif
                        enddo

                        if ( indk2 > 0 ) then

                        incsumbeta1 = (sum(incnumer))**2 + Robs * (sum(ddd))**2
!                        incsumbeta2 = 0.0_r8
!                        do igg = 1, group_size
!                           incsumbeta2 = incsumbeta2 + incdenom(igg) + Robs * ddd(igg)**2
!                        enddo
                        incsumbeta2 = sum(incnumer*incnumer) + Robs * sum(ddd*ddd)

                        localization%incbeta1(ilocind,indvars) =               &
                            localization%incbeta1(ilocind,indvars) + incsumbeta1
                        localization%incbeta2(ilocind,indvars) =               &
                            localization%incbeta2(ilocind,indvars) + incsumbeta2

                        localization%incmeangf(ilocind,indvars) =                 &
                            localization%incmeangf(ilocind,indvars) + (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                              sublocalpha(ig) = sublocnumer(ig)/sublocdenom(ig)
                           enddo
                           elfsumbeta1 = (sum(sublocalpha))**2
                           elfsumbeta2 = sum(sublocalpha*sublocalpha)

                           localization%elfbeta1(ilocind,indvars) =                  &
                               localization%elfbeta1(ilocind,indvars) + elfsumbeta1
                           localization%elfbeta2(ilocind,indvars) =                  &
                               localization%elfbeta2(ilocind,indvars) + elfsumbeta2

                           localization%elfmeangf(ilocind,indvars) =                 &
                               localization%elfmeangf(ilocind,indvars) + (elfsumbeta1/elfsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                           localization%num2(ilocind,ig,indvars) =                            &
                               localization%num2(ilocind,ig,indvars) + 1
!                           localization%sumx(ilocind,ig,ndims,inddim) =                            &
!                               localization%sumx(ilocind,ig,ndims,inddim) + sublocx(ig)
!                           localization%sumy(ilocind,ig,ndims,inddim) =                            &
!                               localization%sumy(ilocind,ig,ndims,inddim) + sublocy(ig)
!
!                           localization%sumx2(ilocind,ig,ndims,inddim) =                           &
!                               localization%sumx2(ilocind,ig,ndims,inddim) + sublocx2(ig)
!                           localization%sumy2(ilocind,ig,ndims,inddim) =                           &
!                               localization%sumy2(ilocind,ig,ndims,inddim) + sublocy2(ig)

                           localization%numerator(ilocind,ig,indvars) =                       &
                               localization%numerator(ilocind,ig,indvars) + sublocnumer(ig)
                           localization%denominator(ilocind,ig,indvars) =                     &
                               localization%denominator(ilocind,ig,indvars) + sublocdenom(ig)

!                           localization%alpha1(ilocind,ig,ndims,inddim) =                          &
!                               localization%alpha1(ilocind,ig,ndims,inddim) + sublocalpha(ig)
!                           localization%alpha2(ilocind,ig,ndims,inddim) =                          &
!                               localization%alpha2(ilocind,ig,ndims,inddim) + sublocalpha(ig)**2

                           localization%correlation(ilocind,ig,indvars) =                     &
                               localization%correlation(ilocind,ig,indvars) + abs(corr(ig))
                           localization%sumcorr(ilocind,ig,indvars) =                         &
                               localization%sumcorr(ilocind,ig,indvars) + corr(ig)
                           localization%sumcorrsq(ilocind,ig,indvars) =                       &
                               localization%sumcorrsq(ilocind,ig,indvars) + corr(ig)**2
                           enddo   ! ig
                        endif   ! indk2 > 0
                        endif   ! ifELF

                        if ( ifELF == 1 .and. ifINFCOR == 1 ) then

                           localization%num_infcor(ilocind,indvars) =                    &
                               localization%num_infcor(ilocind,indvars) + 1

                           sumbeta1 = (sum(reg_coef_infcor))**2
                           sumbeta2 = 0.0_r8
                           do igg = 1, group_size
                              sumbeta2 = sumbeta2 + reg_coef_infcor(igg)**2
                           enddo
                           localization%beta1_infcor(ilocind,indvars) =                  &
                               localization%beta1_infcor(ilocind,indvars) + sumbeta1
                           localization%beta2_infcor(ilocind,indvars) =                  &
                               localization%beta2_infcor(ilocind,indvars) + sumbeta2

                           localization%meangf_infcor(ilocind,indvars) =                 &
                               localization%meangf_infcor(ilocind,indvars) + (sumbeta1/sumbeta2-1.0_r8)/(group_size-1.0_r8)

!                           incsumbeta1 = (sum(incnumer_infcor))**2 + Robs * (sum(ddd_infcor))**2
!                           incsumbeta2 = 0.0_r8
!                           do igg = 1, group_size
!                              incsumbeta2 = incsumbeta2 + incdenom_infcor(igg) + Robs * ddd_infcor(igg)**2
!                           enddo
!                           localization%incbeta1_infcor(ilocind,indvars) =               &
!                               localization%incbeta1_infcor(ilocind,indvars) + incsumbeta1
!                           localization%incbeta2_infcor(ilocind,indvars) =               &
!                               localization%incbeta2_infcor(ilocind,indvars) + incsumbeta2
!   
!                           localization%incmeangf_infcor(ilocind,indvars) =                 &
!                               localization%incmeangf_infcor(ilocind,indvars) + (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           if ( indk2 > 0 ) then

                           incsumbeta1 = (sum(incnumer_infcor))**2 + Robs * (sum(ddd_infcor))**2
!                           incsumbeta2 = 0.0_r8
!                           do igg = 1, group_size
!                              incsumbeta2 = incsumbeta2 + incdenom_infcor(igg) + Robs * ddd_infcor(igg)**2
!                           enddo
                           incsumbeta2 = sum(incnumer_infcor*incnumer_infcor) + Robs * sum(ddd_infcor*ddd_infcor)
                           localization%incbeta1_infcor(ilocind,indvars) =               &
                               localization%incbeta1_infcor(ilocind,indvars) + incsumbeta1
                           localization%incbeta2_infcor(ilocind,indvars) =               &
                               localization%incbeta2_infcor(ilocind,indvars) + incsumbeta2
   
                           localization%incmeangf_infcor(ilocind,indvars) =                 &
                               localization%incmeangf_infcor(ilocind,indvars) + (incsumbeta1/incsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                              sublocalpha_infcor(ig) = sublocnumer_infcor(ig)/sublocdenom_infcor(ig)
                           enddo
                           elfsumbeta1 = (sum(sublocalpha_infcor))**2
                           elfsumbeta2 = sum(sublocalpha_infcor*sublocalpha_infcor)

                           localization%elfbeta1_infcor(ilocind,indvars) =                  &
                               localization%elfbeta1_infcor(ilocind,indvars) + elfsumbeta1
                           localization%elfbeta2_infcor(ilocind,indvars) =                  &
                               localization%elfbeta2_infcor(ilocind,indvars) + elfsumbeta2

                           localization%elfmeangf_infcor(ilocind,indvars) =                 &
                               localization%elfmeangf_infcor(ilocind,indvars) + (elfsumbeta1/elfsumbeta2-1.0_r8)/(group_size-1.0_r8)

                           do ig = 1, group_size
                           localization%num2_infcor(ilocind,ig,indvars) =                            &
                               localization%num2_infcor(ilocind,ig,indvars) + 1

                           localization%numerator_infcor(ilocind,ig,indvars) =                       &
                               localization%numerator_infcor(ilocind,ig,indvars) + sublocnumer_infcor(ig)
                           localization%denominator_infcor(ilocind,ig,indvars) =                     &
                               localization%denominator_infcor(ilocind,ig,indvars) + sublocdenom_infcor(ig)

                           localization%correlation_infcor(ilocind,ig,indvars) =                     &
                               localization%correlation_infcor(ilocind,ig,indvars) + abs(corr(ig))
                           localization%sumcorr_infcor(ilocind,ig,indvars) =                         &
                               localization%sumcorr_infcor(ilocind,ig,indvars) + corr(ig)
                           localization%sumcorrsq_infcor(ilocind,ig,indvars) =                       &
                               localization%sumcorrsq_infcor(ilocind,ig,indvars) + corr(ig)**2

                           enddo   ! ig
                        endif   ! indk2 > 0
                        endif   ! ifELF and ifINFCOR

                     enddo   ! itime

                  enddo  ! ii
               enddo     ! jj

            enddo  ! i
         enddo     ! j


         deallocate(yprior,xprior)
         deallocate(ccc,ddd)
         deallocate(corr,reg_coef)
         deallocate(sublocx,sublocy,sublocx2,sublocy2,sublocnumer,sublocdenom,sublocalpha)
         deallocate(incnumer,incdenom,incalpha)

         deallocate(ccc_infcor,ddd_infcor)
         deallocate(reg_coef_infcor,sublocnumer_infcor,sublocdenom_infcor,sublocalpha_infcor)
         deallocate(incnumer_infcor,incdenom_infcor,incalpha_infcor)

         deallocate(subset_regcoef)

      enddo LoopVarIndex
!!!$OMP END PARALLEL DO


  
! NOTE: the writting in OMP does not work...
!       thus put the writting out of OMP loop
!   do inddim = 1, bgrid%number_of_state_variables
print *, 'start writing'
   do indvars = 1, ndims

      writelenloc = locnum

      ! find the state variable index and its vertical index
      if ( indvars == 1 ) then
           inddim = 1
           kk = 1
      else
           if ( indvars <= 6 ) then
                inddim = 2
           elseif ( indvars <= 11 ) then
                inddim = 3
           else
                inddim = 4
           endif
           indk1 = mod(indvars-1,5)
           if ( indk1 == 0 ) then
                kk = 5
           else
                kk = indk1
           endif
      endif 

      if ( k < 10 ) then
         write(idchar1,'(i1)') k
      else
         write(idchar1,'(i2)') k
      endif
      if ( kk < 10 ) then
         write(idchar2,'(i1)') kk
      else
         write(idchar2,'(i2)') kk
      endif
      var_name_tmp = trim(obs%obstype_metadata(idobs))//'_'//         &
                     trim(bgrid%description(inddim))//'_obslev'//trim(idchar1)//'_varlev'//trim(idchar2)
      iunit = 10 + inddim
      open(iunit, FILE=var_name_tmp, FORM='FORMATTED')

      allocate(numtmp(writelenloc))
      numtmp(1:writelenloc) = localization%num(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) numtmp
      write(iunit, FMT=locreal_format) localization%beta1(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%beta2(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%meangf(1:writelenloc,indvars)

      write(iunit, FMT=locreal_format) localization%incbeta1(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%incbeta2(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%incmeangf(1:writelenloc,indvars)

      write(iunit, FMT=locreal_format) localization%elfbeta1(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%elfbeta2(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%elfmeangf(1:writelenloc,indvars)

      if ( ifELF == 1 ) then
      do ig = 1, group_size
         numtmp(1:writelenloc) = localization%num2(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) numtmp

!      write(iunit, FMT=locreal_format) localization%sumx(1:writelenloc,ig,ndims,inddim)
!      write(iunit, FMT=locreal_format) localization%sumy(1:writelenloc,ig,ndims,inddim)
!      write(iunit, FMT=locreal_format) localization%sumx2(1:writelenloc,ig,ndims,inddim)
!      write(iunit, FMT=locreal_format) localization%sumy2(1:writelenloc,ig,ndims,inddim)
      write(iunit, FMT=locreal_format) localization%numerator(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%denominator(1:writelenloc,ig,indvars)
!      write(iunit, FMT=locreal_format) localization%alpha1(1:writelenloc,ig,ndims,inddim)
!      write(iunit, FMT=locreal_format) localization%alpha2(1:writelenloc,ig,ndims,inddim)
      write(iunit, FMT=locreal_format) localization%correlation(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%sumcorr(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%sumcorrsq(1:writelenloc,ig,indvars)
      enddo
      endif

      if ( ifELF == 1 .and. ifINFCOR == 1 ) then
      numtmp(1:writelenloc) = localization%num_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) numtmp
      write(iunit, FMT=locreal_format) localization%beta1_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%beta2_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%meangf_infcor(1:writelenloc,indvars)

      write(iunit, FMT=locreal_format) localization%incbeta1_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%incbeta2_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%incmeangf_infcor(1:writelenloc,indvars)

      write(iunit, FMT=locreal_format) localization%elfbeta1_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%elfbeta2_infcor(1:writelenloc,indvars)
      write(iunit, FMT=locreal_format) localization%elfmeangf_infcor(1:writelenloc,indvars)

      if ( ifELF == 1 ) then
      do ig = 1, group_size
         numtmp(1:writelenloc) = localization%num2_infcor(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) numtmp

      write(iunit, FMT=locreal_format) localization%numerator_infcor(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%denominator_infcor(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%correlation_infcor(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%sumcorr_infcor(1:writelenloc,ig,indvars)
      write(iunit, FMT=locreal_format) localization%sumcorrsq_infcor(1:writelenloc,ig,indvars)
      enddo
      endif
      endif

      deallocate(numtmp)

      close(iunit)

      ! write 3D variable
      allocate(numtmp(nrho))

      var_name_tmp = trim(obs%obstype_metadata(idobs))//'_'//         &
                     trim(bgrid%description(inddim))//'_obslev'//trim(idchar1)//'_varlev'//trim(idchar2)//'_correlspace'
      iunit = 10 + inddim
      open(iunit, FILE=var_name_tmp, FORM='FORMATTED')

      do ig = 1, group_size
         do igg = 1, locnum
            numtmp(1:nrho) = localization%num_rf(igg,1:nrho,ig,indvars)
            write(iunit, FMT=locreal_format) numtmp
         enddo
         do igg = 1, locnum
            write(iunit, FMT=locreal_format) localization%gbetasq(igg,1:nrho,ig,indvars)
         enddo
         do igg = 1, locnum
            write(iunit, FMT=locreal_format) localization%tgbeta(igg,1:nrho,ig,indvars)
         enddo
         if ( ifELF == 1 .and. ifINFCOR == 1 ) then
         do igg = 1, locnum
            write(iunit, FMT=locreal_format) localization%gbetasq_infcor(igg,1:nrho,ig,indvars)
         enddo
         do igg = 1, locnum
            write(iunit, FMT=locreal_format) localization%tgbeta_infcor(igg,1:nrho,ig,indvars)
         enddo
         endif
      enddo

      close(iunit)

      deallocate(numtmp)

   enddo

print *, 'end writing'

   deallocate(localization%num, localization%beta1, localization%beta2, localization%meangf)
   deallocate(localization%incbeta1, localization%incbeta2, localization%incmeangf)
   deallocate(localization%elfbeta1,localization%elfbeta2,localization%elfmeangf)
   deallocate(localization%num2, localization%numerator, localization%denominator)
   deallocate(localization%sumx, localization%sumy, localization%sumx2, localization%sumy2)
   deallocate(localization%alpha1, localization%alpha2)
   deallocate(localization%correlation)
   deallocate(localization%sumcorr,localization%sumcorrsq)

   deallocate(localization%num_infcor, localization%beta1_infcor, localization%beta2_infcor, localization%meangf_infcor)
   deallocate(localization%incbeta1_infcor, localization%incbeta2_infcor, localization%incmeangf_infcor)
   deallocate(localization%elfbeta1_infcor,localization%elfbeta2_infcor,localization%elfmeangf_infcor)
   deallocate(localization%num2_infcor, localization%numerator_infcor, localization%denominator_infcor)
   deallocate(localization%sumx_infcor, localization%sumy_infcor, localization%sumx2_infcor, localization%sumy2_infcor)
   deallocate(localization%alpha1_infcor, localization%alpha2_infcor)
   deallocate(localization%correlation_infcor)
   deallocate(localization%sumcorr_infcor,localization%sumcorrsq_infcor)

   deallocate(localization%num_rf,localization%gbetasq,localization%tgbeta)
   deallocate(localization%gbetasq_infcor,localization%tgbeta_infcor)

   enddo LoopObsLevels
 
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
!
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


!subroutine compute_model_height(we, wes, sn, sns, bt, bts, hgt, phb, ph, lon, lat, nmhgt, mhgt)
!! compute the model height for each model height type
!  
!integer, intent(in)          :: we, wes, sn, sns, bt, bts
!real(r8), intent(in)         :: hgt(:, :)
!real(r8), intent(in)         :: phb(:, :, :), ph(:, :, :)
!real(r8), intent(in)         :: lon(:, :), lat(:, :)
!integer, intent(in)          :: nmhgt
!real(r8), intent(out)        :: mhgt(:, :, :, :)
!
!integer      :: i,j,k,ind
!real(r8)     :: geop, lattmp
!
!! mhgtindex: 1 - surface mass
!!            2 - surface U
!!            3 - surface T
!!!            4 - 3d U
!!!            5 - 3d V
!!!            6 - 3d W (vertical half level)
!!!            7 - 3d T (vertical full level, mass)
!!            4 - 3d T
!
!do ind = 1, nmhgt
!   if ( ind == 1 ) then
!        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn)
!   elseif ( ind == 2 ) then
!        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn) + 10.0
!   elseif ( ind == 3 ) then
!        mhgt(1:we, 1:sn, 1, ind) = hgt(1:we, 1:sn) + 2.0
!   elseif ( ind == 4 ) then
!        do i = 1, we
!           do j = 1, sn
!              do k = 1, bt
!                 geop = ( (phb(i,j,k)+ph(i,j,k))  + (phb(i,j,k+1)+ph(i,j,k+1)) ) / (2.0*gravity)
!                 mhgt(i,j,k,ind) = compute_geometric_height(geop, lat(i,j))
!              enddo
!           enddo
!        enddo 
!   else
!        print *, 'ERROR in compute_model_height'
!        print *, 'ERROR: do not support model height type', ind
!   endif
!enddo
!
!return
!end subroutine compute_model_height

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



