module model_mod

! This module provides routines to work with COSMO data
! files in the DART framework
!
! Author: Jan D. Keller
!         Meteorological Institute, University of Bonn, Germany
!         2011-09-15
!

  use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                               rad2deg, deg2rad, PI
  
  use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                               print_time, print_date, set_calendar_type,        &
                               operator(*),  operator(+), operator(-),           &
                               operator(>),  operator(<), operator(/),           &
                               operator(/=), operator(<=)

  use   cosmo_data_mod, only : cosmo_meta,cosmo_hcoord,cosmo_non_state_data,     &
                               get_cosmo_info,get_data_from_binary,              &
                               set_vertical_coords,grib_header_type

  use     location_mod, only : location_type, get_dist, query_location,          &
                               get_close_maxdist_init, get_close_type,           &
                               set_location, get_location, horiz_dist_only,      & 
                               vert_is_undef,    VERTISUNDEF,                    &
                               vert_is_surface,  VERTISSURFACE,                  &
                               vert_is_level,    VERTISLEVEL,                    &
                               vert_is_pressure, VERTISPRESSURE,                 &
                               vert_is_height,   VERTISHEIGHT,                   &
                               get_close_obs_init, get_close_obs

  use    utilities_mod, only : register_module, error_handler,                   &
                               E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                               nc_check, do_output, to_upper,                    &
                               find_namelist_in_file, check_namelist_read,       &
                               open_file, file_exist, find_textfile_dims,        &
                               file_to_text

  use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,                            &
                               KIND_V_WIND_COMPONENT,                            &
                               KIND_VERTICAL_VELOCITY,                           &
                               KIND_TEMPERATURE,                                 &
                               KIND_PRESSURE,                                    &
                               KIND_SPECIFIC_HUMIDITY,                           &
                               KIND_CLOUD_LIQUID_WATER,                          &
                               KIND_CLOUD_ICE,                                   &
                               KIND_SURFACE_ELEVATION,                           &
                               KIND_SURFACE_GEOPOTENTIAL

  use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

  use byte_mod, only: to_float1,from_float1,word_to_byte_data,byte_to_word_signed,word_to_byte

  use netcdf 

  implicit none

  public  :: get_model_size
  public  :: static_init_model
  public  :: get_state_meta_data
  public  :: get_model_time_step
  public  :: model_interpolate
  public  :: init_conditions
  public  :: init_time
  public  :: adv_1step
  public  :: end_model
  public  :: nc_write_model_atts
  public  :: nc_write_model_vars
  public  :: pert_model_state
  public  :: get_close_maxdist_init
  public  :: get_close_obs_init
  public  :: get_close_obs
  public  :: ens_mean_for_model

!  public  :: grib_to_sv
!  public  :: sv_to_grib

  private :: set_allowed_state_vector_vars
  private :: ll_to_xyz_vector
  private :: ll_to_xyz_single
  private :: get_enclosing_grid_box
  private :: bilinear_interpolation
  private :: linear_interpolation
!  public :: linear_interpolation
!  private :: data_to_state_vector
  private :: get_vertical_boundaries

  public  :: get_state_time
  public  :: get_state_vector
  public  :: write_grib_file

  INTERFACE sv_to_field
    MODULE PROCEDURE sv_to_field_2d
    MODULE PROCEDURE sv_to_field_3d
  END INTERFACE

  type dart_variable_info
    character(len=16)    :: varname_short
    character(len=256)   :: varname_long
    character(len=32)    :: units
    logical              :: is_present
    integer              :: nx
    integer              :: ny
    integer              :: nz
    real(r8),allocatable :: vertical_level(:)
    integer              :: vertical_coordinate
    integer              :: horizontal_coordinate
    integer,allocatable  :: state_vector_sindex(:) ! starting index in state vector for every vertical level
    integer,allocatable  :: cosmo_state_index(:)   ! index in cosmo state of every vertical level
  end type dart_variable_info

  integer,parameter              :: n_state_vector_vars=8
  integer,parameter              :: n_non_state_vars=1

  character(len=256)             :: string
  logical, save                  :: module_initialized = .FALSE.

  type(cosmo_meta),allocatable   :: cosmo_vars(:)
  type(cosmo_hcoord)             :: cosmo_lonlat(3)
  integer                        :: nvars

  character(len=256)             :: cosmo_filename
  integer                        :: ivctype
  integer                        :: model_dt
  real(r8)                       :: model_perturbation_amplitude
  logical                        :: output_state_vector

  namelist /model_nml/  &
   cosmo_filename,ivctype,model_dt,model_perturbation_amplitude,output_state_vector

  integer                        :: model_size
  type(time_type)                :: model_timestep  ! smallest time to adv model

  integer, parameter             :: n_max_kinds=200

  integer                        :: allowed_state_vector_vars(n_state_vector_vars)
  integer                        :: allowed_non_state_vars(1:n_max_kinds)
  logical                        :: is_allowed_state_vector_var(n_max_kinds)
  logical                        :: is_allowed_non_state_var(n_max_kinds)
  type(dart_variable_info)       :: state_vector_vars(1:n_max_kinds)
  type(cosmo_non_state_data)     :: non_state_data

  real(r8),allocatable           :: state_vector(:)

  type(random_seq_type) :: random_seq

  character(len=256)    :: error_string,error_string2

  type(time_type)                :: cosmo_fc_time
  type(time_type)                :: cosmo_an_time

  type(grib_header_type),allocatable  :: grib_header(:)

contains

  function get_model_size()

    integer :: get_model_size
    
    if ( .not. module_initialized ) call static_init_model
    
    get_model_size = model_size

  end function get_model_size

  subroutine static_init_model()

    integer                       :: iunit,io,ivar,ikind,sv_length,i
    integer                       :: sidx,eidx
    integer,allocatable           :: pp_index(:)
    real(r8),allocatable          :: data(:,:)

    real(r8),parameter            :: g = 9.80665

    if ( module_initialized ) return ! only need to do this once.
    
    module_initialized=.TRUE.

    call set_allowed_state_vector_vars()

    ! read the DART namelist for this model
    call find_namelist_in_file('input.nml', 'model_nml', iunit)
    read(iunit, nml = model_nml, iostat = io)
    call check_namelist_read(iunit, io, 'model_nml')
    print*,TRIM(cosmo_filename)

    call get_cosmo_info(cosmo_filename,cosmo_vars,cosmo_lonlat,grib_header,&
                        is_allowed_state_vector_var,cosmo_fc_time)

    state_vector_vars(:)%is_present=.false.

    model_size=maxval(cosmo_vars(:)%dart_eindex)
    nvars=size(cosmo_vars,1)

    sv_length=0
    do ivar=1,nvars
      ikind=cosmo_vars(ivar)%dart_kind
      if (is_allowed_state_vector_var(ikind)) then
        sv_length=sv_length+cosmo_vars(ivar)%dims(1)*cosmo_vars(ivar)%dims(2)
      end if
    end do

    allocate(state_vector(1:sv_length))

    do ivar=1,nvars
      ikind=cosmo_vars(ivar)%dart_kind

      if (is_allowed_state_vector_var(ikind)) then
        if (.not. state_vector_vars(ikind)%is_present) then
          state_vector_vars(ikind)%is_present=.true.
          state_vector_vars(ikind)%varname_short=cosmo_vars(ivar)%varname_short
          state_vector_vars(ikind)%varname_long=cosmo_vars(ivar)%varname_long
          state_vector_vars(ikind)%units=cosmo_vars(ivar)%units
          state_vector_vars(ikind)%nx=cosmo_vars(ivar)%dims(1)
          state_vector_vars(ikind)%ny=cosmo_vars(ivar)%dims(2)
          state_vector_vars(ikind)%nz=cosmo_vars(ivar)%dims(3)
          state_vector_vars(ikind)%horizontal_coordinate=cosmo_vars(ivar)%hcoord_type
          if (state_vector_vars(ikind)%nz>1) then
            state_vector_vars(ikind)%vertical_coordinate=VERTISLEVEL
          else
            state_vector_vars(ikind)%vertical_coordinate=VERTISSURFACE
          end if

          allocate(state_vector_vars(ikind)%vertical_level(1:state_vector_vars(ikind)%nz))
          allocate(state_vector_vars(ikind)%state_vector_sindex(1:state_vector_vars(ikind)%nz))
          allocate(state_vector_vars(ikind)%cosmo_state_index(1:state_vector_vars(ikind)%nz))

        end if

        state_vector_vars(ikind)%vertical_level(cosmo_vars(ivar)%ilevel)=cosmo_vars(ivar)%dart_level
        state_vector_vars(ikind)%state_vector_sindex(cosmo_vars(ivar)%ilevel)=cosmo_vars(ivar)%dart_sindex
        state_vector_vars(ikind)%cosmo_state_index(cosmo_vars(ivar)%ilevel)=ivar

      end if

      if (is_allowed_non_state_var(ikind)) then
        if (ikind==KIND_SURFACE_ELEVATION) then
          allocate(data(1:cosmo_vars(ivar)%dims(1),1:cosmo_vars(ivar)%dims(2)))
          data=get_data_from_binary(cosmo_filename,grib_header(ivar),cosmo_vars(ivar)%dims(1),cosmo_vars(ivar)%dims(2))
          if (.not. allocated(non_state_data%surface_orography)) then
            allocate(non_state_data%surface_orography(1:cosmo_vars(ivar)%dims(1),1:cosmo_vars(ivar)%dims(2)))
          end if
          non_state_data%surface_orography(:,:)=data(:,:)
          deallocate(data)
        end if
        if ((ikind==KIND_SURFACE_GEOPOTENTIAL).and.(.not. allocated(non_state_data%surface_orography))) then
          allocate(data(1:cosmo_vars(ivar)%dims(1),1:cosmo_vars(ivar)%dims(2)))
          data=get_data_from_binary(cosmo_filename,grib_header(ivar),cosmo_vars(ivar)%dims(1),cosmo_vars(ivar)%dims(2))
          allocate(non_state_data%surface_orography(1:cosmo_vars(ivar)%dims(1),1:cosmo_vars(ivar)%dims(2)))
          non_state_data%surface_orography(:,:)=data(:,:)/g
          deallocate(data)
        end if
      end if
    end do

    setlevel : do ivar=1,nvars
      if (cosmo_vars(ivar)%dart_kind==KIND_U_WIND_COMPONENT) then

        if (state_vector_vars(KIND_PRESSURE)%is_present) then
          allocate(pp_index(1:state_vector_vars(KIND_PRESSURE)%nz))
          pp_index(:)=state_vector_vars(KIND_PRESSURE)%state_vector_sindex(:)
        else
          allocate(pp_index(1:1))
          pp_index(1)=-1
        end if
        call set_vertical_coords(cosmo_filename,grib_header(ivar),non_state_data,ivctype,state_vector,pp_index)

        exit setlevel
      end if
    end do setlevel

    return
  end subroutine static_init_model

  subroutine get_state_meta_data(index_in,location,var_type)

    integer, intent(in)            :: index_in
    type(location_type)            :: location
    integer, optional, intent(out) :: var_type

    integer                        :: ivar,var,varindex,hindex,dims(3)
    real(r8)                       :: lon,lat,vloc

    if (.NOT. module_initialized) CALL static_init_model()
    
    var=-1

    findindex : DO ivar=1,nvars
      IF ((index_in >= cosmo_vars(ivar)%dart_sindex) .AND. (index_in <= cosmo_vars(ivar)%dart_eindex)) THEN
        var=ivar
        varindex=index_in-cosmo_vars(ivar)%dart_sindex+1
        dims=cosmo_vars(ivar)%dims
        hindex=MOD(varindex,(dims(1)*dims(2)))
        vloc=cosmo_vars(ivar)%dart_level
        lon=cosmo_lonlat(cosmo_vars(ivar)%hcoord_type)%lon(hindex)
        lat=cosmo_lonlat(cosmo_vars(ivar)%hcoord_type)%lat(hindex)

        location=set_location(lon,lat,vloc,VERTISLEVEL)
        var_type=cosmo_vars(ivar)%dart_kind
        EXIT findindex
      END IF
    END DO findindex

    IF( var == -1 ) THEN
      write(string,*) 'Problem, cannot find base_offset, index_in is: ', index_in
!      call error_handler(E_ERR,'get_state_meta_data',string,source,revision,revdate)
    ENDIF

  end subroutine get_state_meta_data

  function get_model_time_step()

    type(time_type) :: get_model_time_step
    
    if ( .not. module_initialized ) call static_init_model
    
    model_timestep=set_time(model_dt)
    get_model_time_step = model_timestep
    return

  end function get_model_time_step

  subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

    ! Passed variables
    
    real(r8),            intent(in)  :: x(:)
    type(location_type), intent(in)  :: location
    integer,             intent(in)  :: obs_type
    real(r8),            intent(out) :: interp_val
    integer,             intent(out) :: istatus
    
    ! Local storage
    
    real(r8),allocatable :: xyz_grid(:,:)
    real(r8)             :: point_coords(1:3)
    real(r8)             :: xyz_point(3)
    real(r8)             :: lo1,lo2,la1,la2
    real(r8)             :: m1,m2,n1,n2,xc,yc
    
    integer              :: i,j,hbox(2,2),n,vbound(2),sindex,eindex
    real(r8)             :: hbox_weight(2,2),hbox_val(2,2),hbox_lon(2,2),hbox_lat(2,2)
    real(r8)             :: vbound_weight(2),val1,val2
    real(r8),allocatable :: hgrid_data(:)

    ! Error codes:
    ! istatus = 99 : unknown error
    ! istatus = 10 : observation type is not in state vector
    ! istatus = 15 : observation lies outside the model domain (horizontal)
    ! istatus = 16 : observation lies outside the model domain (vertical)
    ! istatus = 19 : observation vertical coordinate is not supported

    IF ( .not. module_initialized ) call static_init_model
    
    interp_val = MISSING_R8     ! the DART bad value flag
    istatus = 99                ! unknown error
    
    if (state_vector_vars(obs_type)%is_present) then

      ! horizontal interpolation

      n=size(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,1)

      allocate(xyz_grid(1:n,1:3))
      point_coords(1:3)=get_location(location)

      ! calculate the angles (in reference to lon/lat) between the desired location and all horizontal grid points
      xyz_grid=ll_to_xyz_vector(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,&
                        cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat)

      xyz_point=ll_to_xyz_single(point_coords(1),point_coords(2))

      ! Find grid indices of box enclosing the observation location

!      call get_enclosing_grid_box(xyz_point,xyz_grid,n,state_vector_vars(obs_type)%nx,state_vector_vars(obs_type)%ny,hbox,hbox_weight)
      call get_enclosing_grid_box_lonlat(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,&
                                         cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat,&
                                         point_coords(1:2),n,state_vector_vars(obs_type)%nx,state_vector_vars(obs_type)%ny,hbox,hbox_weight)
      
!      print*,cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon(hbox(1,1))
!      print*,cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat(hbox(1,1))

      if (hbox(1,1)==-1) then
        istatus=15
        return
      end if

      ! determine vertical level above and below obsevation
      call get_vertical_boundaries(hbox,hbox_weight,obs_type,query_location(location,'which_vert'),point_coords(3),vbound,vbound_weight,istatus)

      ! check if observation is in vertical domain and vertical coordinate system is supported
      if (vbound(1)==-1) then
        return
      end if
      
      ! Perform a bilinear interpolation from the grid box to the desired location
      ! for the level above and below the observation

      sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(1))

      do i=1,2
        do j=1,2
          hbox_val(i,j)=x(sindex+hbox(i,j)-1)
          hbox_lon(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon(hbox(i,j))
          hbox_lat(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat(hbox(i,j))
        end do
      end do

      call bilinear_interpolation(hbox_val,hbox_weight,hbox_lon,hbox_lat,point_coords,val1,istatus)

      sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(2))
      do i=1,2
        do j=1,2
          hbox_val(i,j)=x(sindex+hbox(i,j)-1)
        end do
      end do

      call bilinear_interpolation(hbox_val,hbox_weight,hbox_lon,hbox_lat,point_coords,val2,istatus)

      ! vertical interpolation of horizontally interpolated values

      interp_val=val1*vbound_weight(1)+val2*vbound_weight(2)
      istatus=0
      
      return

    else
      istatus=10
      return
    end if

  end subroutine model_interpolate

  subroutine init_conditions(x)
    !------------------------------------------------------------------
    !
    ! Returns a model state vector, x, that is some sort of appropriate
    ! initial condition for starting up a long integration of the model.
    ! At present, this is only used if the namelist parameter 
    ! start_from_restart is set to .false. in the program perfect_model_obs.
    ! If this option is not to be used in perfect_model_obs, or if no 
    ! synthetic data experiments using perfect_model_obs are planned, 
    ! this can be a NULL INTERFACE.
    
    real(r8), intent(out) :: x(:)
    
    if ( .not. module_initialized ) call static_init_model
    
    x = 0.0_r8
    
  end subroutine init_conditions

  subroutine init_time(time)
    type(time_type), intent(out) :: time

    time=set_time(0,0)

    return
  end subroutine init_time
  
  subroutine adv_1step(x, time)

    ! As COSMO can only be advanced as a separate executable,
    ! this is a NULL INTERFACE.
    
    real(r8),        intent(inout) :: x(:)
    type(time_type), intent(in)    :: time
    
    if ( .not. module_initialized ) call static_init_model
    
    if (do_output()) then
      call print_time(time,'NULL interface adv_1step (no advance) DART time is')
      call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
    endif

    return
    
  end subroutine adv_1step
  
  subroutine end_model()

    deallocate(cosmo_vars)
    deallocate(state_vector)

    return
  end subroutine end_model

  function nc_write_model_atts( ncFileID ) result (ierr)
    
    integer, intent(in)  :: ncFileID      ! netCDF file identifier
    integer              :: ierr          ! return value of function


    integer              :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
    integer              :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
    integer              :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
    integer              :: LineLenDimID
    integer              :: StateVarVarID,StateVarID,VarID
    integer              :: ikind,ndims,idim,dims(100),nx,ny,nz,i
    character(len=3)     :: ckind

    integer              :: lonVarID, latVarID, ulonVarID, ulatVarID, vlonVarID, vlatVarID
    integer              :: levVarID, wlevVarID

    character(len=128)   :: filename
    real(r8)             :: levs(1:500),wlevs(1:501)
   
    if ( .not. module_initialized ) call static_init_model

    ierr = -1 ! assume things go poorly

    write(filename,*) 'ncFileID', ncFileID

    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file, 
    ! and then put into define mode.
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                       'nc_write_model_atts', 'inquire '//trim(filename))
    call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! We need the dimension ID for the number of copies/ensemble members, and
    ! we might as well check to make sure that Time is the Unlimited dimension. 
    ! Our job is create the 'model size' dimension.
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
     'nc_write_model_atts','inq_dimid NMLlinelen')
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
     'nc_write_model_atts', 'copy dimid '//trim(filename))
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
     'nc_write_model_atts', 'time dimid '//trim(filename))
    
    if ( TimeDimID /= unlimitedDimId ) then
      write(error_string,*)'Time Dimension ID ',TimeDimID, &
       ' should equal Unlimited Dimension ID',unlimitedDimID
!      call error_handler(E_ERR,'nc_write_model_atts', error_string, source, revision, revdate)
    endif

    !-------------------------------------------------------------------------------
    ! Define the model size / state variable dimension / whatever ...
    !-------------------------------------------------------------------------------
    call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
     dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! Write Global Attributes 
    !-------------------------------------------------------------------------------

!    call DATE_AND_TIME(crdate,crtime,crzone,values)
!    write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
!     values(1), values(2), values(3), values(5), values(6), values(7)
    
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,str1    ), &
!                  'nc_write_model_atts', 'creation put '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
!                  'nc_write_model_atts', 'source put '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
!                  'nc_write_model_atts', 'revision put '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
!                  'nc_write_model_atts', 'revdate put '//trim(filename))
!    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'cosmo' ), &
!                  'nc_write_model_atts', 'model put '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! Here is the extensible part. The simplest scenario is to output the state vector,
    ! parsing the state vector into model-specific parts is complicated, and you need
    ! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
    ! complicated part.
    !-------------------------------------------------------------------------------
    
    if ( output_state_vector ) then
      
      !----------------------------------------------------------------------------
      ! Create a variable for the state vector
      !----------------------------------------------------------------------------
      
      ! Define the state vector coordinate variable and some attributes.
      call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                    dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                    'statevariable def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                    'nc_write_model_atts','statevariable long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                    'nc_write_model_atts', 'statevariable units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                    'nc_write_model_atts', 'statevariable valid_range '//trim(filename))
      
      ! Define the actual (3D) state vector, which gets filled as time goes on ... 
      call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                    dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                    'nc_write_model_atts','state def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                    'nc_write_model_atts', 'state long_name '//trim(filename))
                   
      ! Leave define mode so we can fill the coordinate variable.
      call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))
      
      ! Fill the state variable coordinate variable
      call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                    'nc_write_model_atts', 'state put_var '//trim(filename))
      
    else
      
      !----------------------------------------------------------------------------
      ! We need to output the prognostic variables.
      !----------------------------------------------------------------------------
      ! Define the new dimensions IDs
      !----------------------------------------------------------------------------

      
      findnxny : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          nx=state_vector_vars(ikind)%nx
          ny=state_vector_vars(ikind)%ny
          exit findnxny
        end if
      end do findnxny

      findnz : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          if ((state_vector_vars(ikind)%nz>1) .and. (ikind .ne. KIND_VERTICAL_VELOCITY)) then
            nz=state_vector_vars(ikind)%nz
            exit findnz
          end if
        end if
      end do findnz

      ! Standard Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='LON', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=lonVarID),&
                    'nc_write_model_atts', 'LON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitudes of grid'), &
                    'nc_write_model_atts', 'LON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'LON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LON valid_range '//trim(filename))

      ! U Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='ULON', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=ulonVarID),&
                    'nc_write_model_atts', 'ULON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitudes for U-wind'), &
                    'nc_write_model_atts', 'ULON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULON valid_range '//trim(filename))

      ! V Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='VLON', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=vlonVarID),&
                    'nc_write_model_atts', 'VLON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'long_name', 'longitudes for V-wind'), &
                    'nc_write_model_atts', 'VLON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'VLON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLON valid_range '//trim(filename))
      
      ! Standard Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='LAT', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=latVarID),&
                    'nc_write_model_atts', 'LAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitudes of grid'), &
                    'nc_write_model_atts', 'LAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'LAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LAT valid_range '//trim(filename))

      ! U Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='ULAT', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=ulatVarID),&
                    'nc_write_model_atts', 'ULAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitudes for U-wind'), &
                    'nc_write_model_atts', 'ULAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULAT valid_range '//trim(filename))

      ! V Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='VLAT', xtype=nf90_real, &
       dimids=(/ nx*ny /), varid=vlatVarID),&
                    'nc_write_model_atts', 'VLAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'long_name', 'latitudes for V-wind'), &
                    'nc_write_model_atts', 'VLAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'VLAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLAT valid_range '//trim(filename))

      ! Standard Z Levels
      call nc_check(nf90_def_var(ncFileID,name='LEV', xtype=nf90_real, &
       dimids=(/ nz /), varid=levVarID),&
                    'nc_write_model_atts', 'LEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'long_name', 'standard hybrid model levels'), &
                    'nc_write_model_atts', 'LEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'LEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'LEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'LEV valid_range '//trim(filename))

      ! W-wind Z Levels
      call nc_check(nf90_def_var(ncFileID,name='WLEV', xtype=nf90_real, &
       dimids=(/ nz+1 /), varid=wlevVarID),&
                    'nc_write_model_atts', 'WLEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'long_name', 'standard model levels for W-wind'), &
                    'nc_write_model_atts', 'WLEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'WLEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'WLEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'WLEV valid_range '//trim(filename))

      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then

          error_string = trim(filename)//' '//trim(state_vector_vars(ikind)%varname_short)

          dims(1)=lonVarId
          dims(2)=latVarId
          if (ikind==KIND_U_WIND_COMPONENT) then
            dims(1)=ulonVarId
            dims(2)=ulatVarId
          end if
          if (ikind==KIND_V_WIND_COMPONENT) then
            dims(1)=vlonVarId
            dims(2)=vlatVarId
          end if

          idim=3
          if (state_vector_vars(ikind)%nz>1) then
            dims(idim)=levVarId
            if (ikind==KIND_VERTICAL_VELOCITY) then
              wlevs(1:nz+1)=state_vector_vars(ikind)%vertical_level(1:nz+1)
              dims(idim)=wlevVarId
            else
              levs(1:nz)=state_vector_vars(ikind)%vertical_level(1:nz)
            end if
            idim=idim+1
          end if
          ! Put ensemble member dimension here
          dims(idim)=unlimitedDimId
          ndims=idim

          call nc_check(nf90_def_var(ncid=ncFileID, name=trim(state_vector_vars(ikind)%varname_short), xtype=nf90_real, &
                        dimids = dims(1:ndims), varid=VarID),&
                        'nc_write_model_atts', trim(error_string)//' def_var' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(state_vector_vars(ikind)%varname_long)), &
                        'nc_write_model_atts', trim(error_string)//' put_att long_name' )

          write(ckind,'(I)') ikind
          call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(ckind)), &
                        'nc_write_model_atts', trim(error_string)//' put_att dart_kind' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(state_vector_vars(ikind)%units)), &
                        'nc_write_model_atts', trim(error_string)//' put_att units' )

        end if
      end do

      !----------------------------------------------------------------------------
      ! Fill the coordinate variables
      !----------------------------------------------------------------------------
      
      call nc_check(nf90_put_var(ncFileID, lonVarID, cosmo_lonlat(1)%lon ), &
                    'nc_write_model_atts', 'LON put_var '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, latVarID, cosmo_lonlat(1)%lat ), &
                    'nc_write_model_atts', 'LAT put_var '//trim(filename))
      
      call nc_check(nf90_put_var(ncFileID, ulonVarID, cosmo_lonlat(2)%lon ), &
                    'nc_write_model_atts', 'ULON put_var '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, ulatVarID, cosmo_lonlat(2)%lat ), &
                    'nc_write_model_atts', 'ULAT put_var '//trim(filename))
      
      call nc_check(nf90_put_var(ncFileID, vlonVarID, cosmo_lonlat(3)%lon ), &
                    'nc_write_model_atts', 'VLON put_var '//trim(filename))
      call nc_check(nf90_put_var(ncFileID, vlatVarID, cosmo_lonlat(3)%lat ), &
                    'nc_write_model_atts', 'VLAT put_var '//trim(filename))

      call nc_check(nf90_put_var(ncFileID, levVarID, levs(1:nz) ), &
                    'nc_write_model_atts', 'LEV put_var '//trim(filename))

      call nc_check(nf90_put_var(ncFileID, wlevVarID, wlevs(1:nz+1) ), &
                    'nc_write_model_atts', 'WLEV put_var '//trim(filename))
      
    end if

    !-------------------------------------------------------------------------------
    ! Flush the buffer and leave netCDF file open
    !-------------------------------------------------------------------------------
    call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')
    
    ierr = 0 ! If we got here, things went well.
    
    return
  end function nc_write_model_atts
  
  function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)         
    !------------------------------------------------------------------
    ! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
    !
    ! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
    ! return code is always '0 == normal', since the fatal errors stop execution.
    !
    ! For the lorenz_96 model, each state variable is at a separate location.
    ! that's all the model-specific attributes I can think of ...
    !
    ! assim_model_mod:init_diag_output uses information from the location_mod
    !     to define the location dimension and variable ID. All we need to do
    !     is query, verify, and fill ...
    !
    ! Typical sequence for adding new dimensions,variables,attributes:
    ! NF90_OPEN             ! open existing netCDF dataset
    !    NF90_redef         ! put into define mode
    !    NF90_def_dim       ! define additional dimensions (if any)
    !    NF90_def_var       ! define variables: from name, type, and dims
    !    NF90_put_att       ! assign attribute values
    ! NF90_ENDDEF           ! end definitions: leave define mode
    !    NF90_put_var       ! provide values for variable
    ! NF90_CLOSE            ! close: save updated netCDF dataset
    
    integer,                intent(in) :: ncFileID      ! netCDF file identifier
    real(r8), dimension(:), intent(in) :: state_vec
    integer,                intent(in) :: copyindex
    integer,                intent(in) :: timeindex
    integer                            :: ierr          ! return value of function
    
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
    character(len=NF90_MAX_NAME)          :: varname 
    integer :: i,ikind, VarID, ncNdims, dimlen,ndims,vardims(3)
    integer :: TimeDimID, CopyDimID
    
    real(r8), allocatable, dimension(:)       :: data_1d_array
    real(r8), allocatable, dimension(:,:)     :: data_2d_array
    real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
    real(r8), allocatable, dimension(:,:,:,:) :: data_4d_array
    
    character(len=128) :: filename
    
    if ( .not. module_initialized ) call static_init_model
    
    ierr = -1 ! assume things go poorly
    
    !--------------------------------------------------------------------
    ! we only have a netcdf handle here so we do not know the filename
    ! or the fortran unit number.  but construct a string with at least
    ! the netcdf handle, so in case of error we can trace back to see
    ! which netcdf file is involved.
    !--------------------------------------------------------------------
    
    write(filename,*) 'ncFileID', ncFileID
    
    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file, 
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
     'nc_write_model_vars', 'inq_dimid copy '//trim(filename))
    
    call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
     'nc_write_model_vars', 'inq_dimid time '//trim(filename))
    
    if ( output_state_vector ) then
      
      call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
       'nc_write_model_vars', 'state inq_varid '//trim(filename))
      call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
       'nc_write_model_vars', 'state put_var '//trim(filename))
      
    else

      !----------------------------------------------------------------------------
      ! We need to process the prognostic variables.
      !----------------------------------------------------------------------------
      
      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
            
          varname = trim(state_vector_vars(ikind)%varname_short)
          error_string = trim(filename)//' '//trim(varname)

          ! Ensure netCDF variable is conformable with progvar quantity.
          ! The TIME and Copy dimensions are intentionally not queried
          ! by looping over the dimensions stored in the progvar type.
          
          call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
                        'nc_write_model_vars', 'inq_varid '//trim(error_string))
          
          call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
                        'nc_write_model_vars', 'inquire '//trim(error_string))

          mystart(:)=1
          mycount(:)=1
          
          if (state_vector_vars(ikind)%nz==1) then
            ndims=2
          else
            ndims=3
          end if

          vardims(1)=state_vector_vars(ikind)%nx
          vardims(2)=state_vector_vars(ikind)%ny
          vardims(3)=state_vector_vars(ikind)%nz

          DimCheck : do i = 1,ndims
            
            write(error_string,'(a,i2,A)') 'inquire dimension ',i,trim(error_string2)
            call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
             'nc_write_model_vars', trim(error_string))
            
            if ( dimlen /= vardims(i) ) then
              write(error_string,*) trim(error_string2),' dim/dimlen ',i,dimlen,' not ',vardims(i)
              write(error_string2,*)' but it should be.'
!              call error_handler(E_ERR, 'nc_write_model_vars', trim(error_string), &
!               source, revision, revdate, text2=trim(error_string2))
            endif
            
            mycount(i) = dimlen
            
          end do DimCheck
          
          where(dimIDs == CopyDimID) mystart = copyindex
          where(dimIDs == CopyDimID) mycount = 1
          where(dimIDs == TimeDimID) mystart = timeindex
          where(dimIDs == TimeDimID) mycount = 1

          if (ndims==2) then
            allocate(data_2d_array(vardims(1),vardims(2)))
            call sv_to_field(data_2d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(error_string2))
            deallocate(data_2d_array)

          elseif (ndims==3) then

            allocate(data_3d_array(vardims(1),vardims(2),vardims(3)))
            call sv_to_field(data_3d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(error_string2))
            deallocate(data_3d_array)
          else
            write(error_string, *) 'no support for data array of dimension ', ncNdims
!            call error_handler(E_ERR,'sv_to_restart_file', string1, &
!                          source,revision,revdate)
          end if

        end if

      end do

    end if

    return

  end function nc_write_model_vars
  
  subroutine pert_model_state(state, pert_state, interf_provided)
    !------------------------------------------------------------------
    !
    ! Perturbs a model state for generating initial ensembles.
    ! The perturbed state is returned in pert_state.
    ! A model may choose to provide a NULL INTERFACE by returning
    ! .false. for the interf_provided argument. This indicates to
    ! the filter that if it needs to generate perturbed states, it
    ! may do so by adding a perturbation to each model state 
    ! variable independently. The interf_provided argument
    ! should be returned as .true. if the model wants to do its own
    ! perturbing of states.
    
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    real(r8)              :: stddev,mean

    integer               :: ikind,ilevel,i,istart,iend
    logical, save         :: random_seq_init = .false.
    
    if ( .not. module_initialized ) call static_init_model
    
    interf_provided = .true.

    ! Initialize my random number sequence (no seed is submitted here!)
    if(.not. random_seq_init) then
      call init_random_seq(random_seq)
      random_seq_init = .true.
    endif
    
    ! add some uncertainty to every state vector element
    do ikind=1,size(state_vector_vars)
      if (state_vector_vars(ikind)%is_present) then
        do ilevel=1,state_vector_vars(ikind)%nz
          istart=state_vector_vars(ikind)%state_vector_sindex(ilevel)
          iend=istart+(state_vector_vars(ikind)%nx*state_vector_vars(ikind)%ny)-1

          mean=sum(abs(state(istart:iend)))/float(iend-istart+1)
          stddev=sqrt(sum((state(istart:iend)-mean)**2))/float(iend-istart+1)

          do i=istart,iend
            pert_state(i) = random_gaussian(random_seq, state(i),model_perturbation_amplitude*stddev)
          end do
          if ((ikind==KIND_SPECIFIC_HUMIDITY) .or. &
              (ikind==KIND_CLOUD_LIQUID_WATER) .or. &
              (ikind==KIND_CLOUD_ICE)) then
            where (pert_state(istart:iend)<0.)
              pert_state(istart:iend)=0.
            end where
          end if
        end do
      end if
    enddo

    return
    
  end subroutine pert_model_state
  
  subroutine ens_mean_for_model(ens_mean)

    real(r8), dimension(:), intent(in) :: ens_mean
    return

  end subroutine ens_mean_for_model
  
  subroutine set_allowed_state_vector_vars()
    
    is_allowed_state_vector_var(:)=.FALSE.
    is_allowed_non_state_var(:)=.FALSE.

    allowed_state_vector_vars(1)=KIND_U_WIND_COMPONENT
    is_allowed_state_vector_var(KIND_U_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(1)=KIND_U_WIND_COMPONENT
    is_allowed_state_vector_var(KIND_U_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(2)=KIND_V_WIND_COMPONENT
    is_allowed_state_vector_var(KIND_V_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(3)=KIND_VERTICAL_VELOCITY
    is_allowed_state_vector_var(KIND_VERTICAL_VELOCITY)=.TRUE.
    allowed_state_vector_vars(4)=KIND_TEMPERATURE
    is_allowed_state_vector_var(KIND_TEMPERATURE)=.TRUE.
    allowed_state_vector_vars(5)=KIND_PRESSURE
    is_allowed_state_vector_var(KIND_PRESSURE)=.TRUE.
    allowed_state_vector_vars(6)=KIND_SPECIFIC_HUMIDITY
    is_allowed_state_vector_var(KIND_SPECIFIC_HUMIDITY)=.TRUE.
    allowed_state_vector_vars(7)=KIND_CLOUD_LIQUID_WATER
    is_allowed_state_vector_var(KIND_CLOUD_LIQUID_WATER)=.TRUE.
    allowed_state_vector_vars(8)=KIND_CLOUD_ICE
    is_allowed_state_vector_var(KIND_CLOUD_ICE)=.TRUE.

    allowed_non_state_vars(1)=KIND_SURFACE_ELEVATION
    is_allowed_non_state_var(KIND_SURFACE_ELEVATION)=.TRUE.
    allowed_non_state_vars(2)=KIND_SURFACE_GEOPOTENTIAL
    is_allowed_non_state_var(KIND_SURFACE_GEOPOTENTIAL)=.TRUE.

    return

  end subroutine set_allowed_state_vector_vars

  function ll_to_xyz_vector(lon,lat) RESULT (xyz)
    
    ! Passed variables
    
    real(r8),allocatable :: xyz(:,:)      ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat(:),lon(:) ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius

    integer              :: n

    ! define output vector size to be the same as the input vector size
    ! second dimension (3) is x,y,z

    n=SIZE(lat,1)
    ALLOCATE(xyz(1:n,1:3))

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.

    ! caclulate the x,y,z-coordinates

    xyz(1:n,1)=radius*sin(lat(1:n)*deg2rad)*cos(lon(1:n)*deg2rad)
    xyz(1:n,2)=radius*sin(lat(1:n)*deg2rad)*sin(lon(1:n)*deg2rad)
    xyz(1:n,3)=radius*cos(lat(1:n)*deg2rad)
    
    return
  end function ll_to_xyz_vector

  function ll_to_xyz_single(lon,lat) result (xyz)
    
    ! Passed variables
    
    real(r8)             :: xyz(1:3) ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat,lon  ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius

    integer              :: i,j,n

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.

    ! caclulate the x,y,z-coordinates

    xyz(1)=radius*sin(lat*deg2rad)*cos(lon*deg2rad)
    xyz(2)=radius*sin(lat*deg2rad)*sin(lon*deg2rad)
    xyz(3)=radius*cos(lat*deg2rad)
    
    return
  end function ll_to_xyz_single

  subroutine get_enclosing_grid_box(p,g,n,nx,ny,b,bw)
    
    integer,intent(in)   :: n,nx,ny
    real(r8),intent(in)  :: p(1:3),g(1:n,1:3)
    integer,intent(out)  :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)             :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:3),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2)
    integer              :: i,j,minidx(2),boxidx(2),idx1,idx2,xb,yb

    work(2:nx+1,2:ny+1,1:3)=RESHAPE( g, (/ nx,ny,3 /))

    do i=2,nx+1
      work(i,1,1:3)=work(i,2,1:3)-(work(i,3,1:3)-work(i,2,1:3))
      work(i,ny+2,1:3)=work(i,ny+1,1:3)-(work(i,ny,1:3)-work(i,ny+1,1:3))
    end do

    do j=2,ny+1
      work(1,j,1:3)=work(2,j,1:3)-(work(3,j,1:3)-work(2,j,1:3))
      work(nx+2,j,1:3)=work(nx+1,j,1:3)-(work(nx,j,1:3)-work(nx+1,j,1:3))
    end do

    work(1,1,1:3)=work(2,2,1:3)-0.5*(sqrt(2.)*(work(2,2,1:3)-work(1,2,1:3))+sqrt(2.)*(work(2,2,1:3)-work(2,1,1:3)))
    work(1,ny+2,1:3)=work(2,ny+1,1:3)-0.5*(sqrt(2.)*(work(2,ny+1,1:3)-work(1,ny+1,1:3))+sqrt(2.)*(work(2,ny+1,1:3)-work(2,ny+2,1:3)))
    work(nx+2,1,1:3)=work(nx+1,2,1:3)-0.5*(sqrt(2.)*(work(nx+1,2,1:3)-work(nx+2,2,1:3))+sqrt(2.)*(work(nx+1,2,1:3)-work(nx+1,1,1:3)))
    work(nx+2,ny+2,1:3)=work(nx+1,ny+1,1:3)-0.5*(sqrt(2.)*(work(nx+1,ny+1,1:3)-work(nx+2,ny+1,1:3))+sqrt(2.)*(work(nx+1,ny+1,1:3)-work(nx+1,ny+2,1:3)))

    do i=1,nx+2
      do j=1,ny+2
        dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
      end do
    end do

    minidx(:)=minloc(dist)

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    end if


    do i=0,1
      do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))
      end do
    end do 

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5))
    yb=minidx(2)+(2*(boxidx(2)-0.5))

!    print*,minidx
!    print*,xb,yb

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
        do j=1,2
          b(i,j)=((minidx(2)+(j-1)*(boxidx(2)-0.5)*2)*ny)+(minidx(1)+(i-1)*(2*(boxidx(1)-0.5)))
        end do
      end do

      do i=1,2
        do j=1,2
          boxdist(i,j)=dist(mod(b(i,j),ny),b(i,j)/ny)
        end do
      end do

      bw(:,:)=1./boxdist(:,:)
!      bw(:,:)=(((1.-boxdist(:,:))/(1.1*maxval(boxdist)))**2)/((boxdist(:,:)/(1.1*maxval(boxdist)))**2)
      bw=bw/sum(bw)
      b(:,:)=b(:,:)-1
    end if

    return

  end subroutine get_enclosing_grid_box

  subroutine get_enclosing_grid_box_lonlat(lon,lat,p,n,nx,ny,b,bw)
    
    integer,intent(in)   :: n,nx,ny
    real(r8),intent(in)  :: p(1:2),lon(1:n),lat(1:n)
    integer,intent(out)  :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)             :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:2),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2),pw(2)

    integer              :: i,j,minidx(2),boxidx(2),idx1,idx2,xb,yb,bx(2,2),by(2,2)

    work(2:nx+1,2:ny+1,1)=reshape(lon,(/ nx,ny /))*deg2rad
    work(2:nx+1,2:ny+1,2)=reshape(lat,(/ nx,ny /))*deg2rad
    pw=p*deg2rad

    do i=2,nx+1
      work(i,1,1:2)=work(i,2,1:2)-(work(i,3,1:2)-work(i,2,1:2))
      work(i,ny+2,1:2)=work(i,ny+1,1:2)-(work(i,ny,1:2)-work(i,ny+1,1:2))
    end do

    do j=2,ny+1
      work(1,j,1:2)=work(2,j,1:2)-(work(3,j,1:2)-work(2,j,1:2))
      work(nx+2,j,1:2)=work(nx+1,j,1:2)-(work(nx,j,1:2)-work(nx+1,j,1:2))
    end do

    work(1,1,1:2)=work(2,2,1:2)-0.5*(sqrt(2.)*(work(2,2,1:2)-work(1,2,1:2))+sqrt(2.)*(work(2,2,1:2)-work(2,1,1:2)))
    work(1,ny+2,1:2)=work(2,ny+1,1:2)-0.5*(sqrt(2.)*(work(2,ny+1,1:2)-work(1,ny+1,1:2))+sqrt(2.)*(work(2,ny+1,1:2)-work(2,ny+2,1:2)))
    work(nx+2,1,1:2)=work(nx+1,2,1:2)-0.5*(sqrt(2.)*(work(nx+1,2,1:2)-work(nx+2,2,1:2))+sqrt(2.)*(work(nx+1,2,1:2)-work(nx+1,1,1:2)))
    work(nx+2,ny+2,1:2)=work(nx+1,ny+1,1:2)-0.5*(sqrt(2.)*(work(nx+1,ny+1,1:2)-work(nx+2,ny+1,1:2))+sqrt(2.)*(work(nx+1,ny+1,1:2)-work(nx+1,ny+2,1:2)))

    do i=1,nx+2
      do j=1,ny+2
!        dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
        dist(i,j)=6173.*acos(cos(work(i,j,2)-pw(2))-cos(work(i,j,2))*cos(pw(2))*(1-cos(work(i,j,1)-pw(1))))
      end do
    end do

    minidx(:)=minloc(dist)

!    print*,minidx

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    end if

    open(21,file='/daten02/jkeller/testbox.bin',form='unformatted')
    write(21) nx
    write(21) ny

!    print*,mod(b(i,j),nx),(b(i,j)/nx)+1

    do i=0,1
      do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))/4.
!        write(*,'(4(I5))') minidx(1)+i-1,minidx(1)+i,minidx(2)+j-1,minidx(2)+j
        write(21) (minidx(2)+j-1),minidx(1)+i-1,&
                  (minidx(2)+j-1),minidx(1)+i,&
                  (minidx(2)+j),minidx(1)+i-1,&
                  (minidx(2)+j),minidx(1)+i
        write(21) boxdist(i+1,j+1)
      end do
    end do 

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5))
    yb=minidx(2)+(2*(boxidx(2)-0.5))

!    print*,minidx
!    print*,xb,yb

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
        do j=1,2
          bx(i,j)=minidx(1)+(i-1)*(2*(boxidx(1)-0.5))
          by(i,j)=minidx(2)+(j-1)*(2*(boxidx(2)-0.5))
!          b(i,j)=(((minidx(2)+(j-1)*(boxidx(2)-0.5)*2)-1)*nx)+(minidx(1)+(i-1)*(2*(boxidx(1)-0.5)))
        end do
      end do

      do i=1,2
        do j=1,2
          boxdist(i,j)=dist(bx(i,j),by(i,j))
        end do
      end do

      bw(:,:)=1./boxdist(:,:)
!      bw(:,:)=(((1.-boxdist(:,:))/(1.1*maxval(boxdist)))**2)/((boxdist(:,:)/(1.1*maxval(boxdist)))**2)
      bw=bw/sum(bw)
      bx=bx-1
      by=by-1
      b(:,:)=(by-1)*nx+bx
    end if

!    write(*,'(4(F6.3,1X))') lon(b(1,1)),lat(b(1,1)),lon(b(2,1)),lat(b(2,1))
!    write(*,'(4(F6.3,1X))') lon(b(1,2)),lat(b(1,2)),lon(b(2,2)),lat(b(2,2))

    write(21) b
    write(21) minidx-1
    write(21) lon
    write(21) lat
    close(21)
    return

  end subroutine get_enclosing_grid_box_lonlat


!  subroutine bilinear_interpolation(b,p,data,v,istatus,lon,lat)

  subroutine bilinear_interpolation(bv,bw,blo,bla,p,v,istatus)

    ! Passed variables
    
    real(r8),intent(in)  :: bv(2,2),bw(2,2),blo(2,2),bla(2,2)
    real(r8),intent(in)  :: p(3)
    real(r8),intent(out) :: v
    integer,intent(out)  :: istatus
    
    ! Local storage
    
    real(r8)             :: x1,x2,lo1,lo2,la1,la2
    real(r8)             :: m1,m2,n1,n2,xc,yc,d1,d2,d
    
!    write(*,'(3(F8.5,1X))') bv(1,1),blo(1,1),bla(1,1)
!    write(*,'(3(F8.5,1X))') bv(2,1),blo(2,1),bla(2,1)

    call linear_interpolation(p(1),p(2),bv(1,1),blo(1,1),bla(1,1),&
                                        bv(2,1),blo(2,1),bla(2,1),&
                                        x1,lo1,la1)

!    write(*,'(3(F8.5,1X))') x1,lo1,la1

    call linear_interpolation(p(1),p(2),bv(1,2),blo(1,2),bla(1,2),&
                                        bv(2,2),blo(2,2),bla(2,2),&
                                        x2,lo2,la2)

!    write(*,'(3(F8.5,1X))') x2,lo2,la2

    d1=sqrt((lo1-p(1))**2+(la1-p(2))**2)
    d2=sqrt((lo2-p(1))**2+(la2-p(2))**2)
    d=sqrt((lo1-lo2)**2+(la1-la2)**2)

    v=(1.-(d1/d))*x1+(1.-(d2/d))*x2

!    print*,1.-d1/d,1.-d2/d
!    print*,v
    
    return

  end subroutine bilinear_interpolation

  subroutine linear_interpolation(lop,lap,x1,lo1,la1,x2,lo2,la2,x,lo,la)

    real(r8),intent(in)  :: lo1,lo2,la1,la2,x1,x2,lop,lap
    real(r8),intent(out) :: lo,la,x

    real(r8)             :: m1,m2,n1,n2,d1,d2,d

    m1=(la2-la1)/(lo2-lo1)
    if (m1 .ne. 0.) then
      n1=la1-lo1*m1
      m2=-1./m1
      n2=lap-lop*m2
      lo=(n2-n1)/(m1-m2)
      la=lo*m1+n1

      d1=sqrt((lo1-lo)**2+(la1-la)**2)
      d2=sqrt((lo2-lo)**2+(la2-la)**2)
      d=sqrt((lo1-lo2)**2+(la1-la2)**2)
    else
      la=la1
      lo=lop

      d1=sqrt((lo1-lo)**2+(la1-la)**2)
      d2=sqrt((lo2-lo)**2+(la2-la)**2)
      d=sqrt((lo1-lo2)**2+(la1-la2)**2)
    end if

    x=(1.-(d1/d))*x1+(1.-(d2/d))*x2

    return

  end subroutine linear_interpolation

!  SUBROUTINE data_to_state_vector(bdata,nvar,bpos,blen,allowed,v,sv,x)
!
!    ! data_to_state_vector calls subroutines to extract the data fields from the binary data
!    ! and put it into the state vector
!
!    INTEGER(kind=1),INTENT(in)       :: bdata(:)
!    INTEGER,INTENT(in)               :: nvar
!    INTEGER,INTENT(in)               :: bpos(1:nvar,1:4)
!    INTEGER,INTENT(in)               :: blen(1:nvar,1:4)
!    LOGICAL,INTENT(in)               :: allowed(:)
!    TYPE(cosmo_meta),INTENT(in)      :: v(1:nvar)
!    TYPE(dart_variable_info)         :: sv(:)
!    REAL(r8),ALLOCATABLE,INTENT(out) :: x(:)
!
!    INTEGER                          :: index,ivar,ikind,ilevel,n
!    REAL(r8),ALLOCATABLE             :: data(:,:)
!
!    n=0
!    DO ivar=1,nvar
!      ikind=v(ivar)%dart_kind
!      IF (allowed(ikind)) THEN
!        n=n+v(ivar)%dims(1)*v(ivar)%dims(2)
!      END IF
!    END DO
!
!    DO ivar=1,nvar
!      ikind=v(ivar)%dart_kind
!      IF (allowed(ikind)) THEN
!!        CALL read_data(bdata,bpos(ivar,4))
!!        CALL field_to_state_vector(bdata,nvar,bpos,blen,v,x)
!      END IF
!    END DO
!        
!!        ALLOCATE(state_vector_vars(ikind)%state_vector_sindex(1:state_vector_vars(ikind)%nz))
!
!
!  END SUBROUTINE data_to_state_vector
!
  subroutine get_vertical_boundaries(hb,hw,otype,vcs,p,b,w,istatus)

    real(r8),intent(in)  :: hw(2,2),p,vcs
    integer,intent(in)   :: hb(2,2),otype
    integer,intent(out)  :: b(2),istatus
    real(r8),intent(out) :: w(2)

    integer              :: k,nlevel,x1,x2,x3,x4,y1,y2,y3,y4
    real(r8)             :: u,l
    real(r8),allocatable :: klevel(:),hlevel(:),plevel(:)

    b(:)=-1

    if (vcs==-2. .or. vcs==4. ) then
      istatus=19
      return
    end if

    ! surface not implemented, return out of vertical domain
    if (vcs==-1.) then
      istatus=16
      return
    end if

    x1=mod(hb(1,1),size(non_state_data%pfl,1))
    x2=mod(hb(2,1),size(non_state_data%pfl,1))
    x3=mod(hb(1,2),size(non_state_data%pfl,1))
    x4=mod(hb(2,2),size(non_state_data%pfl,1))
    y1=hb(1,1)/size(non_state_data%pfl,1)
    y2=hb(2,1)/size(non_state_data%pfl,1)
    y3=hb(1,2)/size(non_state_data%pfl,1)
    y4=hb(2,2)/size(non_state_data%pfl,1)
    
    if (otype .ne. KIND_VERTICAL_VELOCITY) then
      nlevel=non_state_data%nfl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hfl(x1,y1,:)+&
             hw(2,1)*non_state_data%hfl(x2,y2,:)+&
             hw(1,2)*non_state_data%hfl(x3,y3,:)+&
             hw(2,2)*non_state_data%hfl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%pfl(x1,y1,:)+&
             hw(2,1)*non_state_data%pfl(x2,y2,:)+&
             hw(1,2)*non_state_data%pfl(x3,y3,:)+&
             hw(2,2)*non_state_data%pfl(x4,y4,:)
    else
      nlevel=non_state_data%nhl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hhl(x1,y1,:)+&
             hw(2,1)*non_state_data%hhl(x2,y2,:)+&
             hw(1,2)*non_state_data%hhl(x3,y3,:)+&
             hw(2,2)*non_state_data%hhl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%phl(x1,y1,:)+&
             hw(2,1)*non_state_data%phl(x2,y2,:)+&
             hw(1,2)*non_state_data%phl(x3,y3,:)+&
             hw(2,2)*non_state_data%phl(x4,y4,:)
    end if

    do k=1,nlevel-1
      if (vcs==1.) then
        u=klevel(k+1)
        l=klevel(k)
      end if
      if (vcs==2.) then
        u=plevel(k+1)
        l=plevel(k)
      end if
      if (vcs==3.) then
        u=hlevel(k+1)
        l=hlevel(k)
      end if

      if (u>=p .and. l<=p) then
        b(1)=k
        b(2)=k+1
        w(1)=1.-(p-l)/(u-l)
        w(2)=1.-(u-p)/(u-l)
        return
      end if

    end do

    istatus=16
    return
  end subroutine get_vertical_boundaries

  subroutine sv_to_field_2d(f,x,v)

    real(r8),intent(out)                :: f(:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie

    is=v%state_vector_sindex(1)
    ie=is+v%nx*v%ny-1
    f(:,:) = reshape(x(is:ie),(/ v%nx,v%ny /))

    return

  end subroutine sv_to_field_2d

  subroutine sv_to_field_3d(f,x,v)

    real(r8),intent(out)                :: f(:,:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie,iz

    do iz=1,v%nz
      is=v%state_vector_sindex(iz)
      ie=is+v%nx*v%ny-1
      f(:,:,iz) = reshape(x(is:ie),(/ v%nx,v%ny /))
    end do
    
    return

  end subroutine sv_to_field_3d

!  subroutine grib_to_sv(filename, state_vector, model_time)
!    !------------------------------------------------------------------
!    ! Reads the current time and state variables from a cosmo grib
!    ! file and packs them into a dart state vector.
!    
!    character(len=*), intent(in)    :: filename 
!    real(r8),         intent(inout) :: state_vector(:)
!    type(time_type),  intent(out)   :: model_time
!    
!    ! temp space to hold data while we are reading it
!    integer  :: i, j, k, l, ni, nj, nk, nl, ivar, indx
!    real(r8), allocatable, dimension(:)         :: data_1d_array
!    real(r8), allocatable, dimension(:,:)       :: data_2d_array
!    real(r8), allocatable, dimension(:,:,:)     :: data_3d_array
!    real(r8), allocatable, dimension(:,:,:,:)   :: data_4d_array
!    
!    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
!    character(len=NF90_MAX_NAME) :: varname 
!    integer :: sidx,eidx,cosmo_var
!    
!    if ( .not. module_initialized ) call static_init_model
!    
!    state_vector(:) = MISSING_R8
!    
!    ! Check that the input file exists ... 
!    
!    model_time = get_state_time()
!    
!    STOP
!    if (do_output()) &
!     call print_time(model_time,'time in restart file '//trim(filename))
!    if (do_output()) &
!     call print_date(model_time,'date in restart file '//trim(filename))
!    
!    sidx=1
!
!    do ivar=1,n_state_vector_vars
!      
!      if (state_vector_vars(ivar)%is_present) then
!        varname = trim(state_vector_vars(ivar)%varname_short)
!        error_string = trim(filename)//' '//trim(varname)
!
!        do ilevel=1,state_vector_vars(ivar)%nz
!
!          cosmo_var=state_vector_vars(ivar)%cosmo_state_index(ilevel)
!          len=(state_vector_vars(ivar)%nx*state_vector_vars(ivar)%ny)-1
!          eidx=sidx+len-1
!          data=get_data(bdata,bpos(cosmo_var,1:4),blen(cosmo_var,1:4))
!          x(sidx:eidx)=reshape(data,(/ len /))
!          
!        end do
!      end if
!    end do
! 
!    call nc_check(nf90_close(ncid), &
!     'grib_to_sv','close '//trim(filename))
!    
!  end subroutine grib_to_sv
  
  function get_state_time result (time)
    type(time_type) :: time
    
    if ( .not. module_initialized ) call static_init_model
    time=cosmo_fc_time

    return

  end function get_state_time

  function get_state_vector result (sv)

    real(r8)             :: sv(1:model_size)

    integer              :: ivar,ikind,nx,ny,sidx,eidx
    real(r8),allocatable :: data(:,:)

    if ( .not. module_initialized ) call static_init_model

    call set_allowed_state_vector_vars()

    do ivar=1,nvars
      ikind=cosmo_vars(ivar)%dart_kind
      if (is_allowed_state_vector_var(ikind)) then
        nx=state_vector_vars(ikind)%nx
        ny=state_vector_vars(ikind)%ny
        allocate(data(1:nx,1:ny))
        data=get_data_from_binary(cosmo_filename,grib_header(ivar),nx,ny)
        sidx=cosmo_vars(ivar)%dart_sindex
        eidx=cosmo_vars(ivar)%dart_eindex
        state_vector(sidx:eidx)=reshape(data,(/ (nx*ny) /))
        deallocate(data)
      end if
    end do

    sv(:)=state_vector(:)

    return

  end function get_state_vector

  subroutine write_grib_file(sv,nfile)

    real(r8),intent(in)           :: sv(:)
    character(len=128),intent(in) :: nfile

    integer                       :: irec,ivar,istat,ipos,nrec
    integer                       :: len,hlen,ix,iy,nx,ny,idx,naddbyte
    integer                       :: dval,ibsf,idsf
    integer,allocatable           :: griblen(:)
    real(r8)                      :: bsf,dsf
    integer(kind=1)               :: bin4(4),gribword(4),bin42(4)
    integer(kind=1),allocatable   :: bytearr(:)
    real(r8),allocatable          :: data(:,:)
    real(r8)                      :: ref_value


    OPEN(10,FILE=TRIM(cosmo_filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)

    if ( .not. module_initialized ) call static_init_model
!    sv(:)=state_vector(:)

!    generate byte_array to write
    len=0
    nvars=2
    allocate(griblen(1:nvars))
    DO ivar=1,nvars
      len=len+size(grib_header(ivar)%pds)+size(grib_header(ivar)%gds)+grib_header(ivar)%data_length+8+4
      griblen(ivar)=size(grib_header(ivar)%pds)+size(grib_header(ivar)%gds)+grib_header(ivar)%data_length+8+4
    END DO
    
!    if (MOD(len,4) .NE. 0) len=len+4-MOD(len,4)

    ALLOCATE(bytearr(1:len))
    bytearr(:)=0

    ipos=1

    DO ivar=1,nvars

!temp(ivar)%start_record=irec
!          temp(ivar)%start_offset=ioff

      nx=cosmo_vars(ivar)%dims(1)
      ny=cosmo_vars(ivar)%dims(2)
      allocate(data(1:nx,1:ny))

      idx=cosmo_vars(ivar)%dart_sindex
      
      data(:,:)=reshape(sv(idx:idx+nx*ny),(/ nx,ny /))

      ! write word GRIB
      gribword(1)=ICHAR('G')
      gribword(2)=ICHAR('R')
      gribword(3)=ICHAR('I')
      gribword(4)=ICHAR('B')
      bytearr(ipos:ipos+3)=gribword
      ipos=ipos+4

      call word_to_byte(griblen(ivar),bytearr(ipos:ipos+2),3)
      bytearr(ipos+3)=1
      print*,bytearr(ipos:ipos+3)
      ipos=ipos+4

      ! write PDS
      hlen=size(grib_header(ivar)%pds)
      bytearr(ipos:ipos+hlen-1)=grib_header(ivar)%pds
      ipos=ipos+hlen

      ! write GDS
      hlen=size(grib_header(ivar)%gds)
      bytearr(ipos:ipos+hlen-1)=grib_header(ivar)%gds
      ipos=ipos+hlen

      iy=ipos

      ! write BDS-header
      hlen=size(grib_header(ivar)%bds)
      bytearr(ipos:ipos+hlen-1)=grib_header(ivar)%bds 
      print*,grib_header(ivar)%bds

      ref_value=minval(data)
      bin4(1:4)=from_float1(ref_value)
      bytearr(ipos+6:ipos+9)=bin4
      CALL byte_to_word_signed(bytearr(ipos+4:ipos+5),ibsf,2)
      bsf=FLOAT(ibsf)

      CALL byte_to_word_signed(grib_header(ivar)%pds(27:28),idsf,2)
      dsf=FLOAT(idsf)

      ipos=ipos+hlen

      DO ix=iy,iy+30
        print*,ix,bytearr(ix)
      END DO

      DO iy=1,ny
        DO ix=1,nx
!          print*,data(ix,iy)
          dval=int((data(ix,iy)-ref_value)*((10.**dsf)/(2.**bsf)))
!          print*,word_to_byte_data(dval)
          bytearr(ipos:ipos+1)=word_to_byte_data(dval)
          ipos=ipos+2
        END DO
      END DO

      deallocate(data)

      naddbyte=IAND(grib_header(ivar)%bds(4),15)/8
      ipos=ipos+naddbyte

      ! write word GRIB
      gribword(1)=ICHAR('7')
      gribword(2)=ICHAR('7')
      gribword(3)=ICHAR('7')
      gribword(4)=ICHAR('7')
      bytearr(ipos:ipos+3)=gribword
      ipos=ipos+4

    END DO

    OPEN(11,FILE=TRIM(nfile),FORM='UNFORMATTED')
    print*,len
    WRITE(11) bytearr(1:len)

    ipos=1
    nrec=(len/4)
    DO irec=1,nrec
!      if (ipos .GT. 300) print*,ipos,ipos+3
!      READ(10,rec=irec,iostat=istat) bin4
!      if (ipos .GT. 300) print*,bin4
      bin4=bytearr(ipos:ipos+3)
!      WRITE(11,rec=irec,iostat=istat) bin4
!      if (ipos .GT. 300) print*,bin4
!      if (ipos .GT. 300) print*,''
      ipos=ipos+4
    END DO

    

!    if ((len-nrec*4)>0) WRITE(11,rec=nrec+1,iostat=istat) bytearr(ipos:ipos+(len-nrec*4)-1)

    CLOSE(10)
    CLOSE(11)

!    istat=0
!    irec=1
!    print*,'vergleich'
!    DO WHILE(istat==0)
!      OPEN(10,FILE=TRIM(cosmo_filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
!      OPEN(11,FILE=TRIM(nfile),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
!      READ(10,rec=irec,iostat=istat) bin4
!      READ(11,rec=irec,iostat=istat) bin42
!      print*,irec,bin4-bin42
!      print*,irec,bin4
!      print*,irec,bin42
!      print*,''
!      irec=irec+1
!    ENDDO
!    CLOSE(10)
!    CLOSE(11)

    return

  end subroutine write_grib_file
  
end module model_mod
