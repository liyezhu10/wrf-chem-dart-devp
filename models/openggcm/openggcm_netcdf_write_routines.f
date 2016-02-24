! This is F90 code written with F77 syntax because (for expediency) it will
! simply be appended to the mhd-iono.for file.
! This means all line continuation characters are in column 6
! Appently compiles with long lines, so can go beyond column 72

      module netcdf_mod
      
      use netcdf
      implicit none
      private
      
      public :: wr_netcdf_model_time,
     &          wr_netcdf_ctim_grid,
     &          wr_netcdf_interface_grid,
     &          wr_netcdf_2D,
     &          wr_netcdf_3D
      
      contains 
      
!=======================================================================
      
      subroutine wr_netcdf_model_time(ncid, model_time)
      
      integer, intent(in) :: ncid
      real*8,  intent(in) :: model_time
      
      integer :: VarID
      integer :: io1, io2
      
!     io1 = nf90_redef(ncid) 
!     call nc_check(io1, 'wr_netcdf_model_time', 'redef')

      io1 = nf90_def_var(ncid, 'time', nf90_double, VarID)
      io2 = nf90_put_att(ncid, VarID, 'units', 'seconds since 1966-01-01')

      call nc_check(io1, 'wr_netcdf_model_time', 'def_var time')
      call nc_check(io2, 'wr_netcdf_model_time', 'put_att:time:units')
      
      io1 = nf90_enddef(ncid)
      call nc_check(io1, 'wr_netcdf_model_time', 'enddef')
      
      io1 = nf90_put_var(ncid, VarID, model_time)
      call nc_check(io1, 'wr_netcdf_model_time', 'put_var: model_time')
      
      end subroutine wr_netcdf_model_time
      
!=======================================================================
      
!> ctim grid will have variables:
!> cg_lon
!> cg_lat
!> cg_height
      
      subroutine wr_netcdf_ctim_grid(ncid, nlon, nlat, nheight)
      
      integer, intent(in) :: ncid
      integer, intent(in) :: nlon
      integer, intent(in) :: nlat
      integer, intent(in) :: nheight
      
      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------
      
      real*4, allocatable :: latitude(:)
      real*4, allocatable :: longitude(:)
      real*4, allocatable :: height(:)
      real*4  :: dlat, dlon, dheight
      
      integer :: LonDimID, LatDimID, HeightDimID
      integer :: LonVarID, LatVarID, HeightVarID
      integer :: io, io1, io2, io3
      integer :: ilat, ilon, iheight
      
      !-----------------------------------------------------------------------
      ! calculate the required metadata
      
      allocate(longitude(nlon), latitude(nlat), height(nheight))
      
      dlon = 360.0/real(nlon-1)
      do ilon = 1,nlon
         longitude(ilon) = dlon*real(ilon-1)
      enddo
      
      dlat = 180.0/real(nlat-1)
      do ilat = 1,nlat  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
         latitude(ilat) = dlat*real(ilat-1)
      enddo
      
      dheight = 500.0/real(nheight)   ! TJH FIXME ... bogus height array
      do iheight = 1,nheight
         height(iheight) = dheight*real(iheight-1)
      enddo
      
      !-----------------------------------------------------------------------

      io = nf90_redef(ncid) 
      call nc_check(io, 'wr_netcdf_ctim_grid', 'redef')
      
      io1 = nf90_def_dim(ncid, 'cg_lon', nlon, LonDimID)
      io2 = nf90_def_var(ncid, 'cg_lon', nf90_real, (/ LonDimID /), LonVarID)
      io3 = nf90_put_att(ncid, LonVarID, 'short_name', 'geographic longitude' )
      io  = nf90_put_att(ncid, LonVarID, 'units', 'degrees')
      
      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_lon')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_lon')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_lon short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_lon units')
      
      io1 = nf90_def_dim(ncid, 'cg_lat', nlat, LatDimID)
      io2 = nf90_def_var(ncid, 'cg_lat', nf90_real, (/ LatDimID /), LatVarID)
      io3 = nf90_put_att(ncid, LatVarID, 'short_name', 'geographic latitude')
      io  = nf90_put_att(ncid, LatVarID, 'units', 'degrees')

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_lat')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_lat')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_lat short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_lat units')

      io1 = nf90_def_dim(ncid, 'cg_height', nheight, HeightDimID)
      io2 = nf90_def_var(ncid, 'cg_height', nf90_real, (/ HeightDimID /), HeightVarID)
      io3 = nf90_put_att(ncid, HeightVarID, 'short_name', 'geographic height')
      io  = nf90_put_att(ncid, HeightVarID, 'units', 'kilometers')

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_height')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_height')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_height short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_height units')

      ! leave define mode so we can fill
      
      io = nf90_enddef(ncid)
      call nc_check(io, 'wr_netcdf_ctim_grid', 'enddef')
      
      io1 = nf90_put_var(ncid,    LonVarID, longitude)
      io2 = nf90_put_var(ncid,    LatVarID,  latitude)
      io3 = nf90_put_var(ncid, HeightVarID,    height)

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'put_var: cg_lon')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'put_var: cg_lat')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_var: cg_height')
      
      end subroutine wr_netcdf_ctim_grid
      
!=======================================================================
      
      subroutine wr_netcdf_interface_grid(ncid, nlon, nlat, nheight)
      
      integer, intent(in) :: ncid
      integer, intent(in) :: nlon
      integer, intent(in) :: nlat
      integer, intent(in) :: nheight
      
      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------
      
      real*4, allocatable :: latitude(:)
      real*4, allocatable :: longitude(:)
      real*4, allocatable :: height(:)
      real*4  :: dlat, dlon, dheight
      
      integer :: LonDimID, LatDimID, HeightDimID
      integer :: LonVarID, LatVarID, HeightVarID
      integer :: io, io1, io2, io3
      integer :: ilat, ilon, iheight
      
      !=======================================================================

      write(*,*)'ERROR: wr_netcdf_interface_grid not written, not tested'
      write(*,*)'ERROR: wr_netcdf_interface_grid not written, not tested'
      write(*,*)'ERROR: wr_netcdf_interface_grid not written, not tested'
      write(*,*)'ERROR: wr_netcdf_interface_grid not written, not tested'
      stop

      !-----------------------------------------------------------------------
      ! calculate the required metadata
      
      allocate(longitude(nlon), latitude(nlat), height(nheight))
      
      dlon = 360.0/real(nlon-1)
      do ilon = 1,nlon
         longitude(ilon) = dlon*real(ilon-1)
      enddo
      
      dlat = 180.0/real(nlat-1)
      do ilat = 1,nlat  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
         latitude(ilat) = dlat*real(ilat-1)
      enddo
      
      dheight = 500.0/real(nheight)   ! TJH FIXME ... bogus height array
      do iheight = 1,nheight
         height(iheight) = dheight*real(iheight-1)
      enddo
      
      !-----------------------------------------------------------------------

      io = nf90_redef(ncid) 
      call nc_check(io, 'wr_netcdf_interface_grid', 'redef')
      
      io1 = nf90_def_dim(ncid, 'ig_lon', nlon, LonDimID)
      io2 = nf90_def_var(ncid, 'ig_lon', nf90_real, (/ LonDimID /), LonVarID)
      io3 = nf90_put_att(ncid, LonVarID, 'short_name', 'geographic longitude' )
      
      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_dim ig_lon')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'def_var ig_lon')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_lon short_name')
      
      io1 = nf90_def_dim(ncid, 'ig_lat', nlat, LatDimID)
      io2 = nf90_def_var(ncid, 'ig_lat', nf90_real, (/ LatDimID /), LatVarID)
      io3 = nf90_put_att(ncid, LatVarID, 'short_name', 'geographic latitude')

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_dim ig_lat')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'def_var ig_lat')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_lat short_name')

      io1 = nf90_def_dim(ncid, 'ig_height', nheight, HeightDimID)
      io2 = nf90_def_var(ncid, 'ig_height', nf90_real, (/ HeightDimID /), HeightVarID)
      io3 = nf90_put_att(ncid, LatVarID, 'short_name', 'geographic height')

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_dim ig_height')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'def_var ig_height')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_height short_name')

      ! leave define mode so we can fill
      
      io = nf90_enddef(ncid)
      call nc_check(io, 'wr_netcdf_interface_grid', 'enddef')
      
      io1 = nf90_put_var(ncid,    LonVarID, longitude)
      io2 = nf90_put_var(ncid,    LatVarID,  latitude)
      io3 = nf90_put_var(ncid, HeightVarID,    height)

      call nc_check(io1, 'wr_netcdf_interface_grid', 'put_var: ig_lon')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'put_var: ig_lat')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_var: ig_height')
      
      end subroutine wr_netcdf_interface_grid
      
!=======================================================================
      
      subroutine wr_netcdf_3D(ncid, gridtype, tensor, tensorname, tensorunits, tensorshort)
      
      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      real*4,           intent(in) :: tensor(:,:,:)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort
      
      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------
      
      integer :: LonDimID, LatDimID, HeightDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon, nlat, nheight
      
      ! call nc_check(nf90_create(filename, NF90_CLOBBER, ncid),'wr_netcdf_3D')
      !-----------------------------------------------------------------------
      
      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      else
         write(*,*)'ERROR :: stop'
         stop
      endif
      
      ! Make sure the netcdf dimensions match the incoming variable
      
      io1 = nf90_inquire_dimension(ncid,    LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid,    LatDimID, len=nlat)
      io3 = nf90_inquire_dimension(ncid, HeightDimID, len=nheight)
      
      call nc_check(io1, 'wr_netcdf_3D', 'inquire lon    dimension length')
      call nc_check(io2, 'wr_netcdf_3D', 'inquire lat    dimension length')
      call nc_check(io3, 'wr_netcdf_3D', 'inquire height dimension length')
      
      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension stop'
         stop
      endif 
      
      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension stop'
         stop
      endif 
      
      if (nheight /= size(tensor,3)) then
         write(*,*)'ERROR :: failed height dimension stop'
         stop
      endif 
      
      ! If everything is consistent, go ahead and define it
      
      io = nf90_redef(ncid) 
      call nc_check(io, 'wr_netcdf_3D', 'redef '//trim(tensorname))
      
      ! Now define the variable/tensor/hyperslab
      
      io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ LonDimID, LatDimID, HeightDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_3D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_3D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_3D', 'put_att units'//trim(tensorname))
      
      ! leave define mode so we can fill
      
      io = nf90_enddef(ncid)
      call nc_check(io, 'wr_netcdf_3D', 'enddef')
      
      io = nf90_put_var(ncid,     VarID, tensor)
      call nc_check(io, 'wr_netcdf_3D', 'put_var: tensor')
      
      end subroutine wr_netcdf_3D
      
!=======================================================================
      
      subroutine wr_netcdf_2D(ncid, gridtype, tensor, tensorname, tensorunits, tensorshort)
      
      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      real*4,           intent(in) :: tensor(:,:)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort
      
      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------
      
      integer :: LonDimID, LatDimID, VarID 
      integer :: io, io1, io2, io3
      integer :: nlon, nlat
      
      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID)
      else
         write(*,*)'ERROR :: stop'
         stop
      endif
      
      ! Make sure the netcdf dimensions match the incoming variable
      
      io1 = nf90_inquire_dimension(ncid, LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid, LatDimID, len=nlat)
      
      call nc_check(io1, 'wr_netcdf_2D', 'inquire lon dimension length'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_2D', 'inquire lat dimension length'//trim(tensorname))
      
      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension stop'
         stop
      endif
      
      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension stop'
         stop
      endif
      
      ! call nc_check(nf90_create(filename, NF90_CLOBBER, ncid),'wr_netcdf_2D')
      !-----------------------------------------------------------------------
      io = nf90_redef(ncid) 
      call nc_check(io, 'wr_netcdf_2D', 'redef'//trim(tensorname))
      
      io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ LonDimID, LatDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))
      
      call nc_check(io1, 'wr_netcdf_2D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_2D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_2D', 'put_att units'//trim(tensorname))
      
      ! leave define mode so we can fill
      
      io = nf90_enddef(ncid)
      call nc_check(io, 'wr_netcdf_2D', 'enddef'//trim(tensorname))
      
      io = nf90_put_var(ncid, VarID, tensor)
      call nc_check(io, 'wr_netcdf_2D', 'put_var: '//trim(tensorname))
      
      end subroutine wr_netcdf_2D
      
!=======================================================================
      
      subroutine get_ctim_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      
      integer, intent(in)  :: ncid
      integer, intent(out) :: lonDimID, latDimID
      integer, optional, intent(out) :: heightDimID
      
      integer :: io1, io2, io3
      
      ! We know the dimension names
      
      io1 = nf90_inq_dimid(ncid, 'cg_lon', lonDimID)
      io2 = nf90_inq_dimid(ncid, 'cg_lat', latDimID)
      
      call nc_check(io1, 'get_ctim_dimension_ids', 'inq_dimid cg_lon')
      call nc_check(io2, 'get_ctim_dimension_ids', 'inq_dimid cg_lat')
      
      if (present(heightDimID)) then
         io3 = nf90_inq_dimid(ncid, 'cg_height', heightDimID)
         call nc_check(io3, 'get_ctim_dimension_ids', 'inq_dimid cg_height')
      endif
      
      end subroutine get_ctim_dimension_ids
      
!=======================================================================
      
      subroutine get_interface_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      
      integer, intent(in)  :: ncid
      integer, intent(out) :: lonDimID, latDimID
      integer, optional, intent(out) :: heightDimID
      
      integer :: io1, io2, io3
      
      ! We know the dimension names
      
      io1 = nf90_inq_dimid(ncid, 'ig_lon', lonDimID)
      io2 = nf90_inq_dimid(ncid, 'ig_lat', latDimID)
      
      call nc_check(io1, 'get_interface_dimension_ids', 'inq_dimid ig_lon')
      call nc_check(io2, 'get_interface_dimension_ids', 'inq_dimid ig_lat')
      
      if (present(heightDimID)) then
         io3 = nf90_inq_dimid(ncid, 'ig_height', heightDimID)
         call nc_check(io3, 'get_interface_dimension_ids', 'inq_dimid ig_height')
      endif
      
      end subroutine get_interface_dimension_ids
      
!=======================================================================
      
      subroutine nc_check(istatus, subr_name, context)
      
      integer,          intent(in)           :: istatus
      character(len=*), intent(in)           :: subr_name
      character(len=*), intent(in), optional :: context
      
      character(len=512) :: error_msg
      
      ! if no error, nothing to do here.  we are done.
      if( istatus == nf90_noerr) return
      
      ! something wrong.  construct an error string and call the handler.
      ! context is optional, but is very useful if specified.
      if (present(context) ) then
          error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))
      else
          error_msg = nf90_strerror(istatus)
      endif
      
      write(*,*)'ERROR: '//trim(subr_name)//' '//trim(error_msg)
      stop
      
      end subroutine nc_check
      
!===================================================================
! End of netcdf module
!===================================================================
      
      end module netcdf_mod
      
