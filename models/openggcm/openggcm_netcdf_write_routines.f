
      module netcdf_mod

      use netcdf
      implicit none
      private

      public :: wr_netcdf_model_time,
     &          wr_netcdf_ctim_grid,
     &          wr_netcdf_interface_grid,
     &          wr_netcdf_r4_1D,
     &          wr_netcdf_r4_2D,
     &          wr_netcdf_r8_2D,
     &          wr_netcdf_r4_3D,
     &          wr_netcdf_r8_3D

      contains

!===================================================================

      subroutine death()

      stop

      end subroutine death

!=======================================================================
! This is F90 code written with F77 syntax because (for expediency) it will
! simply be appended to the mhd-iono.for file.
! This means all line continuation characters are in column 6
! Appently compiles with long lines, so can go beyond column 72
!
! Furthermore, the build mechanism of openggcm cannot tolerate having a
! module contained inside the mhd-iono.for   file, so the module baggage
! must be stripped off when appending to mhd-iono.for.
!
! ... delete everything above this line when appending
!=======================================================================

      subroutine wr_netcdf_model_time(ncid, model_time)

      use netcdf
      implicit none

      integer, intent(in) :: ncid
      real*8,  intent(in) :: model_time

      integer :: VarID
      integer :: io1, io2

      io1 = nf90_def_var(ncid, 'time', nf90_double, VarID)
      io2 = nf90_put_att(ncid, VarID, 'units', 'seconds since 1966-01-01')

      call nc_check(io1, 'wr_netcdf_model_time', 'def_var time')
      call nc_check(io2, 'wr_netcdf_model_time', 'put_att:time:units')

      io1 = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io1, 'wr_netcdf_model_time', 'enddef')

      io1 = nf90_put_var(ncid, VarID, model_time)
      call nc_check(io1, 'wr_netcdf_model_time', 'put_var: model_time')

      end subroutine wr_netcdf_model_time

!=======================================================================

      subroutine wr_netcdf_ctim_grid(ncid, nlon, nlonname, nlonunits, nlonshort,
     &                                     nlat, nlatname, nlatunits, nlatshort,
     &                 nheight, height, nheightname, nheightunits, nheightshort)

      use netcdf
      implicit none

      integer, intent(in) :: ncid
      integer, intent(in) :: nlon
      character(len=*), intent(in) :: nlonname, nlonunits, nlonshort
      integer, intent(in) :: nlat
      character(len=*), intent(in) :: nlatname, nlatunits, nlatshort
      integer, intent(in) :: nheight
      real*4,  intent(in) :: height(nheight)
      character(len=*), intent(in) :: nheightname, nheightunits, nheightshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      real*4, allocatable :: longitude(:)
      real*4, allocatable :: latitude(:)
      real*4  :: dlat, dlon

      integer :: LonDimID, LatDimID, HeightDimID
      integer :: LonVarID, LatVarID, HeightVarID
      integer :: io, io1, io2, io3
      integer :: ilat, ilon

      !-----------------------------------------------------------------------
      ! calculate the required metadata

      allocate(longitude(nlon), latitude(nlat))

      dlon = 360.0/real(nlon-1)
      do ilon = 1,nlon
         longitude(ilon) = dlon*real(ilon-1)
      enddo

      dlat = 180.0/real(nlat-1)
      do ilat = 1,nlat           ! from south-to-north
         latitude(ilat) = dlat*real(ilat-1) - 90.0
      enddo

      !-----------------------------------------------------------------------

      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_ctim_grid', 'redef')

      io1 = nf90_def_dim(ncid, trim(nlonname), nlon, LonDimID)
      io2 = nf90_def_var(ncid, trim(nlonname), nf90_real, (/ LonDimID /), LonVarID)
      io3 = nf90_put_att(ncid, LonVarID, 'short_name', trim(nlonshort) )
      io  = nf90_put_att(ncid, LonVarID, 'units', trim(nlonunits))

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_lon')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_lon')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_lon short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_lon units')

      io1 = nf90_def_dim(ncid, trim(nlatname), nlat, LatDimID)
      io2 = nf90_def_var(ncid, trim(nlatname), nf90_real, (/ LatDimID /), LatVarID)
      io3 = nf90_put_att(ncid, LatVarID, 'short_name', trim(nlatshort))
      io  = nf90_put_att(ncid, LatVarID, 'units', trim(nlatunits))

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_lat')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_lat')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_lat short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_lat units')

      io1 = nf90_def_dim(ncid, trim(nheightname), nheight, HeightDimID)
      io2 = nf90_def_var(ncid, trim(nheightname), nf90_real, (/ HeightDimID /), HeightVarID)
      io3 = nf90_put_att(ncid, HeightVarID, 'short_name', trim(nheightshort))
      io  = nf90_put_att(ncid, HeightVarID, 'units', trim(nheightunits))

      call nc_check(io1, 'wr_netcdf_ctim_grid', 'def_dim cg_height')
      call nc_check(io2, 'wr_netcdf_ctim_grid', 'def_var cg_height')
      call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_att cg_height short_name')
      call nc_check(io , 'wr_netcdf_ctim_grid', 'put_att cg_height units')

      io = nf90_enddef(ncid) ! leave define mode so we can fill
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

      use netcdf
      implicit none

      integer, intent(in) :: ncid
      integer, intent(in) :: nlon
      integer, intent(in) :: nlat
      integer, intent(in) :: nheight

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      real*4, allocatable :: longitude(:)
      real*4, allocatable :: latitude(:)
      real*4, allocatable :: height(:)
      real*4  :: dlat, dlon, dheight

      integer :: LonDimID, LatDimID, HeightDimID
      integer :: LonVarID, LatVarID, HeightVarID
      integer :: io, io1, io2, io3
      integer :: ilat, ilon, iheight

      !=======================================================================

      ! calculate the required metadata with BOGUS VALUES at the moment

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

      io1 = nf90_def_dim(ncid, 'ig_lon'   , nlon   , LonDimID)
      io2 = nf90_def_dim(ncid, 'ig_lat'   , nlat   , LatDimID)
      io3 = nf90_def_dim(ncid, 'ig_height', nheight, HeightDimID)

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_dim ig_lon')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'def_dim ig_lat')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'def_dim ig_height')

      io1 = nf90_def_var(ncid, 'ig_lon', nf90_real, (/ LonDimID /), LonVarID)
      io2 = nf90_put_att(ncid, LonVarID, 'short_name', 'magnetic longitude' )
      io3 = nf90_put_att(ncid, LonVarID, 'units', 'degrees' )

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_var ig_lon')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'put_att ig_lon short_name')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_lon units')

      io1 = nf90_def_var(ncid, 'ig_lat', nf90_real, (/ LatDimID /), LatVarID)
      io2 = nf90_put_att(ncid, LatVarID, 'short_name', 'magnetic colatitude')
      io3 = nf90_put_att(ncid, LatVarID, 'units', 'degrees' )

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_var ig_lat')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'put_att ig_lat short_name')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_lat units')

      io1 = nf90_def_var(ncid, 'ig_height', nf90_real, (/ HeightDimID /), HeightVarID)
      io2 = nf90_put_att(ncid, HeightVarID, 'short_name', 'height')
      io3 = nf90_put_att(ncid, HeightVarID, 'units', 'kilometers' )

      call nc_check(io1, 'wr_netcdf_interface_grid', 'def_var ig_height')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'put_att ig_height short_name')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_att ig_height units')

      ! must also put out some information so DART can convert magnetic lon/lat to
      ! geographic lon/lat on demand.

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_interface_grid', 'enddef')

      io1 = nf90_put_var(ncid,    LonVarID, longitude)
      io2 = nf90_put_var(ncid,    LatVarID,  latitude)
      io3 = nf90_put_var(ncid, HeightVarID,    height)

      call nc_check(io1, 'wr_netcdf_interface_grid', 'put_var: ig_lon')
      call nc_check(io2, 'wr_netcdf_interface_grid', 'put_var: ig_lat')
      call nc_check(io3, 'wr_netcdf_interface_grid', 'put_var: ig_height')

      end subroutine wr_netcdf_interface_grid

!=======================================================================

      subroutine wr_netcdf_r4_1D(ncid, gridtype, dim1, tensor, tensorname, tensorunits, tensorshort)

      use netcdf
      implicit none

      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      integer,          intent(in) :: dim1
      real*4,           intent(in) :: tensor(dim1)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      integer :: LonDimID, LatDimID, LevDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon

      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      else
         write(*,*)"ERROR :: unknown grid type '"//trim(gridtype)//"'"
         call death()
      endif

      ! Make sure the netcdf dimensions match the incoming variable

      io1 = nf90_inquire_dimension(ncid, LonDimID, len=nlon)
      call nc_check(io1, 'wr_netcdf_r4_1D', 'inquire lon dimension length'//trim(tensorname))

      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed dimension '//trim(tensorname)
         call death()
      endif

      !-----------------------------------------------------------------------
      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_r4_1D', 'redef'//trim(tensorname))

      io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ LonDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_r4_1D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r4_1D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_r4_1D', 'put_att units'//trim(tensorname))

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_r4_1D', 'enddef'//trim(tensorname))

      io = nf90_put_var(ncid, VarID, tensor)
      call nc_check(io, 'wr_netcdf_r4_1D', 'put_var: '//trim(tensorname))

      end subroutine wr_netcdf_r4_1D

!=======================================================================

      subroutine wr_netcdf_r4_2D(ncid, gridtype, dim1, dim2, tensor, tensorname, tensorunits, tensorshort)

      use netcdf
      implicit none

      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      integer,          intent(in) :: dim1, dim2
      real*4,           intent(in) :: tensor(dim1,dim2)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      integer :: LonDimID, LatDimID, LevDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon, nlat

      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      else
         write(*,*)"ERROR :: unknown grid type '"//trim(gridtype)//"'"
         call death()
      endif

      ! Make sure the netcdf dimensions match the incoming variable

      io1 = nf90_inquire_dimension(ncid, LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid, LatDimID, len=nlat)

      call nc_check(io1, 'wr_netcdf_r4_2D', 'inquire lon dimension length'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r4_2D', 'inquire lat dimension length'//trim(tensorname))

      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlon,' got ',size(tensor,1)
         call death()
      endif

      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlat,' got ',size(tensor,2)
         call death()
      endif

      !-----------------------------------------------------------------------
      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_r4_2D', 'redef'//trim(tensorname))

      io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ LonDimID, LatDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_r4_2D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r4_2D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_r4_2D', 'put_att units'//trim(tensorname))

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_r4_2D', 'enddef'//trim(tensorname))

      io = nf90_put_var(ncid, VarID, tensor)
      call nc_check(io, 'wr_netcdf_r4_2D', 'put_var: '//trim(tensorname))

      end subroutine wr_netcdf_r4_2D

!=======================================================================

      subroutine wr_netcdf_r8_2D(ncid, gridtype, dim1, dim2, tensor, tensorname, tensorunits, tensorshort)

      use netcdf
      implicit none

      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      integer,          intent(in) :: dim1, dim2
      real*8,           intent(in) :: tensor(dim1,dim2)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      integer :: LonDimID, LatDimID, LevDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon, nlat

      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, LevDimID)
      else
         write(*,*)"ERROR :: unknown grid type '"//trim(gridtype)//"'"
         call death()
      endif

      ! Make sure the netcdf dimensions match the incoming variable

      io1 = nf90_inquire_dimension(ncid, LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid, LatDimID, len=nlat)

      call nc_check(io1, 'wr_netcdf_r8_2D', 'inquire lon dimension length'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r8_2D', 'inquire lat dimension length'//trim(tensorname))

      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlon,' got ',size(tensor,1)
         call death()
      endif

      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlat,' got ',size(tensor,2)
         call death()
      endif

      !-----------------------------------------------------------------------
      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_r8_2D', 'redef'//trim(tensorname))

      io1 = nf90_def_var(ncid, tensorname, nf90_double, (/ LonDimID, LatDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_r8_2D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r8_2D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_r8_2D', 'put_att units'//trim(tensorname))

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_r8_2D', 'enddef'//trim(tensorname))

      io = nf90_put_var(ncid, VarID, tensor)
      call nc_check(io, 'wr_netcdf_r8_2D', 'put_var: '//trim(tensorname))

      end subroutine wr_netcdf_r8_2D

!=======================================================================

      subroutine wr_netcdf_r4_3D(ncid, gridtype, dim1, dim2, dim3, tensor, tensorname, tensorunits, tensorshort)

      use netcdf
      implicit none

      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      integer,          intent(in) :: dim1, dim2, dim3
      real*4,           intent(in) :: tensor(dim1,dim2,dim3)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      integer :: LonDimID, LatDimID, HeightDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon, nlat, nheight

      !-----------------------------------------------------------------------

      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      else
         write(*,*)"ERROR :: unknown grid type '"//trim(gridtype)//"'"
         call death()
      endif

      ! Make sure the netcdf dimensions match the incoming variable

      io1 = nf90_inquire_dimension(ncid,    LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid,    LatDimID, len=nlat)
      io3 = nf90_inquire_dimension(ncid, HeightDimID, len=nheight)

      call nc_check(io1, 'wr_netcdf_r4_3D', 'inquire lon    dimension length')
      call nc_check(io2, 'wr_netcdf_r4_3D', 'inquire lat    dimension length')
      call nc_check(io3, 'wr_netcdf_r4_3D', 'inquire height dimension length')

      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlon,' got ',size(tensor,1)
         call death()
      endif

      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlat,' got ',size(tensor,2)
         call death()
      endif

      if (nheight /= size(tensor,3)) then
         write(*,*)'ERROR :: failed height dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nheight,' got ',size(tensor,3)
         call death()
      endif

      ! If everything is consistent, go ahead and define it

      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_r4_3D', 'redef '//trim(tensorname))

      ! Now define the variable/tensor/hyperslab

      io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ LonDimID, LatDimID, HeightDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_r4_3D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r4_3D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_r4_3D', 'put_att units'//trim(tensorname))

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_r4_3D', 'enddef')

      io = nf90_put_var(ncid,     VarID, tensor)
      call nc_check(io, 'wr_netcdf_r4_3D', 'put_var: tensor')

      end subroutine wr_netcdf_r4_3D

!=======================================================================

      subroutine wr_netcdf_r8_3D(ncid, gridtype, dim1, dim2, dim3, tensor, tensorname, tensorunits, tensorshort)

      use netcdf
      implicit none

      integer,          intent(in) :: ncid
      character(len=*), intent(in) :: gridtype ! 'ctim' or 'interface'
      integer,          intent(in) :: dim1, dim2, dim3
      real*8,           intent(in) :: tensor(dim1,dim2,dim3)
      character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

      !-----------------------------------------------------------------------
      ! local storage
      !-----------------------------------------------------------------------

      integer :: LonDimID, LatDimID, HeightDimID, VarID
      integer :: io, io1, io2, io3
      integer :: nlon, nlat, nheight

      !-----------------------------------------------------------------------

      if (trim(gridtype) == 'ctim') then
         call  get_ctim_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      elseif (trim(gridtype) == 'interface') then
         call  get_interface_dimension_ids(ncid, lonDimID, latDimID, heightDimID)
      else
         write(*,*)"ERROR :: unknown grid type '"//trim(gridtype)//"'"
         call death()
      endif

      ! Make sure the netcdf dimensions match the incoming variable

      io1 = nf90_inquire_dimension(ncid,    LonDimID, len=nlon)
      io2 = nf90_inquire_dimension(ncid,    LatDimID, len=nlat)
      io3 = nf90_inquire_dimension(ncid, HeightDimID, len=nheight)

      call nc_check(io1, 'wr_netcdf_r8_3D', 'inquire lon    dimension length')
      call nc_check(io2, 'wr_netcdf_r8_3D', 'inquire lat    dimension length')
      call nc_check(io3, 'wr_netcdf_r8_3D', 'inquire height dimension length')

      if (nlon /= size(tensor,1)) then
         write(*,*)'ERROR :: failed longitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlon,' got ',size(tensor,1)
         call death()
      endif

      if (nlat /= size(tensor,2)) then
         write(*,*)'ERROR :: failed latitude dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nlat,' got ',size(tensor,2)
         call death()
      endif

      if (nheight /= size(tensor,3)) then
         write(*,*)'ERROR :: failed height dimension '//trim(tensorname)
         write(*,*)'ERROR :: expected ',nheight,' got ',size(tensor,3)
         call death()
      endif

      ! If everything is consistent, go ahead and define it

      io = nf90_redef(ncid)
      call nc_check(io, 'wr_netcdf_r8_3D', 'redef '//trim(tensorname))

      ! Now define the variable/tensor/hyperslab

      io1 = nf90_def_var(ncid, tensorname, nf90_double, (/ LonDimID, LatDimID, HeightDimID /), VarID)
      io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
      io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

      call nc_check(io1, 'wr_netcdf_r8_3D', 'def_var:'//trim(tensorname))
      call nc_check(io2, 'wr_netcdf_r8_3D', 'put_att short_name'//trim(tensorname))
      call nc_check(io3, 'wr_netcdf_r8_3D', 'put_att units'//trim(tensorname))

      io = nf90_enddef(ncid) ! leave define mode so we can fill
      call nc_check(io, 'wr_netcdf_r8_3D', 'enddef')

      io = nf90_put_var(ncid,     VarID, tensor)
      call nc_check(io, 'wr_netcdf_r8_3D', 'put_var: tensor')

      end subroutine wr_netcdf_r8_3D

!=======================================================================

      subroutine get_ctim_dimension_ids(ncid, lonDimID, latDimID, heightDimID)

      use netcdf
      implicit none

      integer, intent(in)  :: ncid
      integer, intent(out) :: lonDimID, latDimID
      integer, intent(out) :: heightDimID

      integer :: io1, io2, io3

      ! We know the dimension names

      io1 = nf90_inq_dimid(ncid, 'cg_lon', lonDimID)
      io2 = nf90_inq_dimid(ncid, 'cg_lat', latDimID)
      io3 = nf90_inq_dimid(ncid, 'cg_height', heightDimID)

      call nc_check(io1, 'get_ctim_dimension_ids', 'inq_dimid cg_lon')
      call nc_check(io2, 'get_ctim_dimension_ids', 'inq_dimid cg_lat')
      call nc_check(io3, 'get_ctim_dimension_ids', 'inq_dimid cg_height')

      end subroutine get_ctim_dimension_ids

!=======================================================================

      subroutine get_interface_dimension_ids(ncid, lonDimID, latDimID, heightDimID)

      use netcdf
      implicit none

      integer, intent(in)  :: ncid
      integer, intent(out) :: lonDimID, latDimID
      integer, intent(out) :: heightDimID

      integer :: io1, io2, io3

      ! We know the dimension names

      io1 = nf90_inq_dimid(ncid, 'ig_lon', lonDimID)
      io2 = nf90_inq_dimid(ncid, 'ig_lat', latDimID)
      io3 = nf90_inq_dimid(ncid, 'ig_height', heightDimID)

      call nc_check(io1, 'get_interface_dimension_ids', 'inq_dimid ig_lon')
      call nc_check(io2, 'get_interface_dimension_ids', 'inq_dimid ig_lat')
      call nc_check(io3, 'get_interface_dimension_ids', 'inq_dimid ig_height')

      end subroutine get_interface_dimension_ids

!=======================================================================

      subroutine nc_check(istatus, subr_name, context)

      use netcdf
      implicit none

      integer,          intent(in) :: istatus
      character(len=*), intent(in) :: subr_name
      character(len=*), intent(in) :: context

      character(len=512) :: error_msg

      ! if no error, nothing to do here.  we are done.
      if( istatus == nf90_noerr) return

      ! something wrong.  construct an error string and call the handler.
      error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))

      write(*,*)'ERROR: '//trim(subr_name)//' '//trim(error_msg)
      call death()

      end subroutine nc_check

!===================================================================
! End of netcdf module
!===================================================================

      end module netcdf_mod

