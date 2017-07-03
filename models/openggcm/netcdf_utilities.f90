! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! README ...
! To test this routine I needed to make a routine called 'death()' to match
! what is being used in openggcm. NETCDF_UTILITIES:DEATH() SHOULD NOT BE
! USED WHEN THIS MODULE IS USED WITH OPENGGCM.

module netcdf_utilities

implicit none
private

!----------------------------------------------------------------------------
! Attributes for variable kinds -- compiler flags are not portable
! 'digits12' is reserved for real variables that MUST retain 64 bits of
! precision. DO NOT CHANGE '12' to a smaller number. BAD BAD BAD things happen.
!----------------------------------------------------------------------------

integer, parameter :: i4 = SELECTED_INT_KIND(8)
integer, parameter :: i8 = SELECTED_INT_KIND(17)
integer, parameter :: r4 = SELECTED_REAL_KIND(6,30)
integer, parameter :: c4 = SELECTED_REAL_KIND(6,30)
integer, parameter :: r8 = SELECTED_REAL_KIND(12,100)
integer, parameter :: c8 = SELECTED_REAL_KIND(12,100)
integer, parameter :: digits12 = SELECTED_REAL_KIND(12)

public :: nc_check,                  &
          wr_netcdf_model_time,      &
          wr_netcdf_ctim_grid,       &
          wr_netcdf_interface_grid,  &
          wr_netcdf_oplus_grid,      &
          wr_netcdf,                 &
          rd_netcdf

interface wr_netcdf
   module procedure wr_netcdf_r4_1D
   module procedure wr_netcdf_r4_2D
   module procedure wr_netcdf_r4_3D
   module procedure wr_netcdf_r8_1D
   module procedure wr_netcdf_r8_2D
   module procedure wr_netcdf_r8_3D
end interface wr_netcdf

interface rd_netcdf
   module procedure rd_netcdf_r4_1D
   module procedure rd_netcdf_r4_2D
   module procedure rd_netcdf_r4_3D
   module procedure rd_netcdf_r8_1D
   module procedure rd_netcdf_r8_2D
   module procedure rd_netcdf_r8_3D
end interface rd_netcdf

contains

!=======================================================================

subroutine wr_netcdf_model_time(ncid, model_time)

use netcdf
implicit none

integer,  intent(in) :: ncid
real(r8), intent(in) :: model_time

integer :: VarID
integer :: io1, io2

io1 = nf90_def_var(ncid, 'time', nf90_double, VarID)
io2 = nf90_put_att(ncid, VarID, 'units', 'seconds since 1966-01-01')

call nc_check(io1, 'wr_netcdf_model_time', 'def_var time')
call nc_check(io2, 'wr_netcdf_model_time', 'put_att time units')

io1 = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io1, 'wr_netcdf_model_time', 'enddef')

io1 = nf90_put_var(ncid, VarID, model_time)
call nc_check(io1, 'wr_netcdf_model_time', 'put_var model_time')

end subroutine wr_netcdf_model_time

!=======================================================================

subroutine wr_netcdf_ctim_grid(ncid, nlon, nlonname, nlonunits, nlonshort, &
                                     nlat, nlatname, nlatunits, nlatshort, &
                 nheight, height, nheightname, nheightunits, nheightshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: nlon
character(len=*), intent(in) :: nlonname, nlonunits, nlonshort
integer,          intent(in) :: nlat
character(len=*), intent(in) :: nlatname, nlatunits, nlatshort
integer,          intent(in) :: nheight
real(r8),         intent(in) :: height(nheight)
character(len=*), intent(in) :: nheightname, nheightunits, nheightshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

real(r4), allocatable :: longitude(:)
real(r4), allocatable :: latitude(:)
real(r4)  :: dlat, dlon

integer :: LonDimID, LatDimID, HeightDimID
integer :: LonVarID, LatVarID, HeightVarID
integer :: io, io1, io2, io3
integer :: ilat, ilon

!-----------------------------------------------------------------------
! calculate the required metadata

allocate(longitude(nlon), latitude(nlat))

dlon = 360.0_r4/real(nlon-1,r4)
do ilon = 1,nlon
   longitude(ilon) = dlon*real(ilon-1,r4)
enddo

dlat = 180.0_r4/real(nlat-1,r4)
do ilat = 1,nlat           ! from south-to-north for the geographic grid
   latitude(ilat) = dlat*real(ilat-1,r4) - 90.0_r4
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

call nc_check(io1, 'wr_netcdf_ctim_grid', 'put_var cg_lon')
call nc_check(io2, 'wr_netcdf_ctim_grid', 'put_var cg_lat')
call nc_check(io3, 'wr_netcdf_ctim_grid', 'put_var cg_height')

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

real(r4), allocatable :: longitude(:)
real(r4), allocatable :: latitude(:)
real(r4), allocatable :: height(:)
real(r4)  :: dlat, dlon, dheight

integer :: LonDimID, LatDimID, HeightDimID
integer :: LonVarID, LatVarID, HeightVarID
integer :: io, io1, io2, io3
integer :: ilat, ilon, iheight

!=======================================================================

! calculate the required metadata with BOGUS VALUES at the moment
! for the height array

allocate(longitude(nlon), latitude(nlat), height(nheight))

dlon = 360.0_r4/real(nlon-1,r4)
do ilon = 1,nlon
   longitude(ilon) = dlon*real(ilon-1,r4)
enddo

dlat = 180.0_r4/real(nlat-1,r4)
do ilat = 1,nlat  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   latitude(ilat) = dlat*real(ilat-1,r4)
enddo

dheight = 500.0_r4/real(nheight,r4)   ! TJH FIXME ... bogus height array
do iheight = 1,nheight
   height(iheight) = dheight*real(iheight-1,r4)
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

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_interface_grid', 'enddef')

io1 = nf90_put_var(ncid,    LonVarID, longitude)
io2 = nf90_put_var(ncid,    LatVarID,  latitude)
io3 = nf90_put_var(ncid, HeightVarID,    height)

call nc_check(io1, 'wr_netcdf_interface_grid', 'put_var ig_lon')
call nc_check(io2, 'wr_netcdf_interface_grid', 'put_var ig_lat')
call nc_check(io3, 'wr_netcdf_interface_grid', 'put_var ig_height')

end subroutine wr_netcdf_interface_grid

!=======================================================================
!> Since the lon,lat,height variables are all 3-dimensional, all the
!> dimensions must be declared before any of the variables can be declared.

subroutine wr_netcdf_oplus_grid(ncid, nlon, nlonname, nlonunits, nlonshort, &
                                      nlat, nlatname, nlatunits, nlatshort, &
                                   nheight, nheightname, nheightunits, nheightshort, &
                                 oplus_lon, oplus_lat, oplus_height)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: nlon
character(len=*), intent(in) :: nlonname, nlonunits, nlonshort
integer,          intent(in) :: nlat
character(len=*), intent(in) :: nlatname, nlatunits, nlatshort
integer,          intent(in) :: nheight
character(len=*), intent(in) :: nheightname, nheightunits, nheightshort
real(r8),         intent(in) :: oplus_lon(   nheight,nlat,nlon)
real(r8),         intent(in) :: oplus_lat(   nheight,nlat,nlon)
real(r8),         intent(in) :: oplus_height(nheight,nlat,nlon)

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: LonDimID, LatDimID, HeightDimID
integer :: LonVarID, LatVarID, HeightVarID
integer :: io, io1, io2, io3

!-----------------------------------------------------------------------

io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_oplus_grid', 'redef')

io1 = nf90_def_dim(ncid, trim(nlonname),    nlon,    LonDimID)
io2 = nf90_def_dim(ncid, trim(nlatname),    nlat,    LatDimID)
io3 = nf90_def_dim(ncid, trim(nheightname), nheight, HeightDimID)
call nc_check(io1, 'wr_netcdf_oplus_grid', 'def_dim oplus_lon')
call nc_check(io2, 'wr_netcdf_oplus_grid', 'def_dim oplus_lat')
call nc_check(io3, 'wr_netcdf_oplus_grid', 'def_dim oplus_height')

! The dimensions must be specified in C-like order, i.e. reversed from their Fortran order.

io1 = nf90_def_var(ncid, trim(nlonname), nf90_real, (/ HeightDimID, LatDimID, LonDimID /), LonVarID)
io2 = nf90_put_att(ncid, LonVarID, 'short_name', trim(nlonshort) )
io3 = nf90_put_att(ncid, LonVarID, 'units', trim(nlonunits))

call nc_check(io1, 'wr_netcdf_oplus_grid', 'def_var oplus_lon')
call nc_check(io2, 'wr_netcdf_oplus_grid', 'put_att oplus_lon short_name')
call nc_check(io3, 'wr_netcdf_oplus_grid', 'put_att oplus_lon units')

io1 = nf90_def_var(ncid, trim(nlatname), nf90_real, (/ HeightDimID, LatDimID, LonDimID /), LatVarID)
io2 = nf90_put_att(ncid, LatVarID, 'short_name', trim(nlatshort))
io3 = nf90_put_att(ncid, LatVarID, 'units', trim(nlatunits))

call nc_check(io1, 'wr_netcdf_oplus_grid', 'def_var oplus_lat')
call nc_check(io2, 'wr_netcdf_oplus_grid', 'put_att oplus_lat short_name')
call nc_check(io3, 'wr_netcdf_oplus_grid', 'put_att oplus_lat units')

io1 = nf90_def_var(ncid, trim(nheightname), nf90_real, (/ HeightDimID, LatDimID, LonDimID /), HeightVarID)
io2 = nf90_put_att(ncid, HeightVarID, 'short_name', trim(nheightshort))
io3 = nf90_put_att(ncid, HeightVarID, 'units', trim(nheightunits))

call nc_check(io1, 'wr_netcdf_oplus_grid', 'def_var oplus_height')
call nc_check(io2, 'wr_netcdf_oplus_grid', 'put_att oplus_height short_name')
call nc_check(io3, 'wr_netcdf_oplus_grid', 'put_att oplus_height units')

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_oplus_grid', 'enddef')

io1 = nf90_put_var(ncid,    LonVarID, oplus_lon)
io2 = nf90_put_var(ncid,    LatVarID, oplus_lat)
io3 = nf90_put_var(ncid, HeightVarID, oplus_height)

call nc_check(io1, 'wr_netcdf_oplus_grid', 'put_var oplus_lon')
call nc_check(io2, 'wr_netcdf_oplus_grid', 'put_var oplus_lat')
call nc_check(io3, 'wr_netcdf_oplus_grid', 'put_var oplus_height')

end subroutine wr_netcdf_oplus_grid

!=======================================================================

subroutine wr_netcdf_r4_1D(ncid, dim1, dim1name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1
character(len=*), intent(in) :: dim1name
real(r4),         intent(in) :: tensor(dim1)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length

!-----------------------------------------------------------------------
! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
call nc_check(io1, 'wr_netcdf_r4_1D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
call nc_check(io1, 'wr_netcdf_r4_1D', 'inquire_dimension '//trim(dim1name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

!-----------------------------------------------------------------------
io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r4_1D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ Dim1_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r4_1D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_1D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r4_1D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r4_1D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid, VarID, tensor)
call nc_check(io, 'wr_netcdf_r4_1D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r4_1D

!=======================================================================

subroutine wr_netcdf_r8_1D(ncid, dim1, dim1name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1
character(len=*), intent(in) :: dim1name
real(r8),         intent(in) :: tensor(dim1)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length

!-----------------------------------------------------------------------
! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
call nc_check(io1, 'wr_netcdf_r8_1D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
call nc_check(io1, 'wr_netcdf_r8_1D', 'inquire_dimension '//trim(dim1name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

!-----------------------------------------------------------------------
io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r8_1D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_double, (/ Dim1_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r8_1D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_1D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r8_1D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r8_1D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid, VarID, tensor)
call nc_check(io, 'wr_netcdf_r8_1D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r8_1D


!=======================================================================

subroutine wr_netcdf_r4_2D(ncid, dim1, dim1name, dim2, dim2name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1, dim2
character(len=*), intent(in) :: dim1name, dim2name
real(r4),         intent(in) :: tensor(dim1,dim2)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, Dim2_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length, dim2length

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
io2 = nf90_inq_dimid(ncid, trim(dim2name), Dim2_ID)

call nc_check(io1, 'wr_netcdf_r4_2D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_2D', 'inq_dimid '//trim(dim2name)//' '//trim(tensorname))

! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
io2 = nf90_inquire_dimension(ncid, Dim2_ID, len=dim2length)

call nc_check(io1, 'wr_netcdf_r4_2D', 'inquire dimension length'//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_2D', 'inquire dimension length'//trim(dim2name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

if (dim2length /= size(tensor,2)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim2name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim2length,' expected ',size(tensor,2)
   call death()
endif

io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r4_2D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ Dim1_ID, Dim2_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r4_2D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_2D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r4_2D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r4_2D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid, VarID, tensor)
call nc_check(io, 'wr_netcdf_r4_2D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r4_2D

!=======================================================================

subroutine wr_netcdf_r8_2D(ncid, dim1, dim1name, dim2, dim2name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1, dim2
character(len=*), intent(in) :: dim1name, dim2name
real(r8),         intent(in) :: tensor(dim1,dim2)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, Dim2_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length, dim2length

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
io2 = nf90_inq_dimid(ncid, trim(dim2name), Dim2_ID)

call nc_check(io1, 'wr_netcdf_r8_2D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_2D', 'inq_dimid '//trim(dim2name)//' '//trim(tensorname))

! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
io2 = nf90_inquire_dimension(ncid, Dim2_ID, len=dim2length)

call nc_check(io1, 'wr_netcdf_r8_2D', 'inquire dimension length'//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_2D', 'inquire dimension length'//trim(dim2name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

if (dim2length /= size(tensor,2)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim2name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim2length,' expected ',size(tensor,2)
   call death()
endif

io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r8_2D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_double, (/ Dim1_ID, Dim2_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r8_2D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_2D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r8_2D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r8_2D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid, VarID, tensor)
call nc_check(io, 'wr_netcdf_r8_2D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r8_2D

!=======================================================================

subroutine wr_netcdf_r4_3D(ncid, dim1, dim1name, dim2, dim2name, dim3, dim3name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1, dim2, dim3
character(len=*), intent(in) :: dim1name, dim2name, dim3name
real(r4),         intent(in) :: tensor(dim1,dim2,dim3)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, Dim2_ID, Dim3_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length, dim2length, dim3length

!-----------------------------------------------------------------------

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
io2 = nf90_inq_dimid(ncid, trim(dim2name), Dim2_ID)
io3 = nf90_inq_dimid(ncid, trim(dim3name), Dim3_ID)

call nc_check(io1, 'wr_netcdf_r4_3D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_3D', 'inq_dimid '//trim(dim2name)//' '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r4_3D', 'inq_dimid '//trim(dim3name)//' '//trim(tensorname))

! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
io2 = nf90_inquire_dimension(ncid, Dim2_ID, len=dim2length)
io3 = nf90_inquire_dimension(ncid, Dim3_ID, len=dim3length)

call nc_check(io1, 'wr_netcdf_r4_3D', 'inquire dimension length'//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_3D', 'inquire dimension length'//trim(dim2name)//' '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r4_3D', 'inquire dimension length'//trim(dim3name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

if (dim2length /= size(tensor,2)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim2length,' expected ',size(tensor,2)
   call death()
endif

if (dim3length /= size(tensor,3)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim3length,' expected ',size(tensor,3)
   call death()
endif

! If everything is consistent, go ahead and define it

io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r4_3D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_real, (/ Dim1_ID, Dim2_ID, Dim3_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r4_3D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r4_3D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r4_3D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r4_3D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid,     VarID, tensor)
call nc_check(io, 'wr_netcdf_r4_3D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r4_3D

!=======================================================================

subroutine wr_netcdf_r8_3D(ncid, dim1, dim1name, dim2, dim2name, dim3, dim3name, tensor, tensorname, tensorunits, tensorshort)

use netcdf
implicit none

integer,          intent(in) :: ncid
integer,          intent(in) :: dim1, dim2, dim3
character(len=*), intent(in) :: dim1name, dim2name, dim3name
real(r8),         intent(in) :: tensor(dim1,dim2,dim3)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

integer :: Dim1_ID, Dim2_ID, Dim3_ID, VarID
integer :: io, io1, io2, io3
integer :: dim1length, dim2length, dim3length

!-----------------------------------------------------------------------

io1 = nf90_inq_dimid(ncid, trim(dim1name), Dim1_ID)
io2 = nf90_inq_dimid(ncid, trim(dim2name), Dim2_ID)
io3 = nf90_inq_dimid(ncid, trim(dim3name), Dim3_ID)

call nc_check(io1, 'wr_netcdf_r8_3D', 'inq_dimid '//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_3D', 'inq_dimid '//trim(dim2name)//' '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r8_3D', 'inq_dimid '//trim(dim3name)//' '//trim(tensorname))

! Make sure the netcdf dimensions match the incoming variable

io1 = nf90_inquire_dimension(ncid, Dim1_ID, len=dim1length)
io2 = nf90_inquire_dimension(ncid, Dim2_ID, len=dim2length)
io3 = nf90_inquire_dimension(ncid, Dim3_ID, len=dim3length)

call nc_check(io1, 'wr_netcdf_r8_3D', 'inquire dimension length'//trim(dim1name)//' '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_3D', 'inquire dimension length'//trim(dim2name)//' '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r8_3D', 'inquire dimension length'//trim(dim3name)//' '//trim(tensorname))

if (dim1length /= size(tensor,1)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim1length,' expected ',size(tensor,1)
   call death()
endif

if (dim2length /= size(tensor,2)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim2length,' expected ',size(tensor,2)
   call death()
endif

if (dim3length /= size(tensor,3)) then
   write(*,*)'ERROR :: failed dimension '//trim(dim1name)//' '//trim(tensorname)
   write(*,*)'ERROR :: got ',dim3length,' expected ',size(tensor,3)
   call death()
endif

! If everything is consistent, go ahead and define it

io = nf90_redef(ncid)
call nc_check(io, 'wr_netcdf_r8_3D', 'redef '//trim(tensorname))

io1 = nf90_def_var(ncid, tensorname, nf90_double, (/ Dim1_ID, Dim2_ID, Dim3_ID /), VarID)
io2 = nf90_put_att(ncid, VarID, 'short_name', trim(tensorshort))
io3 = nf90_put_att(ncid, VarID,      'units', trim(tensorunits))

call nc_check(io1, 'wr_netcdf_r8_3D', 'def_var '//trim(tensorname))
call nc_check(io2, 'wr_netcdf_r8_3D', 'put_att short_name '//trim(tensorname))
call nc_check(io3, 'wr_netcdf_r8_3D', 'put_att units '//trim(tensorname))

io = nf90_enddef(ncid) ! leave define mode so we can fill
call nc_check(io, 'wr_netcdf_r8_3D', 'enddef '//trim(tensorname))

io = nf90_put_var(ncid,     VarID, tensor)
call nc_check(io, 'wr_netcdf_r8_3D', 'put_var '//trim(tensorname))

end subroutine wr_netcdf_r8_3D

!=======================================================================
! Start of the read routines
!=======================================================================

subroutine rd_netcdf_r4_1D(ncid, tensorname, dim1, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1
real(r4),         intent(out) :: tensor(dim1)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r4_1D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if (dim1 /= ndim1) then
   write(*,*)'ERROR :: rd_netcdf_r4_1D: failed dimensions '//trim(tensorname),dim1
   write(*,*)'ERROR :: rd_netcdf_r4_1D: expected dimensions ',ndim1
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r4_1D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r4_1D

!=======================================================================

subroutine rd_netcdf_r8_1D(ncid, tensorname, dim1, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1
real(r8),         intent(out) :: tensor(dim1)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r8_1D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if (dim1 /= ndim1) then
   write(*,*)'ERROR :: rd_netcdf_r8_1D: failed dimensions '//trim(tensorname),dim1
   write(*,*)'ERROR :: rd_netcdf_r8_1D: expected dimensions ',ndim1
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r8_1D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r8_1D

!=======================================================================

subroutine rd_netcdf_r4_2D(ncid, tensorname, dim1, dim2, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1, dim2
real(r4),         intent(out) :: tensor(dim1,dim2)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r4_2D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if ((dim1 /= ndim1) .or. (dim2 /= ndim2)) then
   write(*,*)'ERROR :: rd_netcdf_r4_2D: failed dimensions '//trim(tensorname),dim1,dim2
   write(*,*)'ERROR :: rd_netcdf_r4_2D: expected dimensions ',ndim1,ndim2
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r4_2D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r4_2D

!=======================================================================

subroutine rd_netcdf_r8_2D(ncid, tensorname, dim1, dim2, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1, dim2
real(r8),         intent(out) :: tensor(dim1,dim2)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r8_2D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if ((dim1 /= ndim1) .or. (dim2 /= ndim2)) then
   write(*,*)'ERROR :: rd_netcdf_r8_2D: failed dimensions '//trim(tensorname),dim1,dim2
   write(*,*)'ERROR :: rd_netcdf_r8_2D: expected dimensions ',ndim1,ndim2
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r8_2D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r8_2D

!=======================================================================

subroutine rd_netcdf_r4_3D(ncid, tensorname, dim1, dim2, dim3, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1, dim2, dim3
real(r4),         intent(out) :: tensor(dim1,dim2,dim3)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r4_3D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if ((dim1 /= ndim1) .or. (dim2 /= ndim2) .or. (dim3 /= ndim3)) then
   write(*,*)'ERROR :: rd_netcdf_r4_3D: failed dimensions '//trim(tensorname),dim1,dim2,dim3
   write(*,*)'ERROR :: rd_netcdf_r4_3D: expected dimensions ',ndim1,ndim2,ndim3
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r4_3D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r4_3D

!=======================================================================

subroutine rd_netcdf_r8_3D(ncid, tensorname, dim1, dim2, dim3, tensor)

use netcdf
implicit none

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: tensorname
integer,          intent(in)  :: dim1, dim2, dim3
real(r8),         intent(out) :: tensor(dim1,dim2,dim3)

integer :: dim1ID, dim2ID, dim3ID, VarID
integer :: ndim1, ndim2, ndim3
integer :: io

io = nf90_inq_varid(ncid, tensorname, VarID)
call nc_check(io, 'rd_netcdf_r8_3D', 'inq_varid: '//trim(tensorname))

call get_dimensions(ncid, VarID, tensorname, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

if ((dim1 /= ndim1) .or. (dim2 /= ndim2) .or. (dim3 /= ndim3)) then
   write(*,*)'ERROR :: rd_netcdf_r8_3D: failed dimensions '//trim(tensorname),dim1,dim2,dim3
   write(*,*)'ERROR :: rd_netcdf_r8_3D: expected dimensions ',ndim1,ndim2,ndim3
   call death()
endif

io = nf90_get_var(ncid, VarID, tensor)
call nc_check(io, 'rd_netcdf_r8_3D', 'get_var: '//trim(tensorname))

end subroutine rd_netcdf_r8_3D

!=======================================================================

subroutine get_dimensions(ncid, VarID, variablename, dim1ID, ndim1, dim2ID, ndim2, dim3ID, ndim3)

use netcdf
implicit none

integer,          intent(in)  :: ncid
integer,          intent(in)  :: VarID
character(len=*), intent(in)  :: variablename ! for output messages
integer,          intent(out) :: dim1ID, ndim1
integer,          intent(out) :: dim2ID, ndim2
integer,          intent(out) :: dim3ID, ndim3

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens

integer :: io
integer :: numdims, i

character(len=512) :: string1

! set defaults
dimlens(:) = -1
dim1ID     = -1
dim2ID     = -1
dim3ID     = -1

io = nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=numdims)
call nc_check(io, 'get_dimensions', 'inquire variable '//trim(variablename))

if ((numdims > 3) .or. (numdims < 1)) then
   write(*,*)'ERROR: unexpected number of dimensions for '//trim(variablename)
   write(*,*)'ERROR: expected 1 or 2 or 3, got ',numdims
   call death()
endif

do i = 1,numdims
   write(string1,*)trim(variablename)//' dimension ',i
   io = nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i))
   call nc_check(io, 'get_dimensions ', trim(string1))
enddo

ndim1  = dimlens(1)
dim1ID =  dimids(1)

ndim2  = dimlens(2)
dim2ID =  dimids(2)

ndim3  = dimlens(3)
dim3ID =  dimids(3)

end subroutine get_dimensions

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

!=======================================================================
! ... REMOVE the DEATH() routine if you want to use the one in mod_corea.f
!=======================================================================

subroutine death()

stop

end subroutine death

!===================================================================

end module netcdf_utilities

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

