
subroutine wr_netcdf(dim1, dim1name, dim1short, &
                     dim2, dim2name, dim2short, &
                     dim3, dim3name, dim3short, &
                     tensor, tensorname, tensorunits, tensorshort, &
                     run_length_seconds, ncid)

use netcdf
implicit none

integer,          intent(in) :: dim1
character(len=*), intent(in) :: dim1name, dim1short
integer,          intent(in) :: dim2
character(len=*), intent(in) :: dim2name, dim2short
integer,          intent(in) :: dim3
character(len=*), intent(in) :: dim3name, dim3short
real*4,           intent(in) :: tensor(:,:,:)
character(len=*), intent(in) :: tensorname, tensorunits, tensorshort
real*8,           intent(in) :: run_length_seconds
integer,          intent(in) :: ncid

!-----------------------------------------------------------------------
! local storage
!-----------------------------------------------------------------------

real*4, allocatable :: latitude(:)
real*4, allocatable :: longitude(:)
real*8  :: time_offset_seconds
real*4  :: dlat, dlon

integer :: ncid, dim1_ID, dim2_ID, dim3_ID, VarID, 
integer :: LonVarID, LatVarID, LevVarID, TimeVarID
integer :: io, iunit

integer :: nthe, nphi
integer :: ilat, ilon

!=======================================================================

!-----------------------------------------------------------------------
! calculate the required metadata

nthe = -1
nphi = -1

if (trim(dim1name) == 'nthe') nthe = dim1
if (trim(dim2name) == 'nthe') nthe = dim2
if (trim(dim1name) == 'nphi') nphi = dim1
if (trim(dim2name) == 'nphi') nphi = dim2

if (nthe < 1) then
   write(*,*)'ERROR: unable to determine nthe.'
   stop
endif

if (nphi < 1) then
   write(*,*)'ERROR: unable to determine nphi.'
   stop
endif

allocate(latitude(nthe), longitude(nphi))

dlat = 180.0/real(nthe-1)
do ilat = 1,nthe  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   latitude(ilat) = dlat*real(ilat-1)
enddo

dlon = 360.0/real(nphi-1)
do ilon = 1,nphi  ! NORTH-to-SOUTH, NORTH is magnetic latitude 0.0
   longitude(ilon) = dlon*real(ilon-1)
enddo

! call nc_check(nf90_create(filename, NF90_CLOBBER, ncid),'wr_netcdf')
!-----------------------------------------------------------------------

io = nf90_def_var(ncid, time, nf90_double, TimeVarID)
call nc_check(io, 'wr_netcdf', 'def_var time')
io = nf90_put_att(ncid, TimeVarID, 'units', 'seconds since 1966-01-01')
call nc_check(io, 'wr_netcdf', 'put_att:time:units')

io = nf90_def_dim(ncid, dim1name, dim1, dim1_ID)
call nc_check(io, 'wr_netcdf', 'def_dim '//trim(dim1name))
io = nf90_def_var(ncid, dim1name, nf90_real, (/ dim1_ID /), LonVarID)
call nc_check(io, 'wr_netcdf', 'def_var:'//trim(dim1name))
io = nf90_put_att(ncid, LonVarID, 'short_name', dim1short)
call nc_check(io, 'wr_netcdf', 'put_att short_name '//trim(dim1name))

io = nf90_def_dim(ncid, dim2name, dim2, dim2_ID)
call nc_check(io, 'wr_netcdf', 'def_dim '//trim(dim2name))
io = nf90_def_var(ncid, dim2name, nf90_real, (/ dim2_ID /), LatVarID)
call nc_check(io, 'wr_netcdf', 'def_var:'//trim(dim2name))
io = nf90_put_att(ncid, LatVarID, 'short_name', dim2short)
call nc_check(io, 'wr_netcdf', 'put_att short_name '//trim(dim2name))

if (dim3 > 1) then
   io = nf90_def_dim(ncid, dim3name, dim3, dim3_ID)
   call nc_check(io, 'wr_netcdf', 'def_dim '//trim(dim3name))
   io = nf90_def_var(ncid, dim3name, nf90_real, (/ dim3_ID /), LevVarID)
   call nc_check(io, 'wr_netcdf', 'def_var:'//trim(dim3name))
   io = nf90_put_att(ncid, LevVarID, 'short_name', dim3short)
   call nc_check(io, 'wr_netcdf', 'put_att short_name '//trim(dim3name))
endif

! Now define the variable/tensor/hyperslab

if (dim3 > 1) then
   io = nf90_def_var(ncid, tensorname, nf90_real, (/ dim1_ID, dim2_ID, dim3_ID /), VarID)
   call nc_check(io, 'wr_netcdf', 'def_var:'//trim(tensorname))
else
   io = nf90_def_var(ncid, tensorname, nf90_real, (/ dim1_ID, dim2_ID /), VarID)
   call nc_check(io, 'wr_netcdf', 'def_var:'//trim(tensorname))
endif

io = nf90_put_att(ncid, VarID, 'short_name', tensorshort)
call nc_check(io, 'wr_netcdf', 'put_att short_name'//trim(tensorname))
io = nf90_put_att(ncid, VarID,      'units', tensorunits)
call nc_check(io, 'wr_netcdf', 'put_att units'//trim(tensorname))

! leave define mode so we can fill

io = nf90_enddef(ncid)
call nc_check(io, 'wr_netcdf', 'enddef')

io = nf90_put_var(ncid,  LonVarID, longitude)
call nc_check(io, 'wr_netcdf', 'put_var: lon')
io = nf90_put_var(ncid,  LatVarID, latitude)
call nc_check(io, 'wr_netcdf', 'put_var: lat')
io = nf90_put_var(ncid,     VarID, tensor)
call nc_check(io, 'wr_netcdf', 'put_var: tensor')
io = nf90_put_var(ncid, TimeVarID, time_offset_seconds)
call nc_check(io, 'wr_netcdf', 'put_var: time')

! io = nf90_close(ncid)
! call nc_check(io, 'wr_netcdf', 'nf90_close')

contains 

   subroutine nc_check(istatus, subr_name, context)
      integer, intent (in)                   :: istatus
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
      write(*,*)'ERROR: '//trim(subr_name)//' '//trim(error_msg)
      write(*,*)'ERROR: '//trim(subr_name)//' '//trim(error_msg)
      stop

   end subroutine nc_check

end subroutine wr_netcdf
