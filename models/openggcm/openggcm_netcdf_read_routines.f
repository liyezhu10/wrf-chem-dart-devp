! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

      module netcdf_read_mod

      implicit none
      private

      public :: rd_netcdf

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

      subroutine rd_netcdf_r4_1D(ncid, tensorname, dim1, tensor)

      use netcdf
      implicit none

      integer,          intent(in)  :: ncid
      character(len=*), intent(in)  :: tensorname
      integer,          intent(in)  :: dim1
      real*4,           intent(out) :: tensor(dim1)

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
      real*8,           intent(out) :: tensor(dim1)

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
      real*4,           intent(out) :: tensor(dim1,dim2)

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
      real*8,           intent(out) :: tensor(dim1,dim2)

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
      real*4,           intent(out) :: tensor(dim1,dim2,dim3)

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
      real*8,           intent(out) :: tensor(dim1,dim2,dim3)

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
! .... remove everything below here when putting into openggcm
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

      subroutine death()

      stop

      end subroutine death

!===================================================================

      end module netcdf_read_mod


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

