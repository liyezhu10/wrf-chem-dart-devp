! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module HDF5_utilities_mod

use        types_mod, only : i2, i4, r4, r8, MISSING_R8, MISSING_I
use    utilities_mod, only : nc_check, E_MSG, E_ERR, error_handler
use time_manager_mod, only : time_type, operator(>=), set_time, get_time

use HDF5

implicit none
private

public :: H5_CRTDAT, H5_RDWT, h5_get_rank, h5_get_dimensions


! interface hf_get_var
!    module procedure hf_get_int_1d
!    module procedure hf_get_real_1d
! end interface

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2, string3


contains


function h5_get_rank(dspace_id, error) result(rank)

integer(HID_T), intent(in)  :: dspace_id
integer,        intent(out) :: error
integer :: rank

call h5sget_simple_extent_ndims_f(dspace_id, rank, error)

write(*,*)'TJH rank is ',rank

end function h5_get_rank



subroutine h5_get_dimensions(dspace_id, dims, error)

integer(HID_T),   intent(in)  :: dspace_id
integer(HSIZE_T), intent(out) :: dims(:)
integer,          intent(out) :: error

integer(HSIZE_T) :: maxdims(size(dims))

! get the dimensions of the dataspace
call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)

write(*,*)'TJH    dims is ',dims
write(*,*)'TJH maxdims is ',maxdims

end subroutine h5_get_dimensions



subroutine H5_CRTDAT()

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This routine is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! The following example shows how to create an empty dataset.
! It creates a file called 'dsetf.h5', defines the
! dataset dataspace, creates a dataset which is a 4x6 integer array,
! and then closes the dataspace, the dataset, and the file.
!
! This example is used in the HDF5 Tutorial.

character(len=8), parameter :: filename = "dsetf.h5" ! File name
character(len=4), parameter :: dsetname = "dset"     ! Dataset name

integer(HID_T) :: file_id       ! File identifier
integer(HID_T) :: dset_id       ! Dataset identifier
integer(HID_T) :: dspace_id     ! Dataspace identifier

integer(HSIZE_T), dimension(2) :: dims = (/4,6/) ! Dataset dimensions
integer     ::   rank = 2                        ! Dataset rank

integer     ::   error ! Error flag

! Initialize FORTRAN interface.
call h5open_f(error)

! Create a new file using default properties.
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

! Create the dataspace.
call h5screate_simple_f(rank, dims, dspace_id, error)

! Create the dataset with default properties.
call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

! End access to the dataset and release resources used by it.
call h5dclose_f(dset_id, error)

! Terminate access to the data space.
call h5sclose_f(dspace_id, error)

! Close the file.
call h5fclose_f(file_id, error)

! Close FORTRAN interface.
call h5close_f(error)

end subroutine H5_CRTDAT


subroutine H5_RDWT

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This routine is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the COPYING file, which can be found at the root of the source code       *
!   distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
!   If you do not have access to either file, you may request a copy from     *
!   help@hdfgroup.org.                                                        *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! The following example shows how to write and read to/from an existing dataset.
! It opens the file created in the previous example, obtains the dataset
! identifier, writes the data to the dataset in the file,
! then reads the dataset  to memory.
!
! This example is used in the HDF5 Tutorial.

! Initialize the dset_data array.

character(len=8), parameter :: filename = "dsetf.h5" ! File name
character(len=4), parameter :: dsetname = "dset"     ! Dataset name

integer(HID_T) :: file_id       ! File identifier
integer(HID_T) :: dset_id       ! Dataset identifier

integer :: error ! Error flag
integer :: i, j

integer, dimension(4,6) :: dset_data, data_out ! Data buffers
integer(HSIZE_T), dimension(2) :: data_dims

DO i = 1, 4
   DO j = 1, 6
      dset_data(i,j) = (i-1)*6 + j
   END DO
END DO


! Initialize FORTRAN interface.
call h5open_f(error)

! Open an existing file.
call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)

! Open an existing dataset.
call h5dopen_f(file_id, dsetname, dset_id, error)

! Write the dataset.

data_dims(1) = 4
data_dims(2) = 6
call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)

! Read the dataset.
call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, data_dims, error)

! Close the dataset.
call h5dclose_f(dset_id, error)

! Close the file.
call h5fclose_f(file_id, error)

! Close FORTRAN interface.
call h5close_f(error)

end subroutine H5_RDWT



end module HDF5_utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
