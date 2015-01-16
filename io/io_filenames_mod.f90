! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module io_filenames_mod

!> \defgroup io_filenames io_filenames
!> Aim is to store the io filenames
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Any module can set the filenames here, then state_vector_io_mod
!> can read from this module to get the filenames.
!> Maybe this is a bit lazy, but I just want a way for the 
!> different modules to set the filenames
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

use utilities_mod, only : do_nml_file, nmlfileunit, do_nml_term, check_namelist_read, &
                          find_namelist_in_file
use model_mod,     only : construct_file_name_in, construct_file_name_out

implicit none

private

! These should probably be set and get functions rather than 
! direct access

public :: io_filenames_init, restart_files_in, restart_files_out

! How do people name there restart files?
! What about domains?
integer, parameter :: max_num_files = 5000

! public arrays of filenames. Do we need arrays for restarts AND extras?
character(len=2048), allocatable :: restart_files_in(:,:), restart_files_out(:,:,:)

! Namelist options
character(len=512) :: restart_in_stub  = 'input'
character(len=512) :: restart_out_stub = 'output'


! Should probably get num_domains, num_restarts from elsewhere. In here for now
namelist / io_filenames_nml / restart_in_stub, restart_out_stub

contains

!----------------------------------
!> read namelist and set up filename arrays
subroutine io_filenames_init(ens_size, num_domains, inflation_in, inflation_out)

integer, intent(in) :: ens_size ! ensemble size + extras
integer, intent(in) :: num_domains
integer :: iunit, io
integer :: dom, num_files, i
character(len = 129) :: inflation_in(2), inflation_out(2)
character(len = 4)   :: extension

!call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "io_filenames_nml", iunit)
read(iunit, nml = io_filenames_nml, iostat = io)
call check_namelist_read(iunit, io, "io_filenames_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=io_filenames_nml)
if (do_nml_term()) write(     *     , nml=io_filenames_nml)

num_files = ens_size + 6 !> @todo

allocate(restart_files_in(num_files, num_domains))
allocate(restart_files_out(num_files, num_domains, 2)) ! for prior and posterior filenames

do dom = 1, num_domains
   do i = 1, ens_size  ! restarts
      restart_files_in(i, dom)  = construct_file_name_in(restart_in_stub, dom, i)
      restart_files_out(i, dom, 1) = construct_file_name_out(restart_out_stub, dom, i)
      write(extension, '(i4.4)') i
      restart_files_out(i, dom, 2) = 'prior_member' // extension
   enddo
enddo

! input extras
do dom = 1, num_domains
   ! mean -never used
   write(restart_files_in(ens_size + 1, dom), '(A, i2.2, A)') 'mean_d', dom
   ! sd -never used
   write(restart_files_in(ens_size + 2, dom), '(A, i2.2, A)') 'sd_d',   dom
   ! prior inf copy
   write(restart_files_in(ens_size + 3, dom), '(A, A, i2.2, A)') trim(inflation_in(1)), '_mean_d', dom
   ! prior inf sd copy
   write(restart_files_in(ens_size + 4, dom), '(A, A, i2.2, A)') trim(inflation_in(1)), '_sd_d', dom
   ! post inf copy
   write(restart_files_in(ens_size + 5, dom), '(A, A, i2.2, A)') trim(inflation_in(2)), '_mean_d', dom
   ! post inf sd copy
   write(restart_files_in(ens_size + 6, dom), '(A, A, i2.2, A)') trim(inflation_in(2)), '_sd_d', dom
enddo

! output extras
do dom = 1, num_domains
   ! Prior
   ! mean
   write(restart_files_out(ens_size + 1, dom, 1), '(A, i2.2, A)') 'PriorDiag_mean_d', dom, '.nc'
   ! sd
   write(restart_files_out(ens_size + 2, dom, 1), '(A, i2.2, A)') 'PriorDiag_sd_d', dom, '.nc'
   ! prior inf copy
   write(restart_files_out(ens_size + 3, dom, 1), '(A, i2.2, A)') 'PriorDiag_inf_mean_d', dom, '.nc'
   ! prior inf sd copy
   write(restart_files_out(ens_size + 4, dom, 1), '(A, i2.2, A)') 'PriorDiag_inf_sd_d', dom, '.nc'
   ! post inf copy - not used
   write(restart_files_out(ens_size + 5, dom, 1), '(A, A, i2.2, A)') trim(inflation_out(2)), '_mean_d', dom, '.nc'
   ! post inf sd copy - not used
   write(restart_files_out(ens_size + 6, dom, 1), '(A, A, i2.2, A)') trim(inflation_out(2)), '_sd_d', dom, '.nc'

   ! Posterior
   ! mean
   write(restart_files_out(ens_size + 1, dom, 2), '(A, i2.2, A)') 'mean_d', dom, '.nc'
   ! sd
   write(restart_files_out(ens_size + 2, dom, 2), '(A, i2.2, A)') 'sd_d', dom, '.nc'
   ! prior inf copy
   write(restart_files_out(ens_size + 3, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(1)), '_mean_d', dom
   ! prior inf sd copy
   write(restart_files_out(ens_size + 4, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(1)), '_sd_d', dom
   ! post inf copy
   write(restart_files_out(ens_size + 5, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(2)), '_mean_d', dom
   ! post inf sd copy
   write(restart_files_out(ens_size + 6, dom, 2), '(A, A, i2.2, A)') trim(inflation_out(2)), '_sd_d', dom


enddo

end subroutine io_filenames_init

!----------------------------------
!> Destroy module storage
subroutine end_io_filenames()

deallocate(restart_files_in, restart_files_out)

end subroutine end_io_filenames

!----------------------------------
!> @}
end module io_filenames_mod
