! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program clean_forcing

!-------------------------------------------------------------------------------
! /glade/p_old/image/thoar/CAM_DATM/4xdaily
!             CAM_DATM.cpl_0055.ha2x1dx6h.2008.nc
!             CAM_DATM.cpl_0055.ha2x1dx6h.2009.nc
!             CAM_DATM.cpl_0055.ha2x1dx6h.2010.nc
!               a2x6h_Faxa_rainc     rainc
!               a2x6h_Faxa_rainl     rainl
!               a2x6h_Faxa_snowc     snowc
!               a2x6h_Faxa_snowl     snowl
!               a2x6h_Faxa_lwdn      lwdn
!               a2x6h_Faxa_swndr     swndr
!               a2x6h_Faxa_swvdr     swvdr
!               a2x6h_Faxa_swndf     swndf
!               a2x6h_Faxa_swvdf     swvdf
!
!   netcdf CAM_DATM.cpl_0055.ha2x1dx6h.2008 {
!   dimensions:
!      time = UNLIMITED ; // (1460 currently)
!      doma_nx = 144 ;
!      doma_ny = 96 ;
!      a2x6h_nx = 144 ;
!      a2x6h_ny = 96 ;
!      ntb = 2 ;
!
!   float a2x6h_Faxa_lwdn(time, a2x6h_ny, a2x6h_nx) ;
!           a2x6h_Faxa_lwdn:_FillValue = 1.e+30 ;
!           a2x6h_Faxa_lwdn:units = "W m-2" ;
!           a2x6h_Faxa_lwdn:long_name = "Downward longwave heat flux" ;
!           a2x6h_Faxa_lwdn:standard_name = "downwelling_longwave_flux" ;
!           a2x6h_Faxa_lwdn:internal_dname = "a2x6h" ;
!           a2x6h_Faxa_lwdn:cell_methods = "time: mean" ;
!
! swvdf timestep  352 max variance 25982340.000 144  76 max/median    24584.617 144  76
! member 64 is 45594.96 ... others are x.y
!
!-------------------------------------------------------------------------------

use            types_mod, only : r8, i8, MISSING_R8

use             sort_mod, only : index_sort

use       random_seq_mod, only : random_seq_type, &
                                 init_random_seq, &
                                 random_gaussian

use        utilities_mod, only : initialize_utilities, &
                                 finalize_utilities, &
                                 find_namelist_in_file, &
                                 check_namelist_read, &
                                 open_file, close_file, &
                                 nmlfileunit, &
                                 do_nml_file, &
                                 do_nml_term, &
                                 error_handler, E_ERR, E_MSG

use netcdf_utilities_mod, only : nc_check,                       &
                                 nc_open_file_readwrite,         &
                                 nc_close_file,                  &
                                 nc_synchronize_file,            &
                                 nc_begin_define_mode,           &
                                 nc_end_define_mode,             &
                                 nc_get_variable_num_dimensions, &
                                 nc_get_variable_size,           &
                                 nc_get_variable,                &
                                 nc_put_variable

use netcdf

implicit none

!-------------------------------------------------------------------------------
! version controlled file description for error handling, do not edit
!-------------------------------------------------------------------------------

character(len=*), parameter :: source   = '$URL$'
character(len=*), parameter :: revision = '$Revision$'
character(len=*), parameter :: revdate  = '$Date$'
character(len=*), parameter :: routine  = 'clean_forcing'

integer, parameter :: ensemble_size = 80
character(len=*), parameter :: directory = '/glade/p/cisl/dares/thoar/CAM_DATM/4xdaily'

character(len=16) :: varname
character(len=16) :: variables(10) = (/'a2x6h_Faxa_lwdn ', &
                                       'a2x6h_Faxa_swndf', &
                                       'a2x6h_Faxa_swvdf', &
                                       'a2x6h_Faxa_swndr', &
                                       'a2x6h_Faxa_swvdr', &
                                       'null            ', &
                                       'a2x6h_Faxa_rainc', &
                                       'a2x6h_Faxa_rainl', &
                                       'a2x6h_Faxa_snowc', &
                                       'a2x6h_Faxa_snowl'/)

integer  :: indices(ensemble_size)
integer  ::    ncid(ensemble_size)
real(r8) ::  sorted(ensemble_size)
character(len=256) :: input_file(ensemble_size)

real(r8), allocatable :: tensor(:,:,:) ! nx,ny,ensemble_size

type(random_seq_type) :: r

integer :: imember, nT, ny, nx
integer :: ncstart(3)
integer :: nccount(3)
integer :: numdims, dimlens(NF90_MAX_VAR_DIMS)
integer :: ivar, varid

integer            :: iunit, io
character(len=512) :: string1, string2, string3

!-------------------------------------------------------------------------------
! namelist
!-------------------------------------------------------------------------------

integer  :: year = 2008
real(r8) :: criterion = 100.0_r8
real(r8) :: lat = MISSING_R8
real(r8) :: lon = MISSING_R8
integer  :: varnum = 1

namelist /clean_forcing_nml/ year, criterion, lat, lon, varnum

!===============================================================================

! Read the namelist entry
call find_namelist_in_file('input.nml', 'clean_forcing_nml', iunit)
read(iunit, nml = clean_forcing_nml, iostat = io)
call check_namelist_read(iunit, io, 'clean_forcing_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=clean_forcing_nml)
if (do_nml_term()) write(     *     , nml=clean_forcing_nml)

call initialize_utilities(progname=routine)

call init_random_seq(r,1)

100 format(A,'/CAM_DATM.cpl_',i4.4,'.ha2x1dx6h.',i4.4,'.nc')

do imember = 1,ensemble_size
   write(input_file(imember),100) directory, imember, year
   ncid(imember) = nc_open_file_readwrite(input_file(imember),routine)
   write(*,*)'Opened '//trim(input_file(imember))
enddo

do ivar = 1,5

   varname = variables(ivar)

   call nc_get_variable_num_dimensions(ncid(1), varname, numdims)
   call nc_get_variable_size(ncid(1), varname, dimlens)

   if ( .true. ) write(*,*)'dimlens are ',dimlens(1:numdims)
   nx = dimlens(1)
   ny = dimlens(2)
   nT = dimlens(3)
   
   allocate(tensor(ensemble_size,nx,ny))
   
   write(string1,'(A,''_'',i4,''.txt'')') trim(varname), year
   iunit = open_file(string1)

   call purify()

   deallocate(tensor)

enddo

call close_file(iunit)

do imember = 1,ensemble_size
   call nc_close_file(ncid(imember),routine)
   write(*,*)'Closed '//trim(input_file(imember))
enddo

call finalize_utilities(routine)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine read_tensor(itime)
integer, intent(in) :: itime

character(len=*), parameter :: routine = 'read_tensor'
integer :: imember, io

MEMBER : do imember = 1,ensemble_size

   ncstart = (/ 1,  1, itime/)
   nccount = (/nx, ny,     1/)

   io = nf90_inq_varid(ncid(imember), varname, varid)
   call nc_check(io, routine, 'inq_var_id '//trim(varname))

   io = nf90_get_var(ncid(imember), varid, tensor(imember,:,:), &
                      start=ncstart, count=nccount)
   call nc_check(io, routine, 'get_var '//trim(varname))

enddo MEMBER

end subroutine read_tensor


subroutine write_tensor(itime)
integer, intent(in) :: itime

character(len=*), parameter :: routine = 'write_tensor'
integer :: imember, io

MEMBER : do imember = 1,ensemble_size

   ncstart = (/ 1,  1, itime/)
   nccount = (/nx, ny,     1/)

   io = nf90_inq_varid(ncid(imember), varname, varid)
   call nc_check(io, routine, 'inq_var_id '//trim(varname))

   io = nf90_put_var(ncid(imember), varid, tensor(imember,:,:), &
                      start=ncstart, count=nccount)
   call nc_check(io, routine, 'put_var '//trim(varname))

enddo MEMBER

end subroutine write_tensor

!-------------------------------------------------------------------------------
!> find and replace the outliers with a noisy estimate of the median

subroutine purify()

! swvdf timestep  352 max variance 25982340.000 144  76 max/median    24584.617 144  76

! The noise distribution is estimated from the center of the ensemble.
! The correct way is to estimate the distribution and 
! replace the outlier with the expected value for that quantile.
! The correct way is overkill. The really correct way may
! borrow strength from the surrounding gridcells.
! Another way may be to simply design and apply a filter.

integer  :: itime, iy, ix, imember
real(r8) :: mean, variance, new_variance
real(r8) :: q1, q2, q3, iqr, noise
logical  :: suspicious
integer(i8) :: num_outliers

TIMESTEP : do itime = 1,nT

   call read_tensor(itime)

   num_outliers = 0_i8

   LATITUDE : do iy=1,ny
   LONGITUDE : do ix=1,nx

!           do imember = 1,ensemble_size
!              write(iunit,*) ix,iy,imember, tensor(imember,ix,iy)
!           enddo

      mean     = sum( tensor(:,ix,iy))/ensemble_size
      variance = sum((tensor(:,ix,iy) - mean)**2)/(ensemble_size-1)

      ! some gridcells are all identical values
      if (variance < tiny(variance)) cycle LONGITUDE

      suspicious = .true. ! presume everyone is guilty

      CLEAN : do while ( suspicious )

         ! calculate quantiles, sort into ascending order

         call index_sort(tensor(:,ix,iy), indices, ensemble_size)
         sorted = tensor(indices,ix,iy)

         q1  =  sorted(20)
         q2  = (sorted(40) + sorted(41)) / 2.0_r8
         q3  =  sorted(60)
         iqr = q3 - q1

         ! See if this cell has an outlier or not by
         ! replacing maximum with new value and recalculate variance

         sorted(ensemble_size) = q2
         mean         = sum(sorted)/ensemble_size
         new_variance = sum((sorted - mean)**2)/(ensemble_size-1)

!        if (new_variance < tiny(new_variance)) then
!           write(*,*)itime,ix,iy,tensor(:,ix,iy)
!           stop
!        endif

!        want to test "if (variance/new_variance > criterion) then" but
!        new variance could be tiny
         if (variance > new_variance * criterion) then
            write(*,*)'suspicious ',itime,ix,iy,variance,new_variance*criterion
            suspicious = .true.
            variance   = new_variance
         else
!           write(*,*)'clean ',itime,ix,iy,variance/new_variance
            suspicious = .false.
         endif

         if (suspicious) then ! log bad column to unique file
            do imember = 1,ensemble_size
               write(iunit,*)        imember,  tensor(        imember ,ix,iy), &
                         indices(imember), tensor(indices(imember),ix,iy)
            enddo
         else
            exit CLEAN
         endif

         ! If we get this far, the maximum is an outlier.

         num_outliers = num_outliers + 1_i8

!        mean   = sum(sorted(10:70))/61.0_r8
!        stddev = sqrt(sum(sorted(10:70) - mean)**2)/60.0_r8
!        write(*,*)'ix,iy,stddev,iqr = ',ix,iy,stddev,iqr

!        noise = random_gaussian(r, 0.0_r8, 3.0_r8*stddev)
         noise = random_gaussian(r, 0.0_r8, iqr)

         tensor(indices(ensemble_size),ix,iy) = q2 + abs(noise) 

         ! log candidate column to same unique file
         do imember = 1,ensemble_size
            write(iunit,*)        imember,  tensor(        imember ,ix,iy), &
                      indices(imember), tensor(indices(imember),ix,iy)
         enddo

      enddo CLEAN

   enddo LONGITUDE
   enddo LATITUDE

   write(*,*)trim(varname)//' timestep ',itime,' had ',num_outliers, &
             ' values replaced when criterion is ',criterion
   if (num_outliers > 0_i8) call write_tensor(itime)

enddo TIMESTEP

end subroutine purify


end program clean_forcing

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
