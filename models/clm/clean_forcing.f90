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
!              a2x6h_Faxa_rainc     rainc
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
!-------------------------------------------------------------------------------

use            types_mod, only : r8,i8

use             sort_mod, only : index_sort

use        utilities_mod, only : initialize_utilities, &
                                 finalize_utilities, &
                                 find_namelist_in_file, &
                                 check_namelist_read, &
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
character(len=*), parameter :: directory = '/glade/p_old/image/thoar/CAM_DATM/4xdaily'

character(len=16) :: varname
character(len=16) :: variables(10) = (/'a2x6h_Faxa_lwdn ', &
                                       'a2x6h_Faxa_swndr', &
                                       'a2x6h_Faxa_swvdr', &
                                       'a2x6h_Faxa_swndf', &
                                       'a2x6h_Faxa_swvdf', &
                                       'null            ', &
                                       'a2x6h_Faxa_rainc', &
                                       'a2x6h_Faxa_rainl', &
                                       'a2x6h_Faxa_snowc', &
                                       'a2x6h_Faxa_snowl'/)

integer :: indices(ensemble_size)
integer ::    ncid(ensemble_size)
character(len=256) :: input_file(ensemble_size)

real(r8), allocatable :: tensor(:,:,:) ! nx,ny,ensemble_size
real(r8) :: minvalue, maxvalue, q1, q2, q3, iqr
real(r8) :: newhigh, newlow, original_variance, new_variance
logical :: suspicious, newmin, newmax

integer :: imember, itime, iy, ix, nT, ny, nx
integer :: ncstart(3)
integer :: nccount(3)
integer :: numdims, dimlens(NF90_MAX_VAR_DIMS)
integer :: varid

integer            :: iunit, io
character(len=512) :: string1, string2, string3

integer(i8) :: icount

!-------------------------------------------------------------------------------
! namelist
!-------------------------------------------------------------------------------

integer :: year = 2008
real(r8) :: iqr_multiplier = 20.0_r8

namelist /clean_forcing_nml/ year, iqr_multiplier

!======================================================================

! Read the namelist entry
call find_namelist_in_file('input.nml', 'clean_forcing_nml', iunit)
read(iunit, nml = clean_forcing_nml, iostat = io)
call check_namelist_read(iunit, io, 'clean_forcing_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=clean_forcing_nml)
if (do_nml_term()) write(     *     , nml=clean_forcing_nml)

call initialize_utilities(progname=routine)

100 format('/glade/p_old/image/thoar/CAM_DATM/4xdaily/CAM_DATM.cpl_',i4.4,'.ha2x1dx6h.',i4.4,'.nc')

do imember = 1,ensemble_size
   write(input_file(imember),100) imember, year

   ncid(imember) = nc_open_file_readwrite(input_file(imember),routine)
enddo

varname = variables(1)

call nc_get_variable_num_dimensions(ncid(1), varname, numdims)
call nc_get_variable_size(ncid(1), varname, dimlens)

if ( .false. ) write(*,*)'dimlens are ',dimlens(1:numdims)
nx = dimlens(1)
ny = dimlens(2)
nT = dimlens(3)

allocate(tensor(ensemble_size,nx,ny))

! TIMESTEP : do itime = 1,nT
TIMESTEP : do itime = 1169,1169

   call fill_tensor(itime)

   icount = 0_i8

   LATITUDE : do iy=1,ny
   LONGITUDE : do ix=1,nx
!  LATITUDE : do iy=33,33
!  LONGITUDE : do ix=1,nx

      q2 = sum(tensor(:,ix,iy))/ensemble_size
      original_variance = sum((tensor(:,ix,iy) - q2)**2)/(ensemble_size -1)

      suspicious = .true.

      CLEAN : do while ( suspicious )

         newmin = .false.
         newmax = .false.

         call index_sort(tensor(:,ix,iy), indices, ensemble_size)

!        q1  = tensor(indices(20),ix,iy)
         q2  = tensor(indices(40),ix,iy) * 0.5_r8 + &
               tensor(indices(41),ix,iy) * 0.5_r8
!        q3  = tensor(indices(60),ix,iy)

!        iqr = q3 - q1

!        minvalue = tensor(indices( 1),ix,iy)
!        maxvalue = tensor(indices(80),ix,iy)

!        newlow   = tensor(indices(10),ix,iy) * 0.5_r8 + &
!                   tensor(indices(11),ix,iy) * 0.5_r8
!        newhigh  = tensor(indices(70),ix,iy) * 0.5_r8 + &
!                   tensor(indices(71),ix,iy) * 0.5_r8

!        if (minvalue < (q1 - iqr_multiplier*iqr)) then
!           write(*,*)'low ',itime, ix, iy, iqr, minvalue, newlow
!           newmin = .true.
!        endif

!        if (maxvalue > (q3 + iqr_multiplier*iqr)) then
!           write(*,*)'high ',itime, ix, iy, iqr, maxvalue, newhigh
!           newmax = .true.
!        endif

!        if ( newmin .or. newmax ) then
!           suspicious = .true.
!           icount = icount + 1_i8
!        else
!           suspicious = .false.
!        endif

!        if (suspicious) then
!           do imember = 1,ensemble_size
!              write(*,*)        imember,  tensor(        imember ,ix,iy), &
!                        indices(imember), tensor(indices(imember),ix,iy)
!           enddo
!        endif

!        if (newmin) tensor(indices( 1),ix,iy) = newlow
!        if (newmax) tensor(indices(80),ix,iy) = newhigh

         ! replace high value with median and recalculate variance

         tensor(indices(80),ix,iy) = q2

         q2 = sum(tensor(:,ix,iy))/ensemble_size
         new_variance = sum((tensor(:,ix,iy) - q2)**2)/(ensemble_size -1)

         if (original_variance/new_variance > 2.0) then
            write(*,*)itime,iy,ix,original_variance/new_variance
            suspicious = .true.
         else
            suspicious = .false.
         endif

      enddo CLEAN

   enddo LONGITUDE
   enddo LATITUDE

   write(*,*)'timestep ',itime,' had ',icount, &
             ' values replaced when iqr_multiplier is ',iqr_multiplier

enddo TIMESTEP


do imember = 1,ensemble_size
   call nc_close_file(ncid(imember),routine)
enddo

call finalize_utilities(routine)

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine fill_tensor(itime)
integer, intent(in) :: itime

character(len=*), parameter :: routine = 'fill_tensor'
integer :: imember, io

MEMBER : do imember = 1,ensemble_size

   ncstart = (/ 1,  1, itime/)
   nccount = (/nx, ny,     1/)

   io = nf90_inq_varid(ncid(imember), varname, varid)
   call nc_check(io, routine, 'inq_var_id '//trim(varname))

   io = nf90_get_var(ncid(imember), varid, tensor(imember,:,:), &
                      start=ncstart, count=nccount)
   call nc_check(io, routine, 'get_var '//trim(varname))

!     minvalue = minval(tensor(imember,:,:))
!     maxvalue = maxval(tensor(imember,:,:))
!
!     write(*,*)'timestep ',itime,' member ',imember, &
!               ' min,max ', minvalue, maxvalue

enddo MEMBER

end subroutine fill_tensor


end program clean_forcing

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
