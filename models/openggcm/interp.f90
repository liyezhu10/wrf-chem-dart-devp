! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

!> development/test code for interpolation in the CTIM/IPE O+ (and other ions/e-) grid.
!> Jimmy Raeder, UNH, August 2017
!>
!> grid description:
!> 1.  a regular lat/lon grid of footpoints in GEO lat/lon grid
!> 2.  from each footpoint, a field line (dipole only) emerges, which has an 
!>     unevenly spaced number of points, however, the spacing of the points along 
!>     the field line is the same for all field lines hz(iz).
!> 3.  The dipolar field lines are more regular in SM (no azimuth component), but 
!>     that does not help much.
!> 4.  The grid is not completely irregular, but made of 'bricks', i.e., 4 field 
!>     lines making a rectangle at the base will make a brick with the 8 points 
!>     given at hz(iz) and hz(iz+1).  Such grids are common in FEM (Finite Element)
!>     analysis.  Such grids are called hexahedron (hex8) element bricks.
!>
!> possible simplifications
!> 1.  assume field lines are radial: may work ok over polar cap, but rapidly gets 
!>     worse at lower latitudes.  Specifically, at mid latitude the horizontal 
!>     displacement error becomes roughly equal the height (~600 km) --> not acceptable.
!> 2.  transform everything to SM.  Then the azimuthal coordinate is evenly spaced, 
!>     leaving us basically only with an irregular grid in latitude and height.
!>     Simplifies the problem somewhat, but we loose generality.  In the future 
!>     this may change, for example when better coordinates are used, i.e., IGRF/AGGCM
!>     coordinates.  Since there is still irregularity we may not save much 
!>     computationally --> dismissed.
!>
!> Considerations:
!> 1.  We are now having a hex8 grid.  Interpolation functions (speficically Lagrange
!>     interpolation) exist for hex8 grids.  These functions are defined for 
!>     isoparametric coordinates, i.e., the unit brick [-1,1]x[-1,1]x[-1,1].
!>     The mapping from the unit brick into x,y,z space is easy, but not linear, so 
!>     the inverse is not easy to come by.  This is fine for FEM but not for us, 
!>     because we need the inverse mapping (x,y,z)-->(xi,yi,zi), where the latter are
!>     the isoparametric coordinates. --> dismiss
!> 2.  One can instead divide the bricks into tetrahedra (tet4).  For tet4, the 
!>     isoparametric mapping is linear, so the inverse is easy.  The unit tet4 has 
!>     the 4 corners (0,0,0) (1,0,0) (0,1,0) and (0,0,1).  It is trivial to find out
!>     if a point lies within a tet4: map (x,y,z) --> (xi,yi,zi), then the condition
!>     is xi>0 && yi>0 && zi>0 && xi+yi+zi<1.  Interpolation is easily done with 
!>     the Lagrange interpolation functions N_j(xi,yi,zi), often called 
!>     'shape functions' in the FEM literature.
!> 3.  Search in tet4 grids is easy and fast, with some pre-computed values 
!>     (center, list of neighbors), i.e., descending distance, which scales as 
!>     ~(N)**(1/3) where N is the total number of tet4.  Also, if the current search 
!>     point is close to the previous (as in LOS integration), the search terminates
!>     quickly.
!> 4.  tet4 grids are quite universal, even for completely irregular grids, once the 
!>     decomposition is found.
!> 5.  We should do this in (x,y,z) space as opposed to (r,phi,the) spherical 
!>     coordinates since DA is all about distance.
!> 6.  The division of bricks into tetrahedra is not unique, see: 
!>     https://www.ics.uci.edu/~eppstein/projects/tetra/
!>     Division into 5 tet4 leads to non-conformity: the faces and edges of tet4s 
!>     do not line up between bricks.  The interpolation is then not continuous 
!>     across brick faces.  The division into 6 tet4 does not have that problem,
!>     so we use that.
!>
!> In this test code I minimize any efforts re/ memory allocation.  The grid will 
!> be read in a simple as possible (ascii).  The purpose is to demonstrate the 
!> algorithm and test accuracy.

program interp

use        types_mod, only : r4, r8

use    utilities_mod, only : register_module, initialize_utilities,      &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, logfileunit,      &
                             finalize_utilities 

use  openggcm_interp_mod, only : g_oplus_pre, g_oplus_int, g_oplus_int_pnt, nsearch

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

real(r8), allocatable ::  opl_grid_geo_lon(:,:,:) 
real(r8), allocatable ::  opl_grid_geo_lat(:,:,:) 
real(r8), allocatable ::  opl_grid_geo_hgt(:,:,:) 
real(r8), allocatable ::                 u(:,:,:) 

real(r8) :: a(3,3), b(3,3)

integer :: nphi ! number of longitudes
integer :: nthe ! number of latitudes
integer :: nnz  ! number of vertical layers

real(r8) :: pi, rad, re, x1, y1, z1, x2, y2, z2 
real(r8) :: s, x, y, z, rr
real(r8) :: a1, a2, a3, a4, output
real(r8) :: tmpz, tmpp, tmpt

integer :: nn
integer :: iz, jz, ilat, jlat, ilon, jlon, i
integer :: iunit, istat
integer :: good, num

! namelist with default values

character(len=256) :: data_file = '../data/grid-oplus'
integer :: test_case = 4
integer :: debug = 0

namelist / test_interp_nml / &
    data_file, test_case, debug

!----------------------------------------------------------------------
! start
!----------------------------------------------------------------------

! Read the DART namelist for this model

call find_namelist_in_file('input.nml', 'test_interp_nml', iunit)
read(iunit, nml = test_interp_nml, iostat = istat)
call check_namelist_read(iunit, istat, 'test_interp_nml')

!Read some data

open(10, file=data_file, status='old')
read(10,*)nphi, nthe, nnz

write(*,*)'nphi expected to be 20   is ', nphi
write(*,*)'nthe expected to be 91   is ', nthe
write(*,*)'nnz  expected to be 50   is ', nnz

allocate( opl_grid_geo_lon(nnz,nthe,nphi), &
          opl_grid_geo_lat(nnz,nthe,nphi), &
          opl_grid_geo_hgt(nnz,nthe,nphi), &
                         u(nnz,nthe,nphi)  )

do jlon=1, nphi  ! longitude
do jlat=1, nthe  ! latitude
do   jz=1, nnz   ! vertical

   read(10,*)iz,ilat,ilon,a1,a2,a3,a4

   opl_grid_geo_lon(iz,ilat,ilon) = a1  ! longitude, degrees
   opl_grid_geo_lat(iz,ilat,ilon) = a2  ! latitude, degrees
   opl_grid_geo_hgt(iz,ilat,ilon) = a3  ! height, meters
                  u(iz,ilat,ilon) = a4  ! values
enddo
enddo
enddo
close(10)

call g_oplus_pre(nphi, nthe, nnz, opl_grid_geo_lon, &
              opl_grid_geo_lat, opl_grid_geo_hgt, test_case)

!This block looks like these could all be parameters

pi  = 4.0_r8*atan(1.0_r8)
rad = pi/180.0_r8
re  = 6372.0e3

if (test_case == 4) then
      do jlon=1,nphi
      do jlat=1,nthe
      do   jz=1,nnz
         !u(jz,jlat,jlon)=opl_grid_geo_lon(jz,jlat,jlon)
         !u(jz,jlat,jlon)=opl_grid_geo_lat(jz,jlat,jlon)
         u(jz,jlat,jlon)=opl_grid_geo_hgt(jz,jlat,jlon)
      enddo
      enddo
      enddo
      num=1000
      good=0
      do i=1,num
        tmpz=rand(0)*554635.+100000.
        tmpp=rand(0)*2.*pi
        tmpt=rand(0)*pi
        x=(re+tmpz)*cos(tmpp)*sin(tmpt)
        y=(re+tmpz)*sin(tmpp)*sin(tmpt)
        z=(re+tmpz)*cos(tmpt)
        call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
        if (istat.eq.0) good=good+1
        write(*,*)'test4 ',x,y,z,sqrt(x**2+y**2+z**2)-re,output,istat
      enddo
      write(*,*) 'found ',good,' out of ',num 
      x=6522000.
      y=0.
      z=0.
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test4 ',x,y,z,sqrt(x**2+y**2+z**2),output,istat

elseif (test_case == 1) then
   !.... simple test in 1-spaced grid
   !     note that grid is overwritten with a simple box grid
   !     in g_oplus_pre(), look for test.eq.1
   !.... make up a grid with point values

      do jlon=1,nphi
      do jlat=1,nthe
      do   jz=1,nnz
         u(jz,jlat,jlon)=jz+jlat+jlon
      enddo
      enddo
      enddo
      x=11.11
      y=22.22
      z=33.33
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1 ',x,y,z,output,istat,x+y+z
      x=11.11
      y=22.22
      z=33.33
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1 ',x,y,z,output,istat,x+y+z
      x=12.11
      y=23.22
      z=34.33
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1 ',x,y,z,output,istat,x+y+z
      x=12.11
      y=13.22
      z=14.33
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1 ',x,y,z,output,istat,x+y+z
      x=12.00
      y=13.00
      z=14.00
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1 ',x,y,z,output,istat,x+y+z
      do i=1,400
      s=0.01*float(i)
      x=2.22 +1.54*s
      y=11.3 -2.54*s
      z=33.1 -1.74*s
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test1a ',x,y,z,output-(x+y+z),istat,nsearch
      enddo

elseif (test_case == 2) then
   !.... now with the real data, do some profiles, just as one would do 
   !     for LOS integration. looking up, should give an e- profile, 
   !.... typical values should be ~10^11 peaking at 250-400 km

   x1 = (re +  70000.00_r8) * cos(111.55_r8 * rad) * cos(70.33_r8 * rad)
   y1 = (re +  70000.00_r8) * sin(111.55_r8 * rad) * cos(70.33_r8 * rad)
   z1 = (re +  70000.00_r8) * sin( 70.33_r8 * rad)
   x2 = (re + 800000.00_r8) * cos(114.55_r8 * rad) * cos(71.33_r8 * rad)
   y2 = (re + 800000.00_r8) * sin(114.55_r8 * rad) * cos(71.33_r8 * rad)
   z2 = (re + 800000.00_r8) * sin( 71.33_r8 * rad)
     
   nn=400 ! Why 400?
   do i=1,nn
      s  = float(i-1)/float(nn-1)
      x  = x1+s*(x2-x1)
      y  = y1+s*(y2-y1)
      z  = z1+s*(z2-z1)
      rr = (sqrt(x**2+y**2+z**2)-re)/1000.0 
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test2 ',x,y,z,output,istat,nsearch
   enddo

else

   ! ... now with the real data,slant LOS like COSMIC
   ! ... looking up, should give an e- profile, typical 
   ! ... values should be ~10^11 peaking at 250-400 km

   x1 = (re + 3000.0e3) * cos( 11.55_r8 * rad) * cos(80.33_r8 * rad)
   y1 = (re + 3000.0e3) * sin( 11.55_r8 * rad) * cos(80.33_r8 * rad)
   z1 = (re + 3000.0e3) * sin( 80.33_r8 * rad)
   x2 = (re + 3000.0e3) * cos(144.55_r8 * rad) * cos( 6.33_r8 * rad)
   y2 = (re + 3000.0e3) * sin(144.55_r8 * rad) * cos( 6.33_r8 * rad)
   z2 = (re + 3000.0e3) * sin(  6.33_r8 * rad)
     
   nn=400 ! Why 400?
   do i=1,nn
      s  = float(i-1)/float(nn-1)
      x  = x1+s*(x2-x1)
      y  = y1+s*(y2-y1)
      z  = z1+s*(z2-z1)
      rr = (sqrt(x**2+y**2+z**2)-re)/1000.0 
      call g_oplus_int(nphi,nthe,nnz,u,x,y,z,output,istat)
      write(*,*)'test3 ',rr,output,istat,nsearch
   enddo

endif

deallocate(opl_grid_geo_lon, opl_grid_geo_lat, opl_grid_geo_hgt, u)

stop

end program interp

