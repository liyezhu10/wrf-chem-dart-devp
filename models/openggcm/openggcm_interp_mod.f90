! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

module openggcm_interp_mod

use        types_mod, only : r4, r8, i8

use    utilities_mod, only : register_module, initialize_utilities,      &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit,   &
                             do_nml_file, do_nml_term, logfileunit,      &
                             finalize_utilities 

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use   state_structure_mod, only : get_dart_vector_index

implicit none
private

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

integer, parameter :: NTET_MAX = 600000  !... max num of tets
integer, parameter :: NPOI_MAX = 128000  !... max num of points
! these are not known a-priory, but could be adjusted later
integer, parameter :: MTET = 80  !... max num of neighbors for each tet
integer, parameter :: MPOI = 80  !... max num of tets for each point
integer, parameter :: MAXDIM = 100  !... max of each dimension

!..... x,y,z grid point coord
real(r8) :: hx(NPOI_MAX),hy(NPOI_MAX),hz(NPOI_MAX) 

!..... tet mapping matrix (location of center and first corner, 
real(r8) :: tcen(3,NTET_MAX),tlow(3,NTET_MAX),tmap(3,3,NTET_MAX)

!..... 4 corners + neicount + neighbors + myself
integer :: itet(MTET+6,NTET_MAX) 
!..... for each point, count, and which tet4 it belongs to
integer :: jtet(MPOI+1,NPOI_MAX) 
!..... remember the ip,it,iz indices
integer :: kpoi(3,NPOI_MAX)
!..... number of tet4 tetrahedra, search steps
integer :: ntet, nsearch

!..... corner id for grid indices
integer :: ii(MAXDIM,MAXDIM,MAXDIM)

real(r8), parameter :: PI  = 4.0_r8*atan(1.0_r8)
real(r8), parameter :: RAD = PI/180.0_r8
real(r8), parameter :: RE  = 6372.0e3

! Logical to keep track of if we have initialized g_oplus_int
logical, save :: module_initialized = .false.

public :: g_oplus_pre, &
          g_oplus_int, &
          g_oplus_int_pnt, &
          convert_to_cartesian

contains

!-----------------------------------------------------------------------
! Put the public routines first.
! The private routines come after.
!-----------------------------------------------------------------------

!> initialize routine - only do this once

subroutine g_oplus_pre(np, nt, nz, gp, gt, gz, test)

integer,           intent(in)  :: np, nt, nz
real(r8),          intent(out) :: gp(nz,nt,np)
real(r8),          intent(out) :: gt(nz,nt,np)
real(r8),          intent(out) :: gz(nz,nt,np)
integer, optional, intent(in)  :: test

character(len=*), parameter :: routine = 'g_oplus_pre'

real(r8) :: a(3,3), b(3,3), det
real(r8) :: r, p, t

!!!!integer :: ii(np,nt,nz) 

integer :: ip, ip1, it, it1, iz, iz1, k, mytest
integer :: i1, i2, i3, i4, i5, i6, i7, i8
integer :: j, jt, kc, kt, l, ll, maxnei, ni, nj

if ( module_initialized ) then
   return ! only need to do this once.
else
   module_initialized = .true.
endif

if (present(test)) then
   mytest = test
else
   mytest = 1
endif

k = 0

!..... make grid in Cartesian coordinates, store in linear array
do ip = 1,np
do it = 1,nt
do iz = 1,nz

   k = k+1
   r = RE  + gz(iz,it,ip)
   p = RAD * gp(iz,it,ip)
   t = RAD * gt(iz,it,ip)

   !.... making a pole hole to avoid collapsed corners
   if(it ==  1) t = RAD*(-89.99_r8)
   if(it == nt) t = RAD*( 89.99_r8) 

   if (mytest == 1) then
      !..... test with simple unit grid
      hx(k)=ip
      hy(k)=it
      hz(k)=iz
   else  ! This is the conversion to cartesian coords
      hx(k) = r*cos(p)*cos(t)
      hy(k) = r*sin(p)*cos(t)
      hz(k) = r*sin(t)
   endif

   !... save to map grid values later
   kpoi(1,k) = iz
   kpoi(2,k) = it
   kpoi(3,k) = ip 
   ii(ip,it,iz) = k
enddo
enddo
enddo
write(0,*)'test ',hx(1),hy(1),hz(1),hx(k),hy(k),hz(k)

6000  format(3i8,3(1x,f12.5))

!...... create tet4 list
write(0,*)'max points ',k
ntet = 0
JTET = 0
ITET = 0 

PHILOOP:   do ip = 1,np
THETALOOP: do it = 1,nt-1
ZLOOP:     do iz = 1,nz-1

     !..... we do the 6-division here, 5-division is also possible
     !      local brick numbering for iz:  
     ! (ip , it) -> 1  
     ! (ip1, it) -> 2 
     ! (ip1,it1) -> 3 
     ! (ip, it1) -> 4  for iz1 add 4
     ! follow numbering in 'The Finite Element Method in Charged 
     ! Particle Optics', Anjam Khursheed, page 116
     !      NOOOOOO, it's wrong, see below for local numbering

      ip1=ip+1
      if(ip1 > np) ip1=1
      it1=it+1
      iz1=iz+1    !..... periodic in ip

      i1=ii(ip,it,iz)
      i2=ii(ip1,it,iz)
      i3=ii(ip1,it1,iz)
      i4=ii(ip,it1,iz)   !.... global numbers of 8 brick corners
      i5=ii(ip,it,iz1)
      i6=ii(ip1,it,iz1)
      i7=ii(ip1,it1,iz1)
      i8=ii(ip,it1,iz1)

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1 !... global indices of the 4 corners
      itet(2,ntet)=i4
      itet(3,ntet)=i8
      itet(4,ntet)=i7 
      !..... tet4 centers
      tcen(1,ntet)=0.25*(hx(i1)+hx(i4)+hx(i8)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i4)+hy(i8)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i4)+hz(i8)+hz(i7)) 
      do l=1,jtet(1,i1)
         if(jtet(l+1,i1).eq.ntet) goto 94971 
      enddo 

      !........ decompose brick8 into 6 tet4, see https://www.ics.uci.edu/~eppstein/projects/tetra/
      !..... make a list of tets that belong to points P1, ...
      !      those tets are neighbors
      !      needed later for the neighbor list
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)

      if (l .gt. NPOI_MAX)  &
         call error_handler(E_ERR, routine, 'error on NPOI', source, revision, revdate)

      jtet(l+1,i1)=ntet 
94971 continue
      do l=1,jtet(1,i4)
      if(jtet(l+1,i4).eq.ntet) goto 94961 
      enddo 
      jtet(1,i4)=jtet(1,i4)+1
      l=jtet(1,i4)
      if (l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i4)=ntet 
94961 continue
      do l=1,jtet(1,i8)
      if(jtet(l+1,i8).eq.ntet) goto 94951 
      enddo 
      jtet(1,i8)=jtet(1,i8)+1
      l=jtet(1,i8)
      if(l.gt.NPOI_MAX) & 
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i8)=ntet 
94951 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94941 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if (l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94941 continue
      !..... now compute mapping matrix isoparametric -->  physical
      a(1,1)=hx(i4)-hx(i1)
      a(1,2)=hx(i8)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i4)-hy(i1)
      a(2,2)=hy(i8)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i4)-hz(i1)
      a(3,2)=hz(i8)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      !..... invert for mapping physical --> isoparameteric
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      !..... store matrix
      tmap(:,:,ntet)=B 
      !..... first corner
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1
      itet(2,ntet)=i5
      itet(3,ntet)=i8
      itet(4,ntet)=i7 
      tcen(1,ntet)=0.25*(hx(i1)+hx(i5)+hx(i8)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i5)+hy(i8)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i5)+hz(i8)+hz(i7)) 
      do l=1,jtet(1,i1)
         if(jtet(l+1,i1).eq.ntet) goto 94911 
      enddo 
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)

      if (l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)

      jtet(l+1,i1)=ntet !..... add me to the list
94911 continue
      do l=1,jtet(1,i5)
      if(jtet(l+1,i5).eq.ntet) goto 94901 
      enddo 
      jtet(1,i5)=jtet(1,i5)+1
      l=jtet(1,i5)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i5)=ntet 
94901 continue
      do l=1,jtet(1,i8)
      if(jtet(l+1,i8).eq.ntet) goto 94891 
      enddo 
      jtet(1,i8)=jtet(1,i8)+1
      l=jtet(1,i8)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i8)=ntet 
94891 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94881 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94881 continue
      a(1,1)=hx(i5)-hx(i1)
      a(1,2)=hx(i8)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i5)-hy(i1)
      a(2,2)=hy(i8)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i5)-hz(i1)
      a(3,2)=hz(i8)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      tmap(:,:,ntet)=B
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1
      itet(2,ntet)=i5
      itet(3,ntet)=i6
      itet(4,ntet)=i7 
      tcen(1,ntet)=0.25*(hx(i1)+hx(i5)+hx(i6)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i5)+hy(i6)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i5)+hz(i6)+hz(i7)) 
      do l=1,jtet(1,i1)
      if(jtet(l+1,i1).eq.ntet) goto 94851 
      enddo 
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i1)=ntet 
94851 continue
      do l=1,jtet(1,i5)
      if(jtet(l+1,i5).eq.ntet) goto 94841 
      enddo 
      jtet(1,i5)=jtet(1,i5)+1
      l=jtet(1,i5)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i5)=ntet 
94841 continue
      do l=1,jtet(1,i6)
      if(jtet(l+1,i6).eq.ntet) goto 94831 
      enddo 
      jtet(1,i6)=jtet(1,i6)+1
      l=jtet(1,i6)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i6)=ntet 
94831 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94821 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94821 continue
      a(1,1)=hx(i5)-hx(i1)
      a(1,2)=hx(i6)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i5)-hy(i1)
      a(2,2)=hy(i6)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i5)-hz(i1)
      a(3,2)=hz(i6)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      tmap(:,:,ntet)=B 
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1
      itet(2,ntet)=i4
      itet(3,ntet)=i3
      itet(4,ntet)=i7 
      tcen(1,ntet)=0.25*(hx(i1)+hx(i4)+hx(i3)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i4)+hy(i3)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i4)+hz(i3)+hz(i7)) 
      do l=1,jtet(1,i1)
      if(jtet(l+1,i1).eq.ntet) goto 94791 
      enddo 
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i1)=ntet 
94791 continue
      do l=1,jtet(1,i4)
      if(jtet(l+1,i4).eq.ntet) goto 94781 
      enddo 
      jtet(1,i4)=jtet(1,i4)+1
      l=jtet(1,i4)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i4)=ntet 
94781 continue
      do l=1,jtet(1,i3)
      if(jtet(l+1,i3).eq.ntet) goto 94771 
      enddo 
      jtet(1,i3)=jtet(1,i3)+1
      l=jtet(1,i3)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i3)=ntet 
94771 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94761 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94761 continue
      a(1,1)=hx(i4)-hx(i1)
      a(1,2)=hx(i3)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i4)-hy(i1)
      a(2,2)=hy(i3)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i4)-hz(i1)
      a(3,2)=hz(i3)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      tmap(:,:,ntet)=B 
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1
      itet(2,ntet)=i2
      itet(3,ntet)=i3
      itet(4,ntet)=i7 
      tcen(1,ntet)=0.25*(hx(i1)+hx(i2)+hx(i3)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i2)+hy(i3)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i2)+hz(i3)+hz(i7)) 
      do l=1,jtet(1,i1)
      if(jtet(l+1,i1).eq.ntet) goto 94731 
      enddo 
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i1)=ntet 
94731 continue
      do l=1,jtet(1,i2)
      if(jtet(l+1,i2).eq.ntet) goto 94721 
      enddo 
      jtet(1,i2)=jtet(1,i2)+1
      l=jtet(1,i2)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i2)=ntet 
94721 continue
      do l=1,jtet(1,i3)
      if(jtet(l+1,i3).eq.ntet) goto 94711 
      enddo 
      jtet(1,i3)=jtet(1,i3)+1
      l=jtet(1,i3)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i3)=ntet 
94711 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94701 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94701 continue
      a(1,1)=hx(i2)-hx(i1)
      a(1,2)=hx(i3)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i2)-hy(i1)
      a(2,2)=hy(i3)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i2)-hz(i1)
      a(3,2)=hz(i3)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      tmap(:,:,ntet)=B 
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 

      !..... create a new tet, store global indices, compute mapping matrix, etc.
      ntet=ntet+1
      itet(1,ntet)=i1
      itet(2,ntet)=i2
      itet(3,ntet)=i6
      itet(4,ntet)=i7 
      tcen(1,ntet)=0.25*(hx(i1)+hx(i2)+hx(i6)+hx(i7)) 
      tcen(2,ntet)=0.25*(hy(i1)+hy(i2)+hy(i6)+hy(i7)) 
      tcen(3,ntet)=0.25*(hz(i1)+hz(i2)+hz(i6)+hz(i7)) 
      do l=1,jtet(1,i1)
      if(jtet(l+1,i1).eq.ntet) goto 94671 
      enddo 
      jtet(1,i1)=jtet(1,i1)+1
      l=jtet(1,i1)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i1)=ntet 
94671 continue
      do l=1,jtet(1,i2)
      if(jtet(l+1,i2).eq.ntet) goto 94661 
      enddo 
      jtet(1,i2)=jtet(1,i2)+1
      l=jtet(1,i2)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i2)=ntet 
94661 continue
      do l=1,jtet(1,i6)
      if(jtet(l+1,i6).eq.ntet) goto 94651 
      enddo 
      jtet(1,i6)=jtet(1,i6)+1
      l=jtet(1,i6)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i6)=ntet 
94651 continue
      do l=1,jtet(1,i7)
      if(jtet(l+1,i7).eq.ntet) goto 94641 
      enddo 
      jtet(1,i7)=jtet(1,i7)+1
      l=jtet(1,i7)
      if(l.gt.NPOI_MAX) &
         call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
      jtet(l+1,i7)=ntet 
94641 continue
      a(1,1)=hx(i2)-hx(i1)
      a(1,2)=hx(i6)-hx(i1)
      a(1,3)=hx(i7)-hx(i1)
      a(2,1)=hy(i2)-hy(i1)
      a(2,2)=hy(i6)-hy(i1)
      a(2,3)=hy(i7)-hy(i1)
      a(3,1)=hz(i2)-hz(i1)
      a(3,2)=hz(i6)-hz(i1)
      a(3,3)=hz(i7)-hz(i1)
      b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      b(1,2) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
      b(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
      det = b(1,1)*a(1,1) + b(1,2)*a(2,1) + b(1,3)*a(3,1)
      b(1,1) = b(1,1)/det
      b(1,2) = b(1,2)/det
      b(1,3) = b(1,3)/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      tmap(:,:,ntet)=B 
      tlow(1,ntet)=hx(i1)
      tlow(2,ntet)=hy(i1)
      tlow(3,ntet)=hz(i1) 
enddo ZLOOP
enddo THETALOOP
enddo PHILOOP

write(0,*)'g_oplus_pre generated tet4 mesh, ntet= ',ntet

!...... for each tet4, find its neighbors.  Need this for fast search.
!       basically invert jtet: the tets belonging to each of one tet's corner points are neighbors.
!       however, there are dups and triplets if they share an edge or face.
maxnei=0
do kt=1,ntet                !..... loop tets
   do ll=1,4                !..... loop corners
      kc=itet(ll,kt)        !..... the global number of that corner
  !if (abs(hx(kc)).lt.1.or.abs(hy(kc)).lt.1.or.abs(hz(kc)).lt.1) then
  !  write(0,*) kt,hx(kc),hy(kc),hz(kc)
  !endif
      nj=jtet(1,kc)         !..... number of tets belonging to that corner
      do j=1,nj             !..... loop tets to add
         jt=jtet(j+1,kc)    !..... add tet jt to list, but need to look if it's already there
         ni=itet(5,kt)      !..... neighbors already in the list
         do it=1,ni
            if(itet(5+it,kt) .eq. jt) goto 100 !..... jump if tet is already in the list
         enddo 
         itet(5,kt)=itet(5,kt)+1
         ni=itet(5,kt)
         if(ni.gt.MTET+1) &
            call error_handler(E_ERR, routine, 'error on MPOI', source, revision, revdate)
         itet(ni+5,kt)=jt   !..... add to list
         itet(ni+6,kt)=kt   !..... add ourselves at the end, makes easier searching later

         maxnei=max(maxnei,ni) !..... we could use that later to adjust dimensions

100      continue
      enddo 
   enddo
enddo 
write(0,*)'max number of neighbors: ',maxnei

return
end subroutine g_oplus_pre


!-----------------------------------------------------------------------

!> interpolates grid array u to point (x,y,z)
!> hint: try to avoid entering points (x,y,z) that are not in the grid.
!> istat=1 will be properly returned, but searches in vain are more expensive.

subroutine g_oplus_int_pnt(np, nt, nz, u, x, y, z, output, istat)

integer,  intent(in) :: np
integer,  intent(in) :: nt
integer,  intent(in) :: nz
real(r8), intent(in) :: u(nz,nt,np)
real(r8), intent(in) :: x,y,z
real(r8), intent(out) :: output
integer,  intent(out) :: istat

character(len=*), parameter :: routine = 'g_oplus_int_pnt'

!..... the point of the last interpolation
!      most likely place to find this one
!      next likely one of the neighbors,
!      so we'd start searching here
integer, save :: itl = -1, ipl = -1

real(r8) :: d1, dd 

integer :: i, it, ilook, ni, ip, k

integer  :: tphi, ttheta, tz
integer  :: iphi, itheta, iz
integer  :: dphi, dtheta, dz
integer  :: i1, i2, i3, i4
real(r8) :: u1, u2, u3, u4
real(r8) :: xi, xtmp, yi, ytmp, zi, ztmp

! call g_oplus_pre()

!write(0,*)'searching for ',x,y,z
!write(0,*)'  initial ipl ',ipl

nsearch=0
output=0.0
istat=0
if(itl.lt.0) goto 590 !..... never called before, need full search
!...... look in itl first
it=itl
nsearch=1
ilook=1 
xtmp=x-tlow(1,it)
ytmp=y-tlow(2,it)
ztmp=z-tlow(3,it)
xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
if( (xi.ge.0.0) .and. (yi.ge.0.0) .and. (zi.ge.0.0) .and. ((xi+yi+zi).le.1.0) ) goto 700
!..... now look at the neighbors of the last itl
ilook=2
!..... look for point in neighbors of IT, jump to 700 if found
!write(0,*)'          look for point in neighbors of IT'
ni=itet(5,itl)+1 !..... number of neighbors, and itself, at the end
do i=1,ni
   it=itet(5+i,itl)
   nsearch=nsearch+1 
   xtmp=x-tlow(1,it)
   ytmp=y-tlow(2,it)
   ztmp=z-tlow(3,it)
   xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
   yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
   zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
   if(xi.lt.0.0) goto 94971 
   if(yi.lt.0.0) goto 94971 
   if(zi.lt.0.0) goto 94971 
   if((xi+yi+zi).gt.1.0) goto 94971 
   itl=it 
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
   goto 700

94971 continue

enddo

590   continue

!write(0,*)'          no success yet, full search, descending distance'
      !...... no success yet, full search, descending distance
      if(ipl.le.0) ipl=1 !...... arbitrary start at the beginning
!write(0,*)'          ipl ',ipl,hx(ipl),hy(ipl),hz(ipl)
      ilook=3 
      !...... establish first distance
      dd=(x-hx(ipl))**2+(y-hy(ipl))**2+(z-hz(ipl))**2
600   continue
      ip=ipl
      !..... loop over neighbors
      iz=kpoi(1,ip)
      itheta=kpoi(2,ip)
      iphi=kpoi(3,ip)
      do dphi = -1,1
        do dtheta = -1,1
          do dz = -1,1
            ttheta=itheta+dtheta
            if (ttheta < 1 .or. ttheta > nt) cycle
            !if (ttheta < 1) then
              !ttheta=1
              !iphi=mod(iphi+np/2,np)
            !endif
            !if (ttheta > nt) then
              !ttheta=nt
              !iphi=mod(iphi+np/2,np)
            !endif
            tz=iz+dz
            if (tz < 1 .or. tz > nz) cycle
            tphi=iphi+dphi
            if (tphi < 1) tphi=np
            if (tphi > np) tphi=1
            k=ii(tphi,ttheta,tz)
            nsearch=nsearch+1 
            d1=(x-hx(k))**2+(y-hy(k))**2+(z-hz(k))**2
            if(d1.lt.dd)then
              dd=d1
              ip=k
            endif
          enddo
        enddo
      enddo
      if(ip.ne.ipl)then
      ipl=ip !..... we did get closer, so we are not there yet
!write(0,*)'          ipl ',ipl,hx(ipl),hy(ipl),hz(ipl)
      goto 600
      endif 
!write(0,*)'          we are no longer getting closer'
      !...... we are no longer getting closer. We may not yet be in it at this point,
      !       but at least in one of its tets, so final search is over all associated tets.
      !..... look for point in test belonging to IP, jump to 700 if found
      ni=jtet(1,ipl) !..... number of neighbors, and itself, at the end
      do i=1,ni
      it=jtet(1+i,ipl)
      nsearch=nsearch+1 
      xtmp=x-tlow(1,it)
      ytmp=y-tlow(2,it)
      ztmp=z-tlow(3,it)
      xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
      yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
      zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
      if(xi.lt.0.0) goto 94931 
      if(yi.lt.0.0) goto 94931 
      if(zi.lt.0.0) goto 94931 
      if((xi+yi+zi).gt.1.0) goto 94931 
      itl=it 
!write(0,*)'          ipl ',ipl,hx(ipl),hy(ipl),hz(ipl)
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
      goto 700
94931 continue
      enddo
!write(0,*)'          we could not find it'
      !...... we could not find it, so (x,y,z) should not be in the grid
      !       not fatal in many cases, for example out of the domain for LOS integration
      istat=1
      output=0.0
      return 

700   continue
!write(0,*)'          we found it'
      !.... we only get here if we found the tet4 for this point
      istat=0 
      !........ now interpolate
      !  see D. Kenwright and D. Lane, Proc. of the 6th IEEE Visualization Conference, 1070-2385/95, 1995, eq. 12
      !   or any FEM book for basis functions (1. degree Legendre on unit tetrahedron)
      !..... map index  local in tet4 (1,2,3,4) -->  linear (kpoi) --> iz,ip,it and get corner values
      i1=itet(1,itl)
      u1=u(kpoi(1,i1),kpoi(2,i1),kpoi(3,i1)) 
      i2=itet(2,itl)
      u2=u(kpoi(1,i2),kpoi(2,i2),kpoi(3,i2)) 
      i3=itet(3,itl)
      u3=u(kpoi(1,i3),kpoi(2,i3),kpoi(3,i3)) 
      i4=itet(4,itl)
      u4=u(kpoi(1,i4),kpoi(2,i4),kpoi(3,i4)) 
      output=u1+(u2-u1)*xi+(u3-u1)*yi+(u4-u1)*zi
      return
      end

!-----------------------------------------------------------------------


subroutine g_oplus_int(state_handle, ens_size, np, nt, nz, x, y, z, output, istat)

type(ensemble_type), intent(in) :: state_handle !< ensemble handle for data to interpolate in
integer,             intent(in) :: ens_size     !< number of ensembles

integer,  intent(in) :: np
integer,  intent(in) :: nt
integer,  intent(in) :: nz
real(r8), intent(in) :: x,y,z
real(r8), intent(out) :: output(ens_size)
integer,  intent(out) :: istat(ens_size)

character(len=*), parameter :: routine = 'g_oplus_int'

integer  :: domain_id, var_id
integer  :: dim_index(3)
integer(i8)  :: state_index

!..... the tet4 of the last interpolation
!      most likely place to find this one
!      next likely one of the neighbors,
!      so we'd start searching here

integer, save :: itl = -1 

real(r8) :: d1, dd

integer :: i, it, ilook, ni, kt, t

real(r8) :: oplus(ens_size,4)
real(r8) :: xi, xtmp, yi, ytmp, zi, ztmp

!write(0,*)'searching for ',x,y,z
!write(0,*)'  initial itl ',itl

nsearch=0
output=0.0
istat=0
if(itl.lt.0) goto 590 !..... never called before, need full search
!...... look in itl first
it=itl
nsearch=1
ilook=1 
xtmp=x-tlow(1,it)
ytmp=y-tlow(2,it)
ztmp=z-tlow(3,it)
xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
if( (xi.ge.0.0) .and. (yi.ge.0.0) .and. (zi.ge.0.0) .and. ((xi+yi+zi).le.1.0) ) goto 700
!..... now look at the neighbors of the last itl
ilook=2
!..... look for point in neighbors of IT, jump to 700 if found
!write(0,*)'          look for point in neighbors of IT'
ni=itet(5,itl)+1 !..... number of neighbors, and itself, at the end
do i=1,ni
   it=itet(5+i,itl)
   nsearch=nsearch+1 
   xtmp=x-tlow(1,it)
   ytmp=y-tlow(2,it)
   ztmp=z-tlow(3,it)
   xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
   yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
   zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
   if(xi.lt.0.0) goto 94971 
   if(yi.lt.0.0) goto 94971 
   if(zi.lt.0.0) goto 94971 
   if((xi+yi+zi).gt.1.0) goto 94971 
   itl=it 
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
   goto 700

94971 continue

enddo

590   continue

!write(0,*)'          no success yet, full search, descending distance'
      !...... no success yet, full search, descending distance
      if(itl.le.0) itl=1 !...... arbitrary start at the beginning
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
      ilook=3 
      !...... establish first distance
      dd=(x-tcen(1,itl))**2+(y-tcen(2,itl))**2+(z-tcen(3,itl))**2
600   continue
      it=itl
      ni=itet(5,itl)+1 !...... number of neighbors
      !..... loop over neighbors
      do i=1,ni
      kt=itet(5+i,itl)
      nsearch=nsearch+1 
      d1=(x-tcen(1,kt))**2+(y-tcen(2,kt))**2+(z-tcen(3,kt))**2
      if(d1.lt.dd)then
      dd=d1
      it=kt
      endif
      enddo
      if(it.ne.itl)then
      itl=it !..... we did get closer, so we are not there yet
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
      goto 600
      endif 
!write(0,*)'          we are no longer getting closer'
      !...... we are no longer getting closer. We may not yet be in it at this point,
      !       but at least in one of its neighbors, so final search is over all neighbors.
      !..... look for point in neighbors of IT, jump to 700 if found
      ni=itet(5,itl)+1 !..... number of neighbors, and itself, at the end
      do i=1,ni
      it=itet(5+i,itl)
      nsearch=nsearch+1 
      xtmp=x-tlow(1,it)
      ytmp=y-tlow(2,it)
      ztmp=z-tlow(3,it)
      xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
      yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
      zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
      if(xi.lt.0.0) goto 94931 
      if(yi.lt.0.0) goto 94931 
      if(zi.lt.0.0) goto 94931 
      if((xi+yi+zi).gt.1.0) goto 94931 
      itl=it 
!write(0,*)'          itl ',itl,tlow(1,itl),tlow(2,itl),tlow(3,itl)
      goto 700
94931 continue
      enddo
!write(0,*)'          we could not find it'
      !...... we could not find it, so (x,y,z) should not be in the grid
      !       not fatal in many cases, for example out of the domain for LOS integration
      istat=1
      output=0.0
      return 

700   continue
!write(0,*)'          we found it'
      !.... we only get here if we found the tet4 for this point
      istat=0 
      !........ now interpolate
      !  see D. Kenwright and D. Lane, Proc. of the 6th IEEE Visualization Conference, 1070-2385/95, 1995, eq. 12
      !   or any FEM book for basis functions (1. degree Legendre on unit tetrahedron)
      !..... map index  local in tet4 (1,2,3,4) -->  linear (kpoi) --> iz,ip,it and get corner values

      domain_id = 1
      var_id    = 1

      !>@todo CHECK THIS THOROUGHLY

      do t = 1,4

         dim_index(1) = kpoi(1, itet(t,itl))    ! dim_index(1) is the index into height
         dim_index(2) = kpoi(2, itet(t,itl))    ! dim_index(2) is the index into lat
         dim_index(3) = kpoi(3, itet(t,itl))    ! dim_index(3) is the index into lon
  
         state_index = get_dart_vector_index(dim_index(1), dim_index(2), dim_index(3), domain_id, var_id)

         oplus(:,t) = get_state(state_index, state_handle)

      enddo

      output=oplus(:,1) + (oplus(:,2) - oplus(:,1))*xi + &
                          (oplus(:,3) - oplus(:,1))*yi + &
                          (oplus(:,4) - oplus(:,1))*zi

      return
      end

!-----------------------------------------------------------------------
!>todo this routine does not appear to be used ... remove?

subroutine invert3x3(A, Ainv, Adet) 
! https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/269207

real(r8), intent(in)  ::  A(3,3) 
real(r8), intent(out) ::  Ainv(3,3)
real(r8), intent(out) ::  Adet 

Ainv(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
Ainv(1,2) = A(3,2)*A(1,3) - A(1,2)*A(3,3)
Ainv(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)

Adet = Ainv(1,1)*A(1,1) + Ainv(1,2)*A(2,1) + Ainv(1,3)*A(3,1)

!>@todo given the Fortran standard, I am not sure Adet can ever be smaller than 'tiny' ...
if (Adet < tiny(Adet)) &
   call error_handler(E_ERR,'invert3x3','nope',source,revision,revdate)

Ainv(1,1) = Ainv(1,1)/Adet
Ainv(1,2) = Ainv(1,2)/Adet
Ainv(1,3) = Ainv(1,3)/Adet

Ainv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))/Adet
Ainv(2,2) = (A(1,1)*A(3,3) - A(3,1)*A(1,3))/Adet
Ainv(2,3) = (A(2,1)*A(1,3) - A(1,1)*A(2,3))/Adet
Ainv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))/Adet
Ainv(3,2) = (A(3,1)*A(1,2) - A(1,1)*A(3,2))/Adet
Ainv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/Adet

return
end 

!-----------------------------------------------------------------------


subroutine convert_to_cartesian(phi,theta,height,x,y,z)

real(r8), intent(in)  :: phi
real(r8), intent(in)  :: theta
real(r8), intent(in)  :: height
real(r8), intent(out) :: x,y,z

real(r8) :: r,p,t

r = RE  + height
p = RAD * phi
t = RAD * theta

x = r*cos(p)*cos(t)
y = r*sin(p)*cos(t)
z = r*sin(t)

end subroutine convert_to_cartesian

end module openggcm_interp_mod
