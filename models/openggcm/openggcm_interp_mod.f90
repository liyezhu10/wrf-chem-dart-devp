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
!integer, parameter :: MAXDIM = 100  !... max of each dimension

!..... x,y,z grid point coord
real(r8) :: hx(NPOI_MAX),hy(NPOI_MAX),hz(NPOI_MAX) 

!..... tet mapping matrix (location of center and first corner, 
real(r8) :: tcen(3,NTET_MAX),tlow(3,NTET_MAX),tmap(3,3,NTET_MAX)
real(r8) :: tcen_rtp(3,NTET_MAX)

!..... 4 corners + neicount + neighbors + myself
integer :: itet(MTET+6,NTET_MAX) 
!..... 4 face neighbors (1:x<0,y&z>0; 2:y<0,x&z>0; 3:z<0,x&y>0; 4:x,y,z>0)
integer :: ifacetet(4,NTET_MAX)
!..... for each point, count, and which tet4 it belongs to
integer :: jtet(MPOI+1,NPOI_MAX) 
!..... remember the ip,it,iz indices
integer :: kpoi(3,NPOI_MAX)
!..... number of tet4 tetrahedra, search steps
integer :: ntet, nsearch

!!..... corner id for grid indices
!integer :: ii(MAXDIM,MAXDIM,MAXDIM)

real(r8), parameter :: PI  = 4.0_r8*atan(1.0_r8)
real(r8), parameter :: RAD = PI/180.0_r8
real(r8), parameter :: RE  = 6372.0e3_r8

! Logical to keep track of if we have initialized g_oplus_int
logical, save :: module_initialized = .false.

public :: g_oplus_pre, &
          g_oplus_int, &
          nsearch, &
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
real(r8) :: xi, xtmp, yi, ytmp, zi, ztmp

integer :: ii(np,nt,nz) 

integer :: ip, ip1, it, it1, iz, iz1, k, mytest
integer :: i1, i2, i3, i4, i5, i6, i7, i8
integer :: j, jt, kc, kt, l, ll, maxnei, ni, nj
integer :: mm, nn, npnt, cmatch(3), ctot, jtc, ktc, ktclow

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
      !hx(k) = r*cos(p)*cos(t)
      !hy(k) = r*sin(p)*cos(t)
      !hz(k) = r*sin(t)
      hx(k) = r
      hy(k) = t
      hz(k) = p
   endif

   !... save to map grid values later
   kpoi(1,k) = iz
   kpoi(2,k) = it
   kpoi(3,k) = ip 
   ii(ip,it,iz) = k
enddo
enddo
enddo

!...... create tet4 list
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i4)+hx(i8)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i4)+hy(i8)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i4)+hz(i8)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 90.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i5)+hx(i8)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i5)+hy(i8)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i5)+hz(i8)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 90.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i5)+hx(i6)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i5)+hy(i6)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i5)+hz(i6)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 180.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i4)+hx(i3)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i4)+hy(i3)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i4)+hz(i3)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 180.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i2)+hx(i3)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i2)+hy(i3)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i2)+hz(i3)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 270.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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
      tcen_rtp(1,ntet)=0.25_r8*(hx(i1)+hx(i2)+hx(i6)+hx(i7)) 
      tcen_rtp(2,ntet)=0.25_r8*(hy(i1)+hy(i2)+hy(i6)+hy(i7)) 
      tcen_rtp(3,ntet)=0.25_r8*(hz(i1)+hz(i2)+hz(i6)+hz(i7)) 
      if (ip.eq.np) tcen_rtp(3,ntet) = tcen_rtp(3,ntet) + 270.0_r8*RAD
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
      if (a(3,1).lt.-180.0_r8*RAD) a(3,1)=a(3,1)+360.0_r8*RAD
      if (a(3,2).lt.-180.0_r8*RAD) a(3,2)=a(3,2)+360.0_r8*RAD
      if (a(3,3).lt.-180.0_r8*RAD) a(3,3)=a(3,3)+360.0_r8*RAD
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

!...... for each tet4, find its neighbors.  Need this for fast search.
!       basically invert jtet: the tets belonging to each of one tet's corner points are neighbors.
!       however, there are dups and triplets if they share an edge or face.
ifacetet=-1 ! initialize array of facing tets
maxnei=0
do kt=1,ntet                !..... loop tets
   r = tcen_rtp(1,kt)
   t = tcen_rtp(2,kt)
   p = tcen_rtp(3,kt)
   tcen(1,kt) = r*cos(p)*cos(t)
   tcen(2,kt) = r*sin(p)*cos(t)
   tcen(3,kt) = r*sin(t)
   do ll=1,4                !..... loop corners
      kc=itet(ll,kt)        !..... the global number of that corner
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

         !...check if tet and neighbor tet share a face
         if (kt.ne.jt) then
           npnt=0
           do mm=1,4
             do nn=1,4
               ktc=itet(mm,kt)
               jtc=itet(nn,jt)
               if (ktc.eq.jtc) then
                 npnt=npnt+1
                 cmatch(npnt)=ktc  !save common corner
               endif
             enddo
           enddo
           !...4 face neighbors (1:x<0,y&z>0; 2:y<0,x&z>0; 3:z<0,x&y>0; 4:x,y,z>0)
           if (npnt.eq.3) then
             ktclow=itet(1,kt)
  
             ctot=0
             do nn=1,3
               xtmp=hx(cmatch(nn))-hx(ktclow)
               ytmp=hy(cmatch(nn))-hy(ktclow)
               ztmp=hz(cmatch(nn))-hz(ktclow)
               if (ztmp.lt.-180.0_r8*RAD) ztmp=ztmp+360.0_r8*RAD
               xi=xtmp*tmap(1,1,kt)+ytmp*tmap(1,2,kt)+ztmp*tmap(1,3,kt)
               yi=xtmp*tmap(2,1,kt)+ytmp*tmap(2,2,kt)+ztmp*tmap(2,3,kt)
               zi=xtmp*tmap(3,1,kt)+ytmp*tmap(3,2,kt)+ztmp*tmap(3,3,kt)
               if (xi.gt.0.9_r8) then
                 ctot=ctot+2 !(1,0,0)
               else if (yi.gt.0.9_r8) then
                 ctot=ctot+3 !(0,1,0)
               else if (zi.gt.0.9_r8) then
                 ctot=ctot+4 !(0,0,1)
               else
                 ctot=ctot+1 !(0,0,0)
               endif
             enddo
  
             if (ctot.eq.6) then
               ifacetet(3,kt)=jt !(1,2,3)
             else if (ctot.eq.7) then
               ifacetet(2,kt)=jt !(1,2,4)
             else if (ctot.eq.8) then
               ifacetet(1,kt)=jt !(1,3,4)
             else if (ctot.eq.9) then
               ifacetet(4,kt)=jt !(2,3,4)
             else
               call error_handler(E_ERR, routine, 'ctot error', source, revision, revdate)
             endif
           endif
         endif

100      continue
      enddo 
   enddo
enddo 

return
end subroutine g_oplus_pre


!-----------------------------------------------------------------------
!> interpolates grid array u to point (x,y,z)
!> hint: try to avoid entering points (x,y,z) that are not in the grid.
!> istat=1 will be properly returned, but searches in vain are more expensive.
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
  real(r8) :: r, t, p, rho

  integer :: istatus, i, it, ilook, ni, kt, xyzloc

  real(r8) :: oplus(ens_size,4)
  real(r8) :: xi, xtmp, yi, ytmp, zi, ztmp, xyzi

  rho = sqrt(x**2+y**2)
  r = sqrt(rho**2+z**2)
  t = atan(z,rho)
  p = atan(y,x)
  if (p.lt.0.0_r8) p=p+2.0_r8*pi

  nsearch=0
  output=0.0
  istatus=1    !...default to 'not found'

  if (itl.lt.0) then
    !..... never called before
    itl=1 !...... arbitrary start at the beginning
    ilook=3 
  else
    !...... look in itl first
    it=itl
    nsearch=1
    ilook=1 
    xtmp=r-tlow(1,it)
    ytmp=t-tlow(2,it)
    ztmp=p-tlow(3,it)
    xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
    yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
    zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
    if ((xi.ge.0.0_r8).and.(yi.ge.0.0_r8).and.(zi.ge.0.0_r8).and.((xi+yi+zi).le.1.0_r8)) then
      istatus=0
    else 
      !..... now look at the neighbors of the last itl
      ilook=2
      !..... look for point in neighbors of IT
      ni=itet(5,itl)+1 !..... number of neighbors, and itself, at the end
      do i=1,ni
        it=itet(5+i,itl)
        nsearch=nsearch+1 
        xtmp=r-tlow(1,it)
        ytmp=t-tlow(2,it)
        ztmp=p-tlow(3,it)
        xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
        yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
        zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
        if ((xi.ge.0.0_r8).and.(yi.ge.0.0_r8).and.(zi.ge.0.0_r8).and.((xi+yi+zi).le.1.0_r8)) then
          istatus=0
          itl=it 
          exit
        endif
      enddo
    endif
  endif
  
  if (istatus.ne.0) then
    !...... no success yet, full search, descending distance

    !...... establish first distance
    it=itl
    dd=(x-tcen(1,itl))**2+(y-tcen(2,itl))**2+(z-tcen(3,itl))**2

    do
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
        cycle
      else 
        exit
      endif 
    enddo

    !...... we are no longer getting closer. We may not yet be in it at this point,
    !...... look in itl first
    it=itl
    xtmp=r-tlow(1,it)
    ytmp=t-tlow(2,it)
    ztmp=p-tlow(3,it)
    xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
    yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
    zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
    if ((xi.ge.0.0_r8).and.(yi.ge.0.0_r8).and.(zi.ge.0.0_r8).and.((xi+yi+zi).le.1.0_r8)) then
      istatus=0
    else 
      !       but at least in one of its neighbors, so final search is over all neighbors.
      !..... look for point in neighbors of IT
      ni=itet(5,itl)+1 !..... number of neighbors, and itself, at the end
      do i=1,ni
        it=itet(5+i,itl)
        nsearch=nsearch+1 
        xtmp=r-tlow(1,it)
        ytmp=t-tlow(2,it)
        ztmp=p-tlow(3,it)
        xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
        yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
        zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
        if ((xi.ge.0.0_r8).and.(yi.ge.0.0_r8).and.(zi.ge.0.0_r8).and.((xi+yi+zi).le.1.0_r8)) then
          istatus=0
          itl=it 
          exit
        endif
      enddo
    endif
  endif
    
  if (istatus.ne.0) then
    ! could not find it, so follow tets that share face

    ! The line from C to P must intersect one of the faces of T. Find that face:
    !   map P to the unit TET —> P* (you do that anyways to see if P is in T)
    !   now 4 possibilities, according where P lies relative to T.  
    !   Let (x,y,z)=P* depending across which face of T lies, either
    !     1. x<0, y,z>0
    !     2. y<0, x,z>0
    !     3. z<0, x,y>0
    !     4. x,y,z>0 (x+y+z>1, already tested, otherwise it would be inside T)
    ! Now, that the face is found, it belongs to one and only one neighbor TET.
    ! Repeat with each new TET until proper TET found

    it=itl

    !itp1=-100

    do while(it.ge.0)
      xtmp=r-tlow(1,it)
      ytmp=t-tlow(2,it)
      ztmp=p-tlow(3,it)
      xi=xtmp*tmap(1,1,it)+ytmp*tmap(1,2,it)+ztmp*tmap(1,3,it)
      yi=xtmp*tmap(2,1,it)+ytmp*tmap(2,2,it)+ztmp*tmap(2,3,it)
      zi=xtmp*tmap(3,1,it)+ytmp*tmap(3,2,it)+ztmp*tmap(3,3,it)
      !if ((xi.ge.0.0).and.(yi.ge.0.0).and.(zi.ge.0.0).and.((xi+yi+zi).le.1.0)) then
      !more lenient here to deal with rounding, face/edge points
      if ((xi.ge.-1.0d-10).and.(yi.ge.-1.0d-10).and.(zi.ge.-1.0d-10).and.((xi+yi+zi).le.(1.0d0+1.0d-10))) then
        istatus=0
        itl=it 
        exit
      endif

      !itp2=itp1 
      !itp1=it

      xyzi=amin1(xi,yi,zi)
      if (xyzi.lt.0.0) then
        xyzloc=minloc( (/ xi,yi,zi /), 1 )
        select case (xyzloc)
          case(1)
            it=ifacetet(1,it)
          case(2)
            it=ifacetet(2,it)
          case(3)
            it=ifacetet(3,it)
          case default
            call error_handler(E_ERR, routine, 'face direction error', source, revision, revdate)
        end select
      else if ((xi+yi+zi).gt.1) then
        it=ifacetet(4,it)
      else
        call error_handler(E_ERR, routine, 'invalid face condition', source, revision, revdate)
      endif

      nsearch=nsearch+1 
    enddo

  endif

  istat(:)=istatus

  if (istatus.eq.0) then
      !.... we only get here if we found the tet4 for this point
      !........ now interpolate
      !  see D. Kenwright and D. Lane, Proc. of the 6th IEEE Visualization Conference, 1070-2385/95, 1995, eq. 12
      !   or any FEM book for basis functions (1. degree Legendre on unit tetrahedron)
      !..... map index  local in tet4 (1,2,3,4) -->  linear (kpoi) --> iz,ip,it and get corner values

      domain_id = 1
      var_id    = 1

      do it = 1,4

         dim_index(1) = kpoi(1, itet(it,itl))    ! dim_index(1) is the index into height
         dim_index(2) = kpoi(2, itet(it,itl))    ! dim_index(2) is the index into lat
         dim_index(3) = kpoi(3, itet(it,itl))    ! dim_index(3) is the index into lon
  
         state_index = get_dart_vector_index(dim_index(1), dim_index(2), dim_index(3), domain_id, var_id)

         oplus(:,it) = get_state(state_index, state_handle)

      enddo

      output=oplus(:,1) + (oplus(:,2) - oplus(:,1))*xi + &
                          (oplus(:,3) - oplus(:,1))*yi + &
                          (oplus(:,4) - oplus(:,1))*zi
  endif

  !...... if we could not find it, (x,y,z) should not be in the grid
  !       not fatal in many cases, for example out of the domain for LOS integration

  return
end subroutine g_oplus_int


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
