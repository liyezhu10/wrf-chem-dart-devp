SUBROUTINE naturel_loc(ucov,vcov,teta,u,v,t,mode,ps)
   USE exner_hyb_loc_m,   ONLY: exner_hyb_loc
   USE parallel_lmdz 
   USE naturel_mod

!**************************************************************************
!   Tarkeshwar Singh
!   IIT Delhi
!   Email: tarkphysics87@gmail.com
! Perpose: Convert ucov, vcov , teta to natural u, v, t and vice versa
! mode =  0 : (ucov,vcov,teta) to (u,v,t)
! mode =  1 : (u,v,t) to (ucov,vcov,teta)
!**************************************************************************

IMPLICIT NONE
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "comconst.h"
#include "comvert.h"

  INTEGER ::  mode
  REAL :: vcov(ijb_v:ije_v,llm) ! meridional covariant wind
  REAL :: ucov(ijb_u:ije_u,llm) ! zonal covariant wind
  REAL :: teta(ijb_u:ije_u,llm) ! potential temperature (K)
  REAL :: ps(ijb_u:ije_u)       ! surface pressure (Pa)

  REAL :: v(ijb_v:ije_v,llm)    ! meridional wind
  REAL :: u(ijb_u:ije_u,llm)    ! zonal wind
  REAL :: t(ijb_u:ije_u,llm)    ! Temperature (K)

  INTEGER ::  l,ij
  LOGICAL, SAVE :: firstcall=.TRUE.
  INTEGER ijbu,ijeu

     call naturel_allocate
!$OMP BARRIER
     call pression_loc (ijnb_u,ap,bp,ps,preslev)
!$OMP BARRIER
     CALL exner_hyb_loc(ijnb_u,ps,preslev,pkss,pkreff)


!$OMP BARRIER     
      if (pole_nord) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov(ij, l) = 0.
           u   (ij, l) = 0.
           ENDDO
        ENDDO
!$OMP END DO NOWAIT        
      ENDIF

      if (pole_sud) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO  l = 1,llm
          DO  ij = 1,iip1
           ucov( ij +ip1jm, l) = 0.
           u   ( ij +ip1jm, l) = 0.
          ENDDO
        ENDDO
!$OMP END DO NOWAIT      
      ENDIF


      if (mode.eq.0) then
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 

         do l=1,llm
            !--U---
            ijbu=ijb_u
            ijeu=ije_u

            if (pole_nord) ijbu = ijb_u + iip1
            if (pole_sud)  ijeu = ije_u - iip1

            do ij=ijbu,ijeu 
               u(ij,l)=ucov(ij,l)/cu(ij)
            enddo
           !--V---
           do ij=ijb_v, ije_v
              v(ij,l)=vcov(ij,l)/cv(ij)
           enddo
           !--T---
            do ij=ijb_u, ije_u
               t(ij,l)=teta(ij,l)*pkreff(ij,l)/cpp
            enddo
         enddo
!$OMP END DO NOWAIT

      else

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
        do l=1,llm
           !--U---
            ijbu=ijb_u
            ijeu=ije_u

            if (pole_nord) ijbu = ijb_u + iip1
            if (pole_sud)  ijeu = ije_u - iip1

            do ij=ijbu,ijeu 
               ucov(ij,l)=u(ij,l)*cu(ij)
           enddo
           !--V---
           do ij=ijb_v, ije_v
               vcov(ij,l)=v(ij,l)*cv(ij)
           enddo
           !--T---
           do ij= ijb_u, ije_u
               teta(ij,l)=t(ij,l)*cpp/pkreff(ij,l)
           enddo

         enddo !llm
!$OMP END DO NOWAIT

      endif ! mode

      return
      end

