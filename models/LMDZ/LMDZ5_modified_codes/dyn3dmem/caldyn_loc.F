!
! $Id$
!
#undef DEBUG_IO
!#define DEBUG_IO

      SUBROUTINE caldyn_loc
     $ (itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
     $  phi,conser,du,dv,dteta,dp,w,pbaru,pbarv,time )
      USE parallel_lmdz
      USE Write_Field_loc
      USE caldyn_mod, ONLY: vcont, ucont, ang, p, massebx, masseby,
     &                      vorpot, ecin, bern, massebxy, convm

      use control_mod   !! DART-LMDZ
      
      IMPLICIT NONE

!=======================================================================
!
!  Auteur :  P. Le Van
!
!   Objet:
!   ------
!
!   Calcul des tendances dynamiques.
!
! Modif 04/93 F.Forget
!=======================================================================

!-----------------------------------------------------------------------
!   0. Declarations:
!   ----------------

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom.h"

!   Arguments:
!   ----------

      LOGICAL,INTENT(IN) :: conser ! triggers printing some diagnostics ! not used
      INTEGER,INTENT(IN) :: itau ! time step index ! not used
      REAL,INTENT(IN) :: vcov(ijb_v:ije_v,llm) ! covariant meridional wind
      REAL,INTENT(IN) :: ucov(ijb_u:ije_u,llm) ! covariant zonal wind
      REAL,INTENT(IN) :: teta(ijb_u:ije_u,llm) ! potential temperature
      REAL,INTENT(IN) :: ps(ijb_u:ije_u) ! surface pressure
      REAL,INTENT(IN) :: phis(ijb_u:ije_u) ! geopotential at the surface
      REAL,INTENT(IN) :: pk(iip1,jjb_u:jje_u,llm) ! Exner at mid-layer
      REAL,INTENT(IN) :: pkf(ijb_u:ije_u,llm) ! filtered Exner
      REAL,INTENT(IN) :: phi(ijb_u:ije_u,llm) ! geopotential
      REAL,INTENT(OUT) :: masse(ijb_u:ije_u,llm) ! air mass
      REAL,INTENT(OUT) :: dv(ijb_v:ije_v,llm) ! tendency on vcov
      REAL,INTENT(OUT) :: du(ijb_u:ije_u,llm) ! tendency on ucov
      REAL,INTENT(OUT) :: dteta(ijb_u:ije_u,llm) ! tenddency on teta
      REAL,INTENT(OUT) :: dp(ijb_u:ije_u) ! tendency on ps
      REAL,INTENT(OUT) :: w(ijb_u:ije_u,llm) ! vertical velocity
      REAL,INTENT(OUT) :: pbaru(ijb_u:ije_u,llm) ! mass flux in the zonal direction
      REAL,INTENT(OUT) :: pbarv(ijb_v:ije_v,llm) ! mass flux in the meridional direction
      REAL,INTENT(IN) :: time ! current time

!   Local:
!   ------

      INTEGER   ij,l,ijb,ije,ierr


!-----------------------------------------------------------------------
!   Compute dynamical tendencies:
!--------------------------------

      ! compute contravariant winds ucont() and vcont
      CALL covcont_loc  ( llm    , ucov    , vcov , ucont, vcont     )
      ! compute pressure p()
      CALL pression_loc ( ip1jmp1, ap      , bp   ,  ps  , p         )
cym      CALL psextbar (   ps   , psexbarxy                          )
c$OMP BARRIER
      ! compute mass in each atmospheric mesh: masse()
      CALL massdair_loc (    p   , masse                             )
      ! compute X and Y-averages of mass, massebx() and masseby()
      CALL massbar_loc  (   masse, massebx , masseby                 )
      ! compute XY-average of mass, massebxy()
      call massbarxy_loc(   masse, massebxy                          )
      ! compute mass fluxes pbaru() and pbarv()
      CALL flumass_loc  ( massebx, masseby,vcont,ucont,pbaru,pbarv   )
      ! compute dteta() , horizontal converging flux of theta
      CALL dteta1_loc   (   teta , pbaru   , pbarv, dteta            )
      ! compute convm(), horizontal converging flux of mass
      CALL convmas1_loc  (   pbaru, pbarv   , convm                  )
c$OMP BARRIER      
      CALL convmas2_loc  (   convm                      )
c$OMP BARRIER
#ifdef DEBUG_IO
      call WriteField_u('ucont',ucont)
      call WriteField_v('vcont',vcont)
      call WriteField_u('p',p)
      call WriteField_u('masse',masse)
      call WriteField_u('massebx',massebx)
      call WriteField_v('masseby',masseby)
      call WriteField_v('massebxy',massebxy)
      call WriteField_u('pbaru',pbaru)
      call WriteField_v('pbarv',pbarv)
      call WriteField_u('dteta',dteta)
      call WriteField_u('convm',convm)
#endif      

c$OMP BARRIER
c$OMP MASTER
      ijb=ij_begin
      ije=ij_end
      ! compute pressure variation due to mass convergence
      DO ij =ijb, ije
         dp( ij ) = convm( ij,1 ) / airesurg( ij )
      ENDDO
c$OMP END MASTER
c$OMP BARRIER
      
      ! compute vertical velocity w()
      CALL vitvert_loc ( convm  , w                                )
      ! compute potential vorticity vorpot()
      CALL tourpot_loc ( vcov   , ucov  , massebxy  , vorpot       )
      ! compute rotation induced du() and dv()
      CALL dudv1_loc   ( vorpot , pbaru , pbarv     , du     , dv  )

#ifdef DEBUG_IO      
      call WriteField_u('w',w)
      call WriteField_v('vorpot',vorpot)
      call WriteField_u('du',du)
      call WriteField_v('dv',dv)
#endif      
      
      ! compute kinetic energy ecin()
      CALL enercin_loc ( vcov   , ucov  , vcont   , ucont  , ecin  )
      ! compute Bernouilli function bern()
      CALL bernoui_loc ( ip1jmp1, llm   , phi       , ecin   , bern)
      ! compute and add du() and dv() contributions from Bernouilli and pressure
      CALL dudv2_loc   ( teta   , pkf   , bern      , du     , dv  )

#ifdef DEBUG_IO
      call WriteField_u('ecin',ecin)
      call WriteField_u('bern',bern)
      call WriteField_u('du',du)
      call WriteField_v('dv',dv)
      call WriteField_u('pkf',pkf)
#endif
      
      ijb=ij_begin-iip1
      ije=ij_end+iip1
      
      if (pole_nord) ijb=ij_begin
      if (pole_sud) ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)       
      DO l=1,llm
         DO ij=ijb,ije
            ang(ij,l) = ucov(ij,l) + constang(ij)
        ENDDO
      ENDDO
c$OMP END DO

      ! compute vertical advection contributions to du(), dv() and dteta()
      CALL advect_new_loc(ang,vcov,teta,w,massebx,masseby,du,dv,dteta) 

C  WARNING probleme de peridocite de dv sur les PC/linux. Pb d'arrondi 
C          probablement. Observe sur le code compile avec pgf90 3.0-1 
      ijb=ij_begin
      ije=ij_end
      if (pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
      DO l = 1, llm
         DO ij = ijb, ije, iip1
           IF( dv(ij,l).NE.dv(ij+iim,l) )  THEN
c         PRINT *,'!!!ATTENTION!!! probleme de periodicite sur vcov',  
c    ,   ' dans caldyn'
c         PRINT *,' l,  ij = ', l, ij, ij+iim,dv(ij+iim,l),dv(ij,l)
          dv(ij+iim,l) = dv(ij,l)
          endif
         enddo
      enddo
c$OMP END DO NOWAIT     

!********** DART-LMDZ *******************
!! Use only when running with OMP only, will not work with MPI
!$OMP MASTER
      IF( mod(itau,50) .eq. 0 .and. print_sortvar)  THEN
        CALL sortvarc
     .( itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time,vcov )

       ENDIF
!$OMP END MASTER
!********** DART-LMDZ *******************


! Ehouarn: NB: output of control variables not implemented...

      RETURN
      END
