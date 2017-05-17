! 
! $Id$
!
c
c
#define DEBUG_IO
#undef DEBUG_IO


      SUBROUTINE leapfrog_loc(ucov0,vcov0,teta0,ps0,
     &                        masse0,phis0,q0,clesphy0,
     &                        time_0)

       USE misc_mod
       USE parallel_lmdz
       USE times
       USE mod_hallo
       USE Bands
       USE Write_Field
       USE Write_Field_p
       USE vampir
       USE timer_filtre, ONLY : print_filtre_timer
       USE infotrac
       USE guide_loc_mod, ONLY : guide_main
       USE getparam
       USE control_mod
       USE mod_filtreg_p
       USE write_field_loc
       USE allocate_field_mod
       USE call_dissip_mod, ONLY : call_dissip
       USE call_calfis_mod, ONLY : call_calfis
       USE leapfrog_mod
       use exner_hyb_loc_m, only: exner_hyb_loc
       use exner_milieu_loc_m, only: exner_milieu_loc
      IMPLICIT NONE

c      ......   Version  du 10/01/98    ..........

c             avec  coordonnees  verticales hybrides 
c   avec nouveaux operat. dissipation * ( gradiv2,divgrad2,nxgraro2 )

c=======================================================================
c
c   Auteur:  P. Le Van /L. Fairhead/F.Hourdin
c   -------
c
c   Objet:
c   ------
c
c   GCM LMD nouvelle grille
c
c=======================================================================
c
c  ... Dans inigeom , nouveaux calculs pour les elongations  cu , cv
c      et possibilite d'appeler une fonction f(y)  a derivee tangente
c      hyperbolique a la  place de la fonction a derivee sinusoidale.

c  ... Possibilite de choisir le shema pour l'advection de
c        q  , en modifiant iadv dans traceur.def  (10/02) .
c
c      Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron,10/99)
c      Pour Van-Leer iadv=10 
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comdissnew.h"
#include "comvert.h"
#include "comgeom.h"
#include "logic.h"
#include "temps.h"
#include "ener.h"
#include "description.h"
#include "serre.h"
!#include "com_io_dyn.h"
#include "iniprint.h"
#include "academic.h"
!      include "mpif.h"
      
      INTEGER,PARAMETER :: longcles = 20
      REAL,INTENT(IN) :: clesphy0( longcles ) ! not used
      REAL,INTENT(IN) :: time_0 ! not used

c   dynamical variables:
      REAL,INTENT(IN) :: ucov0(ijb_u:ije_u,llm)    ! zonal covariant wind
      REAL,INTENT(IN) :: vcov0(ijb_v:ije_v,llm)    ! meridional covariant wind
      REAL,INTENT(IN) :: teta0(ijb_u:ije_u,llm)    ! potential temperature
      REAL,INTENT(IN) :: q0(ijb_u:ije_u,llm,nqtot) ! advected tracers
      REAL,INTENT(IN) :: ps0(ijb_u:ije_u)          ! surface pressure (Pa)
      REAL,INTENT(IN) :: masse0(ijb_u:ije_u,llm)   ! air mass
      REAL,INTENT(IN) :: phis0(ijb_u:ije_u)        ! geopotentiat at the surface

      real zqmin,zqmax

!      REAL,SAVE,ALLOCATABLE :: p (:,:  )               ! pression aux interfac.des couches
!      REAL,SAVE,ALLOCATABLE :: pks(:)                      ! exner au  sol
!      REAL,SAVE,ALLOCATABLE :: pk(:,:)                   ! exner au milieu des couches
!      REAL,SAVE,ALLOCATABLE :: pkf(:,:)                  ! exner filt.au milieu des couches
!      REAL,SAVE,ALLOCATABLE :: phi(:,:)                  ! geopotentiel
!      REAL,SAVE,ALLOCATABLE :: w(:,:)                    ! vitesse verticale

c variables dynamiques intermediaire pour le transport
!      REAL,SAVE,ALLOCATABLE :: pbaru(:,:),pbarv(:,:) !flux de masse

c   variables dynamiques au pas -1
!      REAL,SAVE,ALLOCATABLE :: vcovm1(:,:),ucovm1(:,:)
!      REAL,SAVE,ALLOCATABLE :: tetam1(:,:),psm1(:)
!      REAL,SAVE,ALLOCATABLE :: massem1(:,:)

c   tendances dynamiques
!      REAL,SAVE,ALLOCATABLE :: dv(:,:),du(:,:)
!      REAL,SAVE,ALLOCATABLE :: dteta(:,:),dp(:)
!      REAL,DIMENSION(:,:,:), ALLOCATABLE, SAVE :: dq

c   tendances de la dissipation
!      REAL,SAVE,ALLOCATABLE :: dvdis(:,:),dudis(:,:)
!      REAL,SAVE,ALLOCATABLE :: dtetadis(:,:)

c   tendances physiques
      REAL,SAVE,ALLOCATABLE :: dvfi(:,:),dufi(:,:)
      REAL,SAVE,ALLOCATABLE :: dtetafi(:,:)
      REAL,SAVE,ALLOCATABLE :: dpfi(:)
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: dqfi

c   variables pour le fichier histoire
      REAL dtav      ! intervalle de temps elementaire

      REAL tppn(iim),tpps(iim),tpn,tps
c
      INTEGER itau,itaufinp1,iav
!      INTEGER  iday ! jour julien
      REAL       time 

      REAL  SSUM 
!      REAL,SAVE,ALLOCATABLE :: finvmaold(:,:)

cym      LOGICAL  lafin
      LOGICAL :: lafin
      INTEGER ij,iq,l
      INTEGER ik

      real time_step, t_wrt, t_ops

! jD_cur: jour julien courant
! jH_cur: heure julienne courante
      REAL :: jD_cur, jH_cur
      INTEGER :: an, mois, jour
      REAL :: secondes

      logical :: physic
      LOGICAL first,callinigrads

      data callinigrads/.true./
      character*10 string10

!      REAL,SAVE,ALLOCATABLE :: flxw(:,:) ! flux de masse verticale

c+jld variables test conservation energie
!      REAL,SAVE,ALLOCATABLE :: ecin(:,:),ecin0(:,:)
C     Tendance de la temp. potentiel d (theta)/ d t due a la 
C     tansformation d'energie cinetique en energie thermique
C     cree par la dissipation
!      REAL,SAVE,ALLOCATABLE :: dtetaecdt(:,:)
!      REAL,SAVE,ALLOCATABLE :: vcont(:,:),ucont(:,:)
!      REAL,SAVE,ALLOCATABLE :: vnat(:,:),unat(:,:)
      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
      CHARACTER*15 ztit
!!      INTEGER   ip_ebil_dyn  ! PRINT level for energy conserv. diag.
!      SAVE      ip_ebil_dyn
!      DATA      ip_ebil_dyn/0/
c-jld 

      character*80 dynhist_file, dynhistave_file
      character(len=*),parameter :: modname="leapfrog_loc"
      character*80 abort_message


      logical,PARAMETER :: dissip_conservative=.TRUE.
 
      INTEGER testita
      PARAMETER (testita = 9)

      logical , parameter :: flag_verif = .false.
      
c declaration liees au parallelisme
      INTEGER :: ierr
      LOGICAL :: FirstCaldyn
      LOGICAL :: FirstPhysic
      INTEGER :: ijb,ije,j,i
      type(Request) :: TestRequest
      type(Request) :: Request_Dissip
      type(Request) :: Request_physic

      INTEGER :: true_itau
      INTEGER :: iapptrac
      INTEGER :: AdjustCount
!      INTEGER :: var_time
      LOGICAL :: ok_start_timer=.FALSE.
      LOGICAL, SAVE :: firstcall=.TRUE.
      TYPE(distrib),SAVE :: new_dist
!
!******** DART-LMDZ ************ 
      REAL qpn,qps,qppn(iim),qpps(iim) 
      REAL  start_time1
!******** DART-LMDZ ************ 


c$OMP MASTER
      ItCount=0
c$OMP END MASTER      
      true_itau=0
      FirstCaldyn=.TRUE.
      FirstPhysic=.TRUE.
      iapptrac=0
      AdjustCount = 0
      lafin=.false.
      
!      if (nday>=0) then
!         itaufin   = nday*day_step
!      else
!         itaufin   = -nday
!      endif

!*********************** DART-LMDZ **************************** 

      if(sub_day) then

         if(itau_dyn == day_step) itau_dyn =0
         itaufin   = itau_dyn + sub_day_step
         itau      = itau_dyn

      else
         itaufin   = nday*day_step
         itau      = 0
      end if

!*********************** DART-LMDZ **************************** 

      itaufinp1 = itaufin +1

!     itau = 0   !! DART-LMDZ

      physic=.true.
      if (iflag_phys==0.or.iflag_phys==2) physic=.false.
      CALL init_nan
      CALL leapfrog_allocate
      ucov=ucov0
      vcov=vcov0
      teta=teta0
      ps=ps0
      masse=masse0
      phis=phis0
      q=q0
      
!      iday = day_ini+itau/day_step
!      time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!         IF(time.GT.1.) THEN
!          time = time-1.
!          iday = iday+1
!         ENDIF

c Allocate variables depending on dynamic variable nqtot
!$OMP MASTER
      if (firstcall) then
!      
!      ALLOCATE(p(ijb_u:ije_u,llmp1))
!      ALLOCATE(pks(ijb_u:ije_u))
!      ALLOCATE(pk(ijb_u:ije_u,llm))
!      ALLOCATE(pkf(ijb_u:ije_u,llm))
!      ALLOCATE(phi(ijb_u:ije_u,llm))
!      ALLOCATE(w(ijb_u:ije_u,llm))
!      ALLOCATE(pbaru(ip1jmp1,llm),pbarv(ip1jm,llm))
!      ALLOCATE(vcovm1(ijb_v:ije_v,llm),ucovm1(ijb_u:ije_u,llm))
!      ALLOCATE(tetam1(ijb_u:ije_u,llm),psm1(ijb_u:ije_u))
!      ALLOCATE(massem1(ijb_u:ije_u,llm))
!      ALLOCATE(dv(ijb_v:ije_v,llm),du(ijb_u:ije_u,llm))
!      ALLOCATE(dteta(ijb_u:ije_u,llm),dp(ijb_u:ije_u))      
!      ALLOCATE(dvdis(ijb_v:ije_v,llm),dudis(ijb_u:ije_u,llm))
!      ALLOCATE(dtetadis(ijb_u:ije_u,llm))
      ALLOCATE(dvfi(ijb_v:ije_v,llm),dufi(ijb_u:ije_u,llm))
      ALLOCATE(dtetafi(ijb_u:ije_u,llm))
      ALLOCATE(dpfi(ijb_u:ije_u))
!      ALLOCATE(dq(ijb_u:ije_u,llm,nqtot))
      ALLOCATE(dqfi(ijb_u:ije_u,llm,nqtot))
!      ALLOCATE(dqfi_tmp(iip1,llm,nqtot))
!      ALLOCATE(finvmaold(ijb_u:ije_u,llm))
!      ALLOCATE(flxw(ijb_u:ije_u,llm))
!      ALLOCATE(ecin(ijb_u:ije_u,llm),ecin0(ijb_u:ije_u,llm))
!      ALLOCATE(dtetaecdt(ijb_u:ije_u,llm))
!      ALLOCATE(vcont(ijb_v:ije_v,llm),ucont(ijb_u:ije_u,llm))
!      ALLOCATE(vnat(ijb_v:ije_v,llm),unat(ijb_u:ije_u,llm))
      endif
!$OMP END MASTER      
!$OMP BARRIER

!                CALL dynredem1_loc("restart.nc",0.0,
!     &                           vcov,ucov,teta,q,masse,ps)


!***************************** DART-LMDZ ***************************************
!! After data assimilation using DART, the updated fields at poles 
!!  (teta, ps and q at poles) should have average values at pole.
!! Also need to recumpute coloum mass from updated surface pressure  
!-------------------------------------------------------------------------------
      IF(assimilation) THEN
 
        if (pole_nord) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l  =  1, llm
            DO ij =  1,iim
             tppn(ij)  = aire(  ij    ) * teta(  ij    ,l)
            ENDDO
             tpn  = SSUM(iim,tppn,1)/apoln

            DO ij = 1, iip1
             teta(  ij    ,l) = tpn
            ENDDO
          ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER               
          DO ij =  1,iim
            tppn(ij)  = aire(  ij    ) * ps (  ij    )
          ENDDO
            tpn  = SSUM(iim,tppn,1)/apoln
  
          DO ij = 1, iip1
            ps(  ij    ) = tpn
          ENDDO
c$OMP END MASTER
        endif
        
        if (pole_sud) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l  =  1, llm
            DO ij =  1,iim
             tpps(ij)  = aire(ij+ip1jm) * teta(ij+ip1jm,l)
            ENDDO
             tps  = SSUM(iim,tpps,1)/apols

            DO ij = 1, iip1
             teta(ij+ip1jm,l) = tps
            ENDDO
          ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER               
          DO ij =  1,iim
            tpps(ij)  = aire(ij+ip1jm) * ps (ij+ip1jm)
          ENDDO
            tps  = SSUM(iim,tpps,1)/apols
  
          DO ij = 1, iip1
            ps(ij+ip1jm) = tps
          ENDDO
c$OMP END MASTER
        endif

! After data assimilation, the surface pressure in start.nc has been changed
! Hence need to recompute the air mass and vertical air pressure from new 
!surface pressure.  Without doing it, the surface pressure will have 
!decline trend during assimilation


c$OMP BARRIER
      CALL pression_loc ( ijnb_u, ap, bp, ps, p )
c$OMP BARRIER
      CALL massdair_loc (     p  , masse         )
c$OMP BARRIER



!!----- q ---------
c$OMP BARRIER
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
          DO l = 1, llm
           DO ij = ij_begin, ij_end
            deltaPress(ij,l) =  p(ij,l) - p(ij,l+1)
           ENDDO
          ENDDO
c$OMP END DO NOWAIT

c$OMP BARRIER
          CALL qminimum_loc( q, 2, deltaPress )
c$OMP BARRIER
      IF (pole_nord) THEN

        DO iq = 1, 2 !! only Water vapour and liquid are state vector in LMDZ-DART
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm
             DO ij = 1, iim
               qppn(ij) = aire(   ij   ) * q(   ij   ,l,iq)
             ENDDO
               qpn  =  SSUM(iim,qppn,1)/apoln
             DO ij = 1, iip1
               q(   ij   ,l,iq)  = qpn
             ENDDO
          ENDDO
c$OMP END DO NOWAIT
        ENDDO

      ENDIF

      IF (pole_sud) THEN
        DO iq = 1, 2 !! only Water vapour and liquid are state vector in LMDZ-DART
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l = 1, llm

             DO ij = 1, iim
               qpps(ij) = aire(ij+ip1jm) * q(ij+ip1jm,l,iq)
             ENDDO
               qps  =  SSUM(iim,qpps,1)/apols

             DO ij = 1, iip1
               q(ij+ip1jm,l,iq)  = qps
             ENDDO

          ENDDO
c$OMP END DO NOWAIT

        ENDDO

       ENDIF
 
      END IF  ! assimilation
!***************************** DART-LMDZ ***************************************



c-----------------------------------------------------------------------
c   On initialise la pression et la fonction d'Exner :
c   --------------------------------------------------

c$OMP MASTER
      dq(:,:,:)=0.
      CALL pression ( ijnb_u, ap, bp, ps, p       )
c$OMP END MASTER
      if (pressure_exner) then
      CALL exner_hyb_loc( ijnb_u, ps, p, pks, pk, pkf)
      else 
        CALL exner_milieu_loc( ijnb_u, ps, p, pks, pk, pkf )
      endif
c-----------------------------------------------------------------------
c   Debut de l'integration temporelle:
c   ----------------------------------
c et du parallelisme !!

   1  CONTINUE ! Matsuno Forward step begins here

      jD_cur = jD_ref + day_ini - day_ref +                             &
     &          itau/day_step
      jH_cur = jH_ref + start_time +                                    &
     &         mod(itau,day_step)/float(day_step) 
      if (jH_cur > 1.0 ) then
        jD_cur = jD_cur +1.
        jH_cur = jH_cur -1.
      endif


#ifdef CPP_IOIPSL
      if (ok_guide) then
        call guide_main(itau,ucov,vcov,teta,q,masse,ps)
!$OMP BARRIER
      endif
#endif


c
c     IF( MOD( itau, 10* day_step ).EQ.0 )  THEN
c       CALL  test_period ( ucov,vcov,teta,q,p,phis )
c       PRINT *,' ----   Test_period apres continue   OK ! -----', itau
c     ENDIF 
c
cym      CALL SCOPY( ijmllm ,vcov , 1, vcovm1 , 1 )
cym      CALL SCOPY( ijp1llm,ucov , 1, ucovm1 , 1 )
cym      CALL SCOPY( ijp1llm,teta , 1, tetam1 , 1 )
cym      CALL SCOPY( ijp1llm,masse, 1, massem1, 1 )
cym      CALL SCOPY( ip1jmp1, ps  , 1,   psm1 , 1 )

       if (FirstCaldyn) then
c$OMP MASTER
         ucovm1=ucov
         vcovm1=vcov
         tetam1= teta
         massem1= masse
         psm1= ps
         
! Ehouarn: finvmaold is actually not used       
!         finvmaold = masse
c$OMP END MASTER
c$OMP BARRIER
!         CALL filtreg_p ( finvmaold ,jjb_u,jje_u,jjb_u,jje_u,jjp1, llm,
!     &                    -2,2, .TRUE., 1 )
       else
! Save fields obtained at previous time step as '...m1'
         ijb=ij_begin
         ije=ij_end

c$OMP MASTER           
         psm1     (ijb:ije) = ps    (ijb:ije)
c$OMP END MASTER

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
         DO l=1,llm      
           ije=ij_end
           ucovm1   (ijb:ije,l) = ucov  (ijb:ije,l)
           tetam1   (ijb:ije,l) = teta  (ijb:ije,l)
           massem1  (ijb:ije,l) = masse (ijb:ije,l)
!           finvmaold(ijb:ije,l)=masse(ijb:ije,l)
                 
           if (pole_sud) ije=ij_end-iip1
           vcovm1(ijb:ije,l) = vcov  (ijb:ije,l)
       

         ENDDO
c$OMP ENDDO  


! Ehouarn: finvmaold not used
!          CALL filtreg_p(finvmaold ,jjb_u,jje_u,jj_begin,jj_end,jjp1, 
!     .                    llm, -2,2, .TRUE., 1 )

       endif ! of if (FirstCaldyn)
       
      forward = .TRUE.
      leapf   = .FALSE.
      dt      =  dtvr

c   ...    P.Le Van .26/04/94  ....

cym      CALL SCOPY   ( ijp1llm,   masse, 1, finvmaold,     1 )
cym      CALL filtreg ( finvmaold ,jjp1, llm, -2,2, .TRUE., 1 )

cym  ne sert a rien
cym      call minmax(ijp1llm,q(:,:,3),zqmin,zqmax)

   2  CONTINUE ! Matsuno backward or leapfrog step begins here

c$OMP MASTER
      ItCount=ItCount+1
      if (MOD(ItCount,1)==1) then
        debug=.true.
      else
        debug=.false.
      endif
c$OMP END MASTER
c-----------------------------------------------------------------------

c   date:
c   -----


c   gestion des appels de la physique et des dissipations:
c   ------------------------------------------------------
c
c   ...    P.Le Van  ( 6/02/95 )  ....

      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE.
      apdiss = .FALSE.

      IF( purmats ) THEN
      ! Purely Matsuno time stepping
         IF( MOD(itau,iconser) .EQ.0.AND.  forward    ) conser = .TRUE.
         IF( MOD(itau,dissip_period ).EQ.0.AND..NOT.forward ) 
     s        apdiss = .TRUE.
         IF( MOD(itau,iphysiq ).EQ.0.AND..NOT.forward 
     s          .and. physic                        ) apphys = .TRUE.
      ELSE
      ! Leapfrog/Matsuno time stepping 
         IF( MOD(itau   ,iconser) .EQ. 0              ) conser = .TRUE.
         IF( MOD(itau+1,dissip_period).EQ.0 .AND. .NOT. forward )
     s        apdiss = .TRUE.
         IF( MOD(itau+1,iphysiq).EQ.0.AND.physic) apphys=.TRUE.
      END IF

! Ehouarn: for Shallow Water case (ie: 1 vertical layer),
!          supress dissipation step
      if (llm.eq.1) then
        apdiss=.false.
      endif

cym    ---> Pour le moment      
cym      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE. ! ie: no output of control variables to stdout in //
      
      if (firstCaldyn) then
c$OMP MASTER
          call Set_Distrib(distrib_caldyn)
c$OMP END MASTER
c$OMP BARRIER
          firstCaldyn=.FALSE.
cym          call InitTime
c$OMP MASTER
          call Init_timer
c$OMP END MASTER
      endif

c$OMP MASTER      
      IF (ok_start_timer) THEN
        CALL InitTime
        ok_start_timer=.FALSE.
      ENDIF      
c$OMP END MASTER      


!ym  PAS D'AJUSTEMENT POUR LE MOMENT     
      if (Adjust) then
        AdjustCount=AdjustCount+1
!        if (iapptrac==iapp_tracvl .and. (forward. OR . leapf)
!     &         .and. itau/iphysiq>2 .and. Adjustcount>30) then
        if (Adjustcount>1) then
           AdjustCount=0
c$OMP MASTER 
           call allgather_timer_average

        if (prt_level > 9) then
        
        print *,'*********************************'
        print *,'******    TIMER CALDYN     ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_caldyn(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_caldyn(i),timer_caldyn,i),
     &            '+-',timer_delta(jj_nb_caldyn(i),timer_caldyn,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER VANLEER    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_vanleer(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_vanleer(i),timer_vanleer,i),
     &            '+-',timer_delta(jj_nb_vanleer(i),timer_vanleer,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER DISSIP    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_dissip(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_dissip(i),timer_dissip,i),
     &             '+-',timer_delta(jj_nb_dissip(i),timer_dissip,i)
        enddo
        
!        if (mpi_rank==0) call WriteBands
        
       endif
       
         call AdjustBands_caldyn(new_dist)
!$OMP END MASTER
!$OMP BARRIER
         CALL leapfrog_switch_caldyn(new_dist)
!$OMP BARRIER


!$OMP MASTER
         distrib_caldyn=new_dist
         CALL set_distrib(distrib_caldyn)
!$OMP END MASTER
!$OMP BARRIER
!         call Register_SwapFieldHallo(ucov,ucov,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(ucovm1,ucovm1,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(vcov,vcov,ip1jm,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(vcovm1,vcovm1,ip1jm,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(teta,teta,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(tetam1,tetam1,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(masse,masse,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(massem1,massem1,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(ps,ps,ip1jmp1,1,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(psm1,psm1,ip1jmp1,1,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(pkf,pkf,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(pk,pk,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(pks,pks,ip1jmp1,1,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(phis,phis,ip1jmp1,1,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(phi,phi,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(finvmaold,finvmaold,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)
! 
!        do j=1,nqtot
!         call Register_SwapFieldHallo(q(:,:,j),q(:,:,j),ip1jmp1,llm,
!     &                                jj_nb_caldyn,0,0,TestRequest)
!        enddo
!
!         call Set_Distrib(distrib_caldyn)
!         call SendRequest(TestRequest)
!         call WaitRequest(TestRequest)
         
!$OMP MASTER
        call AdjustBands_dissip(new_dist)
!$OMP END MASTER
!$OMP BARRIER
        CALL leapfrog_switch_dissip(new_dist)
!$OMP BARRIER
!$OMP MASTER
        distrib_dissip=new_dist
!$OMP END MASTER
!$OMP BARRIER
!        call AdjustBands_physic

c$OMP MASTER  
        if (mpi_rank==0) call WriteBands
c$OMP END MASTER  


      endif
      endif       
     
      
      
c-----------------------------------------------------------------------
c   calcul des tendances dynamiques:
c   --------------------------------
c$OMP BARRIER
c$OMP MASTER
       call VTb(VThallo)
c$OMP END MASTER

       call Register_Hallo_u(ucov,llm,1,1,1,1,TestRequest)
       call Register_Hallo_v(vcov,llm,1,1,1,1,TestRequest)
       call Register_Hallo_u(teta,llm,1,1,1,1,TestRequest)
       call Register_Hallo_u(ps,1,1,2,2,1,TestRequest)
       call Register_Hallo_u(pkf,llm,1,1,1,1,TestRequest)
       call Register_Hallo_u(pk,llm,1,1,1,1,TestRequest)
       call Register_Hallo_u(pks,1,1,1,1,1,TestRequest)
       call Register_Hallo_u(p,llmp1,1,1,1,1,TestRequest)
       
c       do j=1,nqtot
c         call Register_Hallo(q(1,1,j),ip1jmp1,llm,1,1,1,1,
c     *                       TestRequest)
c        enddo

       call SendRequest(TestRequest)
c$OMP BARRIER
       call WaitRequest(TestRequest)

c$OMP MASTER
       call VTe(VThallo)
c$OMP END MASTER
c$OMP BARRIER
      
      if (debug) then        
        call WriteField_u('ucov',ucov)
        call WriteField_v('vcov',vcov)
        call WriteField_u('teta',teta)
        call WriteField_u('ps',ps)
        call WriteField_u('masse',masse)
        call WriteField_u('pk',pk)
        call WriteField_u('pks',pks)
        call WriteField_u('pkf',pkf)
        call WriteField_u('phis',phis)
        do j=1,nqtot
          call WriteField_u('q'//trim(int2str(j)),
     .                q(:,:,j))
        enddo
      endif

      
      True_itau=True_itau+1

c$OMP MASTER
      IF (prt_level>9) THEN
        WRITE(lunout,*)"leapfrog_p: Iteration No",True_itau
      ENDIF


      call start_timer(timer_caldyn)

      ! compute geopotential phi()
      CALL geopot_loc  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )

      
      call VTb(VTcaldyn)
c$OMP END MASTER
!      var_time=time+iday-day_ini

c$OMP BARRIER
!      CALL FTRACE_REGION_BEGIN("caldyn")
      time = jD_cur + jH_cur 
      CALL caldyn_loc 
     $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
     $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, time )

!      CALL FTRACE_REGION_END("caldyn")

c$OMP MASTER
      call VTe(VTcaldyn)
c$OMP END MASTER      

#ifdef DEBUG_IO    
      call WriteField_u('du',du)
      call WriteField_v('dv',dv)
      call WriteField_u('dteta',dteta)
      call WriteField_u('dp',dp)
      call WriteField_u('w',w)
      call WriteField_u('pbaru',pbaru)
      call WriteField_v('pbarv',pbarv)
      call WriteField_u('p',p)
      call WriteField_u('masse',masse)
      call WriteField_u('pk',pk)
#endif
c-----------------------------------------------------------------------
c   calcul des tendances advection des traceurs (dont l'humidite)
c   -------------------------------------------------------------

      
      IF( forward. OR . leapf )  THEN
! Ehouarn: NB: fields sent to advtrac are those at the beginning of the time step
         CALL caladvtrac_loc(q,pbaru,pbarv,
     *        p, masse, dq,  teta,
     .        flxw,pk, iapptrac)

!      do j=1,nqtot
!        call WriteField_u('qadv'//trim(int2str(j)),q(:,:,j))
!      enddo

! Ehouarn: Storage of mass flux for off-line tracers... not implemented...

      ENDIF ! of IF( forward. OR . leapf )


c-----------------------------------------------------------------------
c   integrations dynamique et traceurs:
c   ----------------------------------

c$OMP MASTER 
       call VTb(VTintegre)
c$OMP END MASTER
#ifdef DEBUG_IO    
      if (true_itau>20) then
      call WriteField_u('ucovm1',ucovm1)
      call WriteField_v('vcovm1',vcovm1)
      call WriteField_u('tetam1',tetam1)
      call WriteField_u('psm1',psm1)
      call WriteField_u('ucov_int',ucov)
      call WriteField_v('vcov_int',vcov)
      call WriteField_u('teta_int',teta)
      call WriteField_u('ps_int',ps)
      endif
#endif
c$OMP BARRIER
!       CALL FTRACE_REGION_BEGIN("integrd")

       CALL integrd_loc ( 2,vcovm1,ucovm1,tetam1,psm1,massem1 ,
     $         dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps,masse,phis)
!     $              finvmaold                                    )

!       CALL FTRACE_REGION_END("integrd")
c$OMP BARRIER
#ifdef DEBUG_IO    
      call WriteField_u('ucovm1',ucovm1)
      call WriteField_v('vcovm1',vcovm1)
      call WriteField_u('tetam1',tetam1)
      call WriteField_u('psm1',psm1)
      call WriteField_u('ucov_int',ucov)
      call WriteField_v('vcov_int',vcov)
      call WriteField_u('teta_int',teta)
      call WriteField_u('ps_int',ps)
#endif    
c      do j=1,nqtot
c        call WriteField_p('q'//trim(int2str(j)),
c     .                reshape(q(:,:,j),(/iip1,jmp1,llm/)))
c        call WriteField_p('dq'//trim(int2str(j)),
c     .                reshape(dq(:,:,j),(/iip1,jmp1,llm/)))
c      enddo


c$OMP MASTER 
       call VTe(VTintegre)
c$OMP END MASTER
c .P.Le Van (26/04/94  ajout de  finvpold dans l'appel d'integrd)
c
c-----------------------------------------------------------------------
c   calcul des tendances physiques:
c   -------------------------------
c    ########   P.Le Van ( Modif le  6/02/95 )   ###########
c
       IF( purmats )  THEN
          IF( itau.EQ.itaufin.AND..NOT.forward ) lafin = .TRUE.
       ELSE
          IF( itau+1. EQ. itaufin )              lafin = .TRUE.
       ENDIF

cc$OMP END PARALLEL

c
c
       IF( apphys )  THEN
            
c$OMP MASTER
       PRINT*,"******************************************"
       PRINT*,"TKKK itau,jH_cur =",itau,jH_cur
       PRINT*,"******************************************"
c$OMP END MASTER
 
         CALL call_calfis(itau,lafin,clesphy0,ucov,vcov,teta,masse,ps,  
     &                     phis,q,flxw)
! #ifdef DEBUG_IO    
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('pfi',p)
!         call WriteField_u('pkfi',pk)
!         do j=1,nqtot
!           call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!         enddo
! #endif
! c
! c     .......   Ajout   P.Le Van ( 17/04/96 )   ...........
! c
! cc$OMP PARALLEL DEFAULT(SHARED)
! cc$OMP+         PRIVATE(rdaym_ini,rdayvrai,ijb,ije)

! c$OMP MASTER
!          call suspend_timer(timer_caldyn)

!          write(lunout,*)
!      &   'leapfrog_p: Entree dans la physique : Iteration No ',true_itau
! c$OMP END MASTER

!          CALL pression_loc (  ip1jmp1, ap, bp, ps,  p      )

! c$OMP BARRIER
!          CALL exner_hyb_loc(  ip1jmp1, ps, p,pks, pk, pkf )
! c$OMP BARRIER
!            jD_cur = jD_ref + day_ini - day_ref
!      $        + int (itau * dtvr / daysec) 
!            jH_cur = jH_ref +                                            &
!      &              (itau * dtvr / daysec - int(itau * dtvr / daysec)) 
! !         call ju2ymds(jD_cur+jH_cur, an, mois, jour, secondes)

! c rajout debug
! c       lafin = .true.


! c   Inbterface avec les routines de phylmd (phymars ... )
! c   -----------------------------------------------------

! c+jld

! c  Diagnostique de conservation de l'energie : initialisation
!  
! c-jld
! c$OMP BARRIER
! c$OMP MASTER
!         call VTb(VThallo)
! c$OMP END MASTER

! #ifdef DEBUG_IO    
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('pfi',p)
!         call WriteField_u('pkfi',pk)
! #endif
!         call SetTag(Request_physic,800)
!         
!         call Register_SwapField_u(ucov,ucov,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_v(vcov,vcov,distrib_physic,
!      *                            Request_physic,up=2,down=2)

!         call Register_SwapField_u(teta,teta,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(masse,masse,distrib_physic,
!      *                            Request_physic,up=1,down=2)

!         call Register_SwapField_u(p,p,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(pk,pk,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(phis,phis,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(phi,phi,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(w,w,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call Register_SwapField_u(q,q,distrib_physic,
!      *                            Request_physic,up=2,down=2)

!         call Register_SwapField_u(flxw,flxw,distrib_physic,
!      *                            Request_physic,up=2,down=2)
!         
!         call SendRequest(Request_Physic)
! c$OMP BARRIER
!         call WaitRequest(Request_Physic)       

! c$OMP BARRIER
! c$OMP MASTER
!         call Set_Distrib(distrib_Physic)
!         call VTe(VThallo)
!         
!         call VTb(VTphysiq)
! c$OMP END MASTER
! c$OMP BARRIER

! #ifdef DEBUG_IO    
!       call WriteField_u('ucovfi',ucov)
!       call WriteField_v('vcovfi',vcov)
!       call WriteField_u('tetafi',teta)
!       call WriteField_u('pfi',p)
!       call WriteField_u('pkfi',pk)
!       do j=1,nqtot
!         call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!       enddo
! #endif
!        STOP
! c$OMP BARRIER
! !        CALL FTRACE_REGION_BEGIN("calfis")
!         CALL calfis_loc(lafin ,jD_cur, jH_cur,
!      $               ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,
!      $               du,dv,dteta,dq,
!      $               flxw,
!      $               clesphy0, dufi,dvfi,dtetafi,dqfi,dpfi  )
! !        CALL FTRACE_REGION_END("calfis")
! !        ijb=ij_begin
! !        ije=ij_end  
! !        if ( .not. pole_nord) then
! !c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
! !          DO l=1,llm
! !          dufi_tmp(1:iip1,l)   = dufi(ijb:ijb+iim,l) 
! !          dvfi_tmp(1:iip1,l)   = dvfi(ijb:ijb+iim,l)  
! !          dtetafi_tmp(1:iip1,l)= dtetafi(ijb:ijb+iim,l)  
! !          dqfi_tmp(1:iip1,l,:) = dqfi(ijb:ijb+iim,l,:)  
! !          ENDDO
! !c$OMP END DO NOWAIT
! !
! !c$OMP MASTER
! !          dpfi_tmp(1:iip1)     = dpfi(ijb:ijb+iim)  
! !c$OMP END MASTER
! !        endif ! of if ( .not. pole_nord)

! !c$OMP BARRIER
! !c$OMP MASTER
! !        call Set_Distrib(distrib_physic_bis)

! !        call VTb(VThallo)
! !c$OMP END MASTER
! !c$OMP BARRIER
! ! 
! !        call Register_Hallo_u(dufi,llm,
! !     *                      1,0,0,1,Request_physic)
! !        
! !        call Register_Hallo_v(dvfi,llm,
! !     *                      1,0,0,1,Request_physic)
! !        
! !        call Register_Hallo_u(dtetafi,llm,
! !     *                      1,0,0,1,Request_physic)
! !
! !        call Register_Hallo_u(dpfi,1,
! !     *                      1,0,0,1,Request_physic)
! !
! !        do j=1,nqtot
! !          call Register_Hallo_u(dqfi(ijb_u,1,j),llm,
! !     *                        1,0,0,1,Request_physic)
! !        enddo
! !        
! !        call SendRequest(Request_Physic)
! !c$OMP BARRIER
! !        call WaitRequest(Request_Physic)
! !             
! !c$OMP BARRIER
! !c$OMP MASTER
! !        call VTe(VThallo)
! ! 
! !        call set_Distrib(distrib_Physic)
! !c$OMP END MASTER
! !c$OMP BARRIER        
! !                ijb=ij_begin
! !        if (.not. pole_nord) then
! !        
! !c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
! !          DO l=1,llm
! !            dufi(ijb:ijb+iim,l) = dufi(ijb:ijb+iim,l)+dufi_tmp(1:iip1,l)
! !            dvfi(ijb:ijb+iim,l) = dvfi(ijb:ijb+iim,l)+dvfi_tmp(1:iip1,l) 
! !            dtetafi(ijb:ijb+iim,l) = dtetafi(ijb:ijb+iim,l)
! !     &                              +dtetafi_tmp(1:iip1,l)
! !            dqfi(ijb:ijb+iim,l,:) = dqfi(ijb:ijb+iim,l,:)
! !     &                              + dqfi_tmp(1:iip1,l,:)
! !          ENDDO
! !c$OMP END DO NOWAIT
! !
! !c$OMP MASTER
! !          dpfi(ijb:ijb+iim)   = dpfi(ijb:ijb+iim)+ dpfi_tmp(1:iip1)
! !c$OMP END MASTER
! !          
! !        endif ! of if (.not. pole_nord)

! #ifdef DEBUG_IO           
!         call WriteField_u('dufi',dufi)
!         call WriteField_v('dvfi',dvfi) 
!         call WriteField_u('dtetafi',dtetafi)
!         call WriteField_u('dpfi',dpfi)
!         do j=1,nqtot
!           call WriteField_u('dqfi'//trim(int2str(j)),dqfi(:,:,j))
!        enddo
! #endif

! c$OMP BARRIER

! c      ajout des tendances physiques:
! c      ------------------------------
! #ifdef DEBUG_IO    
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('psfi',ps)
!         do j=1,nqtot
!           call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!        enddo
! #endif

!          IF (ok_strato) THEN
!            CALL top_bound_loc( vcov,ucov,teta,masse,dufi,dvfi,dtetafi)
!          ENDIF

! #ifdef DEBUG_IO           
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('psfi',ps)
!         do j=1,nqtot
!           call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!        enddo
! #endif

!           CALL addfi_loc( dtphys, leapf, forward   ,
!      $                  ucov, vcov, teta , q   ,ps ,
!      $                 dufi, dvfi, dtetafi , dqfi ,dpfi  )

! #ifdef DEBUG_IO    
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('psfi',ps)
!         do j=1,nqtot
!           call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!        enddo
! #endif

! c$OMP BARRIER
! c$OMP MASTER
!         call VTe(VTphysiq)

!         call VTb(VThallo)
! c$OMP END MASTER

!         call SetTag(Request_physic,800)
!         call Register_SwapField_u(ucov,ucov,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_v(vcov,vcov,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(teta,teta,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(masse,masse,
!      *                               distrib_caldyn,Request_physic)

!         call Register_SwapField_u(p,p,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(pk,pk,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(phis,phis,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(phi,phi,
!      *                               distrib_caldyn,Request_physic)
!         
!         call Register_SwapField_u(w,w,
!      *                               distrib_caldyn,Request_physic)

!         call Register_SwapField_u(q,q,
!      *                               distrib_caldyn,Request_physic)
!         
!         call SendRequest(Request_Physic)
! c$OMP BARRIER
!         call WaitRequest(Request_Physic)     

! c$OMP BARRIER
! c$OMP MASTER
!        call VTe(VThallo)
!        call set_distrib(distrib_caldyn)
! c$OMP END MASTER
! c$OMP BARRIER
! c
! c  Diagnostique de conservation de l'energie : difference
!       IF (ip_ebil_dyn.ge.1 ) THEN 
!           ztit='bil phys'
!           CALL diagedyn(ztit,2,1,1,dtphys
!      e  , ucov    , vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
!       ENDIF 

! #ifdef DEBUG_IO    
!         call WriteField_u('ucovfi',ucov)
!         call WriteField_v('vcovfi',vcov)
!         call WriteField_u('tetafi',teta)
!         call WriteField_u('psfi',ps)
!         do j=1,nqtot
!           call WriteField_u('qfi'//trim(int2str(j)),q(:,:,j))
!        enddo
! #endif


! c-jld
c$OMP MASTER
         if (FirstPhysic) then
           ok_start_timer=.TRUE.
           FirstPhysic=.false.
         endif
c$OMP END MASTER
       ENDIF ! of IF( apphys )

      IF(iflag_phys.EQ.2) THEN ! "Newtonian" case
c$OMP MASTER
         if (FirstPhysic) then
           ok_start_timer=.TRUE.
           FirstPhysic=.false.
         endif
c$OMP END MASTER


c   Calcul academique de la physique = Rappel Newtonien + fritcion 
c   --------------------------------------------------------------
cym       teta(:,:)=teta(:,:)
cym     s  -iphysiq*dtvr*(teta(:,:)-tetarappel(:,:))/taurappel
       ijb=ij_begin
       ije=ij_end
!LF       teta(ijb:ije,:)=teta(ijb:ije,:)
!LF     s  -iphysiq*dtvr*(teta(ijb:ije,:)-tetarappel(ijb:ije,:))/taurappel
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
       do l=1,llm
       teta(ijb:ije,l)=teta(ijb:ije,l) -dtvr*
     &        (teta(ijb:ije,l)-tetarappel(ijb:ije,l))*
     &                 (knewt_g+knewt_t(l)*clat4(ijb:ije))        
       enddo
!$OMP END DO

!$OMP MASTER
       if (planet_type.eq."giant") then
         ! add an intrinsic heat flux at the base of the atmosphere
         teta(ijb:ije,1) = teta(ijb:ije,1)
     &        + dtvr * aire(ijb:ije) * ihf / cpp / masse(ijb:ije,1)
       endif
!$OMP END MASTER
!$OMP BARRIER


       call Register_Hallo_u(ucov,llm,0,1,1,0,Request_Physic)
       call Register_Hallo_v(vcov,llm,1,1,1,1,Request_Physic)
       call SendRequest(Request_Physic)
c$OMP BARRIER
       call WaitRequest(Request_Physic)     
c$OMP BARRIER
       call friction_loc(ucov,vcov,dtvr)
!$OMP BARRIER

        ! Sponge layer (if any)
        IF (ok_strato) THEN
          CALL top_bound_loc(vcov,ucov,teta,masse,dtvr)
!$OMP BARRIER
        ENDIF ! of IF (ok_strato) 
      ENDIF ! of IF(iflag_phys.EQ.2)


        CALL pression_loc ( ip1jmp1, ap, bp, ps, p                  )
c$OMP BARRIER
        if (pressure_exner) then
        CALL exner_hyb_loc( ijnb_u, ps, p, pks, pk, pkf )
        else 
          CALL exner_milieu_loc( ijnb_u, ps, p, pks, pk, pkf )
        endif
c$OMP BARRIER
        CALL massdair_loc(p,masse)
c$OMP BARRIER

cc$OMP END PARALLEL

c-----------------------------------------------------------------------
c   dissipation horizontale et verticale  des petites echelles:
c   ----------------------------------------------------------

      IF(apdiss) THEN
      
        CALL call_dissip(ucov,vcov,teta,p,pk,ps)
!cc$OMP  PARALLEL DEFAULT(SHARED) 
!cc$OMP+          PRIVATE(ijb,ije,tppn,tpn,tpps,tps)
!c$OMP MASTER
!        call suspend_timer(timer_caldyn)
!        
!c       print*,'Entree dans la dissipation : Iteration No ',true_itau
!c   calcul de l'energie cinetique avant dissipation
!c       print *,'Passage dans la dissipation'

!        call VTb(VThallo)
!c$OMP END MASTER

!c$OMP BARRIER

!        call Register_SwapField_u(ucov,ucov,distrib_dissip,
!     *                            Request_dissip,up=1,down=1)

!        call Register_SwapField_v(vcov,vcov,distrib_dissip,
!     *                            Request_dissip,up=1,down=1)

!        call Register_SwapField_u(teta,teta,distrib_dissip,
!     *                            Request_dissip)

!        call Register_SwapField_u(p,p,distrib_dissip,
!     *                            Request_dissip)

!        call Register_SwapField_u(pk,pk,distrib_dissip,
!     *                            Request_dissip)

!        call SendRequest(Request_dissip)       
!c$OMP BARRIER
!        call WaitRequest(Request_dissip)       

!c$OMP BARRIER
!c$OMP MASTER
!        call set_distrib(distrib_dissip)
!        call VTe(VThallo)
!        call VTb(VTdissipation)
!        call start_timer(timer_dissip)
!c$OMP END MASTER
!c$OMP BARRIER

!        call covcont_loc(llm,ucov,vcov,ucont,vcont)
!        call enercin_loc(vcov,ucov,vcont,ucont,ecin0)

!c   dissipation

!!        CALL FTRACE_REGION_BEGIN("dissip")
!        CALL dissip_loc(vcov,ucov,teta,p,dvdis,dudis,dtetadis)

!#ifdef DEBUG_IO    
!        call WriteField_u('dudis',dudis)
!        call WriteField_v('dvdis',dvdis)
!        call WriteField_u('dtetadis',dtetadis)
!#endif
! 
!!      CALL FTRACE_REGION_END("dissip")
!         
!        ijb=ij_begin
!        ije=ij_end
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
!        DO l=1,llm
!          ucov(ijb:ije,l)=ucov(ijb:ije,l)+dudis(ijb:ije,l)
!        ENDDO
!c$OMP END DO NOWAIT        
!        if (pole_sud) ije=ije-iip1
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
!        DO l=1,llm
!          vcov(ijb:ije,l)=vcov(ijb:ije,l)+dvdis(ijb:ije,l)
!        ENDDO
!c$OMP END DO NOWAIT        

!c       teta=teta+dtetadis


!c------------------------------------------------------------------------
!        if (dissip_conservative) then
!C       On rajoute la tendance due a la transform. Ec -> E therm. cree
!C       lors de la dissipation
!c$OMP BARRIER
!c$OMP MASTER
!            call suspend_timer(timer_dissip)
!            call VTb(VThallo)
!c$OMP END MASTER
!            call Register_Hallo_u(ucov,llm,1,1,1,1,Request_Dissip)
!            call Register_Hallo_v(vcov,llm,1,1,1,1,Request_Dissip)
!            call SendRequest(Request_Dissip)
!c$OMP BARRIER
!            call WaitRequest(Request_Dissip)
!c$OMP MASTER
!            call VTe(VThallo)
!            call resume_timer(timer_dissip)
!c$OMP END MASTER
!c$OMP BARRIER            
!            call covcont_loc(llm,ucov,vcov,ucont,vcont)
!            call enercin_loc(vcov,ucov,vcont,ucont,ecin)
!            
!            ijb=ij_begin
!            ije=ij_end
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
!            do l=1,llm
!              do ij=ijb,ije
!                dtetaecdt(ij,l)= (ecin0(ij,l)-ecin(ij,l))/ pk(ij,l)
!                dtetadis(ij,l)=dtetadis(ij,l)+dtetaecdt(ij,l)
!              enddo
!            enddo
!c$OMP END DO NOWAIT            
!       endif

!       ijb=ij_begin
!       ije=ij_end
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
!         do l=1,llm
!           do ij=ijb,ije
!              teta(ij,l)=teta(ij,l)+dtetadis(ij,l)
!           enddo
!         enddo
!c$OMP END DO NOWAIT         
!c------------------------------------------------------------------------


!c    .......        P. Le Van (  ajout  le 17/04/96  )   ...........
!c   ...      Calcul de la valeur moyenne, unique de h aux poles  .....
!c

!        ijb=ij_begin
!        ije=ij_end
!         
!        if (pole_nord) then
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
!          DO l  =  1, llm
!            DO ij =  1,iim
!             tppn(ij)  = aire(  ij    ) * teta(  ij    ,l)
!            ENDDO
!             tpn  = SSUM(iim,tppn,1)/apoln

!            DO ij = 1, iip1
!             teta(  ij    ,l) = tpn
!            ENDDO
!          ENDDO
!c$OMP END DO NOWAIT

!c$OMP MASTER               
!          DO ij =  1,iim
!            tppn(ij)  = aire(  ij    ) * ps (  ij    )
!          ENDDO
!            tpn  = SSUM(iim,tppn,1)/apoln
!  
!          DO ij = 1, iip1
!            ps(  ij    ) = tpn
!          ENDDO
!c$OMP END MASTER
!        endif
!        
!        if (pole_sud) then
!c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
!          DO l  =  1, llm
!            DO ij =  1,iim
!             tpps(ij)  = aire(ij+ip1jm) * teta(ij+ip1jm,l)
!            ENDDO
!             tps  = SSUM(iim,tpps,1)/apols

!            DO ij = 1, iip1
!             teta(ij+ip1jm,l) = tps
!            ENDDO
!          ENDDO
!c$OMP END DO NOWAIT

!c$OMP MASTER               
!          DO ij =  1,iim
!            tpps(ij)  = aire(ij+ip1jm) * ps (ij+ip1jm)
!          ENDDO
!            tps  = SSUM(iim,tpps,1)/apols
!  
!          DO ij = 1, iip1
!            ps(ij+ip1jm) = tps
!          ENDDO
!c$OMP END MASTER
!        endif


!c$OMP BARRIER
!c$OMP MASTER
!        call VTe(VTdissipation)

!        call stop_timer(timer_dissip)
!        
!        call VTb(VThallo)
!c$OMP END MASTER
!        call Register_SwapField_u(ucov,ucov,distrib_caldyn,
!     *                            Request_dissip)

!        call Register_SwapField_v(vcov,vcov,distrib_caldyn,
!     *                            Request_dissip)

!        call Register_SwapField_u(teta,teta,distrib_caldyn,
!     *                            Request_dissip)

!        call Register_SwapField_u(p,p,distrib_caldyn,
!     *                            Request_dissip)

!        call Register_SwapField_u(pk,pk,distrib_caldyn,
!     *                            Request_dissip)

!        call SendRequest(Request_dissip)       
!c$OMP BARRIER
!        call WaitRequest(Request_dissip)       

!c$OMP BARRIER
!c$OMP MASTER
!        call set_distrib(distrib_caldyn)
!        call VTe(VThallo)
!        call resume_timer(timer_caldyn)
!c        print *,'fin dissipation'
!c$OMP END MASTER
!c$OMP BARRIER
       END IF ! of IF(apdiss)

cc$OMP END PARALLEL

c ajout debug
c              IF( lafin ) then  
c                abort_message = 'Simulation finished'
c                call abort_gcm(modname,abort_message,0)
c              ENDIF
        
c   ********************************************************************
c   ********************************************************************
c   .... fin de l'integration dynamique  et physique pour le pas itau ..
c   ********************************************************************
c   ********************************************************************

c   preparation du pas d'integration suivant  ......
cym      call WriteField('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
cym      call WriteField('vcov',reshape(vcov,(/iip1,jjm,llm/)))
c$OMP MASTER      
      call stop_timer(timer_caldyn)
c$OMP END MASTER
      IF (itau==itaumax) then
c$OMP MASTER
            call allgather_timer_average
      call barrier
      if (mpi_rank==0) then
        
        print *,'*********************************'
        print *,'******    TIMER CALDYN     ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_caldyn(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_caldyn(i),timer_caldyn,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER VANLEER    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_vanleer(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_vanleer(i),timer_vanleer,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER DISSIP    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_dissip(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_dissip(i),timer_dissip,i)
        enddo
        
        print *,'*********************************'
        print *,'******    TIMER PHYSIC    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_physic(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_physic(i),timer_physic,i)
        enddo
        
      endif  
      CALL barrier
      print *,'Taille du Buffer MPI (REAL*8)',MaxBufferSize
      print *,'Taille du Buffer MPI utilise (REAL*8)',MaxBufferSize_Used
      print *, 'Temps total ecoule sur la parallelisation :',DiffTime()
      print *, 'Temps CPU ecoule sur la parallelisation :',DiffCpuTime()
      CALL print_filtre_timer
c$OMP END MASTER

!------------------------- DART-LMDZ ----------------------------------
c$OMP MASTER
      if(sub_day)then
         start_time1 = time_0 + float(sub_day_step)/float(day_step)
         if(INT(start_time1)==1) start_time1 = 0.0
      else
         start_time1 = time_0
      end if
c$OMP END MASTER
c$OMP BARRIER

!  Write natuarl variales in restart files 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,0,ps)
c$OMP BARRIER


      CALL dynredem1_loc("restart.nc",start_time1,
     .                               vcov,ucov,teta,q,masse,ps)

! Back convert natural to ucov,vcov & teta 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,1,ps)
c$OMP BARRIER
!------------------------- DART-LMDZ ----------------------------------


c$OMP MASTER
      call fin_getparam
        call finalize_parallel
c$OMP END MASTER
c$OMP BARRIER
        RETURN
      ENDIF
      
      IF ( .NOT.purmats ) THEN
c       ........................................................
c       ..............  schema matsuno + leapfrog  ..............
c       ........................................................

            IF(forward. OR. leapf) THEN
              itau= itau + 1
!              iday= day_ini+itau/day_step
!              time= REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!                IF(time.GT.1.) THEN
!                  time = time-1.
!                  iday = iday+1
!                ENDIF
            ENDIF


            IF( itau. EQ. itaufinp1 ) then

              if (flag_verif) then
                write(79,*) 'ucov',ucov
                write(80,*) 'vcov',vcov
                write(81,*) 'teta',teta
                write(82,*) 'ps',ps
                write(83,*) 'q',q
                WRITE(85,*) 'q1 = ',q(:,:,1)
                WRITE(86,*) 'q3 = ',q(:,:,3)
              endif
  

c$OMP MASTER
              call fin_getparam
              call finalize_parallel
c$OMP END MASTER
              abort_message = 'Simulation finished'
              call abort_gcm(modname,abort_message,0)
              RETURN
            ENDIF
c-----------------------------------------------------------------------
c   ecriture du fichier histoire moyenne:
c   -------------------------------------

            IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
c$OMP BARRIER
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF

#ifdef CPP_IOIPSL
             IF (ok_dynzon) THEN 

              CALL bilan_dyn_loc(2,dtvr*iperiod,dtvr*day_step*periodav, 
     ,             ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,q) 

              ENDIF !ok_dynzon

              IF (ok_dyn_ave) THEN
                 CALL writedynav_loc(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
              ENDIF
#endif
            ENDIF

c-----------------------------------------------------------------------
c   ecriture de la bande histoire:
c   ------------------------------

            IF( MOD(itau,iecri).EQ.0) THEN
             ! Ehouarn: output only during LF or Backward Matsuno
	     if (leapf.or.(.not.leapf.and.(.not.forward))) then

c$OMP BARRIER
c$OMP MASTER
              CALL geopot_loc(ip1jmp1,teta,pk,pks,phis,phi)
c$OMP END MASTER
c$OMP BARRIER
       
#ifdef CPP_IOIPSL
             if (ok_dyn_ins) then
	         CALL writehist_loc(itau,vcov,ucov,teta,phi,q,
     &                              masse,ps,phis)
             endif
#endif
            endif ! of if (leapf.or.(.not.leapf.and.(.not.forward)))
           ENDIF ! of IF(MOD(itau,iecri).EQ.0)

            IF(itau.EQ.itaufin) THEN

c$OMP BARRIER

!              if (planet_type.eq."earth") then
! Write an Earth-format restart file

!------------------------- DART-LMDZ ----------------------------------
!  Write natuarl variales in restart files 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,0,ps)
c$OMP BARRIER

c$OMP MASTER
      if(sub_day)then
         start_time1 = time_0 + float(sub_day_step)/float(day_step)
         if(INT(start_time1)==1) start_time1 = 0.0
      else
         start_time1 = time_0
      end if
c$OMP END MASTER
c$OMP BARRIER
                CALL dynredem1_loc("restart.nc",start_time1,
     &                           vcov,ucov,teta,q,masse,ps)
!              endif ! of if (planet_type.eq."earth")

! Back convert natural to ucov,vcov & teta 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,1,ps)
c$OMP BARRIER
!------------------------- DART-LMDZ ----------------------------------


!              CLOSE(99)
            ENDIF ! of IF (itau.EQ.itaufin)

c-----------------------------------------------------------------------
c   gestion de l'integration temporelle:
c   ------------------------------------

            IF( MOD(itau,iperiod).EQ.0 )    THEN
                    GO TO 1
            ELSE IF ( MOD(itau-1,iperiod). EQ. 0 ) THEN

                   IF( forward )  THEN
c      fin du pas forward et debut du pas backward

                      forward = .FALSE.
                        leapf = .FALSE.
                           GO TO 2

                   ELSE
c      fin du pas backward et debut du premier pas leapfrog

                        leapf =  .TRUE.
                        dt  =  2.*dtvr
                        GO TO 2
                   END IF
            ELSE

c      ......   pas leapfrog  .....

                 leapf = .TRUE.
                 dt  = 2.*dtvr
                 GO TO 2
            END IF ! of IF (MOD(itau,iperiod).EQ.0)
                   !    ELSEIF (MOD(itau-1,iperiod).EQ.0)


      ELSE ! of IF (.not.purmats)

c       ........................................................
c       ..............       schema  matsuno        ...............
c       ........................................................
            IF( forward )  THEN

             itau =  itau + 1
!             iday = day_ini+itau/day_step
!             time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!
!                  IF(time.GT.1.) THEN
!                   time = time-1.
!                   iday = iday+1
!                  ENDIF

               forward =  .FALSE.
               IF( itau. EQ. itaufinp1 ) then  
c$OMP MASTER
                 call fin_getparam
                 call finalize_parallel
c$OMP END MASTER
                 abort_message = 'Simulation finished'
                 call abort_gcm(modname,abort_message,0)
                 RETURN
               ENDIF
               GO TO 2

            ELSE ! of IF(forward) i.e. backward step

              IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF

#ifdef CPP_IOIPSL
               IF (ok_dynzon) THEN
               CALL bilan_dyn_loc(2,dtvr*iperiod,dtvr*day_step*periodav,
     ,           ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,q)
               ENDIF
              
               IF (ok_dyn_ave) THEN
                 CALL writedynav_loc(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
               ENDIF
#endif
              ENDIF ! of IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin)


               IF(MOD(itau,iecri         ).EQ.0) THEN

c$OMP BARRIER
c$OMP MASTER
              CALL geopot_loc(ip1jmp1,teta,pk,pks,phis,phi)
c$OMP END MASTER
c$OMP BARRIER


#ifdef CPP_IOIPSL
              if (ok_dyn_ins) then
                 CALL writehist_loc(itau,vcov,ucov,teta,phi,q,
     &                              masse,ps,phis)
              endif ! of if (ok_dyn_ins)
#endif
              ENDIF ! of IF(MOD(itau,iecri).EQ.0)
              

              IF(itau.EQ.itaufin) THEN
!                if (planet_type.eq."earth") then

!------------------------- DART-LMDZ ----------------------------------
!  Write natuarl variales in restart files 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,0,ps)
c$OMP BARRIER

c$OMP MASTER
      if(sub_day)then
         start_time1 = time_0 + float(sub_day_step)/float(day_step)
         if(INT(start_time1)==1) start_time1 = 0.0
      else
         start_time1 = time_0
      end if
c$OMP END MASTER
c$OMP BARRIER
                   CALL dynredem1_loc("restart.nc",start_time1,
     .                               vcov,ucov,teta,q,masse,ps)
!               endif ! of if (planet_type.eq."earth")

! Back convert natural to ucov,vcov & teta 
      if(write_restart_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,1,ps)
!------------------------- DART-LMDZ ----------------------------------

              ENDIF ! of IF(itau.EQ.itaufin)

              forward = .TRUE.
              GO TO  1

            ENDIF ! of IF (forward)

      END IF ! of IF(.not.purmats)
c$OMP MASTER
      call fin_getparam
      call finalize_parallel
c$OMP END MASTER
      abort_message = 'Simulation finished'
      call abort_gcm(modname,abort_message,0)
      RETURN
      END
