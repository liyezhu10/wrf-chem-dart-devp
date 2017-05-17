!
! $Id$
!
c
c
      PROGRAM gcm

#ifdef CPP_IOIPSL
      USE IOIPSL
#endif

      USE mod_const_mpi, ONLY: init_const_mpi
      USE parallel_lmdz
      USE infotrac
      USE mod_interface_dyn_phys
      USE mod_hallo
      USE Bands
      USE getparam
      USE filtreg_mod
      USE control_mod

#ifdef INCA
! Only INCA needs these informations (from the Earth's physics)
      USE indice_sol_mod
#endif

#ifdef CPP_PHYS
      USE mod_grid_phy_lmdz
      USE mod_phys_lmdz_para, ONLY : klon_mpi_para_nb
      USE mod_phys_lmdz_omp_data, ONLY: klon_omp 
      USE dimphy
      USE comgeomphy
#endif
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
c  ... Possibilite de choisir le schema pour l'advection de
c        q  , en modifiant iadv dans traceur.def  (MAF,10/02) .
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
#include "tracstoke.h"

#ifdef INCA
! Only INCA needs these informations (from the Earth's physics)
!#include "indicesol.h"
#endif

      INTEGER         longcles
      PARAMETER     ( longcles = 20 )
      REAL  clesphy0( longcles )
      SAVE  clesphy0



      REAL zdtvr

c   variables dynamiques
      REAL,ALLOCATABLE,SAVE  :: vcov(:,:),ucov(:,:) ! vents covariants
      REAL,ALLOCATABLE,SAVE  :: teta(:,:)     ! temperature potentielle 
      REAL, ALLOCATABLE,SAVE :: q(:,:,:)      ! champs advectes
      REAL,ALLOCATABLE,SAVE  :: ps(:)         ! pression  au sol
c      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
      REAL,ALLOCATABLE,SAVE  :: masse(:,:)    ! masse d'air
      REAL,ALLOCATABLE,SAVE  :: phis(:)       ! geopotentiel au sol
c      REAL phi(ip1jmp1,llm)                  ! geopotentiel
c      REAL w(ip1jmp1,llm)                    ! vitesse verticale

c variables dynamiques intermediaire pour le transport

c   variables pour le fichier histoire
      REAL dtav      ! intervalle de temps elementaire

      REAL time_0

      LOGICAL lafin
c      INTEGER ij,iq,l,i,j
      INTEGER i,j


      real time_step, t_wrt, t_ops


      LOGICAL call_iniphys
      data call_iniphys/.true./

c+jld variables test conservation energie
c      REAL ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
C     Tendance de la temp. potentiel d (theta)/ d t due a la 
C     tansformation d'energie cinetique en energie thermique
C     cree par la dissipation
c      REAL dhecdt(ip1jmp1,llm)
c      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
c      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
c      CHARACTER (len=15) :: ztit
c-jld 


      character (len=80) :: dynhist_file, dynhistave_file
      character (len=20) :: modname
      character (len=80) :: abort_message
! locales pour gestion du temps
      INTEGER :: an, mois, jour
      REAL :: heure


c-----------------------------------------------------------------------
c    variables pour l'initialisation de la physique :
c    ------------------------------------------------
      INTEGER ngridmx
      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )
      REAL zcufi(ngridmx),zcvfi(ngridmx)
      REAL latfi(ngridmx),lonfi(ngridmx)
      REAL airefi(ngridmx)
      SAVE latfi, lonfi, airefi
      
      INTEGER :: ierr


c-----------------------------------------------------------------------
c   Initialisations:
c   ----------------

      abort_message = 'last timestep reached'
      modname = 'gcm'
      descript = 'Run GCM LMDZ'
      lafin    = .FALSE.
      dynhist_file = 'dyn_hist'
      dynhistave_file = 'dyn_hist_ave'



c----------------------------------------------------------------------
c  lecture des fichiers gcm.def ou run.def
c  ---------------------------------------
c
! Ehouarn: dump possibility of using defrun
!#ifdef CPP_IOIPSL
      CALL conf_gcm( 99, .TRUE. , clesphy0 )
      if (mod(iphysiq, iperiod) /= 0) call abort_gcm("conf_gcm",
     s "iphysiq must be a multiple of iperiod", 1)
!#else
!      CALL defrun( 99, .TRUE. , clesphy0 )
!#endif
c
c
c------------------------------------
c   Initialisation partie parallele
c------------------------------------
      CALL init_const_mpi

      call init_parallel
      call ini_getparam("out.def")
      call Read_Distrib

#ifdef CPP_PHYS
        CALL Init_Phys_lmdz(iim,jjp1,llm,mpi_size,distrib_phys)
#endif
      CALL set_bands
#ifdef CPP_PHYS
      CALL Init_interface_dyn_phys
#endif
      CALL barrier

      if (mpi_rank==0) call WriteBands
      call Set_Distrib(distrib_caldyn)

c$OMP PARALLEL
      call Init_Mod_hallo
c$OMP END PARALLEL

#ifdef CPP_PHYS
c$OMP PARALLEL
      call InitComgeomphy
c$OMP END PARALLEL 
#endif

c-----------------------------------------------------------------------
c   Choix du calendrier
c   -------------------

c      calend = 'earth_365d'

#ifdef CPP_IOIPSL
      if (calend == 'earth_360d') then
        call ioconf_calendar('360d')
        write(lunout,*)'CALENDRIER CHOISI: Terrestre a 360 jours/an'
      else if (calend == 'earth_365d') then
        call ioconf_calendar('noleap')
        write(lunout,*)'CALENDRIER CHOISI: Terrestre a 365 jours/an'
      else if (calend == 'earth_366d') then
        call ioconf_calendar('gregorian')
        write(lunout,*)'CALENDRIER CHOISI: Terrestre bissextile'
      else
        abort_message = 'Mauvais choix de calendrier'
        call abort_gcm(modname,abort_message,1)
      endif
#endif

      IF (type_trac == 'inca') THEN
#ifdef INCA
         call init_const_lmdz(
     $        nbtr,anneeref,dayref,
     $        iphysiq,day_step,nday, 
     $        nbsrf, is_oce,is_sic,
     $        is_ter,is_lic)

         call init_inca_para(
     $        iim,jjm+1,llm,klon_glo,mpi_size,
     $        distrib_phys,COMM_LMDZ)
#endif
      END IF

c-----------------------------------------------------------------------
c   Initialisation des traceurs
c   ---------------------------
c  Choix du nombre de traceurs et du schema pour l'advection
c  dans fichier traceur.def, par default ou via INCA
      call infotrac_init

c Allocation de la tableau q : champs advectes   
      ALLOCATE(ucov(ijb_u:ije_u,llm))
      ALLOCATE(vcov(ijb_v:ije_v,llm))
      ALLOCATE(teta(ijb_u:ije_u,llm))
      ALLOCATE(masse(ijb_u:ije_u,llm))
      ALLOCATE(ps(ijb_u:ije_u))
      ALLOCATE(phis(ijb_u:ije_u))
      ALLOCATE(q(ijb_u:ije_u,llm,nqtot))

c-----------------------------------------------------------------------
c   Lecture de l'etat initial :
c   ---------------------------

c  lecture du fichier start.nc
      if (read_start) then
      ! we still need to run iniacademic to initialize some
      ! constants & fields, if we run the 'newtonian' or 'SW' cases:
        if (iflag_phys.ne.1) then
          CALL iniacademic_loc(vcov,ucov,teta,q,masse,ps,phis,time_0)
        endif

!        if (planet_type.eq."earth") then
! Load an Earth-format start file
         CALL dynetat0_loc("start.nc",vcov,ucov,
     &              teta,q,masse,ps,phis, time_0)
!        endif ! of if (planet_type.eq."earth")

c       write(73,*) 'ucov',ucov
c       write(74,*) 'vcov',vcov
c       write(75,*) 'teta',teta
c       write(76,*) 'ps',ps
c       write(77,*) 'q',q

      endif ! of if (read_start)

c le cas echeant, creation d un etat initial
      IF (prt_level > 9) WRITE(lunout,*)
     .              'GCM: AVANT iniacademic AVANT AVANT AVANT AVANT'
      if (.not.read_start) then
         CALL iniacademic_loc(vcov,ucov,teta,q,masse,ps,phis,time_0)
      endif

c-----------------------------------------------------------------------
c   Lecture des parametres de controle pour la simulation :
c   -------------------------------------------------------
c  on recalcule eventuellement le pas de temps

      IF(MOD(day_step,iperiod).NE.0) THEN
        abort_message = 
     .  'Il faut choisir un nb de pas par jour multiple de iperiod'
        call abort_gcm(modname,abort_message,1)
      ENDIF

      IF(MOD(day_step,iphysiq).NE.0) THEN
        abort_message = 
     * 'Il faut choisir un nb de pas par jour multiple de iphysiq'
        call abort_gcm(modname,abort_message,1)
      ENDIF

      zdtvr    = daysec/REAL(day_step)
        IF(dtvr.NE.zdtvr) THEN
         WRITE(lunout,*)
     .    'WARNING!!! changement de pas de temps',dtvr,'>',zdtvr
        ENDIF

C
C on remet le calendrier à zero si demande
c
      IF (start_time /= starttime) then
        WRITE(lunout,*)' GCM: Attention l''heure de depart lue dans le'
     &,' fichier restart ne correspond pas à celle lue dans le run.def'
        IF (raz_date == 1) then
          WRITE(lunout,*)'Je prends l''heure lue dans run.def'
          start_time = starttime
        ELSE
          WRITE(lunout,*)'Je m''arrete'
          CALL abort
        ENDIF
      ENDIF
      IF (raz_date == 1) THEN
        annee_ref = anneeref
        day_ref = dayref
        day_ini = dayref
        itau_dyn = 0
        itau_phy = 0
        time_0 = 0.
        write(lunout,*)
     .   'GCM: On reinitialise a la date lue dans gcm.def'
      ELSE IF (annee_ref .ne. anneeref .or. day_ref .ne. dayref) THEN
        write(lunout,*)
     .  'GCM: Attention les dates initiales lues dans le fichier'
        write(lunout,*)
     .  ' restart ne correspondent pas a celles lues dans '
        write(lunout,*)' gcm.def'
        write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
        write(lunout,*)' day_ref=',day_ref," dayref=",dayref
        write(lunout,*)' Pas de remise a zero'
      ENDIF
c      if (annee_ref .ne. anneeref .or. day_ref .ne. dayref) then
c        write(lunout,*)
c     .  'GCM: Attention les dates initiales lues dans le fichier'
c        write(lunout,*)
c     .  ' restart ne correspondent pas a celles lues dans '
c        write(lunout,*)' gcm.def'
c        write(lunout,*)' annee_ref=',annee_ref," anneeref=",anneeref
c        write(lunout,*)' day_ref=',day_ref," dayref=",dayref
c        if (raz_date .ne. 1) then
c          write(lunout,*)
c     .    'GCM: On garde les dates du fichier restart'
c        else
c          annee_ref = anneeref
c          day_ref = dayref
c          day_ini = dayref
c          itau_dyn = 0
c          itau_phy = 0
c          time_0 = 0.
c          write(lunout,*)
c     .   'GCM: On reinitialise a la date lue dans gcm.def'
c        endif
c      ELSE
c        raz_date = 0
c      endif

#ifdef CPP_IOIPSL
      mois = 1
      heure = 0.
      call ymds2ju(annee_ref, mois, day_ref, heure, jD_ref)
      jH_ref = jD_ref - int(jD_ref)
      jD_ref = int(jD_ref)

      call ioconf_startdate(INT(jD_ref), jH_ref)

      write(lunout,*)'DEBUG'
      write(lunout,*)'annee_ref, mois, day_ref, heure, jD_ref'
      write(lunout,*)annee_ref, mois, day_ref, heure, jD_ref
      call ju2ymds(jD_ref+jH_ref,an, mois, jour, heure)
      write(lunout,*)'jD_ref+jH_ref,an, mois, jour, heure'
      write(lunout,*)jD_ref+jH_ref,an, mois, jour, heure
#else
! Ehouarn: we still need to define JD_ref and JH_ref
! and since we don't know how many days there are in a year
! we set JD_ref to 0 (this should be improved ...)
      jD_ref=0
      jH_ref=0
#endif

      if (iflag_phys.eq.1) then
      ! these initialisations have already been done (via iniacademic)
      ! if running in SW or Newtonian mode
c-----------------------------------------------------------------------
c   Initialisation des constantes dynamiques :
c   ------------------------------------------
        dtvr = zdtvr
        CALL iniconst

c-----------------------------------------------------------------------
c   Initialisation de la geometrie :
c   --------------------------------
        CALL inigeom

c-----------------------------------------------------------------------
c   Initialisation du filtre :
c   --------------------------
        CALL inifilr
      endif ! of if (iflag_phys.eq.1)
c
c-----------------------------------------------------------------------
c   Initialisation de la dissipation :
c   ----------------------------------

      CALL inidissip( lstardis, nitergdiv, nitergrot, niterh   ,
     *                tetagdiv, tetagrot , tetatemp, vert_prof_dissip)

c-----------------------------------------------------------------------
c   Initialisation de la physique :
c   -------------------------------
      IF (call_iniphys.and.(iflag_phys==1.or.iflag_phys>=100)) THEN
         latfi(1)=rlatu(1)
         lonfi(1)=0.
         zcufi(1) = cu(1)
         zcvfi(1) = cv(1)
         DO j=2,jjm
            DO i=1,iim
               latfi((j-2)*iim+1+i)= rlatu(j)
               lonfi((j-2)*iim+1+i)= rlonv(i)
               zcufi((j-2)*iim+1+i) = cu((j-1)*iip1+i)
               zcvfi((j-2)*iim+1+i) = cv((j-1)*iip1+i)
            ENDDO
         ENDDO
         latfi(ngridmx)= rlatu(jjp1)
         lonfi(ngridmx)= 0.
         zcufi(ngridmx) = cu(ip1jm+1)
         zcvfi(ngridmx) = cv(ip1jm-iim)
         CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,aire,airefi)

         WRITE(lunout,*)
     .       'GCM: WARNING!!! vitesse verticale nulle dans la physique'
! Physics:
#ifdef CPP_PHYS
         CALL iniphysiq(ngridmx,llm,daysec,day_ini,dtphys/nsplit_phys, 
     &                latfi,lonfi,airefi,zcufi,zcvfi,rad,g,r,cpp, 
     &                iflag_phys) 
#endif
         call_iniphys=.false.
      ENDIF ! of IF (call_iniphys.and.(iflag_phys==1.or.iflag_phys>=100))


c-----------------------------------------------------------------------
c   Initialisation des dimensions d'INCA :
c   --------------------------------------
      IF (type_trac == 'inca') THEN
!$OMP PARALLEL
#ifdef INCA
         CALL init_inca_dim(klon_omp,llm,iim,jjm,
     $        rlonu,rlatu,rlonv,rlatv)
#endif
!$OMP END PARALLEL
      END IF

c-----------------------------------------------------------------------
c   Initialisation des I/O :
c   ------------------------


      if (nday>=0) then
         day_end = day_ini + nday
      else
         day_end = day_ini - nday/day_step
      endif

      !! ----- DART-LMDZ -----------------------------------
      if(sub_day) then
       day_end = day_ini+INT(time_0+float(sub_day_step)/float(day_step))
      endif
      !! ----- DART-LMDZ -----------------------------------

  
      WRITE(lunout,300)day_ini,day_end
 300  FORMAT('1'/,15x,'run du jour',i7,2x,'au jour',i7//)

#ifdef CPP_IOIPSL
      call ju2ymds(jD_ref + day_ini - day_ref, an, mois, jour, heure)
      write (lunout,301)jour, mois, an
      call ju2ymds(jD_ref + day_end - day_ref, an, mois, jour, heure)
      write (lunout,302)jour, mois, an
 301  FORMAT('1'/,15x,'run du ', i2,'/',i2,'/',i4)
 302  FORMAT('1'/,15x,'    au ', i2,'/',i2,'/',i4)
#endif

!      if (planet_type.eq."earth") then
! Write an Earth-format restart file
        CALL dynredem0_loc("restart.nc", day_end, phis)
!      endif

      ecripar = .TRUE.

#ifdef CPP_IOIPSL
      time_step = zdtvr
      IF (mpi_rank==0) then
        if (ok_dyn_ins) then
          ! initialize output file for instantaneous outputs
          ! t_ops = iecri * daysec ! do operations every t_ops
          t_ops =((1.0*iecri)/day_step) * daysec  
          t_wrt = daysec ! iecri * daysec ! write output every t_wrt
          t_wrt = daysec ! iecri * daysec ! write output every t_wrt
          CALL inithist(day_ref,annee_ref,time_step,
     &                  t_ops,t_wrt)
        endif

      IF (ok_dyn_ave) THEN 
         ! initialize output file for averaged outputs
         t_ops = iperiod * time_step ! do operations every t_ops
         t_wrt = periodav * daysec   ! write output every t_wrt
         CALL initdynav_loc(day_ref,annee_ref,time_step,t_ops,t_wrt)
        END IF
      ENDIF
      dtav = iperiod*dtvr/daysec
#endif
! #endif of #ifdef CPP_IOIPSL

c  Choix des frequences de stokage pour le offline
c      istdyn=day_step/4     ! stockage toutes les 6h=1jour/4
c      istdyn=day_step/12     ! stockage toutes les 2h=1jour/12
      istdyn=day_step/4     ! stockage toutes les 6h=1jour/12
      istphy=istdyn/iphysiq     


c
c-----------------------------------------------------------------------
c   Integration temporelle du modele :
c   ----------------------------------

c       write(78,*) 'ucov',ucov
c       write(78,*) 'vcov',vcov
c       write(78,*) 'teta',teta
c       write(78,*) 'ps',ps
c       write(78,*) 'q',q

c$OMP PARALLEL DEFAULT(SHARED) COPYIN(/temps/,/logici/,/logicl/)

!----------- DART-LMDZ   ------------------ 
       if(start_from_natural) call naturel_loc(ucov,vcov,teta,
     &                                         ucov,vcov,teta,1,ps)
!----------- DART-LMDZ   ------------------ 

!$OMP BARRIER
!___________________
      CALL leapfrog_loc(ucov,vcov,teta,ps,masse,phis,q,clesphy0,
     .              time_0)
c$OMP END PARALLEL

!      OPEN(unit=5487,file='ok_lmdz',status='replace')
!      WRITE(5487,*) 'ok_lmdz'
!      CLOSE(5487)
      END

