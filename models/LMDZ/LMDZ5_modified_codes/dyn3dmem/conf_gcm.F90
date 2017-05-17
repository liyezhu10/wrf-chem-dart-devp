
! $Id$

SUBROUTINE conf_gcm( tapedef, etatinit, clesphy0 )

  USE control_mod
#ifdef CPP_IOIPSL
  use IOIPSL
#else
  ! if not using IOIPSL, we still need to use (a local version of) getin
  use ioipsl_getincom
#endif
  use misc_mod
  use mod_filtre_fft, ONLY : use_filtre_fft
  use mod_filtre_fft_loc, ONLY : use_filtre_fft_loc=>use_filtre_fft
  use mod_hallo, ONLY : use_mpi_alloc
  USE infotrac, ONLY : type_trac
  use assert_m, only: assert

  IMPLICIT NONE
  !-----------------------------------------------------------------------
  !     Auteurs :   L. Fairhead , P. Le Van  .

  !     Arguments :

  !     tapedef   :
  !     etatinit  :     = TRUE   , on ne  compare pas les valeurs des para- 
  !     -metres  du zoom  avec  celles lues sur le fichier start .
  !      clesphy0 :  sortie  .

  LOGICAL etatinit
  INTEGER tapedef

  INTEGER        longcles
  PARAMETER(     longcles = 20 )
  REAL clesphy0( longcles )

  !   Declarations :
  !   --------------
  include "dimensions.h"
  include "paramet.h"
  include "logic.h"
  include "serre.h"
  include "comdissnew.h"
  include "temps.h"
  include "comconst.h"

  ! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
  ! include "clesphys.h"
  include "iniprint.h"

  !   local:
  !   ------

  CHARACTER ch1*72,ch2*72,ch3*72,ch4*12
  REAL clonn,clatt,grossismxx,grossismyy
  REAL dzoomxx,dzoomyy, tauxx,tauyy
  LOGICAL  fxyhypbb, ysinuss
  INTEGER i
  character(len=*),parameter :: modname="conf_gcm"
  character (len=80) :: abort_message
#ifdef CPP_OMP
  integer,external :: OMP_GET_NUM_THREADS
#endif

  !  -------------------------------------------------------------------

  !       .........     Version  du 29/04/97       ..........

  !   Nouveaux parametres nitergdiv,nitergrot,niterh,tetagdiv,tetagrot,
  !      tetatemp   ajoutes  pour la dissipation   .

  !   Autre parametre ajoute en fin de liste de tapedef : ** fxyhypb ** 

  !  Si fxyhypb = .TRUE. , choix de la fonction a derivee tangente hyperb.
  !    Sinon , choix de fxynew  , a derivee sinusoidale  ..

  !   ......  etatinit = . TRUE. si defrun  est appele dans ETAT0_LMD  ou
  !         LIMIT_LMD  pour l'initialisation de start.dat (dic) et
  !                de limit.dat ( dic)                        ...........
  !           Sinon  etatinit = . FALSE .

  !   Donc etatinit = .F.  si on veut comparer les valeurs de  grossismx ,
  !    grossismy,clon,clat, fxyhypb  lues sur  le fichier  start  avec
  !   celles passees  par run.def ,  au debut du gcm, apres l'appel a 
  !    lectba .  
  !   Ces parmetres definissant entre autres la grille et doivent etre
  !   pareils et coherents , sinon il y aura  divergence du gcm .

  !-----------------------------------------------------------------------
  !   initialisations:
  !   ----------------

  !Config  Key  = lunout
  !Config  Desc = unite de fichier pour les impressions
  !Config  Def  = 6
  !Config  Help = unite de fichier pour les impressions 
  !Config         (defaut sortie standard = 6)
  lunout=6
  CALL getin('lunout', lunout)
  IF (lunout /= 5 .and. lunout /= 6) THEN
     OPEN(UNIT=lunout,FILE='lmdz.out_0000',ACTION='write',  &
          STATUS='unknown',FORM='formatted')
  ENDIF

  adjust=.false.
  call getin('adjust',adjust)

#ifdef CPP_OMP
  ! adjust=y not implemented in case of OpenMP threads...
  !$OMP PARALLEL
  if ((OMP_GET_NUM_THREADS()>1).and.adjust) then
     write(lunout,*)'conf_gcm: Error, adjust should be set to n' &
          ,' when running with OpenMP threads'
     abort_message = 'Wrong value for adjust'
     call abort_gcm(modname,abort_message,1)
  endif
  !$OMP END PARALLEL         
#endif

  itaumax=0
  call getin('itaumax',itaumax);
  if (itaumax<=0) itaumax=HUGE(itaumax)

  !Config  Key  = prt_level
  !Config  Desc = niveau d'impressions de débogage
  !Config  Def  = 0
  !Config  Help = Niveau d'impression pour le débogage
  !Config         (0 = minimum d'impression)
  prt_level = 0
  CALL getin('prt_level',prt_level)

  !-----------------------------------------------------------------------
  !  Parametres de controle du run:
  !-----------------------------------------------------------------------
  !Config  Key  = planet_type
  !Config  Desc = planet type ("earth", "mars", "venus", ...)
  !Config  Def  = earth
  !Config  Help = this flag sets the type of atymosphere that is considered
  planet_type="earth"
  CALL getin('planet_type',planet_type)

  !Config  Key  = calend
  !Config  Desc = type de calendrier utilise
  !Config  Def  = earth_360d
  !Config  Help = valeur possible: earth_360d, earth_365d, earth_366d
  !Config         
  calend = 'earth_360d'
  CALL getin('calend', calend)

  !Config  Key  = dayref
  !Config  Desc = Jour de l'etat initial
  !Config  Def  = 1
  !Config  Help = Jour de l'etat initial ( = 350  si 20 Decembre ,
  !Config         par expl. ,comme ici ) ... A completer
  dayref=1
  CALL getin('dayref', dayref)

  !Config  Key  = anneeref
  !Config  Desc = Annee de l'etat initial
  !Config  Def  = 1998
  !Config  Help = Annee de l'etat  initial 
  !Config         (   avec  4  chiffres   ) ... A completer
  anneeref = 1998
  CALL getin('anneeref',anneeref)

  !Config  Key  = raz_date
  !Config  Desc = Remise a zero de la date initiale
  !Config  Def  = 0 (pas de remise a zero)
  !Config  Help = Remise a zero de la date initiale 
  !Config         0 pas de remise a zero, on garde la date du fichier restart
  !Config         1 prise en compte de la date de gcm.def avec remise a zero
  !Config         des compteurs de pas de temps
  raz_date = 0
  CALL getin('raz_date', raz_date)

  !Config  Key  = resetvarc
  !Config  Desc = Reinit des variables de controle
  !Config  Def  = n
  !Config  Help = Reinit des variables de controle
  resetvarc = .false.
  CALL getin('resetvarc',resetvarc)

  !Config  Key  = nday
  !Config  Desc = Nombre de jours d'integration
  !Config  Def  = 10
  !Config  Help = Nombre de jours d'integration
  !Config         ... On pourait aussi permettre des mois ou des annees !
  nday = 10
  CALL getin('nday',nday)

  !Config  Key  = starttime
  !Config  Desc = Heure de depart de la simulation
  !Config  Def  = 0
  !Config  Help = Heure de depart de la simulation
  !Config         en jour
  starttime = 0
  CALL getin('starttime',starttime)

  !Config  Key  = day_step
  !Config  Desc = nombre de pas par jour
  !Config  Def  = 240 
  !Config  Help = nombre de pas par jour (multiple de iperiod) (
  !Config          ici pour  dt = 1 min ) 
  day_step = 240 
  CALL getin('day_step',day_step)

  !Config  Key  = nsplit_phys
  nsplit_phys = 1 
  CALL getin('nsplit_phys',nsplit_phys)

  !Config  Key  = iperiod
  !Config  Desc = periode pour le pas Matsuno
  !Config  Def  = 5
  !Config  Help = periode pour le pas Matsuno (en pas de temps)
  iperiod = 5
  CALL getin('iperiod',iperiod)

  !Config  Key  = iapp_tracvl
  !Config  Desc = frequence du groupement des flux 
  !Config  Def  = iperiod
  !Config  Help = frequence du groupement des flux (en pas de temps) 
  iapp_tracvl = iperiod
  CALL getin('iapp_tracvl',iapp_tracvl)

  !Config  Key  = iconser
  !Config  Desc = periode de sortie des variables de controle
  !Config  Def  = 240  
  !Config  Help = periode de sortie des variables de controle
  !Config         (En pas de temps)
  iconser = 240  
  CALL getin('iconser', iconser)

  !Config  Key  = iecri
  !Config  Desc = periode d'ecriture du fichier histoire
  !Config  Def  = 1
  !Config  Help = periode d'ecriture du fichier histoire (en jour) 
  iecri = 1
  CALL getin('iecri',iecri)

  !Config  Key  = periodav
  !Config  Desc = periode de stockage fichier histmoy
  !Config  Def  = 1
  !Config  Help = periode de stockage fichier histmoy (en jour) 
  periodav = 1.
  CALL getin('periodav',periodav)

  !Config  Key  = output_grads_dyn
  !Config  Desc = output dynamics diagnostics in 'dyn.dat' file
  !Config  Def  = n
  !Config  Help = output dynamics diagnostics in Grads-readable 'dyn.dat' file
  output_grads_dyn=.false.
  CALL getin('output_grads_dyn',output_grads_dyn)

  !Config  Key  = dissip_period
  !Config  Desc = periode de la dissipation 
  !Config  Def  = 0
  !Config  Help = periode de la dissipation 
  !Config  dissip_period=0 => la valeur sera calcule dans inidissip       
  !Config  dissip_period>0 => on prend cette valeur
  dissip_period = 0
  CALL getin('dissip_period',dissip_period)

  !cc  ....   P. Le Van , modif le 29/04/97 .pour la dissipation  ...
  !cc

  !Config  Key  = lstardis
  !Config  Desc = choix de l'operateur de dissipation
  !Config  Def  = y
  !Config  Help = choix de l'operateur de dissipation
  !Config         'y' si on veut star et 'n' si on veut non-start !
  !Config         Moi y en a pas comprendre ! 
  lstardis = .TRUE.
  CALL getin('lstardis',lstardis)

  !Config  Key  = nitergdiv
  !Config  Desc = Nombre d'iteration de gradiv
  !Config  Def  = 1
  !Config  Help = nombre d'iterations de l'operateur de dissipation 
  !Config         gradiv
  nitergdiv = 1
  CALL getin('nitergdiv',nitergdiv)

  !Config  Key  = nitergrot
  !Config  Desc = nombre d'iterations de nxgradrot
  !Config  Def  = 2
  !Config  Help = nombre d'iterations de l'operateur de dissipation  
  !Config         nxgradrot
  nitergrot = 2
  CALL getin('nitergrot',nitergrot)

  !Config  Key  = niterh
  !Config  Desc = nombre d'iterations de divgrad
  !Config  Def  = 2
  !Config  Help = nombre d'iterations de l'operateur de dissipation
  !Config         divgrad
  niterh = 2
  CALL getin('niterh',niterh)

  !Config  Key  = tetagdiv
  !Config  Desc = temps de dissipation pour div
  !Config  Def  = 7200
  !Config  Help = temps de dissipation des plus petites longeur 
  !Config         d'ondes pour u,v (gradiv)
  tetagdiv = 7200.
  CALL getin('tetagdiv',tetagdiv)

  !Config  Key  = tetagrot
  !Config  Desc = temps de dissipation pour grad
  !Config  Def  = 7200
  !Config  Help = temps de dissipation des plus petites longeur 
  !Config         d'ondes pour u,v (nxgradrot)
  tetagrot = 7200.
  CALL getin('tetagrot',tetagrot)

  !Config  Key  = tetatemp 
  !Config  Desc = temps de dissipation pour h
  !Config  Def  = 7200
  !Config  Help =  temps de dissipation des plus petites longeur 
  !Config         d'ondes pour h (divgrad)   
  tetatemp  = 7200.
  CALL getin('tetatemp',tetatemp )

  ! Parametres controlant la variation sur la verticale des constantes de
  ! dissipation.
  ! Pour le moment actifs uniquement dans la version a 39 niveaux
  ! avec ok_strato=y

  dissip_factz=4.
  dissip_deltaz=10.
  dissip_zref=30.
  CALL getin('dissip_factz',dissip_factz )
  CALL getin('dissip_deltaz',dissip_deltaz )
  CALL getin('dissip_zref',dissip_zref )

  ! top_bound sponge: only active if ok_strato=.true. and iflag_top_bound!=0
  !                   iflag_top_bound=0 for no sponge
  !                   iflag_top_bound=1 for sponge over 4 topmost layers
  !                   iflag_top_bound=2 for sponge from top to ~1% of top layer pressure
  iflag_top_bound=1
  CALL getin('iflag_top_bound',iflag_top_bound)

  ! mode_top_bound : fields towards which sponge relaxation will be done:
  !                  mode_top_bound=0: no relaxation
  !                  mode_top_bound=1: u and v relax towards 0
  !                  mode_top_bound=2: u and v relax towards their zonal mean
  !                  mode_top_bound=3: u,v and pot. temp. relax towards their zonal mean
  mode_top_bound=3
  CALL getin('mode_top_bound',mode_top_bound)

  ! top_bound sponge : inverse of charactericstic relaxation time scale for sponge
  tau_top_bound=1.e-5
  CALL getin('tau_top_bound',tau_top_bound)

  !Config  Key  = coefdis
  !Config  Desc = coefficient pour gamdissip
  !Config  Def  = 0
  !Config  Help = coefficient pour gamdissip  
  coefdis = 0.
  CALL getin('coefdis',coefdis)

  !Config  Key  = purmats
  !Config  Desc = Schema d'integration
  !Config  Def  = n
  !Config  Help = Choix du schema d'integration temporel.
  !Config         y = pure Matsuno sinon c'est du Matsuno-leapfrog
  purmats = .FALSE.
  CALL getin('purmats',purmats)

  !Config  Key  = ok_guide
  !Config  Desc = Guidage
  !Config  Def  = n
  !Config  Help = Guidage
  ok_guide = .FALSE.
  CALL getin('ok_guide',ok_guide)

  if (ok_guide .and. adjust) call abort_gcm("conf_gcm", &
       "adjust does not work with ok_guide", 1)

  !Config  Key  =  read_start
  !Config  Desc = Initialize model using a 'start.nc' file
  !Config  Def  = y
  !Config  Help = y: intialize dynamical fields using a 'start.nc' file
  !               n: fields are initialized by 'iniacademic' routine
  read_start= .true.
  CALL getin('read_start',read_start)

  !Config  Key  = iflag_phys
  !Config  Desc = Avec ls physique 
  !Config  Def  = 1
  !Config  Help = Permet de faire tourner le modele sans 
  !Config         physique.
  iflag_phys = 1
  CALL getin('iflag_phys',iflag_phys)

  !Config  Key  =  iphysiq
  !Config  Desc = Periode de la physique
  !Config  Def  = 5
  !Config  Help = Periode de la physique en pas de temps de la dynamique.
  iphysiq = 5
  CALL getin('iphysiq', iphysiq)

  !Config  Key  = ip_ebil_dyn
  !Config  Desc = PRINT level for energy conserv. diag.
  !Config  Def  = 0
  !Config  Help = PRINT level for energy conservation diag. ;
  !               les options suivantes existent :
  !Config         0 pas de print
  !Config         1 pas de print
  !Config         2 print,
  ip_ebil_dyn = 0
  CALL getin('ip_ebil_dyn',ip_ebil_dyn)
!------------ DART-LMDZ ---------------------------
       start_from_natural = .FALSE.
       CALL getin('start_from_natural',start_from_natural)

       write_restart_natural = .FALSE.
       CALL getin('write_restart_natural',write_restart_natural)

       sub_day = .FALSE.
       CALL getin('sub_day',sub_day)

       sub_day_step = 120
       CALL getin('sub_day_step',sub_day_step)

       assimilation = .FALSE.
       CALL getin('assimilation',assimilation)

       print_sortvar = .FALSE.
       CALL getin('print_sortvar',print_sortvar)
  
!------------ DART-LMDZ -----------------

  !cc  ....   P. Le Van , ajout  le 7/03/95 .pour le zoom ...
  !     .........   (  modif  le 17/04/96 )   .........

  test_etatinit: IF (.not. etatinit) then
     !Config  Key  = clon
     !Config  Desc = centre du zoom, longitude
     !Config  Def  = 0
     !Config  Help = longitude en degres du centre 
     !Config         du zoom
     clonn = 0.
     CALL getin('clon',clonn)

     !Config  Key  = clat
     !Config  Desc = centre du zoom, latitude
     !Config  Def  = 0
     !Config  Help = latitude en degres du centre du zoom
     !Config         
     clatt = 0.
     CALL getin('clat',clatt)

     IF( ABS(clat - clatt).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de clat passee par run.def', &
             ' est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     !Config  Key  = grossismx 
     !Config  Desc = zoom en longitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la longitude
     grossismxx = 1.0
     CALL getin('grossismx',grossismxx)

     IF( ABS(grossismx - grossismxx).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de grossismx passee par ', &
             'run.def est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     !Config  Key  = grossismy
     !Config  Desc = zoom en latitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la latitude
     grossismyy = 1.0
     CALL getin('grossismy',grossismyy)

     IF( ABS(grossismy - grossismyy).GE. 0.001 )  THEN
        write(lunout,*)'conf_gcm: La valeur de grossismy passee par ', &
             'run.def est differente de celle lue sur le fichier  start '
        STOP
     ENDIF

     IF( grossismx.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismx < 1 .   *** '
        STOP
     ELSE
        alphax = 1. - 1./ grossismx
     ENDIF

     IF( grossismy.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismy < 1 .   *** '
        STOP
     ELSE
        alphay = 1. - 1./ grossismy
     ENDIF

     write(lunout,*)'conf_gcm: alphax alphay',alphax,alphay

     !    alphax et alphay sont les anciennes formulat. des grossissements

     !Config  Key  = fxyhypb
     !Config  Desc = Fonction  hyperbolique
     !Config  Def  = y
     !Config  Help = Fonction  f(y)  hyperbolique  si = .true.  
     !Config         sinon  sinusoidale
     fxyhypbb = .TRUE.
     CALL getin('fxyhypb',fxyhypbb)

     IF( .NOT.fxyhypb )  THEN
        IF( fxyhypbb )     THEN
           write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
           write(lunout,*)' *** fxyhypb lu sur le fichier start est ', &
                'F alors  qu il est  T  sur  run.def  ***'
           STOP
        ENDIF
     ELSE
        IF( .NOT.fxyhypbb )   THEN
           write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
           write(lunout,*)' ***  fxyhypb lu sur le fichier start est ', &
                'T alors  qu il est  F  sur  run.def  ****  '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = dzoomx
     !Config  Desc = extension en longitude
     !Config  Def  = 0
     !Config  Help = extension en longitude  de la zone du zoom  
     !Config         ( fraction de la zone totale)
     dzoomxx = 0.0
     CALL getin('dzoomx',dzoomxx)

     IF( fxyhypb )  THEN
        IF( ABS(dzoomx - dzoomxx).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de dzoomx passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = dzoomy
     !Config  Desc = extension en latitude
     !Config  Def  = 0
     !Config  Help = extension en latitude de la zone  du zoom  
     !Config         ( fraction de la zone totale)
     dzoomyy = 0.0
     CALL getin('dzoomy',dzoomyy)

     IF( fxyhypb )  THEN
        IF( ABS(dzoomy - dzoomyy).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de dzoomy passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = taux
     !Config  Desc = raideur du zoom en  X
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  X
     tauxx = 3.0
     CALL getin('taux',tauxx)

     IF( fxyhypb )  THEN
        IF( ABS(taux - tauxx).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de taux passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !Config  Key  = tauyy
     !Config  Desc = raideur du zoom en  Y
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  Y
     tauyy = 3.0
     CALL getin('tauy',tauyy)

     IF( fxyhypb )  THEN
        IF( ABS(tauy - tauyy).GE. 0.001 )  THEN
           write(lunout,*)'conf_gcm: La valeur de tauy passee par ', &
                'run.def est differente de celle lue sur le fichier  start '
           STOP
        ENDIF
     ENDIF

     !c
     IF( .NOT.fxyhypb  )  THEN

        !Config  Key  = ysinus
        !Config  IF   = !fxyhypb
        !Config  Desc = Fonction en Sinus
        !Config  Def  = y
        !Config  Help = Fonction  f(y) avec y = Sin(latit.) si = .true. 
        !Config         sinon y = latit.
        ysinuss = .TRUE.
        CALL getin('ysinus',ysinuss)

        IF( .NOT.ysinus )  THEN
           IF( ysinuss )     THEN
              write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
              write(lunout,*)' *** ysinus lu sur le fichier start est F', &
                   ' alors  qu il est  T  sur  run.def  ***'
              STOP
           ENDIF
        ELSE
           IF( .NOT.ysinuss )   THEN
              write(lunout,*)' ********  PBS DANS  CONF_GCM  ******** '
              write(lunout,*)' *** ysinus lu sur le fichier start est T', &
                   ' alors  qu il est  F  sur  run.def  ****  '
              STOP
           ENDIF
        ENDIF
     ENDIF ! of IF( .NOT.fxyhypb  )

     !Config  Key  = offline
     !Config  Desc = Nouvelle eau liquide
     !Config  Def  = n
     !Config  Help = Permet de mettre en route la
     !Config         nouvelle parametrisation de l'eau liquide !
     offline = .FALSE.
     CALL getin('offline',offline)
     IF (offline .AND. adjust) THEN
        WRITE(lunout,*)  &
             'WARNING : option offline does not work with adjust=y :'
        WRITE(lunout,*) 'the files defstoke.nc, fluxstoke.nc ',  &
             'and fluxstokev.nc will not be created'
        WRITE(lunout,*)  &
             'only the file phystoke.nc will still be created ' 
     END IF

     !Config  Key  = type_trac
     !Config  Desc = Choix de couplage avec model de chimie INCA ou REPROBUS
     !Config  Def  = lmdz
     !Config  Help = 
     !Config         'lmdz' = pas de couplage, pur LMDZ
     !Config         'inca' = model de chime INCA 
     !Config         'repr' = model de chime REPROBUS
     type_trac = 'lmdz'
     CALL getin('type_trac',type_trac)

     !Config  Key  = config_inca
     !Config  Desc = Choix de configuration de INCA
     !Config  Def  = none
     !Config  Help = Choix de configuration de INCA :
     !Config         'none' = sans INCA
     !Config         'chem' = INCA avec calcul de chemie
     !Config         'aero' = INCA avec calcul des aerosols 
     config_inca = 'none'
     CALL getin('config_inca',config_inca)

     !Config  Key  = ok_dynzon 
     !Config  Desc = calcul et sortie des transports 
     !Config  Def  = n 
     !Config  Help = Permet de mettre en route le calcul des transports 
     !Config          
     ok_dynzon = .FALSE. 
     CALL getin('ok_dynzon',ok_dynzon) 

     !Config  Key  = ok_dyn_ins
     !Config  Desc = sorties instantanees dans la dynamique
     !Config  Def  = n 
     !Config  Help = 
     !Config          
     ok_dyn_ins = .FALSE. 
     CALL getin('ok_dyn_ins',ok_dyn_ins) 

     !Config  Key  = ok_dyn_ave
     !Config  Desc = sorties moyennes dans la dynamique
     !Config  Def  = n 
     !Config  Help = 
     !Config          
     ok_dyn_ave = .FALSE. 
     CALL getin('ok_dyn_ave',ok_dyn_ave) 

     write(lunout,*)' #########################################'
     write(lunout,*)' Configuration des parametres du gcm: '
     write(lunout,*)' planet_type = ', planet_type
     write(lunout,*)' calend = ', calend
     write(lunout,*)' dayref = ', dayref
     write(lunout,*)' anneeref = ', anneeref
     write(lunout,*)' nday = ', nday
     write(lunout,*)' day_step = ', day_step
     write(lunout,*)' iperiod = ', iperiod
     write(lunout,*)' nsplit_phys = ', nsplit_phys
     write(lunout,*)' iconser = ', iconser
     write(lunout,*)' iecri = ', iecri
     write(lunout,*)' periodav = ', periodav 
     write(lunout,*)' output_grads_dyn = ', output_grads_dyn
     write(lunout,*)' dissip_period = ', dissip_period
     write(lunout,*)' lstardis = ', lstardis
     write(lunout,*)' nitergdiv = ', nitergdiv
     write(lunout,*)' nitergrot = ', nitergrot
     write(lunout,*)' niterh = ', niterh
     write(lunout,*)' tetagdiv = ', tetagdiv
     write(lunout,*)' tetagrot = ', tetagrot
     write(lunout,*)' tetatemp = ', tetatemp
     write(lunout,*)' coefdis = ', coefdis
     write(lunout,*)' purmats = ', purmats
     write(lunout,*)' read_start = ', read_start
     write(lunout,*)' iflag_phys = ', iflag_phys
     write(lunout,*)' iphysiq = ', iphysiq
     write(lunout,*)' clonn = ', clonn 
     write(lunout,*)' clatt = ', clatt
     write(lunout,*)' grossismx = ', grossismx
     write(lunout,*)' grossismy = ', grossismy
     write(lunout,*)' fxyhypbb = ', fxyhypbb
     write(lunout,*)' dzoomxx = ', dzoomxx
     write(lunout,*)' dzoomy = ', dzoomyy
     write(lunout,*)' tauxx = ', tauxx
     write(lunout,*)' tauyy = ', tauyy
     write(lunout,*)' offline = ', offline
     write(lunout,*)' type_trac = ', type_trac
     write(lunout,*)' config_inca = ', config_inca
     write(lunout,*)' ok_dynzon = ', ok_dynzon 
     write(lunout,*)' ok_dyn_ins = ', ok_dyn_ins 
     write(lunout,*)' ok_dyn_ave = ', ok_dyn_ave 
  else
     !Config  Key  = clon
     !Config  Desc = centre du zoom, longitude
     !Config  Def  = 0
     !Config  Help = longitude en degres du centre 
     !Config         du zoom
     clon = 0.
     CALL getin('clon',clon)

     !Config  Key  = clat
     !Config  Desc = centre du zoom, latitude
     !Config  Def  = 0
     !Config  Help = latitude en degres du centre du zoom
     !Config         
     clat = 0.
     CALL getin('clat',clat)

     !Config  Key  = grossismx 
     !Config  Desc = zoom en longitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la longitude
     grossismx = 1.0
     CALL getin('grossismx',grossismx)

     !Config  Key  = grossismy
     !Config  Desc = zoom en latitude
     !Config  Def  = 1.0
     !Config  Help = facteur de grossissement du zoom,
     !Config         selon la latitude
     grossismy = 1.0
     CALL getin('grossismy',grossismy)

     IF( grossismx.LT.1. )  THEN
        write(lunout,*) &
             'conf_gcm: ***  ATTENTION !! grossismx < 1 .   *** '
        STOP
     ELSE
        alphax = 1. - 1./ grossismx
     ENDIF

     IF( grossismy.LT.1. )  THEN
        write(lunout,*) 'conf_gcm: ***ATTENTION !! grossismy < 1 . *** '
        STOP
     ELSE
        alphay = 1. - 1./ grossismy
     ENDIF

     write(lunout,*)'conf_gcm: alphax alphay ',alphax,alphay

     !    alphax et alphay sont les anciennes formulat. des grossissements

     !Config  Key  = fxyhypb
     !Config  Desc = Fonction  hyperbolique
     !Config  Def  = y
     !Config  Help = Fonction  f(y)  hyperbolique  si = .true.  
     !Config         sinon  sinusoidale
     fxyhypb = .TRUE.
     CALL getin('fxyhypb',fxyhypb)

     !Config  Key  = dzoomx
     !Config  Desc = extension en longitude
     !Config  Def  = 0
     !Config  Help = extension en longitude  de la zone du zoom  
     !Config         ( fraction de la zone totale)
     dzoomx = 0.0
     CALL getin('dzoomx',dzoomx)

     !Config  Key  = dzoomy
     !Config  Desc = extension en latitude
     !Config  Def  = 0
     !Config  Help = extension en latitude de la zone  du zoom  
     !Config         ( fraction de la zone totale)
     dzoomy = 0.0
     CALL getin('dzoomy',dzoomy)

     !Config  Key  = taux
     !Config  Desc = raideur du zoom en  X
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  X
     taux = 3.0
     CALL getin('taux',taux)

     !Config  Key  = tauy
     !Config  Desc = raideur du zoom en  Y
     !Config  Def  = 3
     !Config  Help = raideur du zoom en  Y
     tauy = 3.0
     CALL getin('tauy',tauy)

     !Config  Key  = ysinus
     !Config  IF   = !fxyhypb
     !Config  Desc = Fonction en Sinus
     !Config  Def  = y
     !Config  Help = Fonction  f(y) avec y = Sin(latit.) si = .true. 
     !Config         sinon y = latit.
     ysinus = .TRUE.
     CALL getin('ysinus',ysinus)

     !Config  Key  = offline
     !Config  Desc = Nouvelle eau liquide
     !Config  Def  = n
     !Config  Help = Permet de mettre en route la
     !Config         nouvelle parametrisation de l'eau liquide !
     offline = .FALSE.
     CALL getin('offline',offline)
     IF (offline .AND. adjust) THEN
        WRITE(lunout,*)  &
             'WARNING : option offline does not work with adjust=y :'
        WRITE(lunout,*) 'the files defstoke.nc, fluxstoke.nc ',  &
             'and fluxstokev.nc will not be created'
        WRITE(lunout,*)  &
             'only the file phystoke.nc will still be created ' 
     END IF

     !Config  Key  = type_trac
     !Config  Desc = Choix de couplage avec model de chimie INCA ou REPROBUS
     !Config  Def  = lmdz
     !Config  Help = 
     !Config         'lmdz' = pas de couplage, pur LMDZ
     !Config         'inca' = model de chime INCA 
     !Config         'repr' = model de chime REPROBUS
     type_trac = 'lmdz'
     CALL getin('type_trac',type_trac)

     !Config  Key  = config_inca
     !Config  Desc = Choix de configuration de INCA
     !Config  Def  = none
     !Config  Help = Choix de configuration de INCA :
     !Config         'none' = sans INCA
     !Config         'chem' = INCA avec calcul de chemie
     !Config         'aero' = INCA avec calcul des aerosols 
     config_inca = 'none'
     CALL getin('config_inca',config_inca)

     !Config  Key  = ok_dynzon 
     !Config  Desc = sortie des transports zonaux dans la dynamique
     !Config  Def  = n 
     !Config  Help = Permet de mettre en route le calcul des transports 
     !Config          
     ok_dynzon = .FALSE. 
     CALL getin('ok_dynzon',ok_dynzon) 

     !Config  Key  = ok_dyn_ins
     !Config  Desc = sorties instantanees dans la dynamique
     !Config  Def  = n 
     !Config  Help = 
     !Config          
     ok_dyn_ins = .FALSE. 
     CALL getin('ok_dyn_ins',ok_dyn_ins) 

     !Config  Key  = ok_dyn_ave
     !Config  Desc = sorties moyennes dans la dynamique
     !Config  Def  = n 
     !Config  Help = 
     !Config          
     ok_dyn_ave = .FALSE. 
     CALL getin('ok_dyn_ave',ok_dyn_ave) 

     !Config  Key  = use_filtre_fft
     !Config  Desc = flag d'activation des FFT pour le filtre
     !Config  Def  = false
     !Config  Help = permet d'activer l'utilisation des FFT pour effectuer
     !Config         le filtrage aux poles.
     use_filtre_fft=.FALSE.
     CALL getin('use_filtre_fft',use_filtre_fft)
     use_filtre_fft_loc=use_filtre_fft

     IF (use_filtre_fft .AND. grossismx /= 1.0) THEN
        write(lunout,*)'WARNING !!! '
        write(lunout,*)"Le zoom en longitude est incompatible", &
             " avec l'utilisation du filtre FFT ", &
             "---> FFT filter not active"
        use_filtre_fft=.FALSE.
     ENDIF

     !Config  Key  = use_mpi_alloc
     !Config  Desc = Utilise un buffer MPI en mï¿½moire globale
     !Config  Def  = false
     !Config  Help = permet d'activer l'utilisation d'un buffer MPI
     !Config         en mï¿½moire globale a l'aide de la fonction MPI_ALLOC.
     !Config         Cela peut amï¿½liorer la bande passante des transferts MPI
     !Config         d'un facteur 2  
     use_mpi_alloc=.FALSE.
     CALL getin('use_mpi_alloc',use_mpi_alloc)

     !Config key = ok_strato
     !Config  Desc = activation de la version strato
     !Config  Def  = .FALSE.
     !Config  Help = active la version stratosphérique de LMDZ de F. Lott

     ok_strato=.FALSE.
     CALL getin('ok_strato',ok_strato)

     vert_prof_dissip = merge(1, 0, ok_strato .and. llm==39)
     CALL getin('vert_prof_dissip', vert_prof_dissip)
     call assert(vert_prof_dissip == 0 .or. vert_prof_dissip ==  1, &
          "bad value for vert_prof_dissip")

     !Config  Key  = ok_gradsfile
     !Config  Desc = activation des sorties grads du guidage
     !Config  Def  = n
     !Config  Help = active les sorties grads du guidage

     ok_gradsfile = .FALSE.
     CALL getin('ok_gradsfile',ok_gradsfile)

     !Config  Key  = ok_limit
     !Config  Desc = creation des fichiers limit dans create_etat0_limit
     !Config  Def  = y
     !Config  Help = production du fichier limit.nc requise

     ok_limit = .TRUE.
     CALL getin('ok_limit',ok_limit)

     !Config  Key  = ok_etat0
     !Config  Desc = creation des fichiers etat0 dans create_etat0_limit
     !Config  Def  = y
     !Config  Help = production des fichiers start.nc, startphy.nc requise

     ok_etat0 = .TRUE.
     CALL getin('ok_etat0',ok_etat0)

     write(lunout,*)' #########################################'
     write(lunout,*)' Configuration des parametres de cel0' &
          //'_limit: '
     write(lunout,*)' planet_type = ', planet_type
     write(lunout,*)' calend = ', calend
     write(lunout,*)' dayref = ', dayref
     write(lunout,*)' anneeref = ', anneeref
     write(lunout,*)' nday = ', nday
     write(lunout,*)' day_step = ', day_step
     write(lunout,*)' iperiod = ', iperiod
     write(lunout,*)' iconser = ', iconser
     write(lunout,*)' iecri = ', iecri
     write(lunout,*)' periodav = ', periodav 
     write(lunout,*)' output_grads_dyn = ', output_grads_dyn
     write(lunout,*)' dissip_period = ', dissip_period
     write(lunout,*)' lstardis = ', lstardis
     write(lunout,*)' nitergdiv = ', nitergdiv
     write(lunout,*)' nitergrot = ', nitergrot
     write(lunout,*)' niterh = ', niterh
     write(lunout,*)' tetagdiv = ', tetagdiv
     write(lunout,*)' tetagrot = ', tetagrot
     write(lunout,*)' tetatemp = ', tetatemp
     write(lunout,*)' coefdis = ', coefdis
     write(lunout,*)' purmats = ', purmats
     write(lunout,*)' read_start = ', read_start
     write(lunout,*)' iflag_phys = ', iflag_phys
     write(lunout,*)' iphysiq = ', iphysiq
     write(lunout,*)' clon = ', clon
     write(lunout,*)' clat = ', clat
     write(lunout,*)' grossismx = ', grossismx
     write(lunout,*)' grossismy = ', grossismy
     write(lunout,*)' fxyhypb = ', fxyhypb
     write(lunout,*)' dzoomx = ', dzoomx
     write(lunout,*)' dzoomy = ', dzoomy
     write(lunout,*)' taux = ', taux
     write(lunout,*)' tauy = ', tauy
     write(lunout,*)' offline = ', offline
     write(lunout,*)' type_trac = ', type_trac
     write(lunout,*)' config_inca = ', config_inca
     write(lunout,*)' ok_dynzon = ', ok_dynzon
     write(lunout,*)' ok_dyn_ins = ', ok_dyn_ins 
     write(lunout,*)' ok_dyn_ave = ', ok_dyn_ave 
     write(lunout,*)' use_filtre_fft = ', use_filtre_fft
     write(lunout,*)' use_mpi_alloc = ', use_mpi_alloc
     write(lunout,*)' ok_strato = ', ok_strato
     write(lunout,*)' ok_gradsfile = ', ok_gradsfile
     write(lunout,*)' ok_limit = ', ok_limit
     write(lunout,*)' ok_etat0 = ', ok_etat0
  end IF test_etatinit

END SUBROUTINE conf_gcm
