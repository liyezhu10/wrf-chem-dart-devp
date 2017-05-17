! $Id$
!#define IO_DEBUG

SUBROUTINE physiq (nlon,nlev, &
     debut,lafin,jD_cur, jH_cur,pdtphys, &
     paprs,pplay,pphi,pphis,presnivs,clesphy0, &
     u,v,t,qx, &
     flxmass_w, &
     d_u, d_v, d_t, d_qx, d_ps &
     , dudyn)

  USE ioipsl, only: histbeg, histvert, histdef, histend, histsync, &
       histwrite, ju2ymds, ymds2ju, getin
  USE comgeomphy
  USE phys_cal_mod, only: year_len, mth_len, days_elapsed, jh_1jan, year_cur, &
       mth_cur, phys_cal_update
  USE write_field_phy
  USE dimphy
  USE infotrac
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para
  USE iophy
  USE misc_mod, mydebug=>debug
  USE vampir
  USE pbl_surface_mod, ONLY : pbl_surface
  USE change_srf_frac_mod
  USE surface_data,     ONLY : type_ocean, ok_veget, ok_snow
  USE phys_local_var_mod ! Variables internes non sauvegardees de la physique
  USE phys_state_var_mod ! Variables sauvegardees de la physique
  USE phys_output_var_mod ! Variables pour les ecritures des sorties
  USE phys_output_write_mod
  USE fonte_neige_mod, ONLY  : fonte_neige_get_vars
  USE phys_output_mod
  USE phys_output_ctrlout_mod
  USE iophy
  use open_climoz_m, only: open_climoz ! ozone climatology from a file
  use regr_pr_av_m, only: regr_pr_av
  use netcdf95, only: nf95_close
  !IM for NMC files
  !     use netcdf, only: nf90_fill_real
  use netcdf
  use mod_phys_lmdz_mpi_data, only: is_mpi_root
  USE aero_mod
  use ozonecm_m, only: ozonecm ! ozone of J.-F. Royer
  use conf_phys_m, only: conf_phys
  use radlwsw_m, only: radlwsw
  use phyaqua_mod, only: zenang_an
  USE control_mod
#ifdef REPROBUS
  USE CHEM_REP, ONLY : Init_chem_rep_xjour
#endif
  USE indice_sol_mod
  USE phytrac_mod, ONLY : phytrac

#ifdef CPP_RRTM
  USE YOERAD   , ONLY : NRADLP
#endif

  !IM stations CFMIP
  USE CFMIP_point_locations
  use FLOTT_GWD_rando_m, only: FLOTT_GWD_rando

  IMPLICIT none
  !>======================================================================
  !!
  !! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
  !!
  !! Objet: Moniteur general de la physique du modele
  !!AA      Modifications quant aux traceurs :
  !!AA                  -  uniformisation des parametrisations ds phytrac
  !!AA                  -  stockage des moyennes des champs necessaires
  !!AA                     en mode traceur off-line 
  !!======================================================================
  !!   CLEFS CPP POUR LES IO
  !!   =====================
#define histNMC
  !!======================================================================
  !!    modif   ( P. Le Van ,  12/10/98 )
  !!
  !!  Arguments:
  !!
  !! nlon----input-I-nombre de points horizontaux
  !! nlev----input-I-nombre de couches verticales, doit etre egale a klev
  !! debut---input-L-variable logique indiquant le premier passage
  !! lafin---input-L-variable logique indiquant le dernier passage
  !! jD_cur       -R-jour courant a l'appel de la physique (jour julien)
  !! jH_cur       -R-heure courante a l'appel de la physique (jour julien)
  !! pdtphys-input-R-pas d'integration pour la physique (seconde)
  !! paprs---input-R-pression pour chaque inter-couche (en Pa)
  !! pplay---input-R-pression pour le mileu de chaque couche (en Pa)
  !! pphi----input-R-geopotentiel de chaque couche (g z) (reference sol)
  !! pphis---input-R-geopotentiel du sol
  !! presnivs-input_R_pressions approximat. des milieux couches ( en PA)
  !! u-------input-R-vitesse dans la direction X (de O a E) en m/s
  !! v-------input-R-vitesse Y (de S a N) en m/s
  !! t-------input-R-temperature (K)
  !! qx------input-R-humidite specifique (kg/kg) et d'autres traceurs
  !! d_t_dyn-input-R-tendance dynamique pour "t" (K/s)
  !! d_q_dyn-input-R-tendance dynamique pour "q" (kg/kg/s)
  !! flxmass_w -input-R- flux de masse verticale
  !! d_u-----output-R-tendance physique de "u" (m/s/s)
  !! d_v-----output-R-tendance physique de "v" (m/s/s)
  !! d_t-----output-R-tendance physique de "t" (K/s)
  !! d_qx----output-R-tendance physique de "qx" (kg/kg/s)
  !! d_ps----output-R-tendance physique de la pression au sol
  !!======================================================================
  include "dimensions.h"
  integer jjmp1
  parameter (jjmp1=jjm+1-1/jjm)
  integer iip1
  parameter (iip1=iim+1)

  include "regdim.h"
  include "dimsoil.h"
  include "clesphys.h"
  include "temps.h"
  include "iniprint.h"
  include "thermcell.h"
  !======================================================================
  LOGICAL ok_cvl  ! pour activer le nouveau driver pour convection KE
  PARAMETER (ok_cvl=.TRUE.)
  LOGICAL ok_gust ! pour activer l'effet des gust sur flux surface
  PARAMETER (ok_gust=.FALSE.)
  integer iflag_radia     ! active ou non le rayonnement (MPL)
  save iflag_radia
  !$OMP THREADPRIVATE(iflag_radia)
  !======================================================================
  LOGICAL check ! Verifier la conservation du modele en eau
  PARAMETER (check=.FALSE.)
  LOGICAL ok_stratus ! Ajouter artificiellement les stratus
  PARAMETER (ok_stratus=.FALSE.)
  !======================================================================
  REAL amn, amx
  INTEGER igout
  !======================================================================
  ! Clef controlant l'activation du cycle diurne:
  !cc      LOGICAL cycle_diurne
  !cc      PARAMETER (cycle_diurne=.FALSE.)
  !======================================================================
  ! Modele thermique du sol, a activer pour le cycle diurne:
  !cc      LOGICAL soil_model
  !cc      PARAMETER (soil_model=.FALSE.)
  !======================================================================
  ! Dans les versions precedentes, l'eau liquide nuageuse utilisee dans
  ! le calcul du rayonnement est celle apres la precipitation des nuages.
  ! Si cette cle new_oliq est activee, ce sera une valeur moyenne entre
  ! la condensation et la precipitation. Cette cle augmente les impacts
  ! radiatifs des nuages.
  !cc      LOGICAL new_oliq
  !cc      PARAMETER (new_oliq=.FALSE.)
  !======================================================================
  ! Clefs controlant deux parametrisations de l'orographie:
  !c      LOGICAL ok_orodr
  !cc      PARAMETER (ok_orodr=.FALSE.)
  !cc      LOGICAL ok_orolf
  !cc      PARAMETER (ok_orolf=.FALSE.)
  !======================================================================
  LOGICAL ok_journe ! sortir le fichier journalier
  save ok_journe
  !$OMP THREADPRIVATE(ok_journe)
  !
  LOGICAL ok_mensuel ! sortir le fichier mensuel
  save ok_mensuel
  !$OMP THREADPRIVATE(ok_mensuel)
  !
  LOGICAL ok_instan ! sortir le fichier instantane
  save ok_instan
  !$OMP THREADPRIVATE(ok_instan)
  !
  LOGICAL ok_LES ! sortir le fichier LES 
  save ok_LES                            
  !$OMP THREADPRIVATE(ok_LES)                  
  !
  LOGICAL callstats ! sortir le fichier stats 
  save callstats                            
  !$OMP THREADPRIVATE(callstats)                  
  !
  LOGICAL ok_region ! sortir le fichier regional
  PARAMETER (ok_region=.FALSE.)
  !======================================================================
  real seuil_inversion
  save seuil_inversion
  !$OMP THREADPRIVATE(seuil_inversion)
  integer iflag_ratqs
  save iflag_ratqs
  !$OMP THREADPRIVATE(iflag_ratqs)
  real facteur

  REAL wmax_th(klon)
  REAL tau_overturning_th(klon)

  integer lmax_th(klon)
  integer limbas(klon)
  real ratqscth(klon,klev)
  real ratqsdiff(klon,klev)
  real zqsatth(klon,klev)

  !======================================================================
  !
  INTEGER ivap          ! indice de traceurs pour vapeur d'eau
  PARAMETER (ivap=1)
  INTEGER iliq          ! indice de traceurs pour eau liquide
  PARAMETER (iliq=2)
!CR: on ajoute la phase glace
  INTEGER isol          ! indice de traceurs pour eau glace
  PARAMETER (isol=3)
  !
  !
  ! Variables argument:
  !
  INTEGER nlon
  INTEGER nlev
  REAL, intent(in):: jD_cur, jH_cur

  REAL pdtphys
  LOGICAL debut, lafin
  REAL paprs(klon,klev+1)

! ---- DART-LMDZ -----------------------
  REAL paprs1(klon,klev+1)         
  REAL paprs_glo(klon_glo,klev+1) 
  REAL wo_save_glo(klon_glo,klev)
  REAL wo_save_loc(klon,klev)  
! ---- DART-LMDZ -----------------------

  REAL pplay(klon,klev)
  REAL pphi(klon,klev)
  REAL pphis(klon)
  REAL presnivs(klev)
  REAL znivsig(klev)
  real pir

  REAL u(klon,klev)
  REAL v(klon,klev)
  REAL t(klon,klev),thetal(klon,klev)
  ! thetal: ligne suivante a decommenter si vous avez les fichiers     MPL 20130625
  ! fth_fonctions.F90 et parkind1.F90
  ! sinon thetal=theta
  !     REAL fth_thetae,fth_thetav,fth_thetal
  REAL qx(klon,klev,nqtot)
  REAL flxmass_w(klon,klev)
  REAL d_u(klon,klev)
  REAL d_v(klon,klev)
  REAL d_t(klon,klev)
  REAL d_qx(klon,klev,nqtot)
  REAL d_ps(klon)
  ! Variables pour le transport convectif
  real da(klon,klev),phi(klon,klev,klev),mp(klon,klev)
  real wght_cvfd(klon,klev)
  ! Variables pour le lessivage convectif
  ! RomP >>> 
  real phi2(klon,klev,klev)
  real d1a(klon,klev),dam(klon,klev)
  real ev(klon,klev),ep(klon,klev)
  real clw(klon,klev),elij(klon,klev,klev)
  real epmlmMm(klon,klev,klev),eplaMm(klon,klev)
  ! RomP <<<
  !IM definition dynamique o_trac dans phys_output_open
  !      type(ctrl_out) :: o_trac(nqtot)

  ! variables a une pression donnee
  !
  include "declare_STDlev.h"
  !
  !
  include "radopt.h"
  !
  !


  INTEGER debug
  INTEGER n
  !ym      INTEGER npoints
  !ym      PARAMETER(npoints=klon) 
  !
  INTEGER nregISCtot
  PARAMETER(nregISCtot=1) 
  !
  ! imin_debut, nbpti, jmin_debut, nbptj : parametres pour sorties sur 1 region rectangulaire
  ! y compris pour 1 point
  ! imin_debut : indice minimum de i; nbpti : nombre de points en direction i (longitude)
  ! jmin_debut : indice minimum de j; nbptj : nombre de points en direction j (latitude)
  INTEGER imin_debut, nbpti
  INTEGER jmin_debut, nbptj 
  !IM: region='3d' <==> sorties en global
  CHARACTER*3 region
  PARAMETER(region='3d')
  logical ok_hf
  !
  save ok_hf
  !$OMP THREADPRIVATE(ok_hf)

  INTEGER        longcles
  PARAMETER    ( longcles = 20 )
  REAL clesphy0( longcles      )
  !
  ! Variables propres a la physique
  INTEGER itap
  SAVE itap                   ! compteur pour la physique
  !$OMP THREADPRIVATE(itap)
  !
  REAL,save ::  solarlong0
  !$OMP THREADPRIVATE(solarlong0)

  !
  !  Parametres de l'Orographie a l'Echelle Sous-Maille (OESM):
  !
  !IM 141004     REAL zulow(klon),zvlow(klon),zustr(klon), zvstr(klon)
  REAL zulow(klon),zvlow(klon)
  !
  INTEGER igwd,idx(klon),itest(klon)
  !
  !      REAL,allocatable,save :: run_off_lic_0(:)
!!$OMP THREADPRIVATE(run_off_lic_0)
  !ym      SAVE run_off_lic_0
  !KE43
  ! Variables liees a la convection de K. Emanuel (sb):
  !
  REAL bas, top             ! cloud base and top levels
  SAVE bas
  SAVE top
  !$OMP THREADPRIVATE(bas, top)

  !
  !=================================================================================================
  !CR04.12.07: on ajoute les nouvelles variables du nouveau schema de convection avec poches froides
  ! Variables li\'ees \`a la poche froide (jyg)

  REAL mip(klon,klev)  ! mass flux shed by the adiab ascent at each level
  !
  REAL wape_prescr, fip_prescr
  INTEGER it_wape_prescr
  SAVE wape_prescr, fip_prescr, it_wape_prescr
  !$OMP THREADPRIVATE(wape_prescr, fip_prescr, it_wape_prescr)
  !
  ! variables supplementaires de concvl
  REAL Tconv(klon,klev)
  REAL sij(klon,klev,klev)

  real, save :: alp_bl_prescr=0.
  real, save :: ale_bl_prescr=0.

  real, save :: ale_max=1000.
  real, save :: alp_max=2.

  real, save :: wake_s_min_lsp=0.1

  !$OMP THREADPRIVATE(alp_bl_prescr,ale_bl_prescr)
  !$OMP THREADPRIVATE(ale_max,alp_max)
  !$OMP THREADPRIVATE(wake_s_min_lsp)


  real ok_wk_lsp(klon)

  !RC
  ! Variables li\'ees \`a la poche froide (jyg et rr)
  ! Version diagnostique pour l'instant : pas de r\'etroaction sur la convection

  REAL t_wake(klon,klev),q_wake(klon,klev) ! wake pour la convection

  REAL wake_dth(klon,klev)        ! wake : temp pot difference

  REAL wake_d_deltat_gw(klon,klev)! wake : delta T tendency due to Gravity Wave (/s)
  REAL wake_omgbdth(klon,klev)    ! Wake : flux of Delta_Theta transported by LS omega
  REAL wake_dp_omgb(klon,klev)    ! Wake : vertical gradient of large scale omega
  REAL wake_dtKE(klon,klev)       ! Wake : differential heating (wake - unpertubed) CONV
  REAL wake_dqKE(klon,klev)       ! Wake : differential moistening (wake - unpertubed) CONV
  REAL wake_dtPBL(klon,klev)      ! Wake : differential heating (wake - unpertubed) PBL
  REAL wake_dqPBL(klon,klev)      ! Wake : differential moistening (wake - unpertubed) PBL
  REAL wake_ddeltat(klon,klev),wake_ddeltaq(klon,klev)
  REAL wake_dp_deltomg(klon,klev) ! Wake : gradient vertical de wake_omg
  REAL wake_spread(klon,klev)     ! spreading term in wake_delt
  !
  !pourquoi y'a pas de save??
  !
  INTEGER wake_k(klon)            ! Wake sommet
  !
  REAL t_undi(klon,klev)               ! temperature moyenne dans la zone non perturbee
  REAL q_undi(klon,klev)               ! humidite moyenne dans la zone non perturbee
  !
  !jyg<
  !cc      REAL wake_pe(klon)              ! Wake potential energy - WAPE 
  !>jyg

  REAL wake_gfl(klon)             ! Gust Front Length
  REAL wake_dens(klon)
  !
  !
  REAL dt_dwn(klon,klev)
  REAL dq_dwn(klon,klev)
  REAL wdt_PBL(klon,klev)
  REAL udt_PBL(klon,klev)
  REAL wdq_PBL(klon,klev)
  REAL udq_PBL(klon,klev)
  REAL M_dwn(klon,klev)
  REAL M_up(klon,klev)
  REAL dt_a(klon,klev)
  REAL dq_a(klon,klev)
  REAL, dimension(klon) :: www
  REAL, SAVE :: alp_offset
  !$OMP THREADPRIVATE(alp_offset)

!!!
!=================================================================
!         PROVISOIRE : DECOUPLAGE PBL/WAKE
!         --------------------------------
      REAL wake_deltat_sav(klon,klev)
      REAL wake_deltaq_sav(klon,klev)
!=================================================================

  !
  !RR:fin declarations poches froides
  !=======================================================================================================

  REAL ztv(klon,klev),ztva(klon,klev)
  REAL zpspsk(klon,klev)
  REAL ztla(klon,klev),zqla(klon,klev) 
  REAL zthl(klon,klev)

  !cc nrlmd le 10/04/2012

  !--------Stochastic Boundary Layer Triggering: ALE_BL--------
  !---Propri\'et\'es du thermiques au LCL 
  real zlcl_th(klon)                                     ! Altitude du LCL calcul\'e continument (pcon dans thermcell_main.F90)
  real fraca0(klon)                                      ! Fraction des thermiques au LCL
  real w0(klon)                                          ! Vitesse des thermiques au LCL
  real w_conv(klon)                                      ! Vitesse verticale de grande \'echelle au LCL
  real tke0(klon,klev+1)                                 ! TKE au début du pas de temps
  real therm_tke_max0(klon)                              ! TKE dans les thermiques au LCL 
  real env_tke_max0(klon)                                ! TKE dans l'environnement au LCL 

  !---D\'eclenchement stochastique
  integer :: tau_trig(klon)

  !--------Statistical Boundary Layer Closure: ALP_BL--------
  !---Profils de TKE dans et hors du thermique
  real therm_tke_max(klon,klev)                          ! Profil de TKE dans les thermiques
  real env_tke_max(klon,klev)                            ! Profil de TKE dans l'environnement


  !cc fin nrlmd le 10/04/2012

  ! Variables locales pour la couche limite (al1):
  !
  !Al1      REAL pblh(klon)           ! Hauteur de couche limite
  !Al1      SAVE pblh
  !34EK
  !
  ! Variables locales:
  !
  !AA
  !AA  Pour phytrac 
  REAL u1(klon)             ! vents dans la premiere couche U
  REAL v1(klon)             ! vents dans la premiere couche V

  !@$$      LOGICAL offline           ! Controle du stockage ds "physique"
  !@$$      PARAMETER (offline=.false.)
  !@$$      INTEGER physid
  REAL frac_impa(klon,klev) ! fractions d'aerosols lessivees (impaction)
  REAL frac_nucl(klon,klev) ! idem (nucleation)
  ! RomP >>> 
  REAL beta_prec_fisrt(klon,klev) ! taux de conv de l'eau cond (fisrt)
  ! RomP <<<

  REAL          :: calday

  !IM cf FH pour Tiedtke 080604
  REAL rain_tiedtke(klon),snow_tiedtke(klon)
  !
  !IM 050204 END
  REAL devap(klon) ! evaporation et sa derivee
  REAL dsens(klon) ! chaleur sensible et sa derivee

  !
  ! Conditions aux limites
  !
  !
  REAL :: day_since_equinox
  ! Date de l'equinoxe de printemps
  INTEGER, parameter :: mth_eq=3, day_eq=21
  REAL :: jD_eq

  LOGICAL, parameter :: new_orbit = .true.

  !
  INTEGER lmt_pas
  SAVE lmt_pas                ! frequence de mise a jour
  !$OMP THREADPRIVATE(lmt_pas) 
  real zmasse(klon, llm),exner(klon, llm) 
  !     (column-density of mass of air in a cell, in kg m-2)
  real, parameter:: dobson_u = 2.1415e-05 ! Dobson unit, in kg m-2

  !IM sorties
  REAL un_jour
  PARAMETER(un_jour=86400.)
  INTEGER itapm1 !pas de temps de la physique du(es) mois precedents
  SAVE itapm1    !mis a jour le dernier pas de temps du mois en cours
  !$OMP THREADPRIVATE(itapm1)
  !======================================================================
  !
  ! Declaration des procedures appelees
  !
  EXTERNAL angle     ! calculer angle zenithal du soleil
  EXTERNAL alboc     ! calculer l'albedo sur ocean
  EXTERNAL ajsec     ! ajustement sec
  EXTERNAL conlmd    ! convection (schema LMD)
  !KE43
  EXTERNAL conema3  ! convect4.3
  EXTERNAL fisrtilp  ! schema de condensation a grande echelle (pluie)
  !AA 
  ! JBM (3/14) fisrtilp_tr not loaded
  ! EXTERNAL fisrtilp_tr ! schema de condensation a grande echelle (pluie)
  !                          ! stockage des coefficients necessaires au
  !                          ! lessivage OFF-LINE et ON-LINE
  EXTERNAL hgardfou  ! verifier les temperatures
  EXTERNAL nuage     ! calculer les proprietes radiatives
  !C      EXTERNAL o3cm      ! initialiser l'ozone
  EXTERNAL orbite    ! calculer l'orbite terrestre
  EXTERNAL phyetat0  ! lire l'etat initial de la physique
  EXTERNAL phyredem  ! ecrire l'etat de redemarrage de la physique
  EXTERNAL suphel    ! initialiser certaines constantes
  EXTERNAL transp    ! transport total de l'eau et de l'energie
  EXTERNAL ecribina  ! ecrire le fichier binaire global
  EXTERNAL ecribins  ! ecrire le fichier binaire global
  EXTERNAL ecrirega  ! ecrire le fichier binaire regional
  EXTERNAL ecriregs  ! ecrire le fichier binaire regional
  !IM
  EXTERNAL haut2bas  !variables de haut en bas
  EXTERNAL ini_undefSTD  !initialise a 0 une variable a 1 niveau de pression
  EXTERNAL undefSTD      !somme les valeurs definies d'1 var a 1 niveau de pression
  !     EXTERNAL moy_undefSTD  !moyenne d'1 var a 1 niveau de pression
  !     EXTERNAL moyglo_aire   !moyenne globale d'1 var ponderee par l'aire de la maille (moyglo_pondaire)
  !                            !par la masse/airetot (moyglo_pondaima) et la vraie masse (moyglo_pondmass)
  !
  ! Variables locales
  !
  REAL rhcl(klon,klev)    ! humiditi relative ciel clair
  REAL dialiq(klon,klev)  ! eau liquide nuageuse
  REAL diafra(klon,klev)  ! fraction nuageuse
  REAL cldliq(klon,klev)  ! eau liquide nuageuse
  !
  !XXX PB 
  REAL fluxq(klon,klev, nbsrf)   ! flux turbulent d'humidite
  !
  REAL zxfluxt(klon, klev)
  REAL zxfluxq(klon, klev)
  REAL zxfluxu(klon, klev)
  REAL zxfluxv(klon, klev)

  ! Le rayonnement n'est pas calcule tous les pas, il faut donc
  !                      sauvegarder les sorties du rayonnement
  !ym      SAVE  heat,cool,albpla,topsw,toplw,solsw,sollw,sollwdown
  !ym      SAVE  sollwdownclr, toplwdown, toplwdownclr
  !ym      SAVE  topsw0,toplw0,solsw0,sollw0, heat0, cool0
  !
  INTEGER itaprad
  SAVE itaprad
  !$OMP THREADPRIVATE(itaprad)
  !
  REAL conv_q(klon,klev) ! convergence de l'humidite (kg/kg/s)
  REAL conv_t(klon,klev) ! convergence de la temperature(K/s)

  !
  !  REAL zxsnow(klon)
  REAL zxsnow_dummy(klon)
  !
  REAL dist, rmu0(klon), fract(klon)
  REAL zdtime, zlongi
  !
  REAL qcheck
  REAL z_avant(klon), z_apres(klon), z_factor(klon)
  LOGICAL zx_ajustq
  !
  REAL za, zb
  REAL zx_t, zx_qs, zdelta, zcor, zlvdcp, zlsdcp
  real zqsat(klon,klev)
  INTEGER i, k, iq, ig, j, nsrf, ll, l, iiq
  REAL t_coup
  PARAMETER (t_coup=234.0)

  !ym A voir plus tard !!
  !ym      REAL zx_relief(iim,jjmp1)
  !ym      REAL zx_aire(iim,jjmp1)
  !
  ! Grandeurs de sorties
  REAL s_capCL(klon)
  REAL s_oliqCL(klon), s_cteiCL(klon)
  REAL s_trmb1(klon), s_trmb2(klon)
  REAL s_trmb3(klon)
  !KE43
  ! Variables locales pour la convection de K. Emanuel (sb):

  REAL tvp(klon,klev)       ! virtual temp of lifted parcel
  CHARACTER*40 capemaxcels  !max(CAPE)

  REAL rflag(klon)          ! flag fonctionnement de convect
  INTEGER iflagctrl(klon)          ! flag fonctionnement de convect

  ! -- convect43:
  INTEGER ntra              ! nb traceurs pour convect4.3
  REAL dtvpdt1(klon,klev), dtvpdq1(klon,klev)
  REAL dplcldt(klon), dplcldr(klon)
  !?     .     condm_con(klon,klev),conda_con(klon,klev),
  !?     .     mr_con(klon,klev),ep_con(klon,klev)
  !?     .    ,sadiab(klon,klev),wadiab(klon,klev)
  ! --
  !34EK
  !
  ! Variables du changement
  !
  ! con: convection
  ! lsc: condensation a grande echelle (Large-Scale-Condensation)
  ! ajs: ajustement sec
  ! eva: evaporation de l'eau liquide nuageuse
  ! vdf: couche limite (Vertical DiFfusion)

  ! tendance nulles
  REAL, dimension(klon,klev):: du0, dv0, dt0, dq0, dql0, dqi0

  !
  !********************************************************
  !     declarations

  !********************************************************
  !IM 081204 END
  !
  REAL pen_u(klon,klev), pen_d(klon,klev)
  REAL pde_u(klon,klev), pde_d(klon,klev)
  INTEGER kcbot(klon), kctop(klon), kdtop(klon)
  !
  REAL ratqsc(klon,klev)
  real ratqsbas,ratqshaut,tau_ratqs
  save ratqsbas,ratqshaut,tau_ratqs
  !$OMP THREADPRIVATE(ratqsbas,ratqshaut,tau_ratqs)

  ! Parametres lies au nouveau schema de nuages (SB, PDF)
  real fact_cldcon
  real facttemps
  logical ok_newmicro
  save ok_newmicro
  !$OMP THREADPRIVATE(ok_newmicro)
  !real ref_liq_pi(klon,klev), ref_ice_pi(klon,klev)
  save fact_cldcon,facttemps
  !$OMP THREADPRIVATE(fact_cldcon,facttemps)

  integer iflag_cldcon
  save iflag_cldcon
  !$OMP THREADPRIVATE(iflag_cldcon)
  logical ptconv(klon,klev)
  !IM cf. AM 081204 BEG
  logical ptconvth(klon,klev)
  !IM cf. AM 081204 END
  !
  ! Variables liees a l'ecriture de la bande histoire physique
  !
  !======================================================================
  !

  !
  integer itau_w   ! pas de temps ecriture = itap + itau_phy
  !
  !
  ! Variables locales pour effectuer les appels en serie
  !
  !IM RH a 2m (la surface)
  REAL Lheat

  INTEGER        length
  PARAMETER    ( length = 100 )
  REAL tabcntr0( length       )
  !
  INTEGER ndex2d(iim*jjmp1)
  !IM
  !
  !IM AMIP2 BEG
  REAL moyglo, mountor
  !IM 141004 BEG
  REAL zustrdr(klon), zvstrdr(klon)
  REAL zustrli(klon), zvstrli(klon)
  REAL zustrph(klon), zvstrph(klon)
  REAL zustrhi(klon), zvstrhi(klon)
  REAL aam, torsfc
  !IM 141004 END
  !IM 190504 BEG
  INTEGER ij, imp1jmp1
  PARAMETER(imp1jmp1=(iim+1)*jjmp1)
  !ym A voir plus tard
  REAL zx_tmp(imp1jmp1), airedyn(iim+1,jjmp1)
  REAL padyn(iim+1,jjmp1,klev+1)
  REAL dudyn(iim+1,jjmp1,klev)
  REAL rlatdyn(iim+1,jjmp1)
  !IM 190504 END
  LOGICAL ok_msk
  REAL msk(klon)
  !IM 
  REAL airetot, pi
  !ym A voir plus tard
  !ym      REAL zm_wo(jjmp1, klev)
  !IM AMIP2 END
  !
  REAL zx_tmp_fi2d(klon)      ! variable temporaire grille physique
  REAL zx_tmp_fi3d(klon,klev) ! variable temporaire pour champs 3D 
  REAL zx_tmp_2d(iim,jjmp1)
  REAL zx_lon(iim,jjmp1), zx_lat(iim,jjmp1)
  !
  INTEGER nid_day_seri, nid_ctesGCM
  SAVE nid_day_seri, nid_ctesGCM
  !$OMP THREADPRIVATE(nid_day_seri,nid_ctesGCM)
  !
  !IM 280405 BEG
  !  INTEGER nid_bilKPins, nid_bilKPave
  !  SAVE nid_bilKPins, nid_bilKPave
  !  !$OMP THREADPRIVATE(nid_bilKPins, nid_bilKPave)
  !
  REAL ve_lay(klon,klev) ! transport meri. de l'energie a chaque niveau vert.
  REAL vq_lay(klon,klev) ! transport meri. de l'eau a chaque niveau vert.
  REAL ue_lay(klon,klev) ! transport zonal de l'energie a chaque niveau vert.
  REAL uq_lay(klon,klev) ! transport zonal de l'eau a chaque niveau vert.
  !
  INTEGER nhori, nvert
  REAL zsto
  REAL zstophy, zout

  real zjulian
  save zjulian
  !$OMP THREADPRIVATE(zjulian)

  character*20 modname
  character*80 abort_message
  logical, save ::  ok_sync, ok_sync_omp
  !$OMP THREADPRIVATE(ok_sync)
  real date0
  integer idayref

  ! essai writephys
  integer fid_day, fid_mth, fid_ins
  parameter (fid_ins = 1, fid_day = 2, fid_mth = 3) 
  integer prof2d_on, prof3d_on, prof2d_av, prof3d_av
  parameter (prof2d_on = 1, prof3d_on = 2, &
       prof2d_av = 3, prof3d_av = 4)
  !     Variables liees au bilan d'energie et d'enthalpi
  REAL ztsol(klon)
  REAL      d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec
  REAL      d_h_vcol_phy
  REAL      fs_bound, fq_bound
  SAVE      d_h_vcol_phy
  !$OMP THREADPRIVATE(d_h_vcol_phy)
  REAL      zero_v(klon)
  CHARACTER*40 ztit
  INTEGER   ip_ebil  ! PRINT level for energy conserv. diag.
  SAVE      ip_ebil
  DATA      ip_ebil/0/
  !$OMP THREADPRIVATE(ip_ebil)
  INTEGER   if_ebil ! level for energy conserv. dignostics
  SAVE      if_ebil
  !$OMP THREADPRIVATE(if_ebil)
  REAL q2m(klon,nbsrf)  ! humidite a 2m

  !IM: t2m, q2m, ustar, u10m, v10m et t2mincels, t2maxcels
  CHARACTER*40 t2mincels, t2maxcels       !t2m min., t2m max
  CHARACTER*40 tinst, tave, typeval
  REAL cldtaupi(klon,klev)  ! Cloud optical thickness for pre-industrial (pi) aerosols


  ! Aerosol optical properties
  CHARACTER*4, DIMENSION(naero_grp) :: rfname 
  REAL, DIMENSION(klon,klev)     :: mass_solu_aero    ! total mass concentration for all soluble aerosols[ug/m3]
  REAL, DIMENSION(klon,klev)     :: mass_solu_aero_pi ! - " - (pre-industrial value)

  ! Parameters
  LOGICAL ok_ade, ok_aie    ! Apply aerosol (in)direct effects or not
  LOGICAL ok_cdnc          ! ok cloud droplet number concentration (O. Boucher 01-2013)
  REAL bl95_b0, bl95_b1   ! Parameter in Boucher and Lohmann (1995)
  SAVE ok_ade, ok_aie, ok_cdnc, bl95_b0, bl95_b1
  !$OMP THREADPRIVATE(ok_ade, ok_aie, ok_cdnc, bl95_b0, bl95_b1)
  LOGICAL, SAVE :: aerosol_couple ! true  : calcul des aerosols dans INCA
  ! false : lecture des aerosol dans un fichier
  !$OMP THREADPRIVATE(aerosol_couple)    
  INTEGER, SAVE :: flag_aerosol 
  !$OMP THREADPRIVATE(flag_aerosol) 
  LOGICAL, SAVE :: new_aod
  !$OMP THREADPRIVATE(new_aod) 
  !
  !--STRAT AEROSOL
  LOGICAL, SAVE :: flag_aerosol_strat
  !$OMP THREADPRIVATE(flag_aerosol_strat)
  !c-fin STRAT AEROSOL
  !
  ! Declaration des constantes et des fonctions thermodynamiques
  !
  LOGICAL,SAVE :: first=.true.
  !$OMP THREADPRIVATE(first)

  integer, save::  read_climoz ! read ozone climatology
  !     (let it keep the default OpenMP shared attribute)
  !     Allowed values are 0, 1 and 2
  !     0: do not read an ozone climatology
  !     1: read a single ozone climatology that will be used day and night
  !     2: read two ozone climatologies, the average day and night
  !     climatology and the daylight climatology

  integer, save:: ncid_climoz ! NetCDF file containing ozone climatologies
  !     (let it keep the default OpenMP shared attribute)

  real, pointer, save:: press_climoz(:)
  !     (let it keep the default OpenMP shared attribute)
  !     edges of pressure intervals for ozone climatologies, in Pa, in strictly
  !     ascending order

  integer, save:: co3i = 0
  !     time index in NetCDF file of current ozone fields
  !$OMP THREADPRIVATE(co3i) 

  integer ro3i
  !     required time index in NetCDF file for the ozone fields, between 1
  !     and 360

  INTEGER ierr
  include "YOMCST.h"
  include "YOETHF.h"
  include "FCTTRE.h"
  !IM 100106 BEG : pouvoir sortir les ctes de la physique
  include "conema3.h"
  include "fisrtilp.h"
  include "nuage.h"
  include "compbl.h"
  !IM 100106 END : pouvoir sortir les ctes de la physique
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Declarations pour Simulateur COSP
  !============================================================
  real :: mr_ozone(klon,klev)

  !IM sorties fichier 1D paramLMDZ_phy.nc
  REAL :: zx_tmp_0d(1,1)
  INTEGER, PARAMETER :: np=1
  REAL,dimension(klon_glo)        :: rlat_glo
  REAL,dimension(klon_glo)        :: rlon_glo
  REAL gbils(1), gevap(1), gevapt(1), glat(1), gnet0(1), gnet(1)
  REAL grain(1), gtsol(1), gt2m(1), gprw(1)

  !IM stations CFMIP
  INTEGER, SAVE :: nCFMIP
  !$OMP THREADPRIVATE(nCFMIP)
  INTEGER, PARAMETER :: npCFMIP=120
  INTEGER, ALLOCATABLE, SAVE :: tabCFMIP(:)
  REAL, ALLOCATABLE, SAVE :: lonCFMIP(:), latCFMIP(:)
  !$OMP THREADPRIVATE(tabCFMIP, lonCFMIP, latCFMIP) 
  INTEGER, ALLOCATABLE, SAVE :: tabijGCM(:)
  REAL, ALLOCATABLE, SAVE :: lonGCM(:), latGCM(:)
  !$OMP THREADPRIVATE(tabijGCM, lonGCM, latGCM)
  INTEGER, ALLOCATABLE, SAVE :: iGCM(:), jGCM(:)
  !$OMP THREADPRIVATE(iGCM, jGCM)
  logical, dimension(nfiles)            :: phys_out_filestations
  logical, parameter :: lNMC=.FALSE.

  !IM betaCRF
  REAL, SAVE :: pfree, beta_pbl, beta_free
  !$OMP THREADPRIVATE(pfree, beta_pbl, beta_free)
  REAL, SAVE :: lon1_beta,  lon2_beta, lat1_beta, lat2_beta
  !$OMP THREADPRIVATE(lon1_beta,  lon2_beta, lat1_beta, lat2_beta)
  LOGICAL, SAVE :: mskocean_beta
  !$OMP THREADPRIVATE(mskocean_beta)
  REAL, dimension(klon, klev) :: beta         ! facteur sur cldtaurad et cldemirad pour evaluer les retros liees aux CRF
  REAL, dimension(klon, klev) :: cldtaurad    ! epaisseur optique pour radlwsw pour tester "CRF off"
  REAL, dimension(klon, klev) :: cldtaupirad  ! epaisseur optique pour radlwsw pour tester "CRF off"
  REAL, dimension(klon, klev) :: cldemirad    ! emissivite pour radlwsw pour tester "CRF off"
  REAL, dimension(klon, klev) :: cldfrarad    ! fraction nuageuse

  INTEGER :: nbtr_tmp ! Number of tracer inside concvl
  REAL, dimension(klon,klev) :: sh_in ! Specific humidity entering in phytrac
  integer iostat

  REAL zzz

  !======================================================================
  ! Gestion calendrier : mise a jour du module phys_cal_mod
  !
  CALL phys_cal_update(jD_cur,jH_cur)

  !======================================================================
  ! Ecriture eventuelle d'un profil verticale en entree de la physique.
  ! Utilise notamment en 1D mais peut etre active egalement en 3D
  ! en imposant la valeur de igout.
  !======================================================================d
  if (prt_level.ge.1) then
     igout=klon/2+1/klon
     write(lunout,*) 'DEBUT DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
     write(lunout,*) &
          'nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys'
     write(lunout,*) &
          nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur,pdtphys 

     write(lunout,*) 'paprs, play, phi, u, v, t'
     do k=1,klev
        write(lunout,*) paprs(igout,k),pplay(igout,k),pphi(igout,k), &
             u(igout,k),v(igout,k),t(igout,k)
     enddo
     write(lunout,*) 'ovap (g/kg),  oliq (g/kg)'
     do k=1,klev
        write(lunout,*) qx(igout,k,1)*1000,qx(igout,k,2)*1000.
     enddo
  endif

  !======================================================================

  if (first) then 

     !CR:nvelles variables convection/poches froides

     print*, '================================================='
     print*, 'Allocation des variables locales et sauvegardees'
     call phys_local_var_init
     !
     pasphys=pdtphys
     !     appel a la lecture du run.def physique
     call conf_phys(ok_journe, ok_mensuel, &
          ok_instan, ok_hf, &
          ok_LES, &
          callstats, &
          solarlong0,seuil_inversion, &
          fact_cldcon, facttemps,ok_newmicro,iflag_radia, &
          iflag_cldcon,iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs, &
          ok_ade, ok_aie, ok_cdnc, aerosol_couple,  &
          flag_aerosol, flag_aerosol_strat, new_aod, &
          bl95_b0, bl95_b1, &
          !     nv flags pour la convection et les poches froides
          read_climoz, &
          alp_offset)
     call phys_state_var_init(read_climoz)
     call phys_output_var_init
     print*, '================================================='
     !
     dnwd0=0.0
     ftd=0.0
     fqd=0.0
     cin=0.
     !ym Attention pbase pas initialise dans concvl !!!!
     pbase=0
     !IM 180608

     itau_con=0
     first=.false.

  endif  ! first

  !ym => necessaire pour iflag_con != 2    
  pmfd(:,:) = 0.
  pen_u(:,:) = 0.
  pen_d(:,:) = 0.
  pde_d(:,:) = 0.
  pde_u(:,:) = 0.
  aam=0.

  torsfc=0.
  forall (k=1: llm) zmasse(:, k) = (paprs(:, k)-paprs(:, k+1)) / rg


  modname = 'physiq'
  !IM
  IF (ip_ebil_phy.ge.1) THEN
     DO i=1,klon
        zero_v(i)=0.
     END DO
  END IF

  IF (debut) THEN
     CALL suphel ! initialiser constantes et parametres phys.
  ENDIF

  if(prt_level.ge.1) print*,'CONVERGENCE PHYSIQUE THERM 1 '


  !======================================================================
  ! Gestion calendrier : mise a jour du module phys_cal_mod
  !
  !     CALL phys_cal_update(jD_cur,jH_cur)

  !
  ! Si c'est le debut, il faut initialiser plusieurs choses
  !          ********
  !
  IF (debut) THEN
     !rv
     !CRinitialisation de wght_th et lalim_conv pour la definition de la couche alimentation 
     !de la convection a partir des caracteristiques du thermique
     wght_th(:,:)=1.
     lalim_conv(:)=1 
     !RC
     ustar(:,:)=0.
     u10m(:,:)=0.
     v10m(:,:)=0.
     rain_con(:)=0.
     snow_con(:)=0.
     topswai(:)=0.
     topswad(:)=0.
     solswai(:)=0.
     solswad(:)=0.

     wmax_th(:)=0.
     tau_overturning_th(:)=0.

     IF (type_trac == 'inca') THEN
        ! jg : initialisation jusqu'au ces variables sont dans restart
        ccm(:,:,:) = 0.
        tau_aero(:,:,:,:) = 0.
        piz_aero(:,:,:,:) = 0.
        cg_aero(:,:,:,:) = 0.
     END IF

     rnebcon0(:,:) = 0.0
     clwcon0(:,:) = 0.0
     rnebcon(:,:) = 0.0
     clwcon(:,:) = 0.0

     !IM      
     IF (ip_ebil_phy.ge.1) d_h_vcol_phy=0.
     !
     print*,'iflag_coupl,iflag_clos,iflag_wake', &
          iflag_coupl,iflag_clos,iflag_wake
     print*,'CYCLE_DIURNE', cycle_diurne
     !
     IF (iflag_con.EQ.2.AND.iflag_cldcon.GT.-1) THEN
        abort_message = 'Tiedtke needs iflag_cldcon=-2 or -1'
        CALL abort_gcm (modname,abort_message,1)
     ENDIF
     !
     !
     ! Initialiser les compteurs:
     !
     itap    = 0
     itaprad = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !! Un petit travail \`a faire ici.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (iflag_pbl>1) then
        PRINT*, "Using method MELLOR&YAMADA" 
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! FH 2008/05/02 changement lie a la lecture de nbapp_rad dans phylmd plutot que
     ! dyn3d
     ! Attention : la version precedente n'etait pas tres propre.
     ! Il se peut qu'il faille prendre une valeur differente de nbapp_rad
     ! pour obtenir le meme resultat.
     dtime=pdtphys
     radpas = NINT( 86400./dtime/nbapp_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     CALL phyetat0 ("startphy.nc",clesphy0,tabcntr0)
     IF (klon_glo==1) THEN
        coefh=0. ; coefm=0. ; pbl_tke=0.
        coefh(:,2,:)=1.e-2 ; coefm(:,2,:)=1.e-2 ; pbl_tke(:,2,:)=1.e-2
        PRINT*,'FH WARNING : lignes a supprimer'
     ENDIF
     !IM begin
     print*,'physiq: clwcon rnebcon ratqs',clwcon(1,1),rnebcon(1,1) &
          ,ratqs(1,1)
     !IM end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     ! on remet le calendrier a zero
     !
     IF (raz_date .eq. 1) THEN
        itau_phy = 0
     ENDIF

     !IM cf. AM 081204 BEG
     PRINT*,'cycle_diurne3 =',cycle_diurne
     !IM cf. AM 081204 END
     !
     CALL printflag( tabcntr0,radpas,ok_journe, &
          ok_instan, ok_region )
     !
     IF (ABS(dtime-pdtphys).GT.0.001) THEN
        WRITE(lunout,*) 'Pas physique n est pas correct',dtime, &
             pdtphys
        abort_message='Pas physique n est pas correct '
        !           call abort_gcm(modname,abort_message,1)
        dtime=pdtphys
     ENDIF
     IF (nlon .NE. klon) THEN
        WRITE(lunout,*)'nlon et klon ne sont pas coherents', nlon,  &
             klon
        abort_message='nlon et klon ne sont pas coherents'
        call abort_gcm(modname,abort_message,1)
     ENDIF
     IF (nlev .NE. klev) THEN
        WRITE(lunout,*)'nlev et klev ne sont pas coherents', nlev, &
             klev
        abort_message='nlev et klev ne sont pas coherents'
        call abort_gcm(modname,abort_message,1)
     ENDIF
     !
     IF (dtime*REAL(radpas).GT.21600..AND.cycle_diurne) THEN 
        WRITE(lunout,*)'Nbre d appels au rayonnement insuffisant'
        WRITE(lunout,*)"Au minimum 4 appels par jour si cycle diurne"
        abort_message='Nbre d appels au rayonnement insuffisant'
        call abort_gcm(modname,abort_message,1)
     ENDIF
     WRITE(lunout,*)"Clef pour la convection, iflag_con=", iflag_con
     WRITE(lunout,*)"Clef pour le driver de la convection, ok_cvl=", &
          ok_cvl
     !
     !KE43
     ! Initialisation pour la convection de K.E. (sb):
     IF (iflag_con.GE.3) THEN

        WRITE(lunout,*)"*** Convection de Kerry Emanuel 4.3  "
        WRITE(lunout,*) &
             "On va utiliser le melange convectif des traceurs qui"
        WRITE(lunout,*)"est calcule dans convect4.3"
        WRITE(lunout,*)" !!! penser aux logical flags de phytrac"

        DO i = 1, klon
           ema_cbmf(i) = 0.
           ema_pcb(i)  = 0.
           ema_pct(i)  = 0.
           !          ema_workcbmf(i) = 0.
        ENDDO
        !IM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>BEG
        DO i = 1, klon
           ibas_con(i) = 1
           itop_con(i) = 1
        ENDDO
        !IM15/11/02 rajout initialisation ibas_con,itop_con cf. SB =>END
        !===============================================================================
        !CR:04.12.07: initialisations poches froides
        ! Controle de ALE et ALP pour la fermeture convective (jyg)
        if (iflag_wake>=1) then
           CALL ini_wake(0.,0.,it_wape_prescr,wape_prescr,fip_prescr &
                ,alp_bl_prescr, ale_bl_prescr)
           ! 11/09/06 rajout initialisation ALE et ALP du wake et PBL(YU)
           !        print*,'apres ini_wake iflag_cldcon=', iflag_cldcon
        endif

!        do i = 1,klon
!           Ale_bl(i)=0.
!           Alp_bl(i)=0.
!        enddo

        !================================================================================
        !IM stations CFMIP
        nCFMIP=npCFMIP
        OPEN(98,file='npCFMIP_param.data',status='old', &
             form='formatted',iostat=iostat)
        if (iostat == 0) then
           READ(98,*,end=998) nCFMIP
998        CONTINUE
           CLOSE(98)
           CONTINUE
           IF(nCFMIP.GT.npCFMIP) THEN
              print*,'nCFMIP > npCFMIP : augmenter npCFMIP et recompiler'
              call abort_gcm("physiq", "", 1)
           else
              print*,'physiq npCFMIP=',npCFMIP,'nCFMIP=',nCFMIP
           ENDIF

           !
           ALLOCATE(tabCFMIP(nCFMIP))
           ALLOCATE(lonCFMIP(nCFMIP), latCFMIP(nCFMIP))
           ALLOCATE(tabijGCM(nCFMIP))
           ALLOCATE(lonGCM(nCFMIP), latGCM(nCFMIP))
           ALLOCATE(iGCM(nCFMIP), jGCM(nCFMIP))
           !
           ! lecture des nCFMIP stations CFMIP, de leur numero 
           ! et des coordonnees geographiques lonCFMIP, latCFMIP
           !
           CALL read_CFMIP_point_locations(nCFMIP, tabCFMIP,  &
                lonCFMIP, latCFMIP)
           !
           ! identification des
           ! 1) coordonnees lonGCM, latGCM des points CFMIP dans la grille de LMDZ
           ! 2) indices points tabijGCM de la grille physique 1d sur klon points
           ! 3) indices iGCM, jGCM de la grille physique 2d
           !
           CALL LMDZ_CFMIP_point_locations(nCFMIP, lonCFMIP, latCFMIP, &
                tabijGCM, lonGCM, latGCM, iGCM, jGCM)
           !
        else
           ALLOCATE(tabijGCM(0))
           ALLOCATE(lonGCM(0), latGCM(0))
           ALLOCATE(iGCM(0), jGCM(0))
        end if
     else
        ALLOCATE(tabijGCM(0))
        ALLOCATE(lonGCM(0), latGCM(0))
        ALLOCATE(iGCM(0), jGCM(0))
     ENDIF

     DO i=1,klon
        rugoro(i) = f_rugoro * MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
     ENDDO

     !34EK
     IF (ok_orodr) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! FH sans doute a enlever de finitivement ou, si on le garde, l'activer
        ! justement quand ok_orodr = false.
        ! ce rugoro est utilise par la couche limite et fait double emploi
        ! avec les param\'etrisations sp\'ecifiques de Francois Lott.
        !           DO i=1,klon
        !             rugoro(i) = MAX(1.0e-05, zstd(i)*zsig(i)/2.0)
        !           ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (ok_strato) THEN
           CALL SUGWD_strato(klon,klev,paprs,pplay)
        ELSE
           CALL SUGWD(klon,klev,paprs,pplay)
        ENDIF

        DO i=1,klon
           zuthe(i)=0.
           zvthe(i)=0.
           if(zstd(i).gt.10.)then
              zuthe(i)=(1.-zgam(i))*cos(zthe(i))
              zvthe(i)=(1.-zgam(i))*sin(zthe(i))
           endif
        ENDDO
     ENDIF
     !
     !
     lmt_pas = NINT(86400./dtime * 1.0)   ! tous les jours
     WRITE(lunout,*)'La frequence de lecture surface est de ',  &
          lmt_pas
     !
     capemaxcels = 't_max(X)'
     t2mincels = 't_min(X)'
     t2maxcels = 't_max(X)'
     tinst = 'inst(X)'
     tave = 'ave(X)'
     !IM cf. AM 081204 BEG
     write(lunout,*)'AVANT HIST IFLAG_CON=',iflag_con
     !IM cf. AM 081204 END
     !
     !=============================================================
     !   Initialisation des sorties
     !=============================================================

#ifdef CPP_IOIPSL

     !$OMP MASTER
! FH : if ok_sync=.true. , the time axis is written at each time step
! in the output files. Only at the end in the opposite case
     ok_sync_omp=.false.
     CALL getin('ok_sync',ok_sync_omp)
     call phys_output_open(rlon,rlat,nCFMIP,tabijGCM, &
          iGCM,jGCM,lonGCM,latGCM, &
          jjmp1,nlevSTD,clevSTD,rlevSTD, dtime,ok_veget, &
          type_ocean,iflag_pbl,iflag_pbl_split,ok_mensuel,ok_journe, &
          ok_hf,ok_instan,ok_LES,ok_ade,ok_aie,  &
          read_climoz, phys_out_filestations, &
          new_aod, aerosol_couple, &
          flag_aerosol_strat, pdtphys, paprs, pphis,  &
          pplay, lmax_th, ptconv, ptconvth, ivap,  &
          d_t, qx, d_qx, zmasse, ok_sync_omp)
     !$OMP END MASTER
     !$OMP BARRIER
     ok_sync=ok_sync_omp

     freq_outNMC(1) = ecrit_files(7)
     freq_outNMC(2) = ecrit_files(8)
     freq_outNMC(3) = ecrit_files(9)
     WRITE(lunout,*)'OK freq_outNMC(1)=',freq_outNMC(1)
     WRITE(lunout,*)'OK freq_outNMC(2)=',freq_outNMC(2)
     WRITE(lunout,*)'OK freq_outNMC(3)=',freq_outNMC(3)

     include "ini_histday_seri.h"

     include "ini_paramLMDZ_phy.h"

#endif
     ecrit_reg = ecrit_reg * un_jour
     ecrit_tra = ecrit_tra * un_jour

     !XXXPB Positionner date0 pour initialisation de ORCHIDEE
     date0 = jD_ref 
     WRITE(*,*) 'physiq date0 : ',date0
     !
     !
     !
     ! Prescrire l'ozone dans l'atmosphere
     !
     !
     !c         DO i = 1, klon
     !c         DO k = 1, klev
     !c            CALL o3cm (paprs(i,k)/100.,paprs(i,k+1)/100., wo(i,k),20)
     !c         ENDDO
     !c         ENDDO
     !
     IF (type_trac == 'inca') THEN
#ifdef INCA
        CALL VTe(VTphysiq)
        CALL VTb(VTinca)
        calday = REAL(days_elapsed) + jH_cur
        WRITE(lunout,*) 'initial time chemini', days_elapsed, calday

        CALL chemini(  &
             rg, &
             ra, &
             airephy, &
             rlat, &
             rlon, &
             presnivs, &
             calday, &
             klon, &
             nqtot, &
             pdtphys, &
             annee_ref, &
             day_ref,  &
             itau_phy)

        CALL VTe(VTinca)
        CALL VTb(VTphysiq)
#endif
     END IF
     !
     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Nouvelle initialisation pour le rayonnement RRTM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call iniradia(klon,klev,paprs(1,1:klev+1))

     !$omp single
     if (read_climoz >= 1) then
        call open_climoz(ncid_climoz, press_climoz)
     END IF
     !$omp end single
     !
     !IM betaCRF
     pfree=70000. !Pa
     beta_pbl=1.
     beta_free=1.
     lon1_beta=-180.
     lon2_beta=+180.
     lat1_beta=90.
     lat2_beta=-90.
     mskocean_beta=.FALSE.

     OPEN(99,file='beta_crf.data',status='old', &
          form='formatted',err=9999)
     READ(99,*,end=9998) pfree
     READ(99,*,end=9998) beta_pbl
     READ(99,*,end=9998) beta_free
     READ(99,*,end=9998) lon1_beta
     READ(99,*,end=9998) lon2_beta
     READ(99,*,end=9998) lat1_beta
     READ(99,*,end=9998) lat2_beta
     READ(99,*,end=9998) mskocean_beta
9998 Continue
     CLOSE(99)
9999 Continue
     WRITE(*,*)'pfree=',pfree
     WRITE(*,*)'beta_pbl=',beta_pbl
     WRITE(*,*)'beta_free=',beta_free
     WRITE(*,*)'lon1_beta=',lon1_beta
     WRITE(*,*)'lon2_beta=',lon2_beta
     WRITE(*,*)'lat1_beta=',lat1_beta
     WRITE(*,*)'lat2_beta=',lat2_beta
     WRITE(*,*)'mskocean_beta=',mskocean_beta
  ENDIF
  !
  !   ****************     Fin  de   IF ( debut  )   ***************
  !
  !
  ! Incrementer le compteur de la physique
  !
  itap   = itap + 1
  !
  !
  ! Update fraction of the sub-surfaces (pctsrf) and 
  ! initialize, where a new fraction has appeared, all variables depending 
  ! on the surface fraction.
  !
  CALL change_srf_frac(itap, dtime, days_elapsed+1,  &
       pctsrf, falb1, falb2, ftsol, ustar, u10m, v10m, pbl_tke)


  ! Update time and other variables in Reprobus
  IF (type_trac == 'repr') THEN
#ifdef REPROBUS
     CALL Init_chem_rep_xjour(jD_cur-jD_ref+day_ref)
     print*,'xjour equivalent rjourvrai',jD_cur-jD_ref+day_ref
     CALL Rtime(debut)
#endif
  END IF


  ! Tendances bidons pour les processus qui n'affectent pas certaines
  ! variables.
  du0(:,:)=0.
  dv0(:,:)=0.
  dt0 = 0.
  dq0(:,:)=0.
  dql0(:,:)=0.
  dqi0(:,:)=0.
  !
  ! Mettre a zero des variables de sortie (pour securite)
  !
  DO i = 1, klon
     d_ps(i) = 0.0
  ENDDO
  DO k = 1, klev
     DO i = 1, klon
        d_t(i,k) = 0.0
        d_u(i,k) = 0.0
        d_v(i,k) = 0.0
     ENDDO
  ENDDO
  DO iq = 1, nqtot
     DO k = 1, klev
        DO i = 1, klon
           d_qx(i,k,iq) = 0.0
        ENDDO
     ENDDO
  ENDDO
  da(:,:)=0.
  mp(:,:)=0.
  phi(:,:,:)=0.
  ! RomP >>>
  phi2(:,:,:)=0.
  beta_prec_fisrt(:,:)=0.
  beta_prec(:,:)=0.
  epmlmMm(:,:,:)=0.
  eplaMm(:,:)=0.
  d1a(:,:)=0.
  dam(:,:)=0.
  pmflxr=0.
  pmflxs=0.
  ! RomP <<<

  !
  ! Ne pas affecter les valeurs entrees de u, v, h, et q
  !
  DO k = 1, klev
     DO i = 1, klon
        t_seri(i,k)  = t(i,k)
        u_seri(i,k)  = u(i,k)
        v_seri(i,k)  = v(i,k)
        q_seri(i,k)  = qx(i,k,ivap)
        ql_seri(i,k) = qx(i,k,iliq)
!CR: ATTENTION, on rajoute la variable glace
        if (nqo.eq.2) then
           qs_seri(i,k) = 0.
        else if (nqo.eq.3) then
           qs_seri(i,k) = qx(i,k,isol)
        endif
     ENDDO
  ENDDO
  tke0(:,:)=pbl_tke(:,:,is_ave)
!CR:Nombre de traceurs de l'eau: nqo
!  IF (nqtot.GE.3) THEN
   IF (nqtot.GE.(nqo+1)) THEN
!     DO iq = 3, nqtot        
     DO iq = nqo+1, nqtot  
        DO  k = 1, klev
           DO  i = 1, klon
!              tr_seri(i,k,iq-2) = qx(i,k,iq)
              tr_seri(i,k,iq-nqo) = qx(i,k,iq)
           ENDDO
        ENDDO
     ENDDO
  ELSE
     DO k = 1, klev
        DO i = 1, klon
           tr_seri(i,k,1) = 0.0
        ENDDO
     ENDDO
  ENDIF
  !
  DO i = 1, klon
     ztsol(i) = 0.
  ENDDO
  DO nsrf = 1, nbsrf
     DO i = 1, klon
        ztsol(i) = ztsol(i) + ftsol(i,nsrf)*pctsrf(i,nsrf)
     ENDDO
  ENDDO
  !IM
  IF (ip_ebil_phy.ge.1) THEN 
     ztit='after dynamic'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,1,1,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     !     Comme les tendances de la physique sont ajoute dans la dynamique,
     !     on devrait avoir que la variation d'entalpie par la dynamique
     !     est egale a la variation de la physique au pas de temps precedent.
     !     Donc la somme de ces 2 variations devrait etre nulle.
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol+d_h_vcol_phy, d_qt, 0. &
          , fs_bound, fq_bound )
  END IF

  ! Diagnostiquer la tendance dynamique
  !
  IF (ancien_ok) THEN
     DO k = 1, klev
        DO i = 1, klon
           d_u_dyn(i,k) = (u_seri(i,k)-u_ancien(i,k))/dtime
           d_v_dyn(i,k) = (v_seri(i,k)-v_ancien(i,k))/dtime
           d_t_dyn(i,k) = (t_seri(i,k)-t_ancien(i,k))/dtime
           d_q_dyn(i,k) = (q_seri(i,k)-q_ancien(i,k))/dtime
        ENDDO
     ENDDO
!!! RomP >>>   td dyn traceur
     IF (nqtot.GE.3) THEN
        DO iq = 3, nqtot
           DO k = 1, klev
              DO i = 1, klon
                 d_tr_dyn(i,k,iq-2)= &
                      (tr_seri(i,k,iq-2)-tr_ancien(i,k,iq-2))/dtime
                 !         iiq=niadv(iq)
                 !         print*,i,k," d_tr_dyn",d_tr_dyn(i,k,iq-2),"tra:",iq,tname(iiq)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
!!! RomP <<<
  ELSE
     DO k = 1, klev
        DO i = 1, klon
           d_u_dyn(i,k) = 0.0
           d_v_dyn(i,k) = 0.0
           d_t_dyn(i,k) = 0.0
           d_q_dyn(i,k) = 0.0
        ENDDO
     ENDDO
!!! RomP >>>   td dyn traceur
     IF (nqtot.GE.3) THEN
        DO iq = 3, nqtot
           DO k = 1, klev
              DO i = 1, klon
                 d_tr_dyn(i,k,iq-2)= 0.0
              ENDDO
           ENDDO
        ENDDO
     ENDIF
!!! RomP <<<
     ancien_ok = .TRUE.
  ENDIF
  !
  ! Ajouter le geopotentiel du sol:
  !
  DO k = 1, klev
     DO i = 1, klon
        zphi(i,k) = pphi(i,k) + pphis(i)
     ENDDO
  ENDDO
  !
  ! Verifier les temperatures
  !
  !IM BEG
  IF (check) THEN
     amn=MIN(ftsol(1,is_ter),1000.)
     amx=MAX(ftsol(1,is_ter),-1000.)
     DO i=2, klon
        amn=MIN(ftsol(i,is_ter),amn)
        amx=MAX(ftsol(i,is_ter),amx)
     ENDDO
     !
     PRINT*,' debut avant hgardfou min max ftsol',itap,amn,amx
  ENDIF !(check) THEN
  !IM END
  !
  CALL hgardfou(t_seri,ftsol,'debutphy')
  !
  !IM BEG
  IF (check) THEN
     amn=MIN(ftsol(1,is_ter),1000.)
     amx=MAX(ftsol(1,is_ter),-1000.)
     DO i=2, klon
        amn=MIN(ftsol(i,is_ter),amn)
        amx=MAX(ftsol(i,is_ter),amx)
     ENDDO
     !
     PRINT*,' debut apres hgardfou min max ftsol',itap,amn,amx
  ENDIF !(check) THEN
  !IM END
  !
  ! Mettre en action les conditions aux limites (albedo, sst, etc.).
  ! Prescrire l'ozone et calculer l'albedo sur l'ocean.
  !
  if (read_climoz >= 1) then
     ! Ozone from a file
     ! Update required ozone index:
     ro3i = int((days_elapsed + jh_cur - jh_1jan) / year_len * 360.) + 1
     if (ro3i == 361) ro3i = 360
     ! (This should never occur, except perhaps because of roundup
     ! error. See documentation.)
     if (ro3i /= co3i) then
        ! Update ozone field:
        if (read_climoz == 1) then
           call regr_pr_av(ncid_climoz, (/"tro3"/), julien=ro3i, &
                press_in_edg=press_climoz, paprs=paprs, v3=wo)
        else
           ! read_climoz == 2
           call regr_pr_av(ncid_climoz, (/"tro3         ", "tro3_daylight"/), &
                julien=ro3i, press_in_edg=press_climoz, paprs=paprs, v3=wo)
        end if
        ! Convert from mole fraction of ozone to column density of ozone in a
        ! cell, in kDU:
        forall (l = 1: read_climoz) wo(:, :, l) = wo(:, :, l) * rmo3 / rmd &
             * zmasse / dobson_u / 1e3
        ! (By regridding ozone values for LMDZ only once every 360th of
        ! year, we have already neglected the variation of pressure in one
        ! 360th of year. So do not recompute "wo" at each time step even if
        ! "zmasse" changes a little.)
        co3i = ro3i
     end if

!****************************** DART-LMDZ *********************************************
!    ELSEIF (MOD(itap-1,lmt_pas) == 0) THEN  
     ELSEIF (MOD(itap-1 + itau_phy,lmt_pas) == 0) THEN   
     ! Once per day, update ozone from Royer:

     IF (solarlong0<-999.) then
        ! Generic case with evolvoing season
        zzz=real(days_elapsed+1)
     ELSE IF (abs(solarlong0-1000.)<1.e-4) then
        ! Particular case with annual mean insolation
        zzz=real(90) ! could be revisited
        IF (read_climoz/=-1) THEN
           abort_message ='read_climoz=-1 is recommended when solarlong0=1000.'
           CALL abort_gcm (modname,abort_message,1)
        ENDIF
     ELSE
        ! Case where the season is imposed with solarlong0
        zzz=real(90) ! could be revisited
     ENDIF
     wo(:,:,1)=ozonecm(rlat, paprs,read_climoz,rjour=zzz)
      CALL gather(paprs, paprs_glo)
!$OMP MASTER
     IF (is_mpi_root .AND. is_omp_root) THEN 
      open(155,file="stok_paprs.dat",form='unformatted')      
      write(155)zzz  ! save for next sub day cycles to compute wo    
      write(155)read_climoz  ! save for next sub day cycles to compute wo   
      write(155)paprs_glo  ! save for next sub day cycles to compute wo    
      close(155) 
     ENDIF                                                  
!$OMP END MASTER
    elseif (MOD(itap-1,lmt_pas) == 0) THEN  !! For next sub day cycles 
!$OMP MASTER
     IF (is_mpi_root .AND. is_omp_root) THEN 
     open(156,file="stok_paprs.dat",form='unformatted')               
     read(156)zzz                                                    
     read(156)read_climoz                                           
     read(156)paprs_glo                                            
     ENDIF                                                   
!$OMP END MASTER
     call bcast(zzz)
     call bcast(read_climoz)
     CALL scatter(paprs_glo,paprs1) 
      wo(:,:,1)=ozonecm(rlat, paprs1,read_climoz,rjour=zzz) 
!$OMP MASTER
     IF (is_mpi_root) THEN    
     close(156)              
    ENDIF
!$OMP END MASTER
  ENDIF

!****************************** DART-LMDZ *********************************************



  !
  ! Re-evaporer l'eau liquide nuageuse
  !
  DO k = 1, klev  ! re-evaporation de l'eau liquide nuageuse
     DO i = 1, klon
        zlvdcp=RLVTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
        !jyg<
        !      Attention : Arnaud a propose des formules completement differentes
        !                  A verifier !!!
        zlsdcp=RLSTT/RCPD/(1.0+RVTMP2*q_seri(i,k))
        IF (iflag_ice_thermo .EQ. 0) THEN
           zlsdcp=zlvdcp
        ENDIF
        !>jyg
     
        if (iflag_ice_thermo.eq.0) then   
!pas necessaire a priori

        zdelta = MAX(0.,SIGN(1.,RTT-t_seri(i,k)))
        zb = MAX(0.0,ql_seri(i,k))
        za = - MAX(0.0,ql_seri(i,k)) &
             * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
        t_seri(i,k) = t_seri(i,k) + za
        q_seri(i,k) = q_seri(i,k) + zb
        ql_seri(i,k) = 0.0
        d_t_eva(i,k) = za
        d_q_eva(i,k) = zb

        else

!CR: on ré-évapore eau liquide et glace

!        zdelta = MAX(0.,SIGN(1.,RTT-t_seri(i,k)))
!        zb = MAX(0.0,ql_seri(i,k))
!        za = - MAX(0.0,ql_seri(i,k)) &
!             * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
        zb = MAX(0.0,ql_seri(i,k)+qs_seri(i,k))
        za = - MAX(0.0,ql_seri(i,k))*zlvdcp & 
             - MAX(0.0,qs_seri(i,k))*zlsdcp
        t_seri(i,k) = t_seri(i,k) + za
        q_seri(i,k) = q_seri(i,k) + zb
        ql_seri(i,k) = 0.0
!on évapore la glace
        qs_seri(i,k) = 0.0
        d_t_eva(i,k) = za
        d_q_eva(i,k) = zb
        endif

     ENDDO
  ENDDO
  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after reevap'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,1,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
     !
  END IF

  !
  !=========================================================================
  ! Calculs de l'orbite.
  ! Necessaires pour le rayonnement et la surface (calcul de l'albedo).
  ! doit donc etre plac\'e avant radlwsw et pbl_surface

!!!   jyg 17 Sep 2010 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call ymds2ju(year_cur, mth_eq, day_eq,0., jD_eq)
  day_since_equinox = (jD_cur + jH_cur) - jD_eq
  !
  !   choix entre calcul de la longitude solaire vraie ou valeur fixee a 
  !   solarlong0
  if (solarlong0<-999.) then
     if (new_orbit) then
        ! calcul selon la routine utilisee pour les planetes
        call solarlong(day_since_equinox, zlongi, dist)
     else
        ! calcul selon la routine utilisee pour l'AR4
        CALL orbite(REAL(days_elapsed+1),zlongi,dist)
     endif
  else
     zlongi=solarlong0  ! longitude solaire vraie
     dist=1.            ! distance au soleil / moyenne 
  endif
  if(prt_level.ge.1)                                                &
       write(lunout,*)'Longitude solaire ',zlongi,solarlong0,dist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calcul de l'ensoleillement :
  ! ============================
  ! Pour une solarlong0=1000., on calcule un ensoleillement moyen sur
  ! l'annee a partir d'une formule analytique.
  ! Cet ensoleillement est sym\'etrique autour de l'\'equateur et
  ! non nul aux poles.
  IF (abs(solarlong0-1000.)<1.e-4) then
     call zenang_an(cycle_diurne,jH_cur,rlat,rlon,rmu0,fract)
  ELSE
     !  Avec ou sans cycle diurne
     IF (cycle_diurne) THEN
        zdtime=dtime*REAL(radpas) ! pas de temps du rayonnement (s)
        CALL zenang(zlongi,jH_cur,zdtime,rlat,rlon,rmu0,fract)
     ELSE
        CALL angle(zlongi, rlat, fract, rmu0)
     ENDIF
  ENDIF

  ! AI Janv 2014 
  do i = 1, klon
     if (fract(i).le.0.) then
        JrNt(i)=0.
     else
        JrNt(i)=1.
     endif
  enddo

  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Appel au pbl_surface : Planetary Boudary Layer et Surface
  ! Cela implique tous les interactions des sous-surfaces et la partie diffusion 
  ! turbulent du couche limit. 
  ! 
  ! Certains varibales de sorties de pbl_surface sont utiliser que pour 
  ! ecriture des fihiers hist_XXXX.nc, ces sont :
  !   qsol,      zq2m,      s_pblh,  s_lcl,
  !   s_capCL,   s_oliqCL,  s_cteiCL,s_pblT,
  !   s_therm,   s_trmb1,   s_trmb2, s_trmb3,
  !   zxrugs,    zu10m,     zv10m,   fder,
  !   zxqsurf,   rh2m,      zxfluxu, zxfluxv,
  !   frugs,     agesno,    fsollw,  fsolsw,
  !   d_ts,      fevap,     fluxlat, t2m,
  !   wfbils,    wfbilo,    fluxt,   fluxu, fluxv,
  !
  ! Certains ne sont pas utiliser du tout : 
  !   dsens, devap, zxsnow, zxfluxt, zxfluxq, q2m, fluxq
  !

  ! Calcul de l'humidite de saturation au niveau du sol



  if (iflag_pbl/=0) then

!jyg+nrlmd<
      IF (prt_level .ge. 2 .and. mod(iflag_pbl_split,2) .eq. 1) THEN
        print *,'debut du splitting de la PBL'
      ENDIF
!!!
!=================================================================
!         PROVISOIRE : DECOUPLAGE PBL/WAKE
!         --------------------------------
!
!!      wake_deltat_sav(:,:)=wake_deltat(:,:)
!!      wake_deltaq_sav(:,:)=wake_deltaq(:,:)
!!      wake_deltat(:,:)=0.
!!      wake_deltaq(:,:)=0.
!=================================================================
!>jyg+nrlmd
!
     CALL pbl_surface(  &
          dtime,     date0,     itap,    days_elapsed+1, &
          debut,     lafin, &
          rlon,      rlat,      rugoro,  rmu0,      &
          zsig,      sollwdown, pphi,    cldt,      &
          rain_fall, snow_fall, solsw,   sollw,     &
          t_seri,    q_seri,    u_seri,  v_seri,    &
!nrlmd+jyg<
          wake_deltat, wake_deltaq, wake_cstar, wake_s, &
!>nrlmd+jyg
          pplay,     paprs,     pctsrf,             &
          ftsol,falb1,falb2,ustar,u10m,v10m,wstar, &
          sollwdown, cdragh,    cdragm,  u1,    v1, &
          albsol1,   albsol2,   sens,    evap,   &
          albsol3_lic,runoff,   snowhgt,   qsnow, to_ice, sissnow, &
          zxtsol,    zxfluxlat, zt2m,    qsat2m,  &
          d_t_vdf,   d_q_vdf,   d_u_vdf, d_v_vdf, d_t_diss, &
!nrlmd<
  !jyg<
          d_t_vdf_w, d_q_vdf_w, &
          d_t_vdf_x, d_q_vdf_x, &
          sens_x, zxfluxlat_x, sens_w, zxfluxlat_w, &
  !>jyg
          delta_tsurf,wake_dens, &
          cdragh_x,cdragh_w,cdragm_x,cdragm_w, &
          kh,kh_x,kh_w, &
!>nrlmd
          coefh(1:klon,1:klev,1:nbsrf+1),     coefm(1:klon,1:klev,1:nbsrf+1), &
          slab_wfbils,                 &
          qsol,      zq2m,      s_pblh,  s_lcl, &
!jyg<
          s_pblh_x, s_lcl_x, s_pblh_w, s_lcl_w, &
!>jyg
          s_capCL,   s_oliqCL,  s_cteiCL,s_pblT, &
          s_therm,   s_trmb1,   s_trmb2, s_trmb3, &
          zxrugs,    zustar, zu10m,     zv10m,   fder, &
          zxqsurf,   rh2m,      zxfluxu, zxfluxv, &
          frugs,     agesno,    fsollw,  fsolsw, &
          d_ts,      fevap,     fluxlat, t2m, &
          wfbils,    wfbilo,    fluxt,   fluxu,  fluxv, &
          dsens,     devap,     zxsnow, &
          zxfluxt,   zxfluxq,   q2m,     fluxq, pbl_tke, &
!nrlmd+jyg<
          wake_delta_pbl_TKE &
!>nrlmd+jyg
                      )
!
!=================================================================
!         PROVISOIRE : DECOUPLAGE PBL/WAKE
!         --------------------------------
!
!!      wake_deltat(:,:)=wake_deltat_sav(:,:)
!!      wake_deltaq(:,:)=wake_deltaq_sav(:,:)
!=================================================================
!
!  Add turbulent diffusion tendency to the wake difference variables
    wake_deltat(:,:) = wake_deltat(:,:) + (d_t_vdf_w(:,:)-d_t_vdf_x(:,:))
    wake_deltaq(:,:) = wake_deltaq(:,:) + (d_q_vdf_w(:,:)-d_q_vdf_x(:,:))


     !---------------------------------------------------------------------
     ! ajout des tendances de la diffusion turbulente
     IF (klon_glo==1) THEN
        CALL add_pbl_tend &
             (d_u_vdf,d_v_vdf,d_t_vdf+d_t_diss,d_q_vdf,dql0,dqi0,paprs,'vdf')
     ELSE
        CALL add_phys_tend &
             (d_u_vdf,d_v_vdf,d_t_vdf+d_t_diss,d_q_vdf,dql0,dqi0,paprs,'vdf')
     ENDIF
     !--------------------------------------------------------------------

     if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
     endif

     CALL evappot(klon,nbsrf,ftsol,pplay(:,1),cdragh, &
          t_seri(:,1),q_seri(:,1),u_seri(:,1),v_seri(:,1),evap_pot)


     IF (ip_ebil_phy.ge.2) THEN 
        ztit='after surface_main'
        CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
             , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
             , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
        call diagphy(airephy,ztit,ip_ebil_phy &
             , zero_v, zero_v, zero_v, zero_v, sens &
             , evap  , zero_v, zero_v, ztsol &
             , d_h_vcol, d_qt, d_ec &
             , fs_bound, fq_bound )
     END IF

  ENDIF
  ! =================================================================== c
  !   Calcul de Qsat

  DO k = 1, klev
     DO i = 1, klon
        zx_t = t_seri(i,k)
        IF (thermcep) THEN
           zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
           zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
           zx_qs  = MIN(0.5,zx_qs)
           zcor   = 1./(1.-retv*zx_qs)
           zx_qs  = zx_qs*zcor
        ELSE
           IF (zx_t.LT.t_coup) THEN
              zx_qs = qsats(zx_t)/pplay(i,k)
           ELSE
              zx_qs = qsatl(zx_t)/pplay(i,k)
           ENDIF
        ENDIF
        zqsat(i,k)=zx_qs
     ENDDO
  ENDDO

  if (prt_level.ge.1) then
     write(lunout,*) 'L   qsat (g/kg) avant clouds_gno'
     write(lunout,'(i4,f15.4)') (k,1000.*zqsat(igout,k),k=1,klev)
  endif
  !
  ! Appeler la convection (au choix)
  !
  DO k = 1, klev
     DO i = 1, klon
        conv_q(i,k) = d_q_dyn(i,k)  &
             + d_q_vdf(i,k)/dtime
        conv_t(i,k) = d_t_dyn(i,k)  &
             + d_t_vdf(i,k)/dtime
     ENDDO
  ENDDO
  IF (check) THEN
     za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
     WRITE(lunout,*) "avantcon=", za
  ENDIF
  zx_ajustq = .FALSE.
  IF (iflag_con.EQ.2) zx_ajustq=.TRUE.
  IF (zx_ajustq) THEN
     DO i = 1, klon
        z_avant(i) = 0.0
     ENDDO
     DO k = 1, klev
        DO i = 1, klon
           z_avant(i) = z_avant(i) + (q_seri(i,k)+ql_seri(i,k)) &
                *(paprs(i,k)-paprs(i,k+1))/RG
        ENDDO
     ENDDO
  ENDIF

  ! Calcule de vitesse verticale a partir de flux de masse verticale
  DO k = 1, klev
     DO i = 1, klon
        omega(i,k) = RG*flxmass_w(i,k) / airephy(i)
     END DO
  END DO
  if (prt_level.ge.1) write(lunout,*) 'omega(igout, :) = ', &
       omega(igout, :)

  IF (iflag_con.EQ.1) THEN
     abort_message ='reactiver le call conlmd dans physiq.F'
     CALL abort_gcm (modname,abort_message,1)
     !     CALL conlmd (dtime, paprs, pplay, t_seri, q_seri, conv_q,
     !    .             d_t_con, d_q_con,
     !    .             rain_con, snow_con, ibas_con, itop_con)
  ELSE IF (iflag_con.EQ.2) THEN
     CALL conflx(dtime, paprs, pplay, t_seri, q_seri, &
          conv_t, conv_q, -evap, omega, &
          d_t_con, d_q_con, rain_con, snow_con, &
          pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
          kcbot, kctop, kdtop, pmflxr, pmflxs)
     d_u_con = 0.
     d_v_con = 0.

     WHERE (rain_con < 0.) rain_con = 0.
     WHERE (snow_con < 0.) snow_con = 0.
     DO i = 1, klon
        ibas_con(i) = klev+1 - kcbot(i)
        itop_con(i) = klev+1 - kctop(i)
     ENDDO
  ELSE IF (iflag_con.GE.3) THEN
     ! nb of tracers for the KE convection:
     ! MAF la partie traceurs est faite dans phytrac
     ! on met ntra=1 pour limiter les appels mais on peut
     ! supprimer les calculs / ftra.
     ntra = 1

     !=========================================================================
     !ajout pour la parametrisation des poches froides: calcul de
     !t_wake et t_undi: si pas de poches froides, t_wake=t_undi=t_seri
     do k=1,klev
        do i=1,klon
           if (iflag_wake>=1) then
              t_wake(i,k) = t_seri(i,k) &
                   +(1-wake_s(i))*wake_deltat(i,k)
              q_wake(i,k) = q_seri(i,k) &
                   +(1-wake_s(i))*wake_deltaq(i,k)
              t_undi(i,k) = t_seri(i,k) &
                   -wake_s(i)*wake_deltat(i,k)
              q_undi(i,k) = q_seri(i,k) &
                   -wake_s(i)*wake_deltaq(i,k)
           else
              t_wake(i,k) = t_seri(i,k)
              q_wake(i,k) = q_seri(i,k)
              t_undi(i,k) = t_seri(i,k)
              q_undi(i,k) = q_seri(i,k)
           endif
        enddo
     enddo

     ! Calcul de l'energie disponible ALE (J/kg) et de la puissance
     ! disponible ALP (W/m2) pour le soulevement des particules dans
     ! le modele convectif
     !
     do i = 1,klon
        ALE(i) = 0.
        ALP(i) = 0.
     enddo
     !
     !calcul de ale_wake et alp_wake
     if (iflag_wake>=1) then
        if (itap .le. it_wape_prescr) then
           do i = 1,klon
              ale_wake(i) = wape_prescr
              alp_wake(i) = fip_prescr
           enddo
        else
           do i = 1,klon
              !jyg  ALE=WAPE au lieu de ALE = 1/2 Cstar**2
              !cc           ale_wake(i) = 0.5*wake_cstar(i)**2
              ale_wake(i) = wake_pe(i)
              alp_wake(i) = wake_fip(i)
           enddo
        endif
     else
        do i = 1,klon
           ale_wake(i) = 0.
           alp_wake(i) = 0.
        enddo
     endif
     !combinaison avec ale et alp de couche limite: constantes si pas
     !de couplage, valeurs calculees dans le thermique sinon
     if (iflag_coupl.eq.0) then
        if (debut.and.prt_level.gt.9) &
             WRITE(lunout,*)'ALE et ALP imposes'
        do i = 1,klon
           !on ne couple que ale
           !           ALE(i) = max(ale_wake(i),Ale_bl(i))
           ALE(i) = max(ale_wake(i),ale_bl_prescr)
           !on ne couple que alp
           !           ALP(i) = alp_wake(i) + Alp_bl(i)
           ALP(i) = alp_wake(i) + alp_bl_prescr
        enddo
     else
        IF(prt_level>9)WRITE(lunout,*)'ALE et ALP couples au thermique'
        !         do i = 1,klon
        !             ALE(i) = max(ale_wake(i),Ale_bl(i))
        ! avant        ALP(i) = alp_wake(i) + Alp_bl(i)
        !             ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
        !         write(20,*)'ALE',ALE(i),Ale_bl(i),ale_wake(i)
        !         write(21,*)'ALP',ALP(i),Alp_bl(i),alp_wake(i)
        !         enddo

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Modif FH 2010/04/27. Sans doute temporaire.
        ! Deux options pour le alp_offset : constant si >?? 0 ou
        ! proportionnel ??a w si <0
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimation d'une vitesse verticale effective pour ALP
        www(1:klon)=0.
        do k=2,klev-1
           do i=1,klon
              www(i)=max(www(i),-omega(i,k)*RD*t_seri(i,k)/(RG*paprs(i,k)) &
&                    *zw2(i,k)*zw2(i,k))
!             if (paprs(i,k)>pbase(i)) then
! calcul approche de la vitesse verticale en m/s
!                www(i)=max(www(i),-omega(i,k)*RD*temp(i,k)/(RG*paprs(i,k))
!             endif
!   Le 0.1 est en gros H / ps = 1e5 / 1e4
           enddo
        enddo
        do i=1,klon
           if (www(i)>0. .and. ale_bl(i)>0. ) www(i)=www(i)/ale_bl(i)
        enddo


        do i = 1,klon
           ALE(i) = max(ale_wake(i),Ale_bl(i))
           !cc nrlmd le 10/04/2012----------Stochastic triggering--------------
           if (iflag_trig_bl.ge.1) then
              ALE(i) = max(ale_wake(i),Ale_bl_trig(i))
           endif
           !cc fin nrlmd le 10/04/2012
           if (alp_offset>=0.) then
              ALP(i) = alp_wake(i) + Alp_bl(i) + alp_offset ! modif sb
           else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                _                  _
! Ajout d'une composante 3 * A * w w'2 a  w'3  avec w=www : w max sous pbase
!       ou A est la fraction couverte par les ascendances w'
!       on utilise le fait que A * w'3 = ALP
!       et donc A * w'2 ~ ALP / sqrt(ALE)  (on ajoute 0.1 pour les
!       singularites)
             ALP(i)=alp_wake(i)*(1.+3.*www(i)/( sqrt(ale_wake(i))+0.1) ) &
  &                +alp_bl(i)  *(1.+3.*www(i)/( sqrt(ale_bl(i))  +0.1) )
!             ALP(i)=alp_wake(i)+Alp_bl(i)+alp_offset*min(omega(i,6),0.)
!             if (alp(i)<0.) then
!                print*,'ALP ',alp(i),alp_wake(i) &
!                     ,Alp_bl(i),alp_offset*min(omega(i,6),0.)
!             endif
           endif
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     endif
     do i=1,klon
        if (alp(i)>alp_max) then
           IF(prt_level>9)WRITE(lunout,*)                             &
                'WARNING SUPER ALP (seuil=',alp_max, &
                '): i, alp, alp_wake,ale',i,alp(i),alp_wake(i),ale(i)
           alp(i)=alp_max
        endif
        if (ale(i)>ale_max) then
           IF(prt_level>9)WRITE(lunout,*)                             &
                'WARNING SUPER ALE (seuil=',ale_max, &
                '): i, alp, alp_wake,ale',i,ale(i),ale_wake(i),alp(i)
           ale(i)=ale_max
        endif
     enddo

     !fin calcul ale et alp
     !=======================================================================


     ! sb, oct02:
     ! Schema de convection modularise et vectorise:
     ! (driver commun aux versions 3 et 4)
     !
     IF (ok_cvl) THEN ! new driver for convectL

        IF (type_trac == 'repr') THEN
           nbtr_tmp=ntra
        ELSE
           nbtr_tmp=nbtr
        END IF
        !jyg   iflag_con est dans clesphys
        !c          CALL concvl (iflag_con,iflag_clos,
        CALL concvl (iflag_clos, &
             dtime,paprs,pplay,t_undi,q_undi, &
             t_wake,q_wake,wake_s, &
             u_seri,v_seri,tr_seri,nbtr_tmp, &
             ALE,ALP, &
             sig1,w01, &
             d_t_con,d_q_con,d_u_con,d_v_con,d_tr, &
             rain_con, snow_con, ibas_con, itop_con, sigd, &
             ema_cbmf,plcl,plfc,wbeff,upwd,dnwd,dnwd0, &
             Ma,mip,Vprecip,cape,cin,tvp,Tconv,iflagctrl, &
             pbase,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr,qcondc,wd, &
             ! RomP >>>
             !!     .        pmflxr,pmflxs,da,phi,mp,
             !!     .        ftd,fqd,lalim_conv,wght_th)
             pmflxr,pmflxs,da,phi,mp,phi2,d1a,dam,sij,clw,elij, &
             ftd,fqd,lalim_conv,wght_th, &
             ev, ep,epmlmMm,eplaMm, &
             wdtrainA,wdtrainM,wght_cvfd)
        ! RomP <<<

        !IM begin
        !       print*,'physiq: cin pbase dnwd0 ftd fqd ',cin(1),pbase(1),
        !    .dnwd0(1,1),ftd(1,1),fqd(1,1)
        !IM end
        !IM cf. FH
        clwcon0=qcondc
        pmfu(:,:)=upwd(:,:)+dnwd(:,:)

        do i = 1, klon
           if (iflagctrl(i).le.1) itau_con(i)=itau_con(i)+1
        enddo

     ELSE ! ok_cvl

        ! MAF conema3 ne contient pas les traceurs
        CALL conema3 (dtime, &
             paprs,pplay,t_seri,q_seri, &
             u_seri,v_seri,tr_seri,ntra, &
             sig1,w01, &
             d_t_con,d_q_con,d_u_con,d_v_con,d_tr, &
             rain_con, snow_con, ibas_con, itop_con, &
             upwd,dnwd,dnwd0,bas,top, &
             Ma,cape,tvp,rflag, &
             pbase &
             ,bbase,dtvpdt1,dtvpdq1,dplcldt,dplcldr &
             ,clwcon0)

     ENDIF ! ok_cvl

     !
     ! Correction precip
     rain_con = rain_con * cvl_corr
     snow_con = snow_con * cvl_corr
     !

     IF (.NOT. ok_gust) THEN
        do i = 1, klon
           wd(i)=0.0
        enddo
     ENDIF

     ! =================================================================== c
     ! Calcul des proprietes des nuages convectifs
     !

     !   calcul des proprietes des nuages convectifs
     clwcon0(:,:)=fact_cldcon*clwcon0(:,:)
     call clouds_gno &
          (klon,klev,q_seri,zqsat,clwcon0,ptconv,ratqsc,rnebcon0)

     ! =================================================================== c

     DO i = 1, klon
        itop_con(i) = min(max(itop_con(i),1),klev)
        ibas_con(i) = min(max(ibas_con(i),1),itop_con(i))
     ENDDO

     DO i = 1, klon
        ema_pcb(i)  = paprs(i,ibas_con(i))
     ENDDO
     DO i = 1, klon
        ! L'idicage de itop_con peut cacher un pb potentiel
        ! FH sous la dictee de JYG, CR
        ema_pct(i)  = paprs(i,itop_con(i)+1)

        if (itop_con(i).gt.klev-3) then
           if(prt_level >= 9) then
              write(lunout,*)'La convection monte trop haut '
              write(lunout,*)'itop_con(,',i,',)=',itop_con(i)
           endif
        endif
     ENDDO
  ELSE IF (iflag_con.eq.0) THEN
     write(lunout,*) 'On n appelle pas la convection'
     clwcon0=0.
     rnebcon0=0.
     d_t_con=0.
     d_q_con=0.
     d_u_con=0.
     d_v_con=0.
     rain_con=0.
     snow_con=0.
     bas=1
     top=1
  ELSE
     WRITE(lunout,*) "iflag_con non-prevu", iflag_con
     call abort_gcm("physiq", "", 1)
  ENDIF

  !     CALL homogene(paprs, q_seri, d_q_con, u_seri,v_seri,
  !    .              d_u_con, d_v_con)

  CALL add_phys_tend(d_u_con, d_v_con, d_t_con, d_q_con, dql0, dqi0, paprs, &
       'convection')
  !----------------------------------------------------------------------------

  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after convect'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, rain_con, snow_con, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF
  !
  IF (check) THEN
     za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
     WRITE(lunout,*)"aprescon=", za
     zx_t = 0.0
     za = 0.0
     DO i = 1, klon
        za = za + airephy(i)/REAL(klon)
        zx_t = zx_t + (rain_con(i)+ &
             snow_con(i))*airephy(i)/REAL(klon)
     ENDDO
     zx_t = zx_t/za*dtime
     WRITE(lunout,*)"Precip=", zx_t
  ENDIF
  IF (zx_ajustq) THEN
     DO i = 1, klon
        z_apres(i) = 0.0
     ENDDO
     DO k = 1, klev
        DO i = 1, klon
           z_apres(i) = z_apres(i) + (q_seri(i,k)+ql_seri(i,k)) &
                *(paprs(i,k)-paprs(i,k+1))/RG
        ENDDO
     ENDDO
     DO i = 1, klon
        z_factor(i) = (z_avant(i)-(rain_con(i)+snow_con(i))*dtime) &
             /z_apres(i)
     ENDDO
     DO k = 1, klev
        DO i = 1, klon
           IF (z_factor(i).GT.(1.0+1.0E-08) .OR. &
                z_factor(i).LT.(1.0-1.0E-08)) THEN
              q_seri(i,k) = q_seri(i,k) * z_factor(i)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  zx_ajustq=.FALSE.

  !
  !=============================================================================
  !RR:Evolution de la poche froide: on ne fait pas de separation wake/env 
  !pour la couche limite diffuse pour l instant
  !
  !
  !!! nrlmd le 22/03/2011---Si on met les poches hors des thermiques il faut rajouter cette
  !------------------------- tendance calcul�e hors des poches froides
  !
  if (iflag_wake>=1) then
     DO k=1,klev
        DO i=1,klon
           dt_dwn(i,k)  = ftd(i,k) 
           dq_dwn(i,k)  = fqd(i,k) 
           M_dwn(i,k)   = dnwd0(i,k)
           M_up(i,k)    = upwd(i,k)
           dt_a(i,k)    = d_t_con(i,k)/dtime - ftd(i,k) 
           dq_a(i,k)    = d_q_con(i,k)/dtime - fqd(i,k)
        ENDDO
     ENDDO
!nrlmd+jyg<
     DO k=1,klev
        DO i=1,klon
          wdt_PBL(i,k) =  0.
          wdq_PBL(i,k) =  0.
          udt_PBL(i,k) =  0.
          udq_PBL(i,k) =  0.
        ENDDO
     ENDDO
!
     IF (mod(iflag_pbl_split,2) .EQ. 1) THEN
       DO k=1,klev
        DO i=1,klon
       wdt_PBL(i,k) = wdt_PBL(i,k) + d_t_vdf_w(i,k)/dtime
       wdq_PBL(i,k) = wdq_PBL(i,k) + d_q_vdf_w(i,k)/dtime
       udt_PBL(i,k) = udt_PBL(i,k) + d_t_vdf_x(i,k)/dtime
       udq_PBL(i,k) = udq_PBL(i,k) + d_q_vdf_x(i,k)/dtime
!!        dt_dwn(i,k)  = dt_dwn(i,k) + d_t_vdf_w(i,k)/dtime
!!        dq_dwn(i,k)  = dq_dwn(i,k) + d_q_vdf_w(i,k)/dtime
!!        dt_a  (i,k)    = dt_a(i,k) + d_t_vdf_x(i,k)/dtime
!!        dq_a  (i,k)    = dq_a(i,k) + d_q_vdf_x(i,k)/dtime
        ENDDO
       ENDDO
      ENDIF
      IF (mod(iflag_pbl_split/2,2) .EQ. 1) THEN
       DO k=1,klev
        DO i=1,klon
!!        dt_dwn(i,k)  = dt_dwn(i,k) + 0.
!!        dq_dwn(i,k)  = dq_dwn(i,k) + 0.
!!        dt_a(i,k)   = dt_a(i,k)   + d_t_ajs(i,k)/dtime
!!        dq_a(i,k)   = dq_a(i,k)   + d_q_ajs(i,k)/dtime
        udt_PBL(i,k)   = udt_PBL(i,k)   + d_t_ajs(i,k)/dtime
        udq_PBL(i,k)   = udq_PBL(i,k)   + d_q_ajs(i,k)/dtime
        ENDDO
       ENDDO
      ENDIF
!>nrlmd+jyg

     IF (iflag_wake==2) THEN
        ok_wk_lsp(:)=max(sign(1.,wake_s(:)-wake_s_min_lsp),0.)
        DO k = 1,klev
           dt_dwn(:,k)= dt_dwn(:,k)+ &
                ok_wk_lsp(:)*(d_t_eva(:,k)+d_t_lsc(:,k))/dtime
           dq_dwn(:,k)= dq_dwn(:,k)+ &
                ok_wk_lsp(:)*(d_q_eva(:,k)+d_q_lsc(:,k))/dtime
        ENDDO
     ELSEIF (iflag_wake==3) THEN
        ok_wk_lsp(:)=max(sign(1.,wake_s(:)-wake_s_min_lsp),0.)
        DO k = 1,klev
           DO i=1,klon
              IF (rneb(i,k)==0.) THEN
! On ne tient compte des tendances qu'en dehors des nuages (c'est �|  dire
! a priri dans une region ou l'eau se reevapore).
                dt_dwn(i,k)= dt_dwn(i,k)+ &
                ok_wk_lsp(i)*d_t_lsc(i,k)/dtime
                dq_dwn(i,k)= dq_dwn(i,k)+ &
                ok_wk_lsp(i)*d_q_lsc(i,k)/dtime
              ENDIF
           ENDDO
        ENDDO
     ENDIF

     !
     !calcul caracteristiques de la poche froide
     call calWAKE (paprs,pplay,dtime &
          ,t_seri,q_seri,omega &
          ,dt_dwn,dq_dwn,M_dwn,M_up &
          ,dt_a,dq_a,sigd &
          ,wdt_PBL,wdq_PBL &
          ,udt_PBL,udq_PBL &
          ,wake_deltat,wake_deltaq,wake_dth &
          ,wake_h,wake_s,wake_dens &
          ,wake_pe,wake_fip,wake_gfl &
          ,dt_wake,dq_wake &
          ,wake_k, t_undi,q_undi &
          ,wake_omgbdth,wake_dp_omgb &
          ,wake_dtKE,wake_dqKE &
          ,wake_dtPBL,wake_dqPBL &
          ,wake_omg,wake_dp_deltomg &
          ,wake_spread,wake_Cstar,wake_d_deltat_gw &
          ,wake_ddeltat,wake_ddeltaq)
     !
     !-------------------------------------------------------------------------
     ! ajout des tendances des poches froides
     ! Faire rapidement disparaitre l'ancien dt_wake pour garder un d_t_wake
     ! coherent avec les autres d_t_...
     d_t_wake(:,:)=dt_wake(:,:)*dtime
     d_q_wake(:,:)=dq_wake(:,:)*dtime
     CALL add_phys_tend(du0,dv0,d_t_wake,d_q_wake,dql0,dqi0,paprs,'wake')
     !------------------------------------------------------------------------

  endif  ! (iflag_wake>=1)
  !
  !===================================================================
  !JYG
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after wake'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF

  !      print*,'apres callwake iflag_cldcon=', iflag_cldcon
  !
  !===================================================================
  ! Convection seche (thermiques ou ajustement)
  !===================================================================
  !
  call stratocu_if(klon,klev,pctsrf,paprs, pplay,t_seri &
       ,seuil_inversion,weak_inversion,dthmin)



  d_t_ajsb(:,:)=0.
  d_q_ajsb(:,:)=0.
  d_t_ajs(:,:)=0.
  d_u_ajs(:,:)=0.
  d_v_ajs(:,:)=0.
  d_q_ajs(:,:)=0.
  clwcon0th(:,:)=0.
  !
  !      fm_therm(:,:)=0.
  !      entr_therm(:,:)=0.
  !      detr_therm(:,:)=0.
  !
  IF(prt_level>9)WRITE(lunout,*) &
       'AVANT LA CONVECTION SECHE , iflag_thermals=' &
       ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals
  if(iflag_thermals<0) then
     !  Rien
     !  ====
     IF(prt_level>9)WRITE(lunout,*)'pas de convection seche'


  else

     !  Thermiques
     !  ==========
     IF(prt_level>9)WRITE(lunout,*)'JUSTE AVANT , iflag_thermals=' &
          ,iflag_thermals,'   nsplit_thermals=',nsplit_thermals


     !cc nrlmd le 10/04/2012
     DO k=1,klev+1
        DO i=1,klon
           pbl_tke_input(i,k,is_oce)=pbl_tke(i,k,is_oce)
           pbl_tke_input(i,k,is_ter)=pbl_tke(i,k,is_ter)
           pbl_tke_input(i,k,is_lic)=pbl_tke(i,k,is_lic)
           pbl_tke_input(i,k,is_sic)=pbl_tke(i,k,is_sic)
        ENDDO
     ENDDO
     !cc fin nrlmd le 10/04/2012

     if (iflag_thermals>=1) then
!jyg<
         IF (mod(iflag_pbl_split/2,2) .EQ. 1) THEN
!  Appel des thermiques avec les profils exterieurs aux poches
          DO k=1,klev
           DO i=1,klon
            t_therm(i,k) = t_seri(i,k) - wake_s(i)*wake_deltat(i,k)
            q_therm(i,k) = q_seri(i,k) - wake_s(i)*wake_deltaq(i,k)
           ENDDO
          ENDDO
         ELSE
!  Appel des thermiques avec les profils moyens
          DO k=1,klev
           DO i=1,klon
            t_therm(i,k) = t_seri(i,k)
            q_therm(i,k) = q_seri(i,k)
           ENDDO
          ENDDO
         ENDIF
!>jyg
        call calltherm(pdtphys &
             ,pplay,paprs,pphi,weak_inversion &
!!             ,u_seri,v_seri,t_seri,q_seri,zqsat,debut &  !jyg
             ,u_seri,v_seri,t_therm,q_therm,zqsat,debut &  !jyg
             ,d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs &
             ,fm_therm,entr_therm,detr_therm &
             ,zqasc,clwcon0th,lmax_th,ratqscth &
             ,ratqsdiff,zqsatth &
             !on rajoute ale et alp, et les caracteristiques de la couche alim
             ,Ale_bl,Alp_bl,lalim_conv,wght_th, zmax0, f0, zw2,fraca &
             ,ztv,zpspsk,ztla,zthl &
             !cc nrlmd le 10/04/2012
             ,pbl_tke_input,pctsrf,omega,airephy &
             ,zlcl_th,fraca0,w0,w_conv,therm_tke_max0,env_tke_max0 &
             ,n2,s2,ale_bl_stat &
             ,therm_tke_max,env_tke_max &
             ,alp_bl_det,alp_bl_fluct_m,alp_bl_fluct_tke &
             ,alp_bl_conv,alp_bl_stat &
             !cc fin nrlmd le 10/04/2012
             ,zqla,ztva )
!
!jyg<
         IF (mod(iflag_pbl_split/2,2) .EQ. 1) THEN
!  Si les thermiques ne sont presents que hors des poches, la tendance moyenne
!  associ�e doit etre multipliee par la fraction surfacique qu'ils couvrent.
          DO k=1,klev
           DO i=1,klon
!
            wake_deltat(i,k) = wake_deltat(i,k) - d_t_ajs(i,k)
            wake_deltaq(i,k) = wake_deltaq(i,k) - d_q_ajs(i,k)
            t_seri(i,k) = t_therm(i,k) + wake_s(i)*wake_deltat(i,k)
            q_seri(i,k) = q_therm(i,k) + wake_s(i)*wake_deltaq(i,k)
!
            d_u_ajs(i,k) = d_u_ajs(i,k)*(1.-wake_s(i)) 
            d_v_ajs(i,k) = d_v_ajs(i,k)*(1.-wake_s(i)) 
            d_t_ajs(i,k) = d_t_ajs(i,k)*(1.-wake_s(i)) 
            d_q_ajs(i,k) = d_q_ajs(i,k)*(1.-wake_s(i)) 
!
           ENDDO
          ENDDO
         ELSE
          DO k=1,klev
           DO i=1,klon
            t_seri(i,k) = t_therm(i,k)
            q_seri(i,k) = q_therm(i,k)
           ENDDO
          ENDDO
         ENDIF
!>jyg

        !cc nrlmd le 10/04/2012
        !-----------Stochastic triggering-----------
        if (iflag_trig_bl.ge.1) then
           !
           IF (prt_level .GE. 10) THEN
              print *,'cin, ale_bl_stat, alp_bl_stat ', &
                   cin, ale_bl_stat, alp_bl_stat
           ENDIF


           !----Initialisations
           do i=1,klon
              proba_notrig(i)=1.
              random_notrig(i)=1e6*ale_bl_stat(i)-int(1e6*ale_bl_stat(i))
              if ( ale_bl_trig(i) .lt. abs(cin(i))+1.e-10 ) then 
                 tau_trig(i)=tau_trig_shallow
              else
                 tau_trig(i)=tau_trig_deep
              endif
           enddo
           !
           IF (prt_level .GE. 10) THEN
              print *,'random_notrig, tau_trig ', &
                   random_notrig, tau_trig
              print *,'s_trig,s2,n2 ', &
                   s_trig,s2,n2
           ENDIF

           !Option pour re-activer l'ancien calcul de Ale_bl (iflag_trig_bl=2)
           IF (iflag_trig_bl.eq.1) then

              !----Tirage al\'eatoire et calcul de ale_bl_trig
              do i=1,klon
                 if ( (ale_bl_stat(i) .gt. abs(cin(i))+1.e-10) )  then
                    proba_notrig(i)=(1.-exp(-s_trig/s2(i)))** &
                         (n2(i)*dtime/tau_trig(i))
                    !        print *, 'proba_notrig(i) ',proba_notrig(i)
                    if (random_notrig(i) .ge. proba_notrig(i)) then 
                       ale_bl_trig(i)=ale_bl_stat(i)
                    else
                       ale_bl_trig(i)=0.
                    endif
                 else
                    proba_notrig(i)=1.
                    random_notrig(i)=0.
                    ale_bl_trig(i)=0.
                 endif
              enddo

           ELSE IF (iflag_trig_bl.eq.2) then

              do i=1,klon
                 if ( (Ale_bl(i) .gt. abs(cin(i))+1.e-10) )  then
                    proba_notrig(i)=(1.-exp(-s_trig/s2(i)))** &
                         (n2(i)*dtime/tau_trig(i))
                    !        print *, 'proba_notrig(i) ',proba_notrig(i)
                    if (random_notrig(i) .ge. proba_notrig(i)) then 
                       ale_bl_trig(i)=Ale_bl(i)
                    else
                       ale_bl_trig(i)=0.
                    endif
                 else
                    proba_notrig(i)=1.
                    random_notrig(i)=0.
                    ale_bl_trig(i)=0.
                 endif
              enddo

           ENDIF

           !
           IF (prt_level .GE. 10) THEN
              print *,'proba_notrig, ale_bl_trig ', &
                   proba_notrig, ale_bl_trig
           ENDIF

        endif !(iflag_trig_bl)

        !-----------Statistical closure-----------
        if (iflag_clos_bl.eq.1) then 

           do i=1,klon
              !CR: alp probabiliste
              if (ale_bl_trig(i).gt.0.) then
                 alp_bl(i)=alp_bl(i)/(1.-min(proba_notrig(i),0.999))
              endif
           enddo

        else if (iflag_clos_bl.eq.2) then

           !CR: alp calculee dans thermcell_main
           do i=1,klon
              alp_bl(i)=alp_bl_stat(i)
           enddo

        else

           alp_bl_stat(:)=0.

        endif !(iflag_clos_bl)

        IF (prt_level .GE. 10) THEN
           print *,'ale_bl_trig, alp_bl_stat ',ale_bl_trig, alp_bl_stat
        ENDIF

        !cc fin nrlmd le 10/04/2012

        ! ----------------------------------------------------------------------
        ! Transport de la TKE par les panaches thermiques.
        ! FH : 2010/02/01
        !     if (iflag_pbl.eq.10) then
        !     call thermcell_dtke(klon,klev,nbsrf,pdtphys,fm_therm,entr_therm,
        !    s           rg,paprs,pbl_tke)
        !     endif
        ! ----------------------------------------------------------------------
        !IM/FH: 2011/02/23 
        ! Couplage Thermiques/Emanuel seulement si T<0
        if (iflag_coupl==2) then
         IF (prt_level .GE. 10) THEN
           print*,'Couplage Thermiques/Emanuel seulement si T<0'
         ENDIF
           do i=1,klon
              if (t_seri(i,lmax_th(i))>273.) then
                 Ale_bl(i)=0.
              endif
           enddo
        endif

        do i=1,klon
           !           zmax_th(i)=pphi(i,lmax_th(i))/rg
           !CR:04/05/12:correction calcul zmax
           zmax_th(i)=zmax0(i) 
        enddo

     endif


     !  Ajustement sec
     !  ==============

     ! Dans le cas o\`u on active les thermiques, on fait partir l'ajustement
     ! a partir du sommet des thermiques.
     ! Dans le cas contraire, on demarre au niveau 1.

     if (iflag_thermals>=13.or.iflag_thermals<=0) then

        if(iflag_thermals.eq.0) then
           IF(prt_level>9)WRITE(lunout,*)'ajsec'
           limbas(:)=1
        else
           limbas(:)=lmax_th(:)
        endif

        ! Attention : le call ajsec_convV2 n'est maintenu que momentanneement
        ! pour des test de convergence numerique.
        ! Le nouveau ajsec est a priori mieux, meme pour le cas 
        ! iflag_thermals = 0 (l'ancienne version peut faire des tendances
        ! non nulles numeriquement pour des mailles non concernees.

        if (iflag_thermals==0) then
           ! Calling adjustment alone (but not the thermal plume model)
           CALL ajsec_convV2(paprs, pplay, t_seri,q_seri &
                , d_t_ajsb, d_q_ajsb)
        else if (iflag_thermals>0) then
           ! Calling adjustment above the top of thermal plumes
           CALL ajsec(paprs, pplay, t_seri,q_seri,limbas &
                , d_t_ajsb, d_q_ajsb)
        endif

        !-----------------------------------------------------------------------
        ! ajout des tendances de l'ajustement sec ou des thermiques
        CALL add_phys_tend(du0,dv0,d_t_ajsb,d_q_ajsb,dql0,dqi0,paprs,'ajsb')
        d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_ajsb(:,:)
        d_q_ajs(:,:)=d_q_ajs(:,:)+d_q_ajsb(:,:)

        !---------------------------------------------------------------------

     endif

  endif
  !
  !===================================================================
  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after dry_adjust'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF


  !-------------------------------------------------------------------------
  ! Computation of ratqs, the width (normalized) of the subrid scale 
  ! water distribution
  CALL  calcratqs(klon,klev,prt_level,lunout,        &
       iflag_ratqs,iflag_con,iflag_cldcon,pdtphys,  &
       ratqsbas,ratqshaut,tau_ratqs,fact_cldcon,   &
       ptconv,ptconvth,clwcon0th, rnebcon0th,     &
       paprs,pplay,q_seri,zqsat,fm_therm, &
       ratqs,ratqsc)


  !
  ! Appeler le processus de condensation a grande echelle
  ! et le processus de precipitation
  !-------------------------------------------------------------------------
  IF (prt_level .GE.10) THEN
     print *,'itap, ->fisrtilp ',itap
  ENDIF
  !
  CALL fisrtilp(dtime,paprs,pplay, &
       t_seri, q_seri,ptconv,ratqs, &
       d_t_lsc, d_q_lsc, d_ql_lsc, d_qi_lsc, rneb, cldliq, &
       rain_lsc, snow_lsc, &
       pfrac_impa, pfrac_nucl, pfrac_1nucl, &
       frac_impa, frac_nucl, beta_prec_fisrt, &
       prfl, psfl, rhcl,  &
       zqasc, fraca,ztv,zpspsk,ztla,zthl,iflag_cldcon, &
       iflag_ice_thermo)
  !
  WHERE (rain_lsc < 0) rain_lsc = 0.
  WHERE (snow_lsc < 0) snow_lsc = 0.

  CALL add_phys_tend(du0,dv0,d_t_lsc,d_q_lsc,d_ql_lsc,d_qi_lsc,paprs,'lsc')
  !---------------------------------------------------------------------------
  DO k = 1, klev
     DO i = 1, klon
        cldfra(i,k) = rneb(i,k)
!CR: a quoi ca sert? Faut-il ajouter qs_seri?
        IF (.NOT.new_oliq) cldliq(i,k) = ql_seri(i,k)
     ENDDO
  ENDDO
  IF (check) THEN
     za = qcheck(klon,klev,paprs,q_seri,ql_seri,airephy)
     WRITE(lunout,*)"apresilp=", za
     zx_t = 0.0
     za = 0.0
     DO i = 1, klon
        za = za + airephy(i)/REAL(klon)
        zx_t = zx_t + (rain_lsc(i) &
             + snow_lsc(i))*airephy(i)/REAL(klon)
     ENDDO
     zx_t = zx_t/za*dtime
     WRITE(lunout,*)"Precip=", zx_t
  ENDIF
  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after fisrt'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, rain_lsc, snow_lsc, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF

  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  !
  !-------------------------------------------------------------------
  !  PRESCRIPTION DES NUAGES POUR LE RAYONNEMENT
  !-------------------------------------------------------------------

  ! 1. NUAGES CONVECTIFS
  !
  !IM cf FH
  !     IF (iflag_cldcon.eq.-1) THEN ! seulement pour Tiedtke
  IF (iflag_cldcon.le.-1) THEN ! seulement pour Tiedtke
     snow_tiedtke=0.
     !     print*,'avant calcul de la pseudo precip '
     !     print*,'iflag_cldcon',iflag_cldcon
     if (iflag_cldcon.eq.-1) then
        rain_tiedtke=rain_con
     else
        !       print*,'calcul de la pseudo precip '
        rain_tiedtke=0.
        !         print*,'calcul de la pseudo precip 0'
        do k=1,klev
           do i=1,klon
              if (d_q_con(i,k).lt.0.) then
                 rain_tiedtke(i)=rain_tiedtke(i)-d_q_con(i,k)/pdtphys &
                      *(paprs(i,k)-paprs(i,k+1))/rg
              endif
           enddo
        enddo
     endif
     !
     !     call dump2d(iim,jjm,rain_tiedtke(2:klon-1),'PSEUDO PRECIP ')
     !

     ! Nuages diagnostiques pour Tiedtke
     CALL diagcld1(paprs,pplay, &
          !IM cf FH  .             rain_con,snow_con,ibas_con,itop_con,
          rain_tiedtke,snow_tiedtke,ibas_con,itop_con, &
          diafra,dialiq)
     DO k = 1, klev
        DO i = 1, klon
           IF (diafra(i,k).GT.cldfra(i,k)) THEN
              cldliq(i,k) = dialiq(i,k)
              cldfra(i,k) = diafra(i,k)
           ENDIF
        ENDDO
     ENDDO

  ELSE IF (iflag_cldcon.ge.3) THEN
     !  On prend pour les nuages convectifs le max du calcul de la
     !  convection et du calcul du pas de temps precedent diminue d'un facteur
     !  facttemps
     facteur = pdtphys *facttemps
     do k=1,klev
        do i=1,klon
           rnebcon(i,k)=rnebcon(i,k)*facteur
           if (rnebcon0(i,k)*clwcon0(i,k).gt.rnebcon(i,k)*clwcon(i,k)) &
                then
              rnebcon(i,k)=rnebcon0(i,k)
              clwcon(i,k)=clwcon0(i,k)
           endif
        enddo
     enddo

     !
     !jq - introduce the aerosol direct and first indirect radiative forcings
     !jq - Johannes Quaas, 27/11/2003 (quaas@lmd.jussieu.fr)
     IF (flag_aerosol .gt. 0) THEN
        IF (iflag_rrtm .EQ. 0) THEN !--old radiation
           IF (.NOT. aerosol_couple) THEN
              !
              CALL readaerosol_optic( &
                   debut, new_aod, flag_aerosol, itap, jD_cur-jD_ref, &
                   pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
                   mass_solu_aero, mass_solu_aero_pi,  &
                   tau_aero, piz_aero, cg_aero,  &
                   tausum_aero, tau3d_aero)
           ENDIF
         ELSE                       ! RRTM radiation
           IF (aerosol_couple .AND. config_inca == 'aero' ) THEN
            abort_message='config_inca=aero et rrtm=1 impossible'
            call abort_gcm(modname,abort_message,1)
           ELSE
!
#ifdef CPP_RRTM
             CALL readaerosol_optic_rrtm( debut, aerosol_couple, &
             new_aod, flag_aerosol, itap, jD_cur-jD_ref, &
             pdtphys, pplay, paprs, t_seri, rhcl, presnivs,  &
             tr_seri, mass_solu_aero, mass_solu_aero_pi,  &
             tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm,  &
             tausum_aero, tau3d_aero)
#else

              abort_message='You should compile with -rrtm if running with iflag_rrtm=1'
              call abort_gcm(modname,abort_message,1)
#endif
              !
           ENDIF
        ENDIF
     ELSE
        tausum_aero(:,:,:) = 0.
        IF (iflag_rrtm .EQ. 0) THEN !--old radiation
           tau_aero(:,:,:,:) = 0.
           piz_aero(:,:,:,:) = 0.
           cg_aero(:,:,:,:)  = 0.
        ELSE
           tau_aero_sw_rrtm(:,:,:,:)=0.0
           piz_aero_sw_rrtm(:,:,:,:)=0.0
           cg_aero_sw_rrtm(:,:,:,:)=0.0
        ENDIF
     ENDIF
     !
     !--STRAT AEROSOL
     !--updates tausum_aero,tau_aero,piz_aero,cg_aero
     IF (flag_aerosol_strat) THEN
        IF (prt_level .GE.10) THEN
         PRINT *,'appel a readaerosolstrat', mth_cur
        ENDIF
        IF (iflag_rrtm.EQ.0) THEN
           CALL readaerosolstrato(debut)
        ELSE
#ifdef CPP_RRTM
           CALL readaerosolstrato_rrtm(debut)
#else

           abort_message='You should compile with -rrtm if running with iflag_rrtm=1'
           call abort_gcm(modname,abort_message,1)
#endif
        ENDIF
     ENDIF
     !--fin STRAT AEROSOL


     !   On prend la somme des fractions nuageuses et des contenus en eau

     if (iflag_cldcon>=5) then

        do k=1,klev
           ptconvth(:,k)=fm_therm(:,k+1)>0.
        enddo

        if (iflag_coupl==4) then

           ! Dans le cas iflag_coupl==4, on prend la somme des convertures
           ! convectives et lsc dans la partie des thermiques
           ! Le controle par iflag_coupl est peut etre provisoire.
           do k=1,klev
              do i=1,klon
                 if (ptconv(i,k).and.ptconvth(i,k)) then
                    cldliq(i,k)=cldliq(i,k)+rnebcon(i,k)*clwcon(i,k)
                    cldfra(i,k)=min(cldfra(i,k)+rnebcon(i,k),1.)
                 else if (ptconv(i,k)) then
                    cldfra(i,k)=rnebcon(i,k)
                    cldliq(i,k)=rnebcon(i,k)*clwcon(i,k)
                 endif
              enddo
           enddo

        else if (iflag_coupl==5) then
           do k=1,klev
              do i=1,klon
                 cldfra(i,k)=min(cldfra(i,k)+rnebcon(i,k),1.)
                 cldliq(i,k)=cldliq(i,k)+rnebcon(i,k)*clwcon(i,k)
              enddo
           enddo

        else

           ! Si on est sur un point touche par la convection profonde et pas
           ! par les thermiques, on prend la couverture nuageuse et l'eau nuageuse
           ! de la convection profonde.

           !IM/FH: 2011/02/23 
           ! definition des points sur lesquels ls thermiques sont actifs

           do k=1,klev
              do i=1,klon
                 if (ptconv(i,k).and. .not. ptconvth(i,k)) then
                    cldfra(i,k)=rnebcon(i,k)
                    cldliq(i,k)=rnebcon(i,k)*clwcon(i,k)
                 endif
              enddo
           enddo

        endif

     else

        ! Ancienne version
        cldfra(:,:)=min(max(cldfra(:,:),rnebcon(:,:)),1.)
        cldliq(:,:)=cldliq(:,:)+rnebcon(:,:)*clwcon(:,:)
     endif

  ENDIF

  !     plulsc(:)=0.
  !     do k=1,klev,-1
  !        do i=1,klon
  !              zzz=prfl(:,k)+psfl(:,k)
  !           if (.not.ptconvth.zzz.gt.0.)
  !        enddo prfl, psfl,
  !     enddo
  !
  ! 2. NUAGES STARTIFORMES
  !
  IF (ok_stratus) THEN
     CALL diagcld2(paprs,pplay,t_seri,q_seri, diafra,dialiq)
     DO k = 1, klev
        DO i = 1, klon
           IF (diafra(i,k).GT.cldfra(i,k)) THEN
              cldliq(i,k) = dialiq(i,k)
              cldfra(i,k) = diafra(i,k)
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  ! Precipitation totale
  !
  DO i = 1, klon
     rain_fall(i) = rain_con(i) + rain_lsc(i)
     snow_fall(i) = snow_con(i) + snow_lsc(i)
  ENDDO
  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit="after diagcld"
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF
  !
  ! Calculer l'humidite relative pour diagnostique
  !
  DO k = 1, klev
     DO i = 1, klon
        zx_t = t_seri(i,k)
        IF (thermcep) THEN
           if (iflag_ice_thermo.eq.0) then 
           zdelta = MAX(0.,SIGN(1.,rtt-zx_t))
           else
           zdelta = MAX(0.,SIGN(1.,t_glace_min-zx_t))
           endif
           zx_qs  = r2es * FOEEW(zx_t,zdelta)/pplay(i,k)
           zx_qs  = MIN(0.5,zx_qs)
           zcor   = 1./(1.-retv*zx_qs)
           zx_qs  = zx_qs*zcor
        ELSE
           IF (zx_t.LT.t_coup) THEN
              zx_qs = qsats(zx_t)/pplay(i,k)
           ELSE
              zx_qs = qsatl(zx_t)/pplay(i,k)
           ENDIF
        ENDIF
        zx_rh(i,k) = q_seri(i,k)/zx_qs
        zqsat(i,k)=zx_qs
     ENDDO
  ENDDO

  !IM Calcul temp.potentielle a 2m (tpot) et temp. potentielle 
  !   equivalente a 2m (tpote) pour diagnostique
  !
  DO i = 1, klon
     tpot(i)=zt2m(i)*(100000./paprs(i,1))**RKAPPA
     IF (thermcep) THEN
        IF(zt2m(i).LT.RTT) then
           Lheat=RLSTT
        ELSE
           Lheat=RLVTT
        ENDIF
     ELSE
        IF (zt2m(i).LT.RTT) THEN
           Lheat=RLSTT
        ELSE
           Lheat=RLVTT
        ENDIF
     ENDIF
     tpote(i) = tpot(i)*      &
          EXP((Lheat *qsat2m(i))/(RCPD*zt2m(i)))
  ENDDO

  IF (type_trac == 'inca') THEN
#ifdef INCA
     CALL VTe(VTphysiq)
     CALL VTb(VTinca)
     calday = REAL(days_elapsed + 1) + jH_cur

     call chemtime(itap+itau_phy-1, date0, dtime)
     IF (config_inca == 'aero' .OR. config_inca == 'aeNP') THEN
        CALL AEROSOL_METEO_CALC( &
             calday,pdtphys,pplay,paprs,t,pmflxr,pmflxs, &
             prfl,psfl,pctsrf,airephy,rlat,rlon,u10m,v10m)
     END IF

     zxsnow_dummy(:) = 0.0

     CALL chemhook_begin (calday, &
          days_elapsed+1, &
          jH_cur, &
          pctsrf(1,1), &
          rlat, &
          rlon, &
          airephy, &
          paprs, &
          pplay, &
          coefh(1:klon,1:klev,is_ave), &
          pphi, &
          t_seri, &
          u, &
          v, &
          wo(:, :, 1), &
          q_seri, &
          zxtsol, &
          zxsnow_dummy, &
          solsw, &
          albsol1, &
          rain_fall, &
          snow_fall, &
          itop_con, &
          ibas_con, &
          cldfra, &
          iim, &
          jjm, &
          tr_seri, &
          ftsol, &
          paprs, &
          cdragh, &
          cdragm, &
          pctsrf, &
          pdtphys, &
          itap)

     CALL VTe(VTinca)
     CALL VTb(VTphysiq)
#endif 
  END IF !type_trac = inca
  !     
  ! Calculer les parametres optiques des nuages et quelques
  ! parametres pour diagnostiques:
  !

  IF (aerosol_couple) THEN 
     mass_solu_aero(:,:)    = ccm(:,:,1) 
     mass_solu_aero_pi(:,:) = ccm(:,:,2) 
  END IF

  if (ok_newmicro) then
     IF (iflag_rrtm.NE.0) THEN 
#ifdef CPP_RRTM
        IF (ok_cdnc.AND.NRADLP.NE.3) THEN
           abort_message='RRTM choix incoherent NRADLP doit etre egal a 3 pour ok_cdnc' 
           call abort_gcm(modname,abort_message,1)
        endif
#else

        abort_message='You should compile with -rrtm if running with iflag_rrtm=1'
        call abort_gcm(modname,abort_message,1)
#endif
     ENDIF
     CALL newmicro (ok_cdnc, bl95_b0, bl95_b1, &
          paprs, pplay, t_seri, cldliq, cldfra, &
          cldtau, cldemi, cldh, cldl, cldm, cldt, cldq, &
          flwp, fiwp, flwc, fiwc, &
          mass_solu_aero, mass_solu_aero_pi, &
          cldtaupi, re, fl, ref_liq, ref_ice, &
          ref_liq_pi, ref_ice_pi)
  else
     CALL nuage (paprs, pplay, &
          t_seri, cldliq, cldfra, cldtau, cldemi, &
          cldh, cldl, cldm, cldt, cldq, &
          ok_aie, &
          mass_solu_aero, mass_solu_aero_pi, &
          bl95_b0, bl95_b1, &
          cldtaupi, re, fl)
  endif
  !
  !IM betaCRF
  !
  cldtaurad   = cldtau
  cldtaupirad = cldtaupi
  cldemirad   = cldemi

  !
  if(lon1_beta.EQ.-180..AND.lon2_beta.EQ.180..AND. &
       lat1_beta.EQ.90..AND.lat2_beta.EQ.-90.) THEN
     !
     ! global
     !
     DO k=1, klev
        DO i=1, klon
           if (pplay(i,k).GE.pfree) THEN
              beta(i,k) = beta_pbl
           else
              beta(i,k) = beta_free
           endif
           if (mskocean_beta) THEN
              beta(i,k) = beta(i,k) * pctsrf(i,is_oce)
           endif
           cldtaurad(i,k)   = cldtau(i,k) * beta(i,k)
           cldtaupirad(i,k) = cldtaupi(i,k) * beta(i,k)
           cldemirad(i,k)   = cldemi(i,k) * beta(i,k)
           cldfrarad(i,k)   = cldfra(i,k) * beta(i,k)
        ENDDO
     ENDDO
     !
  else
     !
     ! regional
     !
     DO k=1, klev
        DO i=1,klon
           !
           if (rlon(i).ge.lon1_beta.AND.rlon(i).le.lon2_beta.AND. &
                rlat(i).le.lat1_beta.AND.rlat(i).ge.lat2_beta) THEN
              if (pplay(i,k).GE.pfree) THEN
                 beta(i,k) = beta_pbl
              else
                 beta(i,k) = beta_free
              endif
              if (mskocean_beta) THEN
                 beta(i,k) = beta(i,k) * pctsrf(i,is_oce)
              endif
              cldtaurad(i,k)   = cldtau(i,k) * beta(i,k)
              cldtaupirad(i,k) = cldtaupi(i,k) * beta(i,k)
              cldemirad(i,k)   = cldemi(i,k) * beta(i,k)
              cldfrarad(i,k)   = cldfra(i,k) * beta(i,k)
           endif
           !
        ENDDO
     ENDDO
     !
  endif
  !
  ! Appeler le rayonnement mais calculer tout d'abord l'albedo du sol.
  !
  IF (MOD(itaprad,radpas).EQ.0) THEN

     DO i = 1, klon
        albsol1(i) = falb1(i,is_oce) * pctsrf(i,is_oce) &
             + falb1(i,is_lic) * pctsrf(i,is_lic) &
             + falb1(i,is_ter) * pctsrf(i,is_ter) &
             + falb1(i,is_sic) * pctsrf(i,is_sic)
        albsol2(i) = falb2(i,is_oce) * pctsrf(i,is_oce) &
             + falb2(i,is_lic) * pctsrf(i,is_lic) &
             + falb2(i,is_ter) * pctsrf(i,is_ter) &
             + falb2(i,is_sic) * pctsrf(i,is_sic)
     ENDDO

     if (mydebug) then
        call writefield_phy('u_seri',u_seri,llm)
        call writefield_phy('v_seri',v_seri,llm)
        call writefield_phy('t_seri',t_seri,llm)
        call writefield_phy('q_seri',q_seri,llm)
     endif

     IF (aerosol_couple) THEN 
#ifdef INCA
        CALL radlwsw_inca  &
             (kdlon,kflev,dist, rmu0, fract, solaire, &
             paprs, pplay,zxtsol,albsol1, albsol2, t_seri,q_seri, &
             wo(:, :, 1), &
             cldfrarad, cldemirad, cldtaurad, &
             heat,heat0,cool,cool0,radsol,albpla, &
             topsw,toplw,solsw,sollw, &
             sollwdown, &
             topsw0,toplw0,solsw0,sollw0, &
             lwdn0, lwdn, lwup0, lwup,  &
             swdn0, swdn, swup0, swup, &
             ok_ade, ok_aie, &
             tau_aero, piz_aero, cg_aero, &
             topswad_aero, solswad_aero, &
             topswad0_aero, solswad0_aero, &
             topsw_aero, topsw0_aero, &
             solsw_aero, solsw0_aero, &
             cldtaupirad, &
             topswai_aero, solswai_aero)

#endif
     ELSE
        !
        !IM calcul radiatif pour le cas actuel
        !
        RCO2 = RCO2_act
        RCH4 = RCH4_act
        RN2O = RN2O_act
        RCFC11 = RCFC11_act
        RCFC12 = RCFC12_act
        !
        IF (prt_level .GE.10) THEN
           print *,' ->radlwsw, number 1 '
        ENDIF
        !
        CALL radlwsw &
             (dist, rmu0, fract,  &
             paprs, pplay,zxtsol,albsol1, albsol2,  &
             t_seri,q_seri,wo, &
             cldfrarad, cldemirad, cldtaurad, &
             ok_ade.OR.flag_aerosol_strat, ok_aie, flag_aerosol, &
             flag_aerosol_strat, &
             tau_aero, piz_aero, cg_aero, &
             tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm,&     ! Rajoute par OB pour RRTM
             tau_aero_lw_rrtm, & 
             cldtaupirad,new_aod, &
             zqsat, flwc, fiwc, &
             ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
             heat,heat0,cool,cool0,radsol,albpla, &
             topsw,toplw,solsw,sollw, &
             sollwdown, &
             topsw0,toplw0,solsw0,sollw0, &
             lwdn0, lwdn, lwup0, lwup,  &
             swdn0, swdn, swup0, swup, &
             topswad_aero, solswad_aero, &
             topswai_aero, solswai_aero, &
             topswad0_aero, solswad0_aero, &
             topsw_aero, topsw0_aero, &
             solsw_aero, solsw0_aero, &
             topswcf_aero, solswcf_aero, &
             !-C. Kleinschmitt for LW diagnostics
             toplwad_aero, sollwad_aero,&
             toplwai_aero, sollwai_aero, &
             toplwad0_aero, sollwad0_aero,&
             !-end
             ZLWFT0_i, ZFLDN0, ZFLUP0, &
             ZSWFT0_i, ZFSDN0, ZFSUP0)

        !
        !IM 2eme calcul radiatif pour le cas perturbe ou au moins un
        !IM des taux doit etre different du taux actuel
        !IM Par defaut on a les taux perturbes egaux aux taux actuels
        !
        if (ok_4xCO2atm) then
           if (RCO2_per.NE.RCO2_act.OR.RCH4_per.NE.RCH4_act.OR. &
                RN2O_per.NE.RN2O_act.OR.RCFC11_per.NE.RCFC11_act.OR. &
                RCFC12_per.NE.RCFC12_act) THEN
              !
              RCO2 = RCO2_per
              RCH4 = RCH4_per
              RN2O = RN2O_per
              RCFC11 = RCFC11_per
              RCFC12 = RCFC12_per
              !
              IF (prt_level .GE.10) THEN
                 print *,' ->radlwsw, number 2 '
              ENDIF
              !
              CALL radlwsw &
                   (dist, rmu0, fract,  &
                   paprs, pplay,zxtsol,albsol1, albsol2,  &
                   t_seri,q_seri,wo, &
                   cldfra, cldemi, cldtau, &
                   ok_ade.OR.flag_aerosol_strat, ok_aie, flag_aerosol, &
                   flag_aerosol_strat, &
                   tau_aero, piz_aero, cg_aero, &
                   tau_aero_sw_rrtm, piz_aero_sw_rrtm, cg_aero_sw_rrtm,&     ! Rajoute par OB pour RRTM
                   tau_aero_lw_rrtm, &
                   cldtaupi,new_aod, &
                   zqsat, flwc, fiwc, &
                   ref_liq, ref_ice, ref_liq_pi, ref_ice_pi, &
                   heatp,heat0p,coolp,cool0p,radsolp,albplap, &
                   topswp,toplwp,solswp,sollwp, &
                   sollwdownp, &
                   topsw0p,toplw0p,solsw0p,sollw0p, &
                   lwdn0p, lwdnp, lwup0p, lwupp,  &
                   swdn0p, swdnp, swup0p, swupp, &
                   topswad_aerop, solswad_aerop, &
                   topswai_aerop, solswai_aerop, &
                   topswad0_aerop, solswad0_aerop, &
                   topsw_aerop, topsw0_aerop, &
                   solsw_aerop, solsw0_aerop, &
                   topswcf_aerop, solswcf_aerop, &
                   !-C. Kleinschmitt for LW diagnostics
                   toplwad_aerop, sollwad_aerop,&
                   toplwai_aerop, sollwai_aerop, &
                   toplwad0_aerop, sollwad0_aerop,&
                   !-end
                   ZLWFT0_i, ZFLDN0, ZFLUP0, &
                   ZSWFT0_i, ZFSDN0, ZFSUP0)
           endif
        endif
        !
     ENDIF ! aerosol_couple
     itaprad = 0
  ENDIF ! MOD(itaprad,radpas)
  itaprad = itaprad + 1

  IF (iflag_radia.eq.0) THEN
     IF (prt_level.ge.9) THEN
        PRINT *,'--------------------------------------------------'
        PRINT *,'>>>> ATTENTION rayonnement desactive pour ce cas'
        PRINT *,'>>>>           heat et cool mis a zero '
        PRINT *,'--------------------------------------------------'
     END IF
     heat=0.
     cool=0.
     sollw=0.   ! MPL 01032011
     solsw=0.
     radsol=0.
     swup=0.    ! MPL 27102011 pour les fichiers AMMA_profiles et AMMA_scalars
     swup0=0.
     swdn=0.
     swdn0=0.
     lwup=0.
     lwup0=0.
     lwdn=0.
     lwdn0=0.
  END IF

  !
  ! Ajouter la tendance des rayonnements (tous les pas)
  !
  DO k = 1, klev
     DO i = 1, klon
        t_seri(i,k) = t_seri(i,k) &
             + (heat(i,k)-cool(i,k)) * dtime/RDAY
     ENDDO
  ENDDO
  !
  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after rad'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , topsw, toplw, solsw, sollw, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF
  !
  !
  ! Calculer l'hydrologie de la surface
  !
  !      CALL hydrol(dtime,pctsrf,rain_fall, snow_fall, zxevap,
  !     .            agesno, ftsol,fqsurf,fsnow, ruis)
  !

  !
  ! Calculer le bilan du sol et la derive de temperature (couplage)
  !
  DO i = 1, klon
     !         bils(i) = radsol(i) - sens(i) - evap(i)*RLVTT
     ! a la demande de JLD
     bils(i) = radsol(i) - sens(i) + zxfluxlat(i)
  ENDDO
  !
  !moddeblott(jan95)
  ! Appeler le programme de parametrisation de l'orographie
  ! a l'echelle sous-maille:
  !
  IF (prt_level .GE.10) THEN
     print *,' call orography ? ', ok_orodr
  ENDIF
  !
  IF (ok_orodr) THEN
     !
     !  selection des points pour lesquels le shema est actif:
     igwd=0
     DO i=1,klon
        itest(i)=0
        !        IF ((zstd(i).gt.10.0)) THEN
        IF (((zpic(i)-zmea(i)).GT.100.).AND.(zstd(i).GT.10.0)) THEN
           itest(i)=1
           igwd=igwd+1
           idx(igwd)=i
        ENDIF
     ENDDO
     !        igwdim=MAX(1,igwd)
     !
     IF (ok_strato) THEN

        CALL drag_noro_strato(klon,klev,dtime,paprs,pplay, &
             zmea,zstd, zsig, zgam, zthe,zpic,zval, &
             igwd,idx,itest, &
             t_seri, u_seri, v_seri, &
             zulow, zvlow, zustrdr, zvstrdr, &
             d_t_oro, d_u_oro, d_v_oro)

     ELSE
        CALL drag_noro(klon,klev,dtime,paprs,pplay, &
             zmea,zstd, zsig, zgam, zthe,zpic,zval, &
             igwd,idx,itest, &
             t_seri, u_seri, v_seri, &
             zulow, zvlow, zustrdr, zvstrdr, &
             d_t_oro, d_u_oro, d_v_oro)
     ENDIF
     !
     !  ajout des tendances
     !-----------------------------------------------------------------------------------------
     ! ajout des tendances de la trainee de l'orographie
     CALL add_phys_tend(d_u_oro,d_v_oro,d_t_oro,dq0,dql0,dqi0,paprs,'oro')
     !-----------------------------------------------------------------------------------------
     !
  ENDIF ! fin de test sur ok_orodr
  !
  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  IF (ok_orolf) THEN
     !
     !  selection des points pour lesquels le shema est actif:
     igwd=0
     DO i=1,klon
        itest(i)=0
        IF ((zpic(i)-zmea(i)).GT.100.) THEN
           itest(i)=1
           igwd=igwd+1
           idx(igwd)=i
        ENDIF
     ENDDO
     !        igwdim=MAX(1,igwd)
     !
     IF (ok_strato) THEN

        CALL lift_noro_strato(klon,klev,dtime,paprs,pplay, &
             rlat,zmea,zstd,zpic,zgam,zthe,zpic,zval, &
             igwd,idx,itest, &
             t_seri, u_seri, v_seri, &
             zulow, zvlow, zustrli, zvstrli, &
             d_t_lif, d_u_lif, d_v_lif               )

     ELSE
        CALL lift_noro(klon,klev,dtime,paprs,pplay, &
             rlat,zmea,zstd,zpic, &
             itest, &
             t_seri, u_seri, v_seri, &
             zulow, zvlow, zustrli, zvstrli, &
             d_t_lif, d_u_lif, d_v_lif)
     ENDIF
     !   
     !-----------------------------------------------------------------------------------------
     ! ajout des tendances de la portance de l'orographie
     CALL add_phys_tend(d_u_lif,d_v_lif,d_t_lif,dq0,dql0,dqi0,paprs,'lif')
     !-----------------------------------------------------------------------------------------
     !
  ENDIF ! fin de test sur ok_orolf
  !  HINES GWD PARAMETRIZATION

  IF (ok_hines) then

     CALL hines_gwd(klon,klev,dtime,paprs,pplay, &
          rlat,t_seri,u_seri,v_seri, &
          zustrhi,zvstrhi, &
          d_t_hin, d_u_hin, d_v_hin)
     !
     !  ajout des tendances
     CALL add_phys_tend(d_u_hin,d_v_hin,d_t_hin,dq0,dql0,dqi0,paprs,'hin')

  ENDIF

  if (ok_gwd_rando) then
     call FLOTT_GWD_rando(DTIME, pplay, t_seri, u_seri, v_seri, &
          rain_fall + snow_fall, zustr_gwd_rando, zvstr_gwd_rando, &
          du_gwd_rando, dv_gwd_rando)
     CALL add_phys_tend(du_gwd_rando, dv_gwd_rando, dt0, dq0, dql0,dqi0,paprs, &
          'flott_gwd_rando')
  end if

  ! STRESS NECESSAIRES: TOUTE LA PHYSIQUE

  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  DO i = 1, klon
     zustrph(i)=0.
     zvstrph(i)=0.
  ENDDO
  DO k = 1, klev
     DO i = 1, klon
        zustrph(i)=zustrph(i)+(u_seri(i,k)-u(i,k))/dtime* &
             (paprs(i,k)-paprs(i,k+1))/rg
        zvstrph(i)=zvstrph(i)+(v_seri(i,k)-v(i,k))/dtime* &
             (paprs(i,k)-paprs(i,k+1))/rg
     ENDDO
  ENDDO
  !
  !IM calcul composantes axiales du moment angulaire et couple des montagnes
  !
  IF (is_sequential .and. ok_orodr) THEN 
     CALL aaam_bud (27,klon,klev,jD_cur-jD_ref,jH_cur, &
          ra,rg,romega, &
          rlat,rlon,pphis, &
          zustrdr,zustrli,zustrph, &
          zvstrdr,zvstrli,zvstrph, &
          paprs,u,v, &
          aam, torsfc)
  ENDIF
  !IM cf. FLott END
  !IM
  IF (ip_ebil_phy.ge.2) THEN 
     ztit='after orography'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,2,2,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     call diagphy(airephy,ztit,ip_ebil_phy &
          , zero_v, zero_v, zero_v, zero_v, zero_v &
          , zero_v, zero_v, zero_v, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
  END IF

  !DC Calcul de la tendance due au methane
  IF(ok_qch4) THEN
     CALL METHOX(1,klon,klon,klev,q_seri,d_q_ch4,pplay)
  ! ajout de la tendance d'humidite due au methane
     CALL add_phys_tend(du0,dv0,dt0,d_q_ch4*dtime,dql0,'q_ch4')
  END IF
  !
  !
  !====================================================================
  ! Interface Simulateur COSP (Calipso, ISCCP, MISR, ..)
  !====================================================================
  ! Abderrahmane 24.08.09

  IF (ok_cosp) THEN
     ! adeclarer 
#ifdef CPP_COSP
     IF (itap.eq.1.or.MOD(itap,NINT(freq_cosp/dtime)).EQ.0) THEN

      IF (prt_level .GE.10) THEN
        print*,'freq_cosp',freq_cosp
      ENDIF
        mr_ozone=wo(:, :, 1) * dobson_u * 1e3 / zmasse
        !       print*,'Dans physiq.F avant appel cosp ref_liq,ref_ice=',
        !     s        ref_liq,ref_ice
        call phys_cosp(itap,dtime,freq_cosp, &
             ok_mensuelCOSP,ok_journeCOSP,ok_hfCOSP, &
             ecrit_mth,ecrit_day,ecrit_hf, ok_all_xml, &
             klon,klev,rlon,rlat,presnivs,overlap, &
             fract,ref_liq,ref_ice, &
             pctsrf(:,is_ter)+pctsrf(:,is_lic), &
             zu10m,zv10m,pphis, &
             zphi,paprs(:,1:klev),pplay,zxtsol,t_seri, &
             qx(:,:,ivap),zx_rh,cldfra,rnebcon,flwc,fiwc, &
             prfl(:,1:klev),psfl(:,1:klev), &
             pmflxr(:,1:klev),pmflxs(:,1:klev), &
             mr_ozone,cldtau, cldemi)

        !     L          calipso2D,calipso3D,cfadlidar,parasolrefl,atb,betamol,
        !     L          cfaddbze,clcalipso2,dbze,cltlidarradar,
        !     M          clMISR,
        !     R          clisccp2,boxtauisccp,boxptopisccp,tclisccp,ctpisccp,
        !     I          tauisccp,albisccp,meantbisccp,meantbclrisccp)

     ENDIF

#endif
  ENDIF  !ok_cosp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !AA
  !AA Installation de l'interface online-offline pour traceurs
  !AA
  !====================================================================
  !   Calcul  des tendances traceurs
  !====================================================================
  !

  IF (type_trac=='repr') THEN
     sh_in(:,:) = q_seri(:,:)
  ELSE
     sh_in(:,:) = qx(:,:,ivap)
  END IF

  call phytrac ( &
       itap,     days_elapsed+1,    jH_cur,   debut, &
       lafin,    dtime,     u, v,     t, &
       paprs,    pplay,     pmfu,     pmfd, &
       pen_u,    pde_u,     pen_d,    pde_d, &
       cdragh,   coefh(1:klon,1:klev,is_ave),   fm_therm, entr_therm, &
       u1,       v1,        ftsol,    pctsrf, &
       zustar,   zu10m,     zv10m, &
       wstar(:,is_ave),    ale_bl,         ale_wake, &
       rlat,     rlon, &
       frac_impa,frac_nucl, beta_prec_fisrt,beta_prec, &
       presnivs, pphis,     pphi,     albsol1, &
       sh_in,    rhcl,      cldfra,   rneb, &
       diafra,   cldliq,    itop_con, ibas_con, &
       pmflxr,   pmflxs,    prfl,     psfl, &
       da,       phi,       mp,       upwd, &
       phi2,     d1a,       dam,      sij, wght_cvfd, &        !<<RomP+RL
       wdtrainA, wdtrainM,  sigd,     clw,elij, &   !<<RomP
       ev,       ep,        epmlmMm,  eplaMm, &     !<<RomP
       dnwd,     aerosol_couple,      flxmass_w, &
       tau_aero, piz_aero,  cg_aero,  ccm, &
       rfname, &
       d_tr_dyn, &                                 !<<RomP
       tr_seri)

  IF (offline) THEN

     IF (prt_level.ge.9) &
          print*,'Attention on met a 0 les thermiques pour phystoke'
     call phystokenc ( &
          nlon,klev,pdtphys,rlon,rlat, &
          t,pmfu, pmfd, pen_u, pde_u, pen_d, pde_d, &
          fm_therm,entr_therm, &
          cdragh,coefh(1:klon,1:klev,is_ave),u1,v1,ftsol,pctsrf, &
          frac_impa, frac_nucl, &
          pphis,airephy,dtime,itap, &
          qx(:,:,ivap),da,phi,mp,upwd,dnwd)


  ENDIF

  !
  ! Calculer le transport de l'eau et de l'energie (diagnostique)
  !
  CALL transp (paprs,zxtsol, &
       t_seri, q_seri, u_seri, v_seri, zphi, &
       ve, vq, ue, uq)
  !
  !IM global posePB BEG
  IF(1.EQ.0) THEN
     !
     CALL transp_lay (paprs,zxtsol, &
          t_seri, q_seri, u_seri, v_seri, zphi, &
          ve_lay, vq_lay, ue_lay, uq_lay)
     !
  ENDIF !(1.EQ.0) THEN
  !IM global posePB END
  ! Accumuler les variables a stocker dans les fichiers histoire:
  !

  !================================================================
  ! Conversion of kinetic and potential energy into heat, for
  ! parameterisation of subgrid-scale motions
  !================================================================

  d_t_ec(:,:)=0.
  forall (k=1: llm) exner(:, k) = (pplay(:, k)/paprs(:,1))**RKAPPA
  CALL ener_conserv(klon,klev,pdtphys,u,v,t,qx(:,:,ivap), &
       u_seri,v_seri,t_seri,q_seri,pbl_tke(:,:,is_ave)-tke0(:,:), &
       zmasse,exner,d_t_ec)
  t_seri(:,:)=t_seri(:,:)+d_t_ec(:,:)

  !IM
  IF (ip_ebil_phy.ge.1) THEN 
     ztit='after physic'
     CALL diagetpq(airephy,ztit,ip_ebil_phy,1,1,dtime &
          , t_seri,q_seri,ql_seri,qs_seri,u_seri,v_seri,paprs,pplay &
          , d_h_vcol, d_qt, d_qw, d_ql, d_qs, d_ec)
     !     Comme les tendances de la physique sont ajoute dans la dynamique,
     !     on devrait avoir que la variation d'entalpie par la dynamique
     !     est egale a la variation de la physique au pas de temps precedent.
     !     Donc la somme de ces 2 variations devrait etre nulle.

     call diagphy(airephy,ztit,ip_ebil_phy &
          , topsw, toplw, solsw, sollw, sens &
          , evap, rain_fall, snow_fall, ztsol &
          , d_h_vcol, d_qt, d_ec &
          , fs_bound, fq_bound )
     !
     d_h_vcol_phy=d_h_vcol
     !
  END IF
  !
  !=======================================================================
  !   SORTIES
  !=======================================================================
  !
  !IM initialisation + calculs divers diag AMIP2
  !
  include "calcul_divers.h"
  !
  !IM Interpolation sur les niveaux de pression du NMC
  !   -------------------------------------------------
  !
  include "calcul_STDlev.h"
  !
  ! slp sea level pressure
  slp(:) = paprs(:,1)*exp(pphis(:)/(RD*t_seri(:,1)))
  !
  !cc prw = eau precipitable
  DO i = 1, klon
     prw(i) = 0.
     DO k = 1, klev
        prw(i) = prw(i) + &
             q_seri(i,k)*(paprs(i,k)-paprs(i,k+1))/RG
     ENDDO
  ENDDO
  !
  IF (type_trac == 'inca') THEN
#ifdef INCA
     CALL VTe(VTphysiq)
     CALL VTb(VTinca)

     CALL chemhook_end ( &
          dtime, &
          pplay, &
          t_seri, &
          tr_seri, &
          nbtr, &
          paprs, &
          q_seri, &
          airephy, &
          pphi, &
          pphis, &
          zx_rh)

     CALL VTe(VTinca)
     CALL VTb(VTphysiq)
#endif
  END IF


  !
  ! Convertir les incrementations en tendances
  !
  IF (prt_level .GE.10) THEN
     print *,'Convertir les incrementations en tendances '
  ENDIF
  !
  if (mydebug) then
     call writefield_phy('u_seri',u_seri,llm)
     call writefield_phy('v_seri',v_seri,llm)
     call writefield_phy('t_seri',t_seri,llm)
     call writefield_phy('q_seri',q_seri,llm)
  endif

  DO k = 1, klev
     DO i = 1, klon
        d_u(i,k) = ( u_seri(i,k) - u(i,k) ) / dtime
        d_v(i,k) = ( v_seri(i,k) - v(i,k) ) / dtime
        d_t(i,k) = ( t_seri(i,k)-t(i,k) ) / dtime
        d_qx(i,k,ivap) = ( q_seri(i,k) - qx(i,k,ivap) ) / dtime
        d_qx(i,k,iliq) = ( ql_seri(i,k) - qx(i,k,iliq) ) / dtime
!CR: on ajoute le contenu en glace
        if (nqo.eq.3) then
        d_qx(i,k,isol) = ( qs_seri(i,k) - qx(i,k,isol) ) / dtime
        endif
     ENDDO
  ENDDO
  !
!CR: nb de traceurs eau: nqo
!  IF (nqtot.GE.3) THEN
   IF (nqtot.GE.(nqo+1)) THEN
!     DO iq = 3, nqtot
     DO iq = nqo+1, nqtot
        DO  k = 1, klev
           DO  i = 1, klon
!              d_qx(i,k,iq) = ( tr_seri(i,k,iq-2) - qx(i,k,iq) ) / dtime
               d_qx(i,k,iq) = ( tr_seri(i,k,iq-nqo) - qx(i,k,iq) ) / dtime
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  !IM rajout diagnostiques bilan KP pour analyse MJO par Jun-Ichi Yano
  !IM global posePB      include "write_bilKP_ins.h"
  !IM global posePB      include "write_bilKP_ave.h"
  !

  ! Sauvegarder les valeurs de t et q a la fin de la physique:
  !
  DO k = 1, klev
     DO i = 1, klon
        u_ancien(i,k) = u_seri(i,k)
        v_ancien(i,k) = v_seri(i,k)
        t_ancien(i,k) = t_seri(i,k)
        q_ancien(i,k) = q_seri(i,k)
     ENDDO
  ENDDO

!!! RomP >>>
!CR: nb de traceurs eau: nqo
!  IF (nqtot.GE.3) THEN
   IF (nqtot.GE.(nqo+1)) THEN
!     DO iq = 3, nqtot
     DO iq = nqo+1, nqtot
        DO k = 1, klev
           DO i = 1, klon
!              tr_ancien(i,k,iq-2) = tr_seri(i,k,iq-2)
              tr_ancien(i,k,iq-nqo) = tr_seri(i,k,iq-nqo)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
!!! RomP <<<
  !==========================================================================
  ! Sorties des tendances pour un point particulier
  ! a utiliser en 1D, avec igout=1 ou en 3D sur un point particulier
  ! pour le debug
  ! La valeur de igout est attribuee plus haut dans le programme
  !==========================================================================

  if (prt_level.ge.1) then
     write(lunout,*) 'FIN DE PHYSIQ !!!!!!!!!!!!!!!!!!!!'
     write(lunout,*) &
          'nlon,klev,nqtot,debut,lafin,jD_cur, jH_cur, pdtphys pct tlos'
     write(lunout,*) &
          nlon,klev,nqtot,debut,lafin, jD_cur, jH_cur ,pdtphys, &
          pctsrf(igout,is_ter), pctsrf(igout,is_lic),pctsrf(igout,is_oce), &
          pctsrf(igout,is_sic)
     write(lunout,*) 'd_t_dyn,d_t_con,d_t_lsc,d_t_ajsb,d_t_ajs,d_t_eva'
     do k=1,klev
        write(lunout,*) d_t_dyn(igout,k),d_t_con(igout,k), &
             d_t_lsc(igout,k),d_t_ajsb(igout,k),d_t_ajs(igout,k), &
             d_t_eva(igout,k)
     enddo
     write(lunout,*) 'cool,heat'
     do k=1,klev
        write(lunout,*) cool(igout,k),heat(igout,k)
     enddo

     write(lunout,*) 'd_t_oli,d_t_vdf,d_t_oro,d_t_lif,d_t_ec'
     do k=1,klev
        write(lunout,*) d_t_oli(igout,k),d_t_vdf(igout,k), &
             d_t_oro(igout,k),d_t_lif(igout,k),d_t_ec(igout,k)
     enddo

     write(lunout,*) 'd_ps ',d_ps(igout)
     write(lunout,*) 'd_u, d_v, d_t, d_qx1, d_qx2 '
     do k=1,klev
        write(lunout,*) d_u(igout,k),d_v(igout,k),d_t(igout,k), &
             d_qx(igout,k,1),d_qx(igout,k,2)
     enddo
  endif

  !==========================================================================

  !============================================================
  !   Calcul de la temperature potentielle
  !============================================================
  DO k = 1, klev
     DO i = 1, klon
        !JYG/IM theta en debut du pas de temps
        !JYG/IM       theta(i,k)=t(i,k)*(100000./pplay(i,k))**(RD/RCPD)
        !JYG/IM theta en fin de pas de temps de physique
        theta(i,k)=t_seri(i,k)*(100000./pplay(i,k))**(RD/RCPD)
        ! thetal: 2 lignes suivantes a decommenter si vous avez les fichiers     MPL 20130625
        ! fth_fonctions.F90 et parkind1.F90
        ! sinon thetal=theta
        !       thetal(i,k)=fth_thetal(pplay(i,k),t_seri(i,k),q_seri(i,k),
        !    :         ql_seri(i,k))
        thetal(i,k)=theta(i,k)
     ENDDO
  ENDDO
  !

  ! 22.03.04 BEG
  !=============================================================
  !   Ecriture des sorties
  !=============================================================
#ifdef CPP_IOIPSL

  ! Recupere des varibles calcule dans differents modules
  ! pour ecriture dans histxxx.nc 

  ! Get some variables from module fonte_neige_mod
  CALL fonte_neige_get_vars(pctsrf,  &
       zxfqcalving, zxfqfonte, zxffonte)




  !=============================================================
  ! Separation entre thermiques et non thermiques dans les sorties
  ! de fisrtilp
  !=============================================================

  if (iflag_thermals>=1) then
     d_t_lscth=0.
     d_t_lscst=0.
     d_q_lscth=0.
     d_q_lscst=0.
     do k=1,klev
        do i=1,klon
           if (ptconvth(i,k)) then
              d_t_lscth(i,k)=d_t_eva(i,k)+d_t_lsc(i,k)
              d_q_lscth(i,k)=d_q_eva(i,k)+d_q_lsc(i,k)
           else
              d_t_lscst(i,k)=d_t_eva(i,k)+d_t_lsc(i,k)
              d_q_lscst(i,k)=d_q_eva(i,k)+d_q_lsc(i,k)
           endif
        enddo
     enddo

     do i=1,klon
        plul_st(i)=prfl(i,lmax_th(i)+1)+psfl(i,lmax_th(i)+1)
        plul_th(i)=prfl(i,1)+psfl(i,1)
     enddo
  endif


  !On effectue les sorties:

  CALL phys_output_write(itap, pdtphys, paprs, pphis,               &
       pplay, lmax_th, aerosol_couple,                 &
       ok_ade, ok_aie, ivap, new_aod, ok_sync,         &
       ptconv, read_climoz, clevSTD,                   &
       ptconvth, d_t, qx, d_qx, zmasse,                &
       flag_aerosol, flag_aerosol_strat, ok_cdnc)




  include "write_histday_seri.h"

  include "write_paramLMDZ_phy.h"

#endif

  ! 22.03.04 END
  !
  !====================================================================
  ! Si c'est la fin, il faut conserver l'etat de redemarrage
  !====================================================================
  !

  !        -----------------------------------------------------------------
  !        WSTATS: Saving statistics
  !        -----------------------------------------------------------------
  !        ("stats" stores and accumulates 8 key variables in file "stats.nc"
  !        which can later be used to make the statistic files of the run:
  !        "stats")          only possible in 3D runs !


  IF (callstats) THEN

     call wstats(klon,o_psol%name,"Surface pressure","Pa" &
          ,2,paprs(:,1))
     call wstats(klon,o_tsol%name,"Surface temperature","K", &
          2,zxtsol)
     zx_tmp_fi2d(:) = rain_fall(:) + snow_fall(:)
     call wstats(klon,o_precip%name,"Precip Totale liq+sol", &
          "kg/(s*m2)",2,zx_tmp_fi2d)
     zx_tmp_fi2d(:) = rain_lsc(:) + snow_lsc(:)
     call wstats(klon,o_plul%name,"Large-scale Precip", &
          "kg/(s*m2)",2,zx_tmp_fi2d)
     zx_tmp_fi2d(:) = rain_con(:) + snow_con(:)
     call wstats(klon,o_pluc%name,"Convective Precip", &
          "kg/(s*m2)",2,zx_tmp_fi2d)
     call wstats(klon,o_sols%name,"Solar rad. at surf.", &
          "W/m2",2,solsw)
     call wstats(klon,o_soll%name,"IR rad. at surf.", &
          "W/m2",2,sollw)
     zx_tmp_fi2d(:) = topsw(:)-toplw(:)
     call wstats(klon,o_nettop%name,"Net dn radiatif flux at TOA", &
          "W/m2",2,zx_tmp_fi2d)



     call wstats(klon,o_temp%name,"Air temperature","K", &
          3,t_seri)
     call wstats(klon,o_vitu%name,"Zonal wind","m.s-1", &
          3,u_seri)
     call wstats(klon,o_vitv%name,"Meridional wind", &
          "m.s-1",3,v_seri)
     call wstats(klon,o_vitw%name,"Vertical wind", &
          "m.s-1",3,omega)
     call wstats(klon,o_ovap%name,"Specific humidity", "kg/kg", &
          3,q_seri)



     IF(lafin) THEN
        write (*,*) "Writing stats..."
        call mkstats(ierr)
     ENDIF

  ENDIF !if callstats

  IF (lafin) THEN
     itau_phy = itau_phy + itap
     CALL phyredem ("restartphy.nc")
     !         open(97,form="unformatted",file="finbin")
     !         write(97) u_seri,v_seri,t_seri,q_seri
     !         close(97)
     !$OMP MASTER
     if (read_climoz >= 1) then
        if (is_mpi_root) then
           call nf95_close(ncid_climoz)
        end if
        deallocate(press_climoz) ! pointer
     end if
     !$OMP END MASTER
  ENDIF

  !      first=.false.

  RETURN
END SUBROUTINE physiq
FUNCTION qcheck(klon,klev,paprs,q,ql,aire)
  IMPLICIT none
  !
  ! Calculer et imprimer l'eau totale. A utiliser pour verifier
  ! la conservation de l'eau
  !
  include "YOMCST.h"
  INTEGER klon,klev
  REAL paprs(klon,klev+1), q(klon,klev), ql(klon,klev)
  REAL aire(klon)
  REAL qtotal, zx, qcheck
  INTEGER i, k
  !
  zx = 0.0
  DO i = 1, klon
     zx = zx + aire(i)
  ENDDO
  qtotal = 0.0
  DO k = 1, klev
     DO i = 1, klon
        qtotal = qtotal + (q(i,k)+ql(i,k)) * aire(i) &
             *(paprs(i,k)-paprs(i,k+1))/RG
     ENDDO
  ENDDO
  !
  qcheck = qtotal/zx
  !
  RETURN
END FUNCTION qcheck
SUBROUTINE gr_fi_ecrit(nfield,nlon,iim,jjmp1,fi,ecrit)
  IMPLICIT none
  !
  ! Tranformer une variable de la grille physique a
  ! la grille d'ecriture
  !
  INTEGER nfield,nlon,iim,jjmp1, jjm
  REAL fi(nlon,nfield), ecrit(iim*jjmp1,nfield)
  !
  INTEGER i, n, ig
  !
  jjm = jjmp1 - 1
  DO n = 1, nfield
     DO i=1,iim
        ecrit(i,n) = fi(1,n)
        ecrit(i+jjm*iim,n) = fi(nlon,n)
     ENDDO
     DO ig = 1, nlon - 2
        ecrit(iim+ig,n) = fi(1+ig,n)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE gr_fi_ecrit
