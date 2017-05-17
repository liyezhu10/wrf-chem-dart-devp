!
! $Id$
!
c
      SUBROUTINE dynredem0_loc(fichnom,iday_end,phis)
#ifdef CPP_IOIPSL
      USE IOIPSL
#endif
      USE parallel_lmdz
      USE mod_hallo
      USE infotrac
      USE control_mod, only : sub_day  ! DART-LMDZ

      IMPLICIT NONE
c=======================================================================
c Ecriture du fichier de redemarrage sous format NetCDF (initialisation)
c=======================================================================
c   Declarations:
c   -------------
#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comvert.h"
#include "comgeom.h"
#include "temps.h"
#include "ener.h"
#include "logic.h"
#include "netcdf.inc"
#include "description.h"
#include "serre.h"
#include "iniprint.h"

c   Arguments:
c   ----------
      INTEGER iday_end
      REAL phis(ijb_u:ije_u)
      CHARACTER*(*) fichnom

c   Local:
c   ------
      INTEGER iq,l
      INTEGER length
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      INTEGER ierr
      character*20 modname
      character*80 abort_message

c   Variables locales pour NetCDF:
c
      INTEGER dims2(2), dims3(3), dims4(4)
      INTEGER idim_index
      INTEGER idim_rlonu, idim_rlonv, idim_rlatu, idim_rlatv
      INTEGER idim_s, idim_sig
      INTEGER idim_tim
      INTEGER nid,nvarid

      REAL zan0,zjulian,hours
      INTEGER yyears0,jjour0, mmois0
      character*30 unites
      REAL :: phis_glo(ip1jmp1)
      
      CALL Gather_field_u(phis,phis_glo,1)
      
      
c-----------------------------------------------------------------------
      if (mpi_rank==0) then
      
      modname='dynredem0_loc'

#ifdef CPP_IOIPSL
      call ymds2ju(annee_ref, 1, iday_end, 0.0, zjulian)
      call ju2ymds(zjulian, yyears0, mmois0, jjour0, hours)
#else
! set yyears0, mmois0, jjour0 to 0,1,1 (hours is not used)
      yyears0=0
      mmois0=1
      jjour0=1
#endif                

      DO l=1,length
       tab_cntrl(l) = 0.
      ENDDO
       tab_cntrl(1)  =  REAL(iim)
       tab_cntrl(2)  =  REAL(jjm)
       tab_cntrl(3)  =  REAL(llm)
       tab_cntrl(4)  =  REAL(day_ref)
       tab_cntrl(5)  =  REAL(annee_ref)
       tab_cntrl(6)  = rad
       tab_cntrl(7)  = omeg
       tab_cntrl(8)  = g
       tab_cntrl(9)  = cpp
       tab_cntrl(10) = kappa
       tab_cntrl(11) = daysec
       tab_cntrl(12) = dtvr
       tab_cntrl(13) = etot0
       tab_cntrl(14) = ptot0
       tab_cntrl(15) = ztot0
       tab_cntrl(16) = stot0
       tab_cntrl(17) = ang0
       tab_cntrl(18) = pa
       tab_cntrl(19) = preff
c
c    .....    parametres  pour le zoom      ......   

       tab_cntrl(20)  = clon
       tab_cntrl(21)  = clat
       tab_cntrl(22)  = grossismx
       tab_cntrl(23)  = grossismy
c
      IF ( fxyhypb )   THEN
       tab_cntrl(24) = 1.
       tab_cntrl(25) = dzoomx
       tab_cntrl(26) = dzoomy
       tab_cntrl(27) = 0.
       tab_cntrl(28) = taux
       tab_cntrl(29) = tauy
      ELSE
       tab_cntrl(24) = 0.
       tab_cntrl(25) = dzoomx
       tab_cntrl(26) = dzoomy
       tab_cntrl(27) = 0.
       tab_cntrl(28) = 0.
       tab_cntrl(29) = 0.
       IF( ysinus )  tab_cntrl(27) = 1.
      ENDIF

       tab_cntrl(30) =  REAL(iday_end)
       tab_cntrl(31) =  REAL(itau_dyn + itaufin)

       if(sub_day) tab_cntrl(31) = REAL( itaufin)     ! DART-LMDZ 

c start_time: start_time of simulation (not necessarily 0.)
       tab_cntrl(32) = start_time
c
c    .........................................................
c
c Creation du fichier:
c
      ierr = NF_CREATE(fichnom, NF_CLOBBER, nid)
      IF (ierr.NE.NF_NOERR) THEN
         write(lunout,*)"dynredem0: Pb d ouverture du fichier "
     &                  //trim(fichnom)
         write(lunout,*)' ierr = ', ierr
         CALL ABORT_GCM("DYNREDEM0", "", 1)
      ENDIF
c
c Preciser quelques attributs globaux:
c
      ierr = NF_PUT_ATT_TEXT (nid, NF_GLOBAL, "title", 27,
     .                       "Fichier demmarage dynamique")
c
c Definir les dimensions du fichiers:
c
      ierr = NF_DEF_DIM (nid, "index", length, idim_index)
      ierr = NF_DEF_DIM (nid, "rlonu", iip1, idim_rlonu)
      ierr = NF_DEF_DIM (nid, "rlatu", jjp1, idim_rlatu)
      ierr = NF_DEF_DIM (nid, "rlonv", iip1, idim_rlonv)
      ierr = NF_DEF_DIM (nid, "rlatv", jjm, idim_rlatv)
      ierr = NF_DEF_DIM (nid, "sigs", llm, idim_s)
      ierr = NF_DEF_DIM (nid, "sig", llmp1, idim_sig)
      ierr = NF_DEF_DIM (nid, "temps", NF_UNLIMITED, idim_tim)
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
c
c Definir et enregistrer certains champs invariants:
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"controle",NF_DOUBLE,1,idim_index,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"controle",NF_FLOAT,1,idim_index,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Parametres de controle")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonu",NF_DOUBLE,1,idim_rlonu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonu",NF_FLOAT,1,idim_rlonu,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonu)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatu",NF_DOUBLE,1,idim_rlatu,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatu",NF_FLOAT,1,idim_rlatu,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatu)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlonv",NF_DOUBLE,1,idim_rlonv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlonv",NF_FLOAT,1,idim_rlonv,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 23,
     .                       "Longitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlonv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlonv)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"rlatv",NF_DOUBLE,1,idim_rlatv,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"rlatv",NF_FLOAT,1,idim_rlatv,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Latitudes des points V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,rlatv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,rlatv)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_DOUBLE,1,idim_s,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsigs",NF_FLOAT,1,idim_s,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 28,
     .                       "Numero naturel des couches s")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsigs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsigs)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"nivsig",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"nivsig",NF_FLOAT,1,idim_sig,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 32,
     .                       "Numero naturel des couches sigma")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,nivsig)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,nivsig)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ap",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ap",NF_FLOAT,1,idim_sig,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26,
     .                       "Coefficient A pour hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,ap)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,ap)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"bp",NF_DOUBLE,1,idim_sig,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"bp",NF_FLOAT,1,idim_sig,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 26,
     .                       "Coefficient B pour hybride")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,bp)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,bp)
#endif
c
      ierr = NF_REDEF (nid)
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"presnivs",NF_DOUBLE,1,idim_s,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"presnivs",NF_FLOAT,1,idim_s,nvarid)
#endif
cIM 220306 END
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,presnivs)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,presnivs)
#endif
c
c Coefficients de passage cov. <-> contra. <--> naturel
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonu
      dims2(2) = idim_rlatu
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cu",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cu",NF_FLOAT,2,dims2,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour U")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cu)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cu)
#endif
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatv
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"cv",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"cv",NF_FLOAT,2,dims2,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 29,
     .                       "Coefficient de passage pour V")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,cv)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,cv)
#endif
c
c Aire de chaque maille:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"aire",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"aire",NF_FLOAT,2,dims2,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 22,
     .                       "Aires de chaque maille")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,aire)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,aire)
#endif
c
c Geopentiel au sol:
c
      ierr = NF_REDEF (nid)
      dims2(1) = idim_rlonv
      dims2(2) = idim_rlatu
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"phisinit",NF_DOUBLE,2,dims2,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"phisinit",NF_FLOAT,2,dims2,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Geopotentiel au sol")
      ierr = NF_ENDDEF(nid)
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,phis_glo)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,phis_glo)
#endif
c
c Definir les variables pour pouvoir les enregistrer plus tard:
c
      ierr = NF_REDEF (nid) ! entrer dans le mode de definition
c
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"temps",NF_DOUBLE,1,idim_tim,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"temps",NF_FLOAT,1,idim_tim,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 19,
     .                       "Temps de simulation")
      write(unites,200)yyears0,mmois0,jjour0
200   format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "units", 30,
     .                         unites)

c
      dims4(1) = idim_rlonu
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ucov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ucov",NF_FLOAT,4,dims4,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse U")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatv
      dims4(3) = idim_s
      dims4(4) = idim_tim
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"vcov",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"vcov",NF_FLOAT,4,dims4,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 9,
     .                       "Vitesse V")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"teta",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"teta",NF_FLOAT,4,dims4,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 11,
     .                       "Temperature")
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim

      DO iq=1,nqtot
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,tname(iq),NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,tname(iq),NF_FLOAT,4,dims4,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,ttext(iq))
      ENDDO
c
      dims4(1) = idim_rlonv
      dims4(2) = idim_rlatu
      dims4(3) = idim_s
      dims4(4) = idim_tim
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"masse",NF_DOUBLE,4,dims4,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"masse",NF_FLOAT,4,dims4,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 12,
     .                       "C est quoi ?")
c
      dims3(1) = idim_rlonv
      dims3(2) = idim_rlatu
      dims3(3) = idim_tim
cIM 220306 BEG
#ifdef NC_DOUBLE
      ierr = NF_DEF_VAR (nid,"ps",NF_DOUBLE,3,dims3,nvarid)
#else
      ierr = NF_DEF_VAR (nid,"ps",NF_FLOAT,3,dims3,nvarid)
#endif
cIM 220306 END
      ierr = NF_PUT_ATT_TEXT (nid, nvarid, "title", 15,
     .                       "Pression au sol")
c
      ierr = NF_ENDDEF(nid) ! sortir du mode de definition
      ierr = NF_CLOSE(nid) ! fermer le fichier

      write(lunout,*)'dynredem_loc: iim,jjm,llm,iday_end',
     &               iim,jjm,llm,iday_end
      write(lunout,*)'dynredem_loc: rad,omeg,g,cpp,kappa',
     &        rad,omeg,g,cpp,kappa

      endif  ! mpi_rank==0
      RETURN
      END
      SUBROUTINE dynredem1_loc(fichnom,time,
     .                     vcov,ucov,teta,q,masse,ps)
      USE parallel_lmdz
      USE mod_hallo
      USE infotrac
      USE control_mod
      USE dynredem_mod
      IMPLICIT NONE
c=================================================================
c  Ecriture du fichier de redemarrage sous format NetCDF
c=================================================================
#include "dimensions.h"
#include "paramet.h"
#include "description.h"
#include "netcdf.inc"
#include "comvert.h"
#include "comgeom.h"
#include "temps.h"
#include "iniprint.h"

      INTEGER l
      REAL vcov(ijb_v:ije_v,llm),ucov(ijb_u:ije_u,llm) 
      REAL teta(ijb_u:ije_u,llm)                   
      REAL ps(ijb_u:ije_u),masse(ijb_u:ije_u,llm)                   
      REAL q(ijb_u:ije_u,llm,nqtot)
      CHARACTER*(*) fichnom
     
      REAL time
      INTEGER nid, nvarid, nid_trac, nvarid_trac
      REAL trac_tmp(ijb_u:ije_u,llm)      
      INTEGER ierr, ierr_file
      INTEGER iq
      INTEGER length
      PARAMETER (length = 100)
      REAL tab_cntrl(length) ! tableau des parametres du run
      character*20 modname
      character*80 abort_message
c
      INTEGER nb
      SAVE nb
      DATA nb / 0 /
      REAL,SAVE,ALLOCATABLE :: ucov_glo(:,:),vcov_glo(:,:),teta_glo(:,:)
      REAL,SAVE,ALLOCATABLE :: masse_glo(:,:),ps_glo(:),q_glo(:,:)
      LOGICAL,SAVE :: exist_file
      INTEGER,SAVE :: ierr_var
      
!      call Gather_Field(ucov,ip1jmp1,llm,0)
!      call Gather_Field(vcov,ip1jm,llm,0)
!      call Gather_Field(teta,ip1jmp1,llm,0)
!      call Gather_Field(masse,ip1jmp1,llm,0)
!      call Gather_Field(ps,ip1jmp1,1,0)
      
!      do iq=1,nqtot
!        call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
!     enddo
      
!$OMP MASTER      
      if (mpi_rank==0) then
      modname = 'dynredem1_loc'
      ierr = NF_OPEN(fichnom, NF_WRITE, nid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*)"dynredem1: Pb. d ouverture "//trim(fichnom)
         CALL abort_gcm("dynredem1", "", 1)
      ENDIF

c  Ecriture/extension de la coordonnee temps

      nb = nb + 1
      ierr = NF_INQ_VARID(nid, "temps", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         write(lunout,*) NF_STRERROR(ierr)
         abort_message='Variable temps n est pas definie'
         CALL abort_gcm(modname,abort_message,ierr)
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR1_DOUBLE (nid,nvarid,nb,time)
#else
      ierr = NF_PUT_VAR1_REAL (nid,nvarid,nb,time)
#endif
      write(lunout,*) "dynredem1_loc: Enregistrement pour ", nb, time

c
c  Re-ecriture du tableau de controle, itaufin n'est plus defini quand
c  on passe dans dynredem0
      ierr = NF_INQ_VARID (nid, "controle", nvarid)
      IF (ierr .NE. NF_NOERR) THEN
         abort_message="dynredem1: Le champ <controle> est absent"
         ierr = 1
         CALL abort_gcm(modname,abort_message,ierr)
      ENDIF
#ifdef NC_DOUBLE
      ierr = NF_GET_VAR_DOUBLE(nid, nvarid, tab_cntrl)
#else
      ierr = NF_GET_VAR_REAL(nid, nvarid, tab_cntrl)
#endif
       tab_cntrl(31) =  REAL(itau_dyn + itaufin)

       if(sub_day) tab_cntrl(31) = REAL(itaufin) ! DART-LMDZ

#ifdef NC_DOUBLE
      ierr = NF_PUT_VAR_DOUBLE (nid,nvarid,tab_cntrl)
#else
      ierr = NF_PUT_VAR_REAL (nid,nvarid,tab_cntrl)
#endif
      endif
!$OMP END MASTER

!     
      CALL dynredem_write_u(nid,"ucov",ucov,llm)
      CALL dynredem_write_v(nid,"vcov",vcov,llm)
      CALL dynredem_write_u(nid,"teta",teta,llm)
      CALL dynredem_write_u(nid,"masse",masse,llm)
      CALL dynredem_write_u(nid,"ps",ps,1)

      IF (type_trac /= 'inca') THEN
        DO iq=1,nqtot
          CALL dynredem_write_u(nid,tname(iq),q(:,:,iq),llm)
        ENDDO
      ELSE
        
!$OMP MASTER
        INQUIRE(FILE="start_trac.nc", EXIST=exist_file) 
        PRINT *, "EXIST", exist_file
!$OMP END MASTER
!$OMP BARRIER
      
        IF (exist_file) THEN
!$OMP MASTER
          ierr_file = NF_OPEN ("start_trac.nc", NF_NOWRITE,nid_trac)
          IF (ierr_file .NE.NF_NOERR) THEN
            WRITE(6,*)' Pb d''ouverture du fichier start_trac.nc'
            WRITE(6,*)' ierr = ', ierr_file 
          ENDIF
!$OMP END MASTER

          DO iq=1,nqtot

!$OMP MASTER      
            ierr_var = NF_INQ_VARID (nid_trac, tname(iq), nvarid_trac)
!$OMP END MASTER
!$OMP BARRIER
            IF (ierr == NF_NOERR) THEN
              CALL dynredem_read_u(nid_trac,tname(iq),q(:,:,iq),llm) 
            ENDIF
            CALL dynredem_write_u(nid,tname(iq),q(:,:,iq),llm)  
          ENDDO          
          
        ELSE ! pas de fichier start_tract
          DO iq=1,nqtot
            CALL dynredem_write_u(nid,tname(iq),q(:,:,iq),llm)
          ENDDO
        ENDIF 
      ENDIF


!$OMP MASTER
      IF (mpi_rank==0) THEN
        ierr = NF_CLOSE(nid)
      ENDIF ! mpi_rank==0
!$OMP END MASTER
      
      RETURN
      END

