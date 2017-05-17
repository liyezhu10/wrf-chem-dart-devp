!
! $Id$
!
!-------------------------------------------------------------------------------
!
SUBROUTINE etat0_netcdf(ib, masque, phis, letat0)
!
!-------------------------------------------------------------------------------
! Purpose: Creates initial states
!-------------------------------------------------------------------------------
! Note: This routine is designed to work for Earth
!-------------------------------------------------------------------------------
  USE control_mod
#ifdef CPP_EARTH
  USE startvar
  USE ioipsl
  USE dimphy
  USE infotrac
  USE fonte_neige_mod
  USE pbl_surface_mod
  USE phys_state_var_mod
  USE filtreg_mod
  USE regr_lat_time_climoz_m, ONLY: regr_lat_time_climoz
  USE conf_phys_m,            ONLY: conf_phys
! For parameterization of ozone chemistry:
  USE regr_lat_time_coefoz_m, only: regr_lat_time_coefoz
  USE press_coefoz_m, only: press_coefoz
  USE regr_pr_o3_m, only: regr_pr_o3
  USE netcdf, ONLY : NF90_OPEN, NF90_NOWRITE, NF90_CLOSE, NF90_NOERR
  USE indice_sol_mod
  use exner_hyb_m, only: exner_hyb
  use exner_milieu_m, only: exner_milieu
  use test_disvert_m, only: test_disvert
#endif
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Arguments:
#include "dimensions.h"
#include "paramet.h"
#include "iniprint.h"
  LOGICAL,                    INTENT(IN)    :: ib     ! barycentric interpolat.
  REAL, DIMENSION(iip1,jjp1), INTENT(INOUT) :: masque ! land mask
  REAL, DIMENSION(iip1,jjp1), INTENT(OUT)   :: phis   ! geopotentiel au sol
  LOGICAL,                    INTENT(IN)    :: letat0 ! F: masque only required
#ifndef CPP_EARTH
  WRITE(lunout,*)'limit_netcdf: Earth-specific routine, needs Earth physics'
#else
!-------------------------------------------------------------------------------
! Local variables:
#include "comgeom2.h"
#include "comvert.h"
#include "comconst.h"
#include "dimsoil.h"
#include "temps.h"
  REAL,    DIMENSION(klon)                 :: tsol, qsol
  REAL,    DIMENSION(klon)                 :: sn, rugmer, run_off_lic_0
  REAL,    DIMENSION(iip1,jjp1)            :: orog, rugo, psol
  REAL,    DIMENSION(iip1,jjp1,llm+1)      :: p3d
  REAL,    DIMENSION(iip1,jjp1,llm)        :: uvent, t3d, tpot, qsat, qd
  REAL,    DIMENSION(iip1,jjm ,llm)        :: vvent
  REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: q3d
  REAL,    DIMENSION(klon,nbsrf)           :: qsolsrf, snsrf, evap
  REAL,    DIMENSION(klon,nbsrf)           :: frugs, agesno
  REAL,    DIMENSION(klon,nsoilmx,nbsrf)   :: tsoil

!--- Local variables for sea-ice reading:
  INTEGER                                  :: iml_lic, jml_lic, llm_tmp
  INTEGER                                  :: ttm_tmp, iret, fid
  INTEGER, DIMENSION(1)                    :: itaul
  REAL,    DIMENSION(1)                    :: lev
  REAL                                     :: date
  REAL,    DIMENSION(:,:),   ALLOCATABLE   ::  lon_lic,  lat_lic, fraclic
  REAL,    DIMENSION(:),     ALLOCATABLE   :: dlon_lic, dlat_lic
  REAL,    DIMENSION(iip1,jjp1)            :: flic_tmp

!--- Misc
  CHARACTER(LEN=80)                        :: x, fmt
  INTEGER                                  :: i, j, l, ji
  REAL,    DIMENSION(iip1,jjp1,llm)        :: pk, pls, y
  REAL,    DIMENSION(ip1jmp1)              :: pks

#include "comdissnew.h"
#include "serre.h"
#include "clesphys.h"

  REAL,    DIMENSION(iip1,jjp1,llm)        :: masse
  INTEGER :: itau, iday
  REAL    :: xpn, xps, time, phystep
  REAL,    DIMENSION(iim)                  :: xppn, xpps
  REAL,    DIMENSION(ip1jmp1,llm)          :: pbaru, phi, w
  REAL,    DIMENSION(ip1jm  ,llm)          :: pbarv
  REAL,    DIMENSION(klon)                 :: fder

!--- Local variables for ocean mask reading:
  INTEGER :: nid_o2a, iml_omask, jml_omask
  LOGICAL :: couple=.FALSE.
  REAL,    DIMENSION(:,:), ALLOCATABLE ::  lon_omask, lat_omask, ocemask, ocetmp
  REAL,    DIMENSION(:),   ALLOCATABLE :: dlon_omask,dlat_omask
  REAL,    DIMENSION(klon)             :: ocemask_fi
  INTEGER, DIMENSION(klon-2)           :: isst
  REAL,    DIMENSION(iim,jjp1)         :: zx_tmp_2d
  REAL    :: dummy
  LOGICAL :: ok_newmicro, ok_journe, ok_mensuel, ok_instan, ok_hf
  LOGICAL :: ok_LES, ok_ade, ok_aie, ok_cdnc, aerosol_couple, new_aod, callstats
  INTEGER :: iflag_radia, flag_aerosol
  LOGICAL :: flag_aerosol_strat
  REAL    :: bl95_b0, bl95_b1, fact_cldcon, facttemps, ratqsbas, ratqshaut
  REAL    :: tau_ratqs
  INTEGER :: iflag_cldcon, iflag_ratqs, iflag_coupl, iflag_clos, iflag_wake
  INTEGER :: iflag_thermals, nsplit_thermals
  INTEGER :: iflag_thermals_ed, iflag_thermals_optflux
  REAL    :: tau_thermals, solarlong0,  seuil_inversion
  INTEGER :: read_climoz ! read ozone climatology
!  Allowed values are 0, 1 and 2
!     0: do not read an ozone climatology
!     1: read a single ozone climatology that will be used day and night
!     2: read two ozone climatologies, the average day and night
!     climatology and the daylight climatology
!-------------------------------------------------------------------------------
  REAL    :: alp_offset
  logical found

!--- Constants
  pi     = 4. * ATAN(1.)
  rad    = 6371229.
  daysec = 86400.
  omeg   = 2.*pi/daysec
  g      = 9.8
  kappa  = 0.2857143
  cpp    = 1004.70885
  preff  = 101325.
  pa     = 50000.
  jmp1   = jjm + 1

!--- CONSTRUCT A GRID
  CALL conf_phys(  ok_journe, ok_mensuel, ok_instan, ok_hf, ok_LES,     &
                   callstats,                                           &
                   solarlong0,seuil_inversion,                          &
                   fact_cldcon, facttemps,ok_newmicro,iflag_radia,      &
                   iflag_cldcon,                                        &
                   iflag_ratqs,ratqsbas,ratqshaut,tau_ratqs,            &
                   ok_ade, ok_aie, ok_cdnc, aerosol_couple,             &
                   flag_aerosol, flag_aerosol_strat, new_aod,           &
                   bl95_b0, bl95_b1,                                    &
                   read_climoz,                                         &
                   alp_offset)

! co2_ppm0 : initial value of atmospheric CO2 from .def file (co2_ppm value)
  co2_ppm0 = co2_ppm

  dtvr   = daysec/FLOAT(day_step)
  WRITE(lunout,*)'dtvr',dtvr

  CALL iniconst()
  if (pressure_exner) call test_disvert
  CALL inigeom()

!--- Initializations for tracers
  CALL infotrac_init
  ALLOCATE(q3d(iip1,jjp1,llm,nqtot))

  CALL inifilr()
  CALL phys_state_var_init(read_climoz)

  rlat(1) = ASIN(1.0)
  DO j=2,jjm; rlat((j-2)*iim+2:(j-1)*iim+1)=rlatu(j);     END DO
  rlat(klon) = - ASIN(1.0)
  rlat(:)=rlat(:)*(180.0/pi)

  rlon(1) = 0.0
  DO j=2,jjm; rlon((j-2)*iim+2:(j-1)*iim+1)=rlonv(1:iim); END DO
  rlon(klon) = 0.0
  rlon(:)=rlon(:)*(180.0/pi)

! For a coupled simulation, the ocean mask from ocean model is used to compute
! the weights an to insure ocean fractions are the same for atmosphere and ocean
! Otherwise, mask is created using Relief file.

  WRITE(lunout,*)'Essai de lecture masque ocean'
  iret = NF90_OPEN("o2a.nc", NF90_NOWRITE, nid_o2a)
  IF(iret/=NF90_NOERR) THEN
    WRITE(lunout,*)'ATTENTION!! pas de fichier o2a.nc trouve'
    WRITE(lunout,*)'Run force'
    x='masque'
    masque(:,:)=0.0
    CALL startget_phys2d(x, iip1, jjp1, rlonv, rlatu, masque, 0.0, jjm, &
   &              rlonu, rlatv, ib)
    WRITE(lunout,*)'MASQUE construit : Masque'
    if (prt_level >= 1) WRITE(lunout,'(97I1)') nINT(masque)
    CALL gr_dyn_fi(1, iip1, jjp1, klon, masque, zmasq)
    WHERE(   zmasq(:)<EPSFRA) zmasq(:)=0.
    WHERE(1.-zmasq(:)<EPSFRA) zmasq(:)=1.
  ELSE
    WRITE(lunout,*)'ATTENTION!! fichier o2a.nc trouve'
    WRITE(lunout,*)'Run couple'
    couple=.true.
    iret=NF90_CLOSE(nid_o2a)
    CALL flininfo("o2a.nc", iml_omask, jml_omask, llm_tmp, ttm_tmp, nid_o2a)
    IF(iml_omask/=iim .OR.jml_omask/=jjp1) THEN
      WRITE(lunout,*)'Dimensions non compatibles pour masque ocean'
      WRITE(lunout,*)'iim = ',iim,' iml_omask = ',iml_omask
      WRITE(lunout,*)'jjp1 = ',jjp1,' jml_omask = ',jml_omask
      CALL abort_gcm('etat0_netcdf','Dimensions non compatibles pour masque oc&
     &ean',1)
    END IF
    ALLOCATE(   ocemask(iml_omask,jml_omask),   ocetmp(iml_omask,jml_omask))
    ALLOCATE( lon_omask(iml_omask,jml_omask),lat_omask(iml_omask,jml_omask))
    ALLOCATE(dlon_omask(iml_omask),         dlat_omask(jml_omask))
    CALL flinopen("o2a.nc", .FALSE., iml_omask, jml_omask, llm_tmp, lon_omask,&
   &              lat_omask, lev, ttm_tmp, itaul, date, dt, fid)
    CALL flinget(fid, 'OceMask', iml_omask, jml_omask, llm_tmp, ttm_tmp, &
   &              1, 1, ocetmp)
    CALL flinclo(fid)
    dlon_omask(1:iml_omask) = lon_omask(1:iml_omask,1)
    dlat_omask(1:jml_omask) = lat_omask(1,1:jml_omask)
    ocemask = ocetmp
    IF(dlat_omask(1)<dlat_omask(jml_omask)) THEN
      DO j=1,jml_omask
        ocemask(:,j) = ocetmp(:,jml_omask-j+1)
      END DO
    END IF
!
! Ocean mask to physical grid
!*******************************************************************************
    WRITE(lunout,*)'ocemask '
    WRITE(fmt,"(i4,'i1)')")iml_omask ; fmt='('//ADJUSTL(fmt)
    WRITE(lunout,fmt)int(ocemask)
    ocemask_fi(1)=ocemask(1,1)
    DO j=2,jjm; ocemask_fi((j-2)*iim+2:(j-1)*iim+1)=ocemask(1:iim,j); END DO
    ocemask_fi(klon)=ocemask(1,jjp1)
    zmasq=1.-ocemask_fi
  END IF

  CALL gr_fi_dyn(1,klon,iip1,jjp1,zmasq,masque)

  ! The startget calls need to be replaced by a call to restget to get the
  ! values in the restart file
  x = 'relief'; orog(:,:) = 0.0
  CALL startget_phys2d(x,iip1,jjp1,rlonv,rlatu, orog, 0.0,jjm,rlonu,rlatv,ib,&
 &               masque)

  x = 'rugosite'; rugo(:,:) = 0.0
  CALL startget_phys2d(x,iip1,jjp1,rlonv,rlatu, rugo, 0.0,jjm, rlonu,rlatv,ib)
!  WRITE(lunout,'(49I1)') INT(orog(:,:)*10)
!  WRITE(lunout,'(49I1)') INT(rugo(:,:)*10)

! Sub-surfaces initialization
!*******************************************************************************
  pctsrf=0.
  x = 'psol'; psol(:,:) = 0.0
  CALL startget_phys2d(x,iip1,jjp1,rlonv,rlatu,psol,0.0,jjm,rlonu,rlatv,ib)
!  WRITE(lunout,*) 'PSOL :', psol(10,20)
!  WRITE(lunout,*) ap(:), bp(:)

! Mid-levels pressure computation
!*******************************************************************************
  CALL pression(ip1jmp1, ap, bp, psol, p3d)
  if (pressure_exner) then
    CALL exner_hyb(ip1jmp1, psol, p3d, pks, pk)
  else
    CALL exner_milieu(ip1jmp1,psol,p3d, pks,pk)
  endif
  pls(:,:,:)=preff*(pk(:,:,:)/cpp)**(1./kappa)
!  WRITE(lunout,*) 'P3D :', p3d(10,20,:)
!  WRITE(lunout,*) 'PK:',    pk(10,20,:)
!  WRITE(lunout,*) 'PLS :', pls(10,20,:)

  x = 'surfgeo'; phis(:,:) = 0.0
  CALL startget_phys2d(x,iip1,jjp1,rlonv,rlatu,phis, 0.0,jjm, rlonu,rlatv,ib)

  x = 'u';    uvent(:,:,:) = 0.0
  CALL startget_dyn(x,rlonu,rlatu,pls,y,uvent,0.0,  &
 &                  rlonv,rlatv,ib)

  x = 'v';    vvent(:,:,:) = 0.0
  CALL startget_dyn(x, rlonv,rlatv,pls(:, :jjm, :),y(:, :jjm, :),vvent,0.0, &
 &                  rlonu,rlatu(:jjm),ib)

  x = 't';    t3d(:,:,:) = 0.0
  CALL startget_dyn(x,rlonv,rlatu,pls,y,t3d,0.0,    &
 &                  rlonu,rlatv,ib)

  x = 'tpot'; tpot(:,:,:) = 0.0
  CALL startget_dyn(x,rlonv,rlatu,pls,pk,tpot,0.0,  &
 &                  rlonu,rlatv,ib)

  WRITE(lunout,*) 'T3D min,max:',minval(t3d(:,:,:)),maxval(t3d(:,:,:))
  WRITE(lunout,*) 'PLS min,max:',minval(pls(:,:,:)),maxval(pls(:,:,:))

! Humidity at saturation computation
!*******************************************************************************
  WRITE(lunout,*) 'avant q_sat'
  CALL q_sat(llm*jjp1*iip1, t3d, pls, qsat)
  WRITE(lunout,*) 'apres q_sat'
  WRITE(lunout,*) 'QSAT min,max:',minval(qsat(:,:,:)),maxval(qsat(:,:,:))
!  WRITE(lunout,*) 'QSAT :',qsat(10,20,:)

  x = 'q';    qd (:,:,:) = 0.0
  CALL startget_dyn(x,rlonv,rlatu,pls,qsat,qd,0.0, rlonu,rlatv,ib)
  q3d(:,:,:,:) = 0.0 ; q3d(:,:,:,1) = qd(:,:,:)

! Parameterization of ozone chemistry:
! Look for ozone tracer:
  i = 1
  DO
    found = tname(i)=="O3" .OR. tname(i)=="o3"
    if (found .or. i == nqtot) exit
    i = i + 1
  end do
  if (found) then
    call regr_lat_time_coefoz
    call press_coefoz
    call regr_pr_o3(p3d, q3d(:, :, :, i))
!   Convert from mole fraction to mass fraction:
    q3d(:, :, :, i) = q3d(:, :, :, i)  * 48. / 29.
  end if

!--- OZONE CLIMATOLOGY
  IF(read_climoz>=1) CALL regr_lat_time_climoz(read_climoz)

  x = 'tsol'; tsol(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,tsol,0.0,jjm,rlonu,rlatv,ib)

  x = 'qsol';  qsol(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,qsol,0.0,jjm,rlonu,rlatv,ib)

  x = 'snow';  sn(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,sn,0.0,jjm,rlonu,rlatv,ib)

  x = 'rads';  radsol(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,radsol,0.0,jjm,rlonu,rlatv,ib)

  x = 'rugmer'; rugmer(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,rugmer,0.0,jjm,rlonu,rlatv,ib)

  x = 'zmea';  zmea(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zmea,0.0,jjm,rlonu,rlatv,ib)

  x = 'zstd';  zstd(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zstd,0.0,jjm,rlonu,rlatv,ib)

  x = 'zsig';  zsig(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zsig,0.0,jjm,rlonu,rlatv,ib)

  x = 'zgam';  zgam(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zgam,0.0,jjm,rlonu,rlatv,ib)

  x = 'zthe';  zthe(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zthe,0.0,jjm,rlonu,rlatv,ib)

  x = 'zpic';  zpic(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zpic,0.0,jjm,rlonu,rlatv,ib)

  x = 'zval';  zval(:) = 0.0
  CALL startget_phys1d(x,iip1,jjp1,rlonv,rlatu,klon,zval,0.0,jjm,rlonu,rlatv,ib)

!  WRITE(lunout,'(48I3)') 'TSOL :', INT(tsol(2:klon)-273)

! Soil ice file reading for soil fraction and soil ice fraction
!*******************************************************************************
  CALL flininfo("landiceref.nc", iml_lic, jml_lic, llm_tmp, ttm_tmp, fid)
  ALLOCATE( lat_lic(iml_lic,jml_lic),lon_lic(iml_lic, jml_lic))
  ALLOCATE(dlat_lic(jml_lic),       dlon_lic(iml_lic))
  ALLOCATE( fraclic(iml_lic,jml_lic))
  CALL flinopen("landiceref.nc", .FALSE., iml_lic, jml_lic, llm_tmp,  &
 &               lon_lic, lat_lic, lev, ttm_tmp, itaul, date, dt, fid)
  CALL flinget(fid, 'landice', iml_lic, jml_lic, llm_tmp, ttm_tmp, 1,1, fraclic)
  CALL flinclo(fid)

! Interpolation on model T-grid
!*******************************************************************************
  WRITE(lunout,*)'dimensions de landice iml_lic, jml_lic : ',iml_lic,jml_lic
! conversion if coordinates are in degrees
  IF(MAXVAL(lon_lic)>pi) lon_lic=lon_lic*pi/180.
  IF(MAXVAL(lat_lic)>pi) lat_lic=lat_lic*pi/180.
  dlon_lic(:)=lon_lic(:,1)
  dlat_lic(:)=lat_lic(1,:)
  CALL grille_m( iml_lic, jml_lic, dlon_lic, dlat_lic, fraclic, iim,jjp1,   &
 &               rlonv, rlatu, flic_tmp(1:iim,:) )
  flic_tmp(iip1,:)=flic_tmp(1,:)

!--- To the physical grid
  CALL gr_dyn_fi(1, iip1, jjp1, klon, flic_tmp, pctsrf(:,is_lic))

!--- Adequation with soil/sea mask
  WHERE(pctsrf(:,is_lic)<EPSFRA) pctsrf(:,is_lic)=0. 
  WHERE(zmasq(:)<EPSFRA)         pctsrf(:,is_lic)=0.
  pctsrf(:,is_ter)=zmasq(:)
  DO ji=1,klon
    IF(zmasq(ji)>EPSFRA) THEN 
      IF(pctsrf(ji,is_lic)>=zmasq(ji)) THEN
        pctsrf(ji,is_lic)=zmasq(ji)
        pctsrf(ji,is_ter)=0.
      ELSE
        pctsrf(ji,is_ter)=zmasq(ji)-pctsrf(ji,is_lic)
        IF(pctsrf(ji,is_ter)<EPSFRA) THEN
          pctsrf(ji,is_ter)=0.
          pctsrf(ji,is_lic)=zmasq(ji)
        END IF 
      END IF 
    END IF 
  END DO 

! sub-surface ocean and sea ice (sea ice set to zero for start)
!*******************************************************************************
  pctsrf(:,is_oce)=(1.-zmasq(:))
  WHERE(pctsrf(:,is_oce)<EPSFRA) pctsrf(:,is_oce)=0.
  IF(couple) pctsrf(:,is_oce)=ocemask_fi(:)
  isst=0
  WHERE(pctsrf(2:klon-1,is_oce)>0.) isst=1

! It is checked that the sub-surfaces sum is equal to 1
!*******************************************************************************
  ji=COUNT((ABS(SUM(pctsrf(:,:),dim=2))-1.0)>EPSFRA)
  IF(ji/=0) WRITE(lunout,*) 'pb repartition sous maille pour ',ji,' points'
  CALL gr_fi_ecrit(1, klon, iim, jjp1, zmasq, zx_tmp_2d)
!  WRITE(fmt,"(i3,')')")iim; fmt='(i'//ADJUSTL(fmt)
!  WRITE(lunout,*)'zmasq = '
!  WRITE(lunout,TRIM(fmt))NINT(zx_tmp_2d)
  CALL gr_fi_dyn(1, klon, iip1, jjp1, zmasq, masque)
  WRITE(fmt,"(i4,'i1)')")iip1 ; fmt='('//ADJUSTL(fmt)
  WRITE(lunout,*) 'MASQUE construit : Masque'
  if (prt_level >= 1) WRITE(lunout,TRIM(fmt))NINT(masque(:,:))

! Intermediate computation
!*******************************************************************************
  CALL massdair(p3d,masse)
  WRITE(lunout,*)' ALPHAX ',alphax
  DO l=1,llm
    xppn(:)=aire(1:iim,1   )*masse(1:iim,1   ,l)
    xpps(:)=aire(1:iim,jjp1)*masse(1:iim,jjp1,l)
    xpn=SUM(xppn)/apoln
    xps=SUM(xpps)/apols
    masse(:,1   ,l)=xpn
    masse(:,jjp1,l)=xps
  END DO
  q3d(iip1,:,:,:)=q3d(1,:,:,:)
  phis(iip1,:) = phis(1,:)

  IF(.NOT.letat0) RETURN

! Writing
!*******************************************************************************
  CALL inidissip(lstardis, nitergdiv, nitergrot, niterh, tetagdiv, tetagrot, &
       tetatemp, vert_prof_dissip)
  WRITE(lunout,*)'sortie inidissip'
  itau=0
  itau_dyn=0
  itau_phy=0
  iday=dayref+itau/day_step
  time=FLOAT(itau-(iday-dayref)*day_step)/day_step
  IF(time>1.) THEN
   time=time-1
   iday=iday+1
  END IF
  day_ref=dayref
  annee_ref=anneeref

  CALL geopot( ip1jmp1, tpot, pk, pks, phis, phi )
  WRITE(lunout,*)'sortie geopot'

  CALL caldyn0( itau, uvent, vvent, tpot, psol, masse, pk, phis,               &
                phi,  w, pbaru, pbarv, time+iday-dayref)
  WRITE(lunout,*)'sortie caldyn0'     
  CALL dynredem0( "start.nc", dayref, phis)
  WRITE(lunout,*)'sortie dynredem0'

!-------DART-LMDZ------------------
  if(write_restart_natural) call naturel(uvent,vvent,tpot,uvent,vvent,tpot,0,psol)
  CALL dynredem1( "start.nc", 0.0, vvent, uvent, tpot, q3d, masse, psol)
  if(write_restart_natural) call naturel(uvent,vvent,tpot,uvent,vvent,tpot,1,psol) 
!-------DART-LMDZ------------------

  WRITE(lunout,*)'sortie dynredem1' 

! Physical initial state writting
!*******************************************************************************
  WRITE(lunout,*)'phystep ',dtvr,iphysiq,nbapp_rad
  phystep   = dtvr * FLOAT(iphysiq)
  radpas    = NINT (86400./phystep/ FLOAT(nbapp_rad) )
  WRITE(lunout,*)'phystep =', phystep, radpas

! Init: tsol, qsol, sn, evap, tsoil, rain_fall, snow_fall, solsw, sollw, frugs
!*******************************************************************************
  DO i=1,nbsrf; ftsol(:,i) = tsol; END DO
  DO i=1,nbsrf; snsrf(:,i) = sn;   END DO
  falb1(:,is_ter) = 0.08; falb1(:,is_lic) = 0.6
  falb1(:,is_oce) = 0.5;  falb1(:,is_sic) = 0.6
  falb2 = falb1
  evap(:,:) = 0.
  DO i=1,nbsrf; qsolsrf(:,i)=150.; END DO
  DO i=1,nbsrf; DO j=1,nsoilmx; tsoil(:,j,i) = tsol; END DO; END DO
  rain_fall = 0.; snow_fall = 0.
  solsw = 165.;   sollw = -53.
  t_ancien = 273.15
  q_ancien = 0.
  agesno = 0.
  frugs(:,is_oce) = rugmer(:)
  frugs(:,is_ter) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)
  frugs(:,is_lic) = MAX(1.0e-05,zstd(:)*zsig(:)/2.0)
  frugs(:,is_sic) = 0.001
  fder = 0.0
  clwcon = 0.0
  rnebcon = 0.0
  ratqs = 0.0
  run_off_lic_0 = 0.0 
  rugoro = 0.0

! Before phyredem calling, surface modules and values to be saved in startphy.nc
! are initialized
!*******************************************************************************
  dummy = 1.0
  pbl_tke(:,:,:) = 1.e-8 
  zmax0(:) = 40.
  f0(:) = 1.e-5
  sig1(:,:) = 0.
  w01(:,:) = 0.
  wake_deltat(:,:) = 0.
  wake_deltaq(:,:) = 0.
  wake_s(:) = 0.
  wake_cstar(:) = 0.
  wake_fip(:) = 0.
  wake_pe = 0.
  fm_therm = 0.
  entr_therm = 0.
  detr_therm = 0.

  CALL fonte_neige_init(run_off_lic_0)
  CALL pbl_surface_init( qsol, fder, snsrf, qsolsrf, evap, frugs, agesno, tsoil )
  CALL phyredem( "startphy.nc" )

!  WRITE(lunout,*)'CCCCCCCCCCCCCCCCCC REACTIVER SORTIE VISU DANS ETAT0'
!  WRITE(lunout,*)'entree histclo'
  CALL histclo()

#endif 
!#endif of #ifdef CPP_EARTH
  RETURN

END SUBROUTINE etat0_netcdf
!
!-------------------------------------------------------------------------------
