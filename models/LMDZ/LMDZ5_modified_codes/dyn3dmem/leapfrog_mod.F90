MODULE leapfrog_mod

  REAL,POINTER,SAVE :: ucov(:,:) ! zonal covariant wind
  REAL,POINTER,SAVE :: vcov(:,:) ! meridional covariant wind
  REAL,POINTER,SAVE :: teta(:,:) ! potential temperature
  REAL,POINTER,SAVE :: ps(:) ! surface pressure
  REAL,POINTER,SAVE :: masse(:,:) ! air mass
  REAL,POINTER,SAVE :: phis(:) ! geopotential at the surface
  REAL,POINTER,SAVE :: q(:,:,:) ! advected tracers
  REAL,POINTER,SAVE :: p(:,:) ! interlayer pressure
  REAL,POINTER,SAVE :: pks(:) ! Exner at the surface
  REAL,POINTER,SAVE :: pk(:,:) ! Exner at mid-layer
  REAL,POINTER,SAVE :: pkf(:,:) ! filtered Exner
  REAL,POINTER,SAVE :: phi(:,:) ! geopotential
  REAL,POINTER,SAVE :: w(:,:) ! vertical velocity
  REAL,POINTER,SAVE :: pbaru(:,:)
  REAL,POINTER,SAVE :: pbarv(:,:)
  REAL,POINTER,SAVE :: vcovm1(:,:)
  REAL,POINTER,SAVE :: ucovm1(:,:)
  REAL,POINTER,SAVE :: tetam1(:,:)
  REAL,POINTER,SAVE :: psm1(:)
  REAL,POINTER,SAVE :: massem1(:,:)
  REAL,POINTER,SAVE :: dv(:,:)
  REAL,POINTER,SAVE :: du(:,:)
  REAL,POINTER,SAVE :: dteta(:,:)
  REAL,POINTER,SAVE :: dp(:)
  REAL,POINTER,SAVE :: dq(:,:,:)
  REAL,POINTER,SAVE :: finvmaold(:,:)
  REAL,POINTER,SAVE :: flxw(:,:)
  REAL,POINTER,SAVE :: unat(:,:)
  REAL,POINTER,SAVE :: vnat(:,:)
  REAL,POINTER,SAVE :: deltaPress(:,:) !! DART-LMDZ
 

  
CONTAINS

  SUBROUTINE leapfrog_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  USE infotrac
  USE caldyn_mod,ONLY : caldyn_allocate
  USE integrd_mod,ONLY : integrd_allocate
  USE caladvtrac_mod,ONLY : caladvtrac_allocate
  USE call_calfis_mod,ONLY : call_calfis_allocate
  USE call_dissip_mod, ONLY : call_dissip_allocate
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d


    d=>distrib_caldyn
    CALL allocate_u(ucov,llm,d)
    CALL allocate_v(vcov,llm,d)
    CALL allocate_u(teta,llm,d)
    CALL allocate_u(ps,d)
    CALL allocate_u(masse,llm,d)
    CALL allocate_u(phis,d)
    CALL allocate_u(q,llm,nqtot,d)
    CALL allocate_u(p,llmp1,d)
    CALL allocate_u(pks,d)
    CALL allocate_u(pk,llm,d)
    CALL allocate_u(pkf,llm,d)
    CALL allocate_u(phi,llm,d)
    CALL allocate_u(w,llm,d)
    CALL allocate_u(pbaru,llm,d)
    CALL allocate_v(pbarv,llm,d)
    CALL allocate_v(vcovm1,llm,d)
    CALL allocate_u(ucovm1,llm,d)
    CALL allocate_u(tetam1,llm,d)
    CALL allocate_u(psm1,d)
    CALL allocate_u(massem1,llm,d)
    CALL allocate_v(dv,llm,d)
    CALL allocate_u(du,llm,d)
    CALL allocate_u(dteta,llm,d)
    CALL allocate_u(dp,d)
    CALL allocate_u(dq,llm,nqtot,d)
    CALL allocate_u(finvmaold,llm,d)
    CALL allocate_u(flxw,llm,d)
    CALL allocate_u(unat,llm,d)
    CALL allocate_v(vnat,llm,d)
    CALL allocate_u(deltaPress,llm,d)  ! DART-LMDZ
    
    CALL caldyn_allocate
    CALL integrd_allocate
    CALL caladvtrac_allocate
    CALL call_calfis_allocate
    CALL call_dissip_allocate
        
  END SUBROUTINE leapfrog_allocate
  
  SUBROUTINE leapfrog_switch_caldyn(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE caldyn_mod,ONLY : caldyn_switch_caldyn
  USE integrd_mod,ONLY : integrd_switch_caldyn
  USE caladvtrac_mod,ONLY : caladvtrac_switch_caldyn
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL switch_u(ucov,distrib_caldyn,dist)
    CALL switch_v(vcov,distrib_caldyn,dist)
    CALL switch_u(teta,distrib_caldyn,dist)
    CALL switch_u(ps,distrib_caldyn,dist)
    CALL switch_u(masse,distrib_caldyn,dist)
    CALL switch_u(phis,distrib_caldyn,dist,up=halo_max,down=halo_max)
    CALL switch_u(q,distrib_caldyn,dist)
    CALL switch_u(p,distrib_caldyn,dist)
    CALL switch_u(pks,distrib_caldyn,dist)
    CALL switch_u(pk,distrib_caldyn,dist)
    CALL switch_u(pkf,distrib_caldyn,dist)
    CALL switch_u(phi,distrib_caldyn,dist)
    CALL switch_u(w,distrib_caldyn,dist)
    CALL switch_u(pbaru,distrib_caldyn,dist)
    CALL switch_v(pbarv,distrib_caldyn,dist)
    CALL switch_v(vcovm1,distrib_caldyn,dist)
    CALL switch_u(ucovm1,distrib_caldyn,dist)
    CALL switch_u(tetam1,distrib_caldyn,dist)
    CALL switch_u(psm1,distrib_caldyn,dist)
    CALL switch_u(massem1,distrib_caldyn,dist)
    CALL switch_v(dv,distrib_caldyn,dist)
    CALL switch_u(du,distrib_caldyn,dist)
    CALL switch_u(dteta,distrib_caldyn,dist)
    CALL switch_u(dp,distrib_caldyn,dist)
    CALL switch_u(dq,distrib_caldyn,dist)
    CALL switch_u(finvmaold,distrib_caldyn,dist)
    CALL switch_u(flxw,distrib_caldyn,dist)
    CALL switch_u(unat,distrib_caldyn,dist)
    CALL switch_v(vnat,distrib_caldyn,dist)
    CALL switch_u(deltaPress,distrib_caldyn,dist)   ! DART-LMDZ

    
    CALL caldyn_switch_caldyn(dist)
    CALL integrd_switch_caldyn(dist)
    CALL caladvtrac_switch_caldyn(dist)
    
  END SUBROUTINE leapfrog_switch_caldyn
  
  SUBROUTINE leapfrog_switch_dissip(dist)
  USE allocate_field_mod
  USE bands
  USE parallel_lmdz
  USE call_dissip_mod,ONLY : call_dissip_switch_dissip
  IMPLICIT NONE
    TYPE(distrib),INTENT(IN) :: dist

    CALL call_dissip_switch_dissip(dist)
    
  END SUBROUTINE leapfrog_switch_dissip
  
END MODULE leapfrog_mod  







