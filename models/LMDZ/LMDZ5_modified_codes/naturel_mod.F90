!**************************************************************************
!   Tarkeshwar Singh
!   IIT Delhi
!   Email: tarkphysics87@gmail.com
!**************************************************************************

MODULE naturel_mod
 REAL,POINTER,SAVE :: preslev (:,:  )
 REAL,POINTER,SAVE :: pkss(:) 
 REAL,POINTER,SAVE :: pkreff(:,:) 
CONTAINS

  SUBROUTINE naturel_allocate
  USE bands
  USE allocate_field_mod
  USE parallel_lmdz
  USE dimensions_mod
  IMPLICIT NONE
  TYPE(distrib),POINTER :: d

    d=>distrib_caldyn

    CALL allocate_u(pkss,d)
    CALL allocate_u(pkreff,llm,d)
    CALL allocate_u(preslev,llmp1,d)

  END SUBROUTINE naturel_allocate

END MODULE naturel_mod


 
