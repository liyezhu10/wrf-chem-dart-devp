! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! The following observation kinds are currently created by the converters
! in the NASA_Earthdata directory. As NASA provides them
!	soil_moisture_x:long_name = "Volumetric Soil Moisture from X-band" ;
!	soil_moisture_x:units = "percent" ;
! CLM soil moisture in the restart files is another beast:
!       H2OSOI_LIQ:long_name = "liquid water" ;
!       H2OSOI_LIQ:units = "kg/m2" ;
! CLM soil moisture in the history files is another beast:
!       H2OSOI:long_name = "liquid water" ;
!       H2OSOI:units = "mm3/mm3" ;

!SMOS_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMOS_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMAP_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SMAP_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!SSMI_A_SOIL_MOISTURE,         volumetric soil moisture  percent
!SSMI_D_SOIL_MOISTURE,         volumetric soil moisture  percent
!AMSRE_A_SOIL_MOISTURE_X,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_D_SOIL_MOISTURE_X,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_A_SOIL_MOISTURE_C,      volumetric soil moisture  percent (cvrts to a fraction)
!AMSRE_D_SOIL_MOISTURE_C,      volumetric soil moisture  percent (cvrts to a fraction)
!TRMM_SOIL_MOISTURE,           volumetric soil moisture  percent
!WINDSAT_SOIL_MOISTURE_X,      volumetric soil moisture  percent
!WINDSAT_SOIL_MOISTURE_C,      volumetric soil moisture  percent

! BEGIN DART PREPROCESS KIND LIST
!SOIL_TEMPERATURE,             QTY_SOIL_TEMPERATURE,   COMMON_CODE
!LPRM_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,      COMMON_CODE
!SMOS_A_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!SMOS_D_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!SMAP_A_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!SMAP_D_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!SSMI_A_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!SSMI_D_SOIL_MOISTURE,         QTY_SOIL_MOISTURE,      COMMON_CODE
!AMSRE_A_SOIL_MOISTURE_X,      QTY_SOIL_MOISTURE,      COMMON_CODE
!AMSRE_D_SOIL_MOISTURE_X,      QTY_SOIL_MOISTURE,      COMMON_CODE
!AMSRE_A_SOIL_MOISTURE_C,      QTY_SOIL_MOISTURE,      COMMON_CODE
!AMSRE_D_SOIL_MOISTURE_C,      QTY_SOIL_MOISTURE,      COMMON_CODE
!TRMM_SOIL_MOISTURE,           QTY_SOIL_MOISTURE,      COMMON_CODE
!WINDSAT_SOIL_MOISTURE_X,      QTY_SOIL_MOISTURE,      COMMON_CODE
!WINDSAT_SOIL_MOISTURE_C,      QTY_SOIL_MOISTURE,      COMMON_CODE
! END DART PREPROCESS KIND LIST

