#!/bin/ksh -aeux
#########################################################################
#
# Purpose: EMISS PERTURB
#
#########################################################################
#
# CODE VERSIONS
export WRFDA_VER=WRFDAv3.4_dmpar
export WRFDA_TOOLS_VER=WRFDA_TOOLSv3.4
export WRF_VER=WRFv3.6.1_dmpar
export WRFCHEM_VER=WRFCHEMv3.6.1_dmpar
export DART_VER=DART_CHEM_MY_BRANCH_FRAPPE
#
# EXPERIMENT DETAILS:
export NUM_MEMBERS=20
export MEM_START=1
export EXP=real_FRAPPE_CNTL_VARLOC_REV1.a
#
# TIME DATA:
export START_DATE=2014071606
export END_DATE=2014071606
export DATE=${START_DATE}
#
# CYCLING AND BOUNDARY CONDITION TIME DATA
export EMISS_FREQ=6
export EMISS_FREQ=1
#
# DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export HOME_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
export FRAPPE_DIR=/glade/p/FRAPPE
#
export TRUNK_DIR=${HOME_DIR}/TRUNK
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}${WRF_VER}
export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
export ANTHRO_DATA_DIR=${FRAPPE_DIR}/REAL_TIME_DATA/anthro_emissions
export FIRE_DATA_DIR=${FRAPPE_DIR}/REAL_TIME_DATA/fire_emissions
export BIO_DATA_DIR=${FRAPPE_DIR}/REAL_TIME_DATA/bio_emissions
#
# LOOP THROUGH DATES
while [[ ${DATE} -le ${END_DATE} ]] ; do
   export YYYY=$(echo $DATE | cut -c1-4)
   export MM=$(echo $DATE | cut -c5-6)
   export DD=$(echo $DATE | cut -c7-8)
   export HH=$(echo $DATE | cut -c9-10)
   export RUN_DIR=${FRAPPE_DIR}/${EXP}/${DATE}/wrfchem_chem_emiss
#
# COPY RECENTER CODE
   if [[ -e recenter_chem_emiss_RT_CR_CORR.exe ]]; then 
      rm -rf recenter_chem_emiss_RT_CR_CORR.exe; 
   fi
   cp ${DART_DIR}/recenter_chem_emiss_RT_CR_CORR.exe ./.
#
   export NL_PERT_CHEM=.true.
   export NL_PERT_FIRE=.true.
   export NL_PERT_BIO=.true.
   export WRFCHEMI=wrfchemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
   export WRFFIRECHEMI=wrffirechemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
   export WRFBIOCHEMI=wrfbiochemi_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
   export WRFCHEMI_PARENT=${ANTHRO_DATA_DIR}/${WRFCHEMI}
   export WRFFIRECHEMI_PARENT=${FIRE_DATA_DIR}/${WRFFIRECHEMI}
   export WRFBIOCHEMI_PARENT=${BIO_DATA_DIR}/${WRFBIOCHEMI}
#
# CREATE NAMELIST
   rm -rf recenter_chem_emiss_RT_CR_CORR_nml.nl
   cat << EOF > recenter_chem_emiss_RT_CR_CORR_nml.nl
&recenter_chem_emiss_RT_CR__CORR_nml
wrfchemi='${WRFCHEMI}',
wrffirechemi='${WRFFIRECHEMI}',
wrfbiochemi='${WRFBIOCHEMI}',
wrfchemi_parent='${WRFCHEMI_PARENT}',
wrffirechemi_parent='${WRFFIRECHEMI_PARENT}',
wrfbiochemi_parent='${WRFBIOCHEMI_PARENT}',
sw_chem=${NL_PERT_CHEM},
sw_fire=${NL_PERT_FIRE},
sw_bio=${NL_PERT_BIO},
/
EOF
#
# RUN PERTRUBATION CODE    
      ./recenter_chem_emiss_RT_CR_CORR.exe > temp.out
#
# ADVANCE TIME      
   export DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${EMISS_FREQ} 2>/dev/null)
done
