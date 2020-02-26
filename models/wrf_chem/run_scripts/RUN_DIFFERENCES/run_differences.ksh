#!/bin/ksh -aeux
#########################################################################
#
# Purpose: EMISS PERTURB
#
#########################################################################
#
# CODE VERSIONS
# TIME DATA:
export STR_DATE=2014071406
export END_DATE=2014071700
export CYCLE_TIME=6
export MEMBER=1
export RUN_DIR=run_${MEMBER}
if [[ ${MEMBER} -lt 100 ]]; then export RUN_DIR=run_e0${MEMBER}; fi
if [[ ${MEMBER} -lt 10  ]]; then export RUN_DIR=run_e00${MEMBER}; fi
#
export SCRATCH_DIR=/scratch/summit/mizzi
export PROJECTS_DIR=/projects/mizzi
export TRUNK_DIR=/projects/mizzi/TRUNK
#
export WRFDA_VER=WRFDAv3.9.1.1_dmpar
export WRFVAR_DIR=${TRUNK_DIR}/${WRFDA_VER}
export BUILD_DIR=${WRFVAR_DIR}/var/da
export CNTL_DIR=${SCRATCH_DIR}/real_FRAPPE_CONTROL
export RETR_DIR=${SCRATCH_DIR}/real_FRAPPE_RETR_MOP_CO
export RETR_DIR=${SCRATCH_DIR}/real_FRAPPE_CPSR_MOP_CO
export DIFF_DIR=${SCRATCH_DIR}/DIFF_PLOTS
if [[ ! -d ${DIFF_DIR} ]]; then mkdir ${DIFF_DIR}; fi
cd ${DIFF_DIR}
#
# LOOP THROUGH DATES
export DATE=${STR_DATE}
while [[ ${DATE} -le ${END_DATE} ]] ; do
   export CYCLE_STR=${DATE}
   export CYCLE_END=$(${BUILD_DIR}/da_advance_time.exe ${CYCLE_STR} ${CYCLE_TIME} 2>/dev/null)
   export CYCLE_DATE=${CYCLE_STR}
   while [[ ${CYCLE_DATE} -le ${CYCLE_END} ]] ; do   
      export YYYY=$(echo ${CYCLE_DATE} | cut -c1-4)
      export MM=$(echo ${CYCLE_DATE} | cut -c5-6)
      export DD=$(echo ${CYCLE_DATE} | cut -c7-8)
      export HH=$(echo ${CYCLE_DATE} | cut -c9-10)
#
# DIFFERENCE FILES
      export CNTL_FILE=${CNTL_DIR}/${CYCLE_STR}/wrfchem_cycle_cr/${RUN_DIR}/wrfout_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
      export RETR_FILE=${RETR_DIR}/${CYCLE_STR}/wrfchem_cycle_cr/${RUN_DIR}/wrfout_d01_${YYYY}-${MM}-${DD}_${HH}:00:00
      export DIFF_FORE_FILE=Diff_${MM}${DD}${HH}_AFOR.nc
      export DIFF_ANAL_FILE=Diff_${MM}${DD}${HH}_ZANL.nc
      rm -rf diff_file.nc
      ncdiff ${RETR_FILE} ${CNTL_FILE} diff_file.nc
      if [[ ${CYCLE_DATE} -eq ${CYCLE_STR} ]]; then 
         mv diff_file.nc ${DIFF_ANAL_FILE}
      else
         mv diff_file.nc ${DIFF_FORE_FILE}
      fi
#
# NEXT DATE/TIME
      export CYCLE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${CYCLE_DATE} 1 2>/dev/null)
   done
#
# NEXT CYCLE      
   export DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} ${CYCLE_TIME} 2>/dev/null)
done
