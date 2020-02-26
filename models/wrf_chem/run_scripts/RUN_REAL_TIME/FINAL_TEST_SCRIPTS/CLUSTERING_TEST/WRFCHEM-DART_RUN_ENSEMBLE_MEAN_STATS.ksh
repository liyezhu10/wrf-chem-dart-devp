#!/bin/ksh -aeux
#
# Copyright 2019 University Corporation for Atmospheric Research and 
# Colorado Department of Public Health and Environment.
# 
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
# CONDITIONS OF ANY KIND, either express or implied. See the License for the 
# specific language governing permissions and limitations under the License.
# 
# Development of this code utilized the RMACC Summit supercomputer, which is 
# supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
# the University of Colorado Boulder, and Colorado State University. The Summit 
# supercomputer is a joint effort of the University of Colorado Boulder and 
# Colorado State University.
#
#SBATCH --job-name ens_stats
#SBATCH --output ens_stats.log
#SBATCH --account ${ACCOUNT}
#SBATCH --qos normal
#SBATCH --time 00:05:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --partition shas
#
#########################################################################
#
# CALCULATE ENSEMBLE MEAN_INPUT
#
#########################################################################
#
   if ${RUN_ENSEMBLE_MEAN_INPUT}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/ensemble_mean_input ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/ensemble_mean_input
         cd ${RUN_DIR}/${DATE}/ensemble_mean_input
      else
         cd ${RUN_DIR}/${DATE}/ensemble_mean_input
      fi
      rm -rf wrfinput_d${CR_DOMAIN}_mean
      rm -rf wrfbdy_d${CR_DOMAIN}_mean
      rm -rf wrfinput_d${FR_DOMAIN}_mean
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
         if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
            cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfinput_d${CR_DOMAIN}_${KMEM}
            cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfbdy_d${CR_DOMAIN}_${KMEM}
         else
            cp ${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfinput_d${CR_DOMAIN}_${KMEM}
            cp ${UPDATE_BC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfbdy_d${CR_DOMAIN}_${KMEM}
         fi
         let MEM=${MEM}+1
      done
#      cp ${REAL_DIR}/wrfinput_d${FR_DOMAIN}_${START_FILE_DATE} wrfinput_d${FR_DOMAIN}_mean
#
# Calculate ensemble mean
      ncea -n ${NUM_MEMBERS},4,1 wrfinput_d${CR_DOMAIN}_0001 wrfinput_d${CR_DOMAIN}_mean
      ncea -n ${NUM_MEMBERS},4,1 wrfbdy_d${CR_DOMAIN}_0001 wrfbdy_d${CR_DOMAIN}_mean
#
# Calculate ensemble spread
      rm -rf wrfinput_d${CR_DOMAIN}_tmp*
      rm -rf wrfinput_d${CR_DOMAIN}_sprd 
      ncecat -n ${NUM_MEMBERS},4,1 wrfinput_d${CR_DOMAIN}_0001 wrfinput_d${CR_DOMAIN}_tmp1
      ncwa -a record wrfinput_d${CR_DOMAIN}_tmp1 wrfinput_d${CR_DOMAIN}_tmp2
      ncbo --op_typ='-' wrfinput_d${CR_DOMAIN}_tmp1 wrfinput_d${CR_DOMAIN}_tmp2 wrfinput_d${CR_DOMAIN}_tmp3
      ncra -y rmssdn wrfinput_d${CR_DOMAIN}_tmp3 wrfinput_d${CR_DOMAIN}_sprd
      rm -rf wrfinput_d${CR_DOMAIN}_tmp*
      rm -rf wrfinput_d${CR_DOMAIN}_*0*
      rm -rf wrfbdy_d${CR_DOMAIN}_*0*
   fi
#
#########################################################################
#
# CALCULATE ENSEMBLE MEAN_OUTPUT
#
#########################################################################
#
   if ${RUN_ENSEMBLE_MEAN_OUTPUT}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/ensemble_mean_output ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/ensemble_mean_output
         cd ${RUN_DIR}/${DATE}/ensemble_mean_output
      else
         cd ${RUN_DIR}/${DATE}/ensemble_mean_output
      fi
      if [[ ${DATE} -eq ${INITIAL_DATE} ]]; then
         export OUTPUT_DIR=${WRFCHEM_INITIAL_DIR}
      else
         export OUTPUT_DIR=${WRFCHEM_CYCLE_CR_DIR}
      fi
      rm -rf wrfout_d${CR_DOMAIN}_*
      export P_DATE=${DATE}
      export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${FCST_PERIOD} 2>/dev/null)
      while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do
         export P_YYYY=$(echo $P_DATE | cut -c1-4)
         export P_MM=$(echo $P_DATE | cut -c5-6)
         export P_DD=$(echo $P_DATE | cut -c7-8)
         export P_HH=$(echo $P_DATE | cut -c9-10)
         export P_FILE_DATE=${P_YYYY}-${P_MM}-${P_DD}_${P_HH}:00:00
         let MEM=1
         while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
            export CMEM=e${MEM}
            export KMEM=${MEM}
            if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
            if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
            if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
            rm -rf wrfout_d${CR_DOMAIN}_${KMEM}
            ln -sf ${OUTPUT_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${P_FILE_DATE} wrfout_d${CR_DOMAIN}_${KMEM}
            let MEM=${MEM}+1
         done
#         cp ${OUTPUT_DIR}/run_e001/wrfout_d${CR_DOMAIN}_${P_FILE_DATE} wrfout_d${CR_DOMAIN}_${P_DATE}_mean
#
# Calculate ensemble mean
         ncea -n ${NUM_MEMBERS},4,1 wrfout_d${CR_DOMAIN}_0001 wrfout_d${CR_DOMAIN}_${P_DATE}_mean
         export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${HISTORY_INTERVAL_HR} 2>/dev/null)
      done
   fi
