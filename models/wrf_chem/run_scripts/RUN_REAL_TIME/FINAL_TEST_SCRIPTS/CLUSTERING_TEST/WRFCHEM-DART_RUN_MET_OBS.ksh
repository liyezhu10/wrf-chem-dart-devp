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

#########################################################################
#
# RUN PREPBUFR MET OBSERVATIONS
#
#########################################################################
#
# APM: This block needs to be revised so we can convert a single prepbufr
#      file in real time we can use only the obs that are on the current
#      prepbufr file.
#
   if ${RUN_MET_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/prepbufr_met_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/prepbufr_met_obs
         cd ${RUN_DIR}/${DATE}/prepbufr_met_obs
      else
         cd ${RUN_DIR}/${DATE}/prepbufr_met_obs
      fi
#
# GET PREPBUFR FILES
#           
      export L_DATE=${D_YYYY}${D_MM}${D_DD}06
      export E_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +24 2>/dev/null)
      while [[ ${L_DATE} -le ${E_DATE} ]]; do
         export L_YYYY=$(echo $L_DATE | cut -c1-4)
         export L_YY=$(echo $L_DATE | cut -c3-4)
         export L_MM=$(echo $L_DATE | cut -c5-6)
         export L_DD=$(echo $L_DATE | cut -c7-8)
         export L_HH=$(echo $L_DATE | cut -c9-10)
         cp ${EXPERIMENT_PREPBUFR_DIR}/${L_YYYY}${L_MM}${L_DD}${L_HH}/prepbufr.gdas.${L_YYYY}${L_MM}${L_DD}${L_HH}.wo40.be prepqm${L_YY}${L_MM}${L_DD}${L_HH}
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +6 2>/dev/null)
      done
#
# GET DART input.nml
      rm -rf input.nml
      cp ${DART_DIR}/observations/NCEP/prep_bufr/work/input.nml ./.
#
# RUN_PREPBUFR TO ASCII CONVERTER
      ${DART_DIR}/observations/NCEP/prep_bufr/work/prepbufr_RT.csh ${D_YYYY} ${DD_MM} ${DD_DD} ${DD_DD} ${DART_DIR}/observations/NCEP/prep_bufr/exe > index.file 2>&1
#
# RUN ASCII TO OBS_SEQ CONVERTER
      ${HYBRID_SCRIPTS_DIR}/da_create_dart_ncep_ascii_to_obs_input_nml_RT.ksh
      ${DART_DIR}/observations/NCEP/ascii_to_obs/work/create_real_obs > index_create 2>&1
#
      mv obs_seq${D_DATE} obs_seq_prep_${DATE}.out
   fi
