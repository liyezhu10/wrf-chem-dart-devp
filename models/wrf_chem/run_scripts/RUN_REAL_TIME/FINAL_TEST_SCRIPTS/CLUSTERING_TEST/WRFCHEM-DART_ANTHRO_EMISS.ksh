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
# RUN WRFCHEM_CHEMI
#
#########################################################################
#
   if ${RUN_WRFCHEM_CHEMI}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chemi ]]; then
         mkdir ${RUN_DIR}/${DATE}/wrfchem_chemi
         cd ${RUN_DIR}/${DATE}/wrfchem_chemi
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_chemi
      fi
      export L_DATE=${DATE}
      export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
#
      while [[ ${L_DATE} -le ${LE_DATE} ]]; do
         export L_YYYY=$(echo $L_DATE | cut -c1-4)
         export L_MM=$(echo $L_DATE | cut -c5-6)
         export L_DD=$(echo $L_DATE | cut -c7-8)
         export L_HH=$(echo $L_DATE | cut -c9-10)
#
         export FILE_PATH=${EXPERIMENT_WRFCHEMI_DIR}
         cp ${FILE_PATH}/wrfchemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 ./.
         cp ${FILE_PATH}/wrfchemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 ./.
         chmod a+rwx wrfchemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 
         chmod a+rwx wrfchemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00 
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
      done
   fi
