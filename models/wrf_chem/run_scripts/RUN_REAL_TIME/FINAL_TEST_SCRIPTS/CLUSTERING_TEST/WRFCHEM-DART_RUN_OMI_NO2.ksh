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
# RUN OMI NO2 OBSERVATIONS
#
#########################################################################
#
   if ${RUN_OMI_NO2_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/omi_no2_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/omi_no2_obs
         cd ${RUN_DIR}/${DATE}/omi_no2_obs
      else
         cd ${RUN_DIR}/${DATE}/omi_no2_obs
      fi
#
      export FILE=obs_seq_comb_filtered_${DATE}.out
      cp ${EXPERIMENT_OMI_NO2_DIR}/${DATE}/${FILE} ./obs_seq_${DATE}.out
#
      if [[ -s obs_seq_${DATE}.out ]]; then
         cp obs_seq_${DATE}.out obs_seq_omi_no2_${DATE}.out
         rm obs_seq_${DATE}.out
      else
          touch NO_DATA_${DATE}
      fi
   fi
