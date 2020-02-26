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
# RUN PANDA PM25 OBSERVATIONS
#
#########################################################################
#
   if ${RUN_PANDA_PM25_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/panda_pm25_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/panda_pm25_obs
         cd ${RUN_DIR}/${DATE}/panda_pm25_obs
      else
         cd ${RUN_DIR}/${DATE}/panda_pm25_obs
      fi
#
# GET PANDA DATA
      if [[ ! -e panda_station_coordinates.csv  ]]; then
         cp ${EXPERIMENT_PANDA_DIR}/panda_station_coordinates.csv ./.
      fi
      if [[ ! -e panda_stationData.csv  ]]; then
         cp ${EXPERIMENT_PANDA_DIR}/panda_stationData.csv ./.
      fi
#
      export ASIM_MIN_MN=0
      export ASIM_MIN_SS=0
      export ASIM_MAX_MN=0
      export ASIM_MAX_SS=0
#
# RUN_PANDA_PM25_ASCII_TO_DART
      if [[ ${HH} -eq 0 ]]; then
         export L_YYYY=${ASIM_MIN_YYYY}
         export L_MM=${ASIM_MIN_MM}
         export L_DD=${ASIM_MIN_DD}
         export L_HH=24
         export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
      else
         export L_YYYY=${YYYY}
         export L_MM=${MM}
         export L_DD=${DD}
         export L_HH=${HH}
         export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
      fi
      export NL_YEAR=${L_YYYY}
      export NL_MONTH=${L_MM}
      export NL_DAY=${L_DD}
      export NL_HOUR=${L_HH}
#
      export NL_FILENAME_COORD=\'panda_station_coordinates.csv\'
      export NL_FILENAME_DATA=\'panda_stationData.csv\'
      export NL_LAT_MN=${NNL_MIN_LAT}
      export NL_LAT_MX=${NNL_MAX_LAT}
      export NL_LON_MN=${NNL_MIN_LON}
      export NL_LON_MX=${NNL_MAX_LON}
#
# GET EXECUTABLE
      cp ${DART_DIR}/observations/PANDA/work/panda_pm25_ascii_to_obs ./.
      rm -rf create_panda_obs_nml.nl
      rm -rf input.nml
      ${HYBRID_SCRIPTS_DIR}/da_create_dart_panda_input_nml.ksh
      ./panda_pm25_ascii_to_obs
#
# COPY OUTPUT TO ARCHIVE LOCATION
      export PANDA_OUT_FILE=panda_obs_seq
      export PANDA_ARCH_FILE=obs_seq_panda_pm25_${DATE}.out
      if [[ -s ${PANDA_OUT_FILE} ]]; then
         cp ${PANDA_OUT_FILE} ${PANDA_ARCH_FILE}
         rm ${PANDA_OUT_FILE}
      else
         touch NO_DATA_${D_DATE}
      fi     
   fi
