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
#########################################################################
#
# RUN MOPITT CO OBSERVATIONS
#
#########################################################################
#
   if ${RUN_MOPITT_CO_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/mopitt_co_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/mopitt_co_obs
         cd ${RUN_DIR}/${DATE}/mopitt_co_obs
      else
         cd ${RUN_DIR}/${DATE}/mopitt_co_obs
      fi
#
# SET MOPITT PARAMETERS
      export MOPITT_FILE_PRE=MOP02J-
      export MOPITT_FILE_EXT=-L2V10.1.3.beta.hdf   
      export MOP_OUTFILE=\'MOPITT_CO_${D_DATE}'.dat'\'
      rm -rf ${MOP_OUTFILE}
#
#  SET OBS WINDOW
      export BIN_BEG=${ASIM_MN_HH}
      export BIN_END=${ASIM_MX_HH}
      export FLG=0
#
# SET MOPITT INPUT DATA DIR
      if [[ ${BIN_END} -ne 3 ]]; then
         export MOP_INFILE=\'${EXPERIMENT_MOPITT_CO_DIR}/${MOPITT_FILE_PRE}${ASIM_MX_YYYY}${ASIM_MX_MM}${ASIM_MX_DD}${MOPITT_FILE_EXT}\'
      else
         export FLG=1
         export BIN_END=24
         export MOP_INFILE=\'${EXPERIMENT_MOPITT_CO_DIR}/${MOPITT_FILE_PRE}${ASIM_MN_YYYY}${ASIM_MN_MM}${ASIM_MN_DD}${MOPITT_FILE_EXT}\'
      fi
#
# COPY EXECUTABLE
      export FILE=mopitt_extract_no_transform_RT.pro
      rm -rf ${FILE}
      cp ${DART_DIR}/observations/MOPITT_CO/native_to_ascii/${FILE} ./.
      idl << EOF
.compile mopitt_extract_no_transform_RT.pro
mopitt_extract_no_transform_RT, ${MOP_INFILE}, ${MOP_OUTFILE}, ${BIN_BEG}, ${BIN_END}, ${NL_MIN_LON}, ${NL_MAX_LON}, ${NL_MIN_LAT}, ${NL_MAX_LAT}
EOF
#
# GET ADDITIONAL DATA FOR DAY-TO-DAY CROSSOVER
      if [[ ${FLG} -eq 1 ]];  then
         export FLG=0
         export BIN_BEG=0
         export BIN_END=3
         export MOP_INFILE=\'${EXPERIMENT_MOPITT_CO_DIR}/${MOPITT_FILE_PRE}${ASIM_MX_YYYY}${ASIM_MX_MM}${ASIM_MX_DD}${MOPITT_FILE_EXT}\'
#
         idl << EOF
.compile mopitt_extract_no_transform_RT.pro
mopitt_extract_no_transform_RT, ${MOP_INFILE}, ${MOP_OUTFILE}, ${BIN_BEG}, ${BIN_END}, ${NL_MIN_LON}, ${NL_MAX_LON}, ${NL_MIN_LAT}, ${NL_MAX_LAT}
EOF
      fi   
#
# SET NAMELIST TO CONVERT MOPITT ASCII TO OBS_SEQ 
      export NL_YEAR=${D_YYYY}
      export NL_MONTH=${D_MM}
      export NL_DAY=${D_DD}
      export NL_HOUR=${D_HH}
      if [[ ${D_HH} -eq 24 ]]; then
         export NL_BIN_BEG=21.01
         export NL_BIN_END=3.00
      elif [[ ${D_HH} -eq 6 ]]; then
         export NL_BIN_BEG=3.01
         export NL_BIN_END=9.00
      elif [[ ${D_HH} -eq 12 ]]; then
         export NL_BIN_BEG=9.01
         export NL_BIN_END=15.00
      elif [[ ${D_HH} -eq 18 ]]; then
         export NL_BIN_BEG=15.01
         export NL_BIN_END=21.00
      fi
      cp MOPITT_CO_${D_DATE}.dat ${D_DATE}.dat
      export NL_FILEDIR=\'./\' 
      export NL_FILENAME=${D_DATE}.dat
      export NL_MOPITT_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_MOPITT}\'
      export NL_FAC_OBS_ERROR=${NL_FAC_OBS_ERROR_MOPITT}
      export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
      export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
#
# USE MOPITT DATA 
      rm -rf input.nml
      ${HYBRID_SCRIPTS_DIR}/da_create_dart_mopitt_input_nml.ksh
#
# GET EXECUTABLE
      cp ${DART_DIR}/observations/MOPITT_CO/work/mopitt_ascii_to_obs ./.
      ./mopitt_ascii_to_obs > index.html 2>&1
#
# COPY OUTPUT TO ARCHIVE LOCATION
      export MOPITT_FILE=mopitt_obs_seq${D_DATE}
      touch obs_seq_mopitt_co_${DATE}.out
      if [[ -s ${MOPITT_FILE} ]]; then
         cp ${MOPITT_FILE} obs_seq_mopitt_co_${DATE}.out
      else
         touch NO_MOPITT_CO_${DATE}
      fi
   fi
