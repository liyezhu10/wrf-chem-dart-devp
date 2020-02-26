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
# RUN MODIS AOD OBSERVATIONS
#
#########################################################################
#
   if ${RUN_MODIS_AOD_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/modis_aod_obs ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/modis_aod_obs
         cd ${RUN_DIR}/${DATE}/modis_aod_obs
      else
         cd ${RUN_DIR}/${DATE}/modis_aod_obs
      fi
#
# set file prefix for MODIS
# this depends on versions and file times (edit if necessary)
      export FILE_PRE='MYD04_L2.A'
#
# set file suffix for MODIS
# this depends on versions and file times (edit if necessary)
      export FILE_EXT='hdf'
#
      export MODIS_INDIR=${EXPERIMENT_MODIS_AOD_DIR}
      export OUTFILE=modis_aod_ascii_${YYYY}${MM}${DD}${HH}
      (( N_YYYY=${YYYY}+0 ))
      (( N_MM=${MM}+0 ))
      (( N_DD=${DD}+0 ))
      (( N_HH=${HH}+0 ))
      (( N_ASIM_WIN=${ASIM_WINDOW}+0 ))
#
      export FILE=modis_extract_hdf.pro
      rm -rf ${FILE}
      cp ${DART_DIR}/observations/MODIS/native_to_ascii/${FILE} ./.
      idl << EOF
.compile modis_extract_hdf.pro
modis_extract_hdf, "${MODIS_INDIR}", "${OUTFILE}", ${N_YYYY}, ${N_MM}, ${N_DD}, ${N_HH}, ${N_ASIM_WIN}, ${NL_MIN_LON}, ${NL_MAX_LON}, ${NL_MIN_LAT}, ${NL_MAX_LAT}
EOF
#
# convert ASCII to obs_seq file
#      rm -rf input.nml
#      rm -rf modis_asciidata.input
#      rm -rf modis_obs_seq.out
      if [[ -s modis_aod_ascii_${YYYY}${MM}${DD}${HH} ]]; then
         cp modis_aod_ascii_${YYYY}${MM}${DD}${HH} modis_asciidata.input
         cp ${DART_DIR}/observations/MODIS/work/input.nml ./.
         ${DART_DIR}/observations/MODIS/work/modis_ascii_to_obs
         if [[ -s modis_obs_seq.out ]]; then
            export MODIS_FILE=obs_seq_modis_aod_${DATE}.out
            mv modis_obs_seq.out ${MODIS_FILE}
         else
            touch NO_MODIS_AOD_${DATE}
         fi
      else
         touch NO_MODIS_AOD_${DATE}
      fi
   fi
