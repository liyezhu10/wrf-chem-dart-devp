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
# MODIFIED VERSION OF /DART/models/WRF/regression/CONUS-V3/icbc_real.ksh
# TO SETUP AN ENVIRONMENT TO CONVERT OBSERVATIONS TO obs_seq.
#
# SET TIME INFORMATION
  export START_DATE=2014071700
  export END_DATE=2014080418
#  export END_DATE=2008060112
#  export END_DATE=2008063018
  export TIME_INC=6
  export ASIM_WINDOW=3
#
# SYSTEM SPECIFIC SETTINGS
  export PROCS=8
  export OB_TYPE=obs
#
# PATHS
  export WRFDA_VER=WRFDAv3.4_dmpar
  export WRF_VER=WRFv3.4_dmpar
  export DART_VER=DART_CHEM_MY_BRANCH_ALL
#
# INDEPENDENT DIRECTORIES
  export ROOT_DIR=/glade/p/work/mizzi
#  export DATA_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_dat_No_SVD
  export DATA_DIR=/glade/p/frappe/mizzi/REAL_TIME_DATA_DATA/obs_FRAPPE_CO_DnN_dat_No_SVD
  export ASIM_DIR=/glade/scratch/mizzi/MOPITT_to_OBSSEQ
#
# OUTPUT DIR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RAWR_NO_ROT_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RAWR_F50_NO_ROT_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_NO_ROT_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_ME_NO_ROT_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_QOR_NO_ROT_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_QOR_NO_SCALE_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_QOR_SCALE_SUPR
  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_CPSR_SCALE_SUPR_SINGLE_CLUSTER
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_CPSR_NO_SCALE_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_CPSR_SCALE_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_NO_ROT_RJ3_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_CPSR_SCALE_RJ3_SUPR
#  export OBS_MOPITT_OUT_DIR=/glade/p/acd/mizzi/AVE_TEST_DATA/obs_MOPITT_CO_RETR_CPSR_SCALE_RJ0_SUPR
#
# DEPENDENT DIRECTORIES
  export CODE_DIR=${ROOT_DIR}/TRUNK
  export HYBRID_DIR=${ROOT_DIR}/HYBRID_TRUNK
  export WRF_DIR=${CODE_DIR}/${WRF_VER}
  export VAR_DIR=${CODE_DIR}/${WRFDA_VER}
  export BUILD_DIR=${VAR_DIR}/var/build
  export DART_DIR=${CODE_DIR}/${DART_VER}
  export TOOL_DIR=${VAR_DIR}/var/da
  export HYBRID_SCRIPTS_DIR=${HYBRID_DIR}/hybrid_scripts
#
# MAKE ASSIMILATION DIRECTORY AND GO TO IT
  if [[ ! -d ${ASIM_DIR} ]]; then mkdir -p ${ASIM_DIR}; fi
  mkdir -p ${ASIM_DIR}
#
# BEGIN DAY AND TIME LOOP
  export L_DATE=${START_DATE}
  while [[ ${L_DATE} -le ${END_DATE} ]]; do
     export YYYY=$(echo $L_DATE | cut -c1-4)
     export YY=$(echo $L_DATE | cut -c3-4)
     export MM=$(echo $L_DATE | cut -c5-6)
     export DD=$(echo $L_DATE | cut -c7-8)
     export HH=$(echo $L_DATE | cut -c9-10)
     export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} -${ASIM_WINDOW} 2>/dev/null)  
     export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
     export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
     export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
     export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
#
# DART TIME INFO (NO LEADING ZEROS)
     export DT_YYYY=${YYYY}
     export DT_YY=$(echo $L_DATE | cut -c3-4)
     export DT_MM=${MM} 
     export DT_DD=${DD} 
     export DT_HH=${HH} 
     (( DT_MM = ${DT_MM} + 0 ))
     (( DT_DD = ${DT_DD} + 0 ))
     (( DT_HH = ${DT_HH} + 0 ))
#    
     export YEAR_INIT=${DT_YYYY}
     export MONTH_INIT=${DT_MM}
     export DAY_INIT=${DT_DD}
     export HOUR_INIT=${DT_HH}
     export YEAR_END=${DT_YYYY}
     export MONTH_END=${DT_MM}
     export DAY_END=${DT_DD}
     export HOUR_END=${DT_HH}
     export DA_TIME_WINDOW=0
#
# RUN_MOPITT_ASCII_TO_DART
        cd ${ASIM_DIR}
        if [[ ! -d ${ASIM_DIR}/ascii_to_dart/${YYYY}${MM} ]]; then mkdir -p ${ASIM_DIR}/ascii_to_dart/${YYYY}${MM}; fi
        cd ${ASIM_DIR}/ascii_to_dart/${YYYY}${MM}
        if [[ ${HH} -eq 0 ]] then
           export L_YYYY=${PAST_YYYY}
           export L_MM=${PAST_MM}
           export L_DD=${PAST_DD}
           export L_HH=24
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
           export F_DATE=${YYYY}${MM}${DD}${HH}
        else
           export L_YYYY=${YYYY}
           export L_MM=${MM}
           export L_DD=${DD}
           export L_HH=${HH}
           export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
           export F_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
        fi
        export NL_YEAR=${L_YYYY}
        export NL_MONTH=${L_MM}
        export NL_DAY=${L_DD}
        export NL_HOUR=${L_HH}
        if [[ ${L_HH} -eq 24 ]]; then
           NL_BIN_BEG=21.01
           NL_BIN_END=3.00
        elif [[ ${L_HH} -eq 6 ]]; then
           NL_BIN_BEG=3.01
           NL_BIN_END=9.00
        elif [[ ${L_HH} -eq 12 ]]; then
           NL_BIN_BEG=9.01
           NL_BIN_END=15.00
        elif [[ ${L_HH} -eq 18 ]]; then
           NL_BIN_BEG=15.01
           NL_BIN_END=21.00
        fi
        export NL_FILEDIR="'"${ASIM_DIR}/ascii_to_dart/${YYYY}${MM}/"'" 
        export NL_FILENAME="'"${D_DATE}.dat"'" 
#
# USE MOPITT DATA 
        if [[ -e input.nml ]]; then rm -rf input.nml; fi
        ${HYBRID_SCRIPTS_DIR}/da_create_dart_mopitt_input_nml.ksh
#
# GET INTERMEDIATE ASCII DATA
        cp ${DATA_DIR}/${D_DATE}.dat ./.
#
# GET EXECUTABLE
        cp ${DART_DIR}/observations/MOPITT_CO/work/mopitt_ascii_to_obs ./.
        ./mopitt_ascii_to_obs
#
# COPY OUTPUT TO ARCHIVE LOCATION
        export MOPITT_FILE=mopitt_obs_seq${D_DATE}
        if [[ -e ${MOPITT_FILE} ]]; then
           cp ${MOPITT_FILE} ${OBS_MOPITT_OUT_DIR}/obs_seq_mopitt_${D_DATE}
#           mkdir -p ${OBS_MOPITT_OUT_DIR}/${F_DATE}
#           cp ${MOPITT_FILE} ${OBS_MOPITT_OUT_DIR}/${F_DATE}/obs_seq_mopitt_${F_DATE}
        else
           touch ${OBS_MOPITT_OUT_DIR}/NO_DATA_${D_DATE}
#           mkdir -p ${OBS_MOPITT_OUT_DIR}/${F_DATE}
#           touch ${OBS_MOPITT_OUT_DIR}/${F_DATE}/NO_DATA_${F_DATE}
        fi
#
# LOOP TO NEXT DAY AND TIME 
     export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${TIME_INC} 2>/dev/null)  
  done 
exit
