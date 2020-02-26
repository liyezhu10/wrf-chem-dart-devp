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
# RUN IASI CO OBSERVATIONS
#
#########################################################################
#
   if ${RUN_IASI_CO_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/iasi_co_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/iasi_co_obs
         cd ${RUN_DIR}/${DATE}/iasi_co_obs
      else
         cd ${RUN_DIR}/${DATE}/iasi_co_obs
      fi
#
# set file prefix for IASI
# this depends on versions and file times (edit if necessary)
      export FILE_PRE='VERSION2_NCAR_IASI_xxx_1C_M02'
#
# set file suffix for IASI
# this depends on versions and file times (edit if necessary)
      export FILE_EXT='hdf'
#
      export L_PAST_DATE=$($BUILD_DIR/da_advance_time.exe $DATE -24 2>/dev/null)
      export L_PAST_YY=$(echo $L_PAST_DATE | cut -c1-4)
      export L_PAST_MM=$(echo $L_PAST_DATE | cut -c5-6)
      export L_PAST_DD=$(echo $L_PAST_DATE | cut -c7-8)
      export L_PAST_HH=$(echo $L_PAST_DATE | cut -c9-10)
#
      let TEMP_MIN_HH=${ASIM_MIN_HH}
      let TEMP_MAX_HH=${ASIM_MAX_HH}
      (( BIN_BEG_SEC=${TEMP_MIN_HH}*60*60+1 ))
      (( BIN_END_SEC=${TEMP_MAX_HH}*60*60 ))
#
      export NCNT=3
#
# Test for special case
      if [[ ! ${HH} == 00 ]]; then
#
# Normal cases (06Z, 12Z, 18Z)
         export A_DATE_START=$($BUILD_DIR/da_advance_time.exe ${ASIM_MIN_DATE} -${NCNT} 2>/dev/null)
         export A_DATE=${A_DATE_START}
         while [[ ${A_DATE} -le ${ASIM_MAX_DATE} ]]; do 
            if [[ ${A_DATE} == ${A_DATE_START} ]]; then
               export ASIM_OUTFILE=${YYYY}${MM}${DD}${HH}.dat
               rm -rf ${ASIM_OUTFILE}
               touch ${ASIM_OUTFILE}
            fi
            export A_YY=$(echo $A_DATE | cut -c1-4)
            export A_MM=$(echo $A_DATE | cut -c5-6)
            export A_DD=$(echo $A_DATE | cut -c7-8)
            export A_HH=$(echo $A_DATE | cut -c9-10)
            export ICNT=0
            while [[ ${ICNT} -le ${NCNT} ]]; do
               export TEST=$($BUILD_DIR/da_advance_time.exe ${A_DATE} ${ICNT} 2>/dev/null)
               export ND_YY=$(echo $TEST | cut -c1-4)
               export ND_MM=$(echo $TEST | cut -c5-6)
               export ND_DD=$(echo $TEST | cut -c7-8)
               export ND_HH=$(echo $TEST | cut -c9-10)
               export FILE=`ls ${EXPERIMENT_IASI_CO_DIR}/${A_YY}/${A_MM}/${A_DD}/${FILE_PRE}_${A_YY}${A_MM}${A_DD}${A_HH}*Z_${ND_YY}${ND_MM}${ND_DD}${ND_HH}*Z_*`
               if [[ -e ${FILE} ]]; then 
                  export OUTFILE_NM=TEMP_FILE.dat
                  export INFILE=\'${FILE}\'
                  export OUTFILE=\'${OUTFILE_NM}\'
#
# echo what we are processing at the moment
#                  echo ${INFILE}
#                  echo ${OUTFILE}
#                  echo ${BIN_BEG_SEC}
#                  echo ${BIN_END_SEC}
#
# this is the call to an IDL routine to read and write variables
# if already processed (with output), then this can be skipped (do_iasi=0)
# else this needs to be called
                  export FILE=iasi_extract_no_transform_UA.pro
                  rm -rf ${FILE}
                  cp ${DART_DIR}/observations/IASI_CO/native_to_ascii/${FILE} ./.
                  idl << EOF
.compile iasi_extract_no_transform_UA.pro
iasi_extract_no_transform_UA,${INFILE},${OUTFILE},${BIN_BEG_SEC},${BIN_END_SEC}, ${NL_MIN_LON}, ${NL_MAX_LON}, ${NL_MIN_LAT}, ${NL_MAX_LAT}
EOF
#
# cat the output file to the assimlation window file
                  export ASIM_OUTFILE=${YYYY}${MM}${DD}${HH}.dat
                  if [[ -e ${OUTFILE_NM} ]]; then
                     cat ${OUTFILE_NM} >> ${ASIM_OUTFILE}
                     rm -rf ${OUTFILE_NM}
                  fi
               fi
               (( ICNT=${ICNT}+1 ))
            done
#
# go to next hour
            export AA_DATE=${A_DATE}
            export A_DATE=$(${BUILD_DIR}/da_advance_time.exe ${AA_DATE} 1 2>/dev/null)
         done
      else   
#
# Special case (00Z)
         let TEMP_MIN_HH=${ASIM_MIN_HH}
         let TEMP_MAX_HH=${ASIM_MAX_HH}
         (( BIN_BEG_SEC=${TEMP_MIN_HH}*60*60+1 ))
         (( BIN_END_SEC=${TEMP_MAX_HH}*60*60 ))
         export A_DATE_START=$($BUILD_DIR/da_advance_time.exe ${ASIM_MIN_DATE} -${NCNT} 2>/dev/null)
         export A_DATE=${A_DATE_START}
         while [[ ${A_DATE} -le ${ASIM_MAX_DATE} ]]; do 
            if [[ ${A_DATE} == ${A_DATE_START} ]]; then
               export ASIM_OUTFILE=${YYYY}${MM}${DD}${HH}.dat
               rm -rf ${ASIM_OUTFILE}
               touch ${ASIM_OUTFILE}
            fi
            export A_YY=$(echo $A_DATE | cut -c1-4)
            export A_MM=$(echo $A_DATE | cut -c5-6)
            export A_DD=$(echo $A_DATE | cut -c7-8)
            export A_HH=$(echo $A_DATE | cut -c9-10)
#
            if [[ ${PAST_YY} == ${A_YY} && ${PAST_MM} == ${A_MM} && ${PAST_DD} == ${A_DD} ]]; then
               (( BIN_BEG_SEC=${TEMP_MIN_HH}*60*60+1 ))
               (( BIN_END_SEC=24*60*60 ))
            else
               (( BIN_BEG_SEC=1 ))
               (( BIN_END_SEC=${TEMP_MAX_HH}*60*60 ))
            fi 
            export ICNT=0
            while [[ ${ICNT} -le ${NCNT} ]]; do
               export TEST=$($BUILD_DIR/da_advance_time.exe ${A_DATE} ${ICNT} 2>/dev/null)
               export ND_YY=$(echo $TEST | cut -c1-4)
               export ND_MM=$(echo $TEST | cut -c5-6)
               export ND_DD=$(echo $TEST | cut -c7-8)
               export ND_HH=$(echo $TEST | cut -c9-10)
               export FILE=`ls ${EXPERIMENT_IASI_CO_DIR}/${A_YY}/${A_MM}/${A_DD}/${FILE_PRE}_${A_YY}${A_MM}${A_DD}${A_HH}*Z_${ND_YY}${ND_MM}${ND_DD}${ND_HH}*Z_*`
               if [[ -e ${FILE} ]]; then 
                  export OUTFILE_NM=TEMP_FILE.dat
                  export INFILE=\'${FILE}\'
                  export OUTFILE=\'${OUTFILE_NM}\'
#
# echo what we are processing at the moment
#                  echo ${INFILE}
#                  echo ${OUTFILE}
#                  echo ${BIN_BEG_SEC}
#                  echo ${BIN_END_SEC}
#
# this is the call to an IDL routine to read and write variables
# if already processed (with output), then this can be skipped (do_iasi=0)
# else this needs to be called
                  export FILE=iasi_extract_no_transform_UA.pro
                  rm -rf ${FILE}
                  cp ${DART_DIR}/observations/IASI_CO/native_to_ascii/${FILE} ./.
                  idl << EOF
.compile iasi_extract_no_transform_UA.pro
iasi_extract_no_transform_UA,${INFILE},${OUTFILE},${BIN_BEG_SEC},${BIN_END_SEC}, ${NNL_MIN_LON}, ${NNL_MAX_LON}, ${NNL_MIN_LAT}, ${NNL_MAX_LAT}
EOF
#
# cat the output file to the assimlation window file
                  export ASIM_OUTFILE=${YYYY}${MM}${DD}${HH}.dat
                  if [[ -e ${OUTFILE_NM} ]]; then
                     cat ${OUTFILE_NM} >> ${ASIM_OUTFILE}
                     rm -rf ${OUTFILE_NM}
                  fi
               fi
               (( ICNT=${ICNT}+1 ))
            done
#
# go to next hour
            export AA_DATE=${A_DATE}
            export A_DATE=$(${BUILD_DIR}/da_advance_time.exe ${AA_DATE} 1 2>/dev/null)
         done
      fi
#
# convert to obseq file
      export L_PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${ASIM_WINDOW} 2>/dev/null)  
      export L_PAST_YYYY=$(echo $L_PAST_DATE | cut -c1-4)
      export L_PAST_MM=$(echo $L_PAST_DATE | cut -c5-6)
      export L_PAST_DD=$(echo $L_PAST_DATE | cut -c7-8)
      export L_PAST_HH=$(echo $L_PAST_DATE | cut -c9-10)
#
# DART TIME INFO (NO LEADING ZEROS)
      export DT_YYYY=${YYYY}
      export DT_YY=$(echo $DATE | cut -c3-4)
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
# RUN_IASI_ASCII_TO_DART
      if [[ ${HH} -eq 0 ]]; then
         export L_YYYY=${L_PAST_YYYY}
         export L_MM=${L_PAST_MM}
         export L_DD=${L_PAST_DD}
         export L_HH=24
         export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
         export DD_DATE=${YYYY}${MM}${DD}${HH}
      else
         export L_YYYY=${YYYY}
         export L_MM=${MM}
         export L_DD=${DD}
         export L_HH=${HH}
         export D_DATE=${L_YYYY}${L_MM}${L_DD}${L_HH}
         export DD_DATE=${D_DATE}
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
      export NL_FILEDIR=\'./\' 
      export NL_FILENAME=${D_DATE}.dat
      export NL_IASI_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_IASI}\'
      export NL_IASI_O3_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_IASI}\'
      export NL_FAC_OBS_ERROR=${NL_FAC_OBS_ERROR_IASI}
      export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
      export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
#
# USE IASI DATA 
      rm -rf input.nml
      ${HYBRID_SCRIPTS_DIR}/da_create_dart_iasi_input_nml.ksh
#
# GET INTERMEDIATE ASCII DATA
      if [[ ! -e ${D_DATE}.dat ]]; then cp ${DD_DATE}.dat ./${D_DATE}.dat; fi
#
# GET EXECUTABLE
      cp ${DART_DIR}/observations/IASI_CO/work/iasi_ascii_to_obs ./.
      ./iasi_ascii_to_obs > index.html 2>&1
#
# COPY OUTPUT TO ARCHIVE LOCATION
      export IASI_FILE=iasi_obs_seq${D_DATE}
      if [[ -s ${IASI_FILE} ]]; then
         cp ${IASI_FILE} obs_seq_iasi_co_${DATE}.out
      else
         touch NO_DATA_${D_DATE}
      fi
   fi
