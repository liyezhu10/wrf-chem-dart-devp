#!/bin/ksh -x
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
# RUN GEOGRID
#
#########################################################################
#
   if [[ ${RUN_GEOGRID} = "true" ]]; then
      mkdir -p ${RUN_DIR}/geogrid
      cd ${RUN_DIR}/geogrid
#
      cp ${WPS_DIR}/geogrid.exe ./.
      export NL_DX=${DX_CR}
      export NL_DY=${DX_CR}
      export NL_START_DATE=${FILE_DATE}
      export NL_END_DATE=${NEXT_FILE_DATE}
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
      RANDOM=$$
      export TYPE=SERIAL
      export EXECUTE=geogrid.exe
      if [[ ${TYPE} == PARALLEL ]]; then
         mpirun -np ${TASKS} ./${EXECUTE} > index_${RANDOM}.html 2>&1
      elif [[ ${TYPE} == SERIAL ]]; then
        ./${EXECUTE} > index_${RANDOM}.html 2>&1
      fi
   fi
#
#########################################################################
#
# RUN UNGRIB
#
#########################################################################
#
   if [[ ${RUN_UNGRIB} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/ungrib
      cd ${RUN_DIR}/${DATE}/ungrib
      rm -rf GRIBFILE.*
#
      cp ${VTABLE_DIR}/Vtable.${VTABLE_TYPE} Vtable
      cp ${WPS_DIR}/ungrib.exe ./.
#
      export L_FCST_RANGE=${LBC_END}
      export L_START_DATE=${DATE}
      export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
      export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
      export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
      export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
      export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
      export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
      export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
      export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
      export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
      export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
      export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
      export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
      export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
      export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
      export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
      export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
      export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
      export NL_START_DATE=\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\',\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\'
      export NL_END_DATE=\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\',\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\'
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
# UNTAR THE PARENT FORECAST FILES
      FILES=''
      if [[ -e ${EXPERIMENT_GFS_DIR}/${DATE} ]]; then
         if [[ -e ${EXPERIMENT_GFS_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
#            cd ${EXPERIMENT_GFS_DIR}/${DATE}
            tar xvfs ${EXPERIMENT_GFS_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2}
#            cd ${RUN_DIR}/${DATE}/ungrib
         else
            echo 'APM: ERROR - No GRIB files in directory'
            exit
         fi
         sleep 30
#  
         if [[ ${SINGLE_FILE} == false ]]; then
            export CCHH=${HH}00
            (( LBC_ITR=${LBC_START} ))
            while [[ ${LBC_ITR} -le ${LBC_END} ]]; do
               if [[ ${LBC_ITR} -lt 1000 ]]; then export CFTM=${LBC_ITR}; fi
               if [[ ${LBC_ITR} -lt 100  ]]; then export CFTM=0${LBC_ITR}; fi
               if [[ ${LBC_ITR} -lt 10   ]]; then export CFTM=00${LBC_ITR}; fi
               if [[ ${LBC_ITR} -eq 0    ]]; then export CFTM=000; fi
#               export FILE=${EXPERIMENT_GFS_DIR}/${DATE}/${GRIB_PART1}${START_YEAR}${START_MONTH}${START_DAY}_${CCHH}_${CFTM}.grb2
               export FILE=${GRIB_PART1}${START_YEAR}${START_MONTH}${START_DAY}_${CCHH}_${CFTM}.grb2
               FILES="${FILES} ${FILE}"
               (( LBC_ITR=${LBC_ITR}+${LBC_FREQ} ))
            done
         else
            export FILE=${EXPERIMENT_GFS_DIR}/${DATE}/GFS_Global_0p5deg_20080612_1800.grib2
            FILES="${FILES} ${FILE}"
         fi
      fi
#
# LINK GRIB FILES
      ${WPS_DIR}/link_grib.csh $FILES
#
      RANDOM=$$
      export TYPE=SERIAL
      export EXECUTE=ungrib.exe
      if [[ ${TYPE} == PARALLEL ]]; then
         mpirun -np ${TASKS} ./${EXECUTE} > index_${RANDOM}.html 2>&1
      elif [[ ${TYPE} == SERIAL ]]; then
         ./${EXECUTE} > index_${RANDOM}.html 2>&1
      fi
#
# TAR THE PARENT FORECAST FILES
#      rm -rf *.grb2
#      if [[ -e ${EXPERIMENT_GFS_DIR}/${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
#         rm -rf ${EXPERIMENT_GFS_DIR}/${DATE}/${GRIB_PART1}*.grb2
#      else
#         cd ${EXPERIMENT_GFS_DIR}
#         tar -cf ${GRIB_PART1}${DATE}${GRIB_PART2} ${DATE}
#         mv ${GRIB_PART1}${DATE}${GRIB_PART2} ${DATE}/.
#         if [[ -e ${DATE}/${GRIB_PART1}${DATE}${GRIB_PART2} ]]; then
#            rm -rf ${DATE}/${GRIB_PART1}*.grb2
#         else
#            echo 'APM: Failed to created tar file'
#            exit
#         fi
#         cd ${RUN_DIR}/${DATE}/ungrib
#      fi
   fi
#
#########################################################################
#
# RUN METGRID
#
#########################################################################
#
   if [[ ${RUN_METGRID} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/metgrid
      cd ${RUN_DIR}/${DATE}/metgrid
#
      ln -fs ${GEOGRID_DIR}/geo_em.d${CR_DOMAIN}.nc ./.
      ln -fs ${GEOGRID_DIR}/geo_em.d${FR_DOMAIN}.nc ./.
      ln -fs ../ungrib/FILE:* ./.
      ln -fs ${WPS_DIR}/metgrid/METGRID.TBL.${METGRID_TABLE_TYPE} METGRID.TBL
      ln -fs ${WPS_DIR}/metgrid.exe .
#
      export L_FCST_RANGE=${LBC_END}
      export L_START_DATE=${DATE}
      export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
      export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
      export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
      export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
      export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
      export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
      export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
      export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
      export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
      export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
      export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
      export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
      export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
      export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
      export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
      export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
      export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
      export NL_START_DATE=\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\',\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\'
      export NL_END_DATE=\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\',\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\'
      ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
      RANDOM=$$
      export TYPE=SERIAL
      export EXECUTE=metgrid.exe
      if [[ ${TYPE} == PARALLEL ]]; then
         mpirun -np ${TASKS} ./${EXECUTE} > index_${RANDOM}.html 2>&1
      elif [[ ${TYPE} == SERIAL ]]; then
         ./${EXECUTE} > index_${RANDOM}.html 2>&1
      fi
   fi
#
#########################################################################
#
# RUN REAL
#
#########################################################################
#
   if [[ ${RUN_REAL} = "true" ]]; then 
      mkdir -p ${RUN_DIR}/${DATE}/real
      cd ${RUN_DIR}/${DATE}/real
#
      cp ${WRF_DIR}/main/real.exe ./.
      cp ${EXPERIMENT_HIST_IO_DIR}/hist_io_flds_v1 ./.
      cp ${EXPERIMENT_HIST_IO_DIR}/hist_io_flds_v2 ./.
#
# LINK IN THE METGRID FILES
      export P_DATE=${DATE}
      export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_END} 2>/dev/null)
      while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do
         export P_YYYY=$(echo $P_DATE | cut -c1-4)
         export P_MM=$(echo $P_DATE | cut -c5-6)
         export P_DD=$(echo $P_DATE | cut -c7-8)
         export P_HH=$(echo $P_DATE | cut -c9-10)
         export P_FILE_DATE=${P_YYYY}-${P_MM}-${P_DD}_${P_HH}:00:00.nc
         ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d${CR_DOMAIN}.${P_FILE_DATE} ./.
         ln -sf ${RUN_DIR}/${DATE}/metgrid/met_em.d${FR_DOMAIN}.${P_FILE_DATE} ./.
         export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_FREQ} 2>/dev/null) 
      done
#
# LOOP THROUGH BDY TENDENCY TIMES FOR PERTURB_BC
      export P_DATE=${DATE}
      export P_END_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${FCST_PERIOD} 2>/dev/null)
      while [[ ${P_DATE} -le ${P_END_DATE} ]] ; do      
#
# CREATE WRF NAMELIST
         export NL_IOFIELDS_FILENAME=' '
         export NL_IOFIELDS_FILENAME=\'hist_io_flds_v1\',\'hist_io_flds_v2\'
         export L_FCST_RANGE=${FCST_PERIOD}
         export NL_DX=${DX_CR},${DX_FR}
         export NL_DY=${DX_CR},${DX_FR}
         export L_START_DATE=${P_DATE}
         export L_END_DATE=$($BUILD_DIR/da_advance_time.exe ${L_START_DATE} ${L_FCST_RANGE} 2>/dev/null)
         export L_START_YEAR=$(echo $L_START_DATE | cut -c1-4)
         export L_START_MONTH=$(echo $L_START_DATE | cut -c5-6)
         export L_START_DAY=$(echo $L_START_DATE | cut -c7-8)
         export L_START_HOUR=$(echo $L_START_DATE | cut -c9-10)
         export L_END_YEAR=$(echo $L_END_DATE | cut -c1-4)
         export L_END_MONTH=$(echo $L_END_DATE | cut -c5-6)
         export L_END_DAY=$(echo $L_END_DATE | cut -c7-8)
         export L_END_HOUR=$(echo $L_END_DATE | cut -c9-10)
         export NL_START_YEAR=$(echo $L_START_DATE | cut -c1-4),$(echo $L_START_DATE | cut -c1-4)
         export NL_START_MONTH=$(echo $L_START_DATE | cut -c5-6),$(echo $L_START_DATE | cut -c5-6)
         export NL_START_DAY=$(echo $L_START_DATE | cut -c7-8),$(echo $L_START_DATE | cut -c7-8)
         export NL_START_HOUR=$(echo $L_START_DATE | cut -c9-10),$(echo $L_START_DATE | cut -c9-10)
         export NL_END_YEAR=$(echo $L_END_DATE | cut -c1-4),$(echo $L_END_DATE | cut -c1-4)
         export NL_END_MONTH=$(echo $L_END_DATE | cut -c5-6),$(echo $L_END_DATE | cut -c5-6)
         export NL_END_DAY=$(echo $L_END_DATE | cut -c7-8),$(echo $L_END_DATE | cut -c7-8)
         export NL_END_HOUR=$(echo $L_END_DATE | cut -c9-10),$(echo $L_END_DATE | cut -c9-10)
         export NL_START_DATE=\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\',\'${L_START_YEAR}-${L_START_MONTH}-${L_START_DAY}_${L_START_HOUR}:00:00\'
         export NL_END_DATE=\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\',\'${L_END_YEAR}-${L_END_MONTH}-${L_END_DAY}_${L_END_HOUR}:00:00\'
         ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_nested_RT.ksh
#
         RANDOM=$$
         export TYPE=SERIAL
         export EXECUTE=real.exe
         if [[ ${TYPE} == PARALLEL ]]; then
            mpirun -np ${TASKS} ./${EXECUTE} > index_${RANDOM}.html 2>&1
         elif [[ ${TYPE} == SERIAL ]]; then
            ./${EXECUTE} > index_${RANDOM}.html 2>&1
         fi
#
         mv wrfinput_d${CR_DOMAIN} wrfinput_d${CR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
         mv wrfinput_d${FR_DOMAIN} wrfinput_d${FR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
         mv wrfbdy_d${CR_DOMAIN} wrfbdy_d${CR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
#         mv wrfbdy_d${FR_DOMAIN} wrfbdy_d${FR_DOMAIN}_$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} 0 -W 2>/dev/null)
         export P_DATE=$(${BUILD_DIR}/da_advance_time.exe ${P_DATE} ${LBC_FREQ} 2>/dev/null) 
      done
   fi
#
#########################################################################
#
# RUN EXO_COLDENS
#
#########################################################################
#
   if ${RUN_EXO_COLDENS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/exo_coldens ]]; then
         mkdir ${RUN_DIR}/${DATE}/exo_coldens
         cd ${RUN_DIR}/${DATE}/exo_coldens
      else
         cd ${RUN_DIR}/${DATE}/exo_coldens
      fi
#
# LINK NEEDED FILES
      export FILE_CR=wrfinput_d${CR_DOMAIN}
      export FILE_FR=wrfinput_d${FR_DOMAIN}
      rm -rf ${FILE_CR}
      rm -rf ${FILE_FR}
      ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
      ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
      export FILE=exo_coldens.nc
      rm -rf ${FILE}
      ln -sf ${EXPERIMENT_COLDENS_DIR}/${FILE} ${FILE}
      export FILE=exo_coldens.exe
      rm -rf ${FILE}
      ln -sf ${WES_COLDENS_DIR}/work/${FILE} ${FILE}
#
# CREATE INPUT FILE
      export FILE=exo_coldens.inp
      rm -rf ${FILE}
      cat << EOF > ${FILE}
&control
domains = 2,
/
EOF
#
# RUN exo_coldens
      ./exo_coldens.exe < exo_coldens.inp
#
# TEST WHETHER OUTPUT EXISTS
      export FILE_CR=exo_coldens_d${CR_DOMAIN}
      export FILE_FR=exo_coldens_d${FR_DOMAIN}
      if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
         echo EXO_COLDENS FAILED
         exit
      else
         echo EXO_COLDENS SUCCESS
      fi
   fi
#
#########################################################################
#
# RUN SEASONS_WES
#
#########################################################################
#
   if ${RUN_SEASON_WES}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/seasons_wes ]]; then
         mkdir ${RUN_DIR}/${DATE}/seasons_wes
         cd ${RUN_DIR}/${DATE}/seasons_wes
      else
         cd ${RUN_DIR}/${DATE}/seasons_wes
      fi
#
# LINK NEEDED FILES
      export FILE_CR=wrfinput_d${CR_DOMAIN}
      export FILE_FR=wrfinput_d${FR_DOMAIN}
      rm -rf ${FILE_CR}
      rm -rf ${FILE_FR}
      ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
      ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
      export FILE=season_wes_usgs.nc
      rm -rf ${FILE}
      ln -sf ${EXPERIMENT_COLDENS_DIR}/${FILE} ${FILE}
      export FILE=wesely.exe
      rm -rf ${FILE}
      ln -sf ${WES_COLDENS_DIR}/work/${FILE} ${FILE}
#
# CREATE INPUT FILE
      export FILE=wesely.inp
      rm -rf ${FILE}
      cat << EOF > ${FILE}
&control
domains = 2,
/
EOF
#
# RUN wesely
      ./wesely.exe < wesely.inp
#
# TEST WHETHER OUTPUT EXISTS
      export FILE_CR=wrf_season_wes_usgs_d${CR_DOMAIN}.nc
      export FILE_FR=wrf_season_wes_usgs_d${FR_DOMAIN}.nc
      if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
         echo WESELY FAILED
         exit
      else
         echo WESELY SUCCESS
      fi
   fi
