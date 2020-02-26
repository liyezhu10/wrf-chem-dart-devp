#!/bin/ksh -aux
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
export NUM_MEMBERS=30
export WRFDA_VERSION=WRFDAv3.9.1.1_dmpar
export BUILD_DIR=/projects/mizzi/TRUNK/${WRFDA_VERSION}/var/build
#
export EXP=/real_FRAPPE_RETR_MOP_CO
#export EXP=/real_FRAPPE_CPSR_MOP_CO
export EXP=/OUTLIER_4/real_FRAPPE_RETR_IAS_CO
#export EXP=/real_FRAPPE_CPSR_IAS_CO
#export EXP=/real_FRAPPE_RETR_AIR_CO
#
# ARCHIVING DONE
export EXP=/real_FRAPPE_CONTROL
#
export SOURCE_PATH=/scratch/summit/mizzi${EXP}
#
export DATE_STR=2014071606
export DATE_END=2014071818
#
# List of variables to archive
export ARCHIVE_VARS=XLAT,XLONG,ZNU,U,V,W,PH,PHB,T,T_INIT,MU,MUB,P,ALT,PB,PSFC,QVAPOR,T00,P00,CLDFRA,XLAND,PBLH,CLDFRA2,o3,no,no2,no3,so2,sulf,co,hcho,ch4,P25,P10,BC1,BC2,OC1,OC2,DUST_1,DUST_2,DUST_3,DUST_4,DUST_5,SEAS_1,SEAS_2,SEAS_3,SEAS_4,PM2_5_DRY,PM10
#
export CYCLE_PERIOD=6
export FCST_PERIOD=6
export DOMAIN=01
#
# Subset WRF-Chem output files in the following directories
# ensemble_mean_input
# ensemble_mean_output
# wrfchem_initial
# wrfchem_cycle_cr
# ensmean_cycle_fr
#
export L_DATE=${DATE_STR}
while [[ ${L_DATE} -le ${DATE_END} ]] ; do
#
   cd ${SOURCE_PATH}/${L_DATE}
#
   rm -rf ensemble_mean_input
   rm -rf wrfchem_chem_icbc
   rm -rf wrfchem_chem_emiss 
   if [[ -e dart_filter ]]; then
      cd dart_filter
      rm -rf obs_seq.out
      rm -rf prior_inflate_ic_new
      rm -rf prior_inflate_ic_old
      cd ../
   fi
#
#   if [[ -e ensemble_mean_input ]]; then
#      cd ensemble_mean_input
#      export FILE_IN=wrfinput_d${DOMAIN}_mean
#      export FILE_OUT=wrfinput_d${DOMAIN}_mean.archive
#      ncks -v ${ARCHIVE_VARS} ${FILE_IN} ${FILE_OUT} 
#      export FILE_IN=wrfinput_d${DOMAIN}_sprd
#      export FILE_OUT=wrfinput_d${DOMAIN}_sprd.archive
#      ncks -v ${ARCHIVE_VARS} ${FILE_IN} ${FILE_OUT} 
#      cd ../
#   fi
#
   if [[ -e ensemble_mean_output ]]; then
      cd ensemble_mean_output
      let IHR=0
      while [[ ${IHR} -le ${FCST_PERIOD} ]]; do  
         export LL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${IHR} 2>/dev/null)
         export FILE_IN=wrfout_d${DOMAIN}_${LL_DATE}_mean
         export FILE_OUT=wrfout_d${DOMAIN}_${LL_DATE}_mean.archive
         if [[ -e ${FILE_IN} ]]; then
            rm -rf ${FILE_OUT}
            ncks -v ${ARCHIVE_VARS} ${FILE_IN} ${FILE_OUT} 
            RC=$?
            if [[ ${RC} -eq 0 ]];then
               rm -rf ${FILE_IN}
               mv ${FILE_OUT} ${FILE_IN}
            fi
         fi
         let IHR=${IHR}+1
      done
      cd ../
   fi
#
   if [[ -e wrfchem_initial ]]; then
      cd wrfchem_initial
      let IMEM=1
      while [[ ${IMEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${IMEM}
         if [[ ${IMEM} -lt 100 ]]; then export CMEM=e0${IMEM}; fi
         if [[ ${IMEM} -lt 10 ]]; then export CMEM=e00${IMEM}; fi
         cd run_${CMEM}
         rm -rf ubvals_b40.20th.track1_*
         rm -rf wrfbdy_d01
         rm -rf wrfinput_d01
         let IHR=0
         while [[ ${IHR} -le ${FCST_PERIOD} ]]; do  
            export LL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${IHR} 2>/dev/null)
            export YYYY=$(echo $LL_DATE | cut -c1-4)
            export YY=$(echo $LL_DATE | cut -c3-4)
            export MM=$(echo $LL_DATE | cut -c5-6)
            export DD=$(echo $LL_DATE | cut -c7-8)
            export HH=$(echo $LL_DATE | cut -c9-10)
            export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
            export FILE_IN=wrfout_d${DOMAIN}_${FILE_DATE}
            export FILE_OUT=wrfout_d${DOMAIN}_${FILE_DATE}.archive
            if [[ -e ${FILE_IN} ]]; then
               rm -rf ${FILE_OUT}
               ncks -v ${ARCHIVE_VARS} ${FILE_IN} ${FILE_OUT} 
               RC=$?
               if [[ ${RC} -eq 0 ]];then
                  rm -rf ${FILE_IN}
                  mv ${FILE_OUT} ${FILE_IN}
               fi
            fi
            let IHR=${IHR}+1
         done
         cd ../
         let IMEM=${IMEM}+1
      done
      cd ../
   fi
#
   if [[ -e wrfchem_cycle_cr ]]; then
      cd wrfchem_cycle_cr
      let IMEM=1
      while [[ ${IMEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${IMEM}
         if [[ ${IMEM} -lt 100 ]]; then export CMEM=e0${IMEM}; fi
         if [[ ${IMEM} -lt 10 ]]; then export CMEM=e00${IMEM}; fi
         cd run_${CMEM}
         rm -rf ubvals_b40.20th.track1_*
         rm -rf wrfbdy_d01
         rm -rf wrfinput_d01
         let IHR=0
         while [[ ${IHR} -le ${FCST_PERIOD} ]]; do  
            export LL_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${IHR} 2>/dev/null)
            export YYYY=$(echo $LL_DATE | cut -c1-4)
            export YY=$(echo $LL_DATE | cut -c3-4)
            export MM=$(echo $LL_DATE | cut -c5-6)
            export DD=$(echo $LL_DATE | cut -c7-8)
            export HH=$(echo $LL_DATE | cut -c9-10)
            export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
            export FILE_IN=wrfout_d${DOMAIN}_${FILE_DATE}
            export FILE_OUT=wrfout_d${DOMAIN}_${FILE_DATE}.archive
            if [[ -e ${FILE_IN} ]]; then
               rm -rf ${FILE_OUT} 
               ncks -v ${ARCHIVE_VARS} ${FILE_IN} ${FILE_OUT} 
               RC=$?
               if [[ ${RC} -eq 0 ]];then
                  rm -rf ${FILE_IN}
                  mv ${FILE_OUT} ${FILE_IN}
               fi
            fi
            let IHR=${IHR}+1
         done
         cd ../
         let IMEM=${IMEM}+1
      done
      cd ../
   fi
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${CYCLE_PERIOD} 2>/dev/null)
done
