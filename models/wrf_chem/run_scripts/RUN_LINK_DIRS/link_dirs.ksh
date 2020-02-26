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
export WRFDA_VERSION=WRFDAv3.9.1.1_dmpar
export BUILD_DIR=/projects/mizzi/TRUNK/${WRFDA_VERSION}/var/build
#
export EXP_AIRNOW_CONTROL=/real_FRAPPE_AIRNOW_CONTROL
export EXP_AIRNOW_CO=/real_FRAPPE_RETR_AIR_CO
export EXP_AIRNOW_O3=/real_FRAPPE_RETR_AIR_O3
export EXP_AIRNOW_NO2=/real_FRAPPE_RETR_AIR_NO2
export EXP_AIRNOW_SO2=/real_FRAPPE_RETR_AIR_SO2
#
export DESTINATION_PATH=/scratch/summit/mizzi${EXP_AIRNOW_CONTROL}
export SOURCE_PATH=/scratch/summit/mizzi${EXP_AIRNOW_O3}
#
export DATE_STR=2014071412
export DATE_END=2014071412
#
export CYCLE_PERIOD=6
export L_DATE=${DATE_STR}
#
while [[ ${L_DATE} -le ${DATE_END} ]] ; do
   if [[ ! -d ${DESTINATION_PATH}/${L_DATE} ]]; then mkdir -p ${DESTINATION_PATH}/${L_DATE}; fi
   cd ${DESTINATION_PATH}/${L_DATE}
   ln -sf ${SOURCE_PATH}/${L_DATE}/airnow_co_obs airnow_co_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/airnow_no2_obs airnow_no2_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/airnow_o3_obs airnow_o3_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/airnow_so2_obs airnow_so2_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/combine_obs combine_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/exo_coldens exo_coldens
   ln -sf ${SOURCE_PATH}/${L_DATE}/metgrid metgrid
   ln -sf ${SOURCE_PATH}/${L_DATE}/prepbufr_met_obs prepbufr_met_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/preprocess_obs preprocess_obs
   ln -sf ${SOURCE_PATH}/${L_DATE}/real real
   ln -sf ${SOURCE_PATH}/${L_DATE}/seasons_wes seasons_wes
   ln -sf ${SOURCE_PATH}/${L_DATE}/ungrib ungrib
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_bio wrfchem_bio
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_chem_emiss wrfchem_chem_emiss
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_chemi wrfchem_chemi
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_chem_icbc wrfchem_chem_icbc
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_fire wrfchem_fire
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_met_bc wrfchem_met_bc
   ln -sf ${SOURCE_PATH}/${L_DATE}/wrfchem_met_ic wrfchem_met_ic
   cd ${DESTINATION_PATH}
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${CYCLE_PERIOD} 2>/dev/null)
done
