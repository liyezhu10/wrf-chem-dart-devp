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

##########################################################################
#
# Purpose: Set global environment variables for real_time_wrf_chem
#
#########################################################################
#
# RUNTIME DIRECTORIES
   export GEOGRID_DIR=${RUN_DIR}/geogrid
   export METGRID_DIR=${RUN_DIR}/${DATE}/metgrid
   export REAL_DIR=${RUN_DIR}/${DATE}/real
   export WRFCHEM_MET_IC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_ic
   export WRFCHEM_MET_BC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_bc
   export EXO_COLDENS_DIR=${RUN_DIR}/${DATE}/exo_coldens
   export SEASONS_WES_DIR=${RUN_DIR}/${DATE}/seasons_wes
   export WRFCHEM_BIO_DIR=${RUN_DIR}/${DATE}/wrfchem_bio
   export WRFCHEM_FIRE_DIR=${RUN_DIR}/${DATE}/wrfchem_fire
   export WRFCHEM_CHEMI_DIR=${RUN_DIR}/${DATE}/wrfchem_chemi
   export WRFCHEM_CHEM_EMISS_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
   export WRFCHEM_INITIAL_DIR=${RUN_DIR}/${INITIAL_DATE}/wrfchem_initial
   export WRFCHEM_CYCLE_CR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_cr
   export WRFCHEM_CYCLE_FR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_fr
   export WRFCHEM_LAST_CYCLE_CR_DIR=${RUN_DIR}/${PAST_DATE}/wrfchem_cycle_cr
   export PREPBUFR_MET_OBS_DIR=${RUN_DIR}/${DATE}/prepbufr_met_obs
   export MOPITT_CO_OBS_DIR=${RUN_DIR}/${DATE}/mopitt_co_obs
   export IASI_CO_OBS_DIR=${RUN_DIR}/${DATE}/iasi_co_obs
   export IASI_O3_OBS_DIR=${RUN_DIR}/${DATE}/iasi_o3_obs
   export OMI_NO2_OBS_DIR=${RUN_DIR}/${DATE}/omi_no2_obs
   export AIRNOW_CO_OBS_DIR=${RUN_DIR}/${DATE}/airnow_co_obs
   export AIRNOW_O3_OBS_DIR=${RUN_DIR}/${DATE}/airnow_o3_obs
   export AIRNOW_NO2_OBS_DIR=${RUN_DIR}/${DATE}/airnow_no2_obs
   export AIRNOW_SO2_OBS_DIR=${RUN_DIR}/${DATE}/airnow_so2_obs
   export AIRNOW_PM10_OBS_DIR=${RUN_DIR}/${DATE}/airnow_pm10_obs
   export AIRNOW_PM25_OBS_DIR=${RUN_DIR}/${DATE}/airnow_pm25_obs
   export PANDA_CO_OBS_DIR=${RUN_DIR}/${DATE}/panda_co_obs
   export PANDA_O3_OBS_DIR=${RUN_DIR}/${DATE}/panda_o3_obs
   export PANDA_PM25_OBS_DIR=${RUN_DIR}/${DATE}/panda_pm25_obs
   export MODIS_AOD_OBS_DIR=${RUN_DIR}/${DATE}/modis_aod_obs
   export COMBINE_OBS_DIR=${RUN_DIR}/${DATE}/combine_obs
   export PREPROCESS_OBS_DIR=${RUN_DIR}/${DATE}/preprocess_obs
   export WRFCHEM_CHEM_ICBC_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_icbc
   export DART_FILTER_DIR=${RUN_DIR}/${DATE}/dart_filter
   export UPDATE_BC_DIR=${RUN_DIR}/${DATE}/update_bc
   export BAND_DEPTH_DIR=${RUN_DIR}/${DATE}/band_depth
   export ENSEMBLE_MEAN_INPUT_DIR=${RUN_DIR}/${DATE}/ensemble_mean_input
   export ENSEMBLE_MEAN_OUTPUT_DIR=${RUN_DIR}/${DATE}/ensemble_mean_output
   export REAL_TIME_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_REAL_TIME
