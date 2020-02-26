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

##########################################################################
# Purpose: Set global environment variables for real_time_wrf_chem
#
#########################################################################
#
# DEPENDENT INPUT DATA DIRECTORIES:
   export EXPERIMENT_DIR=${SCRATCH_DIR}
   export RUN_DIR=${EXPERIMENT_DIR}/real_FRAPPE_RETR_MOP_CO_CLUSTER
   export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
   export WPS_GEOG_DIR=${INPUT_DATA_DIR}/${WPS_GEOG_VER}
   export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
   export WRFDA_DIR=${TRUNK_DIR}/${WRFDA_VER}
   export BUILD_DIR=${WRFDA_DIR}/var/da
   export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
   export HYBRID_SCRIPTS_DIR=${DART_DIR}/models/wrf_chem/hybrid_scripts
   export ADJUST_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_EMISS_INV
   export WES_COLDENS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_WES_COLDENS
   export MEGAN_BIO_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_MEGAN_BIO
   export FINN_FIRE_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_FINN_FIRE
   export EXPERIMENT_DATA_DIR=${INPUT_DATA_DIR}/FRAPPE_REAL_TIME_DATA
   export MOZBC_DATA_DIR=${EXPERIMENT_DATA_DIR}/mozart_forecasts
   export EXPERIMENT_STATIC_FILES=${EXPERIMENT_DATA_DIR}/static_files
   export EXPERIMENT_WRFCHEMI_DIR=${EXPERIMENT_DATA_DIR}/anthro_emissions
   export EXPERIMENT_WRFFIRECHEMI_DIR=${EXPERIMENT_DATA_DIR}/fire_emissions
   export EXPERIMENT_WRFBIOCHEMI_DIR=${EXPERIMENT_DATA_DIR}/bio_emissions
   export EXPERIMENT_COLDENS_DIR=${EXPERIMENT_DATA_DIR}/wes_coldens
   export EXPERIMENT_PREPBUFR_DIR=${EXPERIMENT_DATA_DIR}/met_obs_prep_data
   export EXPERIMENT_MOPITT_CO_DIR=${EXPERIMENT_DATA_DIR}/mopitt_co_hdf_data
   export EXPERIMENT_IASI_CO_DIR=${EXPERIMENT_DATA_DIR}/iasi_co_hdf_data
   export EXPERIMENT_IASI_O3_DIR=${EXPERIMENT_DATA_DIR}/iasi_o3_hdf_data
   export EXPERIMENT_OMI_NO2_DIR=${EXPERIMENT_DATA_DIR}/omi_no2_obs_seq_data
   export EXPERIMENT_AIRNOW_DIR=${EXPERIMENT_DATA_DIR}/airnow_csv_data
   export EXPERIMENT_MODIS_AOD_DIR=${EXPERIMENT_DATA_DIR}/modis_aod_hdf_data
   export EXPERIMENT_GFS_DIR=${EXPERIMENT_DATA_DIR}/gfs_forecasts
   export EXPERIMENT_DUST_DIR=${EXPERIMENT_DATA_DIR}/dust_fields
   export EXPERIMENT_HIST_IO_DIR=${EXPERIMENT_DATA_DIR}/hist_io_files
   export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
   export BE_DIR=${WRFDA_DIR}/var/run
   export PERT_CHEM_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT
   export PERT_CHEM_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
   export RUN_BAND_DEPTH_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_BAND_DEPTH
