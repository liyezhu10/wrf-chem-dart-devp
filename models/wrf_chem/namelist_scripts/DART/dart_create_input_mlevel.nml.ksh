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
# Purpose: Create DART input.nml 
#
#########################################################################

echo off
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_assim_model_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_assim_tools_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_cov_cutoff_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_dart_to_wrf_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_ensemble_manager_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_filter_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_location_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_model_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_radar_mod_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_diag_mlevel_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_kind_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_selection_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_seq_coverage_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_seq_to_netcdf_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_sequence_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_sequence_tool_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_perfect_model_obs_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_preprocess_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_reg_factor_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_replace_wrf_fields_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_restart_file_tool_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_restart_file_utility_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_schedule_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_smoother_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_utilities_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_wrf_obs_preproc_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_wrf_to_dart_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_MOPITT_CO_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_IASI_CO_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_IASI_O3_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_AIRNOW_PM10_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_AIRNOW_PM25_nml.ksh
${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_obs_def_MODIS_AOD_nml.ksh
echo on
