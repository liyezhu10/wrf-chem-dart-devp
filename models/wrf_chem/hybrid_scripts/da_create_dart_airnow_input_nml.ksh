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
# Purpose: Script to create DART/WRF input.nml for Ave's 
# iasi_ascii_to_obs_seq fortran format conversion 
#
#########################################################################
#
# CREATE DART/WRF NAMELIST FILE
rm -f create_airnow_obs_nml.nl
touch create_airnow_obs_nml.nl
cat > create_airnow_obs_nml.nl << EOF
&create_airnow_obs_nml
   year0=${NL_YEAR}
   month0=${NL_MONTH}
   day0=${NL_DAY}
   hour0=${NL_HOUR}
   beg_year=${ASIM_MIN_YYYY}
   beg_mon=${ASIM_MIN_MM}
   beg_day=${ASIM_MIN_DD}
   beg_hour=${ASIM_MIN_HH}
   beg_min=${ASIM_MIN_MN}
   beg_sec=${ASIM_MIN_SS}
   end_year=${ASIM_MAX_YYYY}
   end_mon=${ASIM_MAX_MM}
   end_day=${ASIM_MAX_DD}
   end_hour=${ASIM_MAX_HH}
   end_min=${ASIM_MAX_MN}
   end_sec=${ASIM_MAX_SS}
   file_in=${NL_FILENAME}
   lat_mn=${NL_LAT_MN}
   lat_mx=${NL_LAT_MX}
   lon_mn=${NL_LON_MN}                                                                                      
   lon_mx=${NL_LON_MX}
   use_log_co=${NL_USE_LOG_CO}
   use_log_o3=${NL_USE_LOG_O3}
   use_log_nox=${NL_USE_LOG_NOX}
   use_log_so2=${NL_USE_LOG_SO2}
   use_log_pm10=${NL_USE_LOG_PM10}
   use_log_pm25=${NL_USE_LOG_PM25}
/
EOF
#
rm -f input.nml
touch input.nml
cat > input.nml << EOF
&obs_sequence_nml
   write_binary_obs_sequence   = .false.
/
&obs_kind_nml
/
&assim_model_nml
   write_binary_restart_files  =.true.
/
&model_nml
/
&location_nml
/
&utilities_nml
   TERMLEVEL                   = 1,
   logfilename                 = 'dart_log.out',
/
&preprocess_nml
   input_obs_kind_mod_file     = '../../obs_kind/DEFAULT_obs_kind_mod.F90',
   output_obs_kind_mod_file    = '../../obs_kind/obs_kind_mod.f90',
   input_obs_def_mod_file      = '../../obs_def/DEFAULT_obs_def_mod.F90',
   output_obs_def_mod_file     = '../../obs_def/obs_def_mod.f90',
   input_files                 = '../../obs_def/obs_def_reanalysis_bufr_mod.f90',
                                 '../../obs_def/obs_def_gps_mod.f90',
                                 '../../obs_def/obs_def_eval_mod.f90'
/
&merge_obs_seq_nml
   num_input_files             = 2,
   filename_seq                = 'obs_seq2008022206',obs_seq2008022212',
   filename_out                = 'obs_seq_ncep_2008022212'
/
EOF
