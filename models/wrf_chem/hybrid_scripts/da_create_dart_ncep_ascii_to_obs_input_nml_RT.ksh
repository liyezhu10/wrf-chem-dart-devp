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
# Purpose: Script to create DART/WRF Namelist 

#########################################################################
#
# CREATE DART/WRF NAMELIST FILE
rm -f input.nml
touch input.nml
cat > input.nml << EOF
&ncepobs_nml
   year       = ${D_YYYY},
   month      = ${D_MM},
   day        = ${D_DD},
   tot_days   = 1,
   max_num    = 800000,
   select_obs = 0,
   daily_file = .false.,
   ObsBase = 'temp_obs.',
   ADPUPA = .true., 
   AIRCFT = .true., 
   AIRCAR = .true., 
   SATEMP = .true., 
   SFCSHP = .true., 
   ADPSFC = .true., 
   SATWND = .true.,
   obs_U  = .true., 
   obs_V  = .true., 
   obs_T  = .true.,
   obs_PS = .true.,
   obs_QV = .true.,
   lon1   = ${NL_MIN_LON},
   lon2   = ${NL_MAX_LON},
   lat1   = ${NL_MIN_LAT},
   lat2   = ${NL_MAX_LAT},
   obs_time = .true.  /

&obs_sequence_nml
   write_binary_obs_sequence = .false.  /

&assim_model_nml
   write_binary_restart_files = .true.  /

&utilities_nml
   TERMLEVEL = 1,
   module_details = .false.,
   logfilename = 'dart_log.out'  /

&model_nml
   /

&location_nml
   /

&obs_kind_nml
   /

&obs_def_gps_nml
    max_gpsro_obs = 100000 /

&preprocess_nml
    input_obs_kind_mod_file = '${DART_DIR}/obs_kind/DEFAULT_obs_kind_mod.F90',
    output_obs_kind_mod_file= '${DART_DIR}/obs_kind/obs_kind_mod.f90',
    input_obs_def_mod_file  = '${DART_DIR}/obs_def/DEFAULT_obs_def_mod.F90',
    output_obs_def_mod_file = '${DART_DIR}/obs_def/obs_def_mod.f90',
    input_files             = '${DART_DIR}/obs_def/obs_def_reanalysis_bufr_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_altimeter_mod.f90',
                              '${DART_DIR}/obs_def/obs_def_gps_mod.f90' /
EOF
