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
# mopitt_ascii_to_obs_seq fortran format conversion 
#
#########################################################################
#
# CREATE DART/WRF NAMELIST FILE
rm -f input.nml
touch input.nml
cat > input.nml << EOF
&create_mopitt_obs_nml
   year                        = ${NL_YEAR}
   month                       = ${NL_MONTH}
   day                         = ${NL_DAY}
   hour                        = ${NL_HOUR}
   bin_beg                     = ${NL_BIN_BEG}
   bin_end                     = ${NL_BIN_END}
   filedir                     = ${NL_FILEDIR}
   filename                    = ${NL_FILENAME}
   MOPITT_CO_retrieval_type    = ${NL_MOPITT_CO_RETRIEVAL_TYPE}
   fac_obs_error               = ${NL_FAC_OBS_ERROR}
   use_cpsr_co_trunc           = ${NL_USE_CPSR_CO_TRUNC}
   cpsr_co_trunc_lim           = ${NL_CPSR_CO_TRUNC_LIM}
   mopitt_co_vloc              = ${NL_MOPITT_CO_VLOC}
   use_log_co                  = ${NL_USE_LOG_CO}
   lon_min                     = ${NNL_MIN_LON}
   lon_max                     = ${NNL_MAX_LON}
   lat_min                     = ${NNL_MIN_LAT}
   lat_max                     = ${NNL_MAX_LAT}
/
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
&obs_def_MOPITT_CO_nml
   MOPITT_CO_retrieval_type    = ${NL_MOPITT_CO_RETRIEVAL_TYPE:-'RETR'},
   use_log_co                  = ${NL_USE_LOG_CO:-.false.},
/ 
&obs_def_IASI_CO_nml
   IASI_CO_retrieval_type      = ${NL_IASI_CO_RETRIEVAL_TYPE:-'RETR'},
   use_log_co                  = ${NL_USE_LOG_CO:-.false.},
/ 
EOF

