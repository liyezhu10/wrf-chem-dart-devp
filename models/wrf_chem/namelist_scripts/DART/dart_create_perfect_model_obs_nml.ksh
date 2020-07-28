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
# Purpose: Create DART &perfect_model_obs_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
&perfect_model_obs_nml
   first_obs_days               = ${NL_FIRST_OBS_DAYS:-148816},
   first_obs_seconds            = ${NL_FIRST_OBS_SECONDS:-75601},
   init_time_days               = ${NL_INIT_TIME_DAYS:--1},
   init_time_seconds            = ${NL_INIT_TIME_SECONDS:--1},
   last_obs_days                = ${NL_LAST_OBS_DAYS:-148817},
   obs_seq_in_file_name         = ${NL_OBS_SEQ_IN_FILE_NAME:-'obs_seq.in'},
   obs_seq_out_file_name        = ${NL_OBS_SEQ_OUT_FILE_NAME:-'obs_seq.out'},
   output_forward_op_errors     = ${NL_OUTPUT_FORWARD_OP_ERRORS:-.false.},
   output_timestamps            = ${NL_OUTPUT_TIMESTAMPS:-.false.},
   print_every_nth_obs          = ${NL_PRINT_EVERY_NTH_OBS:--1},
   trace_execution              = ${NL_TRACE_EXECUTION:-.false.},
   adv_ens_command              = ${NL_ADV_ENS_COMMAND:-'./advance_model.csh'},
   async                        = ${NL_ASYNC:-2},
   last_obs_seconds             = ${NL_LAST_OBS_SECONDS:-10800},
   output_interval              = ${NL_OUTPUT_INTERVAL:-1},
   output_restart               = ${NL_OUTPUT_RESTART:-.true.},
   restart_in_file_name         = ${NL_RESTART_IN_FILE_NAME:-'perfect_ics'},
   restart_out_file_name        = ${NL_RESTART_OUT_FILE_NAME:-'perfect_restart'},
   silence                      = ${NL_SILENCE:-.false.},
   start_from_restart           = ${NL_START_FROM_RESTART:-.true.},
/ 
EOF
#
# Append namelist section to input.nml
if [[ -f input.nml ]]; then
   cat input.nml_temp >> input.nml
   rm input.nml_temp
else
   mv input.nml_temp input.nml
fi
