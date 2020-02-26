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
# Purpose: Create DART &schedule_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &schedule_nml
   calendar                  = ${NL_CALENDAR:-'Gregorian'},  
   first_bin_start           = ${NL_FIRST_BIN_START_YY:-2008}, ${NL_FIRST_BIN_START_MM:-6}, ${NL_FIRST_BIN_START_DD:-1}, ${NL_FIRST_BIN_START_HH:-0}, ${NL_FIRST_BIN_START_MN:-0}, ${NL_FIRST_BIN_START_SS:-0},
   first_bin_end             = ${NL_FIRST_BIN_END_YY:-2008}, ${NL_FIRST_BIN_END_MM:-6}, ${NL_FIRST_BIN_END_DD:-1}, ${NL_FIRST_BIN_END_HH:-0}, ${NL_FIRST_BIN_END_MN:-0}, ${NL_FIRST_BIN_END_SS:-0},
   last_bin_end              = ${NL_LAST_BIN_END_YY:-2008}, ${NL_LAST_BIN_END_MM:-6}, ${NL_LAST_BIN_END_DD:-1}, ${NL_LAST_BIN_END_HH:-0}, ${NL_LAST_BIN_END_MN:-0}, ${NL_LAST_BIN_END_SS:-0},
   bin_interval_days         = ${NL_BIN_INTERVAL_DAYS:-0},
   bin_interval_seconds      = ${NL_BIN_INTERVAL_SECONDS:-21600},   
   max_num_bins              = ${NL_MAX_NUM_BINS:-1000},
   print_table               = ${NL_PRINT_TABLE:-.true.},
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
