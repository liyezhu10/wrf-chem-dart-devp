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
# Purpose: Create DART &wrf_to_dart_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &wrf_to_dart_nml
   dart_restart_name   = ${NL_DART_RESTART_NAME:-'null'},
   print_data_ranges   = ${NL_PRINT_DATA_RANGES:-.false.},
   debug               = ${NL_DEBUG:-.false.},
   add_emiss           = ${NL_ADD_EMISS:-.false.},
   use_log_co          = ${NL_USE_LOG_CO:-.false.},
   use_log_o3          = ${NL_USE_LOG_O3:-.false.},
   use_log_nox         = ${NL_USE_LOG_NOX:-.false.},
   use_log_so2         = ${NL_USE_LOG_SO2:-.false.},
   use_log_pm10        = ${NL_USE_LOG_PM10:-.false.},
   use_log_pm25        = ${NL_USE_LOG_PM25:-.false.},
   use_log_aod         = ${NL_USE_LOG_AOD:-.false.},
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
