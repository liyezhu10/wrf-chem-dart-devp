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
# Purpose: Create DART &obs_selection_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_selection_nml
   filename_seq          = ${NL_FILENAME_SEQ:-'obs_seq.out'},
   filename_seq_list     = ${NL_FILENAME_SEQ_LIST:-'null'},
   filename_out          = ${NL_FILENAME_OUT:-'obs_seq.processed'},
   selections_file       = ${NL_SELECTIONS_FILE:-'obsdef_mask.txt'},
   selections_is_obs_seq = ${NL_SELECTIONS_IS_OBS_SEQ:-.false.},
   print_only            = ${NL_PRINT_ONLY:-.false.},
   calendar              = ${NL_CALENDAR:-'gregorian'}, 
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
