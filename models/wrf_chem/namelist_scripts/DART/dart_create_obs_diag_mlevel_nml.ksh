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
# Purpose: Create DART &obs_diag_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_diag_nml
   mlevel_edges               = ${NL_MLEVEL},
   obs_sequence_name          = ${NL_OBS_SEQUENCE_NAME:-'obs_seq.final'},
   obs_sequence_list          = ${NL_OBS_SEQUENCE_LIST:-'file_list.txt'},
   first_bin_center           = ${NL_FIRST_BIN_CENTER_YY:-2008}, ${NL_FIRST_BIN_CENTER_MM:-6}, ${NL_FIRST_BIN_CENTER_DD:-1}, ${NL_FIRST_BIN_CENTER_HH:-0}, ${NL_FIRST_BIN_CENTER_MN:-0}, ${NL_FIRST_BIN_CENTER_SS:-0},
   last_bin_center            = ${NL_LAST_BIN_CENTER_YY:-2008}, ${NL_LAST_BIN_CENTER_MM:-6}, ${NL_LAST_BIN_CENTER_DD:-1}, ${NL_LAST_BIN_CENTER_HH:-0}, ${NL_LAST_BIN_CENTER_MN:-0}, ${NL_LAST_BIN_CENTER_SS:-0},
   bin_separation             = ${NL_BIN_SEPARATION_YY:-0}, ${NL_BIN_SEPARATION_MM:-0}, ${NL_BIN_SEPARATION_DD:-0}, ${NL_BIN_SEPARATION_HH:-6}, ${NL_BIN_SEPARATION_MN:-0}, ${NL_BIN_SEPARATION_SS:-0},
   bin_width                  = ${NL_BIN_WIDTH_YY:-0}, ${NL_BIN_WIDTH_MM:-0}, ${NL_BIN_WIDTH_DD:-0}, ${NL_BIN_WIDTH_HH:-6}, ${NL_BIN_WIDTH_MN:-0}, ${NL_BIN_WIDTH_SS:-0},
   time_to_skip               = ${NL_TIME_TO_SKIP_YY:-0}, ${NL_TIME_TO_SKIP_MM:-0}, ${NL_TIME_TO_SKIP_DD:-0}, ${NL_TIME_TO_SKIP_HH:-0}, ${NL_TIME_TO_SKIP_MN:-0}, ${NL_TIME_TO_SKIP_SS:-0},
   max_num_bins               = ${NL_MAX_NUM_BINS:-1000},
   Nregions                   = ${NL_NREGIONS:-1},
   lonlim1                    = ${NL_LONLIM1:-235.0},
   lonlim2                    = ${NL_LONLIM2:-295.0},
   latlim1                    = ${NL_LATLIM1:-25.0},
   latlim2                    = ${NL_LATLIM2:-55.0},
   reg_names                  = ${NL_REG_NAMES:-'North America'},
   print_mismatched_locs      = ${NL_PRINT_MISMATCHED_LOCS:-.false.},
   verbose                    = ${NL_VERBOSE:-.false.},
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
