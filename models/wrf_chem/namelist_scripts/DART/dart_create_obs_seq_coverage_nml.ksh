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
# Purpose: Create DART &obs_seq_coverage_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_seq_coverage_nml
   obs_sequence_list = ${NL_OBS_SEQUENCE_LIST:-'obs_coverage_list.txt'},
   obs_sequence_name = ${NL_OBS_SEQEUNCE_NAME:-'null'},
   obs_of_interest   = ${NL_OBS_OF_INTEREST:-'null'},
   textfile_out      = ${NL_TEXTFILE_OUT:-'obsdef_mask.txt'},
   netcdf_out        = ${NL_NETCDF_OUT:-'obsdef_mask.nc'},
   lonlim1           = ${NL_LONLIM1:-0.0},
   lonlim2           = ${NL_LONLIM2:-360.0},
   latlim1           = ${NL_LATLIM1:--90.0},
   latlim2           = ${NL_LATLIM2:-90.0},
   nTmin             = ${NL_NTMIN:-8},
   nTmax             = ${NL_NTMAX:-8},
   verbose           = ${NL_VERBOSE:-.false.},
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
