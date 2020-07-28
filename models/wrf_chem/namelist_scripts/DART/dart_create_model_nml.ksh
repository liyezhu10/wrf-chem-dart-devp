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
# Purpose: Create DART &model_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &model_nml
   add_emiss                   = ${NL_ADD_EMISS:-.false.},
   use_varloc                  = ${NL_USE_VARLOC:-.true.},
   use_indep_chem_assim        = ${NL_USE_INDEP_CHEM_ASSIM:-.false.},
   default_state_variables     = ${NL_DEFAULT_STATE_VARIABLES:-.false.},
   conv_state_variables        = ${NL_CONV_STATE_VARIABLES:-"null"}
   emiss_chemi_variables       = ${NL_EMISS_CHEMI_VARIABLES:-"null"}
   emiss_firechemi_variables   = ${NL_EMISS_FIRECHEMI_VARIABLES:-"null"}
   wrf_state_bounds            = ${NL_WRF_STATE_BOUNDS:-"null"}
   num_domains                 = ${NL_NUM_DOMAINS:-1},
   calendar_type               = ${NL_CALENDAR_TYPE:-3},
   assimilation_period_seconds = ${NL_ASSIMILATION_PERIOD_SECONDS:-21600},
   vert_localization_coord     = ${NL_VERT_LOCALIZATION_COORD:-3},
   center_search_half_length   = ${NL_CENTER_SEARCH_HALF_LENGTH:-500000.0},
   center_spline_grid_scale    = ${NL_CENTER_SPLINE_GRID_SCALE:-10},   
   sfc_elev_max_diff           = ${NL_SFC_ELEV_MAX_DIFF:-100.0},
   circulation_pres_level      = ${NL_CIRCULATION_PRES_LEVEL:-80000.0},
   circulation_radius          = ${NL_CIRCULATION_RADIUS:-108000.0},
   allow_obs_below_vol         = ${NL_ALLOW_OBS_BELOW_VOL:-.false.},
   output_state_vector         = ${NL_OUTPUT_STATE_VECTOR:-.false.},
   use_log_aod                 = ${NL_USE_LOG_AOD:-.false.},
   use_log_co                  = ${NL_USE_LOG_CO:-.false.},
   use_log_nox                 = ${NL_USE_LOG_NOX:-.false.},
   use_log_o3                  = ${NL_USE_LOG_O3:-.false.},
   use_log_pm10                = ${NL_USE_LOG_PM10:-.false.},
   use_log_pm25                = ${NL_USE_LOG_PM25:-.false.},
   use_log_so2                 = ${NL_USE_LOG_SO2:-.false.},
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
