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
# Purpose: Create DART &wrf_obs_preproc_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &wrf_obs_preproc_nml
   obs_boundary             = ${NL_OBS_BOUNDARY:-0.0},
   increase_bdy_error       = ${NL_INCREASE_BDY_ERROR:-.false.},
   maxobsfac                = ${NL_MAXOBSFAC:-2.5},
   obsdistbdy               = ${NL_OBSDISTBDY:-15.0},
   sfc_elevation_check      = ${NL_SFC_ELEVATION_CHECK:-.false.},
   sfc_elevation_tol        = ${NL_SFC_ELEVATION_TOL:-300.0},
   obs_pressure_top         = ${NL_OBS_PRESSURE_TOP:-0.0},
   obs_height_top           = ${NL_OBS_HEIGHT_TOP:-200000000000000.0},
   include_sig_data         = ${NL_INCLUDE_SIG_DATA:-.true.},
   tc_sonde_radii           = ${NL_TC_SONDE_RADII:--1.0},
   superob_aircraft         = ${NL_SUPEROB_ARICRAFT:-.false.},
   aircraft_horiz_int       = ${NL_AIRCRAFT_HORIZ_INT:-36.0},
   aircraft_pres_int        = ${NL_ARICRAFT_PRES_INT:-2500.0},
   superob_sat_winds        = ${NL_SUPEROB_SAT_WINDS:-.false.},
   sat_wind_horiz_int       = ${NL_SAT_WIND_HORIZ_INT:-100.0},
   sat_wind_pres_int        = ${NL_SAT_WIND_PRES_INT:-2500.0},
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
