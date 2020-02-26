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

#-----------------------------------------------------------------------
# Script da_run_update_bc.ksh
#
# Purpose: Update WRF lateral boundary conditions to be consistent with 
# WRFVAR analysis.
#
#-----------------------------------------------------------------------
export NL_LOW_BDY_ONLY=false
export NL_UPDATE_LSM=false

cp -f $OPS_FORC_FILE real_output 
cp -f $DA_OUTPUT_FILE wrfvar_output
cp -f $BDYCDN_IN wrfbdy_d01_input
cp -f $BDYCDN_IN wrfbdy_d01

cat > parame.in << EOF
&control_param
  wrfvar_output_file = 'wrfvar_output'
  wrf_bdy_file       = 'wrfbdy_d01'
  wrf_input          = 'real_output'
  cycling = .${CYCLING}.
  debug   = .true.
  low_bdy_only = .${NL_LOW_BDY_ONLY}. 
  update_lsm = .${NL_UPDATE_LSM}. /
EOF

cp $BUILD_DIR/da_update_bc.exe .
./da_update_bc.exe

export RC=$?
if [[ $RC != 0 ]]; then
   echo "Update_bc failed with error $RC"
   exit $RC
else
   cp wrfbdy_d01 wrfbdy_d01_output
   cp wrfbdy_d01 $BDYCDN_OUT
fi
exit $?
