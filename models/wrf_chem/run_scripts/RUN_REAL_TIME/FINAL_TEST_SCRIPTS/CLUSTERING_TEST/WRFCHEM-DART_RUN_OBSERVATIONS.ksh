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

# GENERATE BIO, FIRE, and ANTHRO EMISSIONS
   rm -rf jobs.file
   cat << EOF > jobs.file
0 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MOPITT_CO.ksh
1 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_IASI_CO.ksh
2 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_IASI_O3.ksh
3 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_OMI_NO2.ksh
4 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_CO.ksh
5 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_O3.ksh
6 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_NO2.ksh
7 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_SO2.ksh
8 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_PM10.ksh
9 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_PM25.ksh
10 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_CO.ksh
11 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_O3.ksh
12 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_PM25.ksh
13 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MODIS_AOD.ksh
14 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MET_OBS.ksh
EOF
#
   srun --multi-prog ./jobs.file
#
