#!/bin/ksh -aeux
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

#==========================================================
# Run the create_ict_aircraft_chem.pro in IDL
#==========================================================
#
export START_DATE=2014071406
export END_DATE=2014071906
export TIME_INC=6
export PROJ_NUMBER=P19010000
export DART_VERSION=DART_CHEM_MY_BRANCH
export WRFDA_VERSION=WRFDAv3.4_dmpar
#
# Experiment
#export EXPERIMENT=real_FRAPPE_CNTL_VARLOC
#export EXPERIMENT=real_FRAPPE_CNTL_NVARLOC
#export EXPERIMENT=real_FRAPPE_COnXX_VARLOC
export EXPERIMENT=real_FRAPPE_COnXX_VARLOC
#
export SCRATCH_DIR=/glade/scratch/mizzi
export ACD_DIR=/glade/p/acd/mizzi
export WORK_DIR=/glade/p/work/mizzi
export FRAPPE_DIR=/glade/p/FRAPPE
#
export TRUNK_DIR=${WORK_DIR}/TRUNK
export DART_DIR=${TRUNK_DIR}/${DART_VERSION}
export BUILD_DIR=${TRUNK_DIR}/${WRFDA_VERSION}/var/build
#
export RUN_DIR=${SCRATCH_DIR}/FRAPPE_DIAGNOSTICS
if [[ ! -e ${RUN_DIR} ]]; then
   mkdir -p ${RUN_DIR}
   cd ${RUN_DIR}
else
   cd ${RUN_DIR}
fi
cp ${DART_DIR}/models/wrf_chem/work/advance_time ./. 
#
# Get forecast files
export L_DATE=${START_DATE}
rm wrfout_d*
while [[ ${L_DATE} -le ${END_DATE} ]]; do
   export L_YY=$(echo $L_DATE | cut -c1-4)
   export L_MM=$(echo $L_DATE | cut -c5-6)
   export L_DD=$(echo $L_DATE | cut -c7-8)
   export L_HH=$(echo $L_DATE | cut -c9-10)
   export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
   cp ${FRAPPE_DIR}/${EXPERIMENT}/${L_DATE}/ensemble_mean/wrfinput_d01_mean wrfout_d01_${L_FILE_DATE}
   export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${TIME_INC} 2>/dev/null)
done
#
# Get idl code
rm -rf run_idl_code.pro
cp ${DART_DIR}/models/wrf_chem/run_frappe_diagnostics/create_ict_aircraft_chem.pro run_idl_code.pro
chmod +x run_idl_code.pro
#
# Create jobs script to run IASI O3 eps to ascii code
rm -rf job.ksh
rm -rf idl_*.err
rm -rf idl_*.out
touch job.ksh
RANDOM=$$
export JOBRND=idl_$RANDOM
cat <<EOFF >job.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 1                                  # number of total (MPI) tasks
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:10                              # wallclock time (minutes)
#BSUB -q geyser
#
idl <<EOF
.run run_idl_code.pro
exit
EOF
EOFF
#
bsub -K < job.ksh
# end of script
