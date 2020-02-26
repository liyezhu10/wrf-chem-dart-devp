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

#
export PROJ_NUMBER=P19010000
#
# compile code
ifort -check none -CB -C plot_Profile_MOPITT.f90 -o plot_Profile_MOPITT.exe -I${NETCDF}/include -L${NETCDF}/lib -lnetcdff -lnetcdf
#
# Create job script 
RANDOM=$$
export JOBRND=jobPRMOP_$RANDOM
rm -rf jobPRMOP_*.ksh
touch ${JOBRND}.ksh
cat << EOF >${JOBRND}.ksh
#!/bin/ksh -aeux
#BSUB -P ${PROJ_NUMBER}
#BSUB -n 32                                  # number of total (MPI) tasks
#BSUB -R "span[ptile=8]"
#BSUB -J ${JOBRND}                          # job name
#BSUB -o ${JOBRND}.out                      # output filename
#BSUB -e ${JOBRND}.err                      # error filename
#BSUB -W 00:29                              # wallclock time (minutes)
#BSUB -q regular
#
rm -rf jobPRMOP_*.out
rm -rf jobPRMOP_*.err
rm -rf jobPRMOP_*.index
./plot_Profile_MOPITT.exe > ${JOBRND}.index 2>&1 
#
export RC=\$?     
if [[ -f JOB_SUCCESS ]]; then rm -rf JOB_SUCCESS; fi     
if [[ -f JOB_FAILED ]]; then rm -rf JOB_FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch JOB_SUCCESS
else
   touch JOB_FAILED 
   exit
fi
EOF
#
# Submit convert file script for each and wait until job completes
bsub -K < ${JOBRND}.ksh 
rm -rf ${JOBRND}.ksh
#
exit
