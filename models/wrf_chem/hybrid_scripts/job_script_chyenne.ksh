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

export JOBID=$1
export PROJECT_NUMBER=$2
export TIME_LIMIT=$3
export CLASS=$4
export NODES=$5
export CORES=$6
export PROCS=$7
export EXECUTE=$8
#
if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
touch job.ksh
cat << EOF > job.ksh
#!/bin/ksh -aeux
#PBS -N ${JOBID}
#PBS -A ${PROJECT_NUMBER}
#PBS -l walltime=${TIME_LIMIT}
#PBS -q ${CLASS}
#PBS -j oe
#PBS -l select=${NODES}:ncpus=${CORES}:mpiprocs=${PROCS}
mpiexec_mpt ./${EXECUTE}  > index.html 2>&1 
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
