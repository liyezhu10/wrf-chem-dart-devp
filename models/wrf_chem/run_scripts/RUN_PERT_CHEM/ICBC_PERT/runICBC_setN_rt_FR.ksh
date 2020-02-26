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

unalias ls
export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
export IC_METDIR=${WRFCHEM_MET_IC_DIR}
export BC_METDIR=${WRFCHEM_MET_BC_DIR}
export CHEMDIR=${RUN_DIR}/${YEAR}${MONTH}${DAY}${HOUR}/wrfchem_chem_icbc
for ((i = 1; i <= ${NUM_MEMBERS}; i += 1))
do
   if [ "$i" -lt "10"  ]; then 
      export IENS="00"${i}
   fi
   if [ "$i" -lt "100"  ]; then
      if [ "$i" -ge "10"  ]; then
         export IENS="0"${i}
      fi
   fi
   if [ "$i" -lt "1000"  ]; then
      if [ "$i" -ge "100"  ]; then
         export IENS=${i}
      fi
   fi
   export WRFINP=wrfinput_d02_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e${IENS}
   echo 'get inp'
   cp ${IC_METDIR}/${WRFINP} ./.
   ls set${i}
   rm -f mozbc.ic.inp.set${i}
   cat mozbc.ic.inp set${i} > mozbc.ic.inp.set${i}
   echo  'run mozbc'
   ./run_mozbc_rt_FR.csh type=ic mozbc_inp=mozbc.ic.inp.set${i} ens=${IENS}
   echo 'put files'
   echo 'OK'
done
