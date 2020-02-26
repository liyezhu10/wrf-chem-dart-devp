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

export YEAR=${YYYY}
export MONTH=${MM}
export DAY=${DD}
export HOUR=${HH}
export METDIR=${REAL_DIR}
export CHEMDIR=${RUN_DIR}/${DATE}/wrfchem_chem_icbc
export WRFINP_CR=wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
export WRFBDY_CR=wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00
cp ${METDIR}/${WRFINP_CR} ./.
cp ${METDIR}/${WRFBDY_CR} ./.
mv ${WRFINP_CR} wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
mv ${WRFBDY_CR} wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000
rm -f mozbc.ic.inp.set00
cat mozbc.ic.inp set00 > mozbc.ic.inp.set00
rm -f mozbc.bc.inp.set00
cat mozbc.bc.inp set00 > mozbc.bc.inp.set00
./run_mozbc_rt_CR.csh type=ic mozbc_inp=mozbc.ic.inp.set00 ens=000
./run_mozbc_rt_CR.csh type=bc mozbc_inp=mozbc.bc.inp.set00 ens=000
mv wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFINP_CR}
mv wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e000 ${WRFBDY_CR}
