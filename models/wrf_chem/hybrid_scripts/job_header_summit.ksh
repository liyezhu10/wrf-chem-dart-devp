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

rm -rf job.header
cat << EOF > job.header
#!/bin/ksh -aeux
#SBATCH --job-name $1
#SBATCH --output $1.log
#SBATCH --account $6
#SBATCH --qos $2
#SBATCH --time $3
#SBATCH --nodes $4
#SBATCH --ntasks $5
#SBATCH --partition shas
#
EOF

rm -rf job.ksh
cat job.header $7 > job.ksh
rm -rf job.header
