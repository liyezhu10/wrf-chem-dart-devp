#!/bin/ksh -aeux
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
export JOBID=$1
export CLASS=$2
export TIME_LIMIT=$3
export NODES=$4
export TASKS=$5
export EXECUTE=$6
export TYPE=$7
export ACCOUNT=$8
#
rm -rf job.ksh
touch job.ksh
#
if [[ ${TYPE} == PARALLEL ]]; then
   cat << EOF > job.ksh
#!/bin/ksh -aeux
#SBATCH --account ${ACCOUNT}
#SBATCH --job-name ${JOBID}
#SBATCH --qos ${CLASS}
#SBATCH --time ${TIME_LIMIT}
#SBATCH --output ${JOBID}.log
#SBATCH --nodes ${NODES}
#SBATCH --ntasks ${TASKS}
#SBATCH --partition shas
mpirun -np \${SLURM_NTASKS} ./${EXECUTE} > index.html 2>&1
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
#
elif [[ ${TYPE} == SERIAL ]]; then
   cat << EOF > job.ksh
#!/bin/ksh -aeux
#SBATCH --account ${ACCOUNT}
#SBATCH --job-name ${JOBID}
#SBATCH --qos ${CLASS}
#SBATCH --time ${TIME_LIMIT}
#SBATCH --output ${JOBID}.log
#SBATCH --nodes ${NODES}
#SBATCH --ntasks ${TASKS}
#SBATCH --partition shas
./${EXECUTE} > index.html 2>&1
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
#
else
   echo 'APM: Error is job script - Not SERIAL or PARALLEL '
   exit
fi
#
# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/hybrid_scripts/job_script_summit.ksh $
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
# $Revision: 13133 $
# $Date: 2019-04-25 15:47:54 -0600 (Thu, 25 Apr 2019) $
