#!/bin/csh

#BSUB -J mm.1_rma
#BSUB -o mm.1_rma.%J.log
#BSUB -e mm.1_rma.%J.err
#BSUB -q regular 
#BSUB -n 512 
#BSUB -R "span[ptile=16]"
#BSUB -P P86850054
#BSUB -W 01:30

# pin tasks to processors
setenv TARGET_CPU_LIST "-1"

module load job_memusage

time mpirun.lsf job_memusage.exe ./matrix_mm_5k

exit 0
