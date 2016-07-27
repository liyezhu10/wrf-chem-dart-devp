#!/bin/csh
#==================================================================
#==================================================================


setenv OMP_NUM_THREADS 16

#setenv XLSMPOPTS "startproc=0:stride=2:stack=128000000"

./wrf_emp_localization >>&! zzz_omp_landsfcaltimeter 


