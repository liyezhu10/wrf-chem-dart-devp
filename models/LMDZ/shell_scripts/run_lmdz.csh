#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# #**************************************************************************
#  Tarkeshwar Singh
#  PhD, IIT Delhi
#  Email: tarkphysics87@gmail.com
#  
# PURPOSE: copy required inputs and launch each LMDZ member
# #**************************************************************************

limit stacksize unlimited
limit stacksize unlimited
limit datasize unlimited

# If running LMDZ in shared memory parallelization (OMP), then set default OMP 
# parameters. Actual values will be  taken from Control_File.csh file
setenv OMP_NUM_THREADS 9
setenv OMP_STACKSIZE 500M

source ../Control_File.csh

## Link LMDZ def files
ln -sf $LMDZ_DIR/*.def .
ln -sf $LMDZ_DIR/$limit_file limit.nc
ln -sf $LMDZ_DIR/$gcm_exe .

## Determine model advanced date to tag LMDZ hist files
../trans_time
set adv_date = `cat times | tail -1`
echo $adv_date
set hh = `echo $adv_date | cut -c12-13`
echo $hh

set ens_member = `cat element`
echo $ens_member

# Default minimum LMDZ restart time is one day. For assimilation, It has been
# modified to run in subday restart mode e.g 6h. In subday restart, It save ozone
# variables at every 0000 UTC of each day and use the same for rest hours of the day.
# if hours is not 0000 UTC, then copy it for given meber from $HOME dir.

mv ../stok_paprs.dat_$ens_member  stok_paprs.dat

## Required for check_lmdz_run.csh script 
echo "************* ens_member = " $ens_member  > time_omp
echo "        ">> time_omp
date >> time_omp

#-----------------------------------Run LMDZ----------------------------------------------
# It has been observed that ensemble member 1 i.e. (task_id=0) become very slow
# compare to others probably  due to running on master node which involves in
# other task as well. It has been solved by assigning one extra (hostfile_lastFreeNode) 
# node and moving this member on extra node.  

if ($1 == 0) then
 mpirun -machinefile ../hostfile_lastFreeNode -n 1 ./$gcm_exe

 # use script check_lmdz_run.csh and check if any other members are also
 # slow compare to rest. It is expected that assign 3-4 more nodes than required.
 # Uncomment follwing and edit to jump particular slow members on extra nodes.

#else if ($1 == 58) then
# mpirun -machinefile ../host1 -n 1 ./$gcm_exe
#else if ($1 == 30) then
# mpirun -machinefile ../host2 -n 1 ./$gcm_exe
#else if ($1 == 46) then
# mpirun -machinefile ../host3 -n 1 ./$gcm_exe
#else if ($1 == 47) then
# mpirun -machinefile ../host4 -n 1 ./$gcm_exe
else
 
 # manually assing host file for each member
set task_id=$1
 @ host_count = ($task_id + 1 )
 cat ../PBS_NODEFILE | head -$host_count | tail -1 > hostfile
 mpirun -machinefile hostfile -n 1   ./$gcm_exe
endif
#-------------------------------------------------------------------------------------------
## Required for check_lmdz_run.csh script 
date >> time_omp
echo "        ">> time_omp
echo "        ">> time_omp

## renane restart and model outfiles
mv restart.nc start.nc
mv restartphy.nc startphy.nc

mv histhf_0000.nc ../histhf_$ens_member.nc_$hh
mv histins_0000.nc ../histins_$ens_member.nc_$hh

mv stok_paprs.dat ../stok_paprs.dat_$ens_member

## Save sortvar output at every model advance in GLOB_$ens_member file. It helps to monitor the model
## stability.
echo "*******************************************************************************" >> ../GLOB_$ens_member
cat times | head -1 >> ../GLOB_$ens_member
grep GLOB gcm.log >> ../GLOB_$ens_member
