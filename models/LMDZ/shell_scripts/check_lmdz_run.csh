#!/bin/csh
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#**************************************************************************
# Tarkeshwar Singh
# PhD, IIT Delhi
# Email: tarkphysics87@gmail.com
#
# PURPOSE: Check running status ( START and FINISH time) of LMDZ ensembles
# Some LMDZ members may take much time to finish due to slow computations on 
# some nodes. Sometimes few members hangs on the nodes.
# If this happen, allocate some extra free nodes in job.csh and put each 
# extra nodes in diffrent hostfiles. All allocated nodes for the job can be 
# seen by using command 'cat $PBS_NODEFILE  >PBS_NODEFILE' in job.csh script.
# e.g., If you are running 60 member and each on diffrent nodes then assign
# some extra nodes e.g. 65 nodes. edit run_lmdz.run to jump slow LMDZ members 
# on thses free nodes. for example: if LMDZ member 5 ( process id 4) is 
# slow then add follwong line
#  ............
#  ..............
#  else if ($1 == 4) then
#  mpirun -machinefile ../hostfile_freeNode -n 1 ./$gcm_exe
#  ............
#  ..............
# (These lines are already added and commented in run_lmdz.csh script. Uncomment 
# and edit it for particular slow member and their manually allocated hostfile)
#**************************************************************************
 
set num_ens = 60

set count = 0 

while ($count <= $num_ens) 
   
 if (-f `pwd`/advance_temp$count/gcm.log) then

   grep 'Simulation finished'  advance_temp$count/gcm.log > /dev/null

   if ( $status != 0) then
    echo "***************************************************"
    echo "LMDZ simulation NOT finished in DIR" advance_temp$count > tmp
    echo "LMDZ simulation NOT finished in DIR" advance_temp$count 
    #set  start_time = `cat advance_temp$count/time_omp | tail -1`
    set  start_time = `grep IST advance_temp$count/time_omp`
    echo  "START time=" $start_time
    echo  "CURR  time=" `date`
    set string = `grep "TKK" advance_temp$count/gcm.log | tail -1`
    echo "LMDZ current "$string[2-5] 
    echo "***************************************************"

   else
   # echo "LMDZ simulation finished in DIR" advance_temp$count
   endif
 else if ($count == 1) then 
  echo "DIR" advance_temp$count "NOT FOUND"
  exit
 endif

 @ count++

end

if (-f tmp) then
 \rm tmp 
 echo ""
   echo "***************************************************"  
   echo "*** Some/all LMDZ simulations NOT finished ****"
   echo "***************************************************"  
else
 echo ""
   echo "***************************************************"  
   echo "*** ALL LMDZ  simulations Finished ****"
   echo "***************************************************"  
 echo ""

endif

