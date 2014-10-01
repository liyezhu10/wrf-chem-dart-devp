#!/bin/csh

set ppn = $1
set instance = $2

touch wrfHydroStillWorking.dum

mpirun -np $ppn ../../wrf_hydro.exe >& wrfHydroOutput.${instance} 

\rm wrfHydroStillWorking.dum

exit 0
