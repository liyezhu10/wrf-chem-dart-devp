#!/bin/ksh -x
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
   export WRFINP=wrfinput_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e${IENS}
   export WRFBDY=wrfbdy_d01_${YEAR}-${MONTH}-${DAY}_${HOUR}:00:00.e${IENS}
   echo 'get inp'
   cp ${IC_METDIR}/${WRFINP} ./.
   cp ${BC_METDIR}/${WRFBDY} ./.
   ls set${i}
   rm -f mozbc.ic.inp.set${i}
   cat mozbc.ic.inp set${i} > mozbc.ic.inp.set${i}
   rm -f mozbc.bc.inp.set${i}
   cat mozbc.bc.inp set${i} > mozbc.bc.inp.set${i}
   echo  'run mozbc'
   ./run_mozbc_rt_CR.csh type=ic mozbc_inp=mozbc.ic.inp.set${i} ens=${IENS}
   ./run_mozbc_rt_CR.csh type=bc mozbc_inp=mozbc.bc.inp.set${i} ens=${IENS}
   echo 'put files'
   echo 'OK'
done
