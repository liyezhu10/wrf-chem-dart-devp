#!/bin/csh -ef

#-- define the location of the emsemble runs ($RUNPATH/$RUNPREFIX_001, ...)
set RUNPATH=/mnt/lustre/lus0/home/dart/dart/rma_openggcm/models/openggcm/test
set RUNPREFIX=target

#-- define frequency of openggcm call
set FREQ=600

#-- determine ensemble size
set inst_num=1
set inst=`printf "%03d" $inst_num`
set run=$RUNPREFIX"_"$inst
set rundir=$RUNPATH/$run
while ( -d $rundir )
  @ inst_num++
  set inst=`printf "%03d" $inst_num`
  set run=$RUNPREFIX"_"$inst
  set rundir=$RUNPATH/$run
end
@ inst_num--
set ensemble_size=$inst_num
if ( $ensemble_size == 0 ) then
  echo "no "$RUNPREFIX"_XXX directories found in "$RUNPATH
  exit
endif
echo "ensemble size="$ensemble_size

#-- give some time for runs to be entered in queue
sleep 10

#-- loop until all runs leave queue 
#-- (be careful, $RUNPREFIX could exist in other runs)
set testrun=$RUNPREFIX"_"
set inqueue=`qstat -a | grep $testrun | wc -l`
while ( $inqueue != 0 )

  #-- look for semaphores, run filter when directed
  foreach inst_num (`seq $ensemble_size`)
    set inst=`printf "%03d" $inst_num`
    set run=$RUNPREFIX"_"$inst
    set rundir=$RUNPATH/$run
    if ( -r $rundir/target/to_dart.semaphore ) then
      echo "Running DART" `date` $run
      cp -v $rundir/target/DATA.ionos2.nc $rundir/target/dart_posterior.nc
      echo $FREQ > $rundir/target/from_dart.semaphore
      cat $rundir/target/from_dart.semaphore
      rm $rundir/target/to_dart.semaphore
    endif
  end

  sleep 1

  set inqueue=`qstat -a | grep $testrun | wc -l`
end
