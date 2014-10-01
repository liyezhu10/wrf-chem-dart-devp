#!/bin/csh

# Because streamflow cannot simply be perturbed, 
# this script takes a perturbed inital state and runs if forward
# with(?) or without perturbed forcing to generate an 
# ensemble of streamflow at an initial time.

# input arguments endYyyy endMm endDD endHH inDir
# the inputDir should be of the form initialEnsemble.yyyymmdd.scheme 
# so that the start time can be grabbed from it.


if ($#argv != 5 & $#argv != 6) then
    echo "Usage: endYyyy endMm endDD endHH inDir(relative to DARTdir) ppn"
    exit 1
endif
set endYyyy = $1
set endMm   = $2
set endDd   = $3
set endHh   = $4
set inDir   = $5

echo endYyyy: $endYyyy
echo endMm: $endMm
echo endDd: $endDd
echo endHh: $endHh
echo inDir: $inDir
if ($#argv == 5) set ppn = $6
if(! $?ppn) set ppn = 1
if($ppn > 16) then 
    echo "Right now only configured to use a single node, setting ppn=1."
    set ppn = 1
endif 

set wd = `pwd`
echo wd: $wd

# parse the inDir for the date
set inDirDate = `echo $inDir | cut -d. -f3`
set startYyyy = `echo $inDirDate | cut -c1-4`
set startMm = `echo $inDirDate | cut -c5-6`
set startDd = `echo $inDirDate | cut -c7-8`

echo inDirDate: $inDirDate
echo startYyyy: $startYyyy
echo startMm: $startMm
echo startDd: $startDd

set startDateSec = `date -d "UTC $startYyyy-$startMm-$startDd" +%s`
set endDateSec = `date -d "UTC $endYyyy-$endMm-$endDd $endHh hour" +%s`
@ newKhour = ($endDateSec - $startDateSec) / 3600
echo $newKhour


set outPath = initialEnsembleWFlow.$startYyyy$startMm$startDd-$endYyyy$endMm$endDd.`echo $inDir | cut -d. -f3`

set num_states = `ls $inDir/restart.hydro.*.nc -1 | wc -l`


echo "Start Date: $startYyyy $startMm $startDd"
echo "End Date: $endYyyy $endMm $endDd"
echo "Number of hours to advance: $newKhour"
echo "Run dir: RUNS.$outPath"
echo "Input dir: $inDir"
echo "Number of ensemble members: $num_states"


# run all ensembles within this dir
\mkdir -p RUNS.$outPath || exit 1
cd RUNS.$outPath
ln -sf ../$inDir .  ## make this symlink so it's clear what these runs started from. 


## generate the random mean of each perturbation
#echo "../randomSeq.rsh uniform $num_states .5 1.5"
echo "../randomSeq.rsh uniform $num_states 1 1"
\cp ../randomSeq.rsh .
#set randomMeans = `./randomSeq.rsh uniform $num_states .5 2`
set randomMeans = `./randomSeq.rsh uniform $num_states 1 1`
echo The random mean precip multipliers:
echo $randomMeans

## get the name lists
\cp ../namelist.hrldas .
\cp ../hydro.namelist .

\ln -sf ../DOMAIN .

## deal with name list time. this requires a restrat file
\ln -sv $inDir/restart.0001.nc  restart.nc   || exit 2
set numadvances = $newKhour
# Have to match the start time in hrldas with the restart files,
# else the forcing data seems to have no effect.
set restartFileTime = `ncdump -v Times restart.nc | tail -2 | head -1 | cut -d'"' -f2`
set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
set restartFileMm = `echo $restartFileTime | cut -d- -f2`
set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`
 ## wont actually use this restart file
\rm restart.nc

ex namelist.hrldas <<ex_end
g;KHOUR ;s;= .*;= $numadvances;
g;START_YEAR;s;= .*;= $restartFileYyyy;
g;START_MONTH;s;= .*;= $restartFileMm;
g;START_DAY;s;= .*;= $restartFileDd;
g;START_HOUR;s;= .*;= $restartFileHh;
wq
ex_end


# Loop through each state /ensemble member
set state_copy = 1
while($state_copy <= $num_states)

    set instance        = `printf "%04d" $state_copy`
    echo "ensemble member: $instance"

    if ( `ls -1 ensemble.*/wrfHydroStillWorking.dum | wc -l` > 10 ) echo 'Waiting for processors...'
    while ( `ls -1 ensemble.*/wrfHydroStillWorking.dum | wc -l` > 10 )
	sleep 1
    end
   
    \mkdir ensemble.${instance}
    cd ensemble.${instance}

    ## get the name lists
    \cp ../namelist.hrldas .
    \cp ../hydro.namelist .
    
    # Get parameters.
    foreach FILE ( ../../PARAMS.gathered/* ) 
	echo $FILE
	\ln -sf $FILE . || exit 2  
    end
    # need the wrfHydro restart files for the individual ensemble members
    \ln -sv ../$inDir/restart.$instance.nc  restart.nc   || exit 2
    \ln -sv ../$inDir/restart.hydro.$instance.nc  restart.hydro.nc   || exit 2


    ## Perturb forcing 
    ## (0. outside this loop (state_copy) define a uniform random sequence of mean precip adjustment)
    ## 1. Edit the path to the forcing dir in the namelist.hrldas.
ex namelist.hrldas <<ex_end
g;INDIR ;s;= .*;= "FORCING.perturbed.spinup/";
wq
ex_end

    ## 2. for all forcing files in the time period
    @ firstForcingFileInd = `\ls -1 ../../FORCING/ | \grep -n $startYyyy$startMm${startDd}00 | \
                             cut -d: -f1`
    set forcingFiles = `\ls -1 ../../FORCING/ | \tail -n+$firstForcingFileInd | \head -$numadvances`
    ## 3. use the mean precip adjustment to define a sequence of normal errors about that mean
    ## through time. (use randomSeq.rsh, which computes text so you dont have to do it in shell)
    #set timePerts = `../../randomSeq.rsh uniform $numadvances \
    #      		"$randomMeans[$state_copy]*(1-.1)" \
    #                   "$randomMeans[$state_copy]*(1+.1)"`
    set timePerts = `../../randomSeq.rsh uniform $numadvances $randomMeans[$state_copy] \
			    $randomMeans[$state_copy]`


    @ countForc=1
    \mkdir FORCING.perturbed.spinup
    foreach fForcing ( $forcingFiles )
	##  multiply the precip in the copy using ncap2 and write a  
	##     copy to ../../FORCING.perturbed.spinup.iEns.endYyyyendMmendDd 
	set iPerturb = $timePerts[$countForc]

	ncap2 -s "RAINRATE=RAINRATE*${iPerturb}" ../../FORCING/$fForcing \
		 FORCING.perturbed.spinup/$fForcing

	@ countForc++
    end 

    cp ../../run_wrfHydro.csh .
    ./run_wrfHydro.csh $ppn $instance &

    ## clean up the forcing
#    \rm -rf FORCING.perturbed.spinup

    
    cd ../
    @ state_copy++
end


exit 0

