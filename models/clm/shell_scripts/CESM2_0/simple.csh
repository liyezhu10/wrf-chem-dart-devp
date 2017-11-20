#!/bin/csh

/glade/p/work/thoar/CESM/clm_dev_branch/cime/scripts/query_config  --compsets clm

set CASE = pmo
set CASEDIR = /glade/p/work/thoar/cases
set COMPSET = I2000Clm50BgcCrop
set COMPSET = 2000_DATM%GSWP3v1_CLM50%BGC-CROP_SICE_SOCN_RTM_SGLC_SWAV
set COMPSET = I2000Clm45Sp
set COMPSET = 2000_DATM%GSWP3v1_CLM45%SP_SICE_SOCN_RTM_SGLC_SWAV 
set CASEROOT = ${CASEDIR}/${CASE}

\rm -rf ${CASEROOT}
\rm -rf /glade/scratch/thoar/${CASE}/bld
\rm -rf /glade/scratch/thoar/${CASE}/run

echo "TJH: Starting create_newcase ..."

/glade/p/work/thoar/CESM/clm_dev_branch/cime/scripts/create_newcase \
    --case ${CASEROOT} \
    --compset ${COMPSET} \
    --mach cheyenne \
    --res f09_g16 \
    --project P86850054 \
    --run-unsupported || exit 1

    --ninst 3 \

cd ${CASEROOT}

echo "TJH: Finished create_newcase ..."
echo "TJH: Starting case.setup     ..."

./xmlchange RUN_STARTDATE=2000-07-01
./xmlchange DOUT_S=FALSE
./xmlchange CALENDAR=GREGORIAN

./case.setup || exit 2

# These seemed to have no effect when run before case.setup
./xmlchange --subgroup case.run --id JOB_QUEUE          --val economy
./xmlchange --subgroup case.run --id JOB_WALLCLOCK_TIME --val 0:20

echo "TJH: Finished case.setup     ..."
echo "TJH: Starting case.build     ..."

./case.build || exit 3

echo "TJH: Finished case.build     ..."

exit 0


# These modify the env_run.xml so they should be 
./xmlchange DATA_ASSIMILATION=TRUE
./xmlchange DATA_ASSIMILATION_CYCLES=1
./xmlchange DATA_ASSIMILATION_SCRIPT=do_nothing.csh
echo "echo Hello World" >! do_nothing.csh
chmod 755 do_nothing.csh

