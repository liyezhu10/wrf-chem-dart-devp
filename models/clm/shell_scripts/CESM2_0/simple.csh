#!/bin/csh

/glade/p/work/thoar/CESM/clm_dev_branch/cime/scripts/query_config  --compsets clm

set CASE = clm50bgc_test_3
set CASEDIR = /glade/p/work/thoar/cases
set COMPSET = I2000Clm45Sp
set COMPSET = I2000Clm50BgcCrop
set COMPSET = 2000_DATM%GSWP3v1_CLM50%BGC-CROP_SICE_SOCN_MOSART_SGLC_SWAV
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
    --ninst 3 \
    --run-unsupported || exit 1


cd ${CASEROOT}

echo "TJH: Finished create_newcase ..."
echo "TJH: Starting case.setup     ..."

./case.setup || exit 2

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

