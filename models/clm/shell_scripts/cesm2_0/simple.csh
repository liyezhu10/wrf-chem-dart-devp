#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This file is just to see if CESM can be bullt on whatever architecture.
# Single instance ... nothing fancy ... no DART ... just something simple.

/glade/p/work/thoar/CESM/clm_dev_branch/cime/scripts/query_config  --compsets clm

set CASE = pmo_2000
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
    --run-unsupported || exit 1

#    --ninst 2 \
#    --multi-driver \

cd ${CASEROOT}

echo "TJH: Finished create_newcase ..."
echo "TJH: Starting case.setup     ..."

./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=2000-01-01
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
echo ""
echo "  cd  ${CASEROOT}"
echo "  ./case.submit"
echo ""
echo "  if that works, try a do_nothing assimilation cycle."

# These modify env_run.xml so they can be executed AFTER a successful clm cycle.
if ( 1 == 2 ) then
   ./xmlchange DATA_ASSIMILATION=TRUE
   ./xmlchange DATA_ASSIMILATION_CYCLES=1
   ./xmlchange DATA_ASSIMILATION_SCRIPT=do_nothing.csh
   echo "echo Hello from `hostname`" >! do_nothing.csh
   chmod 755 do_nothing.csh
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
