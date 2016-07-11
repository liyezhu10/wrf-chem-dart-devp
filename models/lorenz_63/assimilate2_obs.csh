#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: assimilate1_5.csh 10432 2016-06-29 17:48:31Z thoar $

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN LORENZ_63 ASSIMILATION"

# directory that contains the filter
set RUNDIR = "/Users/hendric/DART/pda/models/lorenz_63/work"

#-------------------------------------------------------------------------
# Get the case-specific variables
#-------------------------------------------------------------------------

cd ${RUNDIR}

set n = 1

set restart_file = 'filter_restart'
set mean_file = 'mean.nc'

foreach OBS_FILE (obs_seq.*.out)

   set TIME = `printf %04d ${n}`   

   # create directory to copy restart and inflation files
   # for each time step
   set ADV_DIR = "advance_time_$TIME"
   mkdir $ADV_DIR

   if (  -e   ${OBS_FILE} ) then
      ln -sf ${OBS_FILE} obs_seq.out
   else
      echo "ERROR ... no observation file ${OBS_FILE}"
      echo "ERROR ... no observation file ${OBS_FILE}"
      exit -1
   endif
   
   ./filter
   
   cp $restart_file*.nc $ADV_DIR
   cp $mean_file        $ADV_DIR
   cp obs_seq.final     $ADV_DIR
   
   @ n++

end

exit(0)

