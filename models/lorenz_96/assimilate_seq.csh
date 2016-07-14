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

if ( $#argv > 0 ) then
   set end_n = $1
else
   set end_n = 3
endif

echo "`date` -- BEGIN LORENZ_96 ASSIMILATION"

# directory that contains the filter
set RUNDIR = "/Users/hendric/DART/pda/models/lorenz_96/work/"

#-------------------------------------------------------------------------
# Get the case-specific variables
#-------------------------------------------------------------------------

cd ${RUNDIR}

set n = 1

set restart_file   = 'filter_restart'
set mean_file      = 'mean.nc'
set prior_inf_sd   = 'prior_inflate_restart_sd.nc'
set prior_inf_mean = 'prior_inflate_restart_mean.nc'

# time stamp set to zero for sequential time steping
cp filter_ics_t0 filter_ics

foreach OBS_FILE (../obs/obs_seq.*.out)

   if ("$n" == "1") then
      cp input.first_step.nml  input.nml 
   else
      cp input.from_netcdf.nml input.nml
   endif

   set TIME = `printf %04d ${n}`   

   # create directory to copy restart and inflation files
   # for each time step
   set ADV_DIR = "advance_time_$TIME"
   mkdir $ADV_DIR

   if (  -e   ${OBS_FILE} ) then
      ln -sf ${OBS_FILE} obs_seq.out
   else
      echo "ERROR ... no observation file ${OBS_FILE}"
      exit -1
   endif
   
   ./filter
   
   mv $restart_file*    $ADV_DIR
   mv $mean_file        $ADV_DIR
   mv obs_seq.final     $ADV_DIR
   # mv $prior_inf_sd     $ADV_DIR
   # mv $prior_inf_mean   $ADV_DIR
   
   # copy ouput inflation from filter to inital inflation file
   # for next run.
   cp $prior_inf_sd   prior_inflate_ics_sd
   cp $prior_inf_mean prior_inflate_ics_mean

   ls -1 $ADV_DIR/$restart_file*.nc > restart_file_list.txt

   if ("$n" == "$end_n") then
      break
   endif

   @ n++

end

# todo FIXME : do we need to keep these around or save them to the ADV_DIR?
rm prior_member*


exit(0)

