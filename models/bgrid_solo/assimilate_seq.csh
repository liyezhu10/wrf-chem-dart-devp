#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

if ( $#argv > 0 ) then
   set end_n = $1
else
   set end_n = 2
endif

echo "`date` -- BEGIN B-GRID ASSIMILATION"

# directory that contains the filter
set DARTDIR = /Users/thoar/svn/DART/pda/models/bgrid_solo/work/
set EXPERIMENT = ${DARTDIR}/test_1

#-------------------------------------------------------------------------
# Get the case-specific variables
#-------------------------------------------------------------------------

mkdir -p ${EXPERIMENT}

cd ${EXPERIMENT}

set n = 1

set restart_file   = 'filter_restart'
set mean_file      = 'mean.nc'
set prior_inf_sd   = 'prior_inflate_restart_sd.nc'
set prior_inf_mean = 'prior_inflate_restart_mean.nc'

ls -1 ${DARTDIR}/filter_ics.00* > restart_file_list.txt

cp ${DARTDIR}/input.nml.cycle input.nml
cp ${DARTDIR}/bgrid.nc        .
cp ${DARTDIR}/filter          .

foreach OBS_FILE ( ${DARTDIR}/obs/obs_seq.*.out )

   set TIME = `printf %04d ${n}`   

   # create directory to copy restart and inflation files
   # for each time step
   set ADV_DIR = "advance_time_$TIME"
   mkdir $ADV_DIR

   if (  -e  ${OBS_FILE} ) then
      ln -sf ${OBS_FILE} obs_seq.out
   else
      echo "ERROR ... no observation file ${OBS_FILE}"
      exit -1
   endif
   
   ./filter || exit 1
   
   mv $restart_file*    $ADV_DIR
   mv $mean_file        $ADV_DIR
   mv obs_seq.final     $ADV_DIR
   mv Prior_Diag.nc     $ADV_DIR
   mv Posterior_Diag.nc $ADV_DIR
   if ( -e sd.nc ) then
      mv   sd.nc        $ADV_DIR
   endif

   # copy ouput inflation from filter to inital inflation file
   # for next run.
   if ( -e $prior_inf_sd ) then
      cp   $prior_inf_sd   prior_inflate_ics_sd
      mv   $prior_inf_sd  $ADV_DIR
   endif

   if ( -e $prior_inf_mean ) then
      cp   $prior_inf_mean prior_inflate_ics_mean
      mv   $prior_inf_mean $ADV_DIR
   endif

   ls -1 $ADV_DIR/$restart_file*.nc > restart_file_list.txt

   if ("$n" == "$end_n") then
      break
   endif

   @ n++

end

# todo FIXME : do we need to keep these around or save them to the ADV_DIR?
rm prior_member*


exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

