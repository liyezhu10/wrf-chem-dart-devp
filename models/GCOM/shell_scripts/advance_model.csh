#!/bin/tcsh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance 
#    and copies the necessary files into the temporary directory
# 2) converts the DART output to input expected by the ocean model
# 3) runs the ocean model
# 4) converts the ocean model output to input expected by DART
#
# The error code from the script reflects which block it failed.
#
# Arguments are the 
# 1) process number of caller, 
# 2) the number of state copies belonging to that process, and 
# 3) the name of the filter_control_file for that process

set process = $1
set num_states = $2
set control_file = $3

echo "My process ID is $process"
echo "I am responsible for advancing $num_states states."
echo "I am going to read my filenames from $control_file"

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to 
# run the ocean model.
#-------------------------------------------------------------------------

# Create a unique name for temporary working directory for this process.
# Create a clean temporary directory and go there.

set temp_dir = `printf "advance_temp%04d" $process`

\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.
foreach FILE ( ../input.nml )
   cp -pv $FILE . || exit 1
end

# Try to ensure that the input.nml has the required value for
# dart_to_gcom_nml:advance_time_present for this context.

echo '1'                      >! ex_commands
echo '/dart_to_gcom_nml'      >> ex_commands
echo '/advance_time_present'  >> ex_commands
echo ':s/\.false\./\.true\./' >> ex_commands
echo ':wq'                    >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

# copy the files used by GCOM & DART
foreach FILE ( ../gcom_geometry.nc )
   ln -sfv $FILE . || exit 1
end

echo 'listing now that the table has been set ...'
ls -l

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -n $ensemble_member_line ../$control_file | tail -n 1`
   set input_file      = `head -n $input_file_line      ../$control_file | tail -n 1`
   set output_file     = `head -n $output_file_line     ../$control_file | tail -n 1`

   # make sure we have a clean logfile for this entire advance
   set logfile = `printf "log_advance.%04d.txt"     $ensemble_member`

   echo "control_file is ../$control_file"             >! $logfile
   echo "working on ensemble_member $ensemble_member"  >> $logfile
   echo "input_file  is $input_file"                   >> $logfile
   echo "output_file is $output_file"                  >> $logfile

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by ocean model.
   # We are going to take a gcom netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED DART dart_to_gcom_input_file is 'dart_restart'
   # The dart_to_gcom_nml:advance_time_present = .TRUE. must be set
   # CENTRALDIR will always contain the gcom_restart.nc file.
   # The most recent gcom timestep is inserted into this file.

   set RESTARTFILE = `printf %04d gcom_restart_${ensemble_member}.nc'`

   ln -sfv ../$input_file    dart_restart >>& $logfile || exit 2
   cp -pv  ../${RESTARTFILE} gcom.r.nc    >>& $logfile || exit 2
   ../dart_to_gcom                        >>& $logfile || exit 2

   # Convey the new gcom 'advance_to' time to gcom via the namelist
   cat gcom_in.DART gcom_in.part2 >! gcom_in

   # update tiegcm param.dat values by grabbing the values from namelist_update  
   # and then overwriting whatever is in the param.dat
   # There is a danger that you match multiple things ... F107 and F107A,
   # SOURCE_START and START, for example ... so try to grep whitespace too ...

   set start_year   = `grep " START_YEAR "   namelist_update`
   set start_day    = `grep " START_DAY "    namelist_update`
   set source_start = `grep " SOURCE_START " namelist_update`
   set start        = `grep " START "        namelist_update`
   set secstart     = `grep " SECSTART "     namelist_update`
   set stop         = `grep " STOP "         namelist_update`
   set secstop      = `grep " SECSTOP "      namelist_update`
   set hist         = `grep " HIST "         namelist_update`
   set sechist      = `grep " SECHIST "      namelist_update`
   set save         = `grep " SAVE "         namelist_update`
   set secsave      = `grep " SECSAVE "      namelist_update`
   set f107         = `grep " F107 "         namelist_update`

   # FIXME SAVE and SECSAVE may not be needed ... they do not exist in 
   # the tiegcm.nml anymore. Is this true for all versions in use ...
   #
   # the way to think about the following sed syntax is this:
   # / SearchStringWithWhiteSpaceToMakeUnique  /c\ the_new_contents_of_the_line 

   sed -e "/ START_YEAR /c\ ${start_year}" \
       -e "/ START_DAY /c\ ${start_day}" \
       -e "/ SOURCE_START /c\ ${source_start}" \
       -e "/ START /c\ ${start}" \
       -e "/ STOP /c\ ${stop}" \
       -e "/ HIST /c\ ${hist}" \
       -e "/ SECSTART /c\ ${secstart}" \
       -e "/ SECSTOP /c\ ${secstop}" \
       -e "/ SECHIST /c\ ${sechist}" \
       -e "/ F107 /c\ ${f107}" \
       -e "/ SOURCE /c\ SOURCE = 'tiegcm_restart_p.nc'" \
       -e "/ OUTPUT /c\ OUTPUT = 'tiegcm_restart_p.nc'" \
       -e "/ SECOUT /c\ SECOUT = 'tiegcm_s.nc'"         \
       tiegcm.nml >! tiegcm.nml.updated

   if ( -e tiegcm.nml.updated ) then
      echo "tiegcm.nml updated with new start/stop time for ensemble member $ensemble_member"
      mv tiegcm.nml tiegcm.nml.original
      mv tiegcm.nml.updated tiegcm.nml
   else
      echo "ERROR tiegcm.nml did not update correctly for ensemble member $ensemble_member."
      exit 1
   endif

   #----------------------------------------------------------------------
   # Block 3: Run the ocean model
   # The CCSM version has a pointer file that contains the name of the
   # last restart. The LANL version has no such mechanism, but the 
   # filename can be predicted from the gcom_in namelist information.
   #----------------------------------------------------------------------
   # the value of MPI is inherited

   rm -f ocn.log.*

   ${MPI} ../gcom || exit 3

   grep "Successful completion of gcom run" ocn.log.*
   set gcomstatus = $status
   if ( $gcomstatus != 0 ) then
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully" 
      echo "ERROR - gcom ensemble member $ensemble_member did not complete successfully" 
      exit 3 
   endif
   
   #----------------------------------------------------------------------
   # Block 4: Convert the ocean model output to form needed by DART
   #----------------------------------------------------------------------

   ls -lrt

   # gcom makes a new restart file and updates the pointer file
   # Rename the gcom pointer file contents to contain the ensemble member info

   set RESTARTFILE = `head -1 rpointer.ocn.restart`
   set NEWFILE = `echo $RESTARTFILE | sed -e "s/gcom/gcom.$ensemble_member/"`
   echo "gcom member $ensemble_member made restart file $NEWFILE"

   mv -v $RESTARTFILE $NEWFILE
  
   echo $NEWFILE         >! rpointer.ocn.restart 
   echo "RESTART_FMT=nc" >> rpointer.ocn.restart

   ln -svf ${NEWFILE} gcom.r.nc || exit 4
   
   # gcom_to_dart reads the restart file after the model advance and writes
   # out an updated DART 'initial conditions' file. This initial conditions
   # file contains a header with the valid time of the ensuing model state.
   # The gcom restart files contain the valid time of the model state.

   ../gcom_to_dart || exit 4

   # The (new,updated) DART restart file name is called 'dart_ics'
   # Move the updated files back to 'centraldir'
   mv -v dart_ics ../$output_file || exit 4
   mv -v rpointer.ocn.restart ../rpointer.ocn.${ensemble_member}.restart || exit 4
   mv -v ${NEWFILE} ../${NEWFILE} || exit 4

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# must communicate the time_manager_nml:stop_count 
# cp -pv gcom_in.DART ../gcom_in

# Change back to original directory and get rid of temporary directory
cd ..
# \rm -rf $temp_dir

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
\rm -rf $control_file

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

