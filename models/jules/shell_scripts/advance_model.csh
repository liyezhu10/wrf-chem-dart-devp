#!/bin/csh
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
# 2) converts the DART output to input expected by the model
# 3) runs the model
# 4) converts the model output to input expected by DART
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

# These are information only, the second two can be deleted any time.
echo "Starting advance_model for process $process at "`date`
echo "Responsible for advancing $num_states "
echo "Using a control file named $control_file"

#--------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run JULES.
#--------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = `printf "advance_temp%04d" $process`

# Create a clean temporary directory and go there
\rm -rf  $temp_dir
mkdir -p $temp_dir
cd       $temp_dir

# Get the REQUIRED namelist files to the temp directory
# SR: Copying all the namelist files of JULES to be on
# the safe side, since some of them will be changed by DART.
# The namelist files have the exact names
# of the directories holding the input files. Therefore,
# I do not think that it is necessary to copy the input files.

\cp ../*nml . || exit 1

# Ensure that the input.nml has the required value for
# dart_to_model_nml:advance_time_present for this context.

sed -e "/advance_time_present /c\ advance_time_present = .TRUE." \
       ../input.nml >! input.nml || exit 1

sed -e "/^file/c\ file = 'jules_restart.nc'" \
       ../initial_conditions.nml >! initial_conditions.nml || exit 1

sed -e "/^file/c\ file = '../tile_fractions.dat'" \
       ../ancillaries.nml >! ancillaries.nml || exit 1

sed -e "/^output_dir/c\output_dir = './'" \
       ../output.nml >! output.nml || exit 1

# Check to see if you are running async==2 
# single-threaded JULES executable being advanced without MPI
set MYSTRING = `grep -A 42 filter_nml input.nml | grep async`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set ASYNC = $MYSTRING[2]

if ($ASYNC == 2) then
   set RUN_CMD = ''
else
   echo "ERROR: JULES must be run with async==2"
   echo "ERROR: JULES must be run with async==2"
   exit 2
endif

# Loop through each state
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2
set output_file_line = 3

while($state_copy <= $num_states)

   set ensemble_member = `head -n $ensemble_member_line ../$control_file | tail -n 1`
   set input_file      = `head -n $input_file_line      ../$control_file | tail -n 1`
   set output_file     = `head -n $output_file_line     ../$control_file | tail -n 1`
   set instance        = `printf "%04d" $ensemble_member`

   # make sure we have a clean logfile for this entire advance
   set logfile = `printf "log_advance.%04d.txt"     $ensemble_member`

   echo "control_file is ../$control_file"             >> $logfile
   echo "working on ensemble_member $ensemble_member"  >> $logfile
   echo "input_file  is $input_file"                   >> $logfile
   echo "output_file is $output_file"                  >> $logfile
   echo "Starting dart_to_jules at "`date`             >> $logfile

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by model.
   # Overwrite the appropriate variables of a JULES netCDF restart file.
   # The DART output file (namelist_update) has the 'advance_to' time 
   # which must be communicated to the model ... through the namelist
   #----------------------------------------------------------------------
  
   # The last output file produced by JULES has to be copied as jules_output.nc
   set jls_out = `ls -1  ../*hour*nc | tail -n 1`
   ln -sf $jls_out                 jules_output.nc     >>& $logfile || exit 2
   set jls_rst = `ls -1 ../*dump*nc | tail -n 1`
   ln -sf $jls_rst                 jules_restart.nc    >>& $logfile || exit 2
   
   cp -pv ../$input_file           dart_restart        >>& $logfile || exit 2
   #cp -pv ../restart.$instance.nc  jules_restart.nc    >>& $logfile || exit 2
   #ln -sf ../output.$instance.nc   jules_output.nc     >>& $logfile || exit 2

   # All the required files are created. Now run dart_to_jules
   # This creates a file called DART_time_control.txt that has
   # the information needed to update the JULES namelist to tell
   # JULES how far to advance.
   
   ../dart_to_jules                                    >>& $logfile || exit 2  

   set start_string = `grep "main_run_start"  DART_time_control.txt`
   set  stop_string = `grep "main_run_end"    DART_time_control.txt`

   sed -e "/main_run_start/c\${start_string}" \
         -e "/main_run_end/c\${stop_string}" timesteps.nml >! timesteps.nml.update

   if ( -e  timesteps.nml.update ) then
      echo "timesteps.nml update with new start/stop time for ensemble member $ensemble_member" >> $logfile
      mv -v timesteps.nml timesteps.nml.original
      mv -v timesteps.nml.update timesteps.nml
   else
      echo "ERROR timesteps.nml did not update correctly for ensemble member $ensemble_member." >> $logfile
      exit 2
   endif
   #----------------------------------------------------------------------
   # Block 3: Run the model
   #----------------------------------------------------------------------
   # SR: We are going to run the model in this block. First we need to know 
   # how long we need to run. This information is generated by dart_to_jules
   # and stored in the work directory.     
   ../jules.exe					       >>& $logfile || exit 3
   



   foreach FILE (*.dump.*)
           set BASE=$FILE:r
           set ROOT=$BASE:r
	   set DATESTR=$ROOT:e
	   set SECONDS=$BASE:e
           set MYSTRING=`printf "%05d" $SECONDS`
           mv $FILE $ROOT.$MYSTRING.nc
   end
   
   foreach FILE (*.hour.*)
           set BASE=$FILE:r
           set SECONDS=$BASE:e
           set MYSTRING=`printf "%05d" $SECONDS`
           mv $FILE $BASE.$DATESTR.$MYSTRING.nc
   end

   #----------------------------------------------------------------------
   # Block 4: Convert the model output to form needed by DART
   # At this point, the model has updated the information in jules_restart.nc
   # We need to get that information back into the DART state vector.
   #
   # The updated information needs to be moved into CENTRALDIR in
   # preparation for the next cycle.
   #----------------------------------------------------------------------
  
   # The last output file produced by JULES has to be copied as jules_output.nc
   set jls_out = `ls -1 *hour*nc | tail -n 1`
   ln -sf $jls_out                 jules_output.nc     >>& $logfile || exit 4
   set jls_rst = `ls -1 *dump*nc | tail -n 1`
   ln -sf $jls_rst                 jules_restart.nc    >>& $logfile || exit 4

 
   ../jules_to_dart                                    >>& $logfile || exit 4

   mv *hour*nc ../
   mv *dump*nc ../

   mv -v dart_ics            ../$output_file || exit 4   
 
   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line      = $input_file_line + 3
   @ output_file_line     = $output_file_line + 3
   echo "Finished model_to_dart at "`date`  >> $logfile

end

# Change back to original directory 
cd ..

# After you are assured this script works as expected, you can actually 
# remove the temporary directory. For now ... leave this commented OUT.
#\rm -rf $temp_dir

# Remove the filter_control file to signal completion
\rm -fv $control_file

echo "Finished advance_model for process $process at "`date`

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

