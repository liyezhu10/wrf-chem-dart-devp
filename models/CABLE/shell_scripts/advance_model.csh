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
# 2) converts the DART output to input expected by the CABLE model
# 3) runs CABLE 
# 4) converts the CABLE model output to input expected by DART
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

echo "process      is $process"
echo "num_states   is $num_states"
echo "control_file is $control_file"

#-------------------------------------------------------------------------
# Block 1: populate a run-time directory with the bits needed to run CABLE
#-------------------------------------------------------------------------

# Get unique name for temporary working directory for this process's stuff
set temp_dir = 'advance_temp'${process}
echo "temp_dir is ${temp_dir}"

# Create a clean temporary directory and go there
\rm -rf  ${temp_dir}
mkdir -p ${temp_dir}
cd       ${temp_dir}

# Get the 'changing' namelist files from CENTRALDIR
# Only the namelists in CENTRALDIR have the updated information about
# the state of the model partway through an assimilation experiment.

cp -pv ../cable.nml . || exit 1
cp -pv ../input.nml . || exit 1

# Try to ensure that the input.nml has the required value for
# dart_to_cable_nml:advance_time_present for this context.

echo '1'                      >! ex_commands
echo '/dart_to_cable_nml'     >> ex_commands
echo '/advance_time_present'  >> ex_commands
echo ':s/\.false\./\.true\./' >> ex_commands
echo ':wq'                    >> ex_commands

( ex input.nml < ex_commands ) >& /dev/null
\rm -f ex_commands

# Loop through each state
set           state_copy = 1
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3
while($state_copy <= $num_states)
   
   set ensemble_member = `head -n $ensemble_member_line ../${control_file} | tail -n 1`
   set input_file      = `head -n $input_file_line      ../${control_file} | tail -n 1`
   set output_file     = `head -n $output_file_line     ../${control_file} | tail -n 1`

   echo "advancing ensemble member $ensemble_member ..."

   #----------------------------------------------------------------------
   # Block 2: Convert the DART output file to form needed by CABLE.
   # We are going to take a CABLE netCDF restart file and simply overwrite the
   # appropriate variables. The DART output file also has the 'advance_to'
   # time - which must be communicated to the model ...
   #----------------------------------------------------------------------

   # The EXPECTED DART filter restart_out_file_name is 'dart_restart'
   # The dart_to_cable_nml:advance_time_present = .TRUE. must be set;
   # this is the purpose of the ex_commands earlier.

   ln -sfv ../${input_file} dart_restart || exit 2

   set CABLEFILE = `printf restart_in_gpcc.%04d.nc ${ensemble_member}`
   ln -sfv ../${CABLEFILE} restart_in_gpcc.nc || exit 2
   
   ../dart_to_cable || exit 2

   # Extract just the forcing for the model advance we need.
   # dart_to_cable produces a file "time_control.txt" that has
   # the forcing times we need to advance CABLE. We use them
   # to extract what we need with ncks.
   #
   #cat time_control.txt 
   #     31622400
   #     31708800
   # dart_to_cable:DART   model date 1981 Jan 01 00:00:00
   # dart_to_cable:DART desired date 1981 Jan 02 00:00:00

   set t1 = `head -n 1 time_control.txt`
   set tN = `head -n 2 time_control.txt | tail -n 1`

   # gswpfile%rainf = '/short/xa5/CABLE-AUX/GPCC-CABLE/prcp_hr_1980-19801x1.nc'
   # gswpfile%LWdown= '/short/xa5/CABLE-AUX/GPCC-CABLE/dlwrf_hr_1980-19801x1.nc'
   # gswpfile%SWdown= '/short/xa5/CABLE-AUX/GPCC-CABLE/dswrf_hr_1980-19801x1.nc'
   # gswpfile%PSurf = '/short/xa5/CABLE-AUX/GPCC-CABLE/pres_hr_1980-19801x1.nc'
   # gswpfile%Qair  = '/short/xa5/CABLE-AUX/GPCC-CABLE/shum_hr_1980-19801x1.nc'
   # gswpfile%Tair  = '/short/xa5/CABLE-AUX/GPCC-CABLE/tas_hr_1980-19801x1.nc'
   # gswpfile%wind  = '/short/xa5/CABLE-AUX/GPCC-CABLE/wind_hr_1980-19801x1.nc'

   foreach FILE ( rainf LWdown SWdown PSurf Qair Tair wind )

      set FILESTRING=`grep "gswpfile%$FILE" cable.nml | sed -e "s#[='\\!]# #g"`
      set FILENAME=$FILESTRING[$#FILESTRING]
      set FILETAIL=$FILENAME:t

      # This actually subsets the files and makes local files.

      ncks -O -d time,${t1}.,${tN}. ${FILENAME} subset_${FILETAIL}

      # This changes the namelist to use the local files.
      
      grep -v "gswpfile%$FILE" cable.nml | grep -v '&end' >! cabletemp.nml
      echo "   gswpfile%$FILE = '"subset_${FILETAIL}"'"   >> cabletemp.nml
      echo '&end'                                         >> cabletemp.nml
   
      mv cabletemp.nml cable.nml   

   end

   #----------------------------------------------------------------------
   # Block 3: Run CABLE
   # casafile%cnpipool    = 'output/cnppool1979.csv'
   # filename%restart_in  = './restart_in_gpcc.nc'
   # filename%restart_out = './restart_out.nc'
   #----------------------------------------------------------------------

   mkdir -p output

   set CNIPOOLSTRING=`grep "casafile%cnpipool" cable.nml | sed -e "s#[='\\!]# #g"`
   set CNIPOOLFILE=$CNIPOOLSTRING[$#CNIPOOLSTRING]
  
   ln -sfv /short/xa5/CABLE-EXE/$CNIPOOLFILE output/.

   echo "mpi will run with the following command <${MPICMD}>"

   ${MPICMD} ../cable-mpi

   set cablestatus = $status
   if ( $cablestatus != 0 ) then
      echo "ERROR - ensemble member $ensemble_member did not complete successfully" 
      echo "ERROR - ensemble member $ensemble_member did not complete successfully" 
      exit 3 
   endif

   # FIXME if cable completes correctly, the carbon pool file must be used.
   # this means the casafile%cnpipool variable must point to the new file.
   
   #----------------------------------------------------------------------
   # Block 4: Convert the CABLE model output to form needed by DART
   #----------------------------------------------------------------------

   ls -lrt

   mv -v restart_out.nc restart_in_gpcc.nc

   # cable_to_dart reads the restart file after the model advance and writes
   # out an updated DART 'initial conditions' file. This initial conditions
   # file contains a header with the valid time of the ensuing model state.
   # The cable restart files contain the valid time of the model state.

   ../cable_to_dart || exit 4

   # The (new,updated) DART restart file name is called 'dart_ics'
   # Move the updated files back to 'centraldir'
   mv -v dart_ics            ../${output_file} || exit 4
   mv -v restart_in_gpcc.nc  ../${CABLEFILE}   || exit 4

   # bookkeeping

   @ state_copy++
   @ ensemble_member_line = $ensemble_member_line + 3
   @ input_file_line = $input_file_line + 3
   @ output_file_line = $output_file_line + 3
end

# must communicate the time_manager_nml:stop_count 
# cp -pv cablegc_in.DART ../cable_in

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

