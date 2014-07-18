#!/bin/tcsh -f
#
# Download the AMSRE data with:
#
# wget -t 0 -T 120 -c -r -nH --cut-dirs=4 \
#    ftp://sidads.colorado.edu/pub/DATASETS/nsidc0301_amsre_ease_grid_tbs/global/2003/

# Directories for Ally
# set RAW_OBS_DIR = /gpfsm/dnb52/projects/p31/DATA/DART_DATA/Observations/AMSRE_obs/AMSRE_RAW
# set AMSRE_WORKDIR = /gpfsm/dnb31/atoure/DART/observations/AMSR-E

# Directories for Tim
set   RAW_OBS_DIR = /glade/scratch/thoar/AMSR-E
set AMSRE_WORKDIR = /glade/scratch/thoar/AMSR-E/work
set      DART_DIR = /glade/p/work/thoar/DART/Tb/observations/AMSR-E/work

set nonomatch

# If the work directory does not exist, make it.
if  ( ! -d ${AMSRE_WORKDIR} ) then
  mkdir -p ${AMSRE_WORKDIR}
endif

cd ${AMSRE_WORKDIR}

\cp -v ${DART_DIR}/input.nml        .
\cp -v ${DART_DIR}/ease_grid_to_obs .

# Set DOY equal to the first day-of-year of interest.
# Set NDAYS equal to the number of days to convert.
# Keep in mind it takes a LONG time to convert 1 day of
# all frequencies, polarizations, passes. For the NH alone,
# this takes more than 40 minutes and results in 6.3+ MILLION
# observations - for a single day! 

@ YEAR = 2003
@ DOY = 120
@ NDAYS = 1

# do not change iday.

@ IDAY = 1

while ($IDAY <= $NDAYS)

   set FileBaseDate = `printf ID2r1-AMSRE-ML%04d%03d $YEAR $DOY`

   echo "looking for ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}*[HV]"
   \ls -1 ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}*[HV] >! file_list.txt

   if ( ! -z file_list.txt ) then

       echo "Input files are:"
       cat file_list.txt

       #The format of name of Date_output_file the has to be:yyyy-mm-dd-00000
       # BUT the actual filename date must be THE NEXT DAY for CLM/DART
       # i.e. the data for April 22nd must be in a file named April 23rd
       @ tomorrow = $DOY + 1
       set FILEDATE = `date -d "$YEAR-01-01 +$tomorrow days -1 day" "+%F"`
       echo "FILEDATE is $FILEDATE"
       set Date_output_file = ${FILEDATE}-00000

       # specify sub-dir of processed data,
       # data in the same month are stored in "/YYYYMM/"
       #Get month and day using $DOY and $YEAR
       set YYYYMM=`date -d "$YEAR-01-01 +$DOY days -1 day" "+%Y%m"`
       echo "YYYYMM is $YYYYMM"

       set Output_Obs_dir = "DART_OBS_SEQ/${YYYYMM}"
       echo "The processed data dir is " $Output_Obs_dir

       if  ( ! -d ${RAW_OBS_DIR}/${Output_Obs_dir} ) then
         mkdir -p ${RAW_OBS_DIR}/${Output_Obs_dir}
       endif

       # specify output Output_fileName, format: obs_seq.yyyy-mm-dd-00000.out
       set Output_fileName = "obs_seq.${Date_output_file}.out"
       set Full_Output_fileName = "${RAW_OBS_DIR}/${Output_Obs_dir}/${Output_fileName}"
       echo "Full_Output_fileName = $Full_Output_fileName"

       echo "Starting ease_grid_to_obs for date $Date_output_file"

       # Convert this days worth of obs to a single file 'obs_seq.out'
       # must rename this file to have the date expected by CLM/DART
       ./ease_grid_to_obs

       # link the new observation sequence file to the filename
       # expected by "ease_grid_to_obs"
       \mv -v  obs_seq.out ${Full_Output_fileName}

       # move to next loop for the following day
       echo "Finished ease_grid_to_obs for date $Date_output_file"
   else
      echo "WARNING : No observation files for $YEAR $DOY ... AKA ..."
      echo "WARNING : ${RAW_OBS_DIR}/${YEAR}/${FileBaseDate}*[HV]"
   endif

   @ DOY ++
   @ IDAY ++

end

echo "Finish all data processing"

exit 0
