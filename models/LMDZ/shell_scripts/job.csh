#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#****************************************************************************
# Tarkeshwar Singh
# PhD, IIT Delhi
# Email: tarkphysics87@gmail.com
#
# Purpose:  Copy all dart ics, submit filter run, store restart and outputs files
#****************************************************************************

#************* PBS SETTING **************************************************
### Set the job name
#PBS -N dart_assim 
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
#PBS -M singh@cas.iitd.ac.in
#PBS -P cas
#PBS -l select=65:ncpus=9:mpiprocs=1
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=100:00:00
#### Get environment variables from submitting shell
#PBS -V
## After job starts, must goto working directory. 
## $PBS_O_WORKDIR is the directory from where the job is fired. 

cd $PBS_O_WORKDIR

set num_proc = 120

# Always assign one extra node and save it in hostfile_lastFreeNode
# Required in run_lmdz.csh script
cat $PBS_NODEFILE  >PBS_NODEFILE
cat PBS_NODEFILE | tail -1 > hostfile_lastFreeNode

echo $PBS_JOBID

#****************************************************************************

# Set alias
# #
  if ( ! $?REMOVE ) then
       set REMOVE = 'rm -rf'
  endif
  if ( ! $?COPY ) then
       set COPY = 'cp -fp'
  endif
  if ( ! $?MOVE ) then
       set MOVE = 'mv -f'
  endif
  if ( ! $?LINK ) then
       set LINK = 'ln -fs'
  endif
# remove any pre existing temporary filles and directories
$REMOVE start_* startphy_*  start.nc assim_model_state_* lmdz_out_temp* advance_temp*
#--------------------------------------------------------------------------------------
source Control_File.csh
 #--------------
 # Determine the number of ensemble members from input.nml,
set ENSEMBLESTRING = `grep -A 42 filter_nml input.nml | grep ens_size`
set ensemble_size  = `echo $ENSEMBLESTRING[3] | sed -e "s#,##"`
set num_ens        = ${ensemble_size}
echo "There are ${num_ens} ensemble members."

#-------------------------------------------------------------------------------
set inf_flavor_string = `grep -A 42 filter_nml input.nml | grep inf_flavor`
set inf_flavor     = `echo $inf_flavor_string[3] | sed -e "s#,##"`
echo "inflation flavor = "$inf_flavor

# start.nc is required for static_init_model subroutine in model_mod.f90 
$COPY $DART_ics_DIR/start_1.nc start.nc
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if ( $assim_job_start_from_init == .true. ) then

# follwoing processes are required for anly one time in working directory 
# when job start from initial
  #Copy dart exe
  $COPY $DART_DIR/models/LMDZ/work/dart_to_lmdz .
  $COPY $DART_DIR/models/LMDZ/work/lmdz_to_dart .
  $COPY $DART_DIR/models/LMDZ/work/filter .
  $COPY $DART_DIR/models/LMDZ/work/fill_inflation_restart .
  $COPY $DART_DIR/models/LMDZ/work/trans_time .
  $COPY $DART_DIR/models/LMDZ/work/advance_time .

  #Copy shell scripts 
  $COPY $DART_DIR/models/LMDZ/shell_scripts/advance_model.csh .
  $COPY $DART_DIR/models/LMDZ/shell_scripts/run_lmdz.csh .

  #--------------
  # check sampling_error_correction flag and if .true. then copy final_full.$num_ens file 
  # for sampling error corrections
   set sampling_error_correction_string = `grep -A 42 assim_tools_nml input.nml | grep sampling_error_correction`
   set sampling_error_correction  = `echo $sampling_error_correction_string[3] | sed -e "s#,##"`
   echo "sampling_error_correction is set to ${sampling_error_correction} "

  if ( $sampling_error_correction == .true.) then
     set precomputed_tables_dir = $DART_DIR/system_simulation/final_full_precomputed_tables

     if( -f $DART_DIR/system_simulation/final_full_precomputed_tables/final_full.$num_ens ) then
       $COPY  $precomputed_tables_dir/final_full.$num_ens .
     else
       echo "ERROR - Sampling error correction file final_full.$num_ens does not exists at $precomputed_tables_dir dir"
       exit
     endif

   endif  # sampling_error_correction = .true.
  #--------
  #copy start_#.nc,startphy_#.nc  & filter_ic_old.### from storage and rename -
  set n = 1
  while($n <= ${num_ens})
   $COPY $DART_ics_DIR/start_$n.nc .
   $COPY $DART_ics_DIR/startphy_$n.nc .
   set from = `printf "%s.%04d" $DART_ics_DIR/filter_ic_new $n`
   set to   = `printf "%s.%04d" filter_ic_old        $n`
   echo copying $from to $to
   $COPY $from $to
   @ n++
  end

  #----------------- alter the data timestamp in a DART restart file filter_ic_old --
  # fill required time in restart_file_tool_nml inside input.nml
  #if ($change_timestamp == 1) then
  #  ./restart_file_tool
  #endif
  #---------------
  if( $inf_flavor != 0 ) then
   # Determine the values of inf_sd_initial and inf_sd_initial from input.nml
   set inf_initial_string = `grep -A 42 filter_nml input.nml | grep inf_initial | tail -1`
   set inf_sd_initial_string = `grep -A 42 filter_nml input.nml | grep inf_sd_initial | tail -1`
   set inf_initial  = `echo $inf_initial_string[3] | sed -e "s#,##" `
   set inf_sd_initial  = `echo $inf_sd_initial_string[3] | sed -e "s#,##"`
   echo "inf_initial and inf_sd_initial  is $inf_initial and $inf_sd_initial "
 
   # create inflation initial restart files so that the values inf_initial_from_restart and 
   # inf_sd_initial_from_restart items in the &filter_nml namelist can be .TRUE. 
   # from the beginning. 
   echo $inf_initial $inf_sd_initial | ./fill_inflation_restart
   $COPY inflate_ics prior_inf_ic_old 
  endif # $inf_flavor

endif # assim_job_start_from_init = .true.

#---------------------------------------------------------------------------------
# set start date of assimilation 
set mon1=`printf %02d $start_month`
set day1=`printf %02d $start_day`

# set end date of assimilation 
set mon2=`printf %02d $end_month`
set day2=`printf %02d $end_day`

# convert start and end date in Gregorian day and second (since year 1601)
# so we can  compute total number of days including rolling over month and
# year boundaries. 
## the output of advance time with the -g input is:
#   gregorian_day_number  seconds
# use $var[1] to return just the day number

set start_d=(`echo ${start_year}${mon1}${day1}00 0 -g | ./advance_time`)
set end_d=(`echo ${end_year}${mon2}${day2}00 0 -g | ./advance_time`)

set curday=`echo ${start_year}${mon1}${day1}00 0 | ./advance_time`

 @ totaldays = ( $end_d[1] - $start_d[1] ) + 1

#----------------------------------------------------------------------------------
if ( $assim_job_start_from_init == .false.) then
# If job does not start from initial day in working  directory then copy restarrt files
# from corresponding OUTPUT_YYYYMMMDD directory
  set n = 1
  while($n <= ${num_ens})
   $COPY $OUTPUT_DIR/OUTPUT_${start_year}${mon1}${day1}/start_$n.nc .
   $COPY $OUTPUT_DIR/OUTPUT_${start_year}${mon1}${day1}/startphy_$n.nc .
   set from = `printf "%s.%04d" $OUTPUT_DIR/OUTPUT_${start_year}${mon1}${day1}/filter_ic_new $n`
   set to   = `printf "%s.%04d" filter_ic_old        $n`
   echo copying $from to $to
   $COPY $from $to 
   @ n++
  end

  $COPY $OUTPUT_DIR/OUTPUT_${start_year}${mon1}${day1}/prior_inf_ic_new prior_inf_ic_old

endif #assim_job_start_from_init = .false.
#*****************************************************************************
#--------------- Start assimiliation loop for ecah day obs file obs_seqYYYYMMDD
set day_count = 1

while ( $day_count <= $totaldays )
     
   set  year=`echo $curday | cut -b1-4`
   set month=`echo $curday | cut -b5-6`
   set   day=`echo $curday | cut -b7-8`
   echo starting assimilation for ${year} ${month} ${day}

   # Create dir for current day to restore output and restarts files
    mkdir $OUTPUT_DIR/OUTPUT_$year$month$day
 
   # Store all restart files at give day frequency : $restart_store_freq
   if ( $day_count % $restart_store_freq == 0 & $day_count > 1 ) then
     echo 'Storing Filter restart files in '$OUTPUT_DIR/OUTPUT_$year$month$day' directory'
     $COPY start_*.nc        $OUTPUT_DIR/OUTPUT_$year$month$day
     $COPY startphy_*.nc     $OUTPUT_DIR/OUTPUT_$year$month$day
     $COPY filter_ic_new.*   $OUTPUT_DIR/OUTPUT_$year$month$day
     $COPY prior_inf_ic_new  $OUTPUT_DIR/OUTPUT_$year$month$day
   endif

   # Link input obs_seqYYYYMMDD obs file for assimilation
   $REMOVE obs_seq.out
   if ( -f $OBS_DIR/obs_seq$year$month$day ) then
      $LINK $OBS_DIR/obs_seq$year$month$day obs_seq.out
   else
    echo ''  
    echo ERROR - file obs_seq$year$month$day  does not found in $OBS_DIR
    exit -obs_seq$year$month$day
   endif

   #-------run filter--------
   echo "=================================================================="
   echo "Assimilating Observation file obs_seq$year$month$day "
   echo "=================================================================="
 

    /usr/bin/time -p mpirun  -machinefile $PBS_NODEFILE -n $num_proc  ./filter > filter.log
 
   grep 'filter: End of main filter assimilation loop'  filter.log
   grep 'filter: End of main filter assimilation loop'  filter.log > /dev/null

   #-------copy outputs file in OUTPUT_#### directory ----------------------
   if ($status == 0) then
     echo "Filter Finished for observation obs_seq$year$month$day"
     $MOVE Posterior_Diag.nc $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE Prior_Diag.nc     $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE obs_seq.final     $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE histhf*.nc*       $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE histins*.nc*      $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE filter.log        $OUTPUT_DIR/OUTPUT_$year$month$day
     $MOVE dart_log.out      $OUTPUT_DIR/OUTPUT_$year$month$day
   else
     echo "WARNING - FILTER  stopped abnormally obs_seq$year$month$day"
     echo "========================================="
     exit -obs_seq$year$month$day
   endif

   #-------move prior_inf_ic_new to prior_inf_ic_old for next sequential run --------------------------
   if( $inf_flavor != 0 ) then 
    if (-e  prior_inf_ic_new ) then
     $COPY prior_inf_ic_new prior_inf_ic_old
    else
      echo ERROR - need prior_inf_ic_new to exist
      exit -obs_seq$year$month$day
    endif
   endif #$inf_flavor 

  #---- move filter_ic_new.## to filter_ic_old.## for next restart run -----------------------------
  set n = 1
  while($n <= ${num_ens})
      set from = `printf "%s.%04d" filter_ic_new $n`
      set to   = `printf "%s.%04d" filter_ic_old $n`
      echo copying $from to $to
      $COPY $from $to
      @ n++
  end
  #------ advance current day by 1d
  set curday=`echo $curday +1d | ./advance_time`
  @ day_count++

end  # day_count

$REMOVE  start.nc assim_model_state_* lmdz_out_temp* advance_temp* stok_paprs* hist*.nc 

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

