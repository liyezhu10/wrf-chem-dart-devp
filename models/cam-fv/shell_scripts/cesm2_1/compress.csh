#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#PBS  -N compress.csh
#PBS  -A P86850054
#PBS  -q premium
#PBS  -l select=16:ncpus=36:mpiprocs=36
#PBS  -l walltime=00:20:00
#PBS  -o CESM2_1_80_3node.out
#PBS  -j oe 

# Most often called by assimilate.csh to compress a batch of files. 
# Compression will depend on the file type,
# and may include lossy compression.
# > Execute something like this 2 places:
#      Before archiving a restart set to archive/rest; all large restart files.
#         clm2.r. compress 93.5% of 7.7 Gb in 1:23.
#         .i. files only gzip compress ~12% in 24 seconds
#      Every cycle: 
#         all the cpl history (forcing) files.
#       ? DART output
#            stages of state files  (worth it?  only 15% in 10s)
#               mean, sd  (no inst)
#            obs_seq.final (no inst)   70% of 1 Gb in 35 sec
#            inflations (no inst)  (should be like stages, but not many files)
#      These 2 have different dates! $save_date  vs. $ATM_DATE_EXT
# 
# ? In which directory will assimilate.csh be when this executes?
# ? Where will it find compres.csh?

if ($#argv == 5) then
   # Called from assimilate.csh (or other script).
   set comp_cmd = 'gzip -k'
   set case_name     = $1
   set ymds          = $2
   set ensemble_size = $3
   set sets          = ($4)
   set stages        = ($5)
   set data_dir      = '.'

else if ($#argv == 0) then
   # Edit these and run as a batch job.
   set comp_cmd      = 'gunzip -k'
   set case_name     = CESM2_1_80_3node
   set ymds          = 2010-07-17-64800
   set ensemble_size = 80
   set sets          = (clm2 cpl cam cice)
   set stages        = (preassim output)
   set data_dir      = /glade/scratch/${USER}/${case_name}/run/Compressed2

else 
   echo "Usage: Cannot run interactively."
   echo "       > vi compress.csh {set it up}; qsub compress.csh"
   echo '   OR, in a calling script'
   echo '       ${scr_dir}/compress.csh case_name YYYY-MM-DD-SSSS ensemble_size "sets" "stages"'
   echo '   sets   = 1 or more of {clm2 cpl cam cice hist dart} to compress, separated by spaces'
   echo '   stages = 1 or more of stages {input, preassim, postassim, output} to compress.'
   exit 17
endif

set cmd = `echo $comp_cmd | cut -d' ' -f1`
if ($cmd == 'gzip') then
   set ext = ''
else if ($cmd == 'gunzip') then
   set ext = '.gz'
else
   echo "unrecognized command $cmd.  Don't know which extension to use"
   exit 27
endif

echo "In compress.csh:"
echo "   case_name     = $case_name"
echo "   date          = $ymds"
echo "   ensemble_size = $ensemble_size"
echo "   sets          = $sets"
echo "   stages        = $stages"
echo "   data dir      = $data_dir"

echo "TJH SKIPPING COMPRESSION ..."
echo "TJH SKIPPING COMPRESSION ..."
echo "TJH SKIPPING COMPRESSION ..."
echo "TJH SKIPPING COMPRESSION ..."
exit 0

cd $data_dir

# ==========================
# Fail if ther are leftover diagnostics files.
ls *.eo > /dev/null
if ($status == 0) then
   echo "ERROR; Existing compression diagnostic files: *.eo.  Exiting"
   exit 37
endif

# ==========================
# Command in assimilate.csh to create the command file.
# compress.csh $case $date $ensemble_size "$sets" "$stages"
# ==========================

# Not needed?
setenv MPI_SHEPHERD true

# The default launch_cf.sh has a 1 second sleep between each command.
# That was to prevent some commands from finishing 
# before they had all been started, which crashed the job.
# My jobs will all be much longer than it takes to submit the jobs,
# so use a launch_cf.sh without the sleep:
# setenv RUNCMD "mpiexec_mpt /gpfs/u/home/raeder/Scripts/launch_cf.sh"

setenv RUNCMD "/gpfs/u/home/raeder/Scripts/mpiexec_mpt_debug /gpfs/u/home/raeder/Scripts/launch_cf.sh"
echo "WARNING; $RUNCMD may have a 1 second delay in it.   Remove if all commands take > a few seconds"

# Set up the command file for receiving commands.
rm -f mycmdfile
touch mycmdfile
# This is incremented continuously over all files; components, members, etc.
set task = 0

# - - - - - - - - - 
# Loop over file types to be compressed

# Restart files, ordered by decreasing size.
# cice?  Since I have the processors and it might not take any more wall clock.
foreach comp ( $sets )
switch ($comp)
   # These 3 cases will all execute the indented code below 'cice',
   # to compress restart files.
   case clm2:
   case cpl:
   case cam:
   case cice:
      # Compress the slowest (largest, most extensive compression) first,
      # so that they can start while the faster jobs are being submitted.
      # This is s minor consideration since I removed the sleep command from launch_cf.sh.
      # Loop over instance number.
      echo "comp = $comp"
      set i=1
      while ( $i <= $ensemble_size)
         set file_name = `printf "%s.%s_%04d.r.%s.nc%s" $case_name $comp $i $ymds $ext`
         echo "   $file_name"
   
         # create a command file where each line calls a script with
         # unique arguments, including a unique work directory name.
         if (-f $file_name) then
            # increment task number, which is a running counter of jobs in mycmdfile.
            @ task++
            echo "$comp_cmd $file_name &> compress_${task}.eo " >> mycmdfile
#             echo "(date --rfc-3339=ns; $comp_cmd $file_name ; date --rfc-3339=ns ) &> compress_${task}.eo " >> mycmdfile
         else
            echo "Could not find $file_name"
         endif

         # advance the loop counter
         @ i++
      end
      breaksw

   case hist:
      # Coupler history (forcing) files, ordered by decreasing size 
      # (Check this when they're working right!).
      echo "case_name = $case_name"
      echo "comp = $comp"
#       foreach type ( 'ha' 'hr2x' 'ha2x3h' 'ha2x1h' 'ha2x1hi')
      foreach type ( ha ha2x1d hr2x ha2x3h ha2x1h ha2x1hi)
         # Loop over instance number
         set i=1
         while ( $i <= $ensemble_size)
            set file_name = `printf "%s.cpl_%04d.%s.%s.nc%s" $case_name $i $type $ymds $ext`
      
            if (-f $file_name) then
               @ task++
               if ($type != 'ha') then
#                   echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
                  # The 3 second sleep is needed to prevent compression of small files
                  # from finishing before all the commands have been started.
                  echo "(date --rfc-3339=ns; $comp_cmd $file_name; date --rfc-3339=ns; sleep 3 ) &> compress_${task}.eo"\
                       >> mycmdfile
               else
                  echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
               endif
            else
               echo "$file_name does not exist"
            endif

            @ i++
         end
      end
      breaksw

   case dart:
      # obs_seq.final (no inst)   70% of 1 Gb in 35 sec
      echo "comp = $comp"
      set file_name = ${case_name}.dart.e.cam_obs_seq_final.${ymds}${ext}
      if (-f $file_name) then
         @ task++
         echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
#          echo "$comp_cmd $file_name &> compress_${task}.eo; date --rfc-3339=ns " >> mycmdfile
      endif

      #    ${case_name}.dart.e.cam_${stages}_mean.2010-07-16-00000.nc${ext}
      #    ${case_name}.dart.e.cam_${stages}_sd.2010-07-16-00000.nc${ext}
      foreach stage ($stages)
         foreach stat ( 'mean' 'sd' )
            # Ensemble means and standard deviations for all stages.
            set file_name = ${case_name}.dart.e.cam_${stage}_${stat}.${ymds}.nc${ext}
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            endif
   
            # Inflation files
            #    Don't compress, save non-output inflations?
            # ${case_name}.dart.rh.cam_${stages}_priorinf_mean.2010-07-17-00000.nc${ext}
            # ${case_name}.dart.rh.cam_${stages}_priorinf_sd.2010-07-17-00000.nc${ext}
#             set file_name = ${case_name}.dart.rh.cam_${stage}_priorinf_${stat}.${ymds}.nc${ext}
#             if (-f $file_name) then
#                @ task++
#                echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
#             else
#                echo "Could not find $file_name"
#             endif
         end
         # stages of state files  (worth it?  only 15% in 10s, but I'm "using" the cores anyway)
         # mean, sd  (no inst)
         
         # Loop over instance number
         # CESM2_1_80_3node.cam_0001.e.preassim.2010-07-16-00000.nc${ext}
         set i=1
         while ( $i <= $ensemble_size)
            set file_name = `printf "%s.cam_%04d.e.%s.%s.nc%s" $case_name $i $stage $ymds $ext`
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            else
               echo "Could not find $file_name"
            endif
            @ i++
         end
      end
      breaksw
      
   default:
      breaksw
endsw
end

# - - - - - - - - - 

echo "Before launching mycmdfile"
echo "The MPI libraries used by mpiexec_mpt_debug are:"
/gpfs/u/home/raeder/Scripts/mpiexec_mpt_debug -v
echo ' '

# actually launch the jobs here
date --rfc-3339=ns
$RUNCMD ./mycmdfile 
set mct_status = $status
echo "mct_status = $mct_status"
# Allow output from mycmdfile to be written even when RUNCMD fails
# before error checking.
sleep 4

# How to know when they are done?
# launch_cf.sh returns after all commands are finished.
echo "after launching mycmdfile; did commands finish before mpiexec_mpt returned?"
date --rfc-3339=ns

if ($mct_status == 0) then
   echo "$RUNCMD supposedly worked"
else
   echo "$RUNCMD ./mycmdfile supposedly failed.  Commands will be checked below"
#  RUNCMD may be generating spurious failure messages and codes.
#    exit
endif

# Check the stati?
if ( -f compress_1.eo) then
   grep $cmd *.eo
   # grep failure = compression success = "not 0"
   set gr_stat = $status
   echo "gr_stat when eo file exists = $gr_stat"
else
   # No eo files = failure of something besides g(un)zip.
   set gr_stat = 0
   echo "gr_stat when eo file does not exist = $gr_stat"
endif

if ($gr_stat == 0) then
   echo "compression failed.  See .eo files with non-0 sizes"
   echo "Remove .eo files after failure is resolved."
   exit 190
else
   # Compression worked; clean up the .eo files and mycmdfile
#    rm *.eo  mycmdfile
   mv *.eo  mycmdfile ..
   # How to remove the files?  This, or compile a list as I compress, above.
   # Or remove the -k argument from gzip.
   if ($cmd == 'gzip') then
      foreach f (`ls ${case_name}*.gz`) 
         if (`ls -l $f | cut -d' ' -f 5` > 0) then
            rm $f:r
         else
            echo "$f has size 0.  Not removing original"
         end
      end
   else if ($cmd == 'gunzip') then
      foreach f (`ls ${case_name}*.nc`) 
         if (`ls -l $f | cut -d' ' -f 5` > 0) then
            rm ${f}.gz
         else
            echo "$f has size 0.  Not removing original"
         end
      end
      rm ${case_name}*obs_seq*${ymds}.gz
   endif
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

