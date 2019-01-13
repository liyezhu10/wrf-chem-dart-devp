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
#PBS  -o compress_2010-07-17-21600.out
#PBS  -j oe 

# -------------------------
# Purpose

# (De)Compress files from the forecast or assimilation. 

# -------------------------
# Method

# This script can compress as well as decompress sets of files.
# When called from a script (normally assimilate.csh), it compresses the files.
# When called from the command line, it can uncompress sets of files,
# if this script has been edited to use the right metadata to
# construct the expected data directories.

# It operates on all files in parallel using the "command file" mechanism.
# The file contains a separate task on each line.
# It has syntax appropriate for the launch script (below), 
# which may not be the same as this script.  This must include
# syntax to put stderr and stdout in a single file, e.g. ( &> ) for bash.
# It is dispatched by the local MPI script to perform N simultaneous operations. 

# Compression method can depend on the file type, and may include lossy compression.
# It is most often called by assimilate.csh, but can be run as a batch job.
# Assimilate.csh runs this in 2 places:
# 1) Every cycle: 
#    +  all the cpl history (forcing) files.
#    +  DART output
#       >  stages of state files  
#             mean, sd  (no inst)
#       >  obs_seq.final (no inst)   70% of 1 Gb in 35 sec
#       >  not inflations (no inst)  
# 1) Before archiving a restart set to archive/rest; all large restart files.

# -------------------------

if ($#argv == 5) then
   # Called from assimilate.csh (or other script).
   set comp_cmd      = 'gzip'
   set case_name     = $1
   set ymds          = $2
   set ensemble_size = $3
   set sets          = ($4)
   set stages        = ($5)
   set data_dir      = '.'

else if ($#argv == 0) then
   # Edit these and run as a batch job.
   set comp_cmd      = 'gunzip'
   set case_name     = CESM2_1_80_3node
   set ymds          = 2010-07-17-64800
   set ensemble_size = 80
   set sets          = (hist dart)
   # set sets          = (clm2 cpl cam cice)
   set stages        = (preassim output)
   set data_dir      = /gpfs/fs1/scratch/raeder/${case_name}/run
   # set data_dir      = /gpfs/fs1/scratch/raeder/${case_name}/archive/rest/2010-07-15-00000
else 
   echo "Usage: Cannot run interactively."
   echo "       > vi compress.csh {set it up}; qsub compress.csh"
   echo '   OR, in a calling script'
   echo '       ${scr_dir}/compress.csh case_name YYYY-MM-DD-SSSS ensemble_size "sets" "stages"'
   echo '   sets   = 1 or more of {clm2 cpl cam cice hist dart} to compress, separated by spaces'
   echo '   stages = 1 or more of stages {input, preassim, postassim, output} to compress.'
   exit 17
endif

echo "In compress.csh:"
echo "   case_name     = $case_name"
echo "   date          = $ymds"
echo "   ensemble_size = $ensemble_size"
echo "   sets          = $sets"
echo "   stages        = $stages"
echo "   data dir      = $data_dir"

cd $data_dir

# ==========================
# Fail if there are leftover diagnostics files.
ls *.eo > /dev/null
if ($status == 0) then
   echo "ERROR; Existing compression diagnostic files: *.eo.  Exiting"
   exit 37
endif

# --------------------------
# Environment and commands.

# PBS requires 'setenv MPI_SHEPHERD true' for the command file to work correctly.
setenv MPI_SHEPHERD true
module list

setenv MPICMD    "/gpfs/u/home/raeder/Scripts/mpiexec_mpt_debug"
setenv LAUNCHCMD `pwd`/launch_cf.sh
echo "launch command is $LAUNCHCMD"

cat << EndOfFile >! launch_cf.sh
#!/bin/bash
line=\$(expr \$PMI_RANK + 1)
INSTANCE=\$(sed -n \${line}p \$1)
eval "\$INSTANCE"
EndOfFile
# Doesn't work; ENDOFFILE doess not stop the cat.
# cat << 'ENDOFFILE' >! launch_cf.sh
# #!/bin/bash
# line=$(expr $PMI_RANK + 1)
# INSTANCE=$(sed -n ${line}p $1)
# eval "$INSTANCE"
# ENDOFFILE

chmod u+x launch_cf.sh

set cmd = `echo $comp_cmd | cut -d' ' -f1`
if ($cmd == 'gzip') then
   set ext = ''
else if ($cmd == 'gunzip') then
   set ext = '.gz'
else
   echo "unrecognized command $cmd.  Don't know which extension to use"
   exit 27
endif

# ------------------------------------------------------------------------------
# Create the command file where each line is a separate command, task, operation, ....

rm -f mycmdfile
touch mycmdfile

# 'task' is incremented continuously over all files; components, members, etc.
# 'task' is a running counter of jobs in mycmdfile.
set task = 0

foreach comp ( $sets )
switch ($comp)
   case {clm2,cpl,cam,cice}:
      # Restart files, ordered by decreasing size.
      echo "comp = $comp"
      set i=1
      while ( $i <= $ensemble_size)
         # E.g. CAM6_80mem.cice_0001.r.2010-07-15-00000.nc
         set file_name = `printf "%s.%s_%04d.r.%s.nc%s" $case_name $comp $i $ymds $ext`
         echo "   $file_name"
   
         # If the expected file exists, add the compression command 
         if (-f $file_name) then
            @ task++
            echo "$comp_cmd $file_name &> compress_${task}.eo " >> mycmdfile
         else
            echo 'Could not find "'$file_name'" to compress.'
         endif

         @ i++
      end
      breaksw

   case hist:
      # Coupler history (forcing) files, ordered by decreasing size .
      echo "comp = $comp"
      foreach type ( ha ha2x1d hr2x ha2x3h ha2x1h ha2x1hi )
         set i=1
         while ( $i <= $ensemble_size)
            # E.g. CAM6_80mem.cpl_0001.ha.2010-07-15-00000.nc
            set file_name = `printf "%s.cpl_%04d.%s.%s.nc%s" $case_name $i $type $ymds $ext`
      
            if (-f $file_name) then
               @ task++
               if ($type != 'ha') then
                  echo "($comp_cmd $file_name; sleep 1) &> compress_${task}.eo" >> mycmdfile
                  # The 1 second sleep is needed to prevent compression of small files
                  # from finishing before all the commands have been started.
                  # This can cause mpiexec_mpt to exit with code 137 and no useful messages.
                  # 'date --rfc-3339=ns' can be added before and after the '$comp_cmd $file_name'
                  # to help debug timing problems.
               else
                  echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
               endif
            else
               echo 'Could not find "'$file_name'" to compress.'
            endif

            @ i++
         end
      end
      breaksw

   case dart:
      # Diagnostic files.
      # It is not worthwhile to compress inflation files ... small, not many files
      # It is also not clear that binary observation sequence files compress effectively.
      echo "comp = $comp"

      # obs_seq.final (no inst)   70% of 1 Gb (ascii) in 35 sec
      # E.g. CAM6_80mem.dart.e.cam_obs_seq_final.2010-07-15-00000
      set file_name = ${case_name}.dart.e.cam_obs_seq_final.${ymds}${ext}
      if (-f $file_name) then
         @ task++
         echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
      endif

      foreach stage ($stages)
         foreach stat ( 'mean' 'sd' )
            # E.g. CAM6_80mem.e.cam_output_mean.2010-07-15-00000.nc
            # E.g. CAM6_80mem.e.cam_output_sd.2010-07-15-00000.nc
            set file_name = ${case_name}.dart.e.cam_${stage}_${stat}.${ymds}.nc${ext}
            if (-f $file_name) then
               @ task++
               echo "$comp_cmd $file_name &> compress_${task}.eo" >> mycmdfile
            endif
         end

         set i=1
         while ( $i <= $ensemble_size)
            # E.g. CAM6_80mem.cam_0001.e.preassim.2010-07-15-00000.nc
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

# ------------------------------------------------------------------------------
# Launch the jobs.

date --rfc-3339=ns

$MPICMD -n $task $LAUNCHCMD ./mycmdfile 

set mpt_status = $status
echo "mpt_status = $mpt_status at "
date --rfc-3339=ns


# Check the statuses?
if ( -f compress_1.eo ) then
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
   rm *.eo  mycmdfile
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

