#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to assimilate observations using DART and OpenGGCM.
#
#=============================================================================
# This block of directives constitutes the preamble for the LSF queuing system
# LSF is used on the IMAGe Linux cluster 'coral'
# LSF is used on the IBM   'bluefire'
#
# the normal way to submit to the queue is:    bsub < run_filter
#
# an explanation of the most common directives follows:
# -J Job name (master script job.csh presumes filter_server.xxxx.log)
# -o STDOUT filename
# -e STDERR filename
# -P      account
# -q queue    cheapest == [standby, economy, (regular,debug), premium] == $$$$
# -n number of processors  (really)
##=============================================================================
#
#BSUB -J filter
#BSUB -o filter.%J.log
#BSUB -q regular
#BSUB -n 16
#BSUB -R "span[ptile=2]"
#BSUB -P 86850054
#BSUB -W 2:00
#BSUB -N -u ${USER}@ucar.edu
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
##
## the normal way to submit to the queue is:    qsub run_filter
##
## an explanation of the most common directives follows:
## -N     Job name
## -r n   Declare job non-rerunable
## -e <arg>  filename for standard error
## -o <arg>  filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=2   requests BOTH processors on the node. On both bangkok
##                     and calgary, there is no way to 'share' the processors
##                     on the node with another job, so you might as well use
##                     them both. (ppn == Processors Per Node)
##=============================================================================
#
#PBS -N filter
#PBS -e filter.err
#PBS -o filter.log
#PBS -l nodes=3:ppn=32
#PBS -l walltime=2:00:00

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv MPI         mpirun.lsf

else if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS   todo ... where does xx come from
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv MPI         aprun -np xx

   env | sort

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     POP
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv MPI         csh

endif

#----------------------------------------------------------------------
# Just an echo of the job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"
echo "${JOBNAME} ($JOBID) started      at "`date`
echo

#----------------------------------------------------------------------
#----------------------------------------------------------------------

set RUNPATH=/mnt/lustre/lus0/home/dart/dart/rma_openggcm/models/openggcm/timtest
cd $RUNPATH

set myname = $0          # this is the name of this script

# some systems don't like the -v option to any of the following 

set OSTYPE = `uname -s`
switch ( ${OSTYPE} )
   case IRIX64:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   case AIX:
      setenv REMOVE 'rm -rf'
      setenv   COPY 'cp -p'
      setenv   MOVE 'mv -f'
      breaksw
   default:
      setenv REMOVE 'rm -rvf'
      setenv   COPY 'cp -vp'
      setenv   MOVE 'mv -fv'
      breaksw
endsw

echo "${JOBNAME} ($JOBID) RUNPATH == $RUNPATH"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = /mnt/lustre/lus0/home/dart/dart/rma_openggcm/models/openggcm
set OBSERVATIONDIR = $DARTDIR/observations

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
#-----------------------------------------------------------------------------

# executables
${COPY} ${DARTDIR}/work/filter                     .
${COPY} ${DARTDIR}/work/input.nml                  .

set SAMP_ERR_FILE = ${DARTROOT}/assimilation_code/programs/system_simulation/work/sampling_error_correction_table.Lanai.nc
${COPY} ${SAMP_ERR_FILE}  sampling_error_correction_table.nc

# COMPLETE HACK ...

${COPY} $DARTDIR/data/Data.ionos2.nc input_prior_inflation_mean.nc
${COPY} $DARTDIR/data/Data.ionos2.nc input_prior_inflation_sd.nc

@ ncycle = 0
while ( 1 ) 
   
   find . -name to_dart.semaphore > semaphore_files.txt
   set nfiles = `cat semaphore_files.txt | wc -l`
   if ($nfiles == $ens_size) then

      #------------------------------------------------------------
      # Get the time 

      @ncycle ++  # debug
      if ( $ncycle > 3 ) break

      set semaphore = `head -n 1 semaphore_files.txt`
      set TimeString = `head -n 1 $semaphore`

      echo "DART: Running assimilation for $TimeString at "`date`

      if ( $TimeString == "done" ) exit

      \rm `find . -name to_dart.semaphore`

      #------------------------------------------------------------
      # link the obs

      ${LINK} ${OBSERVATIONDIR}/obs_seq.$TimeString obs_seq.out

      #------------------------------------------------------------
      # input.nml presumes 'input_file_list.txt', 'output_file_list.txt'

      ls -1 $RUNPATH/hoptest_*/target/DATA.ionos2.nc >! input_file_list.txt
      sed -e "s/DATA.ionos2.nc/dart_posterior.nc/" input_file_list.txt >! output_file_list.txt

      ${MPI} ./filter

      #------------------------------------------------------------
      # rename output files (inflation, diagnostics, etc)
     
      ls -lrt

      #------------------------------------------------------------
      # update inflation files 

      #------------------------------------------------------------
      # signal that DART is done
      # remove the 'to_dart' semaphore
      # create the 'from_dart' semaphore

      foreach FILE ( `cat semaphore_files.txt` )
         \rm $FILE
         set FILE = `echo $FILE | sed -e "s/to_dart/from_dart/"`
         echo "600"                          >! $FILE
         echo "Assimilation for $TimeString" >> $FILE
      end
   endif

   sleep 1
   
end

echo "${JOBNAME} ($JOBID) finished at "`date`
echo "Listing contents of RUNPATH before archiving"
ls -l

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
