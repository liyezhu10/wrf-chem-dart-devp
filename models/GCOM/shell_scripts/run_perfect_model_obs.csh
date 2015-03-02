#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Top level script to generate observations and a TRUE state.
#
# This script only processes a single observation file.
# Still fairly complex; requires a raft of
# data files and most of them are in hardcoded locations.
#
# This script is designed to be executed as a batch job. 
# However, if you comment out the model advance part - you can run this
# interactively to check for logic, file motion, syntax errors and the like.
# It is entirely fine to have directives for both PBS and LSF in the same file.
# I guarantee that as soon as you delete one set, you will change machines
# and wish you had not deleted the directives. 
#
# The script moves the necessary files to a temporary directory that is the
# basis for the DART experiment; this will be called CENTRALDIR. 
#
##=============================================================================
## This block of directives constitutes the preamble for the PBS queuing system
##
## the normal way to submit to the queue is:    qsub run_filter
##
## an explanation of the most common directives follows:
## -N <arg>   Job name
## -r n       Declare job non-rerunable
## -e <arg>   filename for standard error
## -o <arg>   filename for standard out
## -q <arg>   Queue name (small, medium, long, verylong)
## -l nodes=xx:ppn=16   request xx nodes and 16 processors on each node. 
##=============================================================================
#
#PBS -N gcom_pmo
#PBS -r n
#PBS -e gcom_pmo.err
#PBS -o gcom_pmo.log
#PBS -q medium
#PBS -l nodes=2:ppn=16
#
#=============================================================================
## This block of directives constitutes the preamble for the LSF queuing system
##
## the normal way to submit to the queue is:    bsub < run_filter
##
## an explanation of the most common directives follows:
## -J <arg>      Job name (master script job.csh presumes filter_server.xxxx.log)
## -o <arg>      output listing filename
## -q <arg>      queue
## -n <arg>      number of processors  (really)
## -P <arg>      account
## -W <arg>      wall-clock hours:minutes required
## -N -u <arg>   mail this user when job finishes
##=============================================================================
#
#BSUB -J gcom_pmo
#BSUB -o gcom_pmo.%J.log
#BSUB -q regular
#BSUB -n 1
#BSUB -P 8685xxxx
#BSUB -W 2:00
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?PBS_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by PBS
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $PBS_O_WORKDIR
   setenv JOBNAME     $PBS_JOBNAME
   setenv JOBID       $PBS_JOBID
   setenv MYQUEUE     $PBS_QUEUE
   setenv MYHOST      $PBS_O_HOST
   setenv RUN_CMD     "mpirun -np 1 -machinefile $PBS_NODEFILE"

else if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv RUN_CMD     mpirun.lsf

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     pmo
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      $HOST
   setenv RUN_CMD     csh

endif

#----------------------------------------------------------------------
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#----------------------------------------------------------------------

set nonomatch  # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following

switch ("`hostname`")

   case ys*:
      # NCAR "yellowstone"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv EXPERIMENT /glade/p/work/${USER}/${JOBNAME}
      setenv CENTRALDIR /glade/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv BASEOBSDIR ${HOME}/work/DART/UCOAM/models/GCOM/work
      setenv    DARTDIR ${HOME}/work/DART/UCOAM/models/GCOM
      setenv   SERUCOAM ${HOME}/UCOAM-Moh/angie/serucoam
   breaksw

   default:
      # SDSU "dulcinea"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv EXPERIMENT /gcemproject/${USER}/${JOBNAME}
      setenv CENTRALDIR /raid/scratch/${USER}/${JOBNAME}/job_${JOBID}
      setenv BASEOBSDIR /home/${USER}/svn/DART/UCOAM/models/GCOM/work
      setenv    DARTDIR /home/${USER}/svn/DART/UCOAM/models/GCOM
      setenv   SERUCOAM /home/mgarcia/UCOAM-Moh/angie/serucoam

   breaksw
endsw

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory.
#----------------------------------------------------------------------

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

#----------------------------------------------------------------------
# Just an echo of the job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $HOST"
echo "${JOBNAME} ($JOBID) started       at "`date`
echo "${JOBNAME} ($JOBID) CENTRALDIR    is $CENTRALDIR"
echo

#=========================================================================
# Block 1: Populate CENTRALDIR with everything needed to run DART and GCOM.
#
# Get the DART executables, scripts, and input files
# The input.nml will be copied from the DART directory and modified appropriately.
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/perfect_model_obs          .  || exit 1
${COPY} ${DARTDIR}/work/dart_to_gcom               .  || exit 1
${COPY} ${DARTDIR}/work/gcom_to_dart               .  || exit 1
${COPY} ${BASEOBSDIR}/obs_seq.in                   .  || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .  || exit 1

#-----------------------------------------------------------------------------
# Get the ucoam executable, control files, and data files.
#-----------------------------------------------------------------------------

mkdir -p   GRID || exit 1
mkdir -p  PARAM || exit 1
mkdir -p OUTPUT || exit 1

${COPY} ${SERUCOAM}/Main                     gcom.serial.exe || exit 1
${COPY} ${SERUCOAM}/GRID/Grid.dat            GRID            || exit 1
${COPY} ${SERUCOAM}/GRID/ProbSize.dat        GRID            || exit 1
${COPY} ${SERUCOAM}/PARAM/param.dat          PARAM           || exit 1
${COPY} ${SERUCOAM}/OUTPUT/gcamIC20secII.nc  OUTPUT/gcom_restart.nc || exit 1

#=========================================================================
# Block 2: Convert 1 ucoam restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
#=========================================================================
#
# DART namelist settings required:
#
# &perfect_model_obs_nml:  start_from_restart       = .true.
# &perfect_model_obs_nml:  async                    = 2
# &perfect_model_obs_nml:  restart_in_file_name     = 'perfect_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name     = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name    = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days           = -1,
# &perfect_model_obs_nml:  init_time_seconds        = -1,
# &perfect_model_obs_nml:  first_obs_days           = -1,
# &perfect_model_obs_nml:  first_obs_seconds        = -1,
# &perfect_model_obs_nml:  last_obs_days            = -1,
# &perfect_model_obs_nml:  last_obs_seconds         = -1,
# &model_nml:              gcom_restart_file        = 'gcom_restart.nc'
# &model_nml:              gcom_geometry_file       = 'gcom_geometry.nc'
# &gcom_to_dart_nml:       gcom_to_dart_output_file = 'perfect_ics'
#=========================================================================

if ( ! -e ${DARTDIR}/work/input.nml ) then
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   exit -2
endif

# Ensure the namelist has the values required by this script.

sed -e "s/ start_from_restart /c\ start_from_restart = .true." \
       "s/ async /c\ async = 2" \
       "s/ restart_in_file_name /c\ restart_in_file_name = 'perfect_ics'" \
       "s/ obs_sequence_in_name /c\ obs_sequence_in_name = 'obs_seq.in'" \
       "s/ obs_sequence_out_name /c\ obs_sequence_out_name = 'obs_seq.perfect'" \
       "s/ gcom_restart_file /c\ gcom_restart_file = 'gcom_restart.nc'" \
       "s/ gcom_geometry_file /c\ gcom_geometry_file = 'gcom_geometry.nc'" \
       "s/ gcom_to_dart_output_file /c\ gcom_to_dart_output_file = 'perfect_ics'" \
       ${DARTDIR}/work/input.nml >! input.nml

echo "`date` -- BEGIN GCOM-TO-DART"

${LINK} gcom_restart.nc gcom_geometry.nc

${RUN_CMD} ./gcom_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'gcom_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'gcom_to_dart' ... ERROR"
   exit -2
endif

echo "`date` -- END GCOM-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#=========================================================================

echo "`date` -- BEGIN ucoam PERFECT_MODEL_OBS"

${RUN_CMD} ./perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

echo "`date` -- END   ucoam PERFECT_MODEL_OBS"

#=========================================================================
# Block 4: Copy/Move the good stuff back to someplace safe.
# CENTRALDIR is usually on a volatile or temporary filesystem.
# EXPERIMENT is usually someplace long-term.

${MOVE} True_State.nc    ${EXPERIMENT}
${MOVE} obs_seq.perfect  ${EXPERIMENT}
${MOVE} dart_log.out     ${EXPERIMENT}
${MOVE} input.nml        ${EXPERIMENT}
${MOVE} *.csh            ${EXPERIMENT}
${MOVE} $myname          ${EXPERIMENT}

echo "Listing contents of CENTRALDIR after archiving (miss anything?)"
ls -l

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

