#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Script to assimilate observations using DART and the ucoam ocean model.
# This presumes two directories exists that contain all the required bits
# for ucoam and for DART.
#
# It is entirely fine to have BOTH directives in the same file.
# I guarantee that as soon as you delete one set, you will change machines
# and wish you had not deleted the directives. 
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
#PBS -N pmo
#PBS -r n
#PBS -e pmo.err
#PBS -o pmo.log
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
#BSUB -J pmo
#BSUB -o pmo.%J.log
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
   setenv MPI         mpirun

else if ($?LSB_QUEUE) then

   #-------------------------------------------------------------------
   # This is used by LSF
   #-------------------------------------------------------------------

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST
   setenv MPI         mpirun.lsf

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     pmo
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
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#----------------------------------------------------------------------

set nonomatch  # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following

switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      setenv   MOVE '/usr/local/bin/mv -fv'
      setenv   COPY '/usr/local/bin/cp -fv --preserve=timestamps'
      setenv   LINK '/usr/local/bin/ln -fvs'
      setenv REMOVE '/usr/local/bin/rm -fr'

      setenv    TEMPDIR /glade/proj3/image/thoar/
      setenv BASEOBSDIR /glade/proj3/image/Observations/WOD09
   breaksw

   case ys*:
      # NCAR "yellowstone"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv    TEMPDIR /glade/scratch/thoar/
      setenv BASEOBSDIR /glade/p/image/Observations/WOD09
   breaksw

   default:
      # SDSU "dulcinea"
      setenv   MOVE 'mv -fv'
      setenv   COPY 'cp -fv --preserve=timestamps'
      setenv   LINK 'ln -fvs'
      setenv REMOVE 'rm -fr'

      setenv    TEMPDIR /gcemproject/thoar/
      setenv BASEOBSDIR /home/thoar/svn/DART/UCOAM/models/GCOM/work
      setenv    DARTDIR /home/thoar/svn/DART/UCOAM/models/GCOM
      setenv   SERUCOAM /home/mgarcia/UCOAM-Moh/angie/serucoam

   breaksw
endsw

#-------------------------------------------------------------------------
# Create temporary working directory for the perfect model and go there
#-------------------------------------------------------------------------

echo "`date` -- BEGIN GENERATE GCOM TRUE STATE"

set temp_dir = pmo_gcom
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name
# of the form "./${MYCASE}.ucoam.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

#set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
#set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
#set OCN_DAY      = `echo $OCN_DATE[3] | bc`
#set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
#set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`
#
#echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
#echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of ucoam.
#-----------------------------------------------------------------------------

# set YYYYMM   = `printf %04d%02d                ${OCN_YEAR} ${OCN_MONTH}`
# set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
# set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

set OBS_FILE = ${BASEOBSDIR}/obs_seq.in

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#
# DART namelist settings required:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name    = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name   = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days          = -1,
# &perfect_model_obs_nml:  init_time_seconds       = -1,
# &perfect_model_obs_nml:  first_obs_days          = -1,
# &perfect_model_obs_nml:  first_obs_seconds       = -1,
# &perfect_model_obs_nml:  last_obs_days           = -1,
# &perfect_model_obs_nml:  last_obs_seconds        = -1,
# &ucoam_to_dart_nml:        ucoam_to_dart_output_file = 'dart_ics'
#=========================================================================

if ( ! -e ${DARTDIR}/work/input.nml ) then
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${DARTDIR}/work/input.nml not found ... ERROR"
   exit -2
endif

sed -e "s#dart_ics#perfect_ics#" < ${DARTDIR}/work/input.nml >! input.nml

${COPY} ${DARTDIR}/work/perfect_model_obs          . || exit -1
${COPY} ${DARTDIR}/work/dart_to_gcom               . || exit -1
${COPY} ${DARTDIR}/work/gcom_to_dart               . || exit -1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh . || exit -1

set OCN_RESTART_FILENAME = gcom_restart.nc
set OCN_GEOMETRY_FILENAME = gcom_geometry.nc

# This would be better in the long run, but for now ...
# ${COPY} ${SERUCOAM}/OUTPUT/$OCN_RESTART_FILENAME   . || exit -1
${COPY} ${SERUCOAM}/OUTPUT/$OCN_GEOMETRY_FILENAME  . || exit -1

${COPY} ${DARTDIR}/data/gcom_restart_1986-01-01.43200.nc $OCN_RESTART_FILENAME || exit -1

#=========================================================================
# Block 2: Convert 1 ucoam restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
# that came from the contents of the pointer file ../rpointer.ocn.restart
#=========================================================================

echo "`date` -- BEGIN GCOM-TO-DART"

./gcom_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'ucoam_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'ucoam_to_dart' ... ERROR"
   exit -3
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

./perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

${MOVE} True_State.nc    ../ucoam_True_State.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.perfect  ../ucoam_obs_seq.${OCN_DATE_EXT}.perfect
${MOVE} dart_log.out     ../ucoam_dart_log.${OCN_DATE_EXT}.out

echo "`date` -- END   ucoam PERFECT_MODEL_OBS"

#=========================================================================
# Block 4: Update the ucoam restart file
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

# Eat the cookie regardless
${REMOVE} ../ucoam_inflation_cookie
${REMOVE} perfect_ics dart_log.nml

echo "`date` -- END   GENERATE ucoam TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

