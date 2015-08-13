#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script is designed to be submitted as a batch job but may be run from 
# the command line (as a single thread) to check for file motion, etc.
# If running interactively, please comment out the part that actually runs filter.
#
#-----------------------------------------------------------------------------
#
#BSUB -J jules_perfect
#BSUB -o jules_perfect.%J.log
#BSUB -P P3507xxxx
#BSUB -q premium
#BSUB -n 1
#BSUB -W 1:00
#BSUB -N -u ${USER}@ucar.edu

#----------------------------------------------------------------------
# Turns out the scripts are a lot more flexible if you don't rely on 
# the queuing-system-specific variables -- so I am converting them to
# 'generic' names and using the generics throughout the remainder.
#----------------------------------------------------------------------

if ($?LSB_HOSTS) then

   setenv ORIGINALDIR $LS_SUBCWD
   setenv JOBNAME     $LSB_JOBNAME
   setenv JOBID       $LSB_JOBID
   setenv MYQUEUE     $LSB_QUEUE
   setenv MYHOST      $LSB_SUB_HOST

else

   #-------------------------------------------------------------------
   # You can run this interactively to check syntax, file motion, etc.
   #-------------------------------------------------------------------

   setenv ORIGINALDIR `pwd`
   setenv JOBNAME     jules_perfect
   setenv JOBID       $$
   setenv MYQUEUE     Interactive
   setenv MYHOST      `hostname`

endif

#----------------------------------------------------------------------
# Just an echo of job attributes
#----------------------------------------------------------------------

echo
echo "${JOBNAME} ($JOBID) submitted   from $ORIGINALDIR"
echo "${JOBNAME} ($JOBID) submitted   from $MYHOST"
echo "${JOBNAME} ($JOBID) running in queue $MYQUEUE"
echo "${JOBNAME} ($JOBID) running       on $MYHOST"
echo "${JOBNAME} ($JOBID) started   at "`date`
echo

#----------------------------------------------------------------------
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#
# Set variables containing various directory names where we will GET things:
# DARTDIR      The location of the DART jules model directory
# JULESDIR     The location of the JULES executable
# ENSEMBLEDIR  The location of the initial ensemble of JULES files
# BASEOBSDIR   The directory containing the (empty) observation sequence file.
# CENTRALDIR   The run-time location for the experiment.
#----------------------------------------------------------------------

set nonomatch # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following

switch ("`hostname`")

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set     DARTDIR = /glade/p/work/${USER}/DART/jules/models/jules
      set    JULESDIR = /glade/p/work/${USER}/DART/jules/models/jules/work
      set ENSEMBLEDIR = /glade/p/image/RDA_strawman/JULES_ensembles
      set  BASEOBSDIR = /glade/p/image/Observations/land/pmo
      set  CENTRALDIR = /glade/scratch/${user}/DART/${JOBNAME}/job_${JOBID}
   breaksw

   default:
      # Bristol "lorax"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set     DARTDIR = /users/ar15645/DART_JULES_SVN/models/jules
      set    JULESDIR = /users/hydroeng/JULES/jules-vn4.2/build/bin
      set ENSEMBLEDIR = /users/ar15645/coupling_simulations/check_restart_files/namelist 
      set  BASEOBSDIR = /users/ar15645/DART_JULES_SVN/models/jules/work
      set  CENTRALDIR = /users/ar15645/run_dart_experiment
   breaksw
endsw

#----------------------------------------------------------------------
# Make a unique, (empty, clean) temporary directory for the experiment.
# DART literature frequently calls this 'CENTRALDIR'
#----------------------------------------------------------------------

mkdir -p ${CENTRALDIR}
cd ${CENTRALDIR}

set myname = $0          # this is the name of this script

#-----------------------------------------------------------------------------
# Get the DART executables, scripts, and input files
# Get the JULES executable, control files, and data files.
# The advance_model.csh needs dart_to_jules (but this script does not).
#-----------------------------------------------------------------------------

${COPY} ${DARTDIR}/work/input.nml   input.nml          || exit 1
${COPY} ${DARTDIR}/work/jules_to_dart                . || exit 1
${COPY} ${DARTDIR}/work/dart_to_jules                . || exit 1
${COPY} ${DARTDIR}/work/perfect_model_obs            . || exit 1
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh   . || exit 1

${COPY} ${JULESDIR}/jules.exe                jules.exe || exit 1
${COPY} ${ENSEMBLEDIR}/*nml                          . || exit 1
${COPY} ${ENSEMBLEDIR}/data/tile_fractions.dat       . || exit 1

#-----------------------------------------------------------------------------
# Get the empty observation sequence file ... or die right away.
# This file will dictate the length of the JULES forecast.
#-----------------------------------------------------------------------------

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
# &perfect_model_obs_nml:  restart_in_file_name    = 'dart_ics'
# &perfect_model_obs_nml:  obs_sequence_in_name    = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name   = 'obs_seq.perfect'
# &perfect_model_obs_nml:  init_time_days          = -1,
# &perfect_model_obs_nml:  init_time_seconds       = -1,
# &perfect_model_obs_nml:  first_obs_days          = -1,
# &perfect_model_obs_nml:  first_obs_seconds       = -1,
# &perfect_model_obs_nml:  last_obs_days           = -1,
# &perfect_model_obs_nml:  last_obs_seconds        = -1,
# &jules_to_dart_nml:      jules_to_dart_output_file = 'dart_ics'
# &model_nml:              jules_restart_filename    = 'jules_restart.nc'
# &model_nml:              jules_history_filename    = 'jules_history.nc'
#=========================================================================

#sed -e "s#dart_ics#perfect_ics#" < input.nml.original >! input.nml

#=========================================================================
# Block 2: Convert 1 jules restart file to a DART initial conditions file.
# At the end of the block, we have a DART initial condition file  perfect_ics
#=========================================================================

# because pmo does not update the JULES state, we can simply link.
# filter will modify the file, so it must be copied. 
${LINK} ${ENSEMBLEDIR}/output/day0-1.dump.20140101.0.nc                   . || exit 1
${LINK} ${ENSEMBLEDIR}/output/day0-1.hour.nc                              . || exit 1

# the advance_model.csh script needs to have an instance number in the filename.
# model_mod:static_init_model() cannot have an instance number in the filename.
# So - we just link the two for this part.

${LINK} day0-1.dump.20140101.0.nc    jules_restart.nc
${LINK} day0-1.hour.nc               jules_output.nc  

echo "`date` -- BEGIN JULES-TO-DART"

./jules_to_dart

if ($status != 0) then
   echo "ERROR ... DART died in 'jules_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'jules_to_dart' ... ERROR"
   exit -3
endif

echo "`date` -- END JULES-TO-DART"

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#=========================================================================

echo "`date` -- BEGIN JULES PERFECT_MODEL_OBS"

./perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

#${MOVE} True_State.nc    ../jules_True_State.${LND_DATE_EXT}.nc
#${MOVE} obs_seq.perfect  ../jules_obs_seq.${LND_DATE_EXT}.perfect
#${MOVE} dart_log.out     ../jules_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   jules PERFECT_MODEL_OBS"

#=========================================================================
# Block 4: Update the jules restart file
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} perfect_ics dart_log.nml

echo "`date` -- END   GENERATE jules TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

