#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
#=========================================================================

echo "`date` -- BEGIN GENERATE CLM TRUE STATE"
pwd

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# As of CESM2.0, the perfect_model.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE

setenv CASEROOT $1
setenv ASSIMILATION_CYCLE $2

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ASSIMILATION_CYCLE + 1

# xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}
setenv CASE           `./xmlquery CASE        --value`
setenv ENSEMBLE_SIZE  `./xmlquery NINST_LND   --value`
setenv EXEROOT        `./xmlquery EXEROOT     --value`
setenv RUNDIR         `./xmlquery RUNDIR      --value`
setenv ARCHIVE        `./xmlquery DOUT_S_ROOT --value`
setenv TOTALPES       `./xmlquery TOTALPES    --value`
setenv STOP_N         `./xmlquery STOP_N      --value`
setenv DATA_ASSIMILATION_CYCLES `./xmlquery DATA_ASSIMILATION_CYCLES --value`
setenv TASKS_PER_NODE `./xmlquery MAX_TASKS_PER_NODE --value`
cd ${RUNDIR}

# string to be replaced by the setup script or by hand once
set BASEOBSDIR = BOGUSBASEOBSDIR

# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      setenv MP_DEBUG_NOTIMEOUT yes

      set  LAUNCHCMD = mpirun.lsf
   breaksw

   case lone*:
      # UT lonestar
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

      set  LAUNCHCMD = mpirun.lsf
   breaksw

   case la*:
      # LBNL "lawrencium"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set  LAUNCHCMD = "mpiexec -n $NTASKS"
   breaksw

   case r*:
      # NCAR "cheyenne"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      setenv MP_DEBUG_NOTIMEOUT yes
      set  LAUNCHCMD = mpiexec_mpt
   breaksw

   case example*:
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      set LAUNCHCMD = "aprun -n $NTASKS"
   breaksw

   default:
      echo "FATAL ERROR: The system-specific environment must be specified."   
      echo "             Add system-specific info to a case statement in"   
      echo "             ${CASEROOT}/perfect_model.csh"   
      exit 1
   breaksw
endsw

#=========================================================================
# Block 1: Determine time of model state ... from file name
# of the form "./${CASE}.clm2.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#=========================================================================

set FILE = `head -n 1 rpointer.lnd`
set FILE = $FILE:r
set LND_DATE_EXT = `echo $FILE:e`
set LND_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set LND_YEAR     = `echo $LND_DATE[1] | bc`
set LND_MONTH    = `echo $LND_DATE[2] | bc`
set LND_DAY      = `echo $LND_DATE[3] | bc`
set LND_SECONDS  = `echo $LND_DATE[4] | bc`
set LND_HOUR     = `echo $LND_DATE[4] / 3600 | bc`

echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_SECONDS (seconds)"
echo "valid time of model is $LND_YEAR $LND_MONTH $LND_DAY $LND_HOUR (hours)"

#=========================================================================
# Block 2: Populate a run-time directory with the input needed to run DART.
#=========================================================================

# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CLM.
#
# The CLM observations are stored in two sets of directories.
# If you are stopping every 24 hours or more, the obs are in directories like YYYYMM.
# In all other situations the observations come from directories like YYYYMM_6H.
# The only ugly part here is if the first advance and subsequent advances are
# not the same length. The observations _may_ come from different directories.
#
# The contents of the file must match the history file contents if one is using
# the obs_def_tower_mod or could be the 'traditional' +/- 12Z ... or both.
# Since the history file contains the previous days' history ... so must the obs file.

if ($STOP_N >= 24) then
   set OBSDIR = `printf %04d%02d    ${LND_YEAR} ${LND_MONTH}`
else
   set OBSDIR = `printf %04d%02d_6H ${LND_YEAR} ${LND_MONTH}`
endif

set OBS_FILE = ${BASEOBSDIR}/${OBSDIR}/obs_seq.${LND_DATE_EXT}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

# Enforce the namelist assumptions used in this script.

ex input.nml <<ex_end
g;read_input_state_from_file ;s;= .*;= .true.;
g;write_output_state_to_file ;s;= .*;= .false.;
g;single_file_in ;s;= .*;= .false.;
g;single_file_out ;s;= .*;= .false.;
g;input_state_files ;s;= .*;= 'clm_restart.nc', 'clm_history.nc', 'clm_vector_history.nc';
g;output_state_files ;s;= .*;= 'clm_restart.nc', 'clm_history.nc', 'clm_vector_history.nc';
g;input_state_file_list ;s;= .*;= '', '';
g;output_state_file_list ;s;= .*;= '', '';
g;obs_seq_in_file_name ;s;= .*;= 'obs_seq.in';
g;obs_seq_out_file_name ;s;= .*;= 'obs_seq.out';
wq
ex_end

# DART/CLM routines all need a clm_restart.nc, clm_history.nc, etc.
# The flux tower forward operator looks for a CLM history file with
# an instance number in the filename.

set      LND_RESTART_FILENAME = ${CASE}.clm2.r.${LND_DATE_EXT}.nc
set      LND_HISTORY_FILENAME = ${CASE}.clm2.h0.${LND_DATE_EXT}.nc
set  LND_VEC_HISTORY_FILENAME = ${CASE}.clm2.h2.${LND_DATE_EXT}.nc
set     OBS1_HISTORY_FILENAME = ${CASE}.clm2.h1.${LND_DATE_EXT}.nc
set     OBS2_HISTORY_FILENAME = ${CASE}.clm2_0001.h1.${LND_DATE_EXT}.nc

${LINK} $LND_RESTART_FILENAME clm_restart.nc
${LINK} $LND_HISTORY_FILENAME clm_history.nc

if (  -e   $OBS1_HISTORY_FILENAME) then
   ${LINK} $OBS1_HISTORY_FILENAME $OBS2_HISTORY_FILENAME
endif
if (  -e   $LND_VEC_HISTORY_FILENAME) then
   ${LINK} $LND_VEC_HISTORY_FILENAME clm_vector_history.nc
endif

#=========================================================================
# Block 3: Advance the model and harvest the synthetic observations.
# output files are:
# True_state.nc   ...... the DART state
# obs_seq.perfect ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
# perfect_restart ...... which we don't need
#=========================================================================

echo "`date` -- BEGIN CLM PERFECT_MODEL_OBS"

${LAUNCHCMD} ${EXEROOT}/perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit -4
endif

# ${MOVE} True_State.nc    clm_True_State.${LND_DATE_EXT}.nc
# ${MOVE} obs_seq.perfect  clm_obs_seq.${LND_DATE_EXT}.perfect
 ${MOVE} dart_log.out     clm_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   CLM PERFECT_MODEL_OBS"

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} perfect_ics dart_log.nml

echo "`date` -- END   GENERATE CLM TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

