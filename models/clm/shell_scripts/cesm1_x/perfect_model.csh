#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN GENERATE CLM TRUE STATE"

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -v --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -vs'
      set REMOVE = '/usr/local/bin/rm -fr'

      set BASEOBSDIR = /glade/proj3/image/Observations/FluxTower
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -v --preserve=timestamps'
      set   LINK = 'ln -vs'
      set REMOVE = 'rm -fr'
      set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
      setenv MP_DEBUG_NOTIMEOUT yes

      set BASEOBSDIR = /glade/p/image/Observations/land
   breaksw

   case lone*:
      # UT lonestar
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -v --preserve=timestamps'
      set   LINK = '/bin/ln -vs'
      set REMOVE = '/bin/rm -fr'

      set BASEOBSDIR = ${WORK}/DART/observations/snow/work/obs_seqs
   breaksw

   case la*:
      # LBNL "lawrencium"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -v --preserve=timestamps'
      set   LINK = 'ln -vs'
      set REMOVE = 'rm -fr'
      set TASKS_PER_NODE = $MAX_TASKS_PER_NODE

      set BASEOBSDIR = /your/observation/directory/here
   breaksw

   default:
      # all others
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -v --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
   breaksw
endsw

#-------------------------------------------------------------------------
# Block 1: Determine time of model state ... from file name
# of the form "./${CASE}.clm2.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

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

#-----------------------------------------------------------------------------
# Block 2: Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CLM.
#
# The CLM observations are stowed in two sets of directories.
# If you are stopping every 24 hours or more, the obs are in directories like YYYYMM.
# In all other situations the observations come from directories like YYYYMM_6H.
# The only ugly part here is if the first advance and subsequent advances are
# not the same length. The observations _may_ come from different directories.
#
# The contents of the file must match the history file contents if one is using
# the obs_def_tower_mod or could be the 'traditional' +/- 12Z ... or both.
# Since the history file contains the previous days' history ... so must the obs file.
#-----------------------------------------------------------------------------

if ($STOP_N >= 24) then
   set OBSDIR = `printf %04d%02d    ${LND_YEAR} ${LND_MONTH}`
else
   set OBSDIR = `printf %04d%02d_6H ${LND_YEAR} ${LND_MONTH}`
endif

set OBS_FILE = ${BASEOBSDIR}/${OBSDIR}/obs_seq.${LND_DATE_EXT}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in || exit 2
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit 2
endif

#-----------------------------------------------------------------------------
# Block 3: Populate a run-time directory with the input needed to run DART.
#
# &model_nml: clm_restart_filename        = 'clm_restart.nc'
# &model_nml: clm_history_filename        = 'clm_history.nc'
# &model_nml: clm_vector_history_filename = 'clm_vector_history.nc' [OPTIONAL]
#-----------------------------------------------------------------------------

${REMOVE} input.nml

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml . || exit 3
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit 3
endif

# DART/CLM routines all need a clm_restart.nc, clm_history.nc, etc.
# The flux tower forward operator looks for a CLM history file with
# an instance number in the filename.

set      LND_RESTART_FILENAME = ${CASE}.clm2.r.${LND_DATE_EXT}.nc
set      LND_HISTORY_FILENAME = ${CASE}.clm2.h0.${LND_DATE_EXT}.nc
set  LND_VEC_HISTORY_FILENAME = ${CASE}.clm2.h2.${LND_DATE_EXT}.nc
set     OBS1_HISTORY_FILENAME = ${CASE}.clm2.h1.${LND_DATE_EXT}.nc
set     OBS2_HISTORY_FILENAME = ${CASE}.clm2_0001.h1.${LND_DATE_EXT}.nc

# remove any potentially pre-existing linked files 
${REMOVE} clm_restart.nc clm_history.nc clm_vector_history.nc ${OBS2_HISTORY_FILENAME}

${LINK} $LND_RESTART_FILENAME clm_restart.nc || exit 3
${LINK} $LND_HISTORY_FILENAME clm_history.nc || exit 3

if (  -e   $OBS1_HISTORY_FILENAME) then
   ${LINK} $OBS1_HISTORY_FILENAME $OBS2_HISTORY_FILENAME || exit 3
endif
if (  -e   $LND_VEC_HISTORY_FILENAME) then
   ${LINK} $LND_VEC_HISTORY_FILENAME clm_vector_history.nc || exit 3
endif

#-----------------------------------------------------------------------------
# Block 4: harvest the synthetic observations.
# output files are:
# true_state.nc   ...... the DART state
# obs_seq.out     ...... the synthetic observations
# dart_log.out    ...... run-time output of all DART routines
#-----------------------------------------------------------------------------

echo "`date` -- BEGIN CLM PERFECT_MODEL_OBS"

${MPI_RUN_COMMAND} ${EXEROOT}/perfect_model_obs

if ($status != 0) then
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   echo "ERROR ... DART died in 'perfect_model_obs' ... ERROR"
   exit 4
endif

${MOVE} true_state.nc    clm_true_state.${LND_DATE_EXT}.nc
${MOVE} obs_seq.out      clm_obs_seq.${LND_DATE_EXT}.perfect
${MOVE} dart_log.out     clm_dart_log.${LND_DATE_EXT}.out

echo "`date` -- END   CLM PERFECT_MODEL_OBS"

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

${REMOVE} dart_log.nml

echo "`date` -- END   GENERATE CLM TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

