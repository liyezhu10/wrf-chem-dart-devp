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

# As of CESM2.0, the assimilate.csh is called by CESM - and has
# two arguments: the CASEROOT and the DATA_ASSIMILATION_CYCLE

echo "`date` -- BEGIN CLM_ASSIMILATE"
pwd

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

setenv CASEROOT $1
setenv ASSIMILATION_CYCLE $2

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ASSIMILATION_CYCLE + 1

# xmlquery must be executed in $CASEROOT.
cd ${CASEROOT}
setenv CASE           `./xmlquery CASE        --value`
setenv ensemble_size  `./xmlquery NINST_LND   --value`
setenv EXEROOT        `./xmlquery EXEROOT     --value`
setenv RUNDIR         `./xmlquery RUNDIR      --value`
setenv archive        `./xmlquery DOUT_S_ROOT --value`
setenv TOTALPES       `./xmlquery TOTALPES    --value`
setenv DATA_ASSIMILATION_CYCLES `./xmlquery DATA_ASSIMILATION_CYCLES --value`
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
      set TASKS_PER_NODE = `echo $LSB_SUB_RES_REQ | sed -ne '/ptile/s#.*\[ptile=\([0-9][0-9]*\)]#\1#p'`
      setenv MP_DEBUG_NOTIMEOUT yes

      set  LAUNCHCMD = mpirun.lsf
   breaksw

   case lone*:
      # UT lonestar
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

      set BASEOBSDIR = ${WORK}/DART/observations/snow/work/obs_seqs
      set  LAUNCHCMD = mpirun.lsf
   breaksw

   case la*:
      # LBNL "lawrencium"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      set TASKS_PER_NODE = $MAX_TASKS_PER_NODE

      set BASEOBSDIR = /your/observation/directory/here
      set  LAUNCHCMD = "mpiexec -n $NTASKS"
   breaksw

   case ch*:
      # NCAR "cheyenne"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'
      # echo "Trying to set TASKS_PER_NODE using PBS_NUM_PPN $PBS_NUM_PPN"
      # Unavailable for some reason:  set TASKS_PER_NODE = $PBS_NUM_PPN
      set TASKS_PER_NODE = 36
      setenv MP_DEBUG_NOTIMEOUT yes
      set  LAUNCHCMD = mpiexec_mpt
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
      set  LAUNCHCMD = "aprun -n $NTASKS"
   breaksw
endsw

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.lnd_0001`
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
# Get observation sequence file ... or die right away.
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
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${COPY} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" input.nml.$$ >! input.nml
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: the lookup table for SAMPLING ERROR CORRECTION was put in the
# RUNDIR by the DART_config script - nothing to do now.
#=========================================================================

#=========================================================================
# Block 3: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
# inf_out_file_name           = 'prior_inflate_restart', 'post_inflate_restart',
# inf_diag_file_name          = 'prior_inflate_diag',    'post_inflate_diag',
#
# NOTICE: the archiving scripts more or less require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
#
# The inflation file is essentially a duplicate of the DART model state ...
# For the purpose of this script, they are the output of a previous assimilation,
# so they should be named something like prior_inflate_restart.YYYY-MM-DD-SSSSS
#
# NOTICE: inf_initial_from_restart and inf_sd_initial_from_restart are somewhat
# problematic. During the bulk of an experiment, these should be TRUE, since
# we want to read existing inflation files. However, the first assimilation
# might need these to be FALSE and then subsequently be set to TRUE.
# There are two ways to handle this:
#
# 1) Create the initial files offline (perhaps with 'fill_inflation_restart')
#    and stage them with the appropriate names in the RUNDIR.
#    You must manually remove the clm_inflation_cookie file
#    from the RUNDIR in this case.
#    - OR -
# 2) create a cookie file called RUNDIR/clm_inflation_cookie
#    The existence of this file will cause this script to set the
#    namelist appropriately. This script will 'eat' the cookie file
#    to prevent this from happening for subsequent executions. If the
#    inflation file does not exist for them, and it needs to, this script
#    should die. The CESM_DART_config script automatically creates a cookie
#    file to support this option.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#=========================================================================

set  MYSTRING = `grep 'inf_flavor' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep 'inf_initial_from_restart' input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
set  POSTE_TF = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      # we are not using an existing inflation file.
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else if ( -e clm_inflation_cookie ) then
      # We want to use an existing inflation file, but this is
      # the first assimilation so there is no existing inflation
      # file. This is the signal we need to to coerce the namelist
      # to have different values for this execution ONLY.
      # Since the local namelist comes from CASEROOT each time, we're golden.

      set PRIOR_TF = FALSE

ex input.nml <<ex_end
g;inf_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
g;inf_sd_initial_from_restart ;s;= .*;= .${PRIOR_TF}., .${POSTE_TF}.,;
wq
ex_end

   else

      # Checking for a prior inflation mean file from the previous assimilation.

      (ls -rt1 clm_output_priorinf_mean.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_mean.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like clm_output_priorinf_mean.YYYY-MM-DD-SSSSS.nc"
         exit -4
      endif

      # Checking for a prior inflation sd file from the previous assimilation.

      (ls -rt1 clm_output_priorinf_sd.* | tail -n 1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest input_priorinf_sd.nc
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation SD file."
         echo "ERROR: expected something like clm_output_priorinf_sd.YYYY-MM-DD-SSSSS.nc"
         exit 2
      endif

   endif
else
   echo "Prior Inflation           not requested for this assimilation."
endif

# Eat the cookie regardless
${REMOVE} clm_inflation_cookie

#=========================================================================
# Block 4: Convert N CLM restart files to DART initial condition files.
# clm_to_dart is serial code, we can do all of these at the same time
# as long as we can have unique namelists for each of them.
#
# At the end of the block, we have DART initial condition files  filter_ics.[1-N]
# that came from pointer files ../rpointer.lnd_[1-N]
#
# REQUIRED DART namelist settings:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
#                        restart_out_file_name   = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &clm_to_dart_nml:      clm_to_dart_output_file = 'dart_ics',
# &dart_to_clm_nml:      dart_to_clm_input_file  = 'dart_restart',
#                        advance_time_present    = .false.
# &model_nml:            clm_restart_filename        = 'clm_restart.nc'
# &model_nml:            clm_history_filename        = 'clm_history.nc'
# &model_nml:            clm_vector_history_filename = 'clm_vector_history.nc'
#
#=========================================================================

echo "`date` -- BEGIN CLM-TO-DART"

ls -1 ${CASE}.clm2_*.r.${LND_DATE_EXT}.nc  >! restart_files.txt
ls -1 ${CASE}.clm2_*.h0.${LND_DATE_EXT}.nc >! history_files.txt
ls -1 ${CASE}.clm2_*.h2.${LND_DATE_EXT}.nc >! history_vector_files.txt

echo "`date` -- END CLM-TO-DART for all ${ensemble_size} members."

#=========================================================================
# Block 5: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                   = 0,
# &filter_nml:           adv_ens_command         = "no_advance_script",
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &filter_nml:           restart_out_file_name   = 'filter_restart'
# &filter_nml:           obs_sequence_in_name    = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name   = 'obs_seq.final'
# &filter_nml:           init_time_days          = -1,
# &filter_nml:           init_time_seconds       = -1,
# &filter_nml:           first_obs_days          = -1,
# &filter_nml:           first_obs_seconds       = -1,
# &filter_nml:           last_obs_days           = -1,
# &filter_nml:           last_obs_seconds        = -1,
# &ensemble_manager_nml: single_restart_file_in  = .false.
# &ensemble_manager_nml: single_restart_file_out = .false.
# &model_nml:            clm_restart_filename        = 'clm_restart.nc'
# &model_nml:            clm_history_filename        = 'clm_history.nc'
# &model_nml:            clm_vector_history_filename = 'clm_vector_history.nc'
#=========================================================================

# clm always needs a clm_restart.nc, clm_history.nc for geometry information, etc.
# it may or may not need a vector-format history file - depends on user input

set     LND_RESTART_FILENAME = ${CASE}.clm2_0001.r.${LND_DATE_EXT}.nc
set     LND_HISTORY_FILENAME = ${CASE}.clm2_0001.h0.${LND_DATE_EXT}.nc
set LND_VEC_HISTORY_FILENAME = ${CASE}.clm2_0001.h2.${LND_DATE_EXT}.nc

${LINK} ${LND_RESTART_FILENAME} clm_restart.nc
${LINK} ${LND_HISTORY_FILENAME} clm_history.nc
if (  -e   ${LND_VEC_HISTORY_FILENAME} ) then
   ${LINK} ${LND_VEC_HISTORY_FILENAME} clm_vector_history.nc
endif

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -7
echo "`date` -- END FILTER"

# ${MOVE} preassim.nc      clm_preassim.${LND_DATE_EXT}.nc
# ${MOVE} analysis.nc      clm_analysis.${LND_DATE_EXT}.nc

${MOVE} obs_seq.final    clm_obs_seq.${LND_DATE_EXT}.final
${MOVE} dart_log.out     clm_dart_log.${LND_DATE_EXT}.out

# Accomodate any possible inflation files.
# Rename file to reflect current date.

foreach FILE ( output_priorinf_mean.* output_priorinf_sd.* )
   if ( -e ${FILE} ) then
      set BASE = $FILE:r
      set EXTENSION = $FILE:e
      ${MOVE} ${FILE} clm_${FILE}.${LND_DATE_EXT}.${EXTENSION}
   else
      echo "No ${FILE} for ${LND_DATE_EXT}"
   endif
end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END CLM_ASSIMILATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

