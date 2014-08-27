#!/bin/csh
# ex: setup_filter.csh forceCopyDartBuildScripts forceCopyDartScripts forceCopyModelParams
#==============================================================================
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# j.mccreight: july 2014
#
# PURPOSE: 
#  The assimilation file structure for wrfHydro is shown below. 
#  This script:
#    1) check for required directories
#    2) stages parameter files given model choice.
#    3) stages DART scripts
#    4) stages DART executable builds
#    5) check for initial restart files and convert to dart format

#  This script does *NOT* check the namelist paths against the file structure. The
#    model should do this.
#  One can choose to "force copy"  for steps 2-4 by passing the keywords forceCopyDartBuilds,
#    or forceCopyDartScripts.
#  The script determines (from the namelist) if noah or noahMP and always destroys and then 
#    creates the directory PARAMS.gathered to be used by the dart model runs.

# wrfHydro assimilation file structure diagram:
# Notes: 
#   "xyz" means a variable name denoting some interesting case and the 
#     potential for multiple of these.
#   -> indicates symlink
#   centralDir is a "DART" dir in a location-centric wrfHydro dir, it needs to have 
#     everything above it in place as shown below, except for RUN.xyz/ and OBSERVATIONS/. 

# location (this is typically defined by a location/domain for wrfHydro)
#      |
#      + DOMAIN.xyz/ 
#      + FORCING.xyz/
#      + PARAMS.noah.xyz/
#      + PARAMS.noahMP.xyz/
#      + PARAMS.hydro.xyz/
#      + OBSERVATIONS/

#      + RUN.xyz/  (example of a non-dart run directory setup for reference)    
#               |
#               + wrf_hydro.exe (-> e.g. wrf_hydro_NoahMP.parallel.exe)   
#               + *.TBL  -> can be symlinks to ../PARAMS.* or copies.
#               + namelist.hrldas.xyz
#               + namelist.hrldas -> namelist.hrldas.xyz
#               + hydro.namelist.xyz
#               + hydro.namelist -> hydro.namelist.xyz
#               + LOGS/
#               + OUTPUT/  this can only be specified for noah, not MP nor hydro
#               + job.xyz.csh

#      + DART.xyz = $centralDir
#               | (The symlinks are only suggestions)
#               +(r) wrf_hydro.exe (-> e.g. wrf_hydro_NoahMP.parallel.exe)   
#               +(r) input.nml -> input.nml.xyz  (these *.xyz are likely all in this folder too)
#               +(r) obs_seq.out -> OBSERVATIONS/obs_seq.out.xyz
#               +(r) namelist.hrldas -> namelist.hrldas.xyz
#               +(r) hydro.namelist -> hydro.namelist.xyz
#               +(r) DOMAIN/ ->. ../DOMAIN.xyz
#               +(r) PARAMS.noah/ -> ../PARAMS.noah.xyz         (either this line required
#               +(r) PARAMS.noahMp/ -> ../PARAMS.noahMP.xyz      or this line, not both)
#               +(r) PARAMS.hydro/ -> ../PARAMS.hydro.xyz
#               +(r) FORCING/ -> ../FORCING.xyz
#                      |
#                      + yyyymmddhh.LDASIN_DOMAIN*

#               +(r) initialEnsemble/ (for time zero)
#                      |
#                      + restart.iEns.nc        (iEns is a 4 character integer, e.g. 0010)
#                      + restart.hydro.iEns.nc
#                      +(o) restart.assimOnly.iEns.nc

#               +(r) OUTPUT/
#                      |
#                      +(a) model_integration.yyyymmddhh-yyyymmddhh.iEns
#                              |
#                              +(a) RESTART.yyyymmddhh.DOMAIN*.iEns.nc
#                              +(a) HYDRO_RST.yyyy-mm-dd_hh:00_DOMAIN*.iEns.nc
#                              +(o) RESTART_assimOnly.yyyymmddhh.iEns.nc
#                              +(a) *.LDASOUT_DOMAIN*.iEns.nc
#                              +(a) *.LSMOUT_DOMAIN*.iEns
#                              +(a) *.RTOUT_DOMAIN*.iEns
#                              +(a) *.CHRTOUT*.iEns
#                              +(a) *.CHANOBS*.iEns
#                              +(a) frxst_pts_out.txt.iEns
#                              +(a) qstrmvol*.iEns
#                              +(a) GW_*.txt

#               +(o) FORCING.perturbed/  (required if perturbing forcing)
#                      |
#                      + yyyymmddhh.iEns.LDASIN_DOMAIN*

#               +(o) OBSERVATIONS -> ../OBSERVATIONS

#               +(c) PARAMS.gathered   (created by this routine, pulling from the other PARAMS.*)
#               +(c) dart_to_wrfHydro  (DART executable builds)
#               +(c) wrfHydro_to_dart
#               +(c) filter
#               +(c) restart_file_tool
#               +(c) run_filter.csh    (scripts)
#               +(c) advance_model.csh

#               +(ca) restart.nc
#               +(ca) restart.hydro.nc
#               +(cao) restart.assimOnly.nc
#               +(ca) restart.iEns.nc
#               +(ca) restart.hydro.iEns.nc
#               +(cao) restart.assimOnly.iEns.nc
#               +(ca) filter_ics.iEns

#               +(x) advance_temp.pid/  (where the actual model runs are performed)
#                      |
#                      + *.TBL  -> can be symlinks to ../../PARAMS.* or copies.

#                      + restart.nc -> ../restart.iEns.nc -> (see this filebelow) (for current iEns)
#                      + restart.hydro.nc -> ../restart.hydro.iEns.nc  -> (see this file below)
#                      + restart.assimOnly.nc -> ../restart.assimOnly.iEns.nc  -> (see file below)

#               | (output spew of restart files and output files)
#               | (these symlinks go to advance_temp.pid, impyling that model restarts are 
#               |  not sacred and will be altered by DART. Analysis increments stored elsewhere?)
#               |  for integDir=../OUTPUT/model_integration.yyyymmddhh-yyyymmddhh.iEns
#               +(a) restart.iEns.nc -> latest integDir/RESTART.yyyymmddhh.DOMAIN*.iEns.nc 
#               +(a) restart.hydro.iEns.nc -> 
#                               latest integDir/HYDRO_RST.yyyy-mm-dd_hh:00_DOMAIN*.iEns.nc
#               +(a) restart.assimOnly.iEns.nc -> 
#                               latest integDir/RESTART.yyyymmddhh.iEns.nc

#               

# (r) = required, (o) = optional, (c) = created/managed by this routine,
# (a) = created in the assimilation,
# (x) = created/destroyed by dart, must not exist after the run or it's an error (advance_model.csh)

# iEns = 4 digit integer indicating the ith ensemble member
#
# The above file structure implies the following for both RUN and DART runs because 
#   advance_temp.pid is where the dart runs actually take place.
# Some of this stuff may be altered within the namelist copies in advance_temp.pid/
# as necessary, see advance_model.csh.
# ---- namelist.hrldas ----
# HRLDAS_CONSTANTS_FILE = "../DOMAIN.xyz/wrfinput_xyz"
# INDIR = "../FORCING.xyz"  (Altered for perturbed forcing by advance_model.csh)
# OUTDIR = "./" - can only be specified for noah, imples files have to be moved to OUTPUT if desired
# RESTART_FILENAME_REQUESTED = "restart.nc"  (required)
# ---- hydro.namelist ----
# GEO_STATIC_FLNM = "../DOMAIN.xyz/geo_em.d03.nc"
# GEO_FINEGRID_FLNM = "../DOMAIN.xyz/Fulldom_hires_hydrofile_ohd_new_basns_w_cal_params_full_domain.nc"
# RESTART_FILE  = 'restart.hydro.nc' (required)
# ---- input.nml ----
# wrfHydro_to_dart_output_file = 'dart_ics'
# dart_to_wrfHydro_input_file = 'dart_restart'

#==============================================================================
## File structure initialization.

# where this script is run and where we look for and run the model
set centralDir = `pwd`
echo 
echo 'centralDir =' $centralDir

## where we look for dart 
set dartDir = /home/jamesmcc/DART/wrfHydro/models/wrfHydro
echo 'dartDir =' $dartDir

#==============================================================================#
# copy this script, why not?
set setupFilterSource = ${dartDir}/shell_scripts/setup_filter.csh 
set setupFilterTarget = setup_filter.csh
if ( -e $setupFilterTarget ) then 
    if ( `diff -q $setupFilterSource $setupFilterTarget` !~ '' ) then 
	cp $setupFilterSource $setupFilterTarget
	echo "Updating this script: TRY AGAIN!"
	exit 1
    endif 
else 
    cp $setupFilterSource $setupFilterTarget
endif 

#==============================================================================
## Name (logical) variables after input arguments, only modifies these defaults 
## if the proper variables are passed. Only need add an "var=0" line for additional input args
## handled this way.
@ forceCopyDartBuilds = 0
@ forceCopyDartScripts = 0
@ forceCopyModelParams = 0
@ forceCopyAll = 0
@ n = 1
while ($n <= $#argv) 
    set $argv[$n] = 1
    @ n++
end

if ( $forceCopyAll ) then 
    @ forceCopyDartBuilds = 1
    @ forceCopyDartScripts = 1 
    @ forceCopyModelParams = 1
endif 

#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'
## mkdir?

#===============================================================================
#Check the overall dir structure


#===============================================================================
# Clean up previous runs
# Destroy all previous output for the sake of clarity.
\rm -rf OUTPUT
mkdir OUTPUT

\rm -rf advance_temp*
\rm -f  assim_model_state_ic.*
\rm -f  assim_model_state_ud.*
\rm -f  dart_log.*
\rm -f  filter_control*
\rm -f  filter_ics.*
\rm -f  restart*.nc
\rm -f  Posterior_Diag.nc
\rm -f  Prior_Diag.nc


#==============================================================================
# Check for required  files and dirs (excluding parameters)
foreach FILE ( obs_seq.out input.nml wrf_hydro.exe namelist.hrldas hydro.namelist DOMAIN FORCING )
    if (! -e ${FILE} ) then
      echo "MISSING file: ${centralDir}/${FILE}"
      exit 2
   endif
end


#===============================================================================
# Enforce some common sense namelist settings. 

## require restart.nc
set MYSTRING = `grep RESTART_FILENAME_REQUESTED namelist.hrldas | grep -v !`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set WRFINPUT = `echo $MYSTRING[$#MYSTRING]`
if ($WRFINPUT !~ 'restart.nc' & $WRFINPUT !~ './restart.nc' ) then
    echo "The restart file (restart.nc) appears improperly specified in namelist.hrldas"
    exit 2
endif 

## require restart.hydro.nc
set MYSTRING = `grep RESTART_FILE hydro.namelist | grep -v !`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set WRFINPUT = `echo $MYSTRING[$#MYSTRING]`
if ($WRFINPUT !~ 'restart.hydro.nc' & $WRFINPUT !~ './restart.hydro.nc' ) then
    echo "The restart file (restart.hydro.nc) appears improperly specified in hydro.namelist"
    exit 2
endif 


# GEO FILES
# handled by symlinks now so none of this is necessary and this is not
# the appropriate place to check. perhaps this code should move to 
# setup_filter.csh as I'm not currently checking the geofiles
# This is how the "wrfintput" file was handled for noah.
#set MYSTRING = `grep HRLDAS_CONSTANTS_FILE namelist.hrldas`
#set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
#set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
#set WRFINPUT = `echo $MYSTRING[$#MYSTRING]`
#\ln -sv $WRFINPUT .
# do the same for geo_* and Full*?


# FORCING 
# Extract the directory containing the forcing files used to create the LDASIN files.
# This script uses logic that requires ALL of the hourly forcing files
# to be concatenated into a single netCDF file that uses 'time' as the unlimited
# dimension. The required time slices are extracted from this single netCDF file
# into the filenames expected by NOAH. Since each ensemble member generally gets
# a unique atmospheric forcing, this helps minimize the number of files. 
#set MYSTRING   = `grep INDIR namelist.hrldas`
#set MYSTRING   = `echo $MYSTRING | sed -e "s#[=,']# #g"`
#set MYSTRING   = `echo $MYSTRING | sed -e 's#"# #g'`
#set FORCINGDIR = `echo $MYSTRING[$#MYSTRING]`
# echo "Master atmospheric netCDF forcing file(s) coming from $FORCINGDIR"



#==============================================================================
# Parameters.
#   are gathered here for the model specfied in the namelist.hrldas
echo 
echo "PARAMS"
set paramsGathered = "PARAMS.gathered"
${REMOVE} $paramsGathered
\mkdir $paramsGathered

# Determine if this is a "noah" or a "noahmp" (lowercase) run based on the input.nml
# in order to copy the appropriate parameter files
# use the not-commented line specifying lsm_model_choice
set lsmModel = `grep -v '!' input.nml | grep -i lsm_model_choice`
set lsmModel = `echo $lsmModel | sed -e "s#[=,':]# #g"`
set lsmModel = `echo $lsmModel[$#lsmModel]`
set lsmModel = `echo $lsmModel | tr '[:upper:]' '[:lower:]'`

#noah params
## (the following three blocks should be written to use a function..)
if ( "$lsmModel" == "noah" ) then 
    set dir = "PARAMS.noah"
    if (! -e ${dir}) then 
	echo "MISSING directory in working dir: ${dir}/"
	exit 3
    endif 
    #foreach FILE ( GENPARM.TBL  LANDUSE.TBL  SOILPARM.TBL  URBPARM.TBL  VEGPARM.TBL )
    foreach FILE ( GENPARM.TBL  SOILPARM.TBL  URBPARM.TBL  VEGPARM.TBL )
	if (! -e ${dir}/${FILE} )  then
	    echo "MISSING file: ${dir}/${FILE}"
	    exit 4
	endif
	echo "Copying ${dir}/${FILE} to $paramsGathered"
        ${COPY} ${dir}/${FILE} $paramsGathered
    end
endif 

## noahMP params
if ( "$lsmModel" == "noahmp" ) then 
    set dir = "PARAMS.noahMP"
    if (! -e ${dir}) then 
	echo "MISSING directory: ${dir}"
	exit 5
    endif 
    foreach FILE ( DISTR_HYDRO_CAL_PARMS.TBL GENPARM.TBL MPTABLE.TBL SOILPARM.TBL URBPARM.TBL VEGPARM.TBL )
	if (! -e ${dir}/${FILE} )  then
	    echo "MISSING file: ${dir}/${FILE}"
	    exit 6
	endif
	echo "Copying ${dir}/${FILE} to $paramsGathered"
        ${COPY} ${dir}/${FILE} $paramsGathered
    end
endif 

## hydro params
## (for now we'll assume the hydro component is always used. If it
##  is not being used, this is still probably harmless.)
set dir = "PARAMS.hydro"
if (! -e ${dir}) then 
    echo "MISSING directory: ${dir}"
    exit 7
endif 
foreach FILE ( CHANPARM.TBL GWBUCKPARM.TBL HYDRO.TBL LAKEPARM.TBL )
    if (! -e ${dir}/${FILE} )  then
	echo "MISSING file: ${dir}/${FILE}"
	exit 7
    endif
    echo "Copying ${dir}/${FILE} to $paramsGathered"
    ${COPY} ${dir}/${FILE} $paramsGathered
end 

#==============================================================================
echo
echo "DART Executable Builds"
foreach FILE ( dart_to_wrfHydro wrfHydro_to_dart filter restart_file_tool )
    if ( -e ${FILE} && ! $forceCopyDartBuilds )  then
	echo "Using existing $FILE"
    else
	if(! -e ${dartDir}/work/${FILE} ) then 
	    echo "MISSING file: ${dartDir}/work/${FILE}"
	    exit 5
	endif
	${COPY} ${dartDir}/work/${FILE} . 
   endif
end

#==============================================================================
echo
echo "DART Scripts"
foreach FILE ( run_filter.csh advance_model.csh )
   if ( -e ${FILE} && ! $forceCopyDartScripts )  then
      echo "Using existing $FILE"
   else
	if(! -e ${dartDir}/shell_scripts/${FILE}) then
	    echo "MISSING file: ${dartDir}/shell_scripts/${FILE}"
	    exit 6
	endif 
	${COPY} ${dartDir}/shell_scripts/${FILE} .
   endif
end

#==============================================================================
# Initial ensemble
# 1) point to a directory full of noah restart files and pick N of them;
#    The tricky part is that the Time variable in those files is all wrong.
# 2) convert them to DART format;
# 3) make sure the time in each DART file is 'identical'
# 4) the original noah restart files are also needed to start up
#    advance_model.csh for the first model advance.
#

# link to the useful script for creating initial ensembles. 
ln -svf ${dartDir}/shell_scripts/mk_initial_ens.csh .
ln -svf ${dartDir}/shell_scripts/randomSeq.rsh .
ln -svf ${dartDir}/shell_scripts/run_initial_ens_fwd.csh .

# where the ensemble of initial conditions are located
set ensembleDir = ${centralDir}/initialEnsemble
echo
echo 'ensembleDir =' $ensembleDir

## if ensembleDir does not exist, create the directory but complain.
if (! -e ${ensembleDir} ) then 
    echo
    echo "CREATING $ensembleDir"
    \mkdir $ensembleDir
    echo "ENSEMBLE initial condition files have not been established (in the proper location)."
    echo "Please populate or link:  ${centralDir}/${ensembleDir}"
    echo "mk_initial_ens.csh might help!"
    echo "run_initial_ens_fwd.csh might also help!"
    exit 1
endif 

set nLsmFiles = `ls -1 ${ensembleDir}/restart.[^ha]* | wc -l`
set lsmFileList = `ls -1 ${ensembleDir}/restart.[^ha]*`

set   nHydroFiles = `ls -1 ${ensembleDir}/*hydro* | wc -l`
set hydroFileList = `ls -1 ${ensembleDir}/*hydro*`

echo $nLsmFiles
echo $nHydroFiles

if ( ! $nLsmFiles | ! $nHydroFiles) then
    echo "No initial restart files provided"
    echo "mk_initial_ens.csh might help!"
    echo "run_initial_ens_fwd.csh might also help!"
    exit 1
endif 

## From the namelist determine if the noAssim restarts are needed. 
## The line could be commented out (default is blank in model_mod.f90) or set to ''.
set assimOnly_active1 = `grep -v '!' input.nml | grep -i assimOnly_netcdf_filename | wc -l`
set assimOnly_active2 = `grep -v '!' input.nml | grep -i assimOnly_netcdf_filename | \
                         cut -d= -f2 | tr -cd '[[:alnum:]]._-' | wc -m`
set assimOnly_active = 0
if ($assimOnly_active1 && $assimOnly_active2) set assimOnly_active = 1

if ($assimOnly_active) then 
    set   nAssimOnlyFiles = `ls -1 ${ensembleDir}/*assimOnly* | wc -l`
    set assimOnlyFileList = `ls -1 ${ensembleDir}/*assimOnly*`
    echo $nAssimOnlyFiles
endif 

if ( $nLsmFiles !~ $nHydroFiles ) then
    echo "The number of LSM and HYDRO restart files do NOT match in $ensembleDir"
    exit 8
endif 

if ($assimOnly_active) then
    if ( $nLsmFiles !~ $nAssimOnlyFiles ) then
	echo "The number of LSM and assimOnly restart files do NOT match in $ensembleDir"
	exit 8
    endif
endif 

echo "The inital emsemble restart files:"
echo $lsmFileList
echo $hydroFileList
if ($assimOnly_active) then 
    echo $assimOnlyFileList
endif 

#===============================================================================
# Set the date in the namelist.hrldas to match that of the first restart file.
# (because there's no error generate - only bogus output where there is no forcing.).
set restartFileTime = `ncdump -v Times $ensembleDir/restart.0001.nc | tail -2 | head -1 | cut -d'"' -f2`
set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
set restartFileMm = `echo $restartFileTime | cut -d- -f2`
set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`

ex namelist.hrldas <<ex_end
g;START_YEAR;s;= .*;= $restartFileYyyy;
g;START_MONTH;s;= .*;= $restartFileMm;
g;START_DAY;s;= .*;= $restartFileDd;
g;START_HOUR;s;= .*;= $restartFileHh;
wq
ex_end

#===============================================================================
# Convert initial condition restart files to to DART initial condition files.
# This is the last step in advance_model, so these are expected at start of filter.
# Also move the inital restarts to the OUTPUT directory, as they will be updated
# by dart_to_wrfHydro
set initEnsOutDir = OUTPUT/initial_restarts
\mkdir -p $initEnsOutDir
@ ifile = 1
@ ensemble_member = 0
while ($ifile <= $nLsmFiles)
    @ ensemble_member = $ensemble_member + 1
    set fext = `printf %04d $ensemble_member`
    
    # make initial conditions for DART
    \ln -svf ${ensembleDir}/restart.$fext.nc restart.nc
    \ln -svf ${ensembleDir}/restart.hydro.$fext.nc restart.hydro.nc
    if ( $assimOnly_active ) \
	\ln -svf ${ensembleDir}/restart.assimOnly.$fext.nc restart.assimOnly.nc

    ./wrfHydro_to_dart                     || exit 9
    ${MOVE} dart_ics filter_ics.$fext      || exit 10 #??where are the naming assumptions explained??
    #${MOVE} dart_ics assim_model_state_ic.$fext      || exit 10
    
    \cp ${ensembleDir}/restart.$fext.nc $initEnsOutDir/restart.$fext.nc 
    \ln -sf  $initEnsOutDir/restart.$fext.nc .
    \cp ${ensembleDir}/restart.hydro.$fext.nc $initEnsOutDir/restart.hydro.$fext.nc 
    \ln -sf  $initEnsOutDir/restart.hydro.$fext.nc .
    if ( $assimOnly_active ) then 
	\cp ${ensembleDir}/restart.assimOnly.$fext.nc $initEnsOutDir/restart.assimOnly.$fext.nc 
	\ln -sf  $initEnsOutDir/restart.assimOnly.$fext.nc .
    endif

   @ ifile = $ifile + 1
end
# The left over restart.nc and restart.hydro.nc are used to initialize the model.

# Since we have some knowledge of the ensemble size, 
# insert reasonable default values in input.nml.

ex input.nml <<ex_end
/filter
g;ens_size ;s;= .*;= $ensemble_member,;
g;num_output_state_members ;s;= .*;= $ensemble_member,;
g;num_output_obs_members ;s;= .*;= $ensemble_member,;
g;single_restart_file_in ;s;= .*;= .false.,;
g;single_restart_file_out ;s;= .*;= .false.,;
wq
ex_end

# DART needs a copy of the NOAH restart file to determine the sizes
# of the state vector components. We are going to LEAVE the final
# ensemble member restart file linked to the expected restart file name.



#==============================================================================
# Finish up.

# Would like to output put the resulting file structure to a log file. 
# set logFile = "setup_filter_log.txt"
#echo "setup_filter.csh" > ${logFile}
#date >> ${logFile}
## This can be massive for forcing dirs with lots of files. I need to make it exclude the forcing dir
## or do something else intelligent for forcing files.
#${dartDir}/shell_scripts/listDirFollowLinks.sh .. >> ${logFile}

echo
echo "centralDir is ${centralDir}"
echo "Configure     ${centralDir}/input.nml"
echo "Configure     ${centralDir}/namelist.hrldas"
echo "Configure     ${centralDir}/hydro.namelist"
echo "execute       ${centralDir}/run_filter.csh"
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

