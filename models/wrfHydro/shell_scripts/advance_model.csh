#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.
# Can be used as-is with most low-order models and the bgrid model which
# can be advanced using the integrate_model executable.
# 
# Arguments are (created by 'filter' or 'perfect_model_obs' and include):
# 1) the process number of caller,
# 2) the number of ensemble members/state copies belonging to that process, and 
# 3) the name of the control_file for that process.
# 
# If this script finishes and the 'control_file' still exists, it is
# an ERROR CONDITION and means one or more of the ensemble members did
# not advance properly. Despite our best attempts to trap on this
# condition, some MPI installations simply hang, some properly terminate.
#
# This script loops over all the entries in the control_file to advance 
# any/all of the ensemble members.  The number of trips through the 
# loop is the second argument to this script. The control_file contains 
# the information about which ensemble members are to be advanced by THIS TASK.
# Sometimes it may be just one ensemble member, sometimes all of them.
# Read DART/doc/html/filter_async_modes.html and the mpi_intro.html
# for an overview.
#
# jlm An example layout of control and input files would be helpful. These links 
# dont provide much (quick) insight on that.
#
# This script has 4 logical 'blocks':
# 1) creates a clean, temporary directory in which to run a model instance
#      and copies the necessary files into the temporary directory
# 2) copies/converts the DART state vector to something the model can ingest
#      (dart_to_wrfHydro)
# 3) runs the model (wrfHydro)
# 4) copies/converts the model output to input expected by DART
#      (wrfHydro_to_dart)

set      process = $1
set   num_states = $2
set control_file = $3

#----------------------------------------------------------------------
# Block 1: copy necessary input files/executables/files common
#          to all model advances to a clean, temporary directory.
#          These will be used by ALL of the ensemble
#          members being advanced by this script.
#----------------------------------------------------------------------

# Create a unique temporary working directory for this process's stuff
# The run-time directory for the entire experiment is called CENTRALDIR;
# we need to provide a safe haven for each TASK ... in 'temp_dir'.
set temp_dir = 'advance_temp'$process

# Create a clean temporary directory and go there
\rm -rf  $temp_dir  || exit 1
\mkdir -p $temp_dir  || exit 1
cd       $temp_dir  || exit 1

# Get the DART input.nml and the NOAH namelist
\cp ../namelist.hrldas .
\cp ../hydro.namelist .
\cp ../input.nml .   ## i'm not seeing this used here but it might be used by filter.

# Get parameters.
foreach FILE ( ../PARAMS.gathered/* ) 
    echo $FILE
#   \cp -v ../$FILE . || exit 2 ## these dont change so link them instead.
    \ln -sf $FILE . || exit 2  
end

# From the namelist determine if the noAssim restarts are needed. 
# The line could be commented out (default is blank in model_mod.f90) or set to ''.
set assimOnly_active1 = `grep -v '!' input.nml | grep -i assimOnly_netcdf_filename | wc -l`
set assimOnly_active2 = `grep -v '!' input.nml | grep -i assimOnly_netcdf_filename | \
                         cut -d= -f2 | tr -cd '[[:alnum:]]._-' | wc -m`
set assimOnly_active = 0
if ($assimOnly_active1 && $assimOnly_active2) set assimOnly_active = 1
# assimOnly variables dont necessarily mean perturbed forcing?
# for now the do mean perturbed forcing...
# eventually may want to list the variables in the file and set flags based on these. 

#-------------------------------------------------------------------------------
# Loop through each state
set state_copy = 1
# These reference line number in the control file.
set ensemble_member_line = 1
set      input_file_line = 2
set     output_file_line = 3

while($state_copy <= $num_states)
   
    set ensemble_member = `\head -$ensemble_member_line ../$control_file | \tail -1`
    set input_file      = `\head -$input_file_line      ../$control_file | \tail -1`
    set output_file     = `\head -$output_file_line     ../$control_file | \tail -1`
    set instance        = `printf "%04d" $ensemble_member`

    #-------------------------------------------------------------------
    # Block 2: dart_to_wrfHydro
    #          * remove scraps from previous advances
    #          * copy/link ensemble-member-specific files
    #          * convey the advance-to-time to the model
    #          * dart_to_wrfHydro: convert the DART state vector to model format
    #-------------------------------------------------------------------

    echo "advance_model.csh block 2 converting ensemble member $instance"

    # clean up from last advance
    # some of these must be copied at some point?? for diagnostics?
    \rm -f  restart.nc  restart.hydro.nc  dart_restart  restart.assimOnly.nc
    \rm -f wrfHydro_advance_information.txt
    \rm -f  HYDRO_RST.*  RESTART.* 
    # if perturbed forcings are used, there will be *LDASIN_DOMAIN* here. see below.
    \rm -f  *.LDASOUT_DOMAIN*  *LDASIN_DOMAIN*
    \rm -f  *.LSMOUT_DOMAIN*  *.RTOUT_DOMAIN*  *.CHRTOUT*  *.CHANOBS*  frxst_pts_out.txt
    \rm -f  qstrmvol*  diag_hydro.*  stderr.txt stdout.txt  GW_*.txt  *.GW_DOMAIN*
 
    # need the wrfHydro restart files for the output of dart_to_wrfHydro
    \ln -sv ../restart.$instance.nc  restart.nc   || exit 2
    \ln -sv ../restart.hydro.$instance.nc  restart.hydro.nc   || exit 2
    if ( $assimOnly_active ) \
	\ln -sv ../restart.assimOnly.$instance.nc  restart.assimOnly.nc   || exit 2

    # the input file is the name of the dart_restart.instance? 
    # There must be reasons for being cryptic.
    \ln -sv ../$input_file           dart_restart || exit 2
    
    # push the assimilation back to the model
    # This modifies
    # restart.nc -> ../restart.$instance.nc
    # restart.hydro.nc -> ../restart.hydro$instance.nc
    ../dart_to_wrfHydro                          || exit 2

    if ( ! -e wrfHydro_advance_information.txt ) then
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
	echo "ERROR: dart_to_noah failed for member $ensemble_member"
        exit 1
    endif

    # This next two parts are based on using one-hour forcing files
    # since the minimum time to advance the model seems to be 1 hour.
    # (kday, khour, but no kminute, for example)
    # dart_to_wrfHydro provides the setting for namelist.hrldas:khour
    # we need to put that value in the local copy of namelist.hrldas

    set numadvancestr = `\grep -i khour wrfHydro_advance_information.txt`
    set numadvancestr = `echo $numadvancestr | sed -e "s#[=,']# #g"`
    set numadvancestr = `echo $numadvancestr | sed -e 's#"# #g'`
    set numadvances   = `echo $numadvancestr[$#numadvancestr]`

# seems most efficient to only write restarts after the desired advance. 
# but this coul/would? wreak havoc on our currently inept mechanism for getting restarts
# (which is to advance a single hour past the last integration period - which
#  may not correspond to the restart frequency).
# plus there is a check for a single restart file below. 
# also not tested for RESTART_FREQUENCY_HOURS

    ## ALSO have to keep the start time in hrldas abreast of the advancing.
    ## else the forcing data seems to have no effect.
    set restartFileTime = `ncdump -v Times restart.nc | tail -2 | head -1 | cut -d'"' -f2`
    set restartFileYyyy = `echo $restartFileTime | cut -d- -f1`
    set restartFileMm = `echo $restartFileTime | cut -d- -f2`
    set restartFileDd = `echo $restartFileTime | cut -d- -f3 | cut -d_ -f1`
    set restartFileHh = `echo $restartFileTime | cut -d_ -f2 | cut -d: -f1`

ex namelist.hrldas <<ex_end
g;KHOUR ;s;= .*;= $numadvances;
g;START_YEAR;s;= .*;= $restartFileYyyy;
g;START_MONTH;s;= .*;= $restartFileMm;
g;START_DAY;s;= .*;= $restartFileDd;
g;START_HOUR;s;= .*;= $restartFileHh;
wq
ex_end

    echo '******************************************************************************'
    grep START_ namelist.hrldas | grep -v !
    echo '******************************************************************************'

    # The forcing has to be for the NEXT "FORCING_TIMESTEP", apparently.
    # FORCING_TIMESTEP is defined in namelist.input At this point, dart_to_wrfHydro
    # has assumptions that the forcing_timestep is one hour.

    # grep -n identifies the (line number): in the outupt, this becomes skipNlines
    set numfilestring = `\grep -ni nfiles wrfHydro_advance_information.txt`
    set numfilestring = `echo $numfilestring | sed -e "s#[=,':]# #g"`
    set numfilestring = `echo $numfilestring | sed -e 's#"# #g'`
    set numfiles      = `echo $numfilestring[$#numfilestring]`
    set skipNlines    = `echo $numfilestring[1]`

    #-------------------------------------------------------------------
    # Block 2.5: PERTURBED FORCING 
    #
    # It may be desirable to keep the perturbed forcings (for diagnostics or smoother?).
    # Here you can exercise that choice. If kept, perturbed forcings will be stored to
    # ../FORCING.perturbed/yyyymmddhh.iEns.LDASIN_DOMAIN* where the final character 
    # is the same as in ../FORCING/yyyymmddhh.LDASIN_DOMAIN*. 
    # If cycling or re-running is performed for previous timesteps after each analysis, 
    # this may become more complex, with the forcing at a given time actually depending 
    # on the time of the analysis. (That is the forcing is changing as the analysis moves
    # into the future). But Im' not going to worry about that yet. 
    #
    # Create the forcings in the current directory and only move if they are to be kept.  
    # The FORGING.perturbed/* filenames should have date and the instance/nnnn should be specified.  
    # then symlinked to remove this.
    ## Alter the namelist.hrldas to point INDIR='./'

    #This hasnt been tested yet. 
    if ($assimOnly_active1 && $assimOnly_active2) set assimOnly_active = 1

    if ( $assimOnly_active ) then 

	## The name of the forcing file is conveniently supplied in 
	## wrfHydro_advance_information.txt. Though ( fix ) the times in the 2 top lines
	## seem wrong.
	
	@ ifile = 1
	while ($ifile <= $numfiles)
	    @ linenum = $skipNlines + $ifile
	    set FNAME = `\head -$linenum wrfHydro_advance_information.txt | tail -1`
	    set FDATE = `echo $FNAME | sed -e "s#[.,']# #g"`
	    set FDATE = `echo $FDATE[1]`

	    set FFILE = `\ls ../FORCING/$FDATE.LDASIN_DOMAIN*`
	    set FFILElocal = `echo $FFILE | \cut -d'/' -f3`

	    echo $FFILE
	    echo $FFILElocal

	    ## get the precip multiplier out of the restart.assimOnly.nc file.
	    ## this will only work for scalar precip! 
	    set precipMult = `ncks -H -v precipMult ../restart.assimOnly.${instance}.nc \
                              | head -1 | cut -d= -f2`

	    ## multiply and put the restulting file in the right place
            ncap2 -s "RAINRATE=RAINRATE*${precipMult}" \
                     $FFILE ../FORCING.perturbed/ensemble.${instance}/$FFILElocal

	    ## change the location of the input in the namelist.hrldas
	    ## maybe this could be done elsewhere, outside loop, but it's lightweight.
ex namelist.hrldas <<ex_end
g;INDIR;s;= .*;= "../FORCING.perturbed/ensemble.${instance}/";
wq
ex_end

	    ## need to update the timestamp in the restart.
	    ## the timestamp should match that of the forcing file??
	    ## jlm - fixme

	    @ ifile = $ifile + 1
	end

	

    endif 
    
    #-------------------------------------------------------------------
    # Block 3: advance the model
    #          In this case, we are saving the run-time messages to
    #          a LOCAL file, which makes debugging easier.
    #          integrate_model is hardcoded to expect input in temp_ic 
    #          and it creates temp_ud as output. 
    #          Your model will likely be different.
    #-------------------------------------------------------------------
    echo "advance the model"
    ## i just want the model to be quiet so I can focus on the DART output
    mpirun -np 2 ../wrf_hydro.exe >& /tmp/jamesmccWfrHydroEnsOutputJunk.$process
    #../wrf_hydro.exe >& /tmp/jamesmccWfrHydroEnsOutputJunk.$process
    \rm -f /tmp/jamesmccWfrHydroEnsOutputJunk.$process

    @ lsm_status = `\ls -1 RESTART*DOMAIN* | wc -l`
    @ hydro_status = `\ls -1 HYDRO_RST* | wc -l`
    @ numadvancesNum = $numadvances

    if ( $lsm_status < 1 || $hydro_status < 1 )  then
	echo "ERROR: wrfHydro died"
	echo "ERROR: wrfHydro died"
	\ls -l
	exit 23
    endif 

    if ( $lsm_status > $numadvancesNum || $hydro_status > $numadvances )  then
	if ( $lsm_status > $numadvancesNum ) then 
	    \ls -l RESTART*DOMAIN*
	    echo "WARNING: wrfHydro created the above RESTART files. only expected # $numadvances" 
	endif 
	if ( $hydro_status > $numadvancesP1 ) then 
	    \ls -l HYDRO_RST*
	    echo "WARNING: wrfHydro created the above HYDRO_RST files. Only expected # $numadvances" 
	endif 
    endif

    #-------------------------------------------------------------------
    # Block 4: wrfHydro_to_dart and managing output files. 
    #          rename files to reflect the ensemble member ID
    #-------------------------------------------------------------------

    # Do this before setting up the next run as there is an unwatned hydro restart file.

    # Determine model integration period of interest (which may contain multiple indiv
    # model integrations) and create a directory in central dir appropriately stamped 
    # to catch the output during this period.. 
    # The timestamps of the first and last LDASOUT files give us this even though we 
    # dont want the last LDASOUT file (it's for one hour beyond the desired integration period
    # thanks to HRLDAS).
    set integStart = `\ls -1 *LDASOUT_DOMAIN* | \head -1 | \cut -d. -f1`
    set integEnd   = `\ls -1 RESTART*_DOMAIN* | \tail -1 | \cut -d. -f2 | cut -d_ -f1`
    set integDir =  OUTPUT/model_integration.${integStart}-${integEnd}.$instance
    set integDirCurrent =  ../${integDir}
    \mkdir $integDirCurrent

    # Fixed HRLDAS to do restarts at the end of the loop, after time advance,  
    # with/after wrfHydro. So I dont have to clean up a bunch of files.

    # Move the output files (*not* restarts)
    foreach outFile ( GW_*.txt frxst_pts_out.txt qstrmvolrt_accum.txt )
	\mv $outFile ${integDirCurrent}/.
    end 
    # these have their own timestamps but tag them with ensId/instance. 
    foreach outFile ( *.LDASOUT_DOMAIN* *.LSMOUT_DOMAIN* *.RTOUT_DOMAIN* *.CHRTOUT* *.CHANOBS* )
	\mv $outFile ${integDirCurrent}/${outFile}.${instance}.nc
    end 

    ## These are the parameters which were used in the model advance. 
    if ( $assimOnly_active ) cp restart.assimOnly.nc $integDirCurrent/.

    # Set the new/latest restart for ingest to dart
    set RESTARTlsm = `\ls -1  RESTART* | \tail -1`
    set RESTARThydro = `\ls -1  HYDRO_RST* | \tail -1`
    # must force overwrite these existing links
    \ln -sf $RESTARTlsm    restart.nc        || exit 4
    \ln -sf $RESTARThydro  restart.hydro.nc  || exit 4

    ../wrfHydro_to_dart              || exit 4

    \mv -v  dart_ics  ../$output_file          || exit 5
    # this breaks the restart.nc and restart.hydro.nc symlinks
    # but they are reset in the next loop
    \mv -v  ${RESTARTlsm}    ${integDirCurrent}/${RESTARTlsm}.$instance.nc   || exit 5
    \mv -v  ${RESTARThydro}  ${integDirCurrent}/${RESTARThydro}.$instance.nc || exit 5
    \ln -sfv ${integDir}/${RESTARTlsm}.$instance.nc   ../restart.$instance.nc
    \ln -sfv ${integDir}/${RESTARThydro}.$instance.nc ../restart.hydro.$instance.nc
    # the linking (vs. cp ing) in these last two lines implies that the model-created 
    # restart files are not sacred, they will be overwritten by dart_to_wrfHydro

    ## increment
    @ state_copy++
    @ ensemble_member_line = $ensemble_member_line + 3
    @ input_file_line = $input_file_line + 3
    @ output_file_line = $output_file_line + 3

end  ## loop over ensemble members for this process.

# Change back to original directory and get rid of temporary directory.
# If all goes well, there should be no need to keep this directory.
# If you are debugging, you may want to keep this directory. 

cd ..
\rm -rf $temp_dir

# MANDATORY - Remove the control_file to signal completion. If it still
# exists in CENTRALDIR after all the ensemble members have been advanced,
# it means one or more of the advances failed and is an ERROR CONDITION.

\rm -rf $control_file

exit 

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

