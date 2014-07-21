#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: run_pmo.csh 6256 2013-06-12 16:19:10Z thoar $
#
# This is an example script for how to run a perfect_model experiment.
# All the required files have been staged in CENTRALDIR and the namelists
# have been customized for the experiment. Since perfect_model_obs is a single-
# threaded program ... and we are only running one instance of NOAH ... we are
# running this on the command line.

#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Check to make sure all the required files have been staged in CENTRALDIR
#==============================================================================

set CENTRALDIR = `pwd`

@ BAIL = 0
foreach FILE ( wrfinput namelist.hrldas Noah_hrldas_beta SOILPARM.TBL \
               VEGPARM.TBL GENPARM.TBL URBPARM.TBL obs_seq.in input.nml \
               perfect_model_obs dart_to_noah noah_to_dart run_pmo.csh \
               advance_model.csh restart.nc perfect_ics)

   if ( ! -e $FILE ) then
      echo "$FILE is needed but not present in CENTRALDIR"
      @ BAIL = 1
   endif

end

if ( $BAIL > 0 ) then
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   echo "FATAL ERROR ... stage the missing file(s) and try again."
   exit 1
endif

./perfect_model_obs

exit 0

# <next few lines under version control, do not edit>
# $URL: https://proxy.subversion.ucar.edu/DAReS/DART/releases/Lanai/models/noah/shell_scripts/run_pmo.csh $
# $Revision: 6256 $
# $Date: 2013-06-12 10:19:10 -0600 (Wed, 12 Jun 2013) $

