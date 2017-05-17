#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

#****************************************************************************
# Tarkeshwar Singh
# PhD, IIT Delhi
# Email: tarkphysics87@gmail.com
#
# Purpose: This script prints output of sortvar.F90 routine from gcm.log files 
#          for each ensemble member. 
#
# It prints GLOBAL 'masse', 'rmsdpdt', 'energie', 'enstrophie', 'entropie', 'rmsv', 'mt.ang'
# For model stability, these parameter should be conversed ans stable with time.
# If these are increasing or decreasing significantly with time during free run or 
# assimilation then adjust values of 'iperiod', 'tetagdiv', 'tetagrot', 'tetatemp' in gcm.def
# It has been been observed that many LMDZ member become unstable during ensemble free run 
# with zoom at 35 km using  360x180x29 grid points. Use following parameter for it in gcm.def
#
# day_step = 3600 
# iperiod  = 5
# tetagdiv = 500.
# tetagrot = 700.
# tetatemp = 700.
#
# During assimilation, iperiod  = 3 is recomended.
#****************************************************************************

set num_ens = 60

set n = 0
while ($n < $num_ens)
echo "*******************************************************************************"
  @ n1 = ( $n + 1 )
   echo Ensemble Member $n1
echo "*******************************************************************************"
   echo "          masse    rmsdpdt       energie  enstrophie  entropie   rmsv    mt.ang"
   cat advance_temp$n/gcm.log | grep GLOB
  @ n ++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

