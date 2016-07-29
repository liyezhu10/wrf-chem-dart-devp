#!/bin/csh
#
# DART software - Copyright 2004 - 2016 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script will concatenate a series of filter_restart.xxxx.nc files
# into a single file. That file will then be coerced to be compatible
# with the True_State.nc that has xxxx timesteps.
# At present, the DART matlab scripts historically use dimension and variable names
# consistent with the 'dart diagnostic' files. The filter_restart files
# have different dimension and variable names.

foreach FILE ( filter_restart*nc )

  echo "Adding the record dimension 'copy' to $FILE"

  # Temporarily create a 'copy' dimension as the record dimension
  # Turn the record dimension into a fixed dimension (ultimately
  # 'time' will be the record dimension.

  ncecat -O -u copy $FILE temp.nc
  ncks -O --fix_rec_dmn=copy temp.nc $FILE

  ncecat -O -u time $FILE temp.nc
  mv temp.nc $FILE

end

ncrcat filter_restart*nc temp.nc

# we have to add the required variables and global attributes

ncks -A -v time,level,t_ps_lon,t_ps_lat,u_v_lon,u_v_lat,copy,CopyMetaData True_State.nc temp.nc

# lastly, change the CopyMetaData to something that works with the DART diagnostics.
# We are going to lie to it and say this is the ensemble mean
#                                  0123456789012345
ncap2 -O -s 'CopyMetaData(0,0:15)="ensemble mean   "' temp.nc Posterior_Diag.nc

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

