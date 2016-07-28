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

foreach FILE ( filter_ics*nc )

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

# There are a host of problems with renaming things in netCDF4 
# Seems the safest is to do it one variable at a time,
# and then change all the dimensions one at a time.
# This is a terrible waste. 

ncrename -v PS,ps     temp.nc
ncrename -v T,t       temp.nc
ncrename -v U,u       temp.nc
ncrename -v V,v       temp.nc

ncrename -d lon,TmpI  temp.nc
ncrename -d lat,TmpJ  temp.nc
ncrename -d slon,VelI temp.nc
ncrename -d slat,VelJ temp.nc

# we have to add the required variables and global attributes

ncks -A -v time,lev,TmpI,TmpJ,VelI,VelJ,copy,CopyMetaData True_State.nc temp.nc
ncatted -a model,global,a,c,'FMS_Bgrid'             temp.nc

# lastly, change the CopyMetaData to something that works with the DART diagnostics.
# We are going to lie to it and say this is the ensemble mean
#                                  0123456789012345
ncap2 -O -s 'CopyMetaData(0,0:15)="ensemble mean   "' temp.nc Posterior_Diag.nc




# netcdf True_State {
# dimensions:
#    metadatalength = 64 ;
#    locationrank = 3 ;
#    copy = 1 ;
#    time = UNLIMITED ; // (120 currently)
#    NMLlinelen = 129 ;
#    NMLnlines = 347 ;
#    StateVariable = 28200 ;
#    TmpI = 60 ;
#    TmpJ = 30 ;
#    lev = 5 ;
#    VelI = 60 ;
#    VelJ = 29 ;
# int copy(copy) ;
# char CopyMetaData(copy, metadatalength) ;
# float ps(time, copy, TmpJ, TmpI) ;
# float  t(time, copy, lev, TmpJ, TmpI) ;
# float  u(time, copy, lev, VelJ, VelI) ;
# float  v(time, copy, lev, VelJ, VelI) ;

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

