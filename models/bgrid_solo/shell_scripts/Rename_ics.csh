#!/bin/csh
#
# DART software - Copyright 2004 - 2016 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script is/was used to rename the dimensions and variables in the
# existing netCDF files to be consistent between the DART diagnostic files
# and the model restart files. This should not be needed going forward,
# as the input files and the model_mod.f90 are now consistent.
#
# This will allow us to use use Restart_to_Diagnostic.csh
# to recreate a 'Posterior_Diag.nc' from the existing filter_restart.nnnn.nc files
# so we can use the DART diagnostic matlab scripts. 

foreach FILE ( True_State.nc )

   echo "Working on $FILE"

   # There are a host of problems with renaming things in netCDF4 
   # Seems the safest is to do it one variable at a time,
   # and then change all the dimensions one at a time.
   # This is a terrible waste. 

   ncrename -v TmpI,t_ps_lon $FILE
   ncrename -v TmpJ,t_ps_lat $FILE
   ncrename -v VelI,u_v_lon  $FILE
   ncrename -v VelJ,u_v_lat  $FILE
   ncrename -v lev,level     $FILE

   ncrename -d TmpI,t_ps_lon $FILE
   ncrename -d TmpJ,t_ps_lat $FILE
   ncrename -d VelI,u_v_lon  $FILE
   ncrename -d VelJ,u_v_lat  $FILE
   ncrename -d lev,level     $FILE

   # we have to add the required variables and global attributes
   # and delete the ugly history thrashings

   ncatted -a long_name,t_ps_lat,m,c,'latitude of temperature and surface pressure' $FILE
   ncatted -a long_name,t_ps_lon,m,c,'longitude of temperature and surface pressure' $FILE
   ncatted -a long_name,u_v_lat,m,c,'latitude of u and v winds' $FILE
   ncatted -a long_name,u_v_lon,m,c,'longitude of u and v winds' $FILE
   ncatted -a valid_range,u_v_lon,o,d,'0.0,360.0' $FILE
   ncatted -a history,global,d,, $FILE

end

foreach FILE ( perfect_ics*.nc )

   echo "Working on $FILE"

   # There are a host of problems with renaming things in netCDF4 
   # Seems the safest is to do it one variable at a time,
   # and then change all the dimensions one at a time.
   # This is a terrible waste. 

   ncrename -v PS,ps        $FILE
   ncrename -v T,t          $FILE
   ncrename -v U,u          $FILE
   ncrename -v V,v          $FILE

   ncrename -d lon,t_ps_lon $FILE
   ncrename -d lat,t_ps_lat $FILE
   ncrename -d slon,u_v_lon $FILE
   ncrename -d slat,u_v_lat $FILE
   ncrename -d lev,level    $FILE

   ncks -A -v level    True_State.nc $FILE
   ncks -A -v t_ps_lon True_State.nc $FILE
   ncks -A -v t_ps_lat True_State.nc $FILE
   ncks -A -v u_v_lon  True_State.nc $FILE
   ncks -A -v u_v_lat  True_State.nc $FILE

   # we have to add the required variables and global attributes
   # and delete the ugly history thrashings

   ncatted -a long_name,t_ps_lat,m,c,'latitude of temperature and surface pressure' $FILE
   ncatted -a long_name,t_ps_lon,m,c,'longitude of temperature and surface pressure' $FILE
   ncatted -a long_name,u_v_lat,m,c,'latitude of u and v winds' $FILE
   ncatted -a long_name,u_v_lon,m,c,'longitude of u and v winds' $FILE
   ncatted -a history_of_appended_files,global,d,, $FILE
   ncatted -a history,global,d,, $FILE

end

foreach FILE ( filter_ics.*.nc )

   echo "Working on $FILE"

   ncwa -a copy,0 $FILE bob.nc

   mv bob.nc $FILE

   # There are a host of problems with renaming things in netCDF4 
   # Seems the safest is to do it one variable at a time,
   # and then change all the dimensions one at a time.
   # This is a terrible waste. 

   ncrename -v PS,ps        $FILE
   ncrename -v T,t          $FILE
   ncrename -v U,u          $FILE
   ncrename -v V,v          $FILE
   ncrename -v lat,t_ps_lat $FILE
   ncrename -v lon,t_ps_lon $FILE
   ncrename -v slat,u_v_lat $FILE
   ncrename -v slon,u_v_lon $FILE

   ncrename -d lon,t_ps_lon $FILE
   ncrename -d lat,t_ps_lat $FILE
   ncrename -d slon,u_v_lon $FILE
   ncrename -d slat,u_v_lat $FILE
   ncrename -d lev,level    $FILE

   ncks -A -v level True_State.nc $FILE

   # we have to add the required variables and global attributes
   # and delete the ugly history thrashings

   ncatted -a long_name,t_ps_lat,m,c,'latitude of temperature and surface pressure' $FILE
   ncatted -a long_name,t_ps_lon,m,c,'longitude of temperature and surface pressure' $FILE
   ncatted -a long_name,u_v_lat,m,c,'latitude of u and v winds' $FILE
   ncatted -a long_name,u_v_lon,m,c,'longitude of u and v winds' $FILE
   ncatted -a valid_range,u_v_lon,o,d,'0.0,360.0' $FILE
   ncatted -a history_of_appended_files,global,d,, $FILE
   ncatted -a history,global,d,, $FILE

end
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

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

