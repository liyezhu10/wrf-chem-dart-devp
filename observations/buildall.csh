#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

set SNAME = $0
set clobber

set startdir=`pwd`

# The NCEP prep_bufr converters are a special case

set FAILURE = 0
cd NCEP/prep_bufr
./install.sh || set FAILURE = 1

if ( $FAILURE != 0 ) then
   echo
   echo "ERROR unsuccessful build in $dir"
   echo "ERROR unsuccessful build in $dir"
   echo "ERROR unsuccessful build in $dir"
   echo
   exit -1
endif

cd $startdir

foreach project ( AIRNOW IASI IASI_CO IASI_O3 MODIS MOPITT_CO NCEP PANDA )

   set dir = $project:h
   set FAILURE = 0

   switch ("$dir")

      case *NCEP*
         cd $dir/ascii_to_obs/work
         echo "building in `pwd`"
         ./quickbuild.csh || set FAILURE = 1
      breaksw
         
      default:
         cd $dir/work
         echo "building in `pwd`"
         ./quickbuild.csh || set FAILURE = 1
      breaksw
         
   endsw

   if ( $FAILURE != 0 ) then
      echo
      echo "ERROR unsuccessful build in $dir"
      echo "ERROR unsuccessful build in $dir"
      echo "ERROR unsuccessful build in $dir"
      echo
      exit -1
   endif

   cd $startdir
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


