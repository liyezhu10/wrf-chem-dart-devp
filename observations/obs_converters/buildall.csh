#!/bin/csh

set SNAME = $0
set clobber

set startdir=`pwd`
cd $startdir

foreach project ( AIRNOW IASI_CO IASI_O3 MODIS MOPITT_CO NCEP/prep_bufr NCEP/ascii_to_obs PANDA )

   echo
   echo "==================================================="
   echo "Building $project support."
   echo "==================================================="
   echo

   set dir = $project
   set FAILURE = 0

   switch ("$dir")

      case NCEP/prep_bufr
         cd $dir
         echo "building in `pwd`"
         ./install.sh || set FAILURE = 1
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
