#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: quickbuild.csh 7387 2015-01-16 22:44:50Z mizzi $
#
# This script compiles all executables in this directory.

#----------------------------------------------------------------------
# Build all the single-threaded megan-bio files
#----------------------------------------------------------------------

cd ../
./make_fire_emis

cd work
rm -rf *.exe *.nc

mv ../fire_emis ./fire_emis.exe

if (! -e fire_emis.exe) then
   echo fire_emis.exe did not buid properly
   exit 1
endif

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/run_scripts/RUN_MEGAN_BIO/work/quickbuild.csh $
# $Revision: 7387 $
# $Date: 2015-01-16 15:44:50 -0700 (Fri, 16 Jan 2015) $

