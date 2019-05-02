#!/bin/csh
#
# DART software - Copyright © 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
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
# $URL$
# $Revision$
# $Date$

