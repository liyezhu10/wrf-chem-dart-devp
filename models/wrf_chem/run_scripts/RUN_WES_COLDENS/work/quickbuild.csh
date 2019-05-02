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
# Build all the single-threaded exo_coldens and season_wes files
#----------------------------------------------------------------------

cd ../
./make_util exo_coldens

./make_util wesely

cd work
rm -rf *.exe *.nc

mv ../exo_coldens ./exo_coldens.exe

if (! -e exo_coldens.exe) then
   echo exo_coldens.exe did not buid properly
   exit 1
endif

mv ../wesely ./wesely.exe

if (! -e wesely.exe) then
   echo wesely.exe did not buid properly
   exit 1
endif

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

