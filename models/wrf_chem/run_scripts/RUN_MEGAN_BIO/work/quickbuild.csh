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
./make_util megan_bio_emiss
./make_util megan_xform
./make_util surfdata_xform

cd work
rm -rf *.exe *.nc

mv ../megan_bio_emiss ./megan_bio_emiss.exe
mv ../megan_xform ./megan_xform.exe
mv ../surfdata_xform ./surfdata_xform.exe

if (! -e megan_bio_emiss.exe) then
   echo megan_bio_emiss.exe did not buid properly
   exit 1
endif
if (! -e megan_xform.exe) then
   echo megan_xform.exe did not buid properly
   exit 1
endif
if (! -e surfdata_xform.exe) then
   echo surfdata_xform.exe did not buid properly
   exit 1
endif

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/run_scripts/RUN_MEGAN_BIO/work/quickbuild.csh $
# $Revision: 7387 $
# $Date: 2015-01-16 15:44:50 -0700 (Fri, 16 Jan 2015) $

