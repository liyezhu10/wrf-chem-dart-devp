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
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------

\rm -f *.o *.mod 

set MODEL = "WRF-Chem/DART run_scripts"

@ n = 0

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   switch ( $TARGET )
   case mkmf_preprocess:
      @ n = $n + 1
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "${MODEL} build number ${n} is ${PROG}" 
      \rm -f ${PROG}
      csh $TARGET || exit $n
      make        || exit $n
\      breaksw
   endsw
end

\rm -f *.o *.mod  input.nml*_default

echo "Success: All ${MODEL} programs compiled."  

#----------------------------------------------------------------------
# Build the MPI-enabled target(s) 
#----------------------------------------------------------------------

foreach TARGET ( mpi_mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mpi_mkmf_##'`

   switch ( $TARGET )
   case mpi_mkmf_preprocess:
      @ n = $n + 1
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "${MODEL} build number ${n} is ${PROG}" 
      \rm -f ${PROG}
      csh $TARGET -mpi
      make        || exit $n
      breaksw
   endsw
end

\rm -f *.o *.mod input.nml*_default

echo "Success: All MPI-enabled ${MODEL} programs compiled."  

#----------------------------------------------------------------------
# Build all the single-threaded mozbc
#----------------------------------------------------------------------

cd ../mozbc-dart
./make_mozbc

cd ../work
cp ../mozbc-dart/mozbc ./mozbc.exe

if (! -e mozbc.exe) then
   echo mozbc did not buid properly
   exit 1
endif

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT/work/quickbuild.csh $
# $Revision: 7387 $
# $Date: 2015-01-16 15:44:50 -0700 (Fri, 16 Jan 2015) $

