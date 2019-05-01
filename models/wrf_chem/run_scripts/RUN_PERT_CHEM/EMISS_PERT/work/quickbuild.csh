#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: quickbuild.csh 13126 2019-04-25 01:59:32Z thoar@ucar.edu $
#
# This script compiles all executables in this directory.

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the 
# resulting source file is used by all the remaining programs, 
# so this MUST be run first.
#----------------------------------------------------------------------

\rm -f preprocess *.o *.mod

set MODEL = "WRF-Chem/DART run scripts"

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
      breaksw
   endsw
end

\rm -f *.o *.mod input.nml*_default

echo "Success: All single task ${MODEL} programs compiled."  

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

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/branches/mizzi/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT/work/quickbuild.csh $
# $Revision: 13126 $
# $Date: 2019-04-24 19:59:32 -0600 (Wed, 24 Apr 2019) $

