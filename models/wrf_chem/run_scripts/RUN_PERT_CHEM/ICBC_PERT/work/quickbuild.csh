#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# This script compiles all executables in this directory.

\rm -f *.o *.mod Makefile

set MODEL = "WRF-Chem/DART run_scripts"

@ n = 0

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   switch ( $TARGET )
   case mkmf_perturb_chem_icbc_CORR_RT_MA_MPI:
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

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

#----------------------------------------------------------------------
# exit if not compiling with MPI support
#----------------------------------------------------------------------

if ( $#argv == 1 && "$1" == "-mpi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
else if ( $#argv == 1 && "$1" == "-nompi" ) then
  echo "Success: All single task DART programs compiled."  
  echo "Script is exiting without building the MPI version of the DART programs."
  exit 0
else
  echo ""
  echo "Success: All single task DART programs compiled."  
  echo "Script now compiling MPI parallel versions of the DART programs."
  echo "Run the quickbuild.csh script with a -nompi argument or"
  echo "edit the quickbuild.csh script and add an exit line"
  echo "to bypass compiling with MPI to run in parallel on multiple cpus."
  echo ""
endif

#----------------------------------------------------------------------
# to disable an MPI version of perturb_chem_icbc_CORR_RT_MA
# call this script with the -nompi argument.
#----------------------------------------------------------------------

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_perturb_chem_icbc_CORR_RT_MA with MPI"

csh  mkmf_perturb_chem_icbc_CORR_RT_MA -mpi
make || exit $n

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

