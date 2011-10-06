#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# Top level script to run a COSMO assimilation experiment.
#

#=============================================================================
#  initial setup
#=============================================================================

set expname = bob
set ncycles = 10

#=============================================================================
#  main cycle loop
#=============================================================================

set icyc = 1
while ( $icyc < $ncycles )

   #==========================================================================
   #  1. set up and run filter to assimilate next set of observations
   #==========================================================================
   
   sed filter.job.template > filter.job
   bsub < filter.job
   
   #==========================================================================
   #  2. convert filter output to model input files
   #==========================================================================
   
   bsub < convert1.job
   
   #==========================================================================
   #  3. run N copies of the model here
   #==========================================================================
   
   sed model.job.template > model.job
   bsub < model.job
   
   #==========================================================================
   #  4. convert model output files to filter input files
   #==========================================================================
   
   bsub < convert2.job
   
   #==========================================================================
   #  5. (optional) archive files here
   #==========================================================================
   
   echo hsi or mss commands here

   @ icyc ++
end

#=============================================================================
#  bottom of main cycle loop
#=============================================================================

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$


