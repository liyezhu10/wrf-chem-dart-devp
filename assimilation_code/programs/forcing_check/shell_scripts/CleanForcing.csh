#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# qsub -I -l select=1:ncpus=1:mpiprocs=1 -l walltime=04:00:00 -q economy -A P86850054
#
# qcmd -q share -l select=1 -l walltime=06:00:00 -- ./CleanForcing.csh
 
#===========================================================================
# Process the DS199.1 files to remove negative values 
#===========================================================================

cat << EOF >! input.nml
# variables in the DS199.1 datasets and their 'long_name' attribute
#    a2x6h_Sa_tbot      "Temperature at the lowest model level" ;
#    a2x6h_Sa_ptem      "Potential temperature at the lowest model level" ;
#    a2x6h_Sa_shum      "Specific humidity at the lowest model level" ;
#    a2x6h_Sa_dens      "Air density at the lowest model level" ;
#    a2x6h_Sa_pbot      "Pressure at the lowest model level" ;
#    a2x6h_Sa_pslv      "Sea level pressure" ;
#    a2x6h_Faxa_lwdn    "Downward longwave heat flux" ;
#    a2x6h_Faxa_rainc   "Convective            precipitation rate" ;
#    a2x6h_Faxa_rainl   "Large-scale (stable)  precipitation rate" ;
#    a2x6h_Faxa_snowc   "Convective            snow rate (water equivalent)" ;
#    a2x6h_Faxa_snowl   "Large-scale (stable)  snow rate (water equivalent)" ;
#    a2x6h_Faxa_swndr   "Direct near-infrared  incident solar radiation" ;
#    a2x6h_Faxa_swvdr   "Direct visible        incident solar radiation" ;
#    a2x6h_Faxa_swndf   "Diffuse near-infrared incident solar radiation" ;
#    a2x6h_Faxa_swvdf   "Diffuse visible       incident solar radiation" ;

&forcing_check_nml
   debug                   = .false.
   fail_on_missing_field   = .false.
   do_all_numeric_fields   = .false.
   fieldnames              = 'a2x6h_Sa_ptem',
                             'a2x6h_Sa_shum',
                             'a2x6h_Sa_dens',
                             'a2x6h_Sa_pbot',
                             'a2x6h_Sa_pslv',
                             'a2x6h_Faxa_lwdn',
                             'a2x6h_Faxa_rainc',
                             'a2x6h_Faxa_rainl',
                             'a2x6h_Faxa_snowc',
                             'a2x6h_Faxa_snowl',
                             'a2x6h_Faxa_swndr',
                             'a2x6h_Faxa_swvdr',
                             'a2x6h_Faxa_swndf',
                             'a2x6h_Faxa_swvdf' 
   fieldlist_file          = ''
   only_report_differences = .false.
   /

&utilities_nml
   /
EOF


#===========================================================================

# set DSSBASE = /glade/p/rda/data/ds199.1/CAM_DATM.cpl
set OURBASE = /glade/p/image/thoar/CAM_DATM/4xdaily/CAM_DATM.cpl

@ YEAR = 1997
while ($YEAR <= 2010)

   @ MEMBER = 1
   while ($MEMBER <= 80)

      set FILENAME = `printf %s_%04d.ha2x1dx6h.%04d.nc $OURBASE $MEMBER $YEAR`
      echo $FILENAME

      if (-e $FILENAME ) then
         echo $FILENAME | ../work/forcing_check
      else
         echo "$ERROR $FILENAME does not exist."
      endif

      @ MEMBER ++
   end

   @ YEAR ++
end

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
