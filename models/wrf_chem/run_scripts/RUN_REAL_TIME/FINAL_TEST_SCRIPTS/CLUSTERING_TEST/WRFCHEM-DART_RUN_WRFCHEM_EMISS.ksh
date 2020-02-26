#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
# GENERATE BIO, FIRE, and ANTHRO EMISSIONS
   rm -rf jobs.file
   cat << EOF > jobs.file
0 ${CLUSTER_DIR}/WRFCHEM-DART_BIO_EMISS.ksh
1 ${CLUSTER_DIR}/WRFCHEM-DART_FIRE_EMISS.ksh
2 ${CLUSTER_DIR}/WRFCHEM-DART_ANTHRO_EMISS.ksh
EOF
#
   srun --multi-prog ./jobs.file
#
