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
0 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MOPITT_CO.ksh
1 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_IASI_CO.ksh
2 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_IASI_O3.ksh
3 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_OMI_NO2.ksh
4 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_CO.ksh
5 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_O3.ksh
6 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_NO2.ksh
7 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_SO2.ksh
8 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_PM10.ksh
9 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_AIRNOW_PM25.ksh
10 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_CO.ksh
11 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_O3.ksh
12 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PANDA_PM25.ksh
13 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MODIS_AOD.ksh
14 ${CLUSTER_DIR}/WRFCHEM-DART_RUN_MET_OBS.ksh
EOF
#
   srun --multi-prog ./jobs.file
#
