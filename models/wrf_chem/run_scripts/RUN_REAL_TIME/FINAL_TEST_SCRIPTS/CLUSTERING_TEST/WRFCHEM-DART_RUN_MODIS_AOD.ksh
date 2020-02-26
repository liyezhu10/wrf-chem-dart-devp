#!/bin/ksh -aeux
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
#########################################################################
#
# RUN MODIS AOD OBSERVATIONS
#
#########################################################################
#
   if ${RUN_MODIS_AOD_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/modis_aod_obs ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/modis_aod_obs
         cd ${RUN_DIR}/${DATE}/modis_aod_obs
      else
         cd ${RUN_DIR}/${DATE}/modis_aod_obs
      fi
#
# set file prefix for MODIS
# this depends on versions and file times (edit if necessary)
      export FILE_PRE='MYD04_L2.A'
#
# set file suffix for MODIS
# this depends on versions and file times (edit if necessary)
      export FILE_EXT='hdf'
#
      export MODIS_INDIR=${EXPERIMENT_MODIS_AOD_DIR}
      export OUTFILE=modis_aod_ascii_${YYYY}${MM}${DD}${HH}
      (( N_YYYY=${YYYY}+0 ))
      (( N_MM=${MM}+0 ))
      (( N_DD=${DD}+0 ))
      (( N_HH=${HH}+0 ))
      (( N_ASIM_WIN=${ASIM_WINDOW}+0 ))
#
      export FILE=modis_extract_hdf.pro
      rm -rf ${FILE}
      cp ${DART_DIR}/observations/MODIS/native_to_ascii/${FILE} ./.
      idl << EOF
.compile modis_extract_hdf.pro
modis_extract_hdf, "${MODIS_INDIR}", "${OUTFILE}", ${N_YYYY}, ${N_MM}, ${N_DD}, ${N_HH}, ${N_ASIM_WIN}, ${NL_MIN_LON}, ${NL_MAX_LON}, ${NL_MIN_LAT}, ${NL_MAX_LAT}
EOF
#
# convert ASCII to obs_seq file
#      rm -rf input.nml
#      rm -rf modis_asciidata.input
#      rm -rf modis_obs_seq.out
      if [[ -s modis_aod_ascii_${YYYY}${MM}${DD}${HH} ]]; then
         cp modis_aod_ascii_${YYYY}${MM}${DD}${HH} modis_asciidata.input
         cp ${DART_DIR}/observations/MODIS/work/input.nml ./.
         ${DART_DIR}/observations/MODIS/work/modis_ascii_to_obs
         if [[ -s modis_obs_seq.out ]]; then
            export MODIS_FILE=obs_seq_modis_aod_${DATE}.out
            mv modis_obs_seq.out ${MODIS_FILE}
         else
            touch NO_MODIS_AOD_${DATE}
         fi
      else
         touch NO_MODIS_AOD_${DATE}
      fi
   fi
