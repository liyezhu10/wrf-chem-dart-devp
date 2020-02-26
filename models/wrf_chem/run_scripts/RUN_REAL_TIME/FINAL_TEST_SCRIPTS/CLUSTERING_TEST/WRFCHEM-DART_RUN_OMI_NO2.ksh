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
# RUN OMI NO2 OBSERVATIONS
#
#########################################################################
#
   if ${RUN_OMI_NO2_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/omi_no2_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/omi_no2_obs
         cd ${RUN_DIR}/${DATE}/omi_no2_obs
      else
         cd ${RUN_DIR}/${DATE}/omi_no2_obs
      fi
#
      export FILE=obs_seq_comb_filtered_${DATE}.out
      cp ${EXPERIMENT_OMI_NO2_DIR}/${DATE}/${FILE} ./obs_seq_${DATE}.out
#
      if [[ -s obs_seq_${DATE}.out ]]; then
         cp obs_seq_${DATE}.out obs_seq_omi_no2_${DATE}.out
         rm obs_seq_${DATE}.out
      else
          touch NO_DATA_${DATE}
      fi
   fi
