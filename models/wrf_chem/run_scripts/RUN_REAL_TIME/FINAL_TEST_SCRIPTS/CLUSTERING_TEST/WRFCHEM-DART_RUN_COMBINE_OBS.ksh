#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
#########################################################################
#
# RUN COMBINE OBSERVATIONS
#
#########################################################################
#
   if ${RUN_COMBINE_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/combine_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/combine_obs
         cd ${RUN_DIR}/${DATE}/combine_obs
      else
         cd ${RUN_DIR}/${DATE}/combine_obs
      fi
#
# GET EXECUTABLES
      cp ${DART_DIR}/models/wrf_chem/work/obs_sequence_tool ./.
      export NUM_FILES=0
#
# GET OBS_SEQ FILES TO COMBINE
# MET OBS
      if [[ -s ${PREPBUFR_MET_OBS_DIR}/obs_seq_prep_${DATE}.out && ${RUN_MET_OBS} ]]; then 
         (( NUM_FILES=${NUM_FILES}+1 ))
         cp ${PREPBUFR_MET_OBS_DIR}/obs_seq_prep_${DATE}.out ./obs_seq_MET_${DATE}.out
         export FILE_LIST[${NUM_FILES}]=obs_seq_MET_${DATE}.out
      fi
#
# MOPITT CO
      if [[ -s ${MOPITT_CO_OBS_DIR}/obs_seq_mopitt_co_${DATE}.out && ${RUN_MOPITT_CO_OBS} ]]; then 
         cp ${MOPITT_CO_OBS_DIR}/obs_seq_mopitt_co_${DATE}.out ./obs_seq_MOP_CO_${DATE}.out
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_MOP_CO_${DATE}.out
      fi
#
# IASI CO
      if [[ -s ${IASI_CO_OBS_DIR}/obs_seq_iasi_co_${DATE}.out && ${RUN_IASI_CO_OBS} ]]; then 
         cp ${IASI_CO_OBS_DIR}/obs_seq_iasi_co_${DATE}.out ./obs_seq_IAS_CO_${DATE}.out
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_IAS_CO_${DATE}.out
      fi
#
# IASI O3
      if [[ -s ${IASI_O3_OBS_DIR}/obs_seq_iasi_o3_${DATE}.out && ${RUN_IASI_O3_OBS} ]]; then 
         cp ${IASI_O3_OBS_DIR}/obs_seq_iasi_o3_${DATE}.out ./obs_seq_IAS_O3_${DATE}.out   
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_IAS_O3_${DATE}.out
      fi
#
# OMI NO2
      if [[ -s ${OMI_NO2_OBS_DIR}/obs_seq_omi_no2_${DATE}.out && ${RUN_OMI_NO2_OBS} ]]; then 
         cp ${OMI_NO2_OBS_DIR}/obs_seq_omi_no2_${DATE}.out ./obs_seq_OMI_NO2_${DATE}.out   
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_OMI_NO2_${DATE}.out
      fi
#
# AIRNOW O3
      if [[ -s ${AIRNOW_O3_OBS_DIR}/obs_seq_airnow_o3_${DATE}.out && ${RUN_AIRNOW_O3_OBS} ]]; then 
         cp ${AIRNOW_O3_OBS_DIR}/obs_seq_airnow_o3_${DATE}.out ./obs_seq_AIR_O3_${DATE}.out   
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_AIR_O3_${DATE}.out
      fi
#
# AIRNOW CO
      if [[ -s ${AIRNOW_CO_OBS_DIR}/obs_seq_airnow_co_${DATE}.out && ${RUN_AIRNOW_CO_OBS} ]]; then 
         cp ${AIRNOW_CO_OBS_DIR}/obs_seq_airnow_co_${DATE}.out ./obs_seq_AIR_CO_${DATE}.out   
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_AIR_CO_${DATE}.out
      fi
#
# MODIS AOD
      if [[ -s ${MODIS_AOD_OBS_DIR}/obs_seq_modis_aod_${DATE}.out && ${RUN_MODIS_AOD_OBS} ]]; then 
         cp ${MODIS_AOD_OBS_DIR}/obs_seq_modis_aod_${DATE}.out ./obs_seq_MOD_AOD_${DATE}.out   
         (( NUM_FILES=${NUM_FILES}+1 ))
         export FILE_LIST[${NUM_FILES}]=obs_seq_MOD_AOD_${DATE}.out
      fi
      export NL_NUM_INPUT_FILES=${NUM_FILES}
#
# All files present
      if [[ ${NL_NUM_INPUT_FILES} -eq 8 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\',\'${FILE_LIST[4]}\',\'${FILE_LIST[5]}\',\'${FILE_LIST[6]}\',\'${FILE_LIST[7]}\',\'${FILE_LIST[8]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 7 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\',\'${FILE_LIST[4]}\',\'${FILE_LIST[5]}\',\'${FILE_LIST[6]}\',\'${FILE_LIST[7]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 6 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\',\'${FILE_LIST[4]}\',\'${FILE_LIST[5]}\',\'${FILE_LIST[6]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 5 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\',\'${FILE_LIST[4]}\',\'${FILE_LIST[5]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 4 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\',\'${FILE_LIST[4]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 3 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\',\'${FILE_LIST[3]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 2 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\',\'${FILE_LIST[2]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 1 ]]; then
         export NL_FILENAME_SEQ=\'${FILE_LIST[1]}\'
      elif [[ ${NL_NUM_INPUT_FILES} -eq 0 ]]; then
         echo APM: ERROR no obs_seq files for FILTER
         exit
      fi
      export NL_FILENAME_OUT="'obs_seq.proc'"
      export NL_FIRST_OBS_DAYS=${ASIM_MIN_DAY_GREG}
      export NL_FIRST_OBS_SECONDS=${ASIM_MIN_SEC_GREG}
      export NL_LAST_OBS_DAYS=${ASIM_MAX_DAY_GREG}
      export NL_LAST_OBS_SECONDS=${ASIM_MAX_SEC_GREG}
      export NL_SYNONYMOUS_COPY_LIST="'NCEP BUFR observation','MOPITT CO observation','IASI CO observation','IASI O3 observation','AIRNOW observation','MODIS observation'"
      export NL_SYNONYMOUS_QC_LIST="'NCEP QC index','MOPITT CO QC index','IASI CO QC index','IASI O3 QC index','AIRNOW QC index','MODIS QC index'"
      rm -rf input.nml
      export NL_MOPITT_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_MOPITT}\'
      export NL_IASI_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_IASI}\'
      export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
      export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
      ${HYBRID_SCRIPTS_DIR}/da_create_dart_input_nml.ksh       
#
      ./obs_sequence_tool
      mv obs_seq.proc obs_seq_comb_${DATE}.out
   fi
#
#########################################################################
#
# RUN PREPROCESS OBSERVATIONS
#
#########################################################################
#
   if ${RUN_PREPROCESS_OBS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/preprocess_obs ]]; then
         mkdir ${RUN_DIR}/${DATE}/preprocess_obs
         cd ${RUN_DIR}/${DATE}/preprocess_obs
      else
         cd ${RUN_DIR}/${DATE}/preprocess_obs
      fi
      if [[ ${DATE} -eq ${INITIAL_DATE} ]]; then
         echo 'This is initial date cannot run PREPROCESS '
         touch CANNOT_RUN_PREPROCESS
      else
#
# GET WRFINPUT TEMPLATE
         cp ${RUN_DIR}/${PAST_DATE}/ensemble_mean_output/wrfout_d${CR_DOMAIN}_${DATE}_mean wrfinput_d${CR_DOMAIN}
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${CR_DOMAIN}_${FILE_DATE}.e001 wrfbiochemi_d${CR_DOMAIN}
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${FILE_DATE}.e001 wrffirechemi_d${CR_DOMAIN}
         cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${FILE_DATE}.e001 wrfchemi_d${CR_DOMAIN}
#
# GET DART UTILITIES
         cp ${DART_DIR}/models/wrf_chem/work/wrf_dart_obs_preprocess ./.
         cp ${DART_DIR}/models/wrf_chem/WRF_DART_utilities/wrf_dart_obs_preprocess.nml ./.
         rm -rf input.nml
         export NL_DEFAULT_STATE_VARIABLES=.false.
         export NL_MOPITT_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_MOPITT}\'
         export NL_IASI_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_IASI}\'
         export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
         export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
         ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
#
# GET INPUT DATA
         rm -rf obs_seq.old
         rm -rf obs_seq.new
         cp ${COMBINE_OBS_DIR}/obs_seq_comb_${DATE}.out obs_seq.old
#        
         ./wrf_dart_obs_preprocess ${DAY_GREG} ${SEC_GREG} > index_html 2>&1 
         export RC=$?     
         if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
         if [[ -f FAILED ]]; then rm -rf FAILED; fi          
         if [[ $RC = 0 ]]; then
            touch SUCCESS
         else
            touch FAILED 
            exit
         fi
         mv obs_seq.new obs_seq_comb_filtered_${DATE}.out 
      fi
   fi
