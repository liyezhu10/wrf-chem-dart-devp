#!/bin/ksh -x
#
# Copyright 2019 University Corporation for Atmospheric Research and 
# Colorado Department of Public Health and Environment.
# 
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
# CONDITIONS OF ANY KIND, either express or implied. See the License for the 
# specific language governing permissions and limitations under the License.
# 
# Development of this code utilized the RMACC Summit supercomputer, which is 
# supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236),
# the University of Colorado Boulder, and Colorado State University. The Summit 
# supercomputer is a joint effort of the University of Colorado Boulder and 
# Colorado State University.

#########################################################################
#
# RUN DART_FILTER
#
#########################################################################
#
   if ${RUN_DART_FILTER}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/dart_filter ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/dart_filter
         cd ${RUN_DIR}/${DATE}/dart_filter
      else
         cd ${RUN_DIR}/${DATE}/dart_filter
      fi
echo $RETRIEVAL_TYPE_MOPITT
echo $RETRIEVAL_TYPE_IASI
exit

#
# Get DART files
      cp ${DART_DIR}/models/wrf_chem/work/filter      ./.
      cp ${DART_DIR}/system_simulation/final_full_precomputed_tables/final_full.${NUM_MEMBERS} ./.
      cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
      cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
#
# Get background forecasts
      if [[ ${DATE} -eq ${FIRST_FILTER_DATE} ]]; then
         export BACKGND_FCST_DIR=${WRFCHEM_INITIAL_DIR}
      else
         export BACKGND_FCST_DIR=${WRFCHEM_LAST_CYCLE_CR_DIR}
      fi
#
# Get observations
      if [[ ${PREPROCESS_OBS_DIR}/obs_seq_comb_filtered_${START_DATE}.out ]]; then      
         cp  ${PREPROCESS_OBS_DIR}/obs_seq_comb_filtered_${START_DATE}.out obs_seq.out
      else
         echo APM ERROR: NO DART OBSERVATIONS
         exit
      fi
#
# Run WRF_TO_DART
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
         cd ${RUN_DIR}/${DATE}/dart_filter
         rm -rf wrk_wrf_${CMEM}
         mkdir wrk_wrf_${CMEM}
         cd ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_${CMEM}
#
# &wrf_to_dart_nml
         export NL_DART_RESTART_NAME="'../filter_ic_old.${KMEM}'"
         export NL_PRINT_DATA_RANGES=.false.
         export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
         export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
         ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
         cp ${DART_DIR}/models/wrf_chem/work/wrf_to_dart ./.
##
## APM: +++ 
## For _ALL and emission inversion use wrfout instead of wrfapm because some needed fields 
## are not in wrfapm.  Also copy in the emissions files
         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${FILE_DATE} ../wrfinput_d${CR_DOMAIN}_${CMEM}
#         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfapm_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfapm_d${CR_DOMAIN}_${FILE_DATE} ../wrfinput_d${CR_DOMAIN}_${CMEM}
## APM: ---
##
         let MEM=${MEM}+1
      done
#
# APM: EMISSIONS
# APM: copy emission files for emission adjustments
      let IMEM=1
      while [[ ${IMEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${IMEM}
         export KMEM=${IMEM}
         if [[ ${IMEM} -lt 1000 ]]; then export KMEM=0${IMEM}; fi
         if [[ ${IMEM} -lt 100 ]]; then export KMEM=00${IMEM}; export CMEM=e0${IMEM}; fi
         if [[ ${IMEM} -lt 10 ]]; then export KMEM=000${IMEM}; export CMEM=e00${IMEM}; fi
         cd ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_${CMEM}      
         export LL_DATE=${DATE}
         export LL_END_DATE=${DATE}
         export LL_YY=`echo ${LL_DATE} | cut -c1-4`
         export LL_MM=`echo ${LL_DATE} | cut -c5-6`
         export LL_DD=`echo ${LL_DATE} | cut -c7-8`
         export LL_HH=`echo ${LL_DATE} | cut -c9-10`
         export LL_FILE_DATE=${LL_YY}-${LL_MM}-${LL_DD}_${LL_HH}:00:00
         if [[ ${LL_DATE} -le ${FIRST_EMISS_INV_DATE} || ${ADD_EMISS} = "false" ]]; then
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${LL_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${LL_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_CO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_NO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_co,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_no,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_oc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_bc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_c2h4,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch2o,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch3oh,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
         else
            cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfchemi_d${CR_DOMAIN}_${LL_FILE_DATE} wrfchemi_d${CR_DOMAIN}
            cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrffirechemi_d${CR_DOMAIN}_${LL_FILE_DATE} wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_CO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_NO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_co,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_no,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_oc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_bc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_c2h4,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch2o,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch3oh,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
         fi
        let IMEM=${IMEM}+1
      done
#
# Create job script 
      TRANDOM=$$
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
         cd ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_${CMEM}
#
         export JOBRND=${TRANDOM}_wrf2drt
         ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${GENERAL_JOB_CLASS} ${GENERAL_TIME_LIMIT} ${GENERAL_NODES} ${GENERAL_TASKS} wrf_to_dart SERIAL
         sbatch job.ksh
         let MEM=${MEM}+1
      done
#
# Wait for wrf_to_dart to complete for each member
      cd ${RUN_DIR}/${DATE}/dart_filter
      ${HYBRID_SCRIPTS_DIR}/da_run_hold_cu.ksh ${TRANDOM}
      ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
      cp ${EXPERIMENT_STATIC_FILES}/ubvals_b40.20th.track1_1996-2005.nc ./.
#
# APM: +++ another wrfapm / wrfout swap for emission inversion
      cp wrk_wrf_e001/wrfinput_d${CR_DOMAIN} ./
      cp wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} ./
      cp wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} ./
      cp ${BACKGND_FCST_DIR}/run_e001/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#      cp ${BACKGND_FCST_DIR}/run_e001/wrfapm_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
# APM: ---
#
# Copy "out" inflation files from prior cycle to "in" inflation files for current cycle
      if ${USE_DART_INFL}; then
         if [[ ${DATE} -eq ${FIRST_DART_INFLATE_DATE} ]]; then
            export NL_INF_INITIAL_FROM_RESTART_PRIOR=.false.
            export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.false.
            export NL_INF_INITIAL_FROM_RESTART_POST=.false.
            export NL_INF_SD_INITIAL_FROM_RESTART_POST=.false.
         else
            export NL_INF_INITIAL_FROM_RESTART_PRIOR=.true.
            export NL_INF_SD_INITIAL_FROM_RESTART_PRIOR=.true.
            export NL_INF_INITIAL_FROM_RESTART_POST=.true.
            export NL_INF_SD_INITIAL_FROM_RESTART_POST=.true.
         fi
         if [[ ${DATE} -ne ${FIRST_DART_INFLATE_DATE} ]]; then
            if [[ ${NL_INF_FLAVOR_PRIOR} != 0 ]]; then
               export INF_OUT_FILE_NAME_PRIOR=${RUN_DIR}/${PAST_DATE}/dart_filter/prior_inflate_ic_new
               cp ${INF_OUT_FILE_NAME_PRIOR} prior_inflate_ic_old
            fi
            if [[ ${NL_INF_FLAVOR_POST} != 0 ]]; then
               export INF_OUT_FILE_NAME_POST=${RUN_DIR}/${PAST_DATE}/dart_filter/post_inflate_ic_new
               cp ${NL_INF_OUT_FILE_NAME_POST} post_infalte_ic_old
            fi 
         fi
      fi
#
# Generate input.nml
      set -A temp `echo ${ASIM_MIN_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
      (( temp[1]=${temp[1]}+1 ))
      export NL_FIRST_OBS_DAYS=${temp[0]}
      export NL_FIRST_OBS_SECONDS=${temp[1]}
      set -A temp `echo ${ASIM_MAX_DATE} 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
      export NL_LAST_OBS_DAYS=${temp[0]}
      export NL_LAST_OBS_SECONDS=${temp[1]}
#
      export NL_NUM_INPUT_FILES=1
      export NL_FILENAME_SEQ="'obs_seq.out'"
      export NL_FILENAME_OUT="'obs_seq.processed'"
      export NL_MOPITT_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_MOPITT}\'
      export NL_IASI_CO_RETRIEVAL_TYPE=\'${RETRIEVAL_TYPE_IASI}\'
      export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
      export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
#
      rm -rf input.nml
      ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
#
# Make filter_apm_nml for special_outlier_threshold
      rm -rf filter_apm.nml
      cat << EOF > filter_apm.nml
&filter_apm_nml
special_outlier_threshold=${NL_SPECIAL_OUTLIER_THRESHOLD}
/
EOF
#
# Run DART_FILTER
# Create job script for this member and run it 
#
      RANDOM=$$
      export JOBRND=${RANDOM}_filter
      ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${FILTER_JOB_CLASS} ${FILTER_TIME_LIMIT} ${FILTER_NODES} ${FILTER_TASKS} filter PARALLEL
      sbatch -W job.ksh
#
# Check whether DART worked properly
     if [[ ! -f Prior_Diag.nc || ! -f Posterior_Diag.nc || ! -f obs_seq.final ]]; then
        echo APM: ERROR in DART FILTER EXIT
        exit
     fi
#
# Run DART_TO_WRF 
      TRANDOM=$$
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
         cd ${RUN_DIR}/${DATE}/dart_filter
         rm -rf wrk_dart_${CMEM}
         mkdir wrk_dart_${CMEM}
         cd ${RUN_DIR}/${DATE}/dart_filter/wrk_dart_${CMEM}
#
# &dart_to_wrf_nml
         export NL_MODEL_ADVANCE_FILE=.false.
         export NL_DART_RESTART_NAME=\'../filter_ic_new.${KMEM}\'
         export NL_USE_LOG_CO=${USE_LOG_CO_LOGIC}
         export NL_USE_LOG_O3=${USE_LOG_O3_LOGIC}
         rm -rf input.nml
         ${DART_DIR}/models/wrf_chem/namelist_scripts/DART/dart_create_input.nml.ksh
         cp ${DART_DIR}/models/wrf_chem/work/dart_to_wrf ./.
         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfout_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
#         cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfapm_d${CR_DOMAIN}_${FILE_DATE} wrfinput_d${CR_DOMAIN}
         export LL_DATE=${DATE}
         export LL_END_DATE=${DATE}
         export LL_YY=`echo ${LL_DATE} | cut -c1-4`
         export LL_MM=`echo ${LL_DATE} | cut -c5-6`
         export LL_DD=`echo ${LL_DATE} | cut -c7-8`
         export LL_HH=`echo ${LL_DATE} | cut -c9-10`
         export LL_FILE_DATE=${LL_YY}-${LL_MM}-${LL_DD}_${LL_HH}:00:00
         if [[ ${LL_DATE} -le ${FIRST_EMISS_INV_DATE} || ${ADD_EMISS} = "false" ]]; then
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${LL_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${LL_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_CO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_NO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_co,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_no,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_oc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_bc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_c2h4,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch2o,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch3oh,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
         else
            cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrfchemi_d${CR_DOMAIN}_${LL_FILE_DATE} wrfchemi_d${CR_DOMAIN}
            cp ${BACKGND_FCST_DIR}/run_${CMEM}/wrffirechemi_d${CR_DOMAIN}_${LL_FILE_DATE} wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_CO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,E_NO,c,c,"XLONG, XLAT" wrfchemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_co,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_no,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_oc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_bc,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_c2h4,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch2o,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
            ncatted -O -a coordinates,ebu_in_ch3oh,c,c,"XLONG, XLAT" wrffirechemi_d${CR_DOMAIN}
         fi
#
         export JOBRND=${TRANDOM}_drt2wrf
         ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${GENERAL_JOB_CLASS} ${GENERAL_TIME_LIMIT} ${GENERAL_NODES} ${GENERAL_TASKS} dart_to_wrf SERIAL
         sbatch job.ksh
#
         let MEM=${MEM}+1
      done
#
# Wait for dart_to_wrf to complete for each member
      cd ${RUN_DIR}/${DATE}/dart_filter
      ${HYBRID_SCRIPTS_DIR}/da_run_hold_cu.ksh ${TRANDOM}
#
# Copy converted output files
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
         cp wrk_dart_${CMEM}/wrfinput_d${CR_DOMAIN} wrfout_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM} 
         let MEM=${MEM}+1
      done
#
# Calculate the ensemble mean and spread adjustments for the emissions 
      if ${ADD_EMISS}; then 
         cd ${RUN_DIR}/${DATE}/dart_filter
         rm -rf perturb_chem_emiss_corr_nml.nl
         rm -rf perturb_emiss_chem_spec_nml.nl
         rm -rf perturb_chem_emiss_INCRs.exe
#
# Copy templates
         if ${NL_PERT_CHEM}; then
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_mean_prior 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_sprd_prior
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_mean_post
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_sprd_post 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_mean_incr 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_sprd_incr 
         fi
         if ${NL_PERT_FIRE}; then
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_mean_prior
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_sprd_prior
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_mean_post
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_sprd_post 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_mean_incr 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_sprd_incr 
         fi
         if ${NL_PERT_BIO}; then
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_mean_prior
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_sprd_prior
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_mean_post
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_sprd_post 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_mean_incr 
            cp ${RUN_DIR}/${DATE}/dart_filter/wrk_wrf_e001/wrfbiochemi_d${CR_DOMAIN} wrfbiochemi_d${CR_DOMAIN}_sprd_incr 
         fi
#        
         cp ${RUN_DIR}/${DATE}/wrfchem_chem_emiss/perturb_chem_emiss_corr_nml.nl ./      
         cp ${RUN_DIR}/${DATE}/wrfchem_chem_emiss/perturb_emiss_chem_spec_nml.nl ./      
         cp ${PERT_CHEM_EMISS_DIR}/work/perturb_chem_emiss_INCRs.exe ./
#
# SERIAL VERSION
         RANDOM=$$
         export JOBRND=${RANDOM}_cr_emiss_pert
         ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${GENERAL_JOB_CLASS} ${GENERAL_TIME_LIMIT} ${GENERAL_NODES} ${GENERAL_TASKS} perturb_chem_emiss_INCRs.exe SERIAL
         sbatch -W job.ksh
      fi
   fi
#
#########################################################################
#
# UPDATE COARSE RESOLUTION BOUNDARY CONDIIONS
#
#########################################################################
#
   if ${RUN_UPDATE_BC}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/update_bc ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/update_bc
         cd ${RUN_DIR}/${DATE}/update_bc
      else
         cd ${RUN_DIR}/${DATE}/update_bc
      fi
#
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
#
         export CYCLING=true
         export OPS_FORC_FILE=${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${FILE_DATE}.${CMEM}
         export BDYCDN_IN=${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${FILE_DATE}.${CMEM}

         cp ${BDYCDN_IN} wrfbdy_d${CR_DOMAIN}_${FILE_DATE}_prior.${CMEM}
         export DA_OUTPUT_FILE=${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM} 
         export BDYCDN_OUT=wrfbdy_d${CR_DOMAIN}_${FILE_DATE}_filt.${CMEM}    
         ${HYBRID_SCRIPTS_DIR}/da_run_update_bc.ksh > index_update_bc 2>&1
#
         let MEM=$MEM+1
      done
   fi
