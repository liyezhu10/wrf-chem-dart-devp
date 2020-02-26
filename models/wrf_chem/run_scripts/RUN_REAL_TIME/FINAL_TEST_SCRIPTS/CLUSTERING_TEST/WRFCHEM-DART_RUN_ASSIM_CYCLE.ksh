#!/bin/ksh -aeux
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
# RUN WRFCHEM_CYCLE_CR
#
#########################################################################
#
   if ${RUN_WRFCHEM_CYCLE_CR}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_cycle_cr ]]; then
         mkdir -p ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
         cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
      fi
#
# Run WRF-Chem for all ensemble members
      TRANDOM=$$
      let IMEM=1
      export L_NUM_MEMBERS=${NUM_MEMBERS}
      if ${RUN_SPECIAL_FORECAST}; then
         export L_NUM_MEMBERS=${NUM_SPECIAL_FORECAST}
      fi
      while [[ ${IMEM} -le ${L_NUM_MEMBERS} ]]; do
         export MEM=${IMEM}
         export NL_TIME_STEP=${NNL_TIME_STEP}
         if ${RUN_SPECIAL_FORECAST}; then
            export MEM=${SPECIAL_FORECAST_MEM[${IMEM}]}
            let NL_TIME_STEP=${NNL_TIME_STEP}*${SPECIAL_FORECAST_FAC}
         fi
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
         export L_RUN_DIR=run_${CMEM}
         cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
         if ${RUN_SPECIAL_FORECAST}; then
            rm -rf ${L_RUN_DIR}
         fi
         if [[ ! -e ${L_RUN_DIR} ]]; then
            mkdir ${L_RUN_DIR}
            cd ${L_RUN_DIR}
         else
            cd ${L_RUN_DIR}
         fi
#
# Get WRF-Chem parameter files
         cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
         cp ${DART_DIR}/models/wrf_chem/work/input.nml ./.
         cp ${WRFCHEM_DIR}/test/em_real/wrf.exe ./.
         cp ${WRFCHEM_DIR}/test/em_real/aerosol.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/aerosol_lat.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/aerosol_lon.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/aerosol_plev.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/bulkdens.asc_s_0_03_0_9 ./.
         cp ${WRFCHEM_DIR}/test/em_real/bulkradii.asc_s_0_03_0_9 ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAM_ABS_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAM_AEROPT_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.A1B ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.A2 ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.RCP4.5 ./.
         cp ${WRFCHEM_DIR}/test/em_real/CAMtr_volume_mixing_ratio.RCP6 ./.
         cp ${WRFCHEM_DIR}/test/em_real/capacity.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/CCN_ACTIVATE.BIN ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_ALB_ICE_DFS_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_ALB_ICE_DRC_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_ASM_ICE_DFS_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_ASM_ICE_DRC_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_DRDSDT0_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_EXT_ICE_DFS_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_EXT_ICE_DRC_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_KAPPA_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/CLM_TAU_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/coeff_p.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/coeff_q.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/constants.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/ETAMPNEW_DATA.expanded_rain ./.
         cp ${WRFCHEM_DIR}/test/em_real/GENPARM.TBL ./.
         cp ${WRFCHEM_DIR}/test/em_real/grib2map.tbl ./.
         cp ${WRFCHEM_DIR}/test/em_real/gribmap.txt ./.
         cp ${WRFCHEM_DIR}/test/em_real/kernels.asc_s_0_03_0_9 ./.
         cp ${WRFCHEM_DIR}/test/em_real/kernels_z.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/LANDUSE.TBL ./.
         cp ${WRFCHEM_DIR}/test/em_real/masses.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/MPTABLE.TBL ./.
         cp ${WRFCHEM_DIR}/test/em_real/ozone.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/ozone_lat.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/ozone_plev.formatted ./.
         cp ${WRFCHEM_DIR}/test/em_real/RRTM_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/RRTMG_LW_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/RRTMG_SW_DATA ./.
         cp ${WRFCHEM_DIR}/test/em_real/SOILPARM.TBL ./.
         cp ${WRFCHEM_DIR}/test/em_real/termvels.asc ./.
         cp ${WRFCHEM_DIR}/test/em_real/tr49t67 ./.
         cp ${WRFCHEM_DIR}/test/em_real/tr49t85 ./.
         cp ${WRFCHEM_DIR}/test/em_real/tr67t85 ./.
         cp ${WRFCHEM_DIR}/test/em_real/URBPARM.TBL ./.
         cp ${WRFCHEM_DIR}/test/em_real/VEGPARM.TBL ./.
         rm -rf hist_io_flds*
         rm -rf clim_p_trop.nc
         rm -rf ubvals_b40.20th.track1_1996-2005.nc
         cp ${EXPERIMENT_HIST_IO_DIR}/hist_io_flds_v1 ./.
         cp ${EXPERIMENT_HIST_IO_DIR}/hist_io_flds_v2 ./.
#
         cp ${EXPERIMENT_STATIC_FILES}/clim_p_trop.nc ./.
         cp ${EXPERIMENT_STATIC_FILES}/ubvals_b40.20th.track1_1996-2005.nc ./.
         cp ${EXO_COLDENS_DIR}/exo_coldens_d${CR_DOMAIN} ./.
         cp ${SEASONS_WES_DIR}/wrf_season_wes_usgs_d${CR_DOMAIN}.nc ./.
#
# Get WRF-Chem emissions files
         export L_DATE=${START_DATE}
         while [[ ${L_DATE} -le ${END_DATE} ]]; do
            export L_YY=`echo ${L_DATE} | cut -c1-4`
            export L_MM=`echo ${L_DATE} | cut -c5-6`
            export L_DD=`echo ${L_DATE} | cut -c7-8`
            export L_HH=`echo ${L_DATE} | cut -c9-10`
            export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
            if [[ ${L_HH} -eq 00 || ${L_HH} -eq 06 || ${L_HH} -eq 12 || ${L_HH} -eq 18 ]]; then
               cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfbiochemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrfbiochemi_d${CR_DOMAIN}_${L_FILE_DATE}
            fi
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
            cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
            if [[ ${L_DATE} -eq ${START_DATE} ]]; then
               cp ${DART_FILTER_DIR}/wrk_dart_${CMEM}/wrfchemi_d${CR_DOMAIN} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
               cp ${DART_FILTER_DIR}/wrk_dart_${CMEM}/wrffirechemi_d${CR_DOMAIN} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
            else
               cp ${WRFCHEM_CHEM_EMISS_DIR}/wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
               cp ${WRFCHEM_CHEM_EMISS_DIR}/wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}.${CMEM} wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
            fi
            export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} +1 2>/dev/null)
         done
#
# Get WR-Chem input and bdy files
         cp ${DART_FILTER_DIR}/wrfout_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfinput_d${CR_DOMAIN}
         cp ${UPDATE_BC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}_filt.${CMEM} wrfbdy_d${CR_DOMAIN}
#
# APM SKIP_DART_TEST
#         cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfinput_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfinput_d${CR_DOMAIN}
#         cp ${WRFCHEM_CHEM_ICBC_DIR}/wrfbdy_d${CR_DOMAIN}_${START_FILE_DATE}.${CMEM} wrfbdy_d${CR_DOMAIN}
#
# Update the other emission files
         if ${ADD_EMISS}; then
            cp wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE} wrfchemi_d${CR_DOMAIN}_prior
            cp wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE} wrfchemi_d${CR_DOMAIN}
            cp wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE} wrffirechemi_d${CR_DOMAIN}_prior
            cp wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE} wrffirechemi_d${CR_DOMAIN}
            cp ${ADJUST_EMISS_DIR}/work/adjust_chem_emiss.exe ./.
            export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${START_DATE} +1 2>/dev/null)
            while [[ ${L_DATE} -le ${END_DATE} ]]; do 
               export L_YY=$(echo $L_DATE | cut -c1-4)
               export L_MM=$(echo $L_DATE | cut -c5-6)
               export L_DD=$(echo $L_DATE | cut -c7-8)
               export L_HH=$(echo $L_DATE | cut -c9-10)
               export L_FILE_DATE=${L_YY}-${L_MM}-${L_DD}_${L_HH}:00:00
#           
               export NL_WRFCHEMI_PRIOR=wrfchemi_d${CR_DOMAIN}_prior
               export NL_WRFCHEMI_POST=wrfchemi_d${CR_DOMAIN}
               export NL_WRFCHEMI_OLD=wrfchemi_d${CR_DOMAIN}_${L_FILE_DATE}
               export NL_WRFCHEMI_NEW=wrfchemi_d${CR_DOMAIN}_new
               cp ${NL_WRFCHEMI_OLD} ${NL_WRFCHEMI_NEW}
#           
               export NL_WRFFIRECHEMI_PRIOR=wrffirechemi_d${CR_DOMAIN}_prior
               export NL_WRFFIRECHEMI_POST=wrffirechemi_d${CR_DOMAIN}
               export NL_WRFFIRECHEMI_OLD=wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
               export NL_WRFFIRECHEMI_NEW=wrffirechemi_d${CR_DOMAIN}_new
               cp ${NL_WRFFIRECHEMI_OLD} ${NL_WRFFIRECHEMI_NEW}
#
# Make adjust_chem_nml for special_outlier_threshold
               rm -rf adjust_chem_emiss.nml
               cat <<  EOF > adjust_chem_emiss.nml
&adjust_chem_emiss
sw_chem=${NL_PERT_CHEM}
sw_fire=${NL_PERT_FIRE}
sw_biog=${NL_PERT_BIO}
fac=${EMISS_DAMP_CYCLE},
facc=${EMISS_DAMP_INTRA_CYCLE},
nx=${NNXP_CR},
ny=${NNYP_CR},
nz=${NNZP_CR},
nz_chemi=${NZ_CHEMI},
nz_firechemi=${NZ_FIRECHEMI},
nchemi_emiss=${NCHEMI_EMISS},
nfirechemi_emiss=${NFIRECHEMI_EMISS},
wrfchemi_prior='${NL_WRFCHEMI_PRIOR}',
wrfchemi_post='${NL_WRFCHEMI_POST}',
wrfchemi_old='${NL_WRFCHEMI_OLD}',
wrfchemi_new='${NL_WRFCHEMI_NEW}',
wrffirechemi_prior='${NL_WRFFIRECHEMI_PRIOR}',
wrffirechemi_post='${NL_WRFFIRECHEMI_POST}',
wrffirechemi_old='${NL_WRFFIRECHEMI_OLD}',
wrffirechemi_new='${NL_WRFFIRECHEMI_NEW}'
/
EOF
               RANDOM=$$
               export EXECUTE=adjust_chem_emiss.exe
               ./${EXECUTE} > index_${RANDOM}.html 2>&1
               cp ${NL_WRFCHEMI_NEW} ${NL_WRFCHEMI_OLD}
               cp ${NL_WRFFIRECHEMI_NEW} ${NL_WRFFIRECHEMI_OLD}
               export L_DATE=`echo ${L_DATE} +1h | ./advance_time`
            done
         fi
#
# Create WRF-Chem namelist.input 
         export NL_MAX_DOM=1
         export NL_IOFIELDS_FILENAME=\'hist_io_flds_v1\',\'hist_io_flds_v2\'
         rm -rf namelist.input
         ${HYBRID_SCRIPTS_DIR}/da_create_wrfchem_namelist_RT.ksh
         let IMEM=${IMEM}+1
      done
#
      let MEM=1
      while [[ ${MEM} -le ${L_NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         export KMEM=${MEM}
         if [[ ${MEM} -lt 1000 ]]; then export KMEM=0${MEM}; fi
         if [[ ${MEM} -lt 100 ]]; then export KMEM=00${MEM}; export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10 ]]; then export KMEM=000${MEM}; export CMEM=e00${MEM}; fi
         export L_RUN_DIR=run_${CMEM}
         cd ${RUN_DIR}/${DATE}/wrfchem_cycle_cr
         if ${RUN_SPECIAL_FORECAST}; then
            rm -rf ${L_RUN_DIR}
         fi
         if [[ ! -e ${L_RUN_DIR} ]]; then
            mkdir ${L_RUN_DIR}
            cd ${L_RUN_DIR}
         else
            cd ${L_RUN_DIR}
         fi
#
# RUN WRF-Chem
         export JOBRND=${TRANDOM}_wrf
         ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${WRFCHEM_JOB_CLASS} ${WRFCHEM_TIME_LIMIT} ${WRFCHEM_NODES} ${WRFCHEM_TASKS} wrf.exe PARALLEL ${ACCOUNT}
         sbatch job.ksh
         let MEM=${MEM}+1
      done
#
# Wait for all WRF-Chem jobs to complete for each member
      ${HYBRID_SCRIPTS_DIR}/da_run_hold_cu.ksh ${TRANDOM}
   fi
