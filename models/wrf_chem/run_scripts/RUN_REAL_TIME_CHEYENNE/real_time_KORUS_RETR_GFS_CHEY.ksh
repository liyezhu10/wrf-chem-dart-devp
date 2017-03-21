#!/bin/ksh -aux
#########################################################################
#
# Purpose: Set global environment variables for real_time_wrf_chem
#
# NOTE: To generate the chemistry perturbations check NL_SW_GENERATE comments
# To use the same purturbastions from one cycle to the next one must copy
# the pert_file_emiss and pert_file_icbc file to the run directory
#
#########################################################################
#
# CYCLE DATE-TIME:
export CYCLE_STR_DATE=2016052500
export CYCLE_STR_DATE=2016052500
export CYCLE_END_DATE=${CYCLE_STR_DATE}
export CYCLE_END_DATE=2016052500
export CYCLE_DATE=${CYCLE_STR_DATE}
export RETRIEVAL_TYPE=RETR
export ADD_EMISS=.true.
export VARLOC=.false.
export INDEP_CHEM_ASIM=.false.
#
# Run fine scale forecast only
export RUN_FINE_SCALE=false
#
# Restart fine scale forecast only
export RUN_FINE_SCALE_RESTART=false
export RESTART_DATE=2014072312
#
if [[ ${RUN_FINE_SCALE_RESTART} = "true" ]]; then
   export RUN_FINE_SCALE=true
fi
#
# Run WRF-Chem for failed forecasts
export RUN_SPECIAL_FORECAST=false
export NUM_SPECIAL_FORECAST=1
export SPECIAL_FORECAST_FAC=1.
export SPECIAL_FORECAST_FAC=2./3.
export SPECIAL_FORECAST_MEM[1]=18
export SPECIAL_FORECAST_MEM[2]=23
export SPECIAL_FORECAST_MEM[3]=3
export SPECIAL_FORECAST_MEM[4]=4
export SPECIAL_FORECAST_MEM[5]=5
export SPECIAL_FORECAST_MEM[6]=6
export SPECIAL_FORECAST_MEM[7]=7
export SPECIAL_FORECAST_MEM[8]=8
export SPECIAL_FORECAST_MEM[9]=9
export SPECIAL_FORECAST_MEM[10]=10
export SPECIAL_FORECAST_MEM[11]=11
export SPECIAL_FORECAST_MEM[12]=12
export SPECIAL_FORECAST_MEM[13]=13
export SPECIAL_FORECAST_MEM[14]=14
export SPECIAL_FORECAST_MEM[15]=15
export SPECIAL_FORECAST_MEM[16]=16
export SPECIAL_FORECAST_MEM[17]=17
export SPECIAL_FORECAST_MEM[18]=18
export SPECIAL_FORECAST_MEM[19]=19
export SPECIAL_FORECAST_MEM[20]=20
export SPECIAL_FORECAST_MEM[21]=21
export SPECIAL_FORECAST_MEM[22]=22
export SPECIAL_FORECAST_MEM[23]=23
export SPECIAL_FORECAST_MEM[24]=24
export SPECIAL_FORECAST_MEM[25]=25
export SPECIAL_FORECAST_MEM[26]=26
export SPECIAL_FORECAST_MEM[27]=27
export SPECIAL_FORECAST_MEM[28]=28
export SPECIAL_FORECAST_MEM[29]=29
export SPECIAL_FORECAST_MEM[30]=30
#
# Run temporal interpolation for missing background files
export RUN_INTERPOLATE=false
#
# for 2014072212 and 2014072218
#export BACK_DATE=2014072206
#export FORW_DATE=2014072300
#let BACK_WT=.3333
#let BACK_WT=.6667
#
# for 20142900
#export BACK_DATE=2014072818
#export FORW_DATE=2014072906
#let BACK_WT=.5000
#
# for 20142912
#export BACK_DATE=2014072906
#export FORW_DATE=2014072918
#let BACK_WT=.5000
#
while [[ ${CYCLE_DATE} -le ${CYCLE_END_DATE} ]]; do
export DATE=${CYCLE_DATE}
export INITIAL_DATE=2016052500
export FIRST_FILTER_DATE=2016052506
export FIRST_EMISS_INV_DATE=2016052512
export FIRST_DART_INFLATE_DATE=2016052512
export CYCLE_PERIOD=6
export HISTORY_INTERVAL_HR=1
(( HISTORY_INTERVAL_MIN = ${HISTORY_INTERVAL_HR} * 60 ))
echo ${HISTORY_INTERVAL_MIN}
export START_IASI_O3_DATA=2016052500
export END_IASI_O3_DATA=2016060000
export NL_DEBUG_LEVEL=200
#
# CODE VERSIONS:
export WPS_VER=WPSv3.6.1_dmpar
export WPS_GEOG_VER=WPSv3.6.1_GEOG_DATA
export WRFDA_VER=WRFDAv3.4_dmpar
export WRF_VER=WRFv3.6.1_dmpar
export WRFCHEM_VER=WRFCHEMv3.6.1_dmpar
export WRFCHEM_VER_GABI=WRFCHEMv3.6.1_dmpar_GABI_REV1.0
export DART_VER=DART_CHEM_REPOSITORY
#
# ROOT DIRECTORIES:
export SCRATCH_DIR=/glade/scratch/mizzi
export WORK_DIR=/glade/p/work/mizzi
export ACD_DIR=/glade/p/acd/mizzi
export FRAPPE_DIR=/glade/p/FRAPPE
export EXPERIMENT_OUTPUT_DIR=${SCRATCH_DIR}
#
# DEPENDENT INPUT DATA DIRECTORIES:
export EXPERIMENT_DIR=${FRAPPE_DIR}
export RUN_DIR=${EXPERIMENT_DIR}/real_KORUS_RETR_VARLOC_CHEYENE
export TRUNK_DIR=${WORK_DIR}/TRUNK
export WPS_DIR_v3p4p1=/glade/p/work/wrfhelp/PRE_COMPILED_CODE/WPSV3.4.1_intel_dmpar
export WPS_DIR_v3p7p1=/glade/p/work/wrfhelp/PRE_COMPILED_CODE/WPSV3.7.1_intel_serial_large-file
export WPS_DIR_v3p8p1=/glade/p/work/wrfhelp/PRE_COMPILED_CODE/WPSV3.8.1_intel_serial_large-file
export WPS_DIR=${TRUNK_DIR}/${WPS_VER}
export WPS_GEOG_DIR=${TRUNK_DIR}/${WPS_GEOG_VER}/geog
export WRFCHEM_DIR=${TRUNK_DIR}/${WRFCHEM_VER}
export WRFCHEM_DIR_GABI=${FRAPPE_DIR}/FRAPPE_TRUNK/${WRFCHEM_VER_GABI}
export WRFDA_DIR=${TRUNK_DIR}/${WRFDA_VER}
export DART_DIR=${TRUNK_DIR}/${DART_VER}
export BUILD_DIR=${WRFDA_DIR}/var/da
export WRF_DIR=${TRUNK_DIR}/${WRF_VER}
export HYBRID_SCRIPTS_DIR=${DART_DIR}/models/wrf_chem/hybrid_scripts
export ADJUST_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_EMISS_INV
export EXPERIMENT_DATA_DIR=${EXPERIMENT_DIR}/KORUS_REAL_TIME_DATA
export MOZBC_DATA_DIR=${EXPERIMENT_DIR}/KORUS_REAL_TIME_DATA/cam_chem_forecasts
export EXPERIMENT_STATIC_FILES=${EXPERIMENT_DATA_DIR}/static_files
export EXPERIMENT_WRFCHEMI_DIR=${EXPERIMENT_DATA_DIR}/anthro_emissions
export EXPERIMENT_WRFFIRECHEMI_DIR=${EXPERIMENT_DATA_DIR}/fire_emissions
export EXPERIMENT_WRFBIOCHEMI_DIR=${EXPERIMENT_DATA_DIR}/bio_emissions
export EXPERIMENT_COLDENS_DIR=${EXPERIMENT_DATA_DIR}/wes_coldens
export EXPERIMENT_PREPBUFR_DIR=${EXPERIMENT_DATA_DIR}/met_obs_prep_data
export EXPERIMENT_MOPITT_CO_DIR=${EXPERIMENT_DATA_DIR}/mopitt_co_hdf_data
export EXPERIMENT_IASI_CO_DIR=${EXPERIMENT_DATA_DIR}/iasi_co_hdf_data
export EXPERIMENT_IASI_O3_DIR=${EXPERIMENT_DATA_DIR}/iasi_o3_hdf_data
export EXPERIMENT_AIRNOW_DIR=${EXPERIMENT_DATA_DIR}/airnow_csv_data
export EXPERIMENT_MODIS_AOD_DIR=${EXPERIMENT_DATA_DIR}/modis_aod_hdf_data
export EXPERIMENT_GFS_DIR=${EXPERIMENT_DATA_DIR}/gfs_forecasts
export EXPERIMENT_DUST_DIR=${EXPERIMENT_DATA_DIR}/dust_fields
export EXPERIMENT_HIST_IO_DIR=${EXPERIMENT_DATA_DIR}/hist_io_files
export VTABLE_DIR=${WPS_DIR}/ungrib/Variable_Tables
export VTABLE_DIR=${WPS_DIR_v3p7p1}/ungrib/Variable_Tables
export BE_DIR=${WRFDA_DIR}/var/run
export PERT_CHEM_INPUT_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/ICBC_PERT
export PERT_CHEM_EMISS_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_PERT_CHEM/EMISS_PERT
#
cp ${DART_DIR}/models/wrf_chem/work/advance_time ./.
cp ${DART_DIR}/models/wrf_chem/work/input_KORUS.nml ./input.nml
export YYYY=$(echo $DATE | cut -c1-4)
export YY=$(echo $DATE | cut -c3-4)
export MM=$(echo $DATE | cut -c5-6)
export DD=$(echo $DATE | cut -c7-8)
export HH=$(echo $DATE | cut -c9-10)
export DATE_SHORT=${YY}${MM}${DD}${HH}
export FILE_DATE=${YYYY}-${MM}-${DD}_${HH}:00:00
export PAST_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -${CYCLE_PERIOD} 2>/dev/null)
export PAST_YYYY=$(echo $PAST_DATE | cut -c1-4)
export PAST_YY=$(echo $PAST_DATE | cut -c3-4)
export PAST_MM=$(echo $PAST_DATE | cut -c5-6)
export PAST_DD=$(echo $PAST_DATE | cut -c7-8)
export PAST_HH=$(echo $PAST_DATE | cut -c9-10)
export PAST_FILE_DATE=${PAST_YYYY}-${PAST_MM}-${PAST_DD}_${PAST_HH}:00:00
export NEXT_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} +${CYCLE_PERIOD} 2>/dev/null)
export NEXT_YYYY=$(echo $NEXT_DATE | cut -c1-4)
export NEXT_YY=$(echo $NEXT_DATE | cut -c3-4)
export NEXT_MM=$(echo $NEXT_DATE | cut -c5-6)
export NEXT_DD=$(echo $NEXT_DATE | cut -c7-8)
export NEXT_HH=$(echo $NEXT_DATE | cut -c9-10)
export NEXT_FILE_DATE=${NEXT_YYYY}-${NEXT_MM}-${NEXT_DD}_${NEXT_HH}:00:00
#
# DART TIME DATA
export DT_YYYY=${YYYY}
export DT_YY=${YY}
export DT_MM=${MM} 
export DT_DD=${DD} 
export DT_HH=${HH} 
let DDT_MM=DT_MM+0
echo ${DDT_MM}
export DT_MM=${DDT_MM}
#
let DDT_DD=DT_DD+0
export DT_DD=${DDT_DD}
#
let DDT_HH=DT_HH+0 
export DT_HH=${DDT_HH}
if [[ ${HH} -eq 0 ]]; then
   cp ${BUILD_DIR}/da_advance_time.exe ./
   ./da_advance_time.exe ${DATE} -1
#
#   export TMP_DATE=$(${BUILD_DIR}/da_advance_time.exe ${DATE} -1 2>/dev/null)
#echo $TMP_DATE
exit
   export TMP_YYYY=$(echo $TMP_DATE | cut -c1-4)
   export TMP_YY=$(echo $TMP_DATE | cut -c3-4)
   export TMP_MM=$(echo $TMP_DATE | cut -c5-6)
   export TMP_DD=$(echo $TMP_DATE | cut -c7-8)
   export TMP_HH=$(echo $TMP_DATE | cut -c9-10)
   export D_YYYY=${TMP_YYYY}
   export D_YY=${TMP_YY}
   export D_MM=${TMP_MM}
   export D_DD=${TMP_DD}
   export D_HH=24
   (( DD_MM = ${D_MM} + 0 ))
   (( DD_DD = ${D_DD} + 0 ))
   (( DD_HH = ${D_HH} + 0 ))
else
   export D_YYYY=${YYYY}
   export D_YY=${YY}
   export D_MM=${MM}
   export D_DD=${DD}
   export D_HH=${HH}
   (( DD_MM = ${D_MM} + 0 ))
   (( DD_DD = ${D_DD} + 0 ))
   (( DD_HH = ${D_HH} + 0 ))
fi
export D_DATE=${D_YYYY}${D_MM}${D_DD}${D_HH}
exit
#
# CALCULATE GREGORIAN TIMES FOR START AND END OF ASSIMILATION WINDOW
set -A GREG_DATA `echo $DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export DAY_GREG=${GREG_DATA[0]}
export SEC_GREG=${GREG_DATA[1]}
set -A GREG_DATA `echo $NEXT_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export NEXT_DAY_GREG=${GREG_DATA[0]}
export NEXT_SEC_GREG=${GREG_DATA[1]}
export ASIM_WINDOW=3
export ASIM_MIN_DATE=$($BUILD_DIR/da_advance_time.exe $DATE -$ASIM_WINDOW 2>/dev/null)
export ASIM_MIN_YYYY=$(echo $ASIM_MIN_DATE | cut -c1-4)
export ASIM_MIN_YY=$(echo $ASIM_MIN_DATE | cut -c3-4)
export ASIM_MIN_MM=$(echo $ASIM_MIN_DATE | cut -c5-6)
export ASIM_MIN_DD=$(echo $ASIM_MIN_DATE | cut -c7-8)
export ASIM_MIN_HH=$(echo $ASIM_MIN_DATE | cut -c9-10)
export ASIM_MAX_DATE=$($BUILD_DIR/da_advance_time.exe $DATE +$ASIM_WINDOW 2>/dev/null)
export ASIM_MAX_YYYY=$(echo $ASIM_MAX_DATE | cut -c1-4)
export ASIM_MAX_YY=$(echo $ASIM_MAX_DATE | cut -c3-4)
export ASIM_MAX_MM=$(echo $ASIM_MAX_DATE | cut -c5-6)
export ASIM_MAX_DD=$(echo $ASIM_MAX_DATE | cut -c7-8)
export ASIM_MAX_HH=$(echo $ASIM_MAX_DATE | cut -c9-10)
set -A temp `echo $ASIM_MIN_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time`
export ASIM_MIN_DAY_GREG=${temp[0]}
export ASIM_MIN_SEC_GREG=${temp[1]}
set -A temp `echo $ASIM_MAX_DATE 0 -g | ${DART_DIR}/models/wrf_chem/work/advance_time` 
export ASIM_MAX_DAY_GREG=${temp[0]}
export ASIM_MAX_SEC_GREG=${temp[1]}
#
# SELECT COMPONENT RUN OPTIONS:
if [[ ${RUN_SPECIAL_FORECAST} = "false" ]]; then
   export RUN_GEOGRID=true
   export RUN_UNGRIB=false
   export RUN_METGRID=false
   export RUN_REAL=false
   export RUN_PERT_WRFCHEM_MET_IC=false
   export RUN_PERT_WRFCHEM_MET_BC=false
   export RUN_EXO_COLDENS=false
   export RUN_SEASON_WES=false
   export RUN_WRFCHEM_BIO=false
   export RUN_WRFCHEM_FIRE=false
   export RUN_WRFCHEM_CHEMI=false
   export RUN_PERT_WRFCHEM_CHEM_ICBC=false
   export RUN_PERT_WRFCHEM_CHEM_EMISS=false
   export RUN_MOPITT_CO_OBS=true
   export RUN_IASI_CO_OBS=false
   export RUN_IASI_O3_OBS=false
   export RUN_AIRNOW_O3_OBS=false
   export RUN_AIRNOW_CO_OBS=false
   export RUN_MODIS_AOD_OBS=false
   export RUN_MET_OBS=true
   export RUN_COMBINE_OBS=true
   export RUN_PREPROCESS_OBS=true
#
   if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
      export RUN_WRFCHEM_INITIAL=true
      export RUN_DART_FILTER=false
      export RUN_UPDATE_BC=false
      export RUN_WRFCHEM_CYCLE_CR=false
      export RUN_WRFCHEM_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_INPUT=false
      export RUN_ENSMEAN_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_OUTPUT=false
   else
      export RUN_WRFCHEM_INITIAL=false
      export RUN_DART_FILTER=true
      export RUN_UPDATE_BC=true
      export RUN_WRFCHEM_CYCLE_CR=true
      export RUN_WRFCHEM_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_INPUT=true
      export RUN_ENSMEAN_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_OUTPUT=true
   fi
else
   export RUN_GEOGRID=false
   export RUN_UNGRIB=false
   export RUN_METGRID=false
   export RUN_REAL=false
   export RUN_PERT_WRFCHEM_MET_IC=false
   export RUN_PERT_WRFCHEM_MET_BC=false
   export RUN_EXO_COLDENS=false
   export RUN_SEASON_WES=false
   export RUN_WRFCHEM_BIO=false
   export RUN_WRFCHEM_FIRE=false
   export RUN_WRFCHEM_CHEMI=false
   export RUN_PERT_WRFCHEM_CHEM_ICBC=false
   export RUN_PERT_WRFCHEM_CHEM_EMISS=false
   export RUN_MOPITT_CO_OBS=false
   export RUN_IASI_CO_OBS=false
   export RUN_IASI_O3_OBS=false
   export RUN_AIRNOW_O3_OBS=false
   export RUN_AIRNOW_CO_OBS=false
   export RUN_MODIS_AOD_OBS=false
   export RUN_MET_OBS=false
   export RUN_COMBINE_OBS=false
   export RUN_PREPROCESS_OBS=false
#
   if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
      export RUN_WRFCHEM_INITIAL=true
      export RUN_DART_FILTER=false
      export RUN_UPDATE_BC=false
      export RUN_WRFCHEM_CYCLE_CR=false
      export RUN_WRFCHEM_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_INPUT=false
      export RUN_ENSMEAN_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_OUTPUT=false
   else
      export RUN_WRFCHEM_INITIAL=false
      export RUN_DART_FILTER=false
      export RUN_UPDATE_BC=false
      export RUN_WRFCHEM_CYCLE_CR=true
      export RUN_WRFCHEM_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_INPUT=true
      export RUN_ENSMEAN_CYCLE_FR=false
      export RUN_ENSEMBLE_MEAN_OUTPUT=true
   fi
fi
if [[ ${RUN_FINE_SCALE} = "true" ]]; then
   export RUN_GEOGRID=false
   export RUN_UNGRIB=false
   export RUN_METGRID=false
   export RUN_REAL=false
   export RUN_PERT_WRFCHEM_MET_IC=false
   export RUN_PERT_WRFCHEM_MET_BC=false
   export RUN_EXO_COLDENS=false
   export RUN_SEASON_WES=false
   export RUN_WRFCHEM_BIO=false
   export RUN_WRFCHEM_FIRE=false
   export RUN_WRFCHEM_CHEMI=false
   export RUN_PERT_WRFCHEM_CHEM_ICBC=false
   export RUN_PERT_WRFCHEM_CHEM_EMISS=false
   export RUN_MOPITT_CO_OBS=false
   export RUN_IASI_CO_OBS=false
   export RUN_IASI_O3_OBS=false
   export RUN_AIRNOW_O3_OBS=false
   export RUN_AIRNOW_CO_OBS=false
   export RUN_MODIS_AOD_OBS=false
   export RUN_MET_OBS=false
   export RUN_COMBINE_OBS=false
   export RUN_PREPROCESS_OBS=false
   export RUN_WRFCHEM_INITIAL=false
   export RUN_DART_FILTER=false
   export RUN_UPDATE_BC=false
   export RUN_WRFCHEM_CYCLE_CR=false
   export RUN_WRFCHEM_CYCLE_FR=false
   export RUN_ENSEMBLE_MEAN_INPUT=false
   export RUN_ENSMEAN_CYCLE_FR=true
   export RUN_ENSEMBLE_MEAN_OUTPUT=false
fi
#
# FORECAST PARAMETERS:
export USE_DART_INFL=true
export FCST_PERIOD=6
(( CYCLE_PERIOD_SEC=${CYCLE_PERIOD}*60*60 ))
export NUM_MEMBERS=10
export MAX_DOMAINS=02
export CR_DOMAIN=01
export FR_DOMAIN=02
export NNXP_CR=139
export NNYP_CR=139
export NNZP_CR=36
export NNXP_FR=230
export NNYP_FR=230
export NNZP_FR=36
(( NNXP_STAG_CR=${NNXP_CR}+1 ))
(( NNYP_STAG_CR=${NNYP_CR}+1 ))
(( NNZP_STAG_CR=${NNZP_CR}+1 ))
(( NNXP_STAG_FR=${NNXP_FR}+1 ))
(( NNYP_STAG_FR=${NNYP_FR}+1 ))
(( NNZP_STAG_FR=${NNZP_FR}+1 ))
export NNZ_CHEM=6
export NNCHEM_SPC=54
export NNFIRE_SPC=31
export NNBIO_SPC=1
export NZ_CHEMI=${NNZ_CHEM}
export NZ_FIRECHEMI=1
export NCHEMI_EMISS=2
export NFIRECHEMI_EMISS=7
export ISTR_CR=1
export JSTR_CR=1
export ISTR_FR=47
export JSTR_FR=32
export DX_CR=15000
export DX_FR=3000
(( LBC_END=2*${FCST_PERIOD} ))
export LBC_FREQ=3
(( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
export LBC_START=0
export START_DATE=${DATE}
export END_DATE=$($BUILD_DIR/da_advance_time.exe ${START_DATE} ${FCST_PERIOD} 2>/dev/null)
export START_YEAR=$(echo $START_DATE | cut -c1-4)
export START_YEAR_SHORT=$(echo $START_DATE | cut -c3-4)
export START_MONTH=$(echo $START_DATE | cut -c5-6)
export START_DAY=$(echo $START_DATE | cut -c7-8)
export START_HOUR=$(echo $START_DATE | cut -c9-10)
export START_FILE_DATE=${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:00:00
export END_YEAR=$(echo $END_DATE | cut -c1-4)
export END_MONTH=$(echo $END_DATE | cut -c5-6)
export END_DAY=$(echo $END_DATE | cut -c7-8)
export END_HOUR=$(echo $END_DATE | cut -c9-10)
export END_FILE_DATE=${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:00:00
#
# LARGE SCALE FORECAST PARAMETERS:
export FG_TYPE=GFS
export GRIB_PART1=gfs_4_
export GRIB_PART2=.g2.tar
#
# COMPUTER PARAMETERS:
export PROJ_NUMBER_ACD=P19010000
export PROJ_NUMBER_NSC=NACD0009
export GENERAL_QUEUE=regular
export GEOGRID_TIME_LIMIT=0:10
export GEOGRID_NODES=1
export GEOGRID_CPUS=8
export GEOGRID_PROCS=8
export GEOGRID_QUEUE=regular
export WRFCHEM_TIME_LIMIT=6:00
export WRFCHEM_TIME_LIMIT_LONG=10:00
export WRFCHEM_NODES=2
export WRFCHEM_CPUS=36
export WRFCHEM_PROCS=36
export WRFCHEM_QUEUE=regular
export WRFDA_TIME_LIMIT=0:10
export WRFDA_NODES=2
export WRFDA_CPUS=36
export WRFDA_PROCS=36
export WRFDA_QUEUE=regular
export FILTER_TIME_LIMIT=3:59
export FILTER_NODES=2
export FILTER_CPUS=36
export FILTER_PROCS=36
export FILTER_JOB_CLASS=regular
export MISC_TIME_LIMIT=0:02
export MISC_NODES=1
export MISC_CPUS=36
export MISC_PROCS=36
export MISC_QUEUE=regular
#
# RUN DIRECTORIES
export GEOGRID_DIR=${RUN_DIR}/geogrid
export METGRID_DIR=${RUN_DIR}/${DATE}/metgrid
export REAL_DIR=${RUN_DIR}/${DATE}/real
export WRFCHEM_MET_IC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_ic
export WRFCHEM_MET_BC_DIR=${RUN_DIR}/${DATE}/wrfchem_met_bc
export EXO_COLDENS_DIR=${RUN_DIR}/${DATE}/exo_coldens
export SEASONS_WES_DIR=${RUN_DIR}/${DATE}/seasons_wes
export WRFCHEM_BIO_DIR=${RUN_DIR}/${DATE}/wrfchem_bio
export WRFCHEM_FIRE_DIR=${RUN_DIR}/${DATE}/wrfchem_fire
export WRFCHEM_CHEMI_DIR=${RUN_DIR}/${DATE}/wrfchem_chemi
export WRFCHEM_CHEM_EMISS_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
export WRFCHEM_INITIAL_DIR=${RUN_DIR}/${INITIAL_DATE}/wrfchem_initial
export WRFCHEM_CYCLE_CR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_cr
export WRFCHEM_CYCLE_FR_DIR=${RUN_DIR}/${DATE}/wrfchem_cycle_fr
export WRFCHEM_LAST_CYCLE_CR_DIR=${RUN_DIR}/${PAST_DATE}/wrfchem_cycle_cr
export PREPBUFR_MET_OBS_DIR=${RUN_DIR}/${DATE}/prepbufr_met_obs
export MOPITT_CO_OBS_DIR=${RUN_DIR}/${DATE}/mopitt_co_obs
export IASI_CO_OBS_DIR=${RUN_DIR}/${DATE}/iasi_co_obs
export IASI_O3_OBS_DIR=${RUN_DIR}/${DATE}/iasi_o3_obs
export AIRNOW_O3_OBS_DIR=${RUN_DIR}/${DATE}/airnow_o3_obs
export AIRNOW_CO_OBS_DIR=${RUN_DIR}/${DATE}/airnow_co_obs
export MODIS_AOD_OBS_DIR=${RUN_DIR}/${DATE}/modis_aod_obs
export COMBINE_OBS_DIR=${RUN_DIR}/${DATE}/combine_obs
export PREPROCESS_OBS_DIR=${RUN_DIR}/${DATE}/preprocess_obs
export WRFCHEM_CHEM_ICBC_DIR=${RUN_DIR}/${DATE}/wrfchem_chem_icbc
export DART_FILTER_DIR=${RUN_DIR}/${DATE}/dart_filter
export UPDATE_BC_DIR=${RUN_DIR}/${DATE}/update_bc
export ENSEMBLE_MEAN_INPUT_DIR=${RUN_DIR}/${DATE}/ensemble_mean_input
export ENSEMBLE_MEAN_OUTPUT_DIR=${RUN_DIR}/${DATE}/ensemble_mean_output
export REAL_TIME_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_REAL_TIME
#
# WPS PARAMETERS:
export SINGLE_FILE=false
export HOR_SCALE=1500
export VTABLE_TYPE=GFS
export METGRID_TABLE_TYPE=ARW
#
# WRF PREPROCESS PARAMETERS
export NL_MIN_LAT=26.
export NL_MAX_LAT=48.
export NL_MIN_LON=114.
export NL_MAX_LON=140.
export NNL_MIN_LAT=${NL_MIN_LAT}
export NNL_MAX_LAT=${NL_MAX_LAT}
export NNL_MIN_LON=${NL_MIN_LON}
if [[ ${NL_MIN_LON} -gt 180. ]]; then ((NNL_MIN_LON=${NL_MIN_LON}-360.)); fi 
export NNL_MAX_LON=${NL_MAX_LON}
if [[ ${NL_MAX_LON} -gt 180. ]]; then ((NNL_MAX_LON=${NL_MAX_LON}-360.)); fi 
export NL_OBS_PRESSURE_TOP=10000.
#
# TARG_LAT=38.92 (91,61)
# TARG_LON=-111.52 = 248.48 (91,61)
#export NL_MIN_LAT=38.32
#export NL_MAX_LAT=39.52
#export NL_MIN_LON=248.08
#export NL_MAX_LON=248.88
#
# PERT CHEM PARAMETERS
export SPREAD_FAC=0.30
export MOZ_SPREAD=${SPREAD_FAC}
export NL_MEAN=1.0
export NL_SPREAD=${SPREAD_FAC}
#
#########################################################################
#
#  NAMELIST PARAMETERS
#
#########################################################################
#
# WPS SHARE NAMELIST:
export NL_WRF_CORE=\'ARW\'
export NL_MAX_DOM=${MAX_DOMAINS}
export NL_IO_FORM_GEOGRID=2
export NL_OPT_OUTPUT_FROM_GEOGRID_PATH=\'${GEOGRID_DIR}\'
export NL_ACTIVE_GRID=".true.",".true."
#
# WPS GEOGRID NAMELIST:
export NL_S_WE=1,1
export NL_E_WE=${NNXP_STAG_CR},${NNXP_STAG_FR}
export NL_S_SN=1,1
export NL_E_SN=${NNYP_STAG_CR},${NNYP_STAG_FR}
export NL_S_VERT=1,1
export NL_E_VERT=${NNZP_STAG_CR},${NNZP_STAG_FR}
export NL_PARENT_ID="0,1"
export NL_PARENT_GRID_RATIO=1,5
export NL_I_PARENT_START=${ISTR_CR},${ISTR_FR}
export NL_J_PARENT_START=${JSTR_CR},${JSTR_FR}
export NL_GEOG_DATA_RES=\'30s\',\'30s\'
export NL_DX=${DX_CR}
export NL_DY=${DX_CR}
export NL_MAP_PROJ=\'mercator\'
export NL_REF_LAT=37.5
export NL_REF_LON=127.0
export NL_STAND_LON=127.0
export NL_TRUELAT1=37.0
export NL_TRUELAT2=60.0
export NL_GEOG_DATA_PATH=\'${WPS_GEOG_DIR}\'
export NL_OPT_GEOGRID_TBL_PATH=\'${WPS_DIR}/geogrid\'
#
# WPS UNGRIB NAMELIST:
export NL_OUT_FORMAT=\'WPS\'
#
# WPS METGRID NAMELIST:
export NL_IO_FORM_METGRID=2
#
#
#########################################################################
#
# CREATE RUN DIRECTORY
#
#########################################################################
#
if [[ ! -e ${RUN_DIR} ]]; then mkdir ${RUN_DIR}; fi
cd ${RUN_DIR}
#
#########################################################################
#
# RUN GEOGRID
#
#########################################################################
#
if [[ ${RUN_GEOGRID} = "true" ]] then
   mkdir -p ${RUN_DIR}/geogrid
   cd ${RUN_DIR}/geogrid
#
   cp ${WPS_DIR}/geogrid.exe ./.
   export NL_DX=${DX_CR}
   export NL_DY=${DX_CR}
#   export NL_START_DATE=${FILE_DATE}
#   export NL_END_DATE=${NEXT_FILE_DATE}
   ${HYBRID_SCRIPTS_DIR}/da_create_wps_namelist_RT.ksh
#
   if [[ -f job.ksh ]]; then rm -rf job.ksh; fi
   touch job.ksh
   RANDOM=$$
   export JOBRND=geogrid_${RANDOM}
   rm -rf *.jerr
   rm -rf *.jout
   cat << EOF >job.ksh
#!/bin/tcsh -aeux
#PBS -N ${JOBRND}
#PBS -A ${PROJ_NUMBER_NSC}
#PBS -l ${GEOGRID_TIME_LIMIT}
#PBS -q ${GEOGRID_QUEUE} 
#PBS -j ${JOBRND}.jout
#PBS -e ${JOBRND}.jerr
#PBS -M mizzi@ucar.edu
#PBS select=${GEOGRID_NODES}:ncpus=${GEOGRID_CPUS}:mpiprocs=${GEOGRID_PROCS}
#
mpiexec_mpt omplace ./geogrid.exe  > index.html 2>&1 
#
export RC=\$?     
if [[ -f SUCCESS ]]; then rm -rf SUCCESS; fi     
if [[ -f FAILED ]]; then rm -rf FAILED; fi          
if [[ \$RC = 0 ]]; then
   touch SUCCESS
else
   touch FAILED 
   exit
fi
EOF
   qsub job.ksh 
fi
exit

export CYCLE_DATE=${NEXT_DATE}
done
#
