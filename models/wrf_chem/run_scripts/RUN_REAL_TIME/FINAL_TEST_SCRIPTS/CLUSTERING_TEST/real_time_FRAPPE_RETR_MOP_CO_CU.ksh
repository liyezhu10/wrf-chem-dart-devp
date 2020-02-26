#!/bin/ksh -aux
##########################################################################
# Purpose: Set global environment variables for real_time_wrf_chem
#
#########################################################################
#
export INITIAL_DATE=2014071400
export FIRST_FILTER_DATE=2014071406
export FIRST_DART_INFLATE_DATE=2014071406
export FIRST_EMISS_INV_DATE=2014071406
#
# START CYCLE DATE-TIME:
export CYCLE_STR_DATE=2014071606
#
# END CYCLE DATE-TIME:
export CYCLE_END_DATE=2014071606
#export CYCLE_END_DATE=${CYCLE_STR_DATE}
#
export CYCLE_DATE=${CYCLE_STR_DATE}
export NL_FAC_OBS_ERROR_MOPITT=1.00
export NL_FAC_OBS_ERROR_IASI=1.00
export RETRIEVAL_TYPE_MOPITT=RAWR
export RETRIEVAL_TYPE_IASI=RAWR
#
# NOTE: the BC temporal adjustment is setup for 6-hr cycling (BCs at 3hr and 6hr).
# km
export NL_HZ_CORR_LNGTH=400.
# m
export NL_VT_CORR_LNGTH=1000.
# hrs
export NL_TM_CORR_LNGTH_IC=24.
# hrs
export NL_TM_CORR_LNGTH_BC=72.
#
export PERT_CHEM_GENER=false
export USE_LOG=false
if [[ ${USE_LOG} == true ]]; then
   export CO_MIN=NULL
   export CO_MAX=NULL
   export O3_MIN=NULL
   export O3_MAX=NULL
   export NO2_MIN=NULL
   export NO2_MAX=NULL
   export SO2_MIN=NULL
   export SO2_MAX=NULL
   export PM10_MIN=NULL
   export PM10_MAX=NULL
   export PM25_MIN=NULL
   export PM25_MAX=NULL
   export USE_LOG_CO_LOGIC=.true.
   export USE_LOG_O3_LOGIC=.true.
   export USE_LOG_NOX_LOGIC=.false.
   export USE_LOG_SO2_LOGIC=.false.
   export USE_LOG_PM10_LOGIC=.false.
   export USE_LOG_PM25_LOGIC=.false.
   export USE_LOG_AOD_LOGIC=.false.
else
   export CO_MIN=1.e-4
   export CO_MAX=NULL
   export O3_MIN=0.
   export O3_MAX=NULL
   export NO2_MIN=0.
   export NO2_MAX=NULL
   export SO2_MIN=0.
   export SO2_MAX=NULL
   export PM10_MIN=0.
   export PM10_MAX=NULL
   export PM25_MIN=0.
   export PM25_MAX=NULL
   export USE_LOG_CO_LOGIC=.false.
   export USE_LOG_O3_LOGIC=.false.
   export USE_LOG_NOX_LOGIC=.false.
   export USE_LOG_SO2_LOGIC=.false.
   export USE_LOG_PM10_LOGIC=.false.
   export USE_LOG_PM25_LOGIC=.false.
   export USE_LOG_AOD_LOGIC=.false.
fi
#
# CPSR Truncation (limit the number of CPSR modes assimilated)
   export NL_USE_CPSR_CO_TRUNC=.false.
   export NL_CPSR_CO_TRUNC_LIM=4
   export NL_USE_CPSR_O3_TRUNC=.false.
   export NL_CPSR_O3_TRUNC_LIM=4
#
# Vertical localizaton flag (0 - retrieval locations; 1 - averaging kernel locations)
   export NL_MOPITT_CO_VLOC=0
   export NL_IASI_CO_VLOC=0
   export NL_IASI_O3_VLOC=0
#
# If VARLOC = true, then INDEP_CHEM_ASIM = false
# If INDEP_CHEM_ASIM = true, then VARLOC = false
# VARLOC and INDEP_CHEM_ASIM may both be false but
# they may not both be true
#
export VARLOC=.false.
export INDEP_CHEM_ASIM=.true.
#
export ADD_EMISS=false
export EMISS_DAMP_CYCLE=1.0
export EMISS_DAMP_INTRA_CYCLE=1.0
#
BAND_ISO_VAL_CO=.09
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
export SPECIAL_FORECAST_FAC=1./2.
export SPECIAL_FORECAST_FAC=1.
export SPECIAL_FORECAST_FAC=2./3.

export SPECIAL_FORECAST_MEM[1]=21
export SPECIAL_FORECAST_MEM[2]=9
export SPECIAL_FORECAST_MEM[3]=10
export SPECIAL_FORECAST_MEM[4]=17
export SPECIAL_FORECAST_MEM[5]=18
#
# Run temporal interpolation for missing background files
# RUN_UNGRIB, RUN_METGRID, and RUN_REAL must all be false for the interpolation and for cycling
# Currently set up for 6 hr forecasts. It can handle up to 24 hr forecasts
export RUN_INTERPOLATE=false
#
# for 2014072212 and 2014072218
#export BACK_DATE=2014072206
#export FORW_DATE=2014072300
#BACK_WT=.3333
#BACK_WT=.6667
#
# for 20142900
#export BACK_DATE=2014072818
#export FORW_DATE=2014072906
# BACK_WT=.5000
#
# for 20142912
#export BACK_DATE=2014072906
#export FORW_DATE=2014072918
# BACK_WT=.5000
#
while [[ ${CYCLE_DATE} -le ${CYCLE_END_DATE} ]]; do
   export DATE=${CYCLE_DATE}
   export CYCLE_PERIOD=6
   export HISTORY_INTERVAL_HR=1
   (( HISTORY_INTERVAL_MIN = ${HISTORY_INTERVAL_HR} * 60 ))
   export HISTORY_INTERVAL_MIN=${HISTORY_INTERVAL_MIN}
   export START_IASI_O3_DATA=2014060100
   export END_IASI_O3_DATA=2014073118
   export NL_DEBUG_LEVEL=200
#
# CODE VERSIONS:
   export WPS_VER=WPSv3.9.1.1_dmpar
   export WPS_GEOG_VER=GEOG_DATA
   export WRFDA_VER=WRFDAv3.9.1.1_dmpar
   export WRF_VER=WRFv3.9.1.1_dmpar
   export WRFCHEM_VER=WRFCHEMv3.9.1.1_dmpar
   export DART_VER=DART_CHEM_REPOSITORY
#
# ROOT DIRECTORIES:
   export SCRATCH_DIR=/scratch/summit/mizzi
   export WORK_DIR=/projects/mizzi
   export INPUT_DATA_DIR=/gpfs/summit/datasets/GEOSChem_met_emis/wrf

   export EXPERIMENT_DIR=${SCRATCH_DIR}
   export RUN_DIR=${EXPERIMENT_DIR}/real_FRAPPE_RETR_MOP_CO_CLUSTER
   export TRUNK_DIR=${WORK_DIR}/TRUNK
   export DART_DIR=${TRUNK_DIR}/${DART_VER}
   export CLUSTER_DIR=${DART_DIR}/models/wrf_chem/run_scripts/RUN_REAL_TIME/FINAL_TEST_SCRIPTS/CLUSTERING_TEST
#
#########################################################################
#
# SETUP DEPENDENT RUN DIRECTORY PATHS
#
#########################################################################
#
   source ${CLUSTER_DIR}/WRFCHEM-DART_DEPENDENT_RUN_DIRS.ksh
#
#########################################################################
#
# SETUP DATE AND TIME VARIABLES
#
#########################################################################
#
   source ${CLUSTER_DIR}/WRFCHEM-DART_DATE_TIME_VARS.ksh
#
# SELECT COMPONENT RUN OPTIONS:
   if [[ ${RUN_SPECIAL_FORECAST} = "false" ]]; then
      export RUN_GEOGRID=false
      export RUN_UNGRIB=false
      export RUN_METGRID=false
      export RUN_REAL=false
      export RUN_EXO_COLDENS=false
      export RUN_SEASON_WES=false
      export RUN_PERT_WRFCHEM_MET_IC=false
      export RUN_PERT_WRFCHEM_MET_BC=false
      export RUN_WRFCHEM_BIO=false
      export RUN_WRFCHEM_FIRE=false
      export RUN_WRFCHEM_CHEMI=false
      export RUN_PERT_WRFCHEM_CHEM_ICBC=false
      export RUN_PERT_WRFCHEM_CHEM_EMISS=false
      export RUN_MOPITT_CO_OBS=false
      export RUN_IASI_CO_OBS=false
      export RUN_IASI_O3_OBS=false
      export RUN_OMI_NO2_OBS=false
      export RUN_AIRNOW_O3_OBS=false
      export RUN_AIRNOW_CO_OBS=false
      export RUN_AIRNOW_NO2_OBS=false
      export RUN_AIRNOW_SO2_OBS=false
      export RUN_AIRNOW_PM10_OBS=false
      export RUN_AIRNOW_PM25_OBS=false
      export RUN_PANDA_CO_OBS=false
      export RUN_PANDA_O3_OBS=false
      export RUN_PANDA_PM25_OBS=false
      export RUN_MODIS_AOD_OBS=false
      export RUN_MET_OBS=false 
      export RUN_COMBINE_OBS=false
      export RUN_PREPROCESS_OBS=false
#
      if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
         export RUN_WRFCHEM_INITIAL=true
         export RUN_DART_FILTER=false
         export RUN_UPDATE_BC=false
         export RUN_ENSEMBLE_MEAN_INPUT=true
         export RUN_WRFCHEM_CYCLE_CR=false
         export RUN_BAND_DEPTH=false
         export RUN_WRFCHEM_CYCLE_FR=false
         export RUN_ENSMEAN_CYCLE_FR=false
         export RUN_ENSEMBLE_MEAN_OUTPUT=true
      else
         export RUN_WRFCHEM_INITIAL=false
         export RUN_DART_FILTER=false
         export RUN_UPDATE_BC=true
         export RUN_ENSEMBLE_MEAN_INPUT=false
         export RUN_WRFCHEM_CYCLE_CR=true
         export RUN_BAND_DEPTH=false
         export RUN_WRFCHEM_CYCLE_FR=false
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
      export RUN_OMI_NO2_OBS=false
      export RUN_AIRNOW_O3_OBS=false
      export RUN_AIRNOW_CO_OBS=false
      export RUN_AIRNOW_NO2_OBS=false
      export RUN_AIRNOW_SO2_OBS=false
      export RUN_AIRNOW_PM10_OBS=false
      export RUN_AIRNOW_PM25_OBS=false
      export RUN_PANDA_CO_OBS=false
      export RUN_PANDA_O3_OBS=false
      export RUN_PANDA_PM25_OBS=false
      export RUN_MODIS_AOD_OBS=false
      export RUN_MET_OBS=false
      export RUN_COMBINE_OBS=false
      export RUN_PREPROCESS_OBS=false
#
      if [[ ${DATE} -eq ${INITIAL_DATE}  ]]; then
         export RUN_WRFCHEM_INITIAL=true
         export RUN_DART_FILTER=false
         export RUN_UPDATE_BC=false
         export RUN_ENSEMBLE_MEAN_INPUT=false
         export RUN_WRFCHEM_CYCLE_CR=false
         export RUN_BAND_DEPTH=false
         export RUN_WRFCHEM_CYCLE_FR=false
         export RUN_ENSMEAN_CYCLE_FR=false
         export RUN_ENSEMBLE_MEAN_OUTPUT=true
      else
         export RUN_WRFCHEM_INITIAL=false
         export RUN_DART_FILTER=false
         export RUN_UPDATE_BC=false
         export RUN_ENSEMBLE_MEAN_INPUT=false
         export RUN_WRFCHEM_CYCLE_CR=true
         export RUN_BAND_DEPTH=false
         export RUN_WRFCHEM_CYCLE_FR=false
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
      export RUN_OMI_NO2_OBS=false
      export RUN_AIRNOW_O3_OBS=false
      export RUN_AIRNOW_CO_OBS=false
      export RUN_AIRNOW_NO2_OBS=false
      export RUN_AIRNOW_SO2_OBS=false
      export RUN_AIRNOW_PM10_OBS=false
      export RUN_AIRNOW_PM25_OBS=false
      export RUN_PANDA_CO_OBS=false
      export RUN_PANDA_O3_OBS=false
      export RUN_PANDA_PM25_OBS=false
      export RUN_MODIS_AOD_OBS=false
      export RUN_MET_OBS=false
      export RUN_COMBINE_OBS=false
      export RUN_PREPROCESS_OBS=false
      export RUN_WRFCHEM_INITIAL=false
      export RUN_DART_FILTER=false
      export RUN_UPDATE_BC=false
      export RUN_ENSEMBLE_MEAN_INPUT=false
      export RUN_WRFCHEM_CYCLE_CR=false
      export RUN_BAND_DEPTH=false
      export RUN_WRFCHEM_CYCLE_FR=false
      export RUN_ENSMEAN_CYCLE_FR=true
      export RUN_ENSEMBLE_MEAN_OUTPUT=false
   fi
#
# FORECAST PARAMETERS:
   export USE_DART_INFL=true
   export FCST_PERIOD=6
   (( CYCLE_PERIOD_SEC=${CYCLE_PERIOD}*60*60 ))
   export CYCLE_PERIOD_SEC=${CYCLE_PERIOD_SEC}
   export NUM_MEMBERS=30
   export MAX_DOMAINS=02
   export CR_DOMAIN=01
   export FR_DOMAIN=02
   export NNXP_CR=179
   export NNYP_CR=139
   export NNZP_CR=36
   export NNXP_FR=320
   export NNYP_FR=290
   export NNZP_FR=36
   (( NNXP_STAG_CR=${NNXP_CR}+1 ))
   (( NNYP_STAG_CR=${NNYP_CR}+1 ))
   (( NNZP_STAG_CR=${NNZP_CR}+1 ))
   (( NNXP_STAG_FR=${NNXP_FR}+1 ))
   (( NNYP_STAG_FR=${NNYP_FR}+1 ))
   (( NNZP_STAG_FR=${NNZP_FR}+1 ))
   export NNXP_STAG_CR=${NNXP_STAG_CR}
   export NNYP_STAG_CR=${NNYP_STAG_CR}
   export NNZP_STAG_CR=${NNZP_STAG_CR}
   export NNXP_STAG_FR=${NNXP_STAG_FR}
   export NNYP_STAG_FR=${NNYP_STAG_FR}
   export NNZP_STAG_FR=${NNZP_STAG_FR}
   export NSPCS=61
   export NNZ_CHEM=11
   export NNCHEM_SPC=49
   export NNFIRE_SPC=31
   export NNBIO_SPC=1
   export NZ_CHEMI=${NNZ_CHEM}
   export NZ_FIRECHEMI=1
   export NCHEMI_EMISS=2
   export NFIRECHEMI_EMISS=7
   export ISTR_CR=1
   export JSTR_CR=1
   export ISTR_FR=86
   export JSTR_FR=35
   export DX_CR=15000
   export DX_FR=3000
   (( LBC_END=2*${FCST_PERIOD} ))
   export LBC_END=${LBC_END}
   export LBC_FREQ=3
   (( INTERVAL_SECONDS=${LBC_FREQ}*60*60 ))
   export INTERVAL_SECONDS=${INTERVAL_SECONDS}
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
   export ACCOUNT=ucb93_summit2
   export SINGLE_JOB_CLASS=normal
   export SINGLE_TIME_LIMIT=00:05:00
   export SINGLE_NODES=1
   export SINGLE_TASKS=1
   export WRFDA_JOB_CLASS=normal
   export WRFDA_TIME_LIMIT=00:30:00
   export WRFDA_NODES=1
   export WRFDA_TASKS=12
   export PERT_JOB_CLASS=normal
   export PERT_TIME_LIMIT=02:30:00
   export PERT_NODES=2-4
   (( PERT_TASKS=${NUM_MEMBERS}+1 ))
   export PERT_TASKS=${PERT_TASKS}
   export GROUP3_JOB_CLASS=normal
   export GROUP3_TIME_LIMIT=00:10:00
   export GROUP3_NODES=1
   export GROUP3_TASKS=3
   export GROUP15_JOB_CLASS=normal
   export GROUP15_TIME_LIMIT=00:05:00
   export GROUP15_NODES=1
   export GROUP15_TASKS=15
   export FILTER_JOB_CLASS=normal
   export FILTER_TIME_LIMIT=03:30:00
   export FILTER_NODES=2-4
   export FILTER_TASKS=48
   export WRFCHEM_JOB_CLASS=normal
   export WRFCHEM_TIME_LIMIT=01:00:00
   export WRFCHEM_NODES=2-4
   export WRFCHEM_TASKS=48
#
#########################################################################
#
# SETUP RUNTIME OUPUT DIRECTORY PATHS
#
#########################################################################
#
   source ${CLUSTER_DIR}/WRFCHEM-DART_RUNTIME_DIRS.ksh
#
# WPS PARAMETERS:
   export SINGLE_FILE=false
   export HOR_SCALE=1500
   export VTABLE_TYPE=GFS
   export METGRID_TABLE_TYPE=ARW
#
# WRF PREPROCESS PARAMETERS
# TARG_LAT=31.56 (33,15) for 072600
# TARG_LON=-120.14 = 239.85 (33,15)
#   export NL_MIN_LAT=27.5
#   export NL_MAX_LAT=38.5
#   export NL_MIN_LON=234.5
#   export NL_MAX_LON=244.5
#
# NL_MIN_LON, NL_MAX_LON = [-180.,190.]
# NL_MIN_LAT, NL_MAX_LAT = [-90.,90.]
# NNL_MIN_LON, NL_MAX_LON = [0.,360.]
# NNL_MIN_LON, NL_MAX_LON = [-90.,90.]
#
   export NL_MIN_LAT=27
   export NL_MAX_LAT=48
   export NL_MIN_LON=-132
   export NL_MAX_LON=-94
#
   export NNL_MIN_LON=${NL_MIN_LON}
   if [[ ${NL_MIN_LON} -lt 0. ]]; then
      (( NNL_MIN_LON=${NL_MIN_LON}+360 ))
   fi
   export NNL_MAX_LON=${NL_MAX_LON}
   if [[ ${NL_MAX_LON} -lt 0. ]]; then
      (( NNL_MAX_LON=${NL_MAX_LON}+360 ))
   fi
   export NNL_MIN_LAT=${NL_MIN_LAT}
   export NNL_MAX_LAT=${NL_MAX_LAT}
   export NL_OBS_PRESSURE_TOP=1000.
#
# PERT CHEM PARAMETERS
   export SPREAD_FAC=0.30
   export NL_SPREAD_CHEMI=${SPREAD_FAC}
   export NL_SPREAD_FIRE=0.00
   export NL_SPREAD_BIOG=0.00
   export NL_PERT_CHEM=true
   export NL_PERT_FIRE=false
   export NL_PERT_BIO=false
#
#########################################################################
#
# SETUP WPS; WRFDA; WRF; WRFCHEM; DART NAMELISTS
#
#########################################################################
#
   source ${CLUSTER_DIR}/WRFCHEM-DART_NAMELISTS.ksh
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
# RUN WPS CLUSTER: GEOGRID; UNBRIB; METGRID; REAL
# RUN WRFCHEM STATIC FILES: EXO_COLDENS; SEASONS_WES  
#
#########################################################################
#
   if ${RUN_GEOGRID} || ${RUN_UNGRIB} || ${RUN_METGRID} || ${RUN_REAL} || ${RUN_EXO_COLDENS} || ${RUN_SEASON_WES}; then
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh wps ${SINGLE_JOB_CLASS} ${SINGLE_TIME_LIMIT} ${SINGLE_NODES} ${SINGLE_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_WPS.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# RUN INTERPOLATION TO GET MISSING BACKGROUND DATA
#
#########################################################################
#
   if ${RUN_INTERPOLATE}; then
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh intrp ${SINGLE_JOB_CLASS} ${SINGLE_TIME_LIMIT} ${SINGLE_NODES} ${SINGLE_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_INTERPOLATE.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# RUN ENSEMBLE GENERATOR METEOROLOGY ICs/BCs: PERT_WRFCHEM_MET_IC; PERT_WRFCHEM_MET_BC
#
#########################################################################
#
   if ${RUN_PERT_WRFCHEM_MET_IC} || ${RUN_PERT_WRFCHEM_MET_BC}; then 
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh pert_met ${WRFDA_JOB_CLASS} ${WRFDA_TIME_LIMIT} ${WRFDA_NODES} ${WRFDA_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PERT_MET_ICBC.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# GENERATE CHEMISTRY ICs/BCs EMISSIONS: 
# WRFCHEM_BIO, WRFCHEM_FIRE; WRFCHEM_ANTHRO
#
#########################################################################
#
   if ${RUN_WRFCHEM_BIO} || ${RUN_WRFCHEM_FIRE} || ${RUN_WRFCHEM_CHEMI}; then 
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh emiss ${GROUP3_JOB_CLASS} ${GROUP3_TIME_LIMIT} ${GROUP3_NODES} ${GROUP3_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_EMISSIONS.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.file
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# RUN ENSEMBLE GENERATOR CHEMISTRY ICs/BCs AND EMISSIONS: 
# PERT_WRFCHEM_CHEM_ICBC; PERT_WRFCHEM_EMISS
#
#########################################################################
#
   if ${RUN_PERT_WRFCHEM_CHEM_ICBC} || ${RUN_PERT_WRFCHEM_CHEM_EMISS}; then 
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh pert_chem ${PERT_JOB_CLASS} ${PERT_TIME_LIMIT} ${PERT_NODES} ${PERT_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_PERT_CHEM_ICBC_EMISS.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# GENERATE OBSERVATIONS FOR DART: 
# MOPITT_CO; IASI_CO; IASI_O3; OMI_NO2; AIRNOW_CO; AIRNOW_O3; AIRNOW_NO2;
# AIRNOW_SO2; AIRNOW_PM10; AIRNOW_PM25; PANDA_CO; PANDA_O3; PANDA_PM25;
# MODIS_AOD; MET_OBS
#
#########################################################################
#
   if ${RUN_MOPITT_CO_OBS} || ${RUN_IASI_CO_OBS} || ${RUN_IASI_O3_OBS} || ${RUN_OMI_NO2_OBS} || ${RUN_AIRNOW_CO_OBS} || ${RUN_AIRNOW_O3_OBS} || ${RUN_AIRNOW_NO2_OBS} || ${RUN_AIRNOW_SO2_OBS} || ${RUN_AIRNOW_PM10_OBS} || ${RUN_AIRNOW_PM25_OBS} || ${RUN_PANDA_CO_OBS} || ${RUN_PANDA_O3_OBS} || ${RUN_PANDA_PM25_OBS} || ${RUN_MODIS_AOD_OBS} || ${RUN_MET_OBS}; then
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh emiss ${GROUP15_JOB_CLASS} ${GROUP15_TIME_LIMIT} ${GROUP15_NODES} ${GROUP15_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_OBSERVATIONS.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.file
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# COMBINE AND PREPROCESS THE OBSERVATIONS: 
#
#########################################################################
#
   if ${RUN_COMBINE_OBS} || ${RUN_PREPROCESS_OBS}; then
      ${HYBRID_SCRIPTS_DIR}/job_header_summit.ksh pert_chem ${SINGLE_JOB_CLASS} ${SINGLE_TIME_LIMIT} ${SINGLE_NODES} ${SINGLE_TASKS} ${ACCOUNT} ${CLUSTER_DIR}/WRFCHEM-DART_RUN_COMBINE_OBS.ksh
      sbatch -W job.ksh
      export RC=$?     
      rm -rf job.ksh
      rm -rf SUCCESS     
      rm -rf FAILED          
      if [[ $RC = 0 ]]; then
         touch SUCCESS
      else
         touch FAILED 
         exit
      fi
   fi
#
#########################################################################
#
# RUN DART_FILTER AND UPDATE_BC: 
#
#########################################################################
#
   if ${RUN_DART_FILTER} || ${RUN_UPDATE_BC}; then
      ${CLUSTER_DIR}/WRFCHEM-DART_RUN_FILTER_UPDATE_BC.ksh
   fi
#
#########################################################################
#
# RUN INITIAL CYCLE: 
#
#########################################################################
#
   if ${RUN_WRFCHEM_INITIAL}; then
      ${CLUSTER_DIR}/WRFCHEM-DART_RUN_INITIAL_CYCLE.ksh
   fi
#
#########################################################################
#
# RUN ASSIMILATION CYCLE: 
#
#########################################################################
#
   if ${RUN_WRFCHEM_CYCLE_CR}; then
      ${CLUSTER_DIR}/WRFCHEM-DART_RUN_ASSIM_CYCLE.ksh
   fi
#
#########################################################################
#
# RUN FINE GRID FORECAST: 
#
#########################################################################
#
   if ${RUN_BAND_DEPTH} || ${RUN_WRFCHEM_CYCLE_FR} || ${RUN_ENSMEAN_CYCLE_FR}; then
      ${CLUSTER_DIR}/WRFCHEM-DART_RUN_FINE_GRID_FORE.ksh
   fi
#
#########################################################################
#
# RUN ENSEMBLE STATISTICS: 
#
#########################################################################
#
   if ${RUN_ENSEMBLE_MEAN_INPUT} || ${RUN_ENSEMBLE_MEAN_OUTPUT}; then
      ${CLUSTER_DIR}/WRFCHEM-DART_RUN_ENSEMBLE_MEAN_STATS.ksh
   fi
#
#########################################################################
#
# GO TO NEXT CYCLE: 
#
#########################################################################
#
   export CYCLE_DATE=${NEXT_DATE}
done
#
