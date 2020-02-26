#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: job_script_summit.ksh 13133 2019-04-25 21:47:54Z nancy@ucar.edu $
#
#########################################################################
#
# RUN WRFCHEM PERTURB ICBC
#
#########################################################################
#
   if ${RUN_PERT_WRFCHEM_CHEM_ICBC}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chem_icbc ]]; then
         mkdir ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
         cd ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_chem_icbc
      fi
#
# PERTURB CHEM ICBC
      export NL_SW_CORR_TM=false
      if [[ ${DATE} -eq ${INITIAL_DATE} ]]; then export NL_SW_CORR_TM=true; fi
      export NL_SW_SEED=true
#
      cp ${METGRID_DIR}/met_em.d${CR_DOMAIN}.*:00:00.nc ./.
#      cp ${METGRID_DIR}/met_em.d${FR_DOMAIN}.*:00:00.nc ./.
      cp ${PERT_CHEM_INPUT_DIR}/work/perturb_chem_icbc_CORR_RT_MA.exe ./.
      cp ${PERT_CHEM_INPUT_DIR}/work/perturb_chem_icbc_CORR_RT_MA_MPI.exe ./.
      cp ${PERT_CHEM_INPUT_DIR}/work/mozbc.exe ./mozbc.exe
      cp ${PERT_CHEM_INPUT_DIR}/runICBC_parent_rt_CR.ksh ./.
      cp ${PERT_CHEM_INPUT_DIR}/runICBC_parent_rt_FR.ksh ./.
      cp ${PERT_CHEM_INPUT_DIR}/run_mozbc_rt_CR.csh ./.
      cp ${PERT_CHEM_INPUT_DIR}/run_mozbc_rt_FR.csh ./.
      cp ${PERT_CHEM_INPUT_DIR}/set00 ./.
#
# SELECT MOZART DATA FILE
#      if [[ ${YYYY} -eq 2014 ]]; then export MOZBC_DATA=/h0003.nc; fi
      if [[ ${YYYY} -eq 2014 ]]; then export MOZBC_DATA=/h0004.nc; fi
#
# CREATE INPUT FILES COARSE DOMAIN
      rm -rf mozbc.ic.inp
      cat << EOF > mozbc.ic.inp
&control
do_bc     = .false.
do_ic     = .true.
domain    = 1
dir_wrf   = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz   = '${MOZBC_DATA_DIR}'
fn_moz    = '${MOZBC_DATA}'
def_missing_var    = .true.
met_file_prefix    = 'met_em'
met_file_suffix    = '.nc'
met_file_separator = '.'
EOF
      rm -rf mozbc.bc.inp
      cat << EOF > mozbc.bc.inp
&control
do_bc     = .true.
do_ic     = .false.
domain    = 1
dir_wrf   = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
dir_moz   = '${MOZBC_DATA_DIR}'
fn_moz    = '${MOZBC_DATA}'
def_missing_var    = .true.
met_file_prefix    = 'met_em'
met_file_suffix    = '.nc'
met_file_separator = '.'
EOF
#
      ./runICBC_parent_rt_CR.ksh
#
# CREATE INPUT FILES FINE DOMAIN
#      rm -rf mozbc.ic.inp
#      cat << EOF > mozbc.ic.inp
#&control
#do_bc     = .false.
#do_ic     = .true.
#domain    = 2
#dir_wrf   = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
#dir_moz   = '${MOZBC_DATA_DIR}'
#fn_moz    = '${MOZBC_DATA}'
#def_missing_var    = .true.
#met_file_prefix    = 'met_em'
#met_file_suffix    = '.nc'
#met_file_separator = '.'
#EOF
#      rm -rf mozbc.bc.inp
#      cat << EOF > mozbc.bc.inp
#&control
#do_bc     = .true.
#do_ic     = .false.
#domain    = 2
#dir_wrf   = '${RUN_DIR}/${DATE}/wrfchem_chem_icbc/'
#dir_moz   = '${MOZBC_DATA_DIR}'
#fn_moz    = '${MOZBC_DATA}'
#def_missing_var    = .true.
#met_file_prefix    = 'met_em'
#met_file_suffix    = '.nc'
#met_file_separator = '.'
#EOF
#
#      ./runICBC_parent_rt_FR.ksh
#
# GENERATE CHEMISTRY IC/BC ENSEMBLE MEMBERS
#
# CREATE NAMELIST
      export WRFINPEN=wrfinput_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
      export WRFBDYEN=wrfbdy_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
      export WRFINPUT_FLD_RW=wrfinput_d${CR_DOMAIN}
      export WRFINPUT_ERR_RW=wrfinput_d${CR_DOMAIN}_err
      export WRFBDY_FLD_RW=wrfbdy_d${CR_DOMAIN}
      mv ${WRFINPEN} ${WRFINPUT_FLD_RW}
      mv ${WRFBDYEN} ${WRFBDY_FLD_RW}
      rm -rf perturb_chem_icbc_corr_nml.nl
      cat << EOF > perturb_chem_icbc_corr_nml.nl
&perturb_chem_icbc_corr_nml
nx=${NNXP_CR},
ny=${NNYP_CR},
nz=${NNZP_CR},
nchem_spcs=${NSPCS},
pert_path_old='${RUN_DIR}/${PAST_DATE}/wrfchem_chem_icbc',
pert_path_new='${RUN_DIR}/${DATE}/wrfchem_chem_icbc',
nnum_mem=${NUM_MEMBERS},
wrfinput_fld_new='${WRFINPUT_FLD_RW}',
wrfinput_err_new='${WRFINPUT_ERR_RW}',
wrfbdy_fld_new='${WRFBDY_FLD_RW}',
sprd_chem=${SPREAD_FAC},
corr_lngth_hz=${NL_HZ_CORR_LNGTH},
corr_lngth_vt=${NL_VT_CORR_LNGTH},
corr_lngth_tm=${NL_TM_CORR_LNGTH_IC},
corr_tm_delt=${CYCLE_PERIOD},
sw_corr_tm=${NL_SW_CORR_TM},
sw_seed=${NL_SW_SEED},
/
EOF
      rm -rf perturb_chem_icbc_spcs_nml.nl
      cat << EOF > perturb_chem_icbc_spcs_nml.nl
&perturb_chem_icbc_spcs_nml
ch_chem_spc='o3','no','no2','no3','nh3','hno3','hno4','n2o5','ho2','h2o2','ch4','co','ch3o2','ch3ooh','hcho','ch3oh','c2h4','ald','ch3cooh','acet','mgly','pan','mpan','macr','mvk','c2h6','c3h6','c3h8','c2h5oh','c10h16','onit','onitr','isopr','isopn','acetol','glyald','hydrald','mek','bigene','open','bigalk','tol','cres','dms','so2','sulf','BC1','BC2','OC1','OC2','SEAS_1','SEAS_2','SEAS_3','SEAS_4','DUST_1','DUST_2','DUST_3','DUST_4','DUST_5','h2','n2o'
/
EOF
#
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         cp ${WRFINPUT_FLD_RW} ${WRFINPUT_FLD_RW}.${CMEM}   
#         cp ${WRFINPUT_FLD_RW} ${WRFINPUT_ERR_RW}.${CMEM}   
         cp ${WRFBDY_FLD_RW} ${WRFBDY_FLD_RW}.${CMEM}   
         let MEM=MEM+1
      done
      cp ${WRFINPUT_FLD_RW} ${WRFINPUT_FLD_RW}_mean
      cp ${WRFINPUT_FLD_RW} ${WRFINPUT_FLD_RW}_sprd
      cp ${WRFINPUT_FLD_RW} ${WRFINPUT_FLD_RW}_frac
#
      export TYPE=PARALLEL
      export EXECUTE=perturb_chem_icbc_CORR_RT_MA_MPI.exe
      if [[ ${TYPE} == PARALLEL ]]; then
         mpirun -np ${PERT_TASKS} ./${EXECUTE} > index.html 2>&1
      elif [[ ${TYPE} == SERIAL ]]; then
         ./${EXECUTE} > index.html 2>&1
      fi
#
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
         mv ${WRFINPUT_FLD_RW}.${CMEM} ${WRFINPEN}.${CMEM}
         mv ${WRFBDY_FLD_RW}.${CMEM} ${WRFBDYEN}.${CMEM}
         let MEM=MEM+1
      done
#
# COMBINE WRFCHEM WITH WRF CR PARENT FILES
      ncks -A ${REAL_DIR}/${WRFINPEN} ${WRFINPEN}
      ncks -A ${EXPERIMENT_DUST_DIR}/EROD_d${CR_DOMAIN} ${WRFINPEN}
      ncks -A ${REAL_DIR}/${WRFBDYEN} ${WRFBDYEN}
#
# COMBINE WRFCHEM WITH WRF FR DOMAIN PARENT FILES
#      export WRFINPEN=wrfinput_d${FR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
#      ncks -A ${REAL_DIR}/${WRFINPEN} ${WRFINPEN}
#      ncks -A ${EXPERIMENT_DUST_DIR}/EROD_d${FR_DOMAIN} ${WRFINPEN}
#
# LOOP THROUGH ALL MEMBERS IN THE ENSEMBLE
      let MEM=1
      while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
         export CMEM=e${MEM}
         if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
         if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# COMBINE WRFCHEM WITH WRF CR DOMAIN
#
         export WRFINPEN=wrfinput_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
         export WRFBDYEN=wrfbdy_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
         ncks -A ${WRFCHEM_MET_IC_DIR}/${WRFINPEN} ${WRFINPEN}
         ncks -A ${EXPERIMENT_DUST_DIR}/EROD_d${CR_DOMAIN} ${WRFINPEN}
         ncks -A ${WRFCHEM_MET_BC_DIR}/${WRFBDYEN} ${WRFBDYEN}
#
# COMBINE WRFCHEM WITH WRF FR DOMAIN
#         export WRFINPEN=wrfinput_d${FR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
#         ncks -A ${WRFCHEM_MET_IC_DIR}/${WRFINPEN} ${WRFINPEN}
#         ncks -A ${EXPERIMENT_DUST_DIR}/EROD_d${FR_DOMAIN} ${WRFINPEN}
#
         let MEM=MEM+1
      done
   fi
#
#########################################################################
#
# RUN WRFCHEM PERTURB EMISSIONS
#
#########################################################################
#
   if ${RUN_PERT_WRFCHEM_CHEM_EMISS}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_chem_emiss ]]; then
         mkdir ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
         cd ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_chem_emiss
      fi
#
# SET PARAMETERS
      export NL_EMISS_TIME=0.0
      export NL_SW_SEED=true
      export NL_SW_CORR_TM=false
      if [[ ${DATE} -eq ${INITIAL_DATE} ]]; then export NL_SW_CORR_TM=true; fi
#
# COPY PERTURBATION CODE
      rm -rf perturb_chem_emiss_CORR_RT_MA.exe
      rm -rf perturb_chem_emiss_CORR_RT_MA_MPI.exe
      cp ${PERT_CHEM_EMISS_DIR}/work/perturb_chem_emiss_CORR_RT_MA.exe ./.
      cp ${PERT_CHEM_EMISS_DIR}/work/perturb_chem_emiss_CORR_RT_MA_MPI.exe ./.
#
      export L_DATE=${DATE}
      export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
      while [[ ${L_DATE} -le ${LE_DATE} ]] ; do
         export L_YYYY=$(echo $L_DATE | cut -c1-4)
         export L_MM=$(echo $L_DATE | cut -c5-6)
         export L_DD=$(echo $L_DATE | cut -c7-8)
         export L_HH=$(echo $L_DATE | cut -c9-10)
#
         export NL_PERT_PATH_PR=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
         export NL_PERT_PATH_PO=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
         if [[ ${L_DATE} -eq ${DATE} || ${L_HH} -eq 00 ]]; then
            export NL_PERT_PATH_PR=${RUN_DIR}/${PAST_DATE}/wrfchem_chem_emiss
            export NL_PERT_PATH_PO=${RUN_DIR}/${DATE}/wrfchem_chem_emiss
         fi
#
# GET COARSE GRID EMISSON FILES
         export WRFCHEMI=wrfchemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFFIRECHEMI=wrffirechemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export WRFBIOCHEMI=wrfbiochemi_d${CR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
# COPY DATA

         cp ${WRFCHEM_CHEMI_DIR}/${WRFCHEMI} ${WRFCHEMI}
         chmod a+rwx ${WRFCHEMI}
         ncatted -O -a coordinates,E_CO,c,c,"XLONG, XLAT" ${WRFCHEMI}
         ncatted -O -a coordinates,E_NO,c,c,"XLONG, XLAT" ${WRFCHEMI}
         cp ${WRFCHEM_FIRE_DIR}/${WRFFIRECHEMI} ${WRFFIRECHEMI}
         chmod a+rwx ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_co,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_no,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_oc,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_bc,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_c2h4,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_ch2o,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         ncatted -O -a coordinates,ebu_in_ch3oh,c,c,"XLONG, XLAT" ${WRFFIRECHEMI}
         if [[ ${L_HH} -eq 00 || ${L_HH} -eq 06 || ${L_HH} -eq 12 || ${L_HH} -eq 18 ]]; then
            cp ${WRFCHEM_BIO_DIR}/${WRFBIOCHEMI} ${WRFBIOCHEMI}
            chmod a+rwx ${WRFBIOCHEMI}
         fi
         let MEM=1
         while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
            export CMEM=e${MEM}
            if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
            if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#
# link ensemble members 
            export WRFINPUT=wrfinput_d${CR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00.${CMEM}
            rm -rf wrfinput_d${CR_DOMAIN}.${CMEM}
            ln -sf ${WRFCHEM_CHEM_ICBC_DIR}/${WRFINPUT} wrfinput_d${CR_DOMAIN}.${CMEM}
#
# wrfchemi
            rm -rf ${WRFCHEMI}.${CMEM}
            cp ${WRFCHEMI} ${WRFCHEMI}.${CMEM}
#
# wrffire
            rm -rf ${WRFFIRECHEMI}.${CMEM}
            cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
#
# wrfbio
            if [[ ${L_HH} -eq 00 || ${L_HH} -eq 06 || ${L_HH} -eq 12 || ${L_HH} -eq 18 ]]; then
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
            fi
            let MEM=MEM+1
         done
         if [[ ${NL_PERT_CHEM} == true ]]; then 
            cp ${WRFCHEMI} ${WRFCHEMI}_mean
            cp ${WRFCHEMI} ${WRFCHEMI}_sprd
            cp ${WRFCHEMI} ${WRFCHEMI}_frac
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then 
            cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}_mean
            cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}_sprd
            cp ${WRFFIRECHEMI} ${WRFFIRECHEMI}_frac
         fi
         if [[ ${NL_PERT_FIRE} == true ]]; then 
            if [[ ${L_HH} -eq 00 || ${L_HH} -eq 06 || ${L_HH} -eq 12 || ${L_HH} -eq 18 ]]; then
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}_mean
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}_sprd
               cp ${WRFBIOCHEMI} ${WRFBIOCHEMI}_frac
            fi
         fi
#
# CREATE NAMELIST
         rm -rf perturb_chem_emiss_corr_nml.nl
         cat << EOF > perturb_chem_emiss_corr_nml.nl
&perturb_chem_emiss_CORR_nml
nx=${NNXP_CR},
ny=${NNYP_CR},
nz=${NNZP_CR},
nz_chem=${NNZ_CHEM},
nchem_spcs=${NNCHEM_SPC},
nfire_spcs=${NNFIRE_SPC},
nbiog_spcs=${NNBIO_SPC},
pert_path_pr='${NL_PERT_PATH_PR}',
pert_path_po='${NL_PERT_PATH_PO}',
nnum_mem=${NUM_MEMBERS},
wrfchemi='${WRFCHEMI}',
wrffirechemi='${WRFFIRECHEMI}',
wrfbiogchemi='${WRFBIOCHEMI}',
sprd_chem=${NL_SPREAD_CHEMI},
sprd_fire=${NL_SPREAD_FIRE},
sprd_biog=${NL_SPREAD_BIOG},
sw_corr_tm=${NL_SW_CORR_TM},
sw_seed=${NL_SW_SEED},
sw_chem=${NL_PERT_CHEM},
sw_fire=${NL_PERT_FIRE},
sw_biog=${NL_PERT_BIO},
corr_lngth_hz=${NL_HZ_CORR_LNGTH},
corr_lngth_vt=${NL_VT_CORR_LNGTH},
corr_lngth_tm=${NL_TM_CORR_LNGTH_BC},
corr_tm_delt=${NL_EMISS_TIME},
/
EOF
         rm -rf perturb_emiss_chem_spec_nml.nl
         cat << EOF > perturb_emiss_chem_spec_nml.nl
&perturb_chem_emiss_spec_nml
ch_chem_spc='E_CO','E_NO','E_NO2','E_BIGALK','E_BIGENE','E_C2H4','E_C2H5OH','E_C2H6','E_C3H6','E_C3H8','E_CH2O','E_CH3CHO','E_CH3COCH3','E_CH3OH','E_MEK','E_SO2','E_TOLUENE','E_NH3','E_ISOP','E_C10H16','E_sulf','E_CO_A','E_CO_BB','E_CO02','E_CO03','E_XNO','E_XNO2','E_BALD','E_C2H2','E_BENZENE','E_XYLENE','E_CRES','E_HONO','E_PM25I','E_PM25J','E_PM_10','E_ECI','E_ECJ','E_ORGI','E_ORGJ','E_SO4I','E_SO4J','E_NO3I','E_NO3J','E_NH4I','E_NH4J','E_PM_25','E_OC','E_BC',
ch_fire_spc='ebu_in_co','ebu_in_no','ebu_in_so2','ebu_in_bigalk','ebu_in_bigene','ebu_in_c2h4','ebu_in_c2h5oh','ebu_in_c2h6','ebu_in_c3h8','ebu_in_c3h6','ebu_in_ch2o','ebu_in_ch3cho','ebu_in_ch3coch3','ebu_in_ch3oh','ebu_in_mek','ebu_in_toluene','ebu_in_nh3','ebu_in_no2','ebu_in_open','ebu_in_c10h16','ebu_in_ch3cooh','ebu_in_cres','ebu_in_glyald','ebu_in_mgly','ebu_in_gly','ebu_in_acetol','ebu_in_isop','ebu_in_macr','ebu_in_mvk','ebu_in_oc','ebu_in_bc',
ch_biog_spc='MSEBIO_ISOP',
/
EOF
#
         export TYPE=PARALLEL
         export EXECUTE=perturb_chem_emiss_CORR_RT_MA_MPI.exe
         if [[ ${TYPE} == PARALLEL ]]; then
            mpirun -np ${PERT_TASKS} ./${EXECUTE} > index.html 2>&1
         elif [[ ${TYPE} == SERIAL ]]; then
            ./${EXECUTE} > index.html 2>&1
         fi
#
# ADD ENSEMBLE MEAN EMISSIONS ADJUSTMENT FROM PREVIOUS CYCLE
         if ${ADD_EMISS}; then
            rm -rf wrfchemi_d${CR_DOMAIN}_mean_incr
            rm -rf wrfchemi_d${CR_DOMAIN}_sprd_incr
            rm -rf wrffirechemi_d${CR_DOMAIN}_mean_incr
            rm -rf wrffirechemi_d${CR_DOMAIN}_sprd_incr
            rm -rf wrfbiochemi_d${CR_DOMAIN}_mean_incr
            rm -rf wrfbiochemi_d${CR_DOMAIN}_sprd_incr
#	    
            if [[ ${DATE} -gt ${FIRST_FILTER_DATE} ]]; then
               if ${NL_PERT_CHEM}; then
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrfchemi_d${CR_DOMAIN}_mean_incr ./         
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrfchemi_d${CR_DOMAIN}_sprd_incr ./         
               fi
               if ${NL_PERT_FIRE}; then
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrffirechemi_d${CR_DOMAIN}_mean_incr ./         
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrffirechemi_d${CR_DOMAIN}_sprd_incr ./         
               fi
               if ${NL_PERT_BIO}; then
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrfbiochemi_d${CR_DOMAIN}_mean_incr ./         
                  cp ${RUN_DIR}/${PAST_DATE}/dart_filter/wrfbiochemi_d${CR_DOMAIN}_sprd_incr ./         
               fi
               cp ${PERT_CHEM_EMISS_DIR}/work/perturb_chem_emiss_ADD_PRIOR_INCRs.exe ./.
#
               export TYPE=SERIAL
               export EXECUTE=perturb_chem_emiss_ADD_PRIOR_INCRs.exe
               if [[ ${TYPE} == PARALLEL ]]; then
                  mpirun -np ${PERT_TASKS} ./${EXECUTE} > index.html 2>&1
               elif [[ ${TYPE} == SERIAL ]]; then
                  ./${EXECUTE} > index.html 2>&1
               fi
            fi
         fi
#
# GET FINE GRID EMISSON FILES FOR THIS MEMBER
#         export WRFCHEMI=wrfchemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#         export WRFFIRECHEMI=wrffirechemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#         export WRFBIOCHEMI=wrfbiochemi_d${FR_DOMAIN}_${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#         export WRFINPUT=wrfinput_d${FR_DOMAIN}_${YYYY}-${MM}-${DD}_${HH}:00:00
#         export WRFINPUT_DIR=${WRFCHEM_CHEM_ICBC_DIR}
#         cp ${WRFINPUT_DIR}/${WRFINPUT} wrfinput_d${FR_DOMAIN}
#   
#         let MEM=1
#         while [[ ${MEM} -le ${NUM_MEMBERS} ]]; do
#            export CMEM=e${MEM}
#            if [[ ${MEM} -lt 100 ]]; then export CMEM=e0${MEM}; fi
#            if [[ ${MEM} -lt 10  ]]; then export CMEM=e00${MEM}; fi
#            if [[ ${NL_PERT_CHEM} == true ]]; then
#               cp ${WRFCHEM_CHEMI_DIR}/${WRFCHEMI} ${WRFCHEMI}.${CMEM}
#            fi
#            if [[ ${NL_PERT_FIRE} == true ]]; then
#               cp ${WRFCHEM_FIRE_DIR}/${WRFFIRECHEMI} ${WRFFIRECHEMI}.${CMEM}
#            fi
#            if [[ ${NL_PERT_BIO} == true ]]; then
#               cp ${WRFCHEM_BIO_DIR}/${WRFBIOCHEMI} ${WRFBIOCHEMI}.${CMEM}
#            fi
#            let MEM=MEM+1
#         done
#
# CREATE NAMELIST
#         rm -rf perturb_chem_emiss_CORR_nml.nl
#         cat << EOF > perturb_chem_emiss_CORR_nml.nl
#&perturb_chem_emiss_CORR_nml
#nx=${NNXP_FR},
#ny=${NNYP_FR},
#nz=${NNZP_FR},
#nz_chem=${NNZ_CHEM},
#nchem_spc=${NNCHEM_SPC},
#nfire_spc=${NNFIRE_SPC},
#nbio_spc=${NNBIO_SPC},
#pert_path='${RUN_DIR}',
#nnum_mem=${NUM_MEMBERS},
#wrfchemi='${WRFCHEMI}',
#wrffirechemi='${WRFFIRECHEMI}',
#wrfbiochemi='${WRFBIOCHEMI}',
#sprd_chem=${NL_SPREAD_CHEMI},
#sprd_fire=${NL_SPREAD_FIRE},
#sprd_biog=${NL_SPREAD_BIOG},
#sw_chem=${NL_PERT_CHEM},
#sw_fire=${NL_PERT_FIRE},
#sw_biog=${NL_PERT_BIO},
#/
#EOF
#         rm -rf perturb_emiss_chem_spec_nml.nl
#         cat << EOF > perturb_emiss_chem_spec_nml.nl
#&perturb_chem_emiss_spec_nml
#ch_chem_spc='E_CO','E_NO','E_NO2','E_BIGALK','E_BIGENE','E_C2H4','E_C2H5OH','E_C2H6','E_C3H6','E_C3H8','E_CH2O','E_CH3CHO','E_CH3COCH3','E_CH3OH','E_MEK','E_SO2','E_TOLUENE','E_NH3','E_ISOP','E_C10H16','E_sulf','E_CO_A','E_CO_BB','E_COO2','E_COO3','E_XNO','E_XNO2','E_PM25I','E_PM25J','E_PM_10','E_ECI','E_ECJ','E_ORGI',E_ORGJ','E_SO4I','E_SO4J','E_NO3I','E_NO3J','E_NH4I','E_NH4J','E_PM_25','E_OC','E_BC','E_BALD','E_C2H2','E_BENZENE','E_XYLENE','E_CRES','E_HONO',
#ch_fire_spc='ebu_in_co','ebu_in_no','ebu_in_so2','ebu_in_bigalk','ebu_in_bigene','ebu_in_c2h4','ebu_in_c2h5oh','ebu_in_c2h6','ebu_in_c3h8','ebu_in_c3h6','ebu_in_ch2o','ebu_in_ch3cho','ebu_in_ch3coch3','ebu_in_ch3oh','ebu_in_mek','ebu_in_toluene','ebu_in_nh3','ebu_in_no2','ebu_in_open','ebu_in_c10h16','ebu_in_ch3cooh','ebu_in_cres','ebu_in_glyald','ebu_in_mgly','ebu_in_gly','ebu_in_acetol','ebu_in_isop','ebu_in_macr','ebu_in_mvk','ebu_in_oc','ebu_in_bc',
#ch_bio_spc='MSEBIO_ISOP',
#/
#EOF
#
#         RANDOM=$$
#         export JOBRND=${RANDOM}_fr_emiss_pert
#         ${HYBRID_SCRIPTS_DIR}/job_script_summit.ksh ${JOBRND} ${GENERAL_JOB_CLASS} ${GENERAL_TIME_LIMIT} ${GENERAL_NODES} ${GENERAL_TASKS} perturb_chem_emiss_CORR_RT_CONST.exe SERIAL
#         sbatch -W job.ksh
#
# SWAP FILE NAMES FOR TEMPORAL DECORRELATION
#          rm -rf pert_chem_emiss_pr
#          mv pert_chem_emiss_po pert_chem_emiss_pr
#
# ADVANCE TIME
         (( NL_EMISS_TIME=${NL_EMISS_TIME} + 1 ))
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
      done
#
#      rm -rf wrfchemi_d${CR_DOMAIN}_tmp*
#      ncecat -n ${NUM_MEMBERS},3,1 wrfchemi_d${CR_DOMAIN}_${FILE_DATE}.e001 wrfchemi_d${CR_DOMAIN}_tmp1
#      ncwa -a record wrfchemi_d${CR_DOMAIN}_tmp1 wrfchemi_d${CR_DOMAIN}_mean
#      ncbo --op_typ='-' wrfchemi_d${CR_DOMAIN}_tmp1 wrfchemi_d${CR_DOMAIN}_mean wrfchemi_d${CR_DOMAIN}_tmp3
#      ncra -y rmssdn wrfchemi_d${CR_DOMAIN}_tmp3 wrfchemi_d${CR_DOMAIN}_sprd
#      rm -rf wrfchemi_d${CR_DOMAIN}_tmp*
   fi
