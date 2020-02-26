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
# RUN WRFCHEM_FIRE
#
#########################################################################
#
   if ${RUN_WRFCHEM_FIRE}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_fire ]]; then
         mkdir ${RUN_DIR}/${DATE}/wrfchem_fire
         cd ${RUN_DIR}/${DATE}/wrfchem_fire
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_fire
      fi
#
# LINK NEEDED FILES
      export FILE_CR=wrfinput_d${CR_DOMAIN}
      export FILE_FR=wrfinput_d${FR_DOMAIN}
      rm -rf ${FILE_CR}
      rm -rf ${FILE_FR}
      ln -sf ${REAL_DIR}/${FILE_CR}_${FILE_DATE} ${FILE_CR}   
      ln -sf ${REAL_DIR}/${FILE_FR}_${FILE_DATE} ${FILE_FR}   
      rm -rf GLOBAL*.txt
      ln -sf ${EXPERIMENT_WRFFIRECHEMI_DIR}/GLOBAL*.txt ./.
      export FILE=fire_emis.exe
      rm -rf ${FILE}
      ln -sf ${FINN_FIRE_DIR}/work/${FILE} ${FILE}
      rm -rf grass_from_img.nc
      rm -rf shrub_from_img.nc
      rm -rf tempfor_from_img.nc
      rm -rf tropfor_from_img.nc
      ln -sf ${EXPERIMENT_WRFFIRECHEMI_DIR}/grass_from_img.nc
      ln -sf ${EXPERIMENT_WRFFIRECHEMI_DIR}/shrub_from_img.nc
      ln -sf ${EXPERIMENT_WRFFIRECHEMI_DIR}/tempfor_from_img.nc
      ln -sf ${EXPERIMENT_WRFFIRECHEMI_DIR}/tropfor_from_img.nc
#
# CREATE INPUT FILE
      export FILE=fire_emis.mozc.inp
      rm -rf ${FILE}
      cat << EOF > ${FILE}
&control
domains = 2,
fire_filename(1) = 'GLOBAL_FINNv15_JULSEP2014_MOZ4_09222014.txt',
start_date = '${FIRE_START_DATE}', 
end_date = '${FIRE_END_DATE}',
fire_directory = './',
wrf_directory = './',
wrf2fire_map = 'co -> CO', 'no -> NO', 'so2 -> SO2', 'bigalk -> BIGALK',
               'bigene -> BIGENE', 'c2h4 -> C2H4', 'c2h5oh -> C2H5OH',
               'c2h6 -> C2H6', 'c3h8 -> C3H8','c3h6 -> C3H6','ch2o -> CH2O', 'ch3cho -> CH3CHO',
               'ch3coch3 -> CH3COCH3','ch3oh -> CH3OH','mek -> MEK','toluene -> TOLUENE',
               'nh3 -> NH3','no2 -> NO2','open -> BIGALD','c10h16 -> C10H16',
               'ch3cooh -> CH3COOH','cres -> CRESOL','glyald -> GLYALD','mgly -> CH3COCHO',
               'gly -> CH3COCHO','acetol -> HYAC','isop -> ISOP','macr -> MACR'
               'mvk -> MVK',
               'oc -> OC;aerosol','bc -> BC;aerosol'
/
EOF
#
      RANDOM=$$
      export TYPE=SERIAL
      export EXECUTE=fire_emis.exe 
      if [[ ${TYPE} == PARALLEL ]]; then
         echo APM ERROR: CANNOT RUN MPI JOBS IN THIS CONTEXT
         exit 
      elif [[ ${TYPE} == SERIAL ]]; then
        ./${EXECUTE} < fire_emis.mozc.inp > index_${RANDOM}.html 2>&1
      fi
#
      export L_DATE=${DATE}
      while [[ ${L_DATE} -le ${END_DATE} ]]; do
         export L_YYYY=$(echo $L_DATE | cut -c1-4)
         export L_MM=$(echo $L_DATE | cut -c5-6)
         export L_DD=$(echo $L_DATE | cut -c7-8)
         export L_HH=$(echo $L_DATE | cut -c9-10)
         export L_FILE_DATE=${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
         export DD_DATE=${L_YYYY}${L_MM}${L_DD}
#
# TEST WHETHER OUTPUT EXISTS
         export FILE_CR=wrffirechemi_d${CR_DOMAIN}_${L_FILE_DATE}
         export FILE_FR=wrffirechemi_d${FR_DOMAIN}_${L_FILE_DATE}
         if [[ ! -e ${FILE_CR} || ! -e ${FILE_CR} ]]; then
            echo WRFFIRE FAILED
            exit
         else
            echo WRFFIRE SUCCESS
         fi
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 1 2>/dev/null)
      done
   fi
