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
# RUN WRFCHEM_BIO
#
#########################################################################
#
   if ${RUN_WRFCHEM_BIO}; then
      if [[ ! -d ${RUN_DIR}/${DATE}/wrfchem_bio ]]; then
         mkdir ${RUN_DIR}/${DATE}/wrfchem_bio
         cd ${RUN_DIR}/${DATE}/wrfchem_bio
      else
         cd ${RUN_DIR}/${DATE}/wrfchem_bio
      fi
#
# LOOP THROUGHT CURRENT AND NEXT DATE
      export L_DATE=${DATE}
      export LE_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} ${FCST_PERIOD} 2>/dev/null)
      while [[ ${L_DATE} -le ${LE_DATE} ]]; do 
         export L_YYYY=$(echo $L_DATE | cut -c1-4)
         export L_MM=$(echo $L_DATE | cut -c5-6)
         export L_DD=$(echo $L_DATE | cut -c7-8)
         export L_HH=$(echo $L_DATE | cut -c9-10)
         export L_FILE_DATE=${L_YYYY}-${L_MM}-${L_DD}_${L_HH}:00:00
#
# LINK NEEDED FILES
         export FILE_CR=wrfinput_d${CR_DOMAIN}
         export FILE_FR=wrfinput_d${FR_DOMAIN}
         rm -rf ${FILE_CR}
         rm -rf ${FILE_FR}
         cp ${REAL_DIR}/${FILE_CR}_${L_FILE_DATE} ${FILE_CR}   
         cp ${REAL_DIR}/${FILE_FR}_${L_FILE_DATE} ${FILE_FR}   
         export FILE_CR=wrfbiochemi_d${CR_DOMAIN}
         export FILE_FR=wrfbiochemi_d${FR_DOMAIN}
         if [[ ${L_DATE} -eq ${DATE} ]]; then
            rm -rf ${FILE_CR}
            rm -rf ${FILE_FR}
         fi
         rm -rf btr*.nc
         rm -rf DSW*.nc
         rm -rf hrb*.nc
         rm -rf iso*.nc
         rm -rf lai*.nc
         rm -rf ntr*.nc
         rm -rf shr*.nc
         rm -rf TAS*.nc
         cp ${EXPERIMENT_WRFBIOCHEMI_DIR}/MEGAN-DATA/*.nc ./.
         export FILE=megan_bio_emiss.exe
         rm -rf ${FILE}
         cp ${MEGAN_BIO_DIR}/work/${FILE} ${FILE}
#
# CREATE INPUT FILE
         export FILE=megan_bio_emiss.inp
         rm -rf ${FILE}
         cat << EOF > ${FILE}
&control
domains = 2,
start_lai_mnth = 1,
end_lai_mnth = 12
/
EOF
#
         RANDOM=$$
         export TYPE=SERIAL
         export EXECUTE=megan_bio_emiss.exe 
         if [[ ${TYPE} == PARALLEL ]]; then
            echo APM ERROR: CANNOT RUN MPI JOBS IN THIS CONTEXT
            exit 
         elif [[ ${TYPE} == SERIAL ]]; then
           ./${EXECUTE} < megan_bio_emiss.inp > index_${RANDOM}.html 2>&1
         fi
#
# TEST WHETHER OUTPUT EXISTS
         export FILE_CR=wrfbiochemi_d${CR_DOMAIN}
         export FILE_FR=wrfbiochemi_d${FR_DOMAIN}
         if [[ ! -e ${FILE_CR} || ! -e ${FILE_FR} ]]; then
            echo WRFCHEM_BIO FAILED
            exit
         else
            echo WRFCHEM_BIO SUCCESS
            mv ${FILE_CR} ${FILE_CR}_${L_FILE_DATE}
            mv ${FILE_FR} ${FILE_FR}_${L_FILE_DATE}
         fi
         export L_DATE=$(${BUILD_DIR}/da_advance_time.exe ${L_DATE} 6 2>/dev/null)
      done
   fi
