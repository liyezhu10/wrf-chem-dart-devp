#!/bin/ksh -aeux
#
# Script to interpolate for a single file
export FPATH=/glade/p/FRAPPE/real_FRAPPE_CNTL_VARLOC_REV1.a
export DATE=2014072912
export P_BACK_DATE=2014073000
export P_DATE=2014073003
export P_FORW_DATE=2014073006
export BACK_WT=0.5
#
export P_YYYY=$(echo $P_DATE | cut -c1-4)
export P_MM=$(echo $P_DATE | cut -c5-6)
export P_DD=$(echo $P_DATE | cut -c7-8)
export P_HH=$(echo $P_DATE | cut -c9-10)
export P_FILE_DATE=${P_YYYY}-${P_MM}-${P_DD}_${P_HH}:00:00
#
export P_BK_YYYY=$(echo $P_BACK_DATE | cut -c1-4)
export P_BK_MM=$(echo $P_BACK_DATE | cut -c5-6)
export P_BK_DD=$(echo $P_BACK_DATE | cut -c7-8)
export P_BK_HH=$(echo $P_BACK_DATE | cut -c9-10)
export P_BACK_FILE_DATE=${P_BK_YYYY}-${P_BK_MM}-${P_BK_DD}_${P_BK_HH}:00:00
export P_FW_YYYY=$(echo $P_FORW_DATE | cut -c1-4)
export P_FW_MM=$(echo $P_FORW_DATE | cut -c5-6)
export P_FW_DD=$(echo $P_FORW_DATE | cut -c7-8)
export P_FW_HH=$(echo $P_FORW_DATE | cut -c9-10)
export P_FORW_FILE_DATE=${P_FW_YYYY}-${P_FW_MM}-${P_FW_DD}_${P_FW_HH}:00:00
#
export BACK_MET_CR=met_em.d01.${P_BACK_FILE_DATE}.nc
export BACK_MET_FR=met_em.d02.${P_BACK_FILE_DATE}.nc
export FORW_MET_CR=met_em.d01.${P_FORW_FILE_DATE}.nc
export FORW_MET_FR=met_em.d02.${P_FORW_FILE_DATE}.nc
export BACK_FILE_CR=BK_met_em.d01.${P_BACK_FILE_DATE}.nc
export BACK_FILE_FR=BK_met_em.d02.${P_BACK_FILE_DATE}.nc
export FORW_FILE_CR=FW_met_em.d01.${P_FORW_FILE_DATE}.nc
export FORW_FILE_FR=FW_met_em.d02.${P_FORW_FILE_DATE}.nc
export OUTFILE_CR=met_em.d01.${P_FILE_DATE}.nc
export OUTFILE_FR=met_em.d02.${P_FILE_DATE}.nc
#
cd ${FPATH}/${DATE}/metgrid
rm -rf ${BACK_FILE_CR}
rm -rf ${BACK_FILE_FR}
rm -rf ${FORW_FILE_CR}
rm -rf ${FORW_FILE_FR}
rm -rf ${OUTFILE_CR}
rm -rf ${OUTFILE_FR}
ln -sf ${BACK_MET_CR} ${BACK_FILE_CR}
ln -sf ${BACK_MET_FR} ${BACK_FILE_FR}
ln -sf ${FORW_MET_CR} ${FORW_FILE_CR}
ln -sf ${FORW_MET_FR} ${FORW_FILE_FR}
#
export TIME_INTERP_DIR1=/glade/p/work/mizzi/TRUNK/DART_CHEM_MY_BRANCH/models/wrf_chem
export TIME_INTERP_DIR2=run_scripts/RUN_TIME_INTERP
export FIX_TIME_FILE=${TIME_INTERP_DIR1}/${TIME_INTERP_DIR2}/fix_time_stamp.exe
cp ${FIX_TIME_FILE} ./.
export NUM_FIX_DATES=1
#
# CREATE NAMELIST
rm -rf time_stamp_nml.nl
cat << EOF > time_stamp_nml.nl
&time_stamp_nml
time_str1='${P_FILE_DATE}'
file_str='${OUTFILE_CR}'
num_dates=${NUM_FIX_DATES}
file_sw=0
/
EOF
ncflint -w ${BACK_WT} ${BACK_FILE_CR} ${FORW_FILE_CR} ${OUTFILE_CR} 
./fix_time_stamp.exe
#
# CREATE NAMELIST
rm -rf time_stamp_nml.nl
cat << EOF > time_stamp_nml.nl
&time_stamp_nml
time_str1='${P_FILE_DATE}'
file_str='${OUTFILE_FR}'
num_dates=${NUM_FIX_DATES}
file_sw=0
/
EOF
ncflint -w ${BACK_WT} ${BACK_FILE_FR} ${FORW_FILE_FR} ${OUTFILE_FR}
./fix_time_stamp.exe
#
