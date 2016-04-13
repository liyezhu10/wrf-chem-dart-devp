#!/bin/bash
#BSUB -J ENSINITIALIZE[1-MEMBER]       # Name of the job.
#BSUB -o LOG/TSSEXPCODEEXPNO_INI_%J_%I.out  # Appends stdout to file %J.out.
#BSUB -e LOG/TSSEXPCODEEXPNO_INI_%J_%I.out  # Appends stderr to file %J.err.
#BSUB -q serial_30min            # queue
############################ LSF ###################################
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
ENSNUM=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
####################################################################
############### DEFINE RELEVANT PARAMETERS #########################
####################################################################
EXPYR=EXPERIMENTYEAR
EXPNUM=EXPNO; EXPDEF=EXPCODE; EXPINFO=${EXPDEF}${EXPNUM};
RUNLENGTH=1; ENSDEF=ENSCODE; ENSINFO=${ENSDEF}${ENSNUM};
TOTENS=MEMBER;
JOBID=$(bjobs | grep -ir TSS${EXPINFO} | awk '{print $1}')
NPROC=NCORES;
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051; RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
FSMHOM=${USRHOM}/FEOM/; FSMPRE=${USRHOM}/FEOM_PREPROC/
MODELHOM=${FSMHOM}/FEOMENS; FSMINI=${FSMPRE}/HINDCAST_IC; 
WRKDIR=/work/ans051/TSS/${EXPINFO}; ENSDIR=${WRKDIR}/${ENSINFO}
SUBMITFILE=${WRKDIR}/FeoM_SBMT_iNiTiaLiZe_${EXPINFO}
OUTPUTFILE=TSS${EXPINFO}_INI_${LSB_JOBID}_${LSB_JOBINDEX}.out 
####################################################################
######### MESH DIRECTORIES DEPENDING ON THE PARTITIONING ###########
####################################################################
if   [ ${NPROC} -eq 256 ]; then 
	MESHDIR=mesh-T2G1.5L110b
elif [ ${NPROC} -eq 1024 ]; then 
	MESHDIR=mesh-T2G1.5L110b.V2
elif [ ${NPROC} -eq 512 ]; then 
	MESHDIR=mesh-T2G1.5L110b.V3
elif [ ${NPROC} -eq 128 ]; then 
	MESHDIR=mesh-T2G1.5L110b.V4
fi
####################################################################
######### CREATE ENSEMBLE MEMBER DIRECTORY #########################
######### COPY INITIAL CONDITIONS ##################################
######### SET FEOM NAMELIST PARAMETERS #############################
####################################################################
if [ ! -d ${ENSDIR} ]; then
	mkdir ${ENSDIR}
fi
cd ${ENSDIR}
if [ ! -f  ${ENSINFO}.clock ]; then 
	cp ${FSMINI}/EXPHC.${EXPYR}.clock ${ENSINFO}.clock
	cp ${FSMINI}/EXPHC.${EXPYR}.forcing.diag.nc .
	cp ${FSMINI}/EXPHC.${EXPYR}.oce.diag.nc .
	cp ${FSMINI}/EXPHC.${EXPYR}.ice.diag.nc .
	cp ${FSMINI}/EXPHC.${EXPYR}.oce.mean.nc .
	cp ${FSMINI}/EXPHC.${EXPYR}.ice.mean.nc .
	cp ${FSMINI}/EXPHC.${EXPYR}.ice.nc .
	if [ ${EXPYR} -eq 2008 ]; then
	ncks -F -d T,7,7 /work/ans051/TSS/ZEN01/${ENSINFO}/${ENSINFO}.${EXPYR}.oce.nc ${ENSINFO}.${EXPYR}.oce.nc
	elif [ ${EXPYR} -eq 2009 ]; then
	cp ${FSMINI}/${ENSINFO}.${EXPYR}.oce.nc .
	fi
	for any in EXPHC.${EXPYR}.*.nc; do
		newfile=`echo ${any} | sed 's/EXPHC/'${ENSINFO}'/'`;
		mv ${any} ${newfile};
	done
	sed -e  "s/EXPNUM/${EXPNUM}/" -e "s/EXPDEF/${EXPDEF}/" -e \
	        "s/ENSNUM/${ENSNUM}/" -e "s/ENSDEF/${ENSDEF}/" -e \
	        "s/TIMESTEP/7200/" -e  "s/MESHDIR/${MESHDIR}/" -e \
	        "s/RUNLENGTH/${RUNLENGTH}/" \
		${MODELHOM}/namelist.config.template > ${ENSDIR}/namelist.config
	sed -e  "s/WNDFORCING/MFS/" -e "s/RADFORCING/MFS/" -e \
		"s/PRCFORCING/NCEP/" -e "s/ROFFORCING/KARA/" -e \
		"s/SSSFORCING/MFS/" \
		${MODELHOM}/namelist.forcing.template > ${ENSDIR}/namelist.forcing
	cp ${MODELHOM}/namelist.diag namelist.diag
	cp ${MODELHOM}/namelist.ice namelist.ice
	cp ${MODELHOM}/namelist.oce namelist.oce
        cp ${MODELHOM}/fesom.x fesom.x
fi
