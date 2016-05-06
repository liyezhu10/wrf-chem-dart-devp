#!/bin/bash
#BSUB -J TSSENSINI[1-ENSEMBLEMEMBERNO]       # Name of the job.
#BSUB -o LOG/TSSENSINI_%J_%I.out  # Appends stdout to file %J.out.
#BSUB -e LOG/TSSENSINI_%J_%I.out  # Appends stderr to file %J.err.
#BSUB -q serial_30min            # queue
############################ LSF ###################################
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
####################################################################
############### DEFINE RELEVANT PARAMETERS #########################
####################################################################
EXPID=EXPCODE
EXPNO=EXPNO 
EXPYR=EXPERIMENTYEAR
MEMNO=ENSEMBLEMEMBERS
NPROC=NCORES
ENSID=ENSEMBLEID
RNLEN=ASSIMILATIONCYCLE
EXPINFO=${EXPID}${EXPNO};
ENSINFO=${ENSID}${ENSNO};
JOBID=$(bjobs | grep -ir TSS${EXPINFO} | awk '{print $1}')
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051 
RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
DRTDIR=${USRHOM}/DART/FEOM/models/FeoM/work
FSMHOM=${USRHOM}/FEOM
FSMPRE=${USRHOM}/FEOM_PREPROC
MODELHOM=${FSMHOM}/FETSSOM.ENS01
FSMINI=${FSMPRE}/ENSEMBLE_IC 
WRKDIR=/work/ans051/TSS/${EXPINFO}
ENSDIR=${WRKDIR}/${ENSINFO}
SUBMITFILE=${WRKDIR}/FeoM_SBMT_iNiTiaLiZe_${EXPINFO}
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
	if [ ${EXPYR} -eq 2008 ]; then
	ncks -F -d T,7,7 /work/ans051/TSS/ZEN01/${ENSINFO}/${ENSINFO}.${EXPYR}.oce.nc ${ENSINFO}.${EXPYR}.oce.nc
	elif [ ${EXPYR} -eq 2009 ]; then
		INIYR=2008
	cp ${FSMINI}/ENSHC.${INIYR}.clock ${ENSINFO}.clock
	ln -sf ${FSMINI}/ENSHC.${INIYR}.forcing.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.diag.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.oce.mean.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.mean.nc .
	ln -sf ${FSMINI}/ENSHC.${INIYR}.ice.nc .

	ln -sf ${FSMINI}/${ENSINFO}.${INIYR}.oce.nc .
	fi
	for File in ENSHC.${INIYR}.*.nc; do
		NFile=$( echo ${File} | sed 's/ENSHC/'${ENSINFO}'/' );
		mv ${File} ${NFile};
	done
	sed -e  "s/EXPNUM/${EXPNO}/" -e "s/EXPDEF/${EXPID}/" -e \
	        "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/" -e \
	        "s/TIMESTEP/7200/" -e  "s/MESHDIR/${MESHDIR}/" -e \
	        "s/RUNLENGTH/${RNLEN}/" \
		${MODELHOM}/namelist.config.template > ${ENSDIR}/namelist.config
	sed -e  "s/WNDFORCING/MFS/" -e "s/RADFORCING/MFS/" -e \
		"s/PRCFORCING/NCEP/" -e "s/ROFFORCING/KARA/" -e \
		"s/SSSFORCING/MFS/" \
		${MODELHOM}/namelist.forcing.template > ${ENSDIR}/namelist.forcing
	cp ${MODELHOM}/namelist.diag namelist.diag
	cp ${MODELHOM}/namelist.ice namelist.ice
	cp ${MODELHOM}/namelist.oce namelist.oce
	cp ${DRTDIR}/FeoM_time ${ENSDIR}/.
	sed -e  "s/ENSNUM/${ENSNO}/" -e "s/ENSDEF/${ENSID}/"  \
	        ${MODELHOM}/namelist.ensemble > ${ENSDIR}/namelist.ensemble
        cp ${MODELHOM}/fesom.x fesom.x
fi
