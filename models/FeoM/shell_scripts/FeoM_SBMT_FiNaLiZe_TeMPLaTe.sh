#!/bin/bash
#BSUB -J TSSENSFINAL         # Name of the job.
#BSUB -o LOG/TSSENSFINAL_%J.out        # Appends std output to file %J.out.
#BSUB -e LOG/TSSENSFINAL_%J.out        # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
# LSF
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=${LSB_JOBID}                   # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
####################################################################
############### DEFINE RELEVANT PARAMETERS #########################
####################################################################
EXPID=EXPERIMENTNAME
EXPNO=EXPNUMBER
EXPINFO=${EXPID}${EXPNO}
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051
RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
FSMHOM=${USRHOM}/FEOM
FSMPRE=${USRHOM}/FEOM_PREPROC
MODELHOM=${FSMHOM}/FETSSOM.ENS01
WRKDIR=/work/ans051/TSS/${EXPINFO}; cd ${WRKDIR}

SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf 

DAYST=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $2}')
RUNYR=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $3}') 
if [ "${RUNYR}" -le "2009" ] && [ "${DAYST}" -lt "3" ]; then
		bsub < ${SBMTFILE}; echo "EXPERIMENT CONTINUES..."
else
	echo "EXPERIMENT FINISHED..."; exit
fi
