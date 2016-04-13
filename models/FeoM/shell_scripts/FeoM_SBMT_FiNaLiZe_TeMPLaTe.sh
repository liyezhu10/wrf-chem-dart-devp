#!/bin/bash
#BSUB -J TSSENSFINALIZE          # Name of the job.
#BSUB -o LOG/TSSENSFINALIZE_%J.out        # Appends std output to file %J.out.
#BSUB -e LOG/TSSENSFINALIZE_%J.out        # Appends std error  to file %J.out.
#BSUB -q serial_30min            # queue
# LSF
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=${LSB_JOBID}                   # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
####################################################################
############### DEFINE RELEVANT PARAMETERS #########################
####################################################################
EXPNUM=EXPNO; EXPDEF=EXPCODE; EXPINFO=${EXPDEF}${EXPNUM};
TOTENS=MEMBER;
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051; RUNDIR=${USRHOM}/DART/lanai/models/feom_ocn/work
FSMHOM=${USRHOM}/FEOM/; FSMPRE=${USRHOM}/FEOM_PREPROC/
MODELHOM=${FSMHOM}/FEOMENS; 
WRKDIR=/work/ans051/TSS/${EXPINFO}; 

SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf

cd ${WRKDIR}

DAYSTEP=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $2}')
RUNYEAR=$(cat ENS01/ENS01.clock | sed -n 2,2p | awk '{print $3}') 
CHECKFILE=${WRKDIR}/submitcheck.lst
if [ "${RUNYEAR}" -le "2009" ] && [ "${DAYSTEP}" -lt "2" ]; then
  ENSCHECK=$( cat ${CHECKFILE} | wc -l )
  if [ ${ENSCHECK} -eq ${TOTENS} ]; then 
	SUBMITCHECK=$( awk '{SUM += $3} END {print SUM}' ${WRKDIR}/submitcheck.lst ) 
	if [ ${SUBMITCHECK} -eq 0 ]; then 
		cd ${WRKDIR}
		rm ${CHECKFILE}.prev; mv ${CHECKFILE} ${CHECKFILE}.prev
		bsub < ${SBMTFILE}; echo "EXPERIMENT CONTINUES..."
	else 
		echo "One of the ensemble members stopped. Ensemble collapsed..."
	fi
  fi
else
	echo "EXPERIMENT FINISHED..."
fi
exit
