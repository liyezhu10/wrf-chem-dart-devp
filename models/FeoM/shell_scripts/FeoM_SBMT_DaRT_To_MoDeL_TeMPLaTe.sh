#!/bin/bash
#BSUB -J ENSDART2M[1-ENSEMBLEMEMBERNO]         # Name of the job.
#BSUB -o LOG/TSSD2M_%J_%I.out  # Appends std output to file %J.out.
#BSUB -e LOG/TSSD2M_%J_%I.out  # Appends std error to file %J.err.
#BSUB -q serial_30min            # queue
############################ LSF ###################################
TEMPLATE=TeMPLaTe; COPY='cp -f'; REMOVE='rm -f'; LINK='ln -sf'
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
ENSNO=$( echo ${LSB_JOBINDEX} | awk '{ printf("%02d\n", $1) }' )
ENSN4=$( echo ${LSB_JOBINDEX} | awk '{ printf("%04d\n", $1) }' )
####################################################################
############### SET EXPERIMENT INFORMATION #########################
####################################################################
EXPID=EXPERIMENTNAME
EXPNO=EXPNUMBER
EXPYR=EXPERIMENTYR
MEMNO=ENSEMBLEMEMBERS
ENSID=ENSNAME
NPROC=NUMBEROFCORES
EXPINFO=${EXPID}${EXPNO}
ENSINFO=${ENSID}${ENSNO}
CYCLE=ASSIMILATIONCYCLE
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051
RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
FSMHOM=${USRHOM}/FEOM
FSMPRE=${USRHOM}/FEOM_PREPROC
MODELHOM=${FSMHOM}/FETSSOM.ENS01
FSMINI=${FSMPRE}/HINDCAST_IC
DRTDIR=${USRHOM}/DART/FEOM/models/FeoM/work
WRKDIR=/work/ans051/TSS/${EXPINFO}
FILDIR=${WRKDIR}/FILTER
ENSDIR=${WRKDIR}/${ENSINFO}; cd ${ENSDIR}
###################################################################
######## SET AND COPY THE EXECUTABLES #############################
###################################################################
EXE=dart_to_model
DART2M=${DRTDIR}/${EXE}; ${COPY} ${DART2M} .
###################################################################
######## DART INPUT/OUTPUT FILES ##################################
###################################################################
TMPNML=${DRTDIR}/input.nml.${TEMPLATE}
OBSSEQ=${DRTDIR}/obs_seq.out
M2DOUT=filter_ics.${ENSN4}
D2MINP=filter_restart.${ENSN4}
###################################################################
###################################################################
###################################################################
ANALYSISFILE=${ENSDIR}/${ENSINFO}.${EXPYR}.oce.nc
ZCYCLE=$( echo "${CYCLE} - 1" | bc )
${LINK} ${FILDIR}/${D2MINP} .
        ./${EXE} 
exit
