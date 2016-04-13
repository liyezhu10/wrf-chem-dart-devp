#!/bin/bash
#BSUB -J ENSFILTER                             # Name of the job.
#BSUB -o LOG/TSSEXPCODEEXPNO_FILTER_%J.out  # Appends std output to file %J.out.
#BSUB -e LOG/TSSEXPCODEEXPNO_FILTER_%J.out  # Appends std error to file %J.err.
#BSUB -q serial_30min            # queue
############################ LSF ###################################
TEMPLATE=TeMPLaTe; COPY='cp -f'; REMOVE='rm -f'; LINK='ln -sf'
JOBDIR=${LS_SUBCWD}                   # directory of this script
JOBIDN=$( echo ${LSB_JOBID} | awk '{ printf("%08d\n", $1) }' ) # job-id
JOBNAM=${LSB_JOBNAME}                 # name of this script
ENSNUM=$( echo 1 | awk '{ printf("%02d\n", $1) }' )
ENSID=$( echo 1 | awk '{ printf("%04d\n", $1) }' )
####################################################################
############### SET EXPERIMENT INFORMATION #########################
####################################################################
EXPNUM=EXPNO; EXPDEF=EXPCODE; EXPINFO=${EXPDEF}${EXPNUM};
EXPYR=EXPERIMENTYEAR; 
ENSDEF=ENSCODE; ENSINFO=${ENSDEF}${ENSNUM}; TOTENS=MEMBER;
ACYCLE=ASSIMILATIONCYCLE
####################################################################
############### DEFINE RELEVANT DIRECTORIES ########################
####################################################################
USRHOM=/users/home/ans051; RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts
FSMHOM=${USRHOM}/FEOM/; FSMPRE=${USRHOM}/FEOM_PREPROC/
MODELHOM=${FSMHOM}/FEOMENS; FSMINI=${FSMPRE}/HINDCAST_IC; 
DRTDIR=${USRHOM}/DART/FEOM/models/FeoM/work
WRKDIR=/work/ans051/TSS/${EXPINFO}; FILDIR=${WRKDIR}/FILTER
ENSDIR=${WRKDIR}/${ENSINFO}; cd ${FILDIR}
###################################################################
######## SET AND COPY THE EXECUTABLES #############################
###################################################################
EXE=filter
FILTER=${DRTDIR}/${EXE}; ${COPY} ${FILTER} .
###################################################################
######## DART INPUT/OUTPUT FILES ##################################
###################################################################
TMPNML=${DRTDIR}/input.nml.${TEMPLATE}
OBSSEQ=${DRTDIR}/obs_seq.out; ${COPY} ${OBSSEQ} .
M2DOUT=filter_ics.${ENSID}
D2MINP=filter_restart.${ENSID}
###################################################################
###################################################################
###################################################################
ZCYCLE=$( echo "${ACYCLE} - 1" | bc )
	sed -e 's;FEOMRSTFILENAME;'${ENSDIR}'/'${ENSINFO}'.'${EXPYR}'.oce.nc;g' -e \
               's;FEOMGRDFILENAME;;g' -e \
               's;INITIALTIME;'${ACYCLE}';g' -e \
               's;BEFOREINIT;'${ZCYCLE}';g' -e \
               's;MODEL2DARTOUTPUT;'${M2DOUT}';g' -e \
               's;DART2MODELINPUT;'${D2MINP}';g' -e \
               's;ENSEMBLENUMBER;'${TOTENS}';g' \
        ${TMPNML} > input.nml
        ./${EXE} 
exit
