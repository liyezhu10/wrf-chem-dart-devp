#!/bin/bash
EXPCODE=INI; EXPNO=10; EXPYR=2009
EXPINFO=${EXPCODE}${EXPNO};
WRKDIR=/work/ans051/TSS/${EXPINFO}
USRHOM=/users/home/ans051; 
RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts

        TMPLFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_TeMPLaTe.lsf
  	SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf
  sed -e "s;EXPERIMENTNUMBER;${EXPNO};g" -e \
         "s;EXPERIMENTNAME;${EXPCODE};g" -e \
         "s;^EXPYR=.*$;EXPYR=${EXPYR};g" \
	 ${TMPLFILE} > ${SBMTFILE}

echo "\
#############################################################\\
#    SOME DISTINGUISHING INFORMATION ABOUT THE EXPERIMENT:   \\
#                                                            \\
# ${EXPINFO}                                                 \\
#      FIRST TRIAL FOR PERFORMING FULL ASSIMILATION LOOP     \\
#      Model started from 2009, January 01                   \\
#                                                            \\
#                                                            \\
#                                                            \\
#############################################################" >> \
${SBMTFILE}
bsub < ${SBMTFILE}
