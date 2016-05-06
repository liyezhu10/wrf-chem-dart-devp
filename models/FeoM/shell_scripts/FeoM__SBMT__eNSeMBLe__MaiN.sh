#!/bin/bash
EXPID=FB0; EXPNO=01; EXPYR=2009
EXPINFO=${EXPID}${EXPNO};
WRKDIR=/work/ans051/TSS/${EXPINFO}
USRHOM=/users/home/ans051; 
RUNDIR=${USRHOM}/DART/FEOM/models/FeoM/shell_scripts

        TMPLFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_TeMPLaTe.lsf
  	SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf
  sed -e "s;^EXPID=.*$;EXPID=${EXPID};g" -e \
         "s;^EXPNO=.*$;EXPNO=${EXPNO};g" -e \
         "s;^EXPYR=.*$;EXPYR=${EXPYR};g" \
	 ${TMPLFILE} > ${SBMTFILE}

echo "\
#############################################################\\
#    SOME DISTINGUISHING INFORMATION ABOUT THE EXPERIMENT:   \\
#                                                            \\
# ${EXPINFO}                                                 \\
#  FIRST TRIAL FOR PERFORMING FULL ASSIMILATION LOOP         \\
#  using the Ferrybox synthetic temperature and salinity     \\
#  data. The observations are genereted from NR001 nature    \\
#  run. Temperature values are perturbed by std=0.1 degC     \\
#  around mean=0. Similarlar, salinity values are            \\
#  perturbed by std=0.04 psu around mean=0.                  \\
#                                                            \\
#  Assimilation cycle is set to 6 hours and observations     \\
#  with +- 3hr will be assimilated in each cycles.           \\
#                                                            \\
#  Model started from 2009, January 01  00:00                \\
#                                                            \\
#                                                            \\
#                                                            \\
#############################################################" >> \
${SBMTFILE}
bsub < ${SBMTFILE}
