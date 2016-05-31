#!/bin/bash
#-- Load Experiment Environment Variables -----------------     
. FeoM_SBMT_ENV_VARS.sh
#----------------------------------------------------------     
     TMPLFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_TeMPLaTe.lsf
     SBMTFILE=${RUNDIR}/FeoM_SBMT_eNS_MeMBeRS_${EXPINFO}.lsf
     cat ${TMPLFILE} > ${SBMTFILE}

echo "\
##########################################################\\
#    SOME DISTINGUISHING INFORMATION ABOUT THE EXPERIMENT:\\
#                                                         \\
# ${EXPINFO} is similar to FB002 but without assimilation \\
#  Here we only evaluate the forecast with all available  \\
#  data.                                                  \\
#                                                         \\
#  The observations are genereted from NR001 nature       \\
#  run. Temperature values are perturbed by std=0.1 degC  \\
#  around mean=0. Similarlar, salinity values are         \\
#  perturbed by std=0.04 psu around mean=0.               \\
#                                                         \\
#  Assimilation cycle is set to 6 hours and observations  \\
#  with +- 3hr will be assimilated in each cycles.        \\
#                                                         \\
#  Synthetic obs. along a cross-section are created       \\
#  for validation. These obs. won't                       \\
#  get assimilated but will appear in obs_seq.final       \\
#  for evaluation of the assimilation experiment.         \\
#                                                         \\
#  Obs. error for temperature is set to 0.5 degC and      \\
#  0.25 for salinity following Grayet etal 2011           \\
#                                                         \\
#  Model started from 2009, January 01  00:00             \\
#                                                         \\
#                                                         \\
#                                                         \\
##########################################################" >> \
${SBMTFILE}
bsub < ${SBMTFILE}
