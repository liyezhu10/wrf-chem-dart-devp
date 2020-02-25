#!/bin/ksh -x
#########################################################################
#
# Purpose: Create DART &obs_def_AIRNOW_PM10_nml 
#
#########################################################################
#
# Generate namelist section
rm -f input.nml_temp
touch input.nml_temp
cat > input.nml_temp << EOF
 &obs_def_AIRNOW_PM10_nml
   use_log_pm10   = ${NL_USE_LOG_PM10:-.false.},
/ 
EOF
#
# Append namelist section to input.nml
if [[ -f input.nml ]]; then
   cat input.nml_temp >> input.nml
   rm input.nml_temp
else
   mv input.nml_temp input.nml
fi
