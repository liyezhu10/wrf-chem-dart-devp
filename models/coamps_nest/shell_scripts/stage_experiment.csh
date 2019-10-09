#!/bin/tcsh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
#==========================================================================
#
# This script prepares a directory for an assimilation experiment.
# The idea is basically that everything you need should be assembled
# by this script and that this should only be run ONCE per experiment.
# After everything is staged in the experiment directory, another script
# can be run to advance the model and perform the assimilation.
#
# The EXPERIMENTDIR should be structured something like:
# .
# ├── <whatever the coamps executable is called>
# ├── state.vars
# ├── convert.nml
# ├── input.nml
# ├── output_priorinf_mean_d01.nc
# ├── output_priorinf_mean_d02.nc
# ├── output_priorinf_sd_d01.nc
# ├── output_priorinf_sd_d02.nc
# ├── obs_seq_2013011000.out
# ├── sampling_error_correction_table.nc
# ├── trans_coamps_to_dart
# ├── trans_dart_to_coamps
# ├── filter
# ├── run_filter.csh
# ├── instance_0001
# │   ├── coamps_0001_2013011000.hdf5
# │   ├── convert.nml
# │   └── state.vars
# ├── instance_0002
# │   ├── coamps_0002_2013011000.hdf5
# │   ├── convert.nml
# │   └── state.vars
# ...
# └── instance_NNNN
#     ├── coamps_NNNN_2013011000.hdf5
#     ├── convert.nml
#     └── state.vars
#
# This script can be run interactively.
#==========================================================================

# DARTDIR declares the root location of the DART code tree.
# COAMPSDIR declares the location of the COAMPS 'test' code tree
#         because that's where I built my COAMPS and I'm using the
#         default forcing/data files.

switch ("`hostname`")
   case eddy
      set DARTDIR = /home/${USER}/git/DART
      set COAMPSDIR = /home/${USER}/COAMPS
      set ENSEMBLEDIR = /home/${USER}/COAMPS_hdf5_files/Ensemble2018
      set EXPERIMENTDIR = /home/${USER}/${USER}_eddy1/COAMPS_cycling_test2
      set COPY = '/bin/cp --preserve=timestamps -v'
      breaksw
   case ch*
      set DARTDIR = /glade/work/${USER}/DART/rma_coamps
      set COAMPSDIR = /glade/work/${USER}/COAMPS
      set ENSEMBLEDIR = /glade/work/${USER}/COAMPS_hdf5_files/Ensemble2018
      set EXPERIMENTDIR = /glade/scratch/${USER}/COAMPS_cycling_test
      set COPY = '/bin/cp --preserve=timestamps -v'
      breaksw
   default
      breaksw
endsw

set ENSEMBLE_SIZE = 20
set CDTG = 2018051700

if (-e ${EXPERIMENTDIR} ) then
   echo "ERROR: ${EXPERIMENTDIR} already exists."
   echo "Intentionally leaving it alone."
   echo "Must provide a new directory name for the experiment."
#  exit 1
endif

#--------------------------------------------------------------------------
# stage everything in the experiment directory
#--------------------------------------------------------------------------

mkdir -p ${EXPERIMENTDIR}
cd ${EXPERIMENTDIR}

# TJH ${COPY} ${COAMPSDIR}/coamps.exe    ${EXPERIMENTDIR}/.     || exit 1

foreach FILE ( \
         ${ENSEMBLEDIR}/README.txt \
         ${ENSEMBLEDIR}/state.vars \
         ${ENSEMBLEDIR}/convert.nml \
         ${ENSEMBLEDIR}/obs_seq_${CDTG}.out \
         ${ENSEMBLEDIR}/output_priorinf*nc \
         ${DARTDIR}/models/coamps_nest/work/input.nml \
         ${DARTDIR}/models/coamps_nest/work/innov_to_obs_seq \
         ${DARTDIR}/models/coamps_nest/work/filter \
         ${DARTDIR}/models/coamps_nest/work/trans_dart_to_coamps \
         ${DARTDIR}/models/coamps_nest/work/trans_coamps_to_dart \
         ${DARTDIR}/models/coamps_nest/shell_scripts/run_filter.csh.template )

         ${COPY} ${FILE} . || exit 2
end

set SAMPDIR = assimilation_code/programs/gen_sampling_err_table/work
set SAMPFILE = sampling_error_correction_table.nc
${COPY}  ${DARTDIR}/${SAMPDIR}/${SAMPFILE}  .

#--------------------------------------------------------------------------
# customize the user input templates with things that will remain constant
# througout the assimilation. We want to leave the *.template files with
# items that will change with each assimilation cycle.
#--------------------------------------------------------------------------

# Enforce some assumptions.
# There are two nests (DART calls a nest a 'domain' - following WRF nomenclature)
# each nest/domain is read from a separate list of (netCDF) files.
# DART will simply update the netCDF files in-place, so the same lists can be used
# for both input and output.

ex input.nml <<ex_end
g;input_state_files ;s;= .*;= '', '';
g;input_state_file_list ;s;= .*;= 'file_list_domain_1.txt', 'file_list_domain_2.txt';
g;output_state_file_list ;s;= .*;= 'file_list_domain_1.txt', 'file_list_domain_2.txt';
g;ens_size ;s;= .*;= ${ENSEMBLE_SIZE};
g;obs_sequence_in_name ;s;= .*;= 'obs_seq.out';
g;obs_sequence_out_name ;s;= .*;= 'obs_seq.final';
g;num_output_obs_members ;s;= .*;= ${ENSEMBLE_SIZE};
g;num_output_state_members ;s;= .*;= ${ENSEMBLE_SIZE};
g;output_mean ;s;= .*;= .TRUE.;
g;output_sd ;s;= .*;= .TRUE.;
g;perturb_from_single_instance ;s;= .*;= .FALSE.;
g;cdtg ;s;= .*;= "${CDTG}";
wq
ex_end

#--------------------------------------------------------------------------
# Set DART and COAMPS (static) input values.
# There are some replacement strings left in the templates
# that will need to be replaced after each assimilation cycle.
# Things like DSTART and the ININAME ...
#--------------------------------------------------------------------------

# prepare run_filter.csh

sed -e "s#EXPERIMENT_DIRECTORY#${EXPERIMENTDIR}#" run_filter.csh.template >! run_filter.csh
chmod u+x run_filter.csh

#--------------------------------------------------------------------------
# Prepare a set of directories - one for each instance of coamps.
# Usually, everything needed to advance coamps gets copied or linked into
# each of these directories.
#--------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ENSEMBLE_SIZE} )

   set dirname = `printf instance_%03d $member`
   mkdir -p $dirname
   cd $dirname

   set COAMPS_RESTART = `printf coamps_%03d_%s.hdf5 $member $CDTG`

   ${COPY} ../input.nml                      . || exit 4
   ${COPY} ../state.vars                     . || exit 4
   ${COPY} ../convert.nml                    . || exit 4
   ${COPY} ${ENSEMBLEDIR}/${COAMPS_RESTART}  . || exit 4

   cd ..

   @ member++
end

#--------------------------------------------------------------------------
# put some instructions in the experiment directory and echo to screen
#--------------------------------------------------------------------------

cat << ENDOFFILE >! README_assumptions.txt

#--------------------------------------------------------------------------
The EXPERIMENT directory is ${EXPERIMENTDIR}
The DART parts came from    ${DARTDIR}
The COAMPS parts came from  ${COAMPSDIR}
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
The following settings MUST be made to input.nml.template for the logic
in run_filter.csh to work correctly. There are many more that you _may_ want
to change or set to impact the performance of the assimilation.

&filter_nml:
   input_state_file_list        = 'file_list_domain_1.txt', 'file_list_domain_2.txt'
   output_state_file_list       = 'file_list_domain_1.txt', 'file_list_domain_2.txt'
   obs_sequence_in_name         = 'obs_seq.out'
   obs_sequence_out_name        = 'obs_seq.final'

&model_nml
   assimilation_period_days     = <WHATEVER YOU WANT>
   assimilation_period_seconds  = 0
   vert_localization_coord      = 3
   cdtg                         = <WHATEVER IS CORRECT>

&ensemble_manager_nml
   layout         = 2
   tasks_per_node = <YOUR HARDWARE CONFIGURATION>

#--------------------------------------------------------------------------

After these files are confirmed, it should be possible to
run the ${EXPERIMENTDIR}/run_filter.csh script.

ENDOFFILE

cat README_assumptions.txt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

