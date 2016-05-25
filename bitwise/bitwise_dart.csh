#!/bin/csh

#===================================================================
# Create test case
#===================================================================

if ($#argv < 2) then
  echo "invoke build_fileter -help for usage."
  exit -1;
endif

echo "arguments $argv"

#set verbose
echo " "

# print help messages
set helpheader = 0
set i = `echo $argv[1]|cut -c2-`
if( $i == "help" || $i == "h") then
  set helpheader = 1
endif

#===================================================================

if ( $helpheader ) then
cat << EOF1
NAME   

      bitwise_test - tests bitwise between rma_trunk and DART trunk

SYNOPSIS 

      create_test -testname full-test-name 
         [-source_trunk trunk_source_directory] 
         [-source_rma rma_source_directory] 
         [-model model_to_build] 
         [-rundir run-root-directory] 
         [-testcase test-case-directory] 
         [-endian which_endian_to_compile_with] 
         [-type run_r4_or_r8] 
         [-quickbuild true_or_false] 
         [-model_to_dart model_to_dart]
         [-dart_to_model dart_to_model]
         [-model_restart template_restart_for_model]
         [-dart_restart dart_restart_name_for_model_to_dart]
         [-trunk_restart output_dart_restart]
         [-obsfile obs_seq.final_file]
         [-out_stub output_netcdf_restart]
         [-help]

OPTIONS 
     -source_trunk trunk_source_directory
          root directory for the DART trunk
     -source_rma rma_source_directory
          root directory for the rma trunk
     -model model_to_build 
          which model you would like to compile (i.e. wrf, POP, cam, ...)
     -rundir run-root-directory 
          where you would like to run your test case
     -testcase test-case-directory
          test case to copy over
     -model_to_dart model_to_dart
          model to dart program
     -dart_to_model dart_to_model
          dart to model program
     -model_restart template_restart_for_model
          template restart for model
     -dart_restart dart_restart_name_for_model_to_dart
          dart restart name for model_to_dart
     -trunk_restart output_dart_restart
          output restarts from DART trunk
     -obsfile obs_seq.final_file
          observation sequence final file to test
     -out_stub output_netcdf_restart
          output restart stubs for netcdf
     -endian which_endian_to_compile_with 
          optional either big or little. default is little
     -type run_r4_or_r8 
          optional either r4 or r8. default is r8
     -quickbuild true_or_false 
          optional either true or false. default is true
     -help

EOF1

exit;
endif

#===================================================================
# DEFUALTS
#===================================================================
set source_rma   = "/glade/p/work/hendric/DART/rma_bitwise"
set source_trunk = "/glade/p/work/hendric/DART/trunk"
set model        = "cam"
set rundir       = "/glade/scratch/hendric/bitwise"
set testcase     = "/glade/p/image/DART_test_cases/wrf/wrf_small"
set endian       = 'little'
set type         = 'r8'
set quickbuild   = 'true'


# for CAM
set model_to_dart = "model_to_dart"
set dart_to_model = "dart_to_model"
set model_restart = "model_restart.nc"
set dart_restart  = "dart_restart"

set trunk_restart = "filter_restart"
set obsfile        = "obs_seq.final"
set out_stub       = "cam_out"

while ( 1 )
  if ( $#argv < 1 ) break;
  set i = $argv[1];
### check the error of the input argument.
      shift argv
      if($#argv <1 ) then
          echo "ERROR in ${0}: Please input the content for $i."
          exit -1
      endif
      set dash = "-"
      if( $argv[1] =~ $dash* ) then
          echo "ERROR in ${0}: wrong argument for $i.";
          exit -1
      endif
### end of check the error for the input argument
  switch ( $i )
    case "-source_rma"
      set source_rma = "$argv[1]"
      breaksw
    case "-source_trunk"
      set source_trunk = "$argv[1]"
      breaksw
    case "-model"
      set model = $argv[1]
      breaksw
    case "-rundir"
      set rundir = "$argv[1]"
      breaksw
    case "-testcase"
      set testcase = "$argv[1]"
      breaksw
    case "-model_to_dart"
      set model_to_dart = $argv[1]
      breaksw
    case "-dart_to_model"
      set dart_to_model = $argv[1]
      breaksw
    case "-model_restart"
      set model_restart = $argv[1]
      breaksw
    case "-dart_restart"
      set dart_restart = $argv[1]
      breaksw
    case "-dart_out_restart"
      set dart_out_restart = $argv[1]
      breaksw
    case "-trunk_restart"
      set trunk_restart = $argv[1]
      breaksw
    case "-out_stub"
      set out_stub = $argv[1]
      breaksw
    case "-endian"
      set endian = $argv[1]
      breaksw
    case "-type"
      set type = $argv[1]
      breaksw
    case "-quickbuild"
      set quickbuild = $argv[1]
      breaksw
    default:
      echo "unknown input, invoked create_test with no arguments for usage"
      exit -1
      breaksw
  endsw
  shift argv
end


echo " "

set basecase = `basename $testcase`

echo "rundir    : "$rundir
echo "testcase  : "$testcase
echo "basecase  : "$basecase
echo " "
echo "source_rma        : "$source_rma  
echo "source_trunk      : "$source_trunk
echo "testcase          : "$testcase    
echo "endian            : "$endian      
echo "type              : "$type      
echo "model             : "$model       
echo "model_to_dart     : "$model_to_dart
echo "dart_to_model     : "$dart_to_model
echo "model_restart     : "$model_restart
echo "dart_restart      : "$dart_restart
echo "dart_out_restart  : "$dart_out_restart
echo "trunk_restart     : "$trunk_restart
echo "obsfile           : "$obsfile
echo "out_stub          : "$out_stub
echo " " 

if (! -e $rundir) then
  echo "rundir does not exist. making new directory $rundir"
  mkdir $rundir
endif

if (! -e "$rundir/$basecase" ) then
   echo "$rundir/$basecase does not exists"
   echo "copying $testcase to :"
   echo " -> rundir $rundir "
   cp -r $testcase $rundir
endif

# copy development mkmfs over with fp-model_precise and debug flags
if ( "$endian" == "big" ) then
   echo "compiling with big endian"
   cp $source_rma/bitwise/mkmf.template.intel.linux.be.dev $source_rma/mkmf/mkmf.template
   cp $source_rma/bitwise/mkmf.template.intel.linux.be.dev $source_trunk/mkmf/mkmf.template
else
   cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_rma/mkmf/mkmf.template
   cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_trunk/mkmf/mkmf.template
endif

# turn on bitwise testing
cd $source_rma/assim_tools
sed -i 's/lanai_bitwise = .false./lanai_bitwise = .true./' assim_tools_mod.f90

svn revert $source_rma/common/types_mod.f90
svn revert $source_trunk/common/types_mod.f90

# run with r4. currently this is only for wrf_reg testcase
if ("$type" == "r4") then
  cd $source_rma/common
  sed -i 's/  integer, parameter :: r8 = S/!!integer, parameter :: r8 = S/'   types_mod.f90
  sed -i 's/!!integer, parameter :: r8 = r4/  integer, parameter :: r8 = r4/' types_mod.f90 
  cd $source_trunk/common
  sed -i 's/  integer, parameter :: r8 = S/!!integer, parameter :: r8 = S/'   types_mod.f90
  sed -i 's/!!integer, parameter :: r8 = r4/  integer, parameter :: r8 = r4/' types_mod.f90 
endif

# TODO FIXME : need to check that quickbuild compiled successfully
if ("$quickbuild" == "true") then
  echo "rma      $model quickbuild : $source_rma/models/$model/work"
  cd $source_rma/models/$model/work
  csh quickbuild.csh -mpi >& build.txt
endif

echo "test_rma $model : linking filter"
ln -sf $source_rma/models/$model/work/filter $rundir/$basecase/test_rma/filter

if ("$quickbuild" == "true") then
  echo "trunk    $model quickbuild : $source_trunk/models/$model/work"
  cd $source_trunk/models/$model/work
  csh quickbuild.csh -mpi >& build.txt
endif

# link files to test directories
echo "test_trunk $model : linking filter"
ln -sf $source_trunk/models/$model/work/filter         $rundir/$basecase/test_trunk/filter

echo "test_trunk $model : linking $dart_to_model"
ln -sf $source_trunk/models/$model/work/$dart_to_model $rundir/$basecase/test_trunk/$dart_to_model

echo "test_trunk $model : linking $model_to_dart"
ln -sf $source_trunk/models/$model/work/$model_to_dart $rundir/$basecase/test_trunk/$model_to_dart

# stage template restarts
cd $rundir/$basecase/test_rma/

echo "test_rma $model : staging template restarts"
ln -sf $source_rma/bitwise/stage_restarts.csh .
csh stage_restarts.csh $model_restart $out_stub "test_rma" 3
# csh stage_restarts.csh

cd $rundir/$basecase/test_trunk/

echo "test_trunk $model : staging template restarts"
ln -sf $source_rma/bitwise/stage_restarts.csh .
csh stage_restarts.csh $model_restart $out_stub "test_trunk" 3
# csh stage_restarts.csh

# convert model restarts to dart restarts
echo "test_trunk $model : converting $model restart files to filter restarts"
ln -sf $source_rma/bitwise/convert_model_restarts_to_dart.csh .
csh convert_model_restarts_to_dart.csh $model_to_dart $dart_out_restart $model_restart $out_stub $trunk_restart 

# submit the jobs and wait for it to finish
echo "submitting filter for test_rma"
cd $rundir/$basecase/test_rma
bsub -K < run_filter.csh

# check that job completed
set last_log = `ls -rt *.log | tail -n 1`
set filter_finish = `grep "Filter done TIME"  $last_log | awk '{print $5}'`
if ("$filter_finish" == "") then
   echo "filter DID NOT FINISH exiting text"
   exit(0)
else
   echo "filter FINISHED"
endif
 
echo "submitting filter for test_trunk"
cd $rundir/$basecase/test_trunk
bsub -K < run_filter.csh

# check that job completed
set last_log = `ls -rt *.log | tail -n 1`
set filter_finish = `grep "Filter done TIME"  $last_log | awk '{print $5}'`
if ("$filter_finish" == "") then
   echo "filter DID NOT FINISH exiting text"
   exit(0)
else
   echo "filter FINISHED"
endif

# TODO FIXME : would like to submit both jobs in parallel and have
# a better way of checking both jobs have finished before testing 
# differences

echo " "
echo "TESTING DIFFERENCES BETWEEN 'test_rma' AND 'test_trunk' "
echo " with trunk restarts   : $trunk_restart.xxxx           "
echo " and obs sequence file : $obsfile                       "
 
cd $rundir/$basecase/test_trunk

if (-f "$trunk_restart.0001") then
   echo " " 
   echo " CONVERTING DART RESTART FILES TO MODEL FILES "
   echo " " 
   ln -sf $source_rma/bitwise/convert_dart_restarts_to_model.csh .
   csh convert_dart_restarts_to_model.csh $dart_to_model $dart_restart $model_restart $out_stub $trunk_restart
   unlink convert_dart_restarts_to_model.csh
   echo " " 
else
  echo " no restarts "
  exit(0)
endif

printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"
printf "|%-33s|%-8s|%-13s|%-13s|\n" " TESTING $basecase" " P/F  " " test_trunk " " test_rma " 
printf "|%-33s|%-8s|%-13s|%-13s|\n" "---------------------------------" \
                                    "--------" \
                                    "-------------" \
                                    "-------------"

cd $rundir/$basecase
 
set restarts = `ls test_trunk/$out_stub.* | xargs -n 1 basename`
set priors   = "" #`ls test_trunk/prior_member.* | xargs -n 1 basename`

set testday  = `ls -l test_trunk/$out_stub.0001.nc | awk '{print $7}'`
set testhour = `ls -l test_trunk/$out_stub.0001.nc | awk '{print $8}' | head -c 2`

if ( ! -e compare_states ) then
   ln -sf $source_rma/utilities/test/work/compare_states compare_states
endif

if ( ! -e input.nml ) then
   cp test_trunk/input.nml input.nml
   cat ~/scripts/compare_states.nml >> input.nml
endif

# compare asci files
foreach ASCI_FILE ( \
  $obsfile )
  # prior_inflate_restart )

  if ( -e test_trunk/$ASCI_FILE ) then
     if ( -e test_rma/$ASCI_FILE ) then
        set RESULT = "FAILED"
        set trunk_time = `ls -l test_trunk/$ASCI_FILE | awk '{print $6 $7","$8}'`
        set rma_time   = `ls -l test_rma/$ASCI_FILE   | awk '{print $6 $7","$8}'`
        printf "|%-33s|" " $ASCI_FILE"
        diff -q test_trunk/$ASCI_FILE test_rma/$ASCI_FILE > out.txt && set RESULT = "PASSED"
        printf "%8s|%13s|%13s|\n" "$RESULT " "$trunk_time " "$rma_time "
     else
        printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " $trunk_time" " NO FILE"
     endif
  else
     printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " NO FILE" ""
  endif
  
end 

# compare .nc files
foreach NC_FILE ( \
  $restarts \
  Prior_Diag.nc \
  Posterior_Diag.nc )
  # $priors \
  # PriorDiag_inf_mean.nc \
  # PriorDiag_inf_sd.nc \
  # PriorDiag_mean.nc \
  # PriorDiag_sd.nc \
  # prior_inflate_restart )
  
  set newfile = "TRUE"
  if (`echo $NC_FILE | grep "$out_stub"` != "") then
    set DAY  = `ls -l test_trunk/$NC_FILE | awk '{print $7}'`
    set HOUR = `ls -l test_trunk/$NC_FILE | awk '{print $8}' | head -c 2`

    if ("$testday" > "$DAY") set newfile = "FALSE"
    if (("$testday" == "$DAY") && ("$testhour" > "$HOUR")) set newfile = "FALSE"
  endif 

  if ($newfile == "TRUE") then 

     if ( -e test_trunk/$NC_FILE ) then
        if ( -e test_rma/$NC_FILE ) then
           set RESULT = "FAILED"
           set trunk_time = `ls -l test_trunk/$NC_FILE | awk '{print $6 $7","$8}'`
           set rma_time = `ls -l test_rma/$NC_FILE | awk '{print $6 $7","$8}'`
           printf "|%-33s|" " $NC_FILE"
           echo test_trunk/$NC_FILE test_rma/$NC_FILE | ./compare_states > out.txt && set RESULT = "PASSED"
           # echo test_trunk/$NC_FILE test_rma/$NC_FILE | ./compare_states && set RESULT = "PASSED"

           # diff -q test_trunk/$NC_FILE $test_rma/$NC_FILE > out.txt && set RESULT = "PASSED"
           printf "%8s|%13s|%13s|\n" "$RESULT " "$trunk_time " "$rma_time "
        else
           printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " $trunk_time" " NO FILE"
        endif
     else
        printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " NO FILE" ""
     endif

  endif
end
printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"

# if ( -e test_rma/*_forward_ope_errors* ) then
#    rm *_forward_ope_errors*
# endif
# 
# if ( -e out.txt ) then
#    rm out.txt
# endif
# 
# if (-e compare_states) then
#    unlink compare_states
# endif
# 
# if (-e input.nml) then
#    unlink input.nml
# endif
# 
# if (-e dart_log.nml) then
#    rm dart_log.*
# endif
