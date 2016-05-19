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

      test_case - bitwise test

SYNOPSIS 

      create_test -testname full-test-name 
         [-compare baseline_name] 
         [-generate baseline_name] 
         [-testroot test-root-directory] 
         [-pes_file PES_file] 
         [-compset_file COMPSET_file] 
         [-testid id] 
         [-inputdataroot input-data-root-directory]
         [-baselineroot baseline_root_directory] 
         [-clean clean_option] 
         [-help]

OPTIONS 

     -baselineroot baseline_root_directory
            Specifies an alternate root directory for baseline datasets
            used for bfb generate/compare testing.  This option is 
            ignored unless baseline generation or comparison is being 
            done.  this will overwrite any env CCSM_BASELINE setting.

            Set full test name including test case, resolution, component set,
            and machines.  Example is ERS.f19_g15.B.bluefire.  USE shortnames
            for best success.  Each testcase can have options appended.  The
            current supported options are
              _D  = debug
              _E  = esmf interfaces
              _P* = pe count setting where * is the pecount (S,M,L,XL,1,etc)
              _R* = regional/single point mode (pts mode) where * is the pt setting (01,02,etc)

     -testroot test-root-directory
            Set the directory where you want the test case created

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


# for CAM
set model_to_dart = "model_to_dart"
set dart_to_model = "dart_to_model"
set model_restart = "model_restart.nc"
set dart_restart  = "dart_restart"

set trunk_restart = "filter_restart"
set obsfile        = "obs_seq.final"
set out_stub       = "cam_out"

echo "arguments $argv"
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
      echo "source_rma"
      set source_rma = "$argv[1]"
      breaksw
    case "-source_trunk"
      echo "source_trunk"
      set source_trunk = "$argv[1]"
      echo "source_trunk after"
      breaksw
    case "-model"
      echo "model"
      set model = $argv[1]
      breaksw
    case "-rundir"
      echo "rundir"
      set rundir = "$argv[1]"
      breaksw
    case "-testcase"
      echo "testcase"
      set testcase = "$argv[1]"
      breaksw
    case "-model_to_dart"
      echo "model_to_dart $argv[1]"
      set model_to_dart = $argv[1]
      breaksw
    case "-dart_to_model"
      echo "dart_to_model $argv[1]"
      set dart_to_model = $argv[1]
      breaksw
    case "-model_restart"
      set model_restart = $argv[1]
      breaksw
    case "-dart_restart"
      set dart_restart = $argv[1]
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
echo "source_rma     : "$source_rma  
echo "source_trunk   : "$source_trunk
echo "model          : "$model       
echo "testcase       : "$testcase    
echo "endian         : "$endian      
echo "model_to_dart  : "$model_to_dart
echo "dart_to_model  : "$dart_to_model
echo "model_restart  : "$model_restart
echo "dart_restart   : "$dart_restart
echo "trunk_restart  : "$trunk_restart
echo "obsfile        : "$obsfile
echo "out_stub       : "$out_stub
echo " " 

if (! -e $rundir) then
  echo "rundir does not exist. making new directory $rundir"
  mkdir $rundir
endif

if (! -e "$rundir/$basecase" ) then
   echo "$basecase does not exists"
   echo "copying $testcase to "
   echo " rundir $rundir "
   cp -r $testcase $rundir
endif

# copy development mkmfs over with fp-model_precise and debug flags
if ( $endian == "big" ) then
   cp $source_rma/bitwise/mkmf.template.intel.linux.be.dev $source_rma/mkmf/mkmf.template
   cp $source_rma/bitwise/mkmf.template.intel.linux.be.dev $source_trunk/mkmf/mkmf.template
else
   cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_rma/mkmf/mkmf.template
   cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_trunk/mkmf/mkmf.template
endif

# turn on bitwise testing
cd $source_rma/assim_tools
sed -i 's/lanai_bitwise = .false./lanai_bitwise = .true./' assim_tools_mod.f90

# TODO FIXME : need to check that quickbuild compiled successfully
echo "quickbuild : $source_rma/models/$model/work"
cd $source_rma/models/$model/work
# csh quickbuild.csh >& build.txt

echo "linking filter to test_rma"
ln -sf $source_rma/models/$model/work/filter $rundir/$basecase/test_rma/filter

echo "quickbuild : $source_trunk/models/$model/work"
cd $source_trunk/models/$model/work
# csh quickbuild.csh >& build.txt

echo "linking filter to test_trunk"
ln -sf $source_trunk/models/$model/work/filter         $rundir/$basecase/test_trunk/filter

echo "linking dart_to_model to test_trunk"
ln -sf $source_trunk/models/$model/work/$dart_to_model $rundir/$basecase/test_trunk/$dart_to_model

echo "linking model_to_dart to test_trunk"
ln -sf $source_trunk/models/$model/work/$model_to_dart $rundir/$basecase/test_trunk/$model_to_dart

cd $rundir/$basecase/test_trunk/

echo "converting cam restart files to filter restarts"
csh convert_restarts_to_dart.csh

echo "copying template restarts"
csh stage_restarts.csh

echo "submitting filter for test_rma"
cd $rundir/$basecase/test_rma
bsub -K < run_filter.csh
echo "submitting filter for test_trunk"
cd $rundir/$basecase/test_trunk
bsub -K < run_filter.csh

# TODO FIXME : need to wait until both jobs finish to continue

echo " "
echo "TESTING DIFFERENCES BETWEEN 'test_rma' AND 'test_trunk' "
echo " with trunk restarts   : $trunk_restart.xxxx           "
echo " and obs sequence file : $obsfile                       "
 
cd $rundir/$basecase/test_trunk

if (-f "$trunk_restart.0001") then
   echo " " 
   echo " CONVERTING DART RESTART FILES TO MODEL FILES "
   echo " " 
   ln -sf $source_rma/bitwise/convert_restarts.csh .
   csh convert_restarts.csh $dart_to_model $dart_restart $model_restart $out_stub $trunk_restart
   unlink convert_restarts.csh
   echo " " 
else
  echo " no restarts "
  exit(0)
endif

printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"
printf "|%-33s|%8s|%13s|%13s|\n" "TESTING $basecase" "P/F  " "test_trunk " "test_rma " 
printf "|%33s|%8s|%13s|%13s|\n" "---------------------------------" \
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
