#!/bin/csh

#===================================================================
# Create test case
#===================================================================

if ($#argv < 2) then
  echo "invoke build_fileter -help for usage."
  exit -1;
endif

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
         [-source_rma1 trunk_source_directory] 
         [-source_rma2 rma_source_directory] 
         [-model model_to_build] 
         [-rundir run-root-directory] 
         [-testcase test-case-directory] 
         [-endian which_endian_to_compile_with] 
         [-type run_r4_or_r8] 
         [-quickbuild true_or_false] 
         [-model_restart template_restart_for_model]
         [-obsfile obs_seq.final_file]
         [-out_stub output_netcdf_restart]
         [-help]

OPTIONS 
     -source_rma1 trunk_source_directory
          root directory for the DART trunk
     -source_rma2 rma_source_directory
          root directory for the rma trunk
     -model model_to_build 
          which model you would like to compile (i.e. wrf, POP, cam, ...)
     -rundir run-root-directory 
          where you would like to run your test case
     -testcase test-case-directory
          test case to copy over
     -model_restart template_restart_for_model
          template restart for model
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
set source_rma2   = "/glade/p/work/hendric/DART/rma_bitwise"
set source_rma1   = "/glade/p/work/hendric/DART/trunk"
set model         = "cam"
set rundir        = "/glade/scratch/hendric/bitwise"
set testcase      = "/glade/p/image/DART_test_cases/wrf/wrf_small"
set endian        = 'little'
set type          = 'r8'
set quickbuild    = 'true'

set model_restart = "model_restart.nc"

set obsfile       = "obs_seq.final"
set out_stub      = "cam_out"

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
    case "-source_rma2"
      set source_rma2 = "$argv[1]"
      breaksw
    case "-source_rma1"
      set source_rma1 = "$argv[1]"
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
    case "-model_restart"
      set model_restart = $argv[1]
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

set basecase  = `basename $testcase`
set branch1   = `basename $source_rma1`
set branch2   = `basename $source_rma2`
set test_rma1 = "test_$branch1" 
set test_rma2 = "test_$branch2" 

echo "rundir    : "$rundir
echo "testcase  : "$testcase
echo "basecase  : "$basecase
echo " "
echo "source_rma1   : "$source_rma1
echo "branch1       : "$branch1
echo "source_rma2   : "$source_rma2  
echo "branch2       : "$branch2
echo "testcase      : "$testcase    
echo "quickbuild    : "$quickbuild      
echo "endian        : "$endian      
echo "type          : "$type      
echo "model         : "$model       
echo "model_restart : "$model_restart
echo "obsfile       : "$obsfile
echo "out_stub      : "$out_stub
echo " " 
echo "test_rma1 : "$test_rma1
echo "test_rma2 : "$test_rma2
echo " " 

if (! -e $rundir) then
  echo "rundir does not exist. making new directory $rundir"
  mkdir $rundir
endif

if (! -e "$rundir/$basecase" ) then
   echo "$rundir/$basecase does not exists"
   echo "copying $testcase to :"
   echo " -> rundir $rundir "
   mkdir $rundir/$basecase
   cp -r $testcase/test_rma $rundir/$basecase/$test_rma1
   cp -r $testcase/test_rma $rundir/$basecase/$test_rma2
   cp -r $testcase/restarts $rundir/$basecase/
endif

# build compare_states for bitwise testing restarts
cd $source_rma1/utilities/test/work
csh mkmf_compare_states >& make.out
make >& make_compare_states.out

foreach BRANCH ($source_rma1 $source_rma2)
   set branch_name = `basename $BRANCH`
   set test_case   = "test_$branch_name" 

   # copy development mkmfs over with fp-model_precise and debug flags
   if ( "$endian" == "big" ) then
      echo "compiling with big endian"
      cp $BRANCH/bitwise/mkmf.template.intel.linux.be.dev $BRANCH/mkmf/mkmf.template
   else
      cp $BRANCH/bitwise/mkmf.template.intel.linux.le.dev $BRANCH/mkmf/mkmf.template
   endif
   
   # turn on bitwise testing
   cd $BRANCH/assim_tools
   sed -i 's/lanai_bitwise = .false./lanai_bitwise = .true./' assim_tools_mod.f90
   
   svn revert $BRANCH/common/types_mod.f90
   
   # run with r4. currently this is only for wrf_reg testcase
   if ("$type" == "r4") then
     echo "running with r4 for $BRANCH"
     cd $BRANCH/common
     sed -i 's/  integer, parameter :: r8 = S/!!integer, parameter :: r8 = S/'   types_mod.f90
     sed -i 's/!!integer, parameter :: r8 = r4/  integer, parameter :: r8 = r4/' types_mod.f90 
   endif
   
   # TODO FIXME : need to check that quickbuild compiled successfully
   if ("$quickbuild" == "true") then
     echo "$BRANCH $model quickbuild : $BRANCH/models/$model/work"
     cd $BRANCH/models/$model/work
     # csh quickbuild.csh -mpi >& build.txt
     csh ~/scripts/make_pf.csh -mpi >& build.txt
   endif
   
   echo "$BRANCH $model : linking filter"
   ln -sf $BRANCH/models/$model/work/filter $rundir/$basecase/$test_case/filter
   
   # stage template restarts
   cd $rundir/$basecase/$test_case/
   
   echo "$test_case $model : staging template restarts"
   ln -sf $BRANCH/bitwise/stage_restarts.csh .
   csh stage_restarts.csh $model_restart $out_stub $test_case 3

   # submit the jobs and wait for it to finish
   echo "submitting filter for $test_case"
   cd $rundir/$basecase/$test_case
   bsub -K < run_filter.csh

   # check that job completed
   set last_log      = `ls -rt *.log | tail -n 1`
   set filter_finish = `grep "Filter done TIME"  $last_log | awk '{print $5}'`
   if ("$filter_finish" == "") then
      echo "filter DID NOT FINISH exiting text"
      exit(0)
   else
      echo "filter FINISHED"
   endif
    
end

# TODO FIXME : would like to submit both jobs in parallel and have
# a better way of checking both jobs have finished before testing 
# differences

echo " "
echo "TESTING DIFFERENCES BETWEEN '$test_rma1' AND '$test_rma2' "
echo " with model restarts   : $out_stub.xxxx.nc              "
echo " and obs sequence file : $obsfile                       "
 
cd $rundir/$basecase/$test_rma1

printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"
printf "|%-33s|%-8s|%-13s|%-13s|\n" " TESTING $basecase" " P/F" " $branch1" " $branch2" 
printf "|%-33s|%-8s|%-13s|%-13s|\n" "---------------------------------" \
                                    "--------" \
                                    "-------------" \
                                    "-------------"

cd $rundir/$basecase
 
set restarts = `ls $test_rma1/$out_stub.* | xargs -n 1 basename`
set priors   = "" #`ls $test_rma1/prior_member.* | xargs -n 1 basename`

set testday  = `ls -l $test_rma1/$out_stub.0001.nc | awk '{print $7}'`
set testhour = `ls -l $test_rma1/$out_stub.0001.nc | awk '{print $8}' | head -c 2`

if ( ! -e compare_states ) then
   ln -sf $source_rma1/utilities/test/work/compare_states compare_states
endif

if ( ! -e input.nml ) then
   cp $test_rma1/input.nml input.nml
   cat ~/scripts/compare_states.nml >> input.nml
endif

# compare asci files
foreach ASCI_FILE ( \
  $obsfile )
  # prior_inflate_restart )

  if ( -e $test_rma1/$ASCI_FILE ) then
     if ( -e $test_rma2/$ASCI_FILE ) then
        set RESULT = "FAILED"
        set rma_time1 = `ls -l $test_rma1/$ASCI_FILE | awk '{print $6 $7","$8}'`
        set rma_time2 = `ls -l $test_rma2/$ASCI_FILE | awk '{print $6 $7","$8}'`
        printf "|%-33s|" " $ASCI_FILE"
        diff -q $test_rma1/$ASCI_FILE $test_rma2/$ASCI_FILE > out.txt && set RESULT = "PASSED"
        printf "%8s|%13s|%13s|\n" "$RESULT " "$rma_time1 " "$rma_time2 "
     else
        printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " $rma_time1" " NO FILE"
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
    set DAY  = `ls -l $test_rma1/$NC_FILE | awk '{print $7}'`
    set HOUR = `ls -l $test_rma1/$NC_FILE | awk '{print $8}' | head -c 2`

    if ("$testday" > "$DAY") then
       set newfile = "FALSE"
    endif

    if (("$testday" == "$DAY") && ("$testhour" > "$HOUR")) then
       set newfile = "FALSE"
    endif

  endif 

  if ($newfile == "TRUE") then 

     if ( -e $test_rma1/$NC_FILE ) then
        if ( -e $test_rma2/$NC_FILE ) then
           set RESULT = "FAILED"
           set rma_time1 = `ls -l $test_rma1/$NC_FILE | awk '{print $6 $7","$8}'`
           set rma_time2 = `ls -l $test_rma2/$NC_FILE | awk '{print $6 $7","$8}'`
           printf "|%-33s|" " $NC_FILE"
           echo $test_rma1/$NC_FILE $test_rma2/$NC_FILE | ./compare_states > out.txt && set RESULT = "PASSED"
           # echo $test_rma1/$NC_FILE $test_rma2/$NC_FILE | ./compare_states && set RESULT = "PASSED"

           # diff -q $test_rma1/$NC_FILE $test_rma/$NC_FILE > out.txt && set RESULT = "PASSED"
           printf "%8s|%13s|%13s|\n" "$RESULT " "$rma_time1 " "$rma_time2 "
        else
           printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " $rma_time1" " NO FILE"
        endif
     else
           printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " NO FILE" ""
     endif

  endif
end
printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"

# if ( -e $test_rma2/*_forward_ope_errors* ) then
#    rm *_forward_ope_errors*
# endif

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
