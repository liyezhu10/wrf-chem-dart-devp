#!/bin/csh

set source_rma   = "/glade/p/work/hendric/DART/rma_bitwise"
set source_trunk = "/glade/p/work/hendric/DART/trunk"
set model       = "wrf"
set rundir      = "/glade/scratch/hendric/bitwise"
set testcase    = "/glade/p/image/DART_test_cases/wrf/wrf_small"

# for WRF
set model_to_dart = "wrf_to_dart"
set dart_to_model = "dart_to_wrf"
set model_restart = "wrfinput_d01"
set dart_restart  = "dart_wrf_vector"

echo " "

echo "rundir    : $rundir"
echo "testcase  : $testcase"

if (! -e $rundir) then
  echo "rundir does not exist. making new directory $rundir"
  mkdir $rundir
endif

set basecase = `basename $testcase`
echo "basecase $basecase"

if (! -e "$rundir/$basecase" ) then
   echo "$basecase does not exists"
   echo "copying $testcase ."
   cp -r $testcase $rundir
endif


# copy development mkmfs over with fp-model_precise and debug flags
cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_rma/mkmf/mkmf.template
cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_trunk/mkmf/mkmf.template

# TODO FIXME : need to check that quickbuild compiled successfully
echo "building $source_rma"
cd $source_rma/models/$model/work
csh quickbuild.csh >& build.txt
echo "linking filter to test_rma"
ln -sf $source_rma/models/$model/work/filter         $rundir/$basecase/test_rma/filter

echo "building $source_trunk"
cd $source_trunk/models/$model/work
csh quickbuild.csh >& build.txt
echo "linking filter to test_trunk"
ln -sf $source_trunk/models/$model/work/filter         $rundir/$basecase/test_trunk/filter
echo "linking dart_to_model to test_trunk"
ln -sf $source_trunk/models/$model/work/$dart_to_model $rundir/$basecase/test_trunk/$dart_to_model

echo "submitting filter for test_rma"
cd $rundir/test_rma
bsub < run_filter.csh
echo "submitting filter for test_rma"
cd $rundir/test_trunk
bsub < run_filter.csh

exit(0)

set trunk_restarts = "filter_ic_new"
set obsfile  = "obs_seq.final"
set out_stub = "wrf_restart"

echo "TESTING DIFFERENCES BETWEEN 'test_rma' AND 'test_trunk' "
echo " with trunk restarts   : $trunk_restarts.xxxx           "
echo " and obs sequence file : $obsfile                       "
 
cd $rundir/$basecase/test_trunk

if (-f "$trunk_restarts.0001") then
   echo "                       " 
   echo " CONVERTING DART RESTART FILES TO MODEL FILES "
   echo "                       " 
   ln -sf ~/scripts/convert_restarts.csh .
   csh convert_restarts.csh $dart_to_model $dart_restart $model_restart $out_stub
   unlink convert_restarts.csh
else
  echo " no restarts "
  exit(0)
endif

# 
# if (-f "$TEST2/filter_restart.0001") then
#    echo "                       " 
#    echo " '$TEST2' has dart restarts "
#    echo " CONVERTING DART RESTART FILES TO MODEL FILES "
#    echo "                       " 
#    cd $TEST2
#    ln -s ~/scripts/convert_restarts.csh .
#    csh convert_restarts.csh $dart_to_model $dart_restart $model_restart $out_stub
#    unlink convert_restarts.csh
#    cd $TEST_DIR
# endif
# 
# 
# echo "                       " 
# echo "TESTING FROM: $TEST_DIR"
# echo "                       " 
# 
# printf "|%13s%13s|\n" "-----------------------------------" \
#                       "-----------------------------------"
# printf "|%33s|%8s"     "TESTING            " "P/F  "
# printf "|%13s|%13s|\n" "$TEST1 " "$TEST2 " 
# printf "|%33s"    "---------------------------------"
# printf "|%8s"     "--------"
# printf "|%13s"    "-------------"
# printf "|%13s|\n" "-------------"
# 
# set RESTARTS = `ls $TEST1/$out_stub.*     | xargs -n 1 basename`
# set PRIORS   = "" #`ls $TEST1/prior_member.* | xargs -n 1 basename`
# 
# set TESTDAY  = `ls -l $TEST1/$out_stub.0001.nc | awk '{print $7}'`
# set TESTHOUR = `ls -l $TEST1/$out_stub.0001.nc | awk '{print $8}' | head -c 2`
# 
# if ( ! -e compare_states ) then
#    ln -s /glade/u/home/hendric/programs/compare_states compare_states
# endif
# 
# if ( ! -e input.nml ) then
#    cp $TEST1/input.nml input.nml
#    cat ~/scripts/compare_states.nml >> input.nml
# endif
# 
# # compare asci files
# foreach ASCI_FILE ( \
#   $OBSFILE )
#   # prior_inflate_restart )
# 
#   if ( -e $TEST1/$ASCI_FILE ) then
#      if ( -e $TEST2/$ASCI_FILE ) then
#         set RESULT = "FAILED"
#         set TIME1FILE = `ls -l $TEST1/$ASCI_FILE | awk '{print $6 $7","$8}'`
#         set TIME2FILE = `ls -l $TEST2/$ASCI_FILE | awk '{print $6 $7","$8}'`
#         printf "|%-33s|" " $ASCI_FILE"
#         diff -q $TEST1/$ASCI_FILE $TEST2/$ASCI_FILE > out.txt && set RESULT = "PASSED"
#         printf "%8s|%13s|%13s|\n" "$RESULT " "$TIME1FILE " "$TIME2FILE "
#      else
#         printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " $TIME1FILE" " NO FILE"
#      endif
#   else
#      printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " NO FILE" ""
#   endif
#   
# end 
# 
# 
# # compare .nc files
# foreach NC_FILE ( \
#   $RESTARTS \
#   Prior_Diag.nc \
#   Posterior_Diag.nc )
#   # $PRIORS \
#   # PriorDiag_inf_mean.nc \
#   # PriorDiag_inf_sd.nc \
#   # PriorDiag_mean.nc \
#   # PriorDiag_sd.nc \
#   # prior_inflate_restart )
#   
#   set NEWFILE = "TRUE"
#   if (`echo $NC_FILE | grep "$out_stub"` != "") then
#     set DAY  = `ls -l $TEST1/$NC_FILE | awk '{print $7}'`
#     set HOUR = `ls -l $TEST1/$NC_FILE | awk '{print $8}' | head -c 2`
# 
#     if ("$TESTDAY" > "$DAY") set NEWFILE = "FALSE"
#     if (("$TESTDAY" == "$DAY") && ("$TESTHOUR" > "$HOUR")) set NEWFILE = "FALSE"
#   endif 
# 
#   if ($NEWFILE == "TRUE") then 
# 
#      if ( -e $TEST1/$NC_FILE ) then
#         if ( -e $TEST2/$NC_FILE ) then
#            set RESULT = "FAILED"
#            set TIME1FILE = `ls -l $TEST1/$NC_FILE | awk '{print $6 $7","$8}'`
#            set TIME2FILE = `ls -l $TEST2/$NC_FILE | awk '{print $6 $7","$8}'`
#            printf "|%-33s|" " $NC_FILE"
#            # echo $TEST1/$NC_FILE $TEST2/$NC_FILE | ./compare_states > out.txt && set RESULT = "PASSED"
#            echo $TEST1/$NC_FILE $TEST2/$NC_FILE | ./compare_states && set RESULT = "PASSED"
# 
#            # diff -q $TEST1/$NC_FILE $TEST2/$NC_FILE > out.txt && set RESULT = "PASSED"
#            printf "%8s|%13s|%13s|\n" "$RESULT " "$TIME1FILE " "$TIME2FILE "
#         else
#            printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " $TIME1FILE" " NO FILE"
#         endif
#      else
#         printf "|%-33s|%-8s|%-13s|%-13s|\n" " $NC_FILE" " DNE" " NO FILE" ""
#      endif
# 
#   endif
# end
# 
# if ( ! -e compare_states ) then
#    ln -s /glade/u/home/hendric/programs/compare_states compare_states
# endif
# printf "|%13s%13s|\n" "-----------------------------------" \
#                       "-----------------------------------"
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
