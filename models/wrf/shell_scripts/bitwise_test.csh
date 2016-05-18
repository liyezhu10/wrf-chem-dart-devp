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

set basecase = `basename $testcase`

echo "rundir    : $rundir"
echo "testcase  : $testcase"
echo "basecase  : $basecase"
echo " " 

if (! -e $rundir) then
  echo "rundir does not exist. making new directory $rundir"
  mkdir $rundir
endif


if (! -e "$rundir/$basecase" ) then
   echo "$basecase does not exists"
   echo "copying $testcase ."
   cp -r $testcase $rundir
endif

# copy development mkmfs over with fp-model_precise and debug flags
cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_rma/mkmf/mkmf.template
cp $source_rma/bitwise/mkmf.template.intel.linux.le.dev $source_trunk/mkmf/mkmf.template

# turn on bitwise testing
cd $source_rma/assim_tools
sed -i 's/lanai_bitwise = .false./lanai_bitwise = .true./' assim_tools_mod.f90


# TODO FIXME : need to check that quickbuild compiled successfully
echo "building $source_rma"
cd $source_rma/models/$model/work
# csh quickbuild.csh >& build.txt
echo "linking filter to test_rma"
ln -sf $source_rma/models/$model/work/filter $rundir/$basecase/test_rma/filter

echo "building $source_trunk"
cd $source_trunk/models/$model/work
# csh quickbuild.csh >& build.txt
echo "linking filter to test_trunk"
ln -sf $source_trunk/models/$model/work/filter         $rundir/$basecase/test_trunk/filter
echo "linking dart_to_model to test_trunk"
ln -sf $source_trunk/models/$model/work/$dart_to_model $rundir/$basecase/test_trunk/$dart_to_model

echo "submitting filter for test_rma"
cd $rundir/$basecase/test_rma
# bsub < run_filter.csh
echo "submitting filter for test_rma"
cd $rundir/$basecase/test_trunk
# bsub < run_filter.csh

set trunk_restarts = "filter_ic_new"
set obsfile        = "obs_seq.final"
set out_stub       = "wrf_restart"

echo " "
echo "TESTING DIFFERENCES BETWEEN 'test_rma' AND 'test_trunk' "
echo " with trunk restarts   : $trunk_restarts.xxxx           "
echo " and obs sequence file : $obsfile                       "
 
cd $rundir/$basecase/test_trunk

if (-f "$trunk_restarts.0001") then
   echo "                       " 
   echo " CONVERTING DART RESTART FILES TO MODEL FILES "
   echo "                       " 
   ln -sf $source_rma/bitwise/convert_restarts.csh .
   # csh convert_restarts.csh $dart_to_model $dart_restart $model_restart $out_stub $trunk_restarts
   unlink convert_restarts.csh
else
  echo " no restarts "
  exit(0)
endif

printf "|%13s%13s|\n" "-----------------------------------" \
                      "-----------------------------------"
printf "|%33s|%8s|%13s|%13s|\n" "TESTING            " "P/F  " "test_trunk " "test_rma " 
printf "|%33s|%8s|%13s|%13s|\n" "---------------------------------" \
                                "--------" \
                                "-------------" \
                                "-------------"

cd $rundir/$basecase
 
set restarts = `ls test_trunk/$out_stub.*     | xargs -n 1 basename`
set priors   = "" #`ls test_trunk/prior_member.* | xargs -n 1 basename`

set testday  = `ls -l test_trunk/$out_stub.0001.nc | awk '{print $7}'`
set testhour = `ls -l test_trunk/$out_stub.0001.nc | awk '{print $8}' | head -c 2`

if ( ! -e compare_states ) then
   ln -s /glade/u/home/hendric/programs/compare_states compare_states
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
        set TIME1FILE = `ls -l test_trunk/$ASCI_FILE | awk '{print $6 $7","$8}'`
        set TIME2FILE = `ls -l test_rma/$ASCI_FILE   | awk '{print $6 $7","$8}'`
        printf "|%-33s|" " $ASCI_FILE"
        diff -q test_trunk/$ASCI_FILE test_rma/$ASCI_FILE > out.txt && set RESULT = "PASSED"
        printf "%8s|%13s|%13s|\n" "$RESULT " "$TIME1FILE " "$TIME2FILE "
     else
        printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " $TIME1FILE" " NO FILE"
     endif
  else
     printf "|%-33s|%-8s|%-13s|%-13s|\n" " $ASCI_FILE" " DNE" " NO FILE" ""
  endif
  
end 


# # compare .nc files
# foreach NC_FILE ( \
#   $restarts \
#   Prior_Diag.nc \
#   Posterior_Diag.nc )
#   # $priors \
#   # PriorDiag_inf_mean.nc \
#   # PriorDiag_inf_sd.nc \
#   # PriorDiag_mean.nc \
#   # PriorDiag_sd.nc \
#   # prior_inflate_restart )
#   
#   set NEWFILE = "TRUE"
#   if (`echo $NC_FILE | grep "$out_stub"` != "") then
#     set DAY  = `ls -l test_trunk/$NC_FILE | awk '{print $7}'`
#     set HOUR = `ls -l test_trunk/$NC_FILE | awk '{print $8}' | head -c 2`
# 
#     if ("$testday" > "$DAY") set NEWFILE = "FALSE"
#     if (("$testday" == "$DAY") && ("$testhour" > "$HOUR")) set NEWFILE = "FALSE"
#   endif 
# 
#   if ($NEWFILE == "TRUE") then 
# 
#      if ( -e test_trunk/$NC_FILE ) then
#         if ( -e test_rma/$NC_FILE ) then
#            set RESULT = "FAILED"
#            set TIME1FILE = `ls -l test_trunk/$NC_FILE | awk '{print $6 $7","$8}'`
#            set TIME2FILE = `ls -l test_rma/$NC_FILE | awk '{print $6 $7","$8}'`
#            printf "|%-33s|" " $NC_FILE"
#            # echo test_trunk/$NC_FILE test_rma/$NC_FILE | ./compare_states > out.txt && set RESULT = "PASSED"
#            echo test_trunk/$NC_FILE test_rma/$NC_FILE | ./compare_states && set RESULT = "PASSED"
# 
#            # diff -q test_trunk/$NC_FILE $test_rma/$NC_FILE > out.txt && set RESULT = "PASSED"
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
