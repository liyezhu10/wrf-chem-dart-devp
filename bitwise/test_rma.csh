#!/bin/csh

set source_rma1  = "/glade/p/work/hendric/DART/merge_cm1/rma_trunk"
set source_rma2  = "/glade/p/work/hendric/DART/clean_CM1"
set rundir       = "/glade/scratch/hendric/bitwiseCM1"
set quickbuild   = "true"

echo "$argv"

while ( 1 )
  if ( $#argv < 1 ) break;
  set i = $argv[1];
  echo "argument = $i"
  ### end of check the error for the input argument
  switch ( $i )
    case "-quickbuild=false"
      set quickbuild = "false"
      breaksw
    case "wrf"
      csh bitwise_rma.csh  -source_rma1    $source_rma1 \
                           -source_rma2    $source_rma2 \
                           -model          wrf \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/wrf/wrf_small" \
                           -model_restart  wrfinput_d01 \
                           -out_stub       wrf_out
      breaksw
    case "wrf_reg"
      csh bitwise_rma.csh  -source_rma1    $source_rma1 \
                           -source_rma2    $source_rma2 \
                           -model          wrf \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/wrf/wrf_regular" \
                           -model_restart  wrfinput_d01 \
                           -out_stub       wrf_out \
                           -type           r4
      breaksw
    case "mpas"
      csh bitwise_rma.csh  -source_rma1    $source_rma1 \
                           -source_rma2    $source_rma2 \
                           -model          mpas_atm \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/mpas/mpas_small" \
                           -model_restart  x1.40962.restart.nc \
                           -out_stub       mpas_out
      breaksw
    case "cam"
      csh bitwise_rma.csh -source_rma1     $source_rma1 \
                           -source_rma2   $source_rma2 \
                           -model          cam \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/cam/cam_fv" \
                           -model_restart  caminput.nc \
                           -out_stub       cam_out
      breaksw
    case "pop"
      csh bitwise_rma.csh  -source_rma1    $source_rma1 \
                           -source_rma2    $source_rma2 \
                           -model          POP \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/pop/pop_gx1v6" \
                           -endian         big \
                           -model_restart  pop.r.nc \
                           -out_stub       pop_out
      breaksw
    case "cm1"
      csh bitwise_rma.csh  -source_rma1    $source_rma1 \
                           -source_rma2    $source_rma2 \
                           -model          cm1 \
                           -quickbuild     $quickbuild \
                           -rundir         $rundir \
                           -testcase       "/glade/p/image/DART_test_cases/cm1" \
                           -model_restart  cm1out_rst_000001.nc \
                           -out_stub       cm1_out
      breaksw
    default:
      echo "unknown input $i, invoked create_test with no arguments for usage"
      exit -1
      breaksw
  endsw
  shift argv
end
