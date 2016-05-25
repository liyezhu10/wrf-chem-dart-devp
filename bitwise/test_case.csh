#!/bin/csh

set source_trunk = "/glade/p/work/hendric/DART/trunk"
set source_rma   = "/glade/p/work/hendric/DART/rma_bitwise"
set rundir       = "/glade/scratch/hendric/bitwiseDartRMA"
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
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             wrf \
                           -quickbuild        $quickbuild \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/wrf/wrf_small" \
                           -model_to_dart     wrf_to_dart \
                           -dart_to_model     dart_to_wrf \
                           -model_restart     wrfinput_d01 \
                           -dart_restart      dart_wrf_vector \
                           -dart_out_restart  dart_wrf_vector \
                           -trunk_restart     filter_ics\
                           -trunk_out_restart filter_ic_new \
                           -out_stub          wrf_out
      breaksw
    case "wrf2dom"
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             wrf \
                           -quickbuild        $quickbuild \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/wrf/wrf_small_2dom" \
                           -wrf2dom           true \
                           -model_to_dart     wrf_to_dart \
                           -dart_to_model     dart_to_wrf \
                           -model_restart     wrfinput \
                           -dart_restart      dart_wrf_vector \
                           -dart_out_restart  dart_wrf_vector \
                           -trunk_restart     filter_ics\
                           -trunk_out_restart filter_ic_new \
                           -out_stub          wrf_out
      breaksw
    case "wrf_reg"
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             wrf \
                           -quickbuild        $quickbuild \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/wrf/wrf_regular" \
                           -model_to_dart     wrf_to_dart \
                           -dart_to_model     dart_to_wrf \
                           -model_restart     wrfinput_d01 \
                           -dart_restart      dart_wrf_vector \
                           -dart_out_restart  dart_wrf_vector \
                           -trunk_restart     filter_ics\
                           -trunk_out_restart filter_ics_new \
                           -out_stub          wrf_out \
                           -type              r4
      breaksw
    case "mpas"
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             mpas_atm \
                           -quickbuild        $quickbuild \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/mpas/mpas_small" \
                           -model_to_dart     model_to_dart \
                           -dart_to_model     dart_to_model \
                           -model_restart     x1.40962.restart.nc \
                           -dart_restart      dart.ic \
                           -dart_out_restart  dart.ud \
                           -trunk_restart     filter_ics\
                           -trunk_out_restart filter_restart \
                           -out_stub          mpas_out
      breaksw
    case "cam"
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             cam \
                           -quickbuild        $quickbuild \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/cam/cam_fv" \
                           -model_to_dart     cam_to_dart \
                           -dart_to_model     dart_to_cam \
                           -model_restart     caminput.nc \
                           -dart_restart      dart_restart \
                           -dart_out_restart  dart_ics \
                           -trunk_restart     filter_ics \
                           -trunk_out_restart filter_restart \
                           -out_stub          cam_out
      breaksw
    case "pop"
      csh bitwise_dart.csh -source_rma        $source_rma \
                           -source_trunk      $source_trunk \
                           -model             POP \
                           -rundir            $rundir \
                           -testcase          "/glade/p/image/DART_test_cases/pop/pop_gx1v6" \
                           -endian            big \
                           -model_to_dart     pop_to_dart \
                           -dart_to_model     dart_to_pop \
                           -model_restart     pop.r.nc \
                           -dart_restart      dart_restart \
                           -dart_out_restart  dart_ics \
                           -trunk_restart     filter_ics \
                           -trunk_out_restart filter_restart \
                           -out_stub          pop_out
      breaksw
    default:
      echo "unknown input, invoked create_test with no arguments for usage"
      exit -1
      breaksw
  endsw
  shift argv
end
