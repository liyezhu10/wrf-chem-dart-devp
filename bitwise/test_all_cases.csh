#!/bin/csh

set source_trunk = "/glade/p/work/hendric/DART/trunk"
set source_rma   = "/glade/p/work/hendric/DART/rma_bitwise"
set rundir       = "/glade/scratch/hendric/bitwise4"

# wrf_small
if ("$1" == "wrf") then
   csh bitwise_test.csh -source_rma     $source_rma \
                        -source_trunk   $source_trunk \
                        -model          wrf \
                        -rundir         $rundir \
                        -testcase       "/glade/p/image/DART_test_cases/wrf/wrf_small" \
                        -model_to_dart  wrf_to_dart \
                        -dart_to_model  dart_to_wrf \
                        -model_restart  wrfinput_d01 \
                        -dart_restart   dart_wrf_vector \
                        -trunk_restart  filter_ic_new \
                        -out_stub       wrf_out
endif

# pop_gx1v6
if ("$1" == "pop") then
   csh bitwise_test.csh -source_rma     $source_rma \
                        -source_trunk   $source_trunk \
                        -model          POP \
                        -rundir         $rundir \
                        -testcase       "/glade/p/image/DART_test_cases/pop/pop_gx1v6__2" \
                        -endian         big \
                        -model_to_dart  pop_to_dart \
                        -dart_to_model  dart_to_pop \
                        -model_restart  pop.r.nc \
                        -dart_restart   dart_restart \
                        -trunk_restart  filter_restart \
                        -out_stub       pop_out.r
endif

## mpas_small
if ("$1" == "mpas") then
   csh bitwise_test.csh -source_rma     $source_rma \
                        -source_trunk   $source_trunk \
                        -model          mpas_atm \
                        -rundir         $rundir \
                        -testcase       "/glade/p/image/DART_test_cases/mpas/mpas_small" \
                        -model_to_dart  model_to_dart \
                        -dart_to_model  dart_to_model \
                        -model_restart  x1.40962.restart.nc \
                        -dart_restart   dart.ic \
                        -trunk_restart  filter_restart \
                        -out_stub       mpas_out
endif

## cam_fv
if ("$1" == "cam") then
   csh bitwise_test.csh -source_rma     $source_rma \
                        -source_trunk   $source_trunk \
                        -model          cam \
                        -rundir         $rundir \
                        -testcase       "/glade/p/image/DART_test_cases/cam/cam_fv" \
                        -model_to_dart  cam_to_dart \
                        -dart_to_model  dart_to_cam \
                        -model_restart  caminput.nc \
                        -dart_restart   dart_restart \
                        -trunk_restart  filter_restart \
                        -out_stub       cam_out
endif
