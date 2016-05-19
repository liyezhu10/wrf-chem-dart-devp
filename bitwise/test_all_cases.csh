#!/bin/csh

## wrf_small
# csh bitwise_test.csh -source_rma     "/glade/p/work/hendric/DART/rma_bitwise" \
#                      -source_trunk   "/glade/p/work/hendric/DART/trunk" \
#                      -model          wrf \
#                      -rundir         "/glade/scratch/hendric/bitwise" \
#                      -testcase       "/glade/p/image/DART_test_cases/wrf/wrf_small" \
#                      -model_to_dart  wrf_to_dart \
#                      -dart_to_model  dart_to_wrf \
#                      -model_restart  wrfinput_d01 \
#                      -dart_restart   dart_wrf_vector \
#                      -trunk_restart  filter_ic_new \
#                      -out_stub       wrf_out

## pop_gx1v6
# csh bitwise_test.csh -source_rma     "/glade/p/work/hendric/DART/rma_bitwise" \
#                      -source_trunk   "/glade/p/work/hendric/DART/trunk" \
#                      -model          POP \
#                      -rundir         "/glade/scratch/hendric/bitwise" \
#                      -testcase       "/glade/p/image/DART_test_cases/pop/pop_gx1v6" \
#                      -model_to_dart  pop_to_dart \
#                      -dart_to_model  dart_to_pop \
#                      -model_restart  pop.r.nc \
#                      -dart_restart   dart_restart \
#                      -trunk_restart  filter_restart \
#                      -out_stub       pop_out.r

## mpas_small
# csh bitwise_test.csh -source_rma     "/glade/p/work/hendric/DART/rma_bitwise" \
#                      -source_trunk   "/glade/p/work/hendric/DART/trunk" \
#                      -model          mpas_atm \
#                      -rundir         "/glade/scratch/hendric/bitwise" \
#                      -testcase       "/glade/p/image/DART_test_cases/pop/pop_gx1v6" \
#                      -model_to_dart  model_to_dart \
#                      -dart_to_model  dart_to_model \
#                      -model_restart  x1.40962.restart.nc \
#                      -dart_restart   dart.ic \
#                      -trunk_restart  filter_restart \
#                      -out_stub       mpas_out

csh bitwise_test.csh -source_rma     "/glade/p/work/hendric/DART/rma_bitwise" \
                     -source_trunk   "/glade/p/work/hendric/DART/trunk" \
                     -model          cam \
                     -rundir         "/glade/scratch/hendric/bitwise" \
                     -testcase       "/glade/p/image/DART_test_cases/cam/cam_fv" \
                     -model_to_dart  cam_to_dart \
                     -dart_to_model  dart_to_cam \
                     -model_restart  caminput.nc \
                     -dart_restart   dart_restart \
                     -trunk_restart  filter_restart \
                     -out_stub       cam_out
