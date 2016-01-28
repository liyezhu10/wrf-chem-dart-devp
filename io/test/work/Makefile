# Makefile created by mkmf $Id$ 

include ../../../mkmf/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys test_read_write_restarts
fixsys: 
	@ (cd ../../../mpi_utilities; ./fixsystem $(FC) )
obs_model_mod.o: ../../../obs_model/obs_model_mod.f90 utilities_mod.o assim_model_mod.o state_vector_io_mod.o obs_sequence_mod.o obs_def_mod.o time_manager_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_model/obs_model_mod.f90
null_win_mod.o: ../../../distributed/null_win_mod.f90 types_mod.o ensemble_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../distributed/null_win_mod.f90
state_structure_mod.o: ../../../io/state_structure_mod.f90 utilities_mod.o obs_kind_mod.o types_mod.o sort_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/state_structure_mod.f90
model_mod.o: ../../../models/clm/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o obs_kind_mod.o ensemble_manager_mod.o distributed_state_mod.o state_structure_mod.o null_mpi_utilities_mod.o random_seq_mod.o dart_time_io_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/clm/model_mod.f90
adaptive_inflate_mod.o: ../../../adaptive_inflate/adaptive_inflate_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../adaptive_inflate/adaptive_inflate_mod.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
assim_tools_mod.o: ../../../assim_tools/assim_tools_mod.f90 types_mod.o utilities_mod.o sort_mod.o random_seq_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o cov_cutoff_mod.o reg_factor_mod.o location_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o adaptive_inflate_mod.o time_manager_mod.o assim_model_mod.o model_mod.o distributed_state_mod.o quality_control_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_tools/assim_tools_mod.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
cov_cutoff_mod.o: ../../../cov_cutoff/cov_cutoff_mod.f90 types_mod.o utilities_mod.o location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../cov_cutoff/cov_cutoff_mod.f90
state_space_diag_mod.o: ../../../io/state_space_diag_mod.f90 types_mod.o time_manager_mod.o location_mod.o ensemble_manager_mod.o assim_model_mod.o model_mod.o adaptive_inflate_mod.o null_mpi_utilities_mod.o utilities_mod.o io_filenames_mod.o state_vector_io_mod.o state_structure_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/state_space_diag_mod.f90
quality_control_mod.o: ../../../quality_control/quality_control_mod.f90 types_mod.o utilities_mod.o location_mod.o obs_sequence_mod.o obs_def_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../quality_control/quality_control_mod.f90
smoother_mod.o: ../../../smoother/smoother_mod.f90 types_mod.o null_mpi_utilities_mod.o utilities_mod.o ensemble_manager_mod.o time_manager_mod.o state_space_diag_mod.o assim_tools_mod.o obs_sequence_mod.o adaptive_inflate_mod.o state_vector_io_mod.o io_filenames_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../smoother/smoother_mod.f90
io_filenames_mod.o: ../../../io/io_filenames_mod.f90 utilities_mod.o model_mod.o state_structure_mod.o ensemble_manager_mod.o copies_on_off_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/io_filenames_mod.f90
state_vector_io_mod.o: ../../../io/state_vector_io_mod.f90 adaptive_inflate_mod.o direct_netcdf_mod.o types_mod.o null_mpi_utilities_mod.o ensemble_manager_mod.o utilities_mod.o assim_model_mod.o time_manager_mod.o io_filenames_mod.o model_mod.o copies_on_off_mod.o state_structure_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/state_vector_io_mod.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o ensemble_manager_mod.o obs_kind_mod.o assim_model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
sort_mod.o: ../../../sort/sort_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../sort/sort_mod.f90
dart_time_io_mod.o: ../../../io/dart_time_io_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/dart_time_io_mod.f90
location_mod.o: ../../../location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_sphere/location_mod.f90
filter_mod.o: ../../../filter/filter_mod.f90 types_mod.o obs_sequence_mod.o obs_def_mod.o time_manager_mod.o utilities_mod.o assim_model_mod.o assim_tools_mod.o obs_model_mod.o ensemble_manager_mod.o adaptive_inflate_mod.o null_mpi_utilities_mod.o smoother_mod.o random_seq_mod.o state_vector_io_mod.o io_filenames_mod.o state_structure_mod.o forward_operator_mod.o quality_control_mod.o state_space_diag_mod.o copies_on_off_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../filter/filter_mod.f90
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o ensemble_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
test_read_write_restarts.o: ../../../io/test/test_read_write_restarts.f90 types_mod.o utilities_mod.o null_mpi_utilities_mod.o assim_model_mod.o state_vector_io_mod.o ensemble_manager_mod.o io_filenames_mod.o time_manager_mod.o filter_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/test/test_read_write_restarts.f90
copies_on_off_mod.o: ../../../io/copies_on_off_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/copies_on_off_mod.f90
direct_netcdf_mod.o: ../../../io/direct_netcdf_mod.f90 types_mod.o null_mpi_utilities_mod.o ensemble_manager_mod.o time_manager_mod.o utilities_mod.o state_structure_mod.o io_filenames_mod.o copies_on_off_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../io/direct_netcdf_mod.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
distributed_state_mod.o: ../../../distributed/distributed_state_mod.f90 null_mpi_utilities_mod.o types_mod.o ensemble_manager_mod.o null_win_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../distributed/distributed_state_mod.f90
reg_factor_mod.o: ../../../reg_factor/reg_factor_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../reg_factor/reg_factor_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
forward_operator_mod.o: ../../../forward_operator/forward_operator_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o null_mpi_utilities_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o ensemble_manager_mod.o distributed_state_mod.o quality_control_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../forward_operator/forward_operator_mod.f90
ensemble_manager_mod.o: ../../../ensemble_manager/ensemble_manager_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o null_mpi_utilities_mod.o sort_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../ensemble_manager/ensemble_manager_mod.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
./reg_factor_mod.f90: ../../../reg_factor/reg_factor_mod.f90
	cp ../../../reg_factor/reg_factor_mod.f90 .
./cov_cutoff_mod.f90: ../../../cov_cutoff/cov_cutoff_mod.f90
	cp ../../../cov_cutoff/cov_cutoff_mod.f90 .
./copies_on_off_mod.f90: ../../../io/copies_on_off_mod.f90
	cp ../../../io/copies_on_off_mod.f90 .
./state_structure_mod.f90: ../../../io/state_structure_mod.f90
	cp ../../../io/state_structure_mod.f90 .
./location_mod.f90: ../../../location/threed_sphere/location_mod.f90
	cp ../../../location/threed_sphere/location_mod.f90 .
./state_space_diag_mod.f90: ../../../io/state_space_diag_mod.f90
	cp ../../../io/state_space_diag_mod.f90 .
./dart_time_io_mod.f90: ../../../io/dart_time_io_mod.f90
	cp ../../../io/dart_time_io_mod.f90 .
./ensemble_manager_mod.f90: ../../../ensemble_manager/ensemble_manager_mod.f90
	cp ../../../ensemble_manager/ensemble_manager_mod.f90 .
./model_mod.f90: ../../../models/clm/model_mod.f90
	cp ../../../models/clm/model_mod.f90 .
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./adaptive_inflate_mod.f90: ../../../adaptive_inflate/adaptive_inflate_mod.f90
	cp ../../../adaptive_inflate/adaptive_inflate_mod.f90 .
./io_filenames_mod.f90: ../../../io/io_filenames_mod.f90
	cp ../../../io/io_filenames_mod.f90 .
./sort_mod.f90: ../../../sort/sort_mod.f90
	cp ../../../sort/sort_mod.f90 .
./distributed_state_mod.f90: ../../../distributed/distributed_state_mod.f90
	cp ../../../distributed/distributed_state_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./obs_model_mod.f90: ../../../obs_model/obs_model_mod.f90
	cp ../../../obs_model/obs_model_mod.f90 .
./null_win_mod.f90: ../../../distributed/null_win_mod.f90
	cp ../../../distributed/null_win_mod.f90 .
./assim_tools_mod.f90: ../../../assim_tools/assim_tools_mod.f90
	cp ../../../assim_tools/assim_tools_mod.f90 .
./quality_control_mod.f90: ../../../quality_control/quality_control_mod.f90
	cp ../../../quality_control/quality_control_mod.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./direct_netcdf_mod.f90: ../../../io/direct_netcdf_mod.f90
	cp ../../../io/direct_netcdf_mod.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./state_vector_io_mod.f90: ../../../io/state_vector_io_mod.f90
	cp ../../../io/state_vector_io_mod.f90 .
./test_read_write_restarts.f90: ../../../io/test/test_read_write_restarts.f90
	cp ../../../io/test/test_read_write_restarts.f90 .
./smoother_mod.f90: ../../../smoother/smoother_mod.f90
	cp ../../../smoother/smoother_mod.f90 .
./forward_operator_mod.f90: ../../../forward_operator/forward_operator_mod.f90
	cp ../../../forward_operator/forward_operator_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./filter_mod.f90: ../../../filter/filter_mod.f90
	cp ../../../filter/filter_mod.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
SRC = ../../../obs_model/obs_model_mod.f90 ../../../distributed/null_win_mod.f90 ../../../io/state_structure_mod.f90 ../../../models/clm/model_mod.f90 ../../../adaptive_inflate/adaptive_inflate_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../random_seq/random_seq_mod.f90 ../../../assim_tools/assim_tools_mod.f90 ../../../utilities/utilities_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../io/state_space_diag_mod.f90 ../../../quality_control/quality_control_mod.f90 ../../../smoother/smoother_mod.f90 ../../../io/io_filenames_mod.f90 ../../../io/state_vector_io_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../sort/sort_mod.f90 ../../../io/dart_time_io_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../filter/filter_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../io/test/test_read_write_restarts.f90 ../../../io/copies_on_off_mod.f90 ../../../io/direct_netcdf_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../distributed/distributed_state_mod.f90 ../../../reg_factor/reg_factor_mod.f90 ../../../common/types_mod.f90 ../../../forward_operator/forward_operator_mod.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90
OBJ = obs_model_mod.o null_win_mod.o state_structure_mod.o model_mod.o adaptive_inflate_mod.o obs_kind_mod.o random_seq_mod.o assim_tools_mod.o utilities_mod.o time_manager_mod.o cov_cutoff_mod.o state_space_diag_mod.o quality_control_mod.o smoother_mod.o io_filenames_mod.o state_vector_io_mod.o obs_def_mod.o sort_mod.o dart_time_io_mod.o location_mod.o filter_mod.o assim_model_mod.o test_read_write_restarts.o copies_on_off_mod.o direct_netcdf_mod.o null_mpi_utilities_mod.o distributed_state_mod.o reg_factor_mod.o types_mod.o forward_operator_mod.o ensemble_manager_mod.o obs_sequence_mod.o
OFF = ../../../obs_sequence/obs_sequence_mod.f90 ../../../reg_factor/reg_factor_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../io/copies_on_off_mod.f90 ../../../io/state_structure_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../io/state_space_diag_mod.f90 ../../../io/dart_time_io_mod.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../models/clm/model_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../adaptive_inflate/adaptive_inflate_mod.f90 ../../../io/io_filenames_mod.f90 ../../../sort/sort_mod.f90 ../../../distributed/distributed_state_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../obs_model/obs_model_mod.f90 ../../../distributed/null_win_mod.f90 ../../../assim_tools/assim_tools_mod.f90 ../../../quality_control/quality_control_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../io/direct_netcdf_mod.f90 ../../../utilities/utilities_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../common/types_mod.f90 ../../../io/state_vector_io_mod.f90 ../../../io/test/test_read_write_restarts.f90 ../../../smoother/smoother_mod.f90 ../../../forward_operator/forward_operator_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../filter/filter_mod.f90 ../../../random_seq/random_seq_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) test_read_write_restarts
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
test_read_write_restarts: $(OBJ)
	$(LD) $(OBJ) -o test_read_write_restarts  $(LDFLAGS)
