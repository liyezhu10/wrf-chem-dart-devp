# Makefile created by mkmf $Id$ 

include ../../../mkmf/mkmf.template


.DEFAULT:
	-touch $@
all: fixsys perfect_model_obs
fixsys: 
	@ (cd ../../../mpi_utilities; ./fixsystem $(FC) )
obs_model_mod.o: ../../../obs_model/obs_model_mod.f90 utilities_mod.o assim_model_mod.o obs_sequence_mod.o obs_def_mod.o time_manager_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_model/obs_model_mod.f90
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
model_mod.o: ../../../models/mpas_ocn/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o xyz_location_mod.o utilities_mod.o obs_kind_mod.o null_mpi_utilities_mod.o random_seq_mod.o get_geometry_mod.o get_reconstruct_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/mpas_ocn/model_mod.f90
xyz_location_mod.o: ../../../location/threed_cartesian/xyz_location_mod.f90 types_mod.o utilities_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_cartesian/xyz_location_mod.f90
adaptive_inflate_mod.o: ../../../adaptive_inflate/adaptive_inflate_mod.f90 types_mod.o time_manager_mod.o utilities_mod.o random_seq_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../adaptive_inflate/adaptive_inflate_mod.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
perfect_model_obs.o: ../../../perfect_model_obs/perfect_model_obs.f90 types_mod.o utilities_mod.o time_manager_mod.o obs_sequence_mod.o obs_def_mod.o obs_model_mod.o assim_model_mod.o null_mpi_utilities_mod.o random_seq_mod.o ensemble_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../perfect_model_obs/perfect_model_obs.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
assim_tools_mod.o: ../../../assim_tools/assim_tools_mod.f90 types_mod.o utilities_mod.o sort_mod.o random_seq_mod.o obs_sequence_mod.o obs_def_mod.o obs_kind_mod.o cov_cutoff_mod.o reg_factor_mod.o location_mod.o ensemble_manager_mod.o null_mpi_utilities_mod.o adaptive_inflate_mod.o time_manager_mod.o assim_model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_tools/assim_tools_mod.f90
reg_factor_mod.o: ../../../reg_factor/reg_factor_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../reg_factor/reg_factor_mod.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
cov_cutoff_mod.o: ../../../cov_cutoff/cov_cutoff_mod.f90 types_mod.o utilities_mod.o location_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../cov_cutoff/cov_cutoff_mod.f90
get_geometry_mod.o: ../../../models/mpas_ocn/get_geometry_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/mpas_ocn/get_geometry_mod.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o assim_model_mod.o obs_kind_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
ensemble_manager_mod.o: ../../../ensemble_manager/ensemble_manager_mod.f90 types_mod.o utilities_mod.o assim_model_mod.o time_manager_mod.o random_seq_mod.o null_mpi_utilities_mod.o sort_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../ensemble_manager/ensemble_manager_mod.f90
sort_mod.o: ../../../sort/sort_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../sort/sort_mod.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
get_coeff_mod.o: ../../../models/mpas_ocn/get_coeff_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/mpas_ocn/get_coeff_mod.f90
location_mod.o: ../../../location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_sphere/location_mod.f90
get_reconstruct_mod.o: ../../../models/mpas_ocn/get_reconstruct_mod.f90 types_mod.o get_coeff_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/mpas_ocn/get_reconstruct_mod.f90
./perfect_model_obs.f90: ../../../perfect_model_obs/perfect_model_obs.f90
	cp ../../../perfect_model_obs/perfect_model_obs.f90 .
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
./cov_cutoff_mod.f90: ../../../cov_cutoff/cov_cutoff_mod.f90
	cp ../../../cov_cutoff/cov_cutoff_mod.f90 .
./reg_factor_mod.f90: ../../../reg_factor/reg_factor_mod.f90
	cp ../../../reg_factor/reg_factor_mod.f90 .
./xyz_location_mod.f90: ../../../location/threed_cartesian/xyz_location_mod.f90
	cp ../../../location/threed_cartesian/xyz_location_mod.f90 .
./get_coeff_mod.f90: ../../../models/mpas_ocn/get_coeff_mod.f90
	cp ../../../models/mpas_ocn/get_coeff_mod.f90 .
./location_mod.f90: ../../../location/threed_sphere/location_mod.f90
	cp ../../../location/threed_sphere/location_mod.f90 .
./get_geometry_mod.f90: ../../../models/mpas_ocn/get_geometry_mod.f90
	cp ../../../models/mpas_ocn/get_geometry_mod.f90 .
./ensemble_manager_mod.f90: ../../../ensemble_manager/ensemble_manager_mod.f90
	cp ../../../ensemble_manager/ensemble_manager_mod.f90 .
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./model_mod.f90: ../../../models/mpas_ocn/model_mod.f90
	cp ../../../models/mpas_ocn/model_mod.f90 .
./adaptive_inflate_mod.f90: ../../../adaptive_inflate/adaptive_inflate_mod.f90
	cp ../../../adaptive_inflate/adaptive_inflate_mod.f90 .
./sort_mod.f90: ../../../sort/sort_mod.f90
	cp ../../../sort/sort_mod.f90 .
./obs_model_mod.f90: ../../../obs_model/obs_model_mod.f90
	cp ../../../obs_model/obs_model_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./assim_tools_mod.f90: ../../../assim_tools/assim_tools_mod.f90
	cp ../../../assim_tools/assim_tools_mod.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./get_reconstruct_mod.f90: ../../../models/mpas_ocn/get_reconstruct_mod.f90
	cp ../../../models/mpas_ocn/get_reconstruct_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
SRC = ../../../obs_model/obs_model_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../models/mpas_ocn/model_mod.f90 ../../../location/threed_cartesian/xyz_location_mod.f90 ../../../adaptive_inflate/adaptive_inflate_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../perfect_model_obs/perfect_model_obs.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../random_seq/random_seq_mod.f90 ../../../assim_tools/assim_tools_mod.f90 ../../../reg_factor/reg_factor_mod.f90 ../../../utilities/utilities_mod.f90 ../../../common/types_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../models/mpas_ocn/get_geometry_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../sort/sort_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../models/mpas_ocn/get_coeff_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../models/mpas_ocn/get_reconstruct_mod.f90
OBJ = obs_model_mod.o assim_model_mod.o model_mod.o xyz_location_mod.o adaptive_inflate_mod.o obs_kind_mod.o perfect_model_obs.o null_mpi_utilities_mod.o random_seq_mod.o assim_tools_mod.o reg_factor_mod.o utilities_mod.o types_mod.o time_manager_mod.o cov_cutoff_mod.o get_geometry_mod.o obs_def_mod.o ensemble_manager_mod.o sort_mod.o obs_sequence_mod.o get_coeff_mod.o location_mod.o get_reconstruct_mod.o
OFF = ../../../perfect_model_obs/perfect_model_obs.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../cov_cutoff/cov_cutoff_mod.f90 ../../../reg_factor/reg_factor_mod.f90 ../../../location/threed_cartesian/xyz_location_mod.f90 ../../../models/mpas_ocn/get_coeff_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../models/mpas_ocn/get_geometry_mod.f90 ../../../ensemble_manager/ensemble_manager_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../models/mpas_ocn/model_mod.f90 ../../../adaptive_inflate/adaptive_inflate_mod.f90 ../../../sort/sort_mod.f90 ../../../obs_model/obs_model_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../assim_tools/assim_tools_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../utilities/utilities_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../common/types_mod.f90 ../../../models/mpas_ocn/get_reconstruct_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../random_seq/random_seq_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) perfect_model_obs
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
perfect_model_obs: $(OBJ)
	$(LD) $(OBJ) -o perfect_model_obs  $(LDFLAGS)
