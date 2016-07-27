# Makefile created by mkmf $Id: mkmf 4945 2011-06-02 22:29:30Z thoar $ 

include ../../../mkmf/mkmf.template

CPPDEFS = -Duse_netCDF

.DEFAULT:
	-touch $@
all: wrf_emp_localization
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
meteor_mod.o: ../../../observations/MADIS/meteor_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/MADIS/meteor_mod.f90
model_mod.o: ../../../models/wrf/model_mod.f90 types_mod.o time_manager_mod.o location_mod.o utilities_mod.o null_mpi_utilities_mod.o random_seq_mod.o obs_kind_mod.o module_map_utils.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/wrf/model_mod.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
wrf_emp_localization.o: wrf_emp_localization.f90 types_mod.o utilities_mod.o model_mod.o module_map_utils.o null_mpi_utilities_mod.o meteor_mod.o ncep_obs_err_mod.o dewpoint_obs_err_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	wrf_emp_localization.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
module_map_utils.o: ../../../models/wrf/module_map_utils.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/wrf/module_map_utils.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o random_nr_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o time_manager_mod.o assim_model_mod.o obs_kind_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
ncep_obs_err_mod.o: ../../../observations/obs_error/ncep_obs_err_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/obs_error/ncep_obs_err_mod.f90
dewpoint_obs_err_mod.o: ../../../observations/obs_error/dewpoint_obs_err_mod.f90 types_mod.o meteor_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/obs_error/dewpoint_obs_err_mod.f90
random_nr_mod.o: ../../../random_nr/random_nr_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_nr/random_nr_mod.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
location_mod.o: ../../../location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_sphere/location_mod.f90
./random_nr_mod.f90: ../../../random_nr/random_nr_mod.f90
	cp ../../../random_nr/random_nr_mod.f90 .
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
./location_mod.f90: ../../../location/threed_sphere/location_mod.f90
	cp ../../../location/threed_sphere/location_mod.f90 .
./dewpoint_obs_err_mod.f90: ../../../observations/obs_error/dewpoint_obs_err_mod.f90
	cp ../../../observations/obs_error/dewpoint_obs_err_mod.f90 .
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./ncep_obs_err_mod.f90: ../../../observations/obs_error/ncep_obs_err_mod.f90
	cp ../../../observations/obs_error/ncep_obs_err_mod.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./meteor_mod.f90: ../../../observations/MADIS/meteor_mod.f90
	cp ../../../observations/MADIS/meteor_mod.f90 .
./module_map_utils.f90: ../../../models/wrf/module_map_utils.f90
	cp ../../../models/wrf/module_map_utils.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./model_mod.f90: ../../../models/wrf/model_mod.f90
	cp ../../../models/wrf/model_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
SRC = ../../../assim_model/assim_model_mod.f90 ../../../observations/MADIS/meteor_mod.f90 ../../../models/wrf/model_mod.f90 ../../../obs_kind/obs_kind_mod.f90 wrf_emp_localization.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../models/wrf/module_map_utils.f90 ../../../random_seq/random_seq_mod.f90 ../../../common/types_mod.f90 ../../../utilities/utilities_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../observations/obs_error/ncep_obs_err_mod.f90 ../../../observations/obs_error/dewpoint_obs_err_mod.f90 ../../../random_nr/random_nr_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../location/threed_sphere/location_mod.f90
OBJ = assim_model_mod.o meteor_mod.o model_mod.o obs_kind_mod.o wrf_emp_localization.o null_mpi_utilities_mod.o module_map_utils.o random_seq_mod.o types_mod.o utilities_mod.o time_manager_mod.o obs_def_mod.o ncep_obs_err_mod.o dewpoint_obs_err_mod.o random_nr_mod.o obs_sequence_mod.o location_mod.o
OFF = ../../../random_nr/random_nr_mod.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../location/threed_sphere/location_mod.f90 ../../../observations/obs_error/dewpoint_obs_err_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../observations/obs_error/ncep_obs_err_mod.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../observations/MADIS/meteor_mod.f90 ../../../models/wrf/module_map_utils.f90 ../../../time_manager/time_manager_mod.f90 ../../../utilities/utilities_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../common/types_mod.f90 ../../../models/wrf/model_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../random_seq/random_seq_mod.f90
clean: neat
	-rm -f .cppdefs $(OBJ) wrf_emp_localization
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
wrf_emp_localization: $(OBJ)
	$(LD) $(OBJ) -o wrf_emp_localization  $(LDFLAGS)
