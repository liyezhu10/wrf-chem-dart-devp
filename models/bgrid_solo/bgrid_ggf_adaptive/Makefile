# Makefile created by mkmf $Id: mkmf 4945 2011-06-02 22:29:30Z thoar $ 

include ../../../mkmf/mkmf.template

CPPDEFS = -Duse_netCDF

.DEFAULT:
	-touch $@
all: bgrid_ggf
bgrid_cold_start.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90 types_mod.o bgrid_horiz.o bgrid_halo.o topography.o utilities_mod.o fms.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90
bgrid_conserve_energy.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o time_manager_mod.o fms.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90
meteor_mod.o: ../../../observations/MADIS/meteor_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/MADIS/meteor_mod.f90
vert_advection.o: ../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90 types_mod.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90
bgrid_masks.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90 types_mod.o bgrid_horiz.o bgrid_halo.o bgrid_vert.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90
model_mod.o: ../../../models/bgrid_solo/model_mod.f90 bgrid_core_driver.o bgrid_prog_var.o bgrid_horiz.o time_manager_mod.o fms.o bgrid_change_grid.o bgrid_vert.o bgrid_halo.o hs_forcing.o location_mod.o random_seq_mod.o types_mod.o utilities_mod.o obs_kind_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/model_mod.f90
obs_kind_mod.o: ../../../obs_kind/obs_kind_mod.f90 utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_kind/obs_kind_mod.f90
bgrid_horiz.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90 types_mod.o constants.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90
atmosphere.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90 types_mod.o bgrid_core_driver.o bgrid_prog_var.o bgrid_horiz.o time_manager_mod.o utilities_mod.o fms.o bgrid_change_grid.o bgrid_vert.o bgrid_halo.o hs_forcing.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90
fms_io.o: ../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90 utilities_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90
atmos_tracer_driver.o: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90 types_mod.o fms.o time_manager_mod.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o atmos_radon.o atmos_carbon_aerosol.o atmos_sulfur_hex.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90
random_seq_mod.o: ../../../random_seq/random_seq_mod.f90 types_mod.o utilities_mod.o random_nr_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_seq/random_seq_mod.f90
tracer_manager.o: ../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90 types_mod.o fms.o field_manager.o time_manager_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90
bgrid_advection.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_polar_filter.o bgrid_halo.o bgrid_change_grid.o vert_advection.o fms.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90
utilities_mod.o: ../../../utilities/utilities_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../utilities/utilities_mod.f90
time_manager_mod.o: ../../../time_manager/time_manager_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../time_manager/time_manager_mod.f90
topography.o: ../../../models/bgrid_solo/fms_src/shared/topography/topography.f90 types_mod.o gaussian_topog.o horiz_interp.o utilities_mod.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/topography/topography.f90
bgrid_sponge.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90 types_mod.o bgrid_horiz.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o bgrid_halo.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90
bgrid_horiz_adjust.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_halo.o bgrid_change_grid.o constants.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90
horiz_interp.o: ../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90 types_mod.o fms.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90
bgrid_vert.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90 types_mod.o constants.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90
obs_def_mod.o: ../../../obs_def/obs_def_mod.f90 types_mod.o utilities_mod.o location_mod.o time_manager_mod.o assim_model_mod.o obs_kind_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_def/obs_def_mod.f90
fft.o: ../../../models/bgrid_solo/fms_src/shared/fft/fft.F90 types_mod.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/fft/fft.F90
atmos_sulfur_hex.o: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90 types_mod.o fms.o time_manager_mod.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90
bgrid_prog_var.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_halo.o bgrid_cold_start.o fms.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90
bgrid_core_driver.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90 types_mod.o bgrid_core.o bgrid_horiz.o bgrid_vert.o bgrid_prog_var.o bgrid_halo.o bgrid_diagnostics.o bgrid_integrals.o bgrid_conserve_energy.o field_manager.o tracer_manager.o time_manager_mod.o utilities_mod.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90
location_mod.o: ../../../location/threed_sphere/location_mod.f90 types_mod.o utilities_mod.o random_seq_mod.o obs_kind_mod.o null_mpi_utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../location/threed_sphere/location_mod.f90
bgrid_ggf.o: bgrid_ggf.f90 types_mod.o utilities_mod.o null_mpi_utilities_mod.o meteor_mod.o ncep_obs_err_mod.o dewpoint_obs_err_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	bgrid_ggf.f90
bgrid_horiz_diff.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_halo.o bgrid_change_grid.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90
gaussian_topog.o: ../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90 types_mod.o fms.o constants.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90
bgrid_change_grid.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90 types_mod.o bgrid_horiz.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90
assim_model_mod.o: ../../../assim_model/assim_model_mod.f90 types_mod.o location_mod.o time_manager_mod.o utilities_mod.o model_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../assim_model/assim_model_mod.f90
hs_forcing.o: ../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90 types_mod.o constants.o utilities_mod.o fms.o time_manager_mod.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90
constants.o: ../../../models/bgrid_solo/fms_src/shared/constants/constants.f90 types_mod.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/constants/constants.f90
bgrid_vert_adjust.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90 types_mod.o bgrid_vert.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90
bgrid_integrals.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90 types_mod.o bgrid_change_grid.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o time_manager_mod.o utilities_mod.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90
bgrid_diagnostics.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90 types_mod.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_prog_var.o bgrid_change_grid.o time_manager_mod.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90
fms.o: ../../../models/bgrid_solo/fms_src/shared/fms/fms.f90 types_mod.o fms_io.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/fms/fms.f90
null_mpi_utilities_mod.o: ../../../mpi_utilities/null_mpi_utilities_mod.f90 types_mod.o utilities_mod.o time_manager_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../mpi_utilities/null_mpi_utilities_mod.f90
field_manager.o: ../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90 ../../../models/bgrid_solo/fms_src/shared/field_manager/parse.inc types_mod.o fms.o utilities_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I../../../models/bgrid_solo/fms_src/shared/field_manager	../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90
atmos_radon.o: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90 types_mod.o fms.o time_manager_mod.o tracer_manager.o field_manager.o atmos_tracer_utilities.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90
bgrid_polar_filter.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90 types_mod.o bgrid_horiz.o fft.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90
types_mod.o: ../../../common/types_mod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../common/types_mod.f90
atmos_carbon_aerosol.o: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90 types_mod.o fms.o time_manager_mod.o tracer_manager.o field_manager.o atmos_tracer_utilities.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90
atmos_tracer_utilities.o: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90 types_mod.o fms.o time_manager_mod.o tracer_manager.o field_manager.o constants.o horiz_interp.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90
fft99.o: ../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90 types_mod.o constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90
bgrid_core.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90 types_mod.o bgrid_prog_var.o bgrid_horiz.o bgrid_vert.o bgrid_masks.o bgrid_advection.o bgrid_horiz_diff.o bgrid_horiz_adjust.o bgrid_vert_adjust.o bgrid_polar_filter.o bgrid_halo.o bgrid_sponge.o fms.o constants.o field_manager.o tracer_manager.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90
dewpoint_obs_err_mod.o: ../../../observations/obs_error/dewpoint_obs_err_mod.f90 types_mod.o meteor_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/obs_error/dewpoint_obs_err_mod.f90
ncep_obs_err_mod.o: ../../../observations/obs_error/ncep_obs_err_mod.f90 types_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../observations/obs_error/ncep_obs_err_mod.f90
bgrid_halo.o: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90 types_mod.o bgrid_horiz.o fms.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90
obs_sequence_mod.o: ../../../obs_sequence/obs_sequence_mod.f90 types_mod.o location_mod.o obs_def_mod.o obs_kind_mod.o time_manager_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../obs_sequence/obs_sequence_mod.f90
random_nr_mod.o: ../../../random_nr/random_nr_mod.f90 types_mod.o utilities_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	../../../random_nr/random_nr_mod.f90
./bgrid_vert_adjust.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90 .
./atmosphere.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90 .
./atmos_carbon_aerosol.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90 .
./bgrid_sponge.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90 .
./fms_io.F90: ../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90
	cp ../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90 .
./parse.inc: ../../../models/bgrid_solo/fms_src/shared/field_manager/parse.inc
	cp ../../../models/bgrid_solo/fms_src/shared/field_manager/parse.inc .
./atmos_radon.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90 .
./bgrid_halo.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90 .
./fft99.f90: ../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90
	cp ../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90 .
./bgrid_polar_filter.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90 .
./bgrid_masks.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90 .
./model_mod.f90: ../../../models/bgrid_solo/model_mod.f90
	cp ../../../models/bgrid_solo/model_mod.f90 .
./assim_model_mod.f90: ../../../assim_model/assim_model_mod.f90
	cp ../../../assim_model/assim_model_mod.f90 .
./atmos_tracer_utilities.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90 .
./bgrid_conserve_energy.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90 .
./vert_advection.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90 .
./bgrid_integrals.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90 .
./time_manager_mod.f90: ../../../time_manager/time_manager_mod.f90
	cp ../../../time_manager/time_manager_mod.f90 .
./atmos_tracer_driver.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90 .
./atmos_sulfur_hex.f90: ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90 .
./obs_kind_mod.f90: ../../../obs_kind/obs_kind_mod.f90
	cp ../../../obs_kind/obs_kind_mod.f90 .
./types_mod.f90: ../../../common/types_mod.f90
	cp ../../../common/types_mod.f90 .
./obs_def_mod.f90: ../../../obs_def/obs_def_mod.f90
	cp ../../../obs_def/obs_def_mod.f90 .
./random_seq_mod.f90: ../../../random_seq/random_seq_mod.f90
	cp ../../../random_seq/random_seq_mod.f90 .
./random_nr_mod.f90: ../../../random_nr/random_nr_mod.f90
	cp ../../../random_nr/random_nr_mod.f90 .
./bgrid_prog_var.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90 .
./obs_sequence_mod.f90: ../../../obs_sequence/obs_sequence_mod.f90
	cp ../../../obs_sequence/obs_sequence_mod.f90 .
./constants.f90: ../../../models/bgrid_solo/fms_src/shared/constants/constants.f90
	cp ../../../models/bgrid_solo/fms_src/shared/constants/constants.f90 .
./hs_forcing.f90: ../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90 .
./tracer_manager.F90: ../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90
	cp ../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90 .
./bgrid_core_driver.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90 .
./bgrid_horiz_diff.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90 .
./location_mod.f90: ../../../location/threed_sphere/location_mod.f90
	cp ../../../location/threed_sphere/location_mod.f90 .
./bgrid_diagnostics.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90 .
./field_manager.F90: ../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90
	cp ../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90 .
./fms.f90: ../../../models/bgrid_solo/fms_src/shared/fms/fms.f90
	cp ../../../models/bgrid_solo/fms_src/shared/fms/fms.f90 .
./bgrid_cold_start.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90 .
./dewpoint_obs_err_mod.f90: ../../../observations/obs_error/dewpoint_obs_err_mod.f90
	cp ../../../observations/obs_error/dewpoint_obs_err_mod.f90 .
./ncep_obs_err_mod.f90: ../../../observations/obs_error/ncep_obs_err_mod.f90
	cp ../../../observations/obs_error/ncep_obs_err_mod.f90 .
./topography.f90: ../../../models/bgrid_solo/fms_src/shared/topography/topography.f90
	cp ../../../models/bgrid_solo/fms_src/shared/topography/topography.f90 .
./fft.F90: ../../../models/bgrid_solo/fms_src/shared/fft/fft.F90
	cp ../../../models/bgrid_solo/fms_src/shared/fft/fft.F90 .
./bgrid_change_grid.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90 .
./null_mpi_utilities_mod.f90: ../../../mpi_utilities/null_mpi_utilities_mod.f90
	cp ../../../mpi_utilities/null_mpi_utilities_mod.f90 .
./gaussian_topog.f90: ../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90
	cp ../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90 .
./meteor_mod.f90: ../../../observations/MADIS/meteor_mod.f90
	cp ../../../observations/MADIS/meteor_mod.f90 .
./bgrid_core.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90 .
./utilities_mod.f90: ../../../utilities/utilities_mod.f90
	cp ../../../utilities/utilities_mod.f90 .
./bgrid_advection.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90 .
./bgrid_horiz_adjust.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90 .
./bgrid_horiz.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90 .
./horiz_interp.f90: ../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90
	cp ../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90 .
./bgrid_vert.f90: ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90
	cp ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90 .
SRC = ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90 ../../../observations/MADIS/meteor_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90 ../../../models/bgrid_solo/model_mod.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90 ../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90 ../../../random_seq/random_seq_mod.f90 ../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90 ../../../utilities/utilities_mod.f90 ../../../time_manager/time_manager_mod.f90 ../../../models/bgrid_solo/fms_src/shared/topography/topography.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90 ../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90 ../../../obs_def/obs_def_mod.f90 ../../../models/bgrid_solo/fms_src/shared/fft/fft.F90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90 ../../../location/threed_sphere/location_mod.f90 bgrid_ggf.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90 ../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90 ../../../assim_model/assim_model_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90 ../../../models/bgrid_solo/fms_src/shared/constants/constants.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90 ../../../models/bgrid_solo/fms_src/shared/fms/fms.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90 ../../../common/types_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90 ../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90 ../../../observations/obs_error/dewpoint_obs_err_mod.f90 ../../../observations/obs_error/ncep_obs_err_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../random_nr/random_nr_mod.f90 ../../../models/bgrid_solo/fms_src/shared/field_manager/parse.inc
OBJ = bgrid_cold_start.o bgrid_conserve_energy.o meteor_mod.o vert_advection.o bgrid_masks.o model_mod.o obs_kind_mod.o bgrid_horiz.o atmosphere.o fms_io.o atmos_tracer_driver.o random_seq_mod.o tracer_manager.o bgrid_advection.o utilities_mod.o time_manager_mod.o topography.o bgrid_sponge.o bgrid_horiz_adjust.o horiz_interp.o bgrid_vert.o obs_def_mod.o fft.o atmos_sulfur_hex.o bgrid_prog_var.o bgrid_core_driver.o location_mod.o bgrid_ggf.o bgrid_horiz_diff.o gaussian_topog.o bgrid_change_grid.o assim_model_mod.o hs_forcing.o constants.o bgrid_vert_adjust.o bgrid_integrals.o bgrid_diagnostics.o fms.o null_mpi_utilities_mod.o field_manager.o atmos_radon.o bgrid_polar_filter.o types_mod.o atmos_carbon_aerosol.o atmos_tracer_utilities.o fft99.o bgrid_core.o dewpoint_obs_err_mod.o ncep_obs_err_mod.o bgrid_halo.o obs_sequence_mod.o random_nr_mod.o
OFF = ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_vert_adjust.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/driver/solo/atmosphere.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_carbon_aerosol.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_sponge.f90 ../../../models/bgrid_solo/fms_src/shared/fms/fms_io.F90 ../../../models/bgrid_solo/fms_src/shared/field_manager/parse.inc ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_radon.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_halo.f90 ../../../models/bgrid_solo/fms_src/shared/fft/fft99.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_polar_filter.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_masks.f90 ../../../models/bgrid_solo/model_mod.f90 ../../../assim_model/assim_model_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_utilities.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_conserve_energy.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/vert_advection/vert_advection.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_integrals.f90 ../../../time_manager/time_manager_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_tracer_driver.f90 ../../../models/bgrid_solo/fms_src/atmos_shared/tracer_driver/atmos_sulfur_hex.f90 ../../../obs_kind/obs_kind_mod.f90 ../../../common/types_mod.f90 ../../../obs_def/obs_def_mod.f90 ../../../random_seq/random_seq_mod.f90 ../../../random_nr/random_nr_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_prog_var.f90 ../../../obs_sequence/obs_sequence_mod.f90 ../../../models/bgrid_solo/fms_src/shared/constants/constants.f90 ../../../models/bgrid_solo/fms_src/atmos_param/hs_forcing/hs_forcing.f90 ../../../models/bgrid_solo/fms_src/shared/tracer_manager/tracer_manager.F90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core_driver.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_diff.f90 ../../../location/threed_sphere/location_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_diagnostics.f90 ../../../models/bgrid_solo/fms_src/shared/field_manager/field_manager.F90 ../../../models/bgrid_solo/fms_src/shared/fms/fms.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_cold_start.f90 ../../../observations/obs_error/dewpoint_obs_err_mod.f90 ../../../observations/obs_error/ncep_obs_err_mod.f90 ../../../models/bgrid_solo/fms_src/shared/topography/topography.f90 ../../../models/bgrid_solo/fms_src/shared/fft/fft.F90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_change_grid.f90 ../../../mpi_utilities/null_mpi_utilities_mod.f90 ../../../models/bgrid_solo/fms_src/shared/topography/gaussian_topog.f90 ../../../observations/MADIS/meteor_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_core.f90 ../../../utilities/utilities_mod.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_advection.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/model/bgrid_horiz_adjust.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_horiz.f90 ../../../models/bgrid_solo/fms_src/shared/horiz_interp/horiz_interp.f90 ../../../models/bgrid_solo/fms_src/atmos_bgrid/tools/bgrid_vert.f90
clean: neat
	-rm -f .cppdefs $(OBJ) bgrid_ggf
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
bgrid_ggf: $(OBJ)
	$(LD) $(OBJ) -o bgrid_ggf  $(LDFLAGS)
