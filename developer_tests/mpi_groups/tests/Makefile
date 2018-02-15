# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# Makefile for various mpi test and timing programs
#
# This Makefile directly 'include's the mkmf.template file used by the 
#  mkmf system.  If the mkmf.template is set wrong for this system, 
#  these builds will fail.  See ../../mkmf/mkmf.template.* for various
#  compiler and platform settings.

include ../../../build_templates/mkmf.template

# the include defines the following makefile variables:
# MPIFC = 
# MPILD =
# FC =
# LD =
# NETCDF =
# INCS = 
# FFLAGS += -I$(INCS)  
# LIBS =
# LDFLAGS += -L$(NETCDF)/lib -lnetcdf

# ----------------------------------------------------------------------------

SRC = ftest_groups.f90

ALLSRC = $(SRC) Makefile README runme_*

EXE = ftest_groups
ALLEXE = $(EXE)

# ----------------------------------------------------------------------------

# default target - build and test only what is necessary for dart
all:  $(EXE)


# all executables - if you are having problems, this might help diagnose
# what parts of the install are giving problems
everything: $(ALLEXE)

# ----------------------------------------------------------------------------

# f90 test program that calls MPI functions
ftest_groups: ftest_groups.f90
	$(MPIFC) ftest_groups.f90 $(LDFLAGS) -o ftest_groups

# ----------------------------------------------------------------------------

sanity:
	make ftest_groups.f90
	./ftest_groups
# ----------------------------------------------------------------------------

clean:
	rm -f *.o *.mod $(ALLEXE) out_* *testdata.nc

tar:
	tar -zcvf ~/comm.tar.gz $(ALLSRC)

untar:
	tar -zxvf ~/comm.tar.gz

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
