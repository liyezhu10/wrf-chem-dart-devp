#!/bin/csh 

ifort readhdf5.f90 \
      -g -C -check noarg_temp_created -fpe0 -fp-model precise \
      -ftrapuv -traceback -warn declarations,uncalled,unused \
      -I/Users/thoar/intel_16.0.0/include \
      -L/Users/thoar/intel_16.0.0/lib \
      -lhdf5_hl -lhdf5 -lhdf5hl_fortran -lhdf5_fortran

