#       Makefile fuer LMRET, Juni 2000
# ====================================================
# |  METEOROLOGISCHES INSTITUT DER UNIVERSITAET BONN |
# |  Diplomarbeit von Felix Ament, "Initialisierung  |
# |  von Wolken aus Meteosat-Messungen"              |
# ====================================================



# !!!!!!!!!!! Dieser Teil muss nicht veraendert werden !!!!!!!!!!!!!!!!!!!!!!
# Name des Compilers
f95     = ifort

CC      = $(f95) -C=array
LD      = $(f95) -C=array

NETCDF = /home/jkeller/netcdf

CPPFLAGS =  -traceback -I$(NETCDF)/include
LDFLAGS = -L$(NETCDF)/lib -lnetcdf -lnetcdff

# Bestimmt das *.o aus *.cc erzeugt wird
.SUFFIXES: .o .f90                  
# Regel zum Compilieren:
.f90.o: 
	 $(CC) -c $(CPPFLAGS) $*.f90   
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ***********************  M A K E - R E G E L N ***************************

# Namen der benoetigten Objektfiles:
OBJS = byte_mod.o grib_info_mod.o cosmo_data_mod.o model_mod.o
# Regel zum Erzeugen des Programms

model_mod.o:
	$(CC) -c $(CPPFLAGS) model_mod.f90

location_mod.o:
	$(CC) -c $(CPPFLAGS) ../../location/threed_sphere/location_mod.f90

types_mod.o:
	$(CC) -c $(CPPFLAGS) ../../common/types_mod.f90

utilities_mod.o:
	$(CC) -c $(CPPFLAGS) ../../utilities/utilities_mod.f90

obs_kind_mod.o:
	$(CC) -c $(CPPFLAGS) ../../obs_kind/obs_kind_mod.f90

test:  types_mod.o utilities_mod.o obs_kind_mod.o byte_mod.o grib_info_mod.o cosmo_data_mod.o location_mod.o model_mod.o test.o
	$(LD) ${LDFLAGS} types_mod.o utilities_mod.o obs_kind_mod.o byte_mod.o grib_info_mod.o cosmo_data_mod.o location_mod.o model_mod.o test.o -o test

# =============================================================

clean:
	rm -f *.o *~ *.mod

# *******************  A B H A E N G I G K E I T E N ********************

#src_cini.o: gradsoutput.o data_parameters.o 