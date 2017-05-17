!
! $Id $
!

MODULE control_mod

! LF 01/2010
! Remplacement du fichier et common control

  IMPLICIT NONE

  REAL,SAVE :: periodav
  REAL,SAVE :: starttime
  INTEGER,SAVE :: nday ! # of days to run
  INTEGER,SAVE :: day_step ! # of dynamical time steps per day
  INTEGER,SAVE :: iperiod ! make a Matsuno step before avery iperiod-1 LF steps
  INTEGER,SAVE :: iapp_tracvl ! apply (cumulated) traceur advection every
                              ! iapp_tracvl dynamical steps
  INTEGER,SAVE :: nsplit_phys ! number of sub-cycle steps in call to physics
  INTEGER,SAVE :: iconser
  INTEGER,SAVE :: iecri
  INTEGER,SAVE :: dissip_period ! apply dissipation every dissip_period
                                ! dynamical step
  INTEGER,SAVE :: iphysiq ! call physics every iphysiq dynamical steps
  INTEGER,SAVE :: iecrimoy
  INTEGER,SAVE :: dayref
  INTEGER,SAVE :: anneeref ! reference year #
  INTEGER,SAVE :: raz_date
  INTEGER,SAVE :: ip_ebil_dyn
  LOGICAL,SAVE :: offline
  CHARACTER(len=4),SAVE :: config_inca
  CHARACTER(len=10),SAVE :: planet_type ! planet type ('earth','mars',...)
  LOGICAL,SAVE :: output_grads_dyn ! output dynamics diagnostics in
                                   ! binary grads file 'dyn.dat' (y/n)
  LOGICAL,SAVE :: ok_dynzon  ! output zonal transports in dynzon.nc file
  LOGICAL,SAVE ::  ok_dyn_ins ! output instantaneous values of fields
                              ! in the dynamics in NetCDF files dyn_hist*nc
  LOGICAL,SAVE :: ok_dyn_ave ! output averaged values of fields in the dynamics
                             ! in NetCDF files dyn_hist*ave.nc
  LOGICAL,SAVE :: resetvarc  ! allows to reset the variables in sortvarc

!--DART-LMDZ---
  INTEGER :: sub_day_step
  LOGICAL :: assimilation, write_restart_natural, start_from_natural,sub_day
  LOGICAL :: noaa_sst, noaa_sic, variable_dissip, print_sortvar
!--DART-LMDZ---


END MODULE
