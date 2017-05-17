#!/bin/csh

# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
# 
#**************************************************************************
#   Tarkeshwar Singh
#   PhD, IIT Delhi
#   Email: tarkphysics87@gmail.com
#   
#  PURPOSE: Set all control variables for assimilation
#  It will be used by job.csh and run_lmdz.csh
# ****************************************************************************

# Set start date of assimilation
set start_year=2010
set start_month=5
set start_day=17

# Set end date of assimilation
set end_year=2010
set end_month=5
set end_day=30

# if assimilation start from initial then set .true. else set .false.
set  assim_job_start_from_init = .true.

# store DART restart file at every $restart_store_freq day
set restart_store_freq = 1

# set the LMDZ exe, ics and bcs data  names
# These should be exist in $LMDZ_DIR
set gcm_exe    = gcm_360x180x39_phylmd_para_mem.e
set ce0l_exe   = ce0l_360x180x39_phylmd_para.e
set limit_file = limit.nc_360x180x39
# requirs for creating perturbed ensemble members
set start_file = start.nc_360x180x39
set startphy_file = startphy.nc_360x180x39
#
# Set path for OUTPUT storage
set OUTPUT_DIR   = `pwd` 
# Path of DART software
set DART_DIR          = /home/cas/phd/asz118162/DART/kodiak/ 
# Path of DART initial start files
set DART_ics_DIR      = /scratch/cas/phd/asz118162/WORK/DART/monsoon_360x180x39/ENSEMBLES
# Define path where LMDZ gcm.e , ce0l.e , limit.nc and all LMDZ control input files (def files) exists
set LMDZ_DIR          = /scratch/cas/phd/asz118162/WORK/DART/monsoon_360x180x39/LMDZ_init 
# Path of obseravtion obs_seqYYYYMMDD files
set OBS_DIR           = /scratch/cas/phd/asz118162/DATA/OBS/DART/NCEP+GPS_2010_MJJAS 
