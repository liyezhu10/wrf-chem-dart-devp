#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
# 
#***********************************************************************************************
#  Tarkeshwar Singh
#  PhD, IIT Delhi
#  Email: tarkphysics87@gmail.com
#   
#
# PORPOSE    : Run ncl scripts over each Posterior.nc/Prior.nc files resulting from sequential 
#              assimilations to extract
#
#              (1) Mean state vectors (T,U,V,Q,ps,CDLIQ)  on A-grid and press levs from C-grid and model levs
#              (2) State Vectors inflation on A-grid and model levs
#              (3) State Vectors  spread on A-grid and press levs from C-grid and model levs

# requires NCL & CDO
# INPUTS FILES: List of Prior_Diag.nc/Posterior_Diag resulting from sequential restart assimilations

#*************************************************************************************************
# List of Posterior.nc/Prior.nc files
set files_name_list = flist_posterior 

# Path of ncl scripts extract_state_vectors_ens_mean_ml_to_pl.ncl , extract_state_vectors_spread.ncl
# and extract_state_vectors_spread.ncl
set ncl_script_path = /home/cas/phd/asz118162/DART/kodiak/models/LMDZ/ncl_scripts 

# Suffix name to add on output file
set suffix          = Posterior

# Inflation index in Prior/Posterior_Diag.nc files, counted from zero
set inf_index       = 2    

#*************************************************************************************************
rm State_vector_spread_ml_*.nc State_vector_inf_ml_*.nc State_vector_ens_mean_pl_*.nc 

sed -e s/copy=3/copy=$inf_index/g $ncl_script_path/extract_state_vector_inf.ncl > extract_state_vector_inf.ncl  
#----------
set n = 1

foreach file (`cat $files_name_list`) 

    ln -sf $file input.nc

    echo "***************************************************************************************"
    echo $file
    echo "***************************************************************************************"

    echo "************ Extracting State_vector ens mean on press levs  and A-grid ************"
    ncl $ncl_script_path/extract_state_vectors_ens_mean_ml_to_pl.ncl
    mv State_vector_ens_mean_pl.nc State_vector_ens_mean_pl_$n.nc

    echo "************ Extracting State_vector inflation on model levs  and A-grid ************ "
    ncl extract_state_vectors_inf.ncl
    mv State_vector_inf_ml.nc State_vector_inf_ml_$n.nc 

    echo "************ Extracting State_vector spread on model levs  and A-grid  ************ "
    ncl $ncl_script_path/extract_state_vectors_spread_ml_to_pl.ncl
    mv State_vector_spread_ml.nc State_vector_spread_pl_$n.nc

    @ n++
end
#************************************************************************************************

cdo mergetime  State_vector_ens_mean_pl_*.nc ${suffix}_State_vec_ens_mean_pl.nc
cdo mergetime  State_vector_inf_ml_*.nc      ${suffix}_State_vec_inf_ml.nc
cdo mergetime  State_vector_spread_pl_*.nc   ${suffix}_State_vec_spread_pl.nc

rm State_vector_spread_pl_*.nc State_vector_inf_ml_*.nc State_vector_ens_mean_pl_*.nc 
