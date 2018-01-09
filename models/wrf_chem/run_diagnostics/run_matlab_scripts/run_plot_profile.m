%
fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_CNTL/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_CO/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_O3/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_MOP_CO/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_IAS_CO/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_IAS_O3/obs_diag_output.nc';
%
copystring    = 'ens_mean';
%copystring    = 'observation';
%copystring    = 'bias';
%copystring    = 'rmse';
%obsnamevar     = 'AIRNOW_CO';
%obsnamevar     = 'AIRNOW_O3';
obsnamevar     = 'MOPITT_CO_RETRIEVAL';
%obsnamevar     = 'IASI_CO_RETRIEVAL';
%
  plot = plot_profile(fname,copystring,'obsname',obsnamevar);

