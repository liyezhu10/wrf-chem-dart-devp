fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0/obs_epoch_005.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0_loc0p5/obs_epoch_005.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0_ph-loc100/obs_epoch_005.nc';
region        = [0 360 -90 90 -Inf Inf];
ObsTypeString = 'MOPITT_CO_RETRIEVAL';
%ObsTypeString = 'IASI_CO_RETRIEVAL';
%ObsTypeString = 'IASI_O3_RETRIEVAL';
CopyString    = 'NCEP BUFR observation';
QCString      = 'DART quality control';
maxgoodQC     = 2;
verbose       = 1;
twoup         = 0;
plot          = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                      QCString, maxgoodQC, verbose, twoup);
