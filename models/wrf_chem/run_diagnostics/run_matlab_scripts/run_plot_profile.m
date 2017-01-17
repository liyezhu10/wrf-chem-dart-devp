%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/IASnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0/obs_diag_output.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0_loc0p5/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0_ph-loc100/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/XXXnIAS_Exp_2_MgDA_20M_100km_XXnO3_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0_indep/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0_no_ia_co/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/XXXnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_CNTL_VARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_CNTL_NVARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_COnXX_VARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_COnXX_NVARLOC/obs_diag_output.nc';
%
copystring1    = 'ens_mean';
copystring2    = 'observation';
copystring3    = 'bias';
copystring4    = 'rmse';
%
%plot1 = plot_profile(fname,copystring1);
%plot2 = plot_profile(fname,copystring2);
%plot3 = plot_profile(fname,copystring3);
plot4 = plot_profile(fname,copystring4);
