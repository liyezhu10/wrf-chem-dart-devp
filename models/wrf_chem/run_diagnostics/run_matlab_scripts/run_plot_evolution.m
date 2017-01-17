%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/XXXnIAS_Exp_2_MgDA_20M_100km_XXnO3_p10p30_f1p0/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0_indep/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0_no_ia_co/obs_diag_output.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/XXXnIAS_Exp_2_MgDA_20M_100km_COnO3_p10p30_f1p0/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_CNTL_VARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_CNTL_NVARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_COnXX_VARLOC/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/real_FRAPPE_COnXX_NVARLOC/obs_diag_output.nc';
%
npar=1;
%copystring    = {'rmse'};
copystring    = {'totalspread'};
copystring    = {'spread'};
%
nvar=1;
%varname      = {'MOPITT_CO_RETRIEVAL'};
varname      = {'IASI_CO_RETRIEVAL'};
%varname      = {'IASI_O3_RETRIEVAL'};
for ipar=1:npar
for ivar=1:nvar
%	   plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar},'range',[lbnd,ubnd]);
	   plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar});
end
end
