%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_RtDA_40M_p30p30_sp4/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_40M_p30p30_sp4/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_40M_p30p30_sp4/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p20p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p30p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p30p30/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_20M_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_p10p30/obs_diag_output.nc';
%
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bar_1_p10p30/obs_diag_output.nc';
fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_loc_a_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_bar_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_loc_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_3_MgDA_20M_100km_p10p00/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_p10p30/obs_diag_output.nc';
%fname         = '/glade/p/acd/mizzi/DART_OBS_DIAG/MOPCOMB_Exp_2_MgDA_20M_100km_NoRot_p10p30/obs_diag_output.nc';
%
npar=1;
copystring    = {'total spread'};
nvar=4;
varname      = {'RADIOSONDE_TEMPERATURE','RADIOSONDE_U_WIND_COMPONENT','RADIOSONDE_V_WIND_COMPONENT','MOPITT_CO_RETRIEVAL'};
%
for ipar=1:npar
for ivar=1:nvar
    plot = plot_rmse_xxx_evolution(fname,copystring{ipar},varname{ivar});
end
end



