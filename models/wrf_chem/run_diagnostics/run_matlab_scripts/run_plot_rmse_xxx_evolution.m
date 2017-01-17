path='/glade/p/acd/mizzi/DART_OBS_DIAG';
%
% 2008 case files 
%
%obs_diag_file='/obs_diag_output_MOPITT_CO.nc';
%exp='/MOPnXXX_Exp_1_MgDA_20M_100km_COnXX_RETR_NO_ROT_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RAWR_NO_ROT_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RAWR_F50_NO_ROT_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RAWR_NO_ROT_EINV_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RAWR_F50_NO_ROT_EINV_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_ME_NO_ROT_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_QOR_SCALE_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_CPSR_SCALE_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_CPSR_SCALE_BLOC_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_VL1k_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_VL2k_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_VL4k_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_VL8k_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_NO_ROT_RJ3_SUPR';
%exp='/MOPnXXX_Exp_2_MgDA_20M_100km_COnXX_RETR_CPSR_SCALE_RJ3_SUPR';
%
obs_diag_file='/obs_diag_output_MOPITT_CO.nc';
exp='/MOPnIASnMOD_Exp_2_MgDA_20M_100km_COnXXnAOD_CPSR_Joint';
%exp='/MOPnIASnMOD_Exp_2_MgDA_20M_100km_COnXXnAOD_CPSR_Indep';
%
%obs_diag_file='/obs_diag_output_IASI_CO.nc';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_NO_ROT_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_F10_NO_ROT_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_F25_NO_ROT_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_F50_NO_ROT_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_QOR_SCALE_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_CPSR_SCALE_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_F25_CPSR_SCALE_SUPR';
%exp='/XXXnIAS_Exp_2_MgDA_20M_100km_COnXX_RAWR_F50_CPSR_SCALE_SUPR';
%
% FRAPPE files
%
%obs_diag_file='/obs_diag_output_MOPITT_CO.nc';
%exp='/real_FRAPPE_CNTL_VARLOC_REV1.a';
%exp='/real_FRAPPE_CNTL_VARLOC_REV3.a';
%exp='/real_FRAPPE_COnXX_VARLOC_REV1.a';
%exp='/real_FRAPPE_COnXX_VARLOC_REV2.a';
%exp='/real_FRAPPE_COnXX_VARLOC_REV3.a';
%
fname=strcat(path,exp,obs_diag_file);
npar=1;
copystring    = {'totalspread'};
%copystring    = {'spread'};
nvar=1;
obsname      = {'MOPITT_CO_RETRIEVAL'};
%obsname      = {'IASI_CO_RETRIEVAL'};
%obsname      = {'RADIOSONDE_U_WIND_COMPONENT'};
%obsname      = {'RADIOSONDE_V_WIND_COMPONENT'};
%obsname      = {'RADIOSONDE_TEMPERATURE'};
%obsname      = {'IASI_O3_RETRIEVAL'};
lbnd=0.;
ubnd=0.4;
%ubnd=0.3;
%ubnd=1.5;
%%ubnd=3.0;
%
for ipar=1:npar
for ivar=1:nvar
plot = plot_rmse_xxx_evolution(fname,copystring{ipar},'obsname',obsname{ivar});
%plot = plot_rmse_xxx_evolution(fname,copystring{ipar},'obsname',obsname{ivar},'range',[lbnd,ubnd]);
end
end
