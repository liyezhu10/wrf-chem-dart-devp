41;326;0c%
path='/glade2/scratch2/mizzi/DART_OBS_DIAG';
%
%exp         = '/real_FRAPPE_RETR_CNTL/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_AIR_CO/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_AIR_O3/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_MOP_CO/obs_diag_output.nc';
exp         = '/real_FRAPPE_RETR_IAS_CO/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_IAS_O3/obs_diag_output.nc';
fname=strcat(path,exp);
npar=1;
copystring    = {'totalspread'};
%copystring    = {'spread'};
nvar=1;
%obsname      = {'AIRNOW_CO'};
%obsname      = {'AIRNOW_O3'};
%obsname      = {'MOPITT_CO_RETRIEVAL'};
obsname      = {'IASI_CO_RETRIEVAL'};
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
