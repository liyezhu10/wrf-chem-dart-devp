%
path='/scratch/summit/mizzi/DART_OBS_DIAG';
%
%exp         = '/real_FRAPPE_AIRNOW_CONTROL/obs_diag_output.nc';
exp         = '/real_FRAPPE_RETR_AIR_CO/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_AIR_O3/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_AIR_NO2/obs_diag_output.nc';
%exp         = '/real_FRAPPE_RETR_AIR_SO2/obs_diag_output.nc';

%
fname=strcat(path,exp);
%
npar=1;
copystring    = {'totalspread'};
%copystring    = {'spread'};
nvar=1;
obsname      = {'AIRNOW_CO'};
%obsname      = {'AIRNOW_O3'};
%obsname      = {'AIRNOW_NO2'};
%obsname      = {'AIRNOW_SO2'};
lbnd=0.;
ubnd=0.5;
ubnd=4.0;
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
