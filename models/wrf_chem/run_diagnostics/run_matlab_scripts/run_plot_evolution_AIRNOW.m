%
fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_CONTROL/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_CO/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_O3/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_NO2/obs_diag_output.nc';
%fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_AIR_SO2/obs_diag_output.nc';
%
npar=2;
copystring    = {'rmse','totalspread'};
%
nvar=1;
%varname      = {'AIRNOW_CO'};
varname      = {'AIRNOW_O3'};
%varname      = {'AIRNOW_NO2'};
%varname      = {'AIRNOW_SO2'};
for ipar=1:npar
   for ivar=1:nvar
%     plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar},'range',[lbnd,ubnd]);
      plot = plot_evolution(fname,copystring{ipar},'obsname',varname{ivar});
   end
end
