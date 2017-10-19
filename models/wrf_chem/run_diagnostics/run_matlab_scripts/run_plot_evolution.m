%
fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_MOP_CO/obs_diag_output.nc';
%
npar=2;
copystring    = {'rmse','totalspread'};
%
nvar=3;
varname      = {'MOPITT_CO_RETRIEVAL','IASI_O3_RETRIEVAL','AIRNOW_O3'};
for ipar=1:npar
   for ivar=1:nvar
%     plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar},'range',[lbnd,ubnd]);
      plot = plot_evolution(fname,copystring{ipar},'obsname',varname{ivar});
   end
end
