%
fname         = '/glade2/scratch2/mizzi/DART_OBS_DIAG/real_FRAPPE_RETR_MOP_CO/obs_diag_output.nc';
%
npar=1;
%copystring    = {'rmse'};
%copystring    = {'totalspread'};
copystring    = {'spread'};
%
nvar=1;
varname      = {'MOPITT_CO_RETRIEVAL'};
for ipar=1:npar
   for ivar=1:nvar
%     plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar},'range',[lbnd,ubnd]);
      plot = plot_evolution(fname,copystring{ipar},'varname',varname{ivar});
   end
end
