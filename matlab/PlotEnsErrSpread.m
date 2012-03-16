function PlotEnsErrSpread( pinfo )
%% PlotEnsErrSpread     Creates summary plots of error and spread 
%
% PlotEnsErrSpread is intended to be called by 'plot_ens_err_spread'.
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: EnsErrSpread( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models 
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
% var             name of netCDF variable of interest
% var_inds        indices of variables of interest 
%
% Example 0   (9var  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 1 2 3 4 5 6 7 8 9 ];
% PlotEnsErrSpread(pinfo)
%
% Example 1   (Lorenz_96  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'state';
% pinfo.var_inds   = [ 3 4 36 39 22 ];
% PlotEnsErrSpread(pinfo)
%
% Example 2   (Lorenz_96_2scale  model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'X';
% pinfo.var_inds   = [ 3 18 27 ];
% PlotEnsErrSpread(pinfo)
%
% Example 3 (FMS BGrid model)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotEnsErrSpread(pinfo)

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

pinfo  = CheckModelCompatibility(pinfo);
Ntimes = pinfo.truth_time(2);
tend   = pinfo.truth_time(1) + Ntimes - 1;
times  = pinfo.dates(pinfo.truth_time(1):tend);

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

switch lower(pinfo.model)

   case '9var'

      % Get the appropriate copies
      truth      = get_state_copy(pinfo.truth_file, pinfo.var, truth_index, ...
                                  pinfo.truth_time(1), pinfo.truth_time(2)) ;
      ens_mean   = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
      ens_spread = get_state_copy(pinfo.diagn_file, pinfo.var, ens_spread_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3

            ivar = (i - 1)*3 + j;

            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));  
            errTotal    = sum(err)/Ntimes;
            spreadTotal = sum(ens_spread(:,ivar))/Ntimes;
            string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
            string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

            fprintf('%s model Variable %d\n',pinfo.model,ivar)

            subplot(3, 1, j);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread', pinfo.model, ivar);
               title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',Ntimes))
               ylabel('distance')
         end
      end

   case {'lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
	 'lorenz_04','forced_lorenz_96','ikeda','simple_advection'} 
      % Get the appropriate copies
      truth      = get_state_copy(pinfo.truth_file, pinfo.var, truth_index, ...
                                  pinfo.truth_time(1), pinfo.truth_time(2)) ;
      ens_mean   = get_state_copy(pinfo.diagn_file, pinfo.var, ens_mean_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;
      ens_spread = get_state_copy(pinfo.diagn_file, pinfo.var, ens_spread_index, ...
                                  pinfo.diagn_time(1), pinfo.diagn_time(2)) ;

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
            iplot = iplot + 1;
            err         = total_err(ens_mean(:,ivar) , truth(:,ivar));  
            errTotal    = sum(err)/Ntimes;
            spreadTotal = sum(ens_spread(:,ivar))/Ntimes;
            string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
            string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

            subplot(length(pinfo.var_inds), 1, iplot);
               plot(times,err, 'b', ...
                    times,ens_spread(:, ivar), 'r');
               s1 = sprintf('%s model Var %d Ensemble Error Spread', pinfo.model, ivar);
               title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
               legend(string1,string2,0)
               legend boxoff
               xlabel(sprintf('model time (%d timesteps)',Ntimes))
               ylabel('distance')
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','cam'}

      clf;

      truth      = GetCopy('fname',pinfo.truth_file, 'varname', pinfo.var, ...
                      'copyindex', truth_index, 'levelindex', pinfo.levelindex, ...
                      'latindex', pinfo.latindex, 'lonindex', pinfo.lonindex, ...
                      'tindex1', pinfo.truth_time(1), 'tcount',pinfo.truth_time(2)) ;

      ens_mean   = GetCopy('fname',pinfo.diagn_file, 'varname', pinfo.var, ...
                      'copyindex', ens_mean_index, 'levelindex', pinfo.levelindex, ...
                      'latindex', pinfo.latindex, 'lonindex', pinfo.lonindex, ...
                      'tindex1', pinfo.diagn_time(1), 'tcount',pinfo.diagn_time(2)) ;

      ens_spread = GetCopy('fname',pinfo.diagn_file, 'varname', pinfo.var, ...
                      'copyindex', ens_spread_index, 'levelindex', pinfo.levelindex, ...
                      'latindex', pinfo.latindex, 'lonindex', pinfo.lonindex, ...
                      'tindex1', pinfo.diagn_time(1), 'tcount',pinfo.diagn_time(2)) ;

      subplot(2,1,1)
         PlotLocator(pinfo);

      subplot(2,1,2)
         err         = total_err(ens_mean, truth);  
         errTotal    = sum(err)/Ntimes;
         spreadTotal = sum(ens_spread)/Ntimes;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

         plot(times,err, 'b', times,ens_spread, 'r');

         s1 = sprintf('Ensemble Mean Error, Ensemble Spread %s ''%s''',pinfo.model,pinfo.var);
         s2 = sprintf('level %d lat %.2f lon %.2f', ...
                       pinfo.level, pinfo.latitude, pinfo.longitude);
         title({s1, s2, pinfo.diagn_file},'interpreter','none','fontweight','bold');

         legend(string1,string2,0);
         legend boxoff
         xlabel(sprintf('model time (%d timesteps)',Ntimes));
         ylabel('distance');

   case {'mpas_atm'}

      clf;

      truth      = GetCopy('fname',pinfo.truth_file, 'varname', pinfo.var, ...
                      'copyindex', truth_index, 'levelindex', pinfo.levelindex, ...
                      'cellindex', pinfo.cellindex, ...
                      'tindex1', pinfo.truth_time(1), 'tcount',pinfo.truth_time(2)) ;

      ens_mean   = GetCopy('fname',pinfo.diagn_file, 'varname', pinfo.var, ...
                      'copyindex', ens_mean_index, 'levelindex', pinfo.levelindex, ...
                      'cellindex', pinfo.cellindex,  ...
                      'tindex1', pinfo.diagn_time(1), 'tcount',pinfo.diagn_time(2)) ;

      ens_spread = GetCopy('fname',pinfo.diagn_file, 'varname', pinfo.var, ...
                      'copyindex', ens_spread_index, 'levelindex', pinfo.levelindex, ...
                      'cellindex', pinfo.cellindex, ...
                      'tindex1', pinfo.diagn_time(1), 'tcount',pinfo.diagn_time(2)) ;

      subplot(2,1,1)
         PlotLocator(pinfo);

      subplot(2,1,2)
         err         = total_err(ens_mean, truth);  
         errTotal    = sum(err)/Ntimes;
         spreadTotal = sum(ens_spread)/Ntimes;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread = ' num2str(spreadTotal)];

         plot(times,err, 'b', times,ens_spread, 'r');

         s1 = sprintf('Ensemble Mean Error, Ensemble Spread %s ''%s''',pinfo.model,pinfo.var);
         s2 = sprintf('level number %d lat %.2f lon %.2f', ...
                       pinfo.level, pinfo.latCell(pinfo.cellindex), pinfo.lonCell(pinfo.cellindex));
         title({s1, s2, pinfo.diagn_file},'interpreter','none','fontweight','bold');

         legend(string1,string2,0);
         legend boxoff
         xlabel(sprintf('model time (%d timesteps)',Ntimes));
         ylabel('distance');

   otherwise
      error('model %s unknown.',pinfo.model)
end

%======================================================================
% Subfunctions
%======================================================================



function var = GetCopy(varargin)
% Gets a time-series of a single specified copy of a prognostic variable
% at a particular 3D location (level, lat, lon)

for i = 1:2:nargin,
   eval(sprintf('pinfo.%s = varargin{i+1};',varargin{i}))
end

[start, count] = GetNCindices(pinfo,'fname',pinfo.varname);
var            = nc_varget(pinfo.fname, pinfo.varname, start, count);



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90])
   continents;
   axis image
   grid on
