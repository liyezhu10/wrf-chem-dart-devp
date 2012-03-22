function PlotBins(pinfo)
%% PlotBins Plots rank histograms of ensemble mean
%
% PlotBins is intended to be called by 'plot_bins'.
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotBins(pinfo)
% 
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copy tagged 'ensemble mean'
% var_inds        indices of state variables of interest
%
% Example 1 (Lorenz_96  model)
%%-------------------------------------------------------- 
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var_inds   = [3 4 36 39 22];
% PlotBins( pinfo );
%
% Example 2 (FMS BGrid model)
%%-------------------------------------------------------- 
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Prior_Diag.nc';
% pinfo.var        = 'u';
% pinfo.level      = 3;
% pinfo.latitude   = 23.5;
% pinfo.longitude  = 45.67;
% PlotBins( pinfo );

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% Get the state for the truth
truth_index = get_copy_index(pinfo.truth_file,'true state');

switch lower(pinfo.model)

   case '9var'

      % Use three different figures with three subplots each
      for i = 1:3
         figure(i); clf
         for j = 1:3
            ivar = (i - 1)*3 + j;
            truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                             ivar, pinfo.truth_time(1), pinfo.truth_time(2));
            ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
            bins  = rank_hist(ens, truth);
            subplot(3, 1, j);
            bar(bins);
            title(sprintf('%s Variable %d for %s', ...
                  pinfo.model,ivar,pinfo.diagn_file), ...
                  'interpreter','none','fontweight','bold')
            xlabel('rank')
            ylabel('occurrence')
            ax = axis;
            ax(1) = 0.5;
            ax(2) = length(bins)+0.5;
            axis(ax)
            axis tight
         end
      end

   case {'lorenz_63','lorenz_84','lorenz_96','lorenz_04','forced_lorenz_96','ikeda'}

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
         iplot = iplot + 1;
         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %d for %s', ...
               pinfo.model,ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
         xlabel('rank')
         ylabel('occurrence')
         ax = axis;
         ax(1) = 0.5;
         ax(2) = length(bins)+0.5;
         axis(ax)
         axis tight
      end

   case {'lorenz_96_2scale','simple_advection'}

      clf; iplot = 0;
      for ivar = pinfo.var_inds,
         iplot = iplot + 1;

         truth = get_var_series(pinfo.truth_file, pinfo.var, truth_index, ...
                               ivar, pinfo.truth_time(1), pinfo.truth_time(2));
         ens   = get_ens_series(pinfo.diagn_file, pinfo.var, ivar, ...
                                   pinfo.diagn_time(1), pinfo.diagn_time(2));
         bins  = rank_hist(ens, truth);
         subplot(length(pinfo.var_inds), 1, iplot);
         bar(bins);
         title(sprintf('%s Variable %s %d for %s', ...
               pinfo.model,pinfo.var, ivar,pinfo.diagn_file), ...
               'interpreter','none','fontweight','bold')
         xlabel('rank')
         ylabel('occurrence')
         ax = axis;
         ax(1) = 0.5;
         ax(2) = length(bins)+0.5;
         axis(ax)
         axis tight
      end

   case {'fms_bgrid','pe2lyr','mitgcm_ocean','cam','mpas_atm'}

      % It is intended that all 3D models have all the required information
      % set in the corresponding Get<model>Info.m script.

      clf;

      truth = GetCopy(pinfo);
      ens   = GetEns( pinfo);

      subplot(2,1,1)
         PlotLocator(pinfo)

      subplot(2,1,2)
      bins  = rank_hist(ens, truth);
      bar(bins);
      title({ ...
        sprintf('%s ''%s'' for %s ', pinfo.model, pinfo.var, pinfo.diagn_file), ...
        sprintf('level %d lat %.2f lon %.2f',pinfo.level, pinfo.latitude, ...
                 pinfo.longitude)}, 'interpreter','none','fontweight','bold')
      xlabel('rank')
      ylabel('occurrence')
      ax = axis;
      ax(1) = 0.5;
      ax(2) = length(bins)+0.5;
      axis(ax)
      axis tight

   otherwise

      error('model %s unknown',pinfo.model)

end


%======================================================================
% Subfunctions
%======================================================================


function var = GetCopy(pinfo)
% Gets a time-series of a single specified 'true' copy of a prognostic variable 
% at a particular 3D location (level, lat, lon)

pinfo.tindex1  = pinfo.truth_time(1);
pinfo.tcount   = pinfo.truth_time(2);
[start, count] = GetNCindices(pinfo,'truth',pinfo.var);

var = nc_varget(pinfo.truth_file, pinfo.var, start, count);



function var = GetEns(pinfo)
% Gets a time-series of all copies of a prognostic variable 
% at a particular location (level, gridcell).

pinfo.tindex1    = pinfo.diagn_time(1);
pinfo.tcount     = pinfo.diagn_time(2);
pinfo.copyindex1 = pinfo.ensemble_indices(1);
pinfo.copycount  = length(pinfo.ensemble_indices);
[start, count]   = GetNCindices(pinfo,'diagn',pinfo.var);

var = nc_varget(pinfo.diagn_file, pinfo.var, start, count);
% All the singleton dimensions get squeezed out.
% If there is only one time, 'bob' is copy-X-time,
% if there is more than one time, bob is time-X-copy
if (pinfo.diagn_time(2) == 1)
    var = var';
end



function PlotLocator(pinfo)
   plot(pinfo.longitude,pinfo.latitude,'pb','MarkerSize',12,'MarkerFaceColor','b');
   axis([0 360 -90 90]);
   continents;
   
