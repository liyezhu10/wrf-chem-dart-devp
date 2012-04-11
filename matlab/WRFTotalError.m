function WRFTotalError( pinfo )
%% -------------------------------------------------------------------
% Plot the total area-weighted error for each variable.
%---------------------------------------------------------------------

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/mpas/matlab/WRFTotalError.m $
% $Id: WRFTotalError.m 5686 2012-04-11 00:01:13Z thoar $
% $Revision: 5686 $
% $Date: 2012-04-10 18:01:13 -0600 (Tue, 10 Apr 2012) $

% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

%----------------------------------------------------------------------
%
%----------------------------------------------------------------------

for ivar=1:pinfo.num_state_vars,

   varname = pinfo.vars{ivar};

   rmse = zeros(pinfo.time_series_length,1);
   sprd = zeros(pinfo.time_series_length,1);

   for itime=1:pinfo.time_series_length,
   
      fprintf('Processing %s timestep %d of %d ...\n', ...
                varname, itime, pinfo.time_series_length)
   
      truth  = get_hyperslab('fname',pinfo.truth_file, 'varname',varname, ...
                   'copyindex',truth_index, 'timeindex',pinfo.truth_time(1)+itime-1);
      ens    = get_hyperslab('fname',pinfo.diagn_file, 'varname',varname, ...
                   'copyindex',ens_mean_index, 'timeindex',pinfo.diagn_time(1)+itime-1);
      spread = get_hyperslab('fname',pinfo.diagn_file, 'varname',varname, ...
                   'copyindex',ens_spread_index, 'timeindex',pinfo.diagn_time(1)+itime-1);
   
      %% Calculate the mean squared error for each level. By construction,
      %  the WRF grid is close enough to equal area for each grid cell.
      %  The simple mean of each level is OK.

      msqe = mean( (truth(:) - ens(:)) .^2 ); % mean squared error
      msqs = mean(       spread(:)     .^2 ); % mean squared spread

      %% Take the square root of the mean of all levels
      rmse(itime) = sqrt(msqe);
      sprd(itime) = sqrt(msqs);
      
   end % loop over time

   %-------------------------------------------------------------------
   % Each variable in its own figure window
   %-------------------------------------------------------------------
   figure(ivar); clf;
      varunits = nc_attget(pinfo.truth_file, pinfo.vars{ivar}, 'units');

      plot(pinfo.time,rmse,'-', pinfo.time,sprd,'--')

      s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(rmse));
      s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(sprd));

      h = legend(s); legend(h,'boxoff')
      grid on;
      xdates(pinfo.time)
      ylabel(sprintf('global-area-weighted rmse (%s)',varunits))
      s1 = sprintf('%s %s Ensemble Mean', pinfo.model,pinfo.vars{ivar});
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end % loop around variables 

clear truth ens spread err XY_spread


function xdates(dates)
if (length(get(gca,'XTick')) > 6)
   datetick('x','mm.dd.HH','keeplimits'); % 'mm/dd'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('month/day/HH - %s start',monstr);
else
   datetick('x',31,'keeplimits'); %'yyyy-mm-dd HH:MM:SS'
   monstr = datestr(dates(1),31);
   xlabelstring = sprintf('%s start',monstr);
end
xlabel(xlabelstring)
