function CAMTotalError( pinfo )
%% -------------------------------------------------------------------
% Plot the total area-weighted error for each variable.
%---------------------------------------------------------------------

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/mpas/matlab/CAMTotalError.m $
% $Id: CAMTotalError.m 5655 2012-04-05 23:17:16Z thoar $
% $Revision: 5655 $
% $Date: 2012-04-05 17:17:16 -0600 (Thu, 05 Apr 2012) $

% Since the models are "compatible", get the info from either one.
lons       = nc_varget(pinfo.truth_file, 'lon');
gw         = nc_varget(pinfo.truth_file, 'gw');
num_lons   = length(lons);

% make a matrix of weights for each horizontal slice
% ensure the gaussian weights sum to unity.
[~,weights] = meshgrid(ones(1,num_lons),gw);
weights     = weights/sum(weights(:));
wts         = weights(:);

% lili, interpolate from gw on unstaggered point to gw on staggered point (US)
[n1 n2] = size(gw);
if ( n1 == 1 )
    ngw = n2;
    gw_stag = zeros(1,ngw-1);
else
    ngw = n1;
    gw_stag = zeros(ngw-1,1);
end
for i = 1:ngw-1
    gw_stag(i) = (gw(i)+gw(i+1))/2.0;
end   
[~,weights_stag] = meshgrid(ones(1,num_lons),gw_stag);
weights_stag     = weights_stag/sum(weights_stag(:));
wts_stag         = weights_stag(:);

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

      if (length(size(truth)) == 2)
         nlev = 1;
      elseif (length(size(truth)) == 3)
         nlev = size(truth,3);
      else
         error('Dang, this cannot happen in CAM.')
      end

      %% Calculate the weighted mean squared error for each level.

      msqe_Z = zeros(nlev,1);
      sprd_Z = zeros(nlev,1);

      for ilevel=1:nlev,

         slabS2E   = (truth(:,:,ilevel) - ens(:,:,ilevel)).^2;  % OK even if 2D iff ilevel = 1
         if ( strcmpi(varname,'US') )
            XY_err    = sum(slabS2E(:) .* wts_stag);
         else
            XY_err    = sum(slabS2E(:) .* wts);
         end
         slabS2E   = spread(:,:,ilevel).^2;
         if ( strcmpi(varname,'US') )
            XY_spread = sum(slabS2E(:) .* wts_stag);
         else
            XY_spread = sum(slabS2E(:) .* wts);
         end

         msqe_Z(ilevel) = XY_err;
         sprd_Z(ilevel) = XY_spread;

      end % loop over levels

      %% Take the square root of the mean of all levels
      rmse(itime) = sqrt(mean(msqe_Z));
      sprd(itime) = sqrt(mean(sprd_Z));

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
