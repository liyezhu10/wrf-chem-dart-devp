function PlotTotalErr( pinfo )
%% PlotTotalErr Plots summary plots of global error and spread
%
% The error is a measure of the distance between the 
% 'ensemble mean' and the 'true_state'. The distance of the 
% spread is also plotted.
%
% PlotTotalErr is intended to be called by 'plot_total_err'
% The only input argument is a structure with model-dependent
% components.
%
% USAGE: PlotTotalErr( pinfo )
%
% STRUCTURE COMPONENTS FOR low-order models
% truth_file      name of netCDF DART file with copy tagged 'true state'
% diagn_file      name of netCDF DART file with copies tagged 'ensemble mean'
%                 and 'ensemble spread'
%
% Example 1   (Lorenz_63, Lorenz_96, 9var ...)
%%--------------------------------------------------------
% pinfo.truth_file = 'True_State.nc';
% pinfo.diagn_file = 'Posterior_Diag.nc';
% PlotTotalErr( pinfo )

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

% adds the time overlap info to the structure.  
pstruct = CheckModelCompatibility(pinfo);
pinfo   = CombineStructs(pinfo, pstruct) 

model   = nc_attget(pinfo.truth_file, nc_global, 'model');

% Get the netcdf variable indices for desired "copies"
% Get the indices for the true state, ensemble mean and spread
% The metadata is queried to determine which "copy" is appropriate.
truth_index      = get_copy_index(pinfo.truth_file, 'true state' );
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% The times are set before from the model compatibility call
num_times = pinfo.truth_time(2);
times     = nc_varget(pinfo.truth_file,'time', pinfo.truth_time(1)-1, num_times);

switch lower(model)

   case {'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_04','ikeda'}

      %% Get the appropriate netcdf variables
      truth  = get_state_copy(pinfo.truth_file, 'state',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      ens    = get_state_copy(pinfo.diagn_file, 'state',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      spread = get_state_copy(pinfo.diagn_file, 'state',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      num_vars = size(spread,2);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      s1 = sprintf('%s Total Error over all %d variables', model, num_vars);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case 'lorenz_96_2scale'

      %% Simply going to append X,Y together and treat as above.

      % Get the appropriate netcdf variables
      tim    = get_state_copy(pinfo.truth_file, 'X',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      tom    = get_state_copy(pinfo.truth_file, 'Y',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      truth  = [tim tom];
      tim    = get_state_copy(pinfo.diagn_file, 'X',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      tom    = get_state_copy(pinfo.diagn_file, 'Y',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      ens    = [tim tom]; 
      tim    = get_state_copy(pinfo.diagn_file, 'X',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      tom    = get_state_copy(pinfo.diagn_file, 'Y',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      spread = [tim tom]; clear tim tom
      num_vars = size(spread,2);

      % Also need to compute the spread; zero truth for this and
      % compute distance from 0
      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf;
      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      s1 = sprintf('%s Total Error over all %d variables',model,num_vars);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case 'forced_lorenz_96'

      %% This model has the state variables replicated, so there is a difference
      % between num_state_vars and the length of the state variable.
      %forcing           = nc_attget(pinfo.truth_file, nc_global, 'model_forcing');
      %delta_t           = nc_attget(pinfo.truth_file, nc_global, 'model_delta_t');
      %time_step_days    = nc_attget(pinfo.truth_file, nc_global, 'model_time_step_days');
      %time_step_seconds = nc_attget(pinfo.truth_file, nc_global, 'model_time_step_seconds');
      num_model_vars    = nc_attget(pinfo.truth_file, nc_global, 'model_num_state_vars');

      % Get the appropriate netcdf variables

      Whole_truth  = get_state_copy(pinfo.truth_file, 'state',     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
      Whole_ens    = get_state_copy(pinfo.diagn_file, 'state',  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      Whole_spread = get_state_copy(pinfo.diagn_file, 'state',ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
      num_vars = size(Whole_spread,2);

      %--------------------------------------------------------------------------
      % Treat the traditional state variable independent of the forcing variables
      %--------------------------------------------------------------------------

      ind1 = 1;                 % ASSUMPTION: traditional state is first
      indN = num_model_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      clf; subplot(2,1,1);

      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      s1 = sprintf('%s Total Error over statevars %d to %d', model, ind1, indN);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

      %--------------------------------------------------------------------------
      % Now for the forcing
      %--------------------------------------------------------------------------

      ind1 = num_model_vars + 1;
      indN = num_vars;

      truth  = Whole_truth(  :, ind1:indN );
      ens    = Whole_ens(    :, ind1:indN );
      spread = Whole_spread( :, ind1:indN );

      % Compute the spread; zero truth for this and compute distance from 0

      err        = total_err(truth, ens);
      err_spread = total_err(zeros(size(spread)), spread);
      errTotal   = sum(err)/num_times;
      spreadTotal= sum(err_spread)/num_times;
      string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
      string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

      subplot(2,1,2)

      plot(times,err, 'b', times,err_spread, 'r');
      legend(string1,string2,0)
      legend boxoff
      s1 = sprintf('%s Total Error over statevars %d to %d', model, ind1, indN);
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')
      xlabel(sprintf('model time (%d timesteps)',num_times))
      ylabel('Total Error')

   case {'simple_advection'}

      %% if the 'state' variable exists ... then
      % 'concentration','source', and 'wind' do not.

      if ( nc_isvar(pinfo.truth_file, 'state') )
         varlist = {'state'};
      else
         varlist = {'concentration','source','wind','mean_source','source_phase'};
      end

      for ivar = 1:length(varlist)

         % Get the appropriate netcdf variables
         truth  = get_state_copy(pinfo.truth_file, varlist{ivar},     truth_index, ...
                              pinfo.truth_time(1), pinfo.truth_time(2));
         ens    = get_state_copy(pinfo.diagn_file, varlist{ivar},  ens_mean_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
         spread = get_state_copy(pinfo.diagn_file, varlist{ivar},ens_spread_index, ...
                              pinfo.diagn_time(1), pinfo.diagn_time(2));
         num_vars = size(spread,2);

         % Also need to compute the spread; zero truth for this and
         % compute distance from 0
         err        = total_err(truth, ens);
         err_spread = total_err(zeros(size(spread)), spread);
         errTotal   = sum(err)/num_times;
         spreadTotal= sum(err_spread)/num_times;
         string1 = ['time-mean Ensemble Mean Total Error = ' num2str(errTotal)];
         string2 = ['time-mean Ensemble Spread Total Error = ' num2str(spreadTotal)];

         figure(ivar); clf(ivar);
         plot(times,err, 'b', times,err_spread, 'r');
         legend(string1,string2,0)
         legend boxoff
         string1 = sprintf('%s Total Error over all %d variables', model, num_vars);
         string2 = sprintf('''%s'' %s', varlist{ivar}, pinfo.diagn_file);
         title({string1,string2},'interpreter','none','fontweight','bold')
         xlabel(sprintf('model time (%d timesteps)',num_times))
         ylabel('Total Error')

      end

   case 'fms_bgrid'
       %% GFDL bgrid model

      BgridTotalError( pinfo )

   case 'pe2lyr'

       %% primitive equation 2 layer model
       
      Pe2lyrTotalError( pinfo )

   case 'pbl_1d'
       %% planetary boundary layer column model

      PBL1DTotalError( pinfo )

   case 'mitgcm_ocean'
       %% MIT general circulation ocean model

      MITGCMOceanTotalError( pinfo )

   case 'mpas_atm'
      %% unstructured grid atmosphere model
       
      MPAS_ATMTotalError( pinfo )

   otherwise

      fprintf('unsupported model %s -- doing nothing\n',model)

end

%=======================================================================
% helper functions
%=======================================================================

function PBL1DTotalError ( pinfo )

%% Get some standard plotting arrays
num_times = pinfo.truth_time(2);
% z_level  = nc_varget(pinfo.truth_file, 'z_level');
%sl_level  = nc_varget(pinfo.truth_file,'sl_level');
times     = nc_varget(pinfo.truth_file, 'time',pinfo.truth_time(1)-1, num_times);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% U variable

count = [num_times 1];

start  = [pinfo.truth_time(1)-1  truth_index-1];
truth  = nc_varget(pinfo.truth_file,'U', start, count);

start  = [pinfo.diagn_time(1)-1  ens_mean_index-1];
ens    = nc_varget(pinfo.diagn_file,'U', start, count);

start  = [pinfo.diagn_time(1)-1  ens_spread_index-1];
spread = nc_varget(pinfo.diagn_file,'U', start, count);

err        = total_err(              truth,    ens);
err_spread = total_err(zeros(size(spread)), spread);

y_error  = squeeze(mean(err,2));           % mean over all levels
y_spread = squeeze(mean(err_spread,2));    % mean over all levels

plot(times,y_error,'r-',times,y_spread,'g-')
title('PBL_1d mean error of U over time ... all levels.')
xlabel('days')
ylabel('mean (all levels) total error')
axis([-Inf Inf 0 Inf])



function BgridTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

model     = nc_attget(pinfo.truth_file, nc_global, 'model');
timeunits = nc_attget(pinfo.truth_file, 'time',    'units');

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
tlons     = nc_varget(pinfo.truth_file, 'TmpI'); num_tlons  = length(tlons );
tlats     = nc_varget(pinfo.truth_file, 'TmpJ'); num_tlats  = length(tlats );
vlons     = nc_varget(pinfo.truth_file, 'VelI'); num_vlons  = length(vlons );
vlats     = nc_varget(pinfo.truth_file, 'VelJ'); num_vlats  = length(vlats );
levels    = nc_varget(pinfo.truth_file, 'lev' ); num_levels = length(levels);
times     = nc_varget(pinfo.truth_file, 'time', tstart-1, tcount); 
num_times = length(times);

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_levels);
spread_final = zeros(num_times, pinfo.num_state_vars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = reshape(SphereWeights(tlons, tlats),1,num_tlats*num_tlons); % Temperature Grid
vwts = reshape(SphereWeights(vlons, vlats),1,num_vlats*num_vlons); % Velocity    Grid

% Can we afford to get the whole thing at once ???
%----------------------------------------------------------------------
% surface pressure has only one level.
% Get2D   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing surface pressure ...')

ivar   = 1;
ilevel = 1;

truth      = Get2D(pinfo.truth_file, 'ps',      truth_index, tstart, tcount);
ens        = Get2D(pinfo.diagn_file, 'ps',   ens_mean_index, dstart, dcount);
spread     = Get2D(pinfo.diagn_file, 'ps', ens_spread_index, dstart, dcount);

err        = total_err(              truth,    ens, twts );
err_spread = total_err(zeros(size(spread)), spread, twts );

         rms(:,ivar,ilevel) = err;         % spatial mean
spread_final(:,ivar,ilevel) = err_spread;  % spatial mean

clear truth ens spread err err_spread  

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels
for ivar=2:pinfo.num_state_vars,

   fprintf('Processing level %d of %d ...\n',ilevel,num_levels)
   %-------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tcount);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dcount);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dcount);

   switch lower(pinfo.vars{ivar})
      case {'t'}
         err        = total_err(              truth,    ens, twts);
         err_spread = total_err(zeros(size(spread)), spread, twts);
      otherwise
         err        = total_err(              truth,    ens, vwts);
         err_spread = total_err(zeros(size(spread)), spread, vwts);
   end

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;
 
      varunits = nc_attget(pinfo.truth_file, pinfo.vars{ivar}, 'units');

      s1 = sprintf('%s %s Ensemble Mean', model,pinfo.vars{ivar});

      switch lower(pinfo.vars{ivar})
         case {'ps'}
            plot(times,          rms(:, ivar, 1), '-'); hold on;
            plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_levels,
               s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      h = legend(s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end



function MITGCMOceanTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has only prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

model     = nc_attget(pinfo.truth_file, nc_global, 'model');
timeunits = nc_attget(pinfo.truth_file, 'time',    'units');

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
XC        = nc_varget(pinfo.truth_file, 'XC'); num_XC  = length(XC );
YC        = nc_varget(pinfo.truth_file, 'YC'); num_YC  = length(YC );
XG        = nc_varget(pinfo.truth_file, 'XG');
YG        = nc_varget(pinfo.truth_file, 'YG');
ZC        = nc_varget(pinfo.truth_file, 'ZC'); num_ZC  = length(ZC);
times     = nc_varget(pinfo.truth_file, 'time', tstart-1, tcount); 
num_times = length(times);

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_ZC);
spread_final = zeros(num_times, pinfo.num_state_vars, num_ZC);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = reshape(SphereWeights(XC, YC),1,num_YC*num_XC);   % Grid Centers
vwts = reshape(SphereWeights(XG, YG),1,num_YC*num_XC);   % Grid edGes

%----------------------------------------------------------------------
% temperature ...  num_times x num_levels x num_lats x num_lons
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_ZC,     % Loop through all levels
for ivar=1:pinfo.num_state_vars-1,  % knowing ssh is last

   fprintf('Processing level %d of %d ...\n',ilevel,num_ZC)
   %-------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tcount);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dcount);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dcount);

   landmask = find(isfinite(truth(1,:))); % presume all members have same mask
   vweights = vwts(landmask);
   tweights = twts(landmask);

   mytruth  =  truth(:,landmask);
   myens    =    ens(:,landmask);
   myspread = spread(:,landmask);

   switch lower(pinfo.vars{ivar})
      case {'v'}
         err        = total_err(              mytruth,    myens, vweights);
         err_spread = total_err(zeros(size(myspread)), myspread, vweights);
      otherwise
         err        = total_err(              mytruth,    myens, tweights);
         err_spread = total_err(zeros(size(myspread)), myspread, tweights);
   end

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%----------------------------------------------------------------------
% surface pressure has only one level.
% Get2D   returns a   num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

disp('Processing sea surface height ...')

ivar   = find(strcmp(pinfo.vars,'SSH'));
ilevel = 1;

truth      = Get2D(pinfo.truth_file, 'SSH',      truth_index, tstart, tcount);
ens        = Get2D(pinfo.diagn_file, 'SSH',   ens_mean_index, dstart, dcount);
spread     = Get2D(pinfo.diagn_file, 'SSH', ens_spread_index, dstart, dcount);

landmask   = find(isfinite(truth(1,:))); % presume all members have same mask
tweights   = twts(landmask);

mytruth    =  truth(:,landmask);
myens      =    ens(:,landmask);
myspread   = spread(:,landmask);

err        = total_err(              mytruth,    myens, tweights );
err_spread = total_err(zeros(size(myspread)), myspread, tweights );

         rms(:,ivar,ilevel) = err;         % spatial mean
spread_final(:,ivar,ilevel) = err_spread;  % spatial mean

clear truth ens spread err err_spread  

%----------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;
      varunits = nc_attget(pinfo.truth_file,pinfo.vars{ivar},'units');

      s1 = sprintf('%s %s Ensemble Mean', model,pinfo.vars{ivar});

      switch lower(pinfo.vars{ivar})
         case {'ssh'}
            plot(times,          rms(:, ivar, 1), '-'); hold on;
            plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_ZC,
               s{i       } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_ZC} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      h = legend(s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end



function Pe2lyrTotalError( pinfo )
%% -------------------------------------------------------------------
% netcdf has no state vector, it has prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

model     = nc_attget(pinfo.truth_file, nc_global, 'model');
timeunits = nc_attget(pinfo.truth_file, 'time',    'units');

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
tlons     = nc_varget(pinfo.truth_file, 'lon' ); num_tlons  = length(tlons );
tlats     = nc_varget(pinfo.truth_file, 'lat' ); num_tlats  = length(tlats );
levels    = nc_varget(pinfo.truth_file, 'lev' ); num_levels = length(levels);
times     = nc_varget(pinfo.truth_file, 'time', tstart-1, tcount); 
num_times = length(times );

% Initialize storage for error averaging
rms          = zeros(num_times, pinfo.num_state_vars, num_levels);
spread_final = zeros(num_times, pinfo.num_state_vars, num_levels);

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Calculate weights for area-averaging.
twts = reshape(SphereWeights(tlons, tlats),1,num_tlats*num_tlons);

%----------------------------------------------------------------------
% Trying not to assume we can get the whole 3D array at once.
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ilevel = 1:num_levels,     % Loop through all levels
for ivar=1:pinfo.num_state_vars,

   fprintf('Processing level %d of %d ...\n',ilevel,num_levels)
   %%------------------------------------------------------------------
   % all vars organized    num_times x num_levels x num_lats x num_lons
   %-------------------------------------------------------------------

   truth  = GetLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tcount);
   ens    = GetLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dcount);
   spread = GetLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dcount);

   err        = total_err(              truth,    ens, twts);
   err_spread = total_err(zeros(size(spread)), spread, twts);

            rms(:, ivar, ilevel) = err;
   spread_final(:, ivar, ilevel) = err_spread;

end
end % End of level loop

clear truth ens spread err err_spread

%%---------------------------------------------------------------------
% Each variable in its own figure window
%----------------------------------------------------------------------
for ivar=1:pinfo.num_state_vars,

   figure(ivar); clf;

      varunits = nc_attget(pinfo.truth_file, pinfo.vars{ivar}, 'units');

      s1 = sprintf('%s %s Ensemble Mean for %s', model,pinfo.vars{ivar});

      switch lower(pinfo.vars{ivar})
         case {'ps'}
            plot(times,          rms(:, ivar, 1), '-'); hold on;
            plot(times, spread_final(:, ivar, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, ivar, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, ivar, 1)));
         otherwise
            plot(times, squeeze(         rms(:, ivar, :)),'-'); hold on;
            plot(times, squeeze(spread_final(:, ivar, :)),'--');

            for i = 1:num_levels,
               s{i           } = sprintf('level %d error  %.3f', i,mean(         rms(:, ivar, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, ivar, i)));
            end
      end

      h = legend(s); legend(h,'boxoff')
      grid on;
      xlabel(sprintf('time (%s) %d timesteps',timeunits,num_times))
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end


function MPAS_ATMTotalError( pinfo )
%% -------------------------------------------------------------------
% Assume netcdf has only prognostic variables.
% We are going to plot the total error (over a horizontal slice) 
% for each variable and annotate an area-weighted total.
%---------------------------------------------------------------------

model     = nc_attget(pinfo.truth_file, nc_global, 'model');
timeunits = nc_attget(pinfo.truth_file, 'time',    'units');

tstart = pinfo.truth_time(1);
tcount = pinfo.truth_time(2);
dstart = pinfo.diagn_time(1);
dcount = pinfo.diagn_time(2);

% Since the models are "compatible", get the info from either one.
area       = nc_varget(pinfo.truth_file, 'areaCell');
times      = nc_varget(pinfo.truth_file, 'time', tstart-1, tcount);
timeunits  = nc_attget(pinfo.truth_file, 'time', 'units');
timebase   = sscanf(timeunits,'%*s%*s%d%*c%d%*c%d'); % YYYY MM DD
timeorigin = datenum(timebase(1),timebase(2),timebase(3));
dates      = times + timeorigin;
num_times  = length(dates); 

% Get the indices for the true state, ensemble mean and spread                  
% The metadata is queried to determine which "copy" is appropriate.             
truth_index      = get_copy_index(pinfo.truth_file, 'true state');
ens_mean_index   = get_copy_index(pinfo.diagn_file, 'ensemble mean');
ens_spread_index = get_copy_index(pinfo.diagn_file, 'ensemble spread');

% Generate the area-based weights. Since we've got a hybrid vertical
% coordinate system, this is a simplification that is not perfect.
weights = area ./ sum(area(:));

%----------------------------------------------------------------------
% GetLevel returns a    num_times x num_lats*num_lons   2Darray.
%----------------------------------------------------------------------

for ivar=1:pinfo.num_state_vars,
    
   varunits = nc_attget(pinfo.truth_file, pinfo.vars{ivar}, 'units');
   s1 = sprintf('%s %s Ensemble Mean', model, pinfo.vars{ivar});
   fprintf('%s ...\n',s1)

   % Initialize storage for error averaging
   num_levels   = FindNLevels(pinfo.truth_file, pinfo.vars{ivar});
   rms          = zeros(num_times, num_levels);
   spread_final = zeros(num_times, num_levels);

   %% Loop over all levels for this variable
   for ilevel = 1:num_levels,

      %-------------------------------------------------------------------
      % all vars organized    num_times x num_levels x num_lats x num_lons
      %-------------------------------------------------------------------

      truth  = GetUnstructuredLevel(pinfo.truth_file, pinfo.vars{ivar},      truth_index, ilevel, tstart, tcount);
      ens    = GetUnstructuredLevel(pinfo.diagn_file, pinfo.vars{ivar},   ens_mean_index, ilevel, dstart, dcount);
      spread = GetUnstructuredLevel(pinfo.diagn_file, pinfo.vars{ivar}, ens_spread_index, ilevel, dstart, dcount);
   
      err        = total_err(              truth,    ens, weights);
      err_spread = total_err(zeros(size(spread)), spread, weights);
   
               rms(:, ilevel) = err;
      spread_final(:, ilevel) = err_spread;

   end
   
   %% Plot the time-evolution for all levels on one figure.

   figure(ivar); clf;

      switch lower(pinfo.vars{ivar})
         case {'surface_pressure'}
            plot(dates,          rms(:, 1), '-'); hold on;
            plot(dates, spread_final(:, 1), '--');

            s{1} = sprintf('time-mean Ensemble Mean error  = %f', mean(         rms(:, 1)));
            s{2} = sprintf('time-mean Ensemble Spread = %f',      mean(spread_final(:, 1)));
            
          otherwise
            plot(dates,          rms,'-'); hold on;
            plot(dates, spread_final,'--');

            for i = 1:num_levels,
               s{i       }    = sprintf('level %d error  %.3f', i,mean(         rms(:, i)));
               s{i+num_levels} = sprintf('level %d spread %.3f', i,mean(spread_final(:, i)));
            end
      end

      h = legend(s); legend(h,'boxoff')
      grid on;
      if (num_times > 6) 
         datetick('x','mm.dd.HH','keeplimits','keepticks'); % 'mm/dd'
         monstr = datestr(dates(1),31);
         xlabelstring = sprintf('month/day/HH - %s start',monstr);
      else
         datetick('x',31,'keeplimits','keepticks'); %'yyyy-mm-dd HH:MM:SS'
         monstr = datestr(dates(1),31);
         xlabelstring = sprintf('%s start',monstr);
      end
          
      xlabel(xlabelstring)
      ylabel(sprintf('global-area-weighted distance (%s)',varunits))
      title({s1,pinfo.diagn_file},'interpreter','none','fontweight','bold')

end



function slice = Get2D(fname, varname, copyindex, tstartind, tcount )
%% this gets a bit funky if tcount == 1; 
% it automatically squeezes the singleton dimension,
% leaving the dims in the wrong slots.  copyindex will always
% squeeze out, so try to figure out if you get a 3d or 2d return
% back from getnc() and then deal with time explicitly.
% same goes below in GetLevel().

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end
ted = nc_varget(fname, varname, start, count);
n   = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);




function slice = GetLevel(fname, varname, copyindex, ilevel, tstartind, tcount)
%%
myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = ilevel;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end

ted = nc_varget(fname, varname, start, count);
n   = ndims(ted);
if (n == 2)
   [ny,nx] = size(ted);
   nt = 1;
else
   [nt,ny,nx] = size(ted);
end
slice      = reshape(ted,[nt ny*nx]);



function slab = GetUnstructuredLevel(fname, varname, copyindex, ilevel, tstartind, tcount)
%% Unstructured levels are inherently (ntime, Ncell)
%  The crazy nc_varget() function will squeeze out singleton dimensions.
%  If time is a singleton, the array gets returned (Ncell, 1) which is not
%  what we want.

myinfo.diagn_file = fname;
myinfo.copyindex  = copyindex;
myinfo.levelindex = ilevel;
[start, count]    = GetNCindices(myinfo,'diagn',varname);

varinfo = nc_getvarinfo(fname,varname);

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'time'}
         start(i) = tstartind - 1;
         count(i) = tcount;
         break
      otherwise
   end
end

slab = nc_varget(fname, varname, start, count);

if (tcount == 1)
    slab = slab';
end



function wts = SphereWeights(lons, lats)
%% SphereWeights creates weights based on area ...
%
% lons,lats must be 1D arrays (in degrees)

nlons = length(lons);
nlats = length(lats);

if ( numel(lons) ~= nlons )
   disp('longitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end
if ( numel(lats) ~= nlats )
   disp('latitude array is of higher dimension than anticipated.')
   error('Must be a vector.')
end

rads    = zeros(nlats,1);               % Ensure lats is a column vector,
rads(:) = pi*lats/180.0;                % and convert to radians.
wts     = cos( rads ) * ones(1,nlons);  % Results in a [nlat-x-nlon] matrix.
wts     = wts ./ sum(wts(:));           % Normalize to unity.


function nlevels = FindNLevels(fname, varname)
%%
varinfo = nc_getvarinfo(fname,varname);
nlevels = 1; % implicitly surface-only variable

for i = 1:length(varinfo.Dimension)
   switch( lower(varinfo.Dimension{i}))
      case{'nvertlevels','nvertlevelsp1'}
         nlevels = varinfo.Size(i);
         break
      otherwise
   end
end

