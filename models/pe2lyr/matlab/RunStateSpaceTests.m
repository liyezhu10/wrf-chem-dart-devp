function RunStateSpaceTests(dummy)
%% RunStateSpaceTests.m

%------------------------------------------------------------
% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL:
% https://subversion.ucar.edu/DAReS/DART/branches/mpas/models/cam/matlab/RunStateSpaceTests.m $
% $Id$
% $Revision$
% $Date$
%------------------------------------------------------------

if (nargin() > 0)
   interactive = 1;
else
   interactive = 0;
end

if (interactive)
 fprintf('Starting %s\n','plot_bins');
 plot_bins
 fprintf('Finished %s ... pausing, hit any key\n','plot_bins'); pause

 fprintf('Starting %s\n','plot_ens_err_spread');
 plot_ens_err_spread
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_err_spread'); pause

 fprintf('Starting %s\n','plot_ens_time_series');
 plot_ens_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_time_series'); pause

 fprintf('Starting %s\n','plot_ens_mean_time_series');
 plot_ens_mean_time_series
 fprintf('Finished %s ... pausing, hit any key\n','plot_ens_mean_time_series'); pause
end

 clear pinfo; close all; 
 
 truth_file = 'True_State.nc';
 diagn_file = 'Prior_Diag.nc';
 vars1 = CheckModel(diagn_file);
 vars1 = rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.var        = 'u';
 pinfo.levelindex = 1;
 pinfo.lonindex   = 69;
 pinfo.latindex   = 14;
 pinfo.level      = 1;
 pinfo.longitude  = 255.0;
 pinfo.latitude   = 38.966;
 clear vars1 vars2
 
 fprintf('Starting %s\n','PlotBins');
 close all; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 fprintf('Starting %s\n','PlotEnsErrSpread');
 close all; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsErrSpread'); pause

 fprintf('Starting %s\n','PlotEnsTimeSeries');
 pinfo.fname = pinfo.diagn_file;
 close all; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsTimeSeries'); pause

 fprintf('Starting %s\n','PlotEnsMeanTimeSeries');
 close all; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%------------------------------------------------------------
%plot_correl
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_correl');
 clear; close all; plot_correl
 fprintf('Finished %s pausing, hit any key\n','plot_correl'); pause
end

 clear pinfo;
 pinfo                    = CheckModel('Prior_Diag.nc');
 pinfo.time               = nc_varget(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
 pinfo.base_var           = 'v';
 pinfo.comp_var           = 'u';
 pinfo.base_tmeind        =   2;
 pinfo.base_lvlind        =   1;
 pinfo.base_lonind        =  69;
 pinfo.base_latind        =  14;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           =     1;
 pinfo.base_lon           = 255.0;
 pinfo.base_lat           =  38.96;
 pinfo.comp_lvlind        =    2;
 pinfo.comp_lvl           =    2;
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 
 fprintf('Starting %s\n','PlotCorrel');
 PlotCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotCorrel'); pause

%------------------------------------------------------------
%plot_phase_space
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_phase_space');
 clear; close all; plot_phase_space
 fprintf('Finished %s pausing, hit any key\n','plot_phase_space'); pause
end

 clear pinfo; clf
 pinfo              = CheckModel('Prior_Diag.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.var1name     = 'u';
 pinfo.var2name     = 'v';
 pinfo.var3name     = 'z';
 pinfo.var1_lvlind  = 1;
 pinfo.var2_lvlind  = 1;
 pinfo.var3_lvlind  = 1;
 pinfo.var1_lvl     = 1;
 pinfo.var2_lvl     = 1;
 pinfo.var3_lvl     = 1;
 pinfo.var1_latind  = 14;
 pinfo.var2_latind  = 14;
 pinfo.var3_latind  = 14;
 pinfo.var1_lat     = 38.96;
 pinfo.var2_lat     = 38.96;
 pinfo.var3_lat     = 38.96;
 pinfo.var1_lonind  = 69;
 pinfo.var2_lonind  = 69;
 pinfo.var3_lonind  = 69;
 pinfo.var1_lon     = 255;
 pinfo.var2_lon     = 255;
 pinfo.var3_lon     = 255;
 pinfo.ens_mem      = 'ensemble mean';
 pinfo.ltype        = 'k-';

 fprintf('Starting %s\n','PlotPhaseSpace');
 PlotPhaseSpace(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotPhaseSpace'); pause

%------------------------------------------------------------
%plot_reg_factor
%------------------------------------------------------------
% plot_reg_factor

%------------------------------------------------------------
%plot_sawtooth
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_sawtooth');
 clear; close all; plot_sawtooth
 fprintf('Finished %s pausing, hit any key\n','plot_sawtooth'); pause
end

 clear pinfo; clf

 truth_file     = 'True_State.nc';
 prior_file     = 'Prior_Diag.nc';
 posterior_file = 'Posterior_Diag.nc';
 pinfo = CheckModelCompatibility(prior_file,posterior_file);
 pinfo.prior_time     = pinfo.truth_time;
 pinfo.posterior_time = pinfo.diagn_time;
 pinfo.truth_file     = truth_file;
 pinfo.prior_file     = prior_file;
 pinfo.posterior_file = posterior_file;
 pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(prior_file);
 pinfo.var_names  = 'z';
 pinfo.levelindex  = 1;
 pinfo.latindex    = 14;
 pinfo.lonindex    = 69;
 pinfo.level       = 1;
 pinfo.latitude    = 38.96;
 pinfo.longitude   = 255;
 pinfo.copies      = 3;
 pinfo.copyindices = [7 12 17];

 fprintf('Starting %s\n','PlotSawtooth');
 PlotSawtooth(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotSawtooth'); pause

%------------------------------------------------------------
%plot_smoother_err
%------------------------------------------------------------
% plot_smoother_err

%------------------------------------------------------------
%plot_total_err
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_total_err');
 clear; close all; plot_total_err
 fprintf('Finished %s pausing, hit any key\n','plot_total_err'); pause
end

 clear pinfo; clf

 truth_file = 'True_State.nc';
 diagn_file = 'Prior_Diag.nc';
 vars1 = CheckModel(diagn_file);
 vars1 = rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.levels = nc_varget(pinfo.diagn_file,'lev');
 pinfo.lons   = nc_varget(pinfo.diagn_file,'lon');
 pinfo.lats   = nc_varget(pinfo.diagn_file,'lat');

 fprintf('Starting %s\n','PlotTotalErr');
 PlotTotalErr(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotTotalErr'); pause

%------------------------------------------------------------
%plot_var_var_correl
%------------------------------------------------------------

figure(1)
if (interactive)
 fprintf('Starting %s\n','plot_var_var_correl');
 clear; close all; plot_var_var_correl
 fprintf('Finished %s pausing, hit any key\n','plot_var_var_correl'); pause
end

 clear pinfo; clf
 diagn_file = 'Prior_Diag.nc';
 pinfo = CheckModel(diagn_file);
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.base_var           = 'u';
 pinfo.comp_var           = 'v';
 pinfo.base_tmeind        = 2;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           = 1;
 pinfo.base_lvlind        = 1;
 pinfo.base_lat           = 38.966;
 pinfo.base_latind        = 14;
 pinfo.base_lon           = 255.0;
 pinfo.base_lonind        = 69;
 pinfo.comp_lvl           = 1;
 pinfo.comp_lvlind        = 1;
 pinfo.comp_lat           = -38.966;
 pinfo.comp_latind        = 35;
 pinfo.comp_lon           = pinfo.base_lon;
 pinfo.comp_lonind        = pinfo.base_lonind;

 fprintf('Starting %s\n','PlotVarVarCorrel');
 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - correlation evolution
%------------------------------------------------------------
if (interactive)
 fprintf('Starting %s\n','plot_jeff_correl');
 clear; close all; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 clf

 fprintf('Starting %s\n','PlotJeffCorrel');
 PlotJeffCorrel(pinfo)
 fprintf('Finished %s\n','PlotJeffCorrel')

