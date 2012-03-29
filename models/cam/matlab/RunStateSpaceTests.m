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
 plot_bins
 fprintf('Finished %s pausing, hit any key\n','plot_bins'); pause
 plot_ens_err_spread
 fprintf('Finished %s pausing, hit any key\n','plot_ens_err_spread'); pause
 plot_ens_time_series
 fprintf('Finished %s pausing, hit any key\n','plot_ens_time_series'); pause
 plot_ens_mean_time_series
 fprintf('Finished %s pausing, hit any key\n','plot_ens_mean_time_series'); pause
end

 clear pinfo; close all; 
 
 truth_file = './True_State.nc';
 diagn_file = './Prior_Diag.nc';
 vars1 = CheckModel(diagn_file);
 rmfield(vars1,{'time','time_series_length','fname'});
 vars2 = CheckModelCompatibility(truth_file,diagn_file);
 pinfo = CombineStructs(vars1,vars2);
 pinfo.var = 'T';
 pinfo.levelindex = 1;
 pinfo.lonindex   = 182;
 pinfo.latindex   = 93;
 pinfo.level      = 3.5;
 pinfo.longitude  = 254.5;
 pinfo.latitude   = 39.92;
 clear vars1 vars2
 
 close all; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 close all; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsErrSpread'); pause

 close all; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsTimeSeries'); pause

 close all; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%------------------------------------------------------------
%plot_correl
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_correl
 fprintf('Finished %s pausing, hit any key\n','plot_correl'); pause
end

 clear pinfo;
 pinfo                    = CheckModel('./Prior_Diag.nc');
 pinfo.time               = nc_varget(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
 pinfo.base_var           = 'T';
 pinfo.comp_var           = 'T';
 pinfo.base_tmeind        =   3;
 pinfo.base_lvlind        =   5;
 pinfo.base_lonind        = 182;
 pinfo.base_latind        =  93;
 pinfo.base_time          = pinfo.time(pinfo.base_tmeind);
 pinfo.base_lvl           =  37.2;
 pinfo.base_lon           = 254.5;
 pinfo.base_lat           =  39.9;
 pinfo.comp_lvlind        =  pinfo.base_lvlind;
 pinfo.comp_lvl           =  pinfo.base_lvl;
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 
 PlotCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotCorrel'); pause

%------------------------------------------------------------
%plot_phase_space
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_phase_space
 fprintf('Finished %s pausing, hit any key\n','plot_phase_space'); pause
end

 clear pinfo; clf
 pinfo              = CheckModel('./Prior_Diag.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);
 pinfo.var1name     = 'T';
 pinfo.var2name     = 'U';
 pinfo.var3name     = 'V';
 pinfo.var1_lvlind  = 1;
 pinfo.var2_lvlind  = 1;
 pinfo.var3_lvlind  = 1;
 pinfo.var1_lvl     = 3.5;
 pinfo.var2_lvl     = 3.5;
 pinfo.var3_lvl     = 3.5;
 pinfo.var1_latind  = 93;
 pinfo.var2_latind  = 93;
 pinfo.var3_latind  = 93;
 pinfo.var1_lat     = 39.92;
 pinfo.var2_lat     = 39.92;
 pinfo.var3_lat     = 39.92;
 pinfo.var1_lonind  = 182;
 pinfo.var2_lonind  = 182;
 pinfo.var3_lonind  = 182;
 pinfo.var1_lon     = 254.53;
 pinfo.var2_lon     = 254.53;
 pinfo.var3_lon     = 254.53;
 pinfo.ens_mem      = 'ensemble mean';
 pinfo.ltype        = 'k-';

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
 pinfo.var_names  = 'PS';
 pinfo.levelindex  = 1;
 pinfo.latindex    = 93;
 pinfo.lonindex    = 182;
 pinfo.level       = 1;
 pinfo.latitude    = 39.92;
 pinfo.longitude   = 254.53;
 pinfo.copies      = 3;
 pinfo.copyindices = [1 10 50];

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
 clear; close all; plot_total_err
 fprintf('Finished %s pausing, hit any key\n','plot_total_err'); pause
end

 clear pinfo; clf
 pinfo.model              = 'cam';
 pinfo.def_var            = 'state';
 pinfo.num_state_vars     = 9;
 pinfo.min_state_var      = 1;
 pinfo.max_state_var      = 9;
 pinfo.def_state_vars     = [1 2 3 4 5 6 7 8 9];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.truth_time         = [1 1000];
 pinfo.diagn_time         = [1 1000];
 pinfo.time_series_length = 1000;
 pinfo.time               = nc_varget(pinfo.truth_file,'time');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.diagn_file);

 PlotTotalErr(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotTotalErr'); pause

%------------------------------------------------------------
%plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_var_var_correl
 fprintf('Finished %s pausing, hit any key\n','plot_var_var_correl'); pause
end

 clear pinfo; clf
 pinfo.fname              = 'Prior_Diag.nc';
 pinfo.model              = 'cam';
 pinfo.base_var           = 'state';
 pinfo.state_var          = 'state';
 pinfo.base_var_index     = 4;
 pinfo.base_time          = 235;
 pinfo.state_var_index    = 8;
 pinfo.time_series_length = 1000;
 pinfo.time               = nc_varget(pinfo.fname,'time');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - virtually identical to plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; close all; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 clear pinfo; clf
 pinfo.fname              = 'Prior_Diag.nc';
 pinfo.model              = 'cam';
 pinfo.base_var           = 'state';
 pinfo.state_var          = 'state';
 pinfo.base_var_index     = 3;
 pinfo.base_time          = 300;
 pinfo.state_var_index    = 2;
 pinfo.time_series_length = 1000;
 pinfo.time               = nc_varget(pinfo.fname,'time');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

 PlotJeffCorrel(pinfo)

