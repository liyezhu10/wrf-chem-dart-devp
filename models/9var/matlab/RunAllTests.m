function RunAllTests(dummy)
%% RunAllTests.m

%------------------------------------------------------------
% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
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
 pinfo          = CheckModelCompatibility('True_State.nc','Prior_Diag.nc');
 pinfo.var      = 'state';
 pinfo.var_inds = [1 2 3 4 5 6 7 8 9];
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.diagn_file);
 
 close all; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 pinfo.time_series_length = 1000;
 pinfo.time = nc_varget(pinfo.truth_file,'time');

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
 pinfo                    = CheckModel('Prior_Diag.nc');
 pinfo.base_var           = 'state';
 pinfo.base_var_index     = 4;
 pinfo.base_time          = 34;
 pinfo.time               = nc_varget(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
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
 pinfo.fname    = 'True_State.nc';
 pinfo.model    = '9var';
 pinfo.var1name = 'state';
 pinfo.var2name = 'state';
 pinfo.var3name = 'state';
 pinfo.var1ind  = 1;
 pinfo.var2ind  = 2;
 pinfo.var3ind  = 3;
 pinfo.ens_mem  = 'true state';
 pinfo.ltype    = 'k-';
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

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

 clear pinfo; close all
 pinfo.model              = '9var';
 pinfo.def_var            = 'state';
 pinfo.num_state_vars     = 9;
 pinfo.min_state_var      = 1;
 pinfo.max_state_var      = 9;
 pinfo.def_state_vars     = [1 2 3 4 5 6 7 8 9];
 pinfo.prior_file         = 'Prior_Diag.nc';
 pinfo.posterior_file     = 'Posterior_Diag.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.diagn_time         = [-1 -1];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.truth_time         = [-1 -1];
 pinfo.var                = 'state';
 pinfo.var_inds           = [1 2 3 4 5 6 7 8 9];
 pinfo.copyindices        = [7 12 17];
 pinfo.time               = nc_varget(pinfo.truth_file,'time');
 pinfo.time_series_length = length(pinfo.time);
 pinfo.prior_time         = [1 pinfo.time_series_length];
 pinfo.posterior_time     = [1 pinfo.time_series_length];
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.prior_file);

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
 pinfo.model              = '9var';
 pinfo.def_var            = 'state';
 pinfo.num_state_vars     = 9;
 pinfo.min_state_var      = 1;
 pinfo.max_state_var      = 9;
 pinfo.def_state_vars     = [1 2 3 4 5 6 7 8 9];
 pinfo.truth_file         = 'True_State.nc';
 pinfo.diagn_file         = 'Prior_Diag.nc';
 pinfo.time               = nc_varget(pinfo.truth_file,'time');
 pinfo.time_series_length = length(pinfo.time);
 pinfo.truth_time         = [1 pinfo.time_series_length];
 pinfo.diagn_time         = [1 pinfo.time_series_length];
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
 pinfo.model              = '9var';
 pinfo.base_var           = 'state';
 pinfo.state_var          = 'state';
 pinfo.base_var_index     = 4;
 pinfo.base_time          = 235;
 pinfo.state_var_index    = 8;
 pinfo.time               = nc_varget(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
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
 pinfo.model              = '9var';
 pinfo.base_var           = 'state';
 pinfo.state_var          = 'state';
 pinfo.base_var_index     = 3;
 pinfo.state_var_index    = 2;
 pinfo.time               = nc_varget(pinfo.fname,'time');
 pinfo.time_series_length = length(pinfo.time);
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.fname);

 PlotJeffCorrel(pinfo)

