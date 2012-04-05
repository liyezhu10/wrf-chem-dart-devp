function RunAllTests(dummy)
%% RunAllTests.m

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

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
 pinfo = CheckModelCompatibility('True_State.nc','Prior_Diag.nc');
[pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.diagn_file);
 pinfo.var            = 'Y';
 pinfo.var_inds       = [100 200 300];
 
 clf; PlotBins(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotBins'); pause

 clf; PlotEnsErrSpread(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsErrSpread'); pause

 clf; PlotEnsTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsTimeSeries'); pause

 clf; PlotEnsMeanTimeSeries(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotEnsMeanTimeSeries'); pause

%------------------------------------------------------------
%plot_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_correl
 fprintf('Finished %s pausing, hit any key\n','plot_correl'); pause
end

 clear pinfo; clf
 pinfo = CheckModel('Prior_Diag.nc');
 pinfo.base_var           = 'X';
 pinfo.base_var_index     = 10;
 pinfo.base_time          = 125;
 
 PlotCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotCorrel'); pause

%------------------------------------------------------------
%plot_phase_space
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_phase_space
 fprintf('Finished %s pausing, hit any key\n','plot_phase_space'); pause
end

 clear pinfo; clf
 pinfo.fname    = 'Prior_Diag.nc';
 pinfo.model    = 'Lorenz_96_2scale';
 pinfo.var1name = 'X';
 pinfo.var2name = 'Y';
 pinfo.var3name = 'X';
 pinfo.var1ind  = 30;
 pinfo.var2ind  = 300;
 pinfo.var3ind  = 10;
 pinfo.ens_mem  = 'ensemble member10';
 pinfo.ltype    = 'k-';

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
 pinfo    = CheckModelCompatibility('Prior_Diag.nc','Posterior_Diag.nc');
 pinfo.prior_time     = pinfo.truth_time;
 pinfo.prior_file     = pinfo.truth_file;
 pinfo.posterior_time = pinfo.diagn_time;
 pinfo.posterior_file = pinfo.diagn_file;
 pinfo.truth_file     = 'True_State.nc';
 pinfo = rmfield(pinfo,{'diagn_file','truth_time','diagn_time'});
 [pinfo.num_ens_members, pinfo.ensemble_indices] = get_ensemble_indices(pinfo.prior_file);
 pinfo.var            = 'X';
 pinfo.var_inds       = [1 12 24];
 pinfo.copyindices    = [7 12 17];

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

 clear pinfo; close all
 pinfo    = CheckModelCompatibility('True_State.nc','Prior_Diag.nc');
 pinfo.def_state_vars     = [1 12 24];

 PlotTotalErr(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotTotalErr'); pause

%------------------------------------------------------------
%plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_var_var_correl
 fprintf('Finished %s pausing, hit any key\n','plot_var_var_correl'); pause
end

 clear pinfo; clf
 pinfo  = CheckModel('Prior_Diag.nc');
 pinfo.base_var        = 'X';
 pinfo.state_var       = 'Y';
 pinfo.base_var_index  = 18;
 pinfo.base_time       = 128;
 pinfo.state_var_index = 180;

 PlotVarVarCorrel(pinfo)
 fprintf('Finished %s pausing, hit any key\n','PlotVarVarCorrel'); pause

%------------------------------------------------------------
%plot_jeff_correl - identical inputs as plot_var_var_correl
%------------------------------------------------------------
if (interactive)
 clear; clf; plot_jeff_correl
 fprintf('Finished %s pausing, hit any key\n','plot_jeff_correl'); pause
end

 PlotJeffCorrel(pinfo)

