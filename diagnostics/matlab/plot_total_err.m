function plot_total_err(diagn_file)
%% DART:plot_total_err - summary plots of global error and spread
%
% The true state ('true_state.nc') is REQUIRED.
%
% Example 1. Uses default name 'preassim.nc' and will prompt you 
%            for a new filename if 'preassim.nc' does not exist.)
% plot_total_err
%
% Example 2.
% diagn_file = 'preassim.nc';
% plot_total_err(diagn_file)

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

if (nargin == 0)
    disp('Input name of ensemble trajectory file:')
    diagn_file = input('<cr> for preassim.nc\n','s');
    if isempty(diagn_file)
        diagn_file = 'preassim.nc';
    end
elseif (nargin == 1)
    % all good - nothing to do
else
    error('Must supply exactly one filename or none.')
end

truth_file = 'true_state.nc';

if ( exist(truth_file,'file') ~= 2 ), error('%s does not exist.',truth_file); end
if ( exist(diagn_file,'file') ~= 2 ), error('%s does not exist.',diagn_file); end

fprintf('Comparing %s and \n          %s\n', truth_file, diagn_file)

CheckModel(diagn_file); % make sure model is supported - no need for anything else.
pinfo   = CheckModelCompatibility(truth_file,diagn_file);

switch lower(pinfo.model)
    case{'9var','lorenz_63','lorenz_84','lorenz_96','lorenz_96_2scale', ...
         'lorenz_04', 'forced_lorenz_96','ikeda','simple_advection'}

    case{'fms_bgrid'}
        pinfo = GetBgridInfo(pinfo, diagn_file, 'PlotTotalErr');

    case {'wrf'}
      pinfo = GetWRFInfo(pinfo, diagn_file, 'PlotTotalErr');

    case {'cam'}
      pinfo = GetCamInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'pe2lyr'}
      pinfo = GetPe2lyrInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'mitgcm_ocean'}
      pinfo = GetMITgcm_oceanInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'mpas_atm'}
      pinfo = GetMPAS_ATMInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'sqg'}
      pinfo = GetSqgInfo(pinfo, diagn_file, 'PlotTotalErr');

   case {'pop'}
      pinfo = GetPOPInfo(pinfo, diagn_file, 'PlotTotalErr');

   otherwise

      error('%s not implemented yet', pinfo.model)
end

PlotTotalErr( pinfo );
clear pinfo


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
