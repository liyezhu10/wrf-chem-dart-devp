function h = plot_feom_state(fname, varargin) 

% fname = '/glade/scratch/thoar/FILTER/Prior_Diag.nc';
% Simple Example:
% h = plot_feom_state(fname);
%
% Example:
% h = plot_feom_state(fname,'varname','temp','timestep', 1, 'copystring','ensemble spread') 

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

default_varname    = 'salt';
default_range      = [-Inf Inf];
default_timestep   =  1;
default_copystring = 'ensemble mean';
p = inputParser;

addRequired(  p, 'fname'     , @ischar);
addParamValue(p, 'varname'   , default_varname   , @ischar);
addParamValue(p, 'range'     , default_range     , @isnumeric);
addParamValue(p, 'time'      , default_timestep  , @isnumeric);
addParamValue(p, 'copystring', default_copystring, @ischar);
parse(p, fname, varargin{:});

% if you want to echo the input
fprintf('fname      : %s\n',     p.Results.fname)
fprintf('varname    : %s\n',     p.Results.varname)
fprintf('range      : %f %f \n', p.Results.range)
fprintf('timestep   : %d\n',     p.Results.time)
fprintf('copystring : %s\n',     p.Results.copystring)

if (exist(fname,'file') ~= 2)
   error('%s does not exist.',fname)
end

lons   = ncread(fname,'longitudes');
lats   = ncread(fname,'latitudes');
depths = ncread(fname,'depths');
s = zeros(size(lons)) + 3;   % markersize

copyindex = get_copy_index(fname, p.Results.copystring);

x = get_hyperslab('fname',fname, ...
                  'varname',p.Results.varname, ...
                  'copyindex',copyindex, ...
                  'tindex',p.Results.time);

h = scatter3(lons,lats,depths,s,x,'filled');

colorbar;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

