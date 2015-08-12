%% Check_netCDF checks to make sure DART can read and reproduce JULES variables.

%% DART software - Copyright 2004 - 2015 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

% Check to see if the DART diagnostic netCDF files faithfully recreate the
% contents of the Jules restart and output files.

jules_restart = '../work/jules_restart.nc';
jules_output  = '../work/jules_output.nc';
dart_file     = 'check_me_out.nc';

% The DART files have the extra 'copy' and are guaranteed to have a 'time'

% No time dimension for these - easy.
restart_variables = {'canopy', 'cs', 'gs', 'snow_tile', 'sthuf', 't_soil', 'tstar_tile', 'rho_snow', 'snow_depth'};

% Must pick off last time 
output_variables = { 'time', 'latitude', 'longitude', 'precip', 'latent_heat', 'smcl', 'soil_wet', 'esoil'};

%% Work on the restart variables first
%  netcdf jules_restart {
%  dimensions:
%          land = 679 ;
%          tile = 9 ;
%          scpool = 1 ;
%          soil = 4 ;
%  variables:
%          float canopy(    tile, land) ;
%          float cs(      scpool, land) ;
%          float gs(land) ;
%          float snow_tile( tile, land) ;
%          float sthuf(     soil, land) ;
%          float t_soil(    soil, land) ;
%          float tstar_tile(tile, land) ;
%          float rho_snow(  tile, land) ;
%          float snow_depth(tile, land) ;

for ivar=1:length(restart_variables)

   % ncread automatically squeezes out singleton dimensions.

   jules = ncread(jules_restart,restart_variables{ivar});
   dart  = ncread(dart_file,restart_variables{ivar});

   clf;
   plot(jules)
   hold on
   plot(dart)
   legend('jules','dart')
   title(restart_variables{ivar})

   diff = jules - dart;
   fprintf('Maximum %s difference %f \n', restart_variables{ivar}, max(diff(:)))
   fprintf('Minimum %s difference %f \n', restart_variables{ivar}, min(diff(:)))

   disp('Pausing, hit any key to continue ...')
   pause;

end


% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

