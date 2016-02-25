function plot_grid(fname)
%% plot_grid ... plots the ULAT,ULON and TLAT,TLON variables from a netcdf file.
% 
% fname = 'DATA.ionos2.nc';
% plot_grid(fname)
%

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

cg_lon     = nc_varget(fname,'cg_lon');
cg_lat     = nc_varget(fname,'cg_lat');
[cgx, cgy] = meshgrid(cg_lon, cg_lat);

ig_lon     = nc_varget(fname,'ig_lon');
ig_lat     = nc_varget(fname,'ig_lat');
[igx, igy] = meshgrid(ig_lon, ig_lat);

geo_lon = nc_varget(fname,'geo_lon');
geo_lat = nc_varget(fname,'geo_lat');

figure(1); clf; orient landscape

h1 = plot(cgx,cgy,'ko'); hold on;
h2 = plot(igx,igy,'rx');
h3 = plot(geo_lon,geo_lat,'gd');
% legend([h1,h2,h3],'geographic','magnetic','mag2geo')
title({'geographic is black o','magnetic is red x', 'mag2geo is green d'})
xlabel('longitude')
ylabel('latitude')

hold off;

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

