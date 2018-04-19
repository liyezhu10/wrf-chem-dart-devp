function h4 = plot_geo_grid(fname)
%% plot_grid ... plots the ULAT,ULON and TLAT,TLON variables from a netcdf file.
% 
% fname = '../data/DATA.ionos2.nc';
% plot_grid(fname)
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%
%	double time(time) ;
%		time:units = "seconds since 1966-01-01" ;
%	float geo_lon(geo_lon, geo_lat, geo_height) ;
%		geo_lon:short_name = "oplus geographic longitude" ;
%		geo_lon:units = "degrees" ;
%	float geo_lat(geo_lon, geo_lat, geo_height) ;
%		geo_lat:short_name = "oplus geographic latitude" ;
%		geo_lat:units = "degrees" ;
%	float geo_height(geo_lon, geo_lat, geo_height) ;
%		geo_height:short_name = "oplus height" ;
%		geo_height:units = "kilometers" ;
%	double oplus(geo_lon, geo_lat, geo_height) ;
%		oplus:short_name = "O+ number density" ;
%		oplus:units = "particles/cm**3" ;
%

lonmat = ncread(fname,'geo_lon');
latmat = ncread(fname,'geo_lat');
hgtmat = ncread(fname,'geo_height');

figure(1); clf; orient landscape
lon = lonmat(1,:,:);
lat = latmat(1,:,:);

h1 = plot(lon(:),lat(:),'ko');
title({'geographic is black o'})
xlabel('longitude')
ylabel('latitude')
set(gca,'FontSize',20)

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

