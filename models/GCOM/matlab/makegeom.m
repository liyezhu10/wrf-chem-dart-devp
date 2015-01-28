% Just making the coordinate variables in the shapes I believe we need them.

%% DART software - Copyright 2004 - 2013 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

%dimensions:
         nx = 97 ;
         ny = 33 ;
         nz = 33 ;
%        time = UNLIMITED ; // (1 currently)
         sx = 96 ;
         sy = 32 ;
         sz = 32 ;
%variables:
%        double ulon(sz, sy, nx) ;
%                ulon:units = "kilometers_east" ;
%                ulon:long_name = "x-coordinate of the grid nodes" ;
%        double ulat(sz, sy, nx) ;
%                ulat:units = "kilometers_east" ;
%                ulat:long_name = "x-coordinate of the grid nodes" ;
%        double ulev(sz, sy, nx) ;
%                ulev:units = "kilometers_east" ;
%                ulev:long_name = "x-coordinate of the grid nodes" ;

ulon = zeros(nx,sy,sz);
ulat = zeros(nx,sy,sz);
ulev = zeros(nx,sy,sz);

 nccreate('grid.nc','ulon',...
          'Dimensions',{'nx' nx 'sy' sy 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','ulat',...
          'Dimensions',{'nx' nx 'sy' sy 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','ulev',...
          'Dimensions',{'nx' nx 'sy' sy 'sz' sz},...
          'Format','classic');

ncwrite('grid.nc','ulon', ulon);
ncwrite('grid.nc','ulat', ulat);
ncwrite('grid.nc','ulev', ulev);

vlon = zeros(sx,ny,sz);
vlat = zeros(sx,ny,sz);
vlev = zeros(sx,ny,sz);

 nccreate('grid.nc','vlon',...
          'Dimensions',{'sx' sx 'ny' ny 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','vlat',...
          'Dimensions',{'sx' sx 'ny' ny 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','vlev',...
          'Dimensions',{'sx' sx 'ny' ny 'sz' sz},...
          'Format','classic');

ncwrite('grid.nc','vlon', vlon);
ncwrite('grid.nc','vlat', vlat);
ncwrite('grid.nc','vlev', vlev);

wlon = zeros(sx,sy,nz);
wlat = zeros(sx,sy,nz);
wlev = zeros(sx,sy,nz);

 nccreate('grid.nc','wlon',...
          'Dimensions',{'sx' sx 'sy' sy 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','wlat',...
          'Dimensions',{'sx' sx 'sy' sy 'nz' nz},...
          'Format','classic');

 nccreate('grid.nc','wlev',...
          'Dimensions',{'sx' sx 'sy' sy 'nz' nz},...
          'Format','classic');

ncwrite('grid.nc','wlon', wlon);
ncwrite('grid.nc','wlat', wlat);
ncwrite('grid.nc','wlev', wlev);

lon = zeros(sx,sy,sz);
lat = zeros(sx,sy,sz);
lev = zeros(sx,sy,sz);

 nccreate('grid.nc','lon',...
          'Dimensions',{'sx' sx 'sy' sy 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','lat',...
          'Dimensions',{'sx' sx 'sy' sy 'sz' sz},...
          'Format','classic');

 nccreate('grid.nc','lev',...
          'Dimensions',{'sx' sx 'sy' sy 'sz' sz},...
          'Format','classic');

 ncwrite('grid.nc','lon', lon);
 ncwrite('grid.nc','lat', lat);
 ncwrite('grid.nc','lev', lev);

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

