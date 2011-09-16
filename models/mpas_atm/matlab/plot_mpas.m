function plot_mpas(fname)
%% plot_mpas 
% 
% fname = '../data/mpas_output.2010-10-23_03:00:00.nc';
% plot_mpas(fname)
%

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

x   = nc_varget(fname,'xCell');
y   = nc_varget(fname,'yCell');
z   = nc_varget(fname,'zCell');
c   = nc_varget(fname,'surface_pressure');
tri = nc_varget(fname,'cellsOnVertex');
h   = trisurf(tri,x,y,z,c);

set(h,'LineStyle','none')
axis equal
axis off
colorbar

shading interp

