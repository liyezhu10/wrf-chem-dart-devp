% Generate approximately evenly distributed points on sphere using
% Golden Section spiral algorithm 
%    http://www.softimageblog.com/archives/115

% Lili 
% Atmosphere pressure to height (link: http://meteorologytraining.tpub.com/14269/css/14269_75.htm)
%         pressure (hPa)  height (m)
pres_hgt = [ 1000.0        111.0;...
              925.0        766.0;...
              850.0       1457.0;...            
              700.0       3012.0;...
              500.0       5574.0;...
              400.0       7185.0;...
              300.0       9164.0;...
              250.0      10363.0;...
              200.0      11784.0;...
              150.0      13608.0;...
              100.0      16180.0;...
               70.0      18442.0;...
               50.0      20576.0;...
               30.0      23849.0;...
               20.0      26481.0;...
               10.0      31055.0;];


%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL: https://subversion.ucar.edu/DAReS/DART/branches/development/observations/even_sphere/even_sphere.m $
% $Id: even_sphere.m 5839 2012-08-08 15:40:18Z nancy $
% $Revision: 5839 $
% $Date: 2012-08-08 09:40:18 -0600 (Wed, 08 Aug 2012) $

close all; clear;

% Generate obs_sequence input for this problem
% Mandatory levels are : 1000 mb, 850 mb, 700 mb, 500 mb, 400 mb, 
% 300 mb, 200 mb, 150 mb, 100 mb, 50 mb, 30 mb, 20 mb, 10 mb, 7 mb, 5 

% For diagnostic test need a null data value and a null qc value
diagnostic_test = false;

% Input is in hectopascals
%levels = [1000 850 700 500 400 300 200 150 100 50 30 20 10 7 5];
levels = [800 750 700 650 600 550 500 450 400 350 300 250 200 150];
% equvilent height of levels based on pres_hgt
%        [2   2.5 3   3.6 4.2 5   5.7 6.2 7.2 8    9  10  11.8 13.2] (unit:km)
% Hui comments: create synthetic obs from 2km to 13km around 500m apart,
%               could be sparser when height increases

% Date information is overwritten by create_fixed_network_sequence
year = 2008;
month = 10;
day = 1;
hour = 0;

% Number of roughly evenly distributed points in horizontal
n = 600;
x(1:n) = 0;
y(1:n) = 0;
z(1:n) = 0;

% Total number of observations at single time is levels*n*3
num_obs = size(levels, 2) * n * 1;
 
inc = pi * (3 - sqrt(5));
off = 2 / n;
for k = 1:n
   y(k) = (k-1) * off - 1 + (off / 2);
   r = sqrt(1 - y(k)*y(k));
   phi = (k-1) * inc; 
   x(k) = cos(phi) * r;
   z(k) = sin(phi) * r;
end 

% Now convert to latitude and longitude
for k = 1:n
   lon(k) = atan2(y(k), x(k)) + pi;
   lat(k) = asin(z(k));
   % Input is in degrees of lat and longitude
   dlon(k) = rad2deg(lon(k));
   dlat(k) = rad2deg(lat(k));
end

% Need to generate output for driving create_obs_sequence
fid = fopen('even_create_input_gps_pressure', 'w');

% Output the total number of observations
fprintf(fid, '%6i\n', num_obs);
% 0 copies of data
if(diagnostic_test) 
   fprintf(fid, '%2i\n', 1);
else
   fprintf(fid, '%2i\n', 0);
end

% 0 QC fields
if(diagnostic_test) 
   fprintf(fid, '%2i\n', 1);
else
   fprintf(fid, '%2i\n', 0);
end

% Metadata for data copy
if(diagnostic_test)
   fprintf(fid, '"observations"\n');
end
% Metadata for qc
if(diagnostic_test)
   fprintf(fid, '"Quality Control"\n');
end


% Loop to create each observation
for hloc = 1:n
   for vloc = 1:size(levels, 2)
      for field = 1:1
         % 0 indicates that there is another observation; 
         fprintf(fid, '%2i\n', 0);
         % Specify obs kind by string
         if(field == 1)
            fprintf(fid, 'GPSRO_REFRACTIVITY\n');
         end
         % Select local operator
         fprintf(fid, '%2i\n', 1);
         % Select pressure as the vertical coordinate
         fprintf(fid, '%2i\n', 2);
         % The vertical pressure level
         fprintf(fid, '%5i\n', levels(vloc));
         % Lon and lat in degrees next
         fprintf(fid, '%6.2f\n', dlon(hloc));
         fprintf(fid, '%6.2f\n', dlat(hloc));
         % Now the date and time
         fprintf(fid, '%5i %3i %3i %3i %2i %2i \n', year, month, day, hour, 0, 0);
         % Finally, the error variance, 1 for temperature, 4 for winds
         if(field == 1)
            fprintf(fid, '%2i\n', 0);
         else
            fprintf(fid, '%2i\n', 0);
         end
     
         % Need value and qc for testing
         if(diagnostic_test)
            fprintf(fid, '%2i\n', 1);
            fprintf(fid, '%2i\n', 0);
         end
      end   
   end
end

% File name
fprintf(fid, 'set_def.out\n');

% Done with output, close file
fclose(fid);



% Some plotting to visualize this observation set
plot(lon, lat, '*');
hold on

% Plot an approximate CAM T85 grid for comparison
% 256x128 points on A-grid
% Plot latitude lines first
del_lat = pi / 128;
for i = 1:128
   glat = del_lat * i - pi/2;
   a = [0 2*pi];
   b = [glat glat];
   plot(a, b, 'k');
end

del_lon = 2*pi / 256;
for i = 1:256
   glon = del_lon*i;
   a = [glon glon];
   b = [-pi/2,  pi/2];
   plot(a, b, 'k');
end

