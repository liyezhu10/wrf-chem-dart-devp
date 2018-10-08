function plotit_reg(lon, lat, dat)
 
nlon = length(lon);
nlat = length(lat);

%% look for missing r8s
missing_val = -888888;
inds = find(dat==missing_val);
if (~isempty(inds)), dat(inds)  = nan; end

my_data = reshape(dat, nlon, nlat);

whos;

% plots the Y axis reversed from what we want.
% how to fix?
%scatter(lon, lat, 1, my_data);

imagesc(lon, lat, my_data);
colorbar;

end
