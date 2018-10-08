function plotit_irr(lon, lat, dat, msize)
 
nlon = length(lon(1,:));
nlat = length(lat(:,1));


% look for missing r8s
missing_val = -888888;
inds = find(my_data==missing_val);
% choose nan or 0 for plotting, whichever gives
% the easiest-to-interpret results.
if (~isempty(inds)), my_data(inds)  = nan; end
%f (~isempty(inds)), my_data(inds)  = 0.; end

my_data = reshape(dat, nlon, nlat);

scatter(lon, lat, msize, my_data,'filled');
%plot3(lon, lat, my_data);
%imagesc(lon, lat, my_data);
colorbar;

end
