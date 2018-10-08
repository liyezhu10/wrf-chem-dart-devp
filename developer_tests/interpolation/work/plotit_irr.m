function plotit_irr(lon, lat, dat, msize)
 
nlon = length(lon(1,:));
nlat = length(lat(:,1));

my_data = reshape(dat, nlon, nlat);

% look for missing r8s
missing_val = -888888;
inds = find(my_data==missing_val);
if (~isempty(inds)), my_data(inds)  = 0.; end
%if (~isempty(inds)), my_data(inds)  = nan; end

scatter(lon, lat, msize, my_data,'filled');
%plot3(lon, lat, my_data);
%imagesc(lon, lat, my_data);
colorbar;

end
