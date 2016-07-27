
nslat = 191;
nlon  = 288;

load RADIOSONDE_TEMPERATURE_T;

ind = 1;
sumx(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 2;
sumy(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 3;
sumx2(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 4;
sumy2(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 5;
numer(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 6;
denom(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 7;
alpha1(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 8;
alpha2(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);
ind = 9;
corr(1:nslat,:) = RADIOSONDE_TEMPERATURE_T((ind-1)*nslat+1:(ind-1)*nslat+nslat,:);

jeff=numer./denom;

