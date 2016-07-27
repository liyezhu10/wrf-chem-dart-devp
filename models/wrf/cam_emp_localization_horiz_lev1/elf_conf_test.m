
tempvar = load('zzz_104');

tmpx=tempvar(:,1);  % x_truth - x_prior_mean
tmpy=tempvar(:,2);  % delta_y_mean

%% regression analysis
p=polyfit(tmpy,tmpx,1);
f=polyval(tmpy,p);

nsample = size(tmpx,1);
rcoef = p(1);
df = nsample-2;          % degree of freedom
var_tmpx=var(tmpx);
var_tmpy=var(tmpy);
se_b = sqrt(var_tmpx*(nsample-1)/(nsample-2))/sqrt(var_tmpy*(nsample-1));
tsig = rcoef/se_b;


%% pseudo-group alpha
nsample = size(tmpx,1);
%rand_sample = 1000;
%ngroup = floor(nsample/rand_sample);
ngroup = 8;
rand_sample = floor(nsample/ngroup);

rindex(1:rand_sample,1:ngroup) = 0;
tmpx_group(1:rand_sample,1:ngroup) = 0.0;
tmpy_group(1:rand_sample,1:ngroup) = 0.0;
alpha_group(1:ngroup) = 0.0;
ssxx(1:ngroup) = 0.0;
ssyy(1:ngroup) = 0.0;
ssxy(1:ngroup) = 0.0;
r2(1:ngroup) = 0.0;
se_b(1:ngroup) = 0.0;

%for i = 1:ngroup
%    rindex(1:rand_sample,i) = randi(nsample,rand_sample,1);
%end

rindex_tot=randsample(nsample,nsample);
for i = 1:ngroup
    rindex(1:rand_sample,i) = rindex_tot((i-1)*rand_sample+1:i*rand_sample);
end

for i = 1:ngroup
    for j = 1:rand_sample
        tmpx_group(j,i) = tmpx(rindex(j,i));
        tmpy_group(j,i) = tmpy(rindex(j,i));
    end
end
for i = 1:ngroup
    p_tmp = polyfit(tmpy_group(:,i),tmpx_group(:,i),1);
    alpha_group(i) = p_tmp(1);
    
    ssxx(i) = sum((tmpy_group(:,i)-mean(tmpy_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
    ssyy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpx_group(:,i)-mean(tmpx_group(:,i))));
    ssxy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
    r2(i) = ssxy(i)^2/(ssxx(i)*ssyy(i));
    se_b(i) = sqrt((ssyy(i)-ssxy(i)^2/ssxx(i))/(rand_sample-2))/sqrt(ssxx(i));
end
sum_alpha  = 0.0;
sum_alpha2 = 0.0;
for i = 1:ngroup
    sum_alpha  = sum_alpha + alpha_group(i);
    sum_alpha2 = sum_alpha2 + alpha_group(i)^2;
end
beta = (sum_alpha^2/sum_alpha2 - 1)/(ngroup-1);



%% random number from tmpx and tmpy

%meanx=mean(tmpx);
%spdx=sqrt(var(tmpx));
%
%meany=mean(tmpy);
%spdy=sqrt(var(tmpy));
%
%sample_size = 10000;
%
%randnx=randn(sample_size,1);
%randnx=spdx*randnx+meanx;
%
%randny=randn(sample_size,1);
%randny=spdy*randny+meany;
%
%p=polyfit(randny,randnx,1);
%f=polyval(p,randny);
%
%plot(randny,randnx,'bo');
%hold on;
%plot(randny,f,'r-');
%hold off






