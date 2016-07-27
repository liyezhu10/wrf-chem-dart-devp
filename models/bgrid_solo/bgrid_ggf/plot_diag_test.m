
len_dist   = 632;
dlen       = 0.005;
x_dist     = 0:dlen:(len_dist-1)*dlen;
bin        = 10;
len_dist_bin = floor((len_dist-1)/bin)+1;
x_dist_bin = 0:dlen*bin:(len_dist_bin-1)*dlen*bin;
write_dist = 450;


vartmp = load('PRESSURE_ps_obslev1_varlev1_cutoff0400');

meangf=vartmp(4,:)./vartmp(1,:);
ggf=(vartmp(2,:)./vartmp(3,:)-1.0);
plot(meangf);

addpath /glade/user/lililei/work/BGrid/emp_loc/GGF/matlab_dir
cutoff = 0.4;
for i = 1:len_dist;
gc(i)=comp_GC(x_dist(i),cutoff);
end
for i  = 1:len_dist_bin;
gc_bin(i)=comp_GC(x_dist_bin(i),cutoff);
end

figure(1);
plot(x_dist,meangf,'LineWidth',2.0);
hold on;
plot(x_dist,ggf,'r-','LineWidth',2.0);
plot(x_dist,gc,'k-')
grid on;


[nx ny]=size(vartmp);
vartmp_bin = zeros(nx,len_dist_bin);

for i = 1:len_dist_bin

    if ( i == 1 ) 
         inds = 1;
         inde = 1;
    else
         inds = 1+(i-2)*bin+1;
         inde = 1+(i-1)*bin;
    end

    tutu_tmp = squeeze(vartmp(:,inds:inde));
    vartmp_bin(:,i) = nansum(tutu_tmp,2);

end
meangf_bin = vartmp_bin(4,:)./vartmp_bin(1,:);
ggf_bin = (vartmp_bin(2,:)./vartmp_bin(3,:)-1.0)/(2.0-1.0);
figure(2);
plot(x_dist_bin,meangf_bin,'LineWidth',2.0);
hold on;
plot(x_dist_bin,ggf_bin,'r-','LineWidth',2.0);
plot(x_dist_bin,gc_bin,'k-');
grid on;






