
totfile = 1;
% inputname = {'U10_U'};
len_dist = 1600;
dlen     = 0.002;
% len_ml   = 49;

% for plotting use
earth_radius = 6370;              % unit km
cen_lat      = 40.8491;           % center latitude
cen_lat_rad  = cen_lat*pi/180.0;
R_cen_lat    = earth_radius*cos(cen_lat_rad);
% x_ml=[1:1:49];

% process data: bin to larger subsets
bin = 10;
len_dist_bin = len_dist/bin;

obstypename = {'RADIOSONDE_PRESSURE',...
               'RADIOSONDE_TEMPERATURE',...
               'RADIOSONDE_U_WIND_COMPONENT',...
               'RADIOSONDE_V_WIND_COMPONENT'
              };
nobstype = size(obstypename,2);

statename = {'PS','T','US','VS','Q','CLDICE','CLDLIQ'};
nstatetype = size(statename,2);

timename={'2008090112','2008090200','2008090212','2008090300','2010122918','2010122921','2010123000',...
          '2010123003','2010123006','2010123009','2010123012','2010123015','2010123018','2010123021','2010123100',...
          '2010123103','2010123106','2010123109','2010123112','2010123115','2010123118','2010123121','2011010100',...
          '2011010103','2011010106','2011010109','2011010112','2011010115','2011010118','2011010121','2011010200',...
          '2011010203','2011010206','2011010209','2011010112','2011010215','2011010218','2011010221','2011010300',...
          '2011010303','2011010306','2011010309','2011010312','2011010315','2011010318','2011010321'};



for idobs = 4:4
for inddim = 4:4
    
%% load in data
for ifile = 1:totfile
%    string = sprintf('load %s_%s_d0%d_%s',obstypename{idobs},statename{inddim},id,timename{ifile});
%    eval(string);
%    string = sprintf('vartmp = %s_%s_d0%d_%s;',obstypename{idobs},statename{inddim},id,timename{ifile});
%    eval(string);

%     string = sprintf('load %s_%s_%s',obstypename{idobs},statename{inddim},timename{ifile});
%     eval(string);
%     string = sprintf('vartmp = %s_%s_%s;',obstypename{idobs},statename{inddim},timename{ifile});
%     eval(string);

    string = sprintf('load %s_%s',obstypename{idobs},statename{inddim});
    eval(string);
    string = sprintf('vartmp = %s_%s;',obstypename{idobs},statename{inddim});
    eval(string);

    var_dist.num(1:len_dist,ifile)=vartmp(1,1:len_dist);
    var_dist.numer(1:len_dist,ifile)=vartmp(2,1:len_dist);
    var_dist.denom(1:len_dist,ifile)=vartmp(3,1:len_dist);
    var_dist.alpha(1:len_dist,ifile)=vartmp(4,1:len_dist);
    var_dist.alpha2(1:len_dist,ifile)=vartmp(5,1:len_dist);
    var_dist.corr(1:len_dist,ifile)=vartmp(6,1:len_dist);
    
end


%% sum the data in time
tot_dist.num    = sum(var_dist.num,2);
tot_dist.numer  = sum(var_dist.numer,2);
tot_dist.denom  = sum(var_dist.denom,2);
tot_dist.alpha  = sum(var_dist.alpha,2);
tot_dist.alpha2 = sum(var_dist.alpha2,2);
tot_dist.corr   = sum(var_dist.corr,2);


%% bin the data
var_dist_bin.num(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.numer(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.denom(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.alpha(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.alpha2(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.corr(1:len_dist_bin,1:totfile) = 0.0;

tot_dist_bin.num(1:len_dist_bin) = 0.0;
tot_dist_bin.numer(1:len_dist_bin) = 0.0;
tot_dist_bin.denom(1:len_dist_bin) = 0.0;
tot_dist_bin.alpha(1:len_dist_bin) = 0.0;
tot_dist_bin.alpha2(1:len_dist_bin) = 0.0;
tot_dist_bin.corr(1:len_dist_bin) = 0.0;

for k = 1:totfile
for i = 1:len_dist_bin
    startbin = (i-1)*bin + 1;
    endbin   = i*bin;
    for j = startbin:endbin
        var_dist_bin.num(i,k)    = var_dist_bin.num(i,k) + var_dist.num(j,k);
        var_dist_bin.numer(i,k)  = var_dist_bin.numer(i,k) + var_dist.numer(j,k);
        var_dist_bin.denom(i,k)  = var_dist_bin.denom(i,k) + var_dist.denom(j,k);
        var_dist_bin.alpha(i,k)  = var_dist_bin.alpha(i,k) + var_dist.alpha(j,k);
        var_dist_bin.alpha2(i,k) = var_dist_bin.alpha2(i,k) + var_dist.alpha2(j,k);
        var_dist_bin.corr(i,k)   = var_dist_bin.corr(i,k) + var_dist.corr(j,k);
    end
end
end

for i = 1:len_dist_bin
    startbin = (i-1)*bin + 1;
    endbin   = i*bin;
    for j = startbin:endbin
        tot_dist_bin.num(i)    = tot_dist_bin.num(i) + tot_dist.num(j);
        tot_dist_bin.numer(i)  = tot_dist_bin.numer(i) + tot_dist.numer(j);
        tot_dist_bin.denom(i)  = tot_dist_bin.denom(i) + tot_dist.denom(j);
        tot_dist_bin.alpha(i)  = tot_dist_bin.alpha(i) + tot_dist.alpha(j);
        tot_dist_bin.alpha2(i) = tot_dist_bin.alpha2(i) + tot_dist.alpha2(j);
        tot_dist_bin.corr(i)   = tot_dist_bin.corr(i) + tot_dist.corr(j);
    end
    
end


%% do the computations for original data

% .jeff is the alpha defined in Jeff's EOL paper
% .beta is the confidence value of alpha
% .avgcorr is the average of absolute correlation

% for physical distance computations
% compute alpha and beta for total samples (sum in all analysis times)
for i = 1:len_dist
    tot_dist.jeff(i)=tot_dist.numer(i)/tot_dist.denom(i);
    tot_dist.beta(i)=(tot_dist.alpha(i)^2/tot_dist.alpha2(i)-1.0)/(tot_dist.num(i)-1.0);
    tot_dist.avgcorr(i)=tot_dist.corr(i)/tot_dist.num(i);
end

% compute alpha and beta for each analysis time
for j = 1:totfile
    for i = 1:len_dist
        var_dist.jeff(i,j)=var_dist.numer(i,j)/var_dist.denom(i,j);
        var_dist.beta(i,j)=(var_dist.alpha(i,j)^2/var_dist.alpha2(i,j)-1.0)/(var_dist.num(i,j)-1.0);
        var_dist.avgcorr(i,j)=var_dist.corr(i,j)/var_dist.num(i,j);
    end
end

% % compute the mean(.jeffmean) and confidence value(.jeffbeta) of the alpha at each analysis time
% % also compute the mean of (avgcorr) at each analysis time
% for i = 1:len_dist
%     temp(1:totfile)=var_dist.jeff(i,1:totfile);
%     var_dist.jeffmean(i)=mean(temp);
%     var_dist.jeffvar(i)=var(temp);
%     
%     temp1=0.0;
%     temp2=0.0;
%     for j=1:totfile
%         temp1=temp1+var_dist.jeff(i,j);
%         temp2=temp2+var_dist.jeff(i,j)^2;
%     end
%     var_dist.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
%     
%     temp(1:totfile)=var_dist.avgcorr(i,1:totfile);
%     var_dist.avgcorrmean(i)=mean(temp);
%     
% end


%% do the computations for binned data

% .jeff is the alpha defined in Jeff's EOL paper
% .beta is the confidence value of alpha
% .avgcorr is the average of absolute correlation

% for physical distance computations
% compute alpha and beta for total samples (sum in all analysis times)
for i = 1:len_dist_bin
    tot_dist_bin.jeff(i)=tot_dist_bin.numer(i)/tot_dist_bin.denom(i);
    tot_dist_bin.beta(i)=(tot_dist_bin.alpha(i)^2/tot_dist_bin.alpha2(i)-1.0)/(tot_dist_bin.num(i)-1.0);
    tot_dist_bin.avgcorr(i)=tot_dist_bin.corr(i)/tot_dist_bin.num(i);
end

% compute alpha and beta for each analysis time
for j = 1:totfile
    for i = 1:len_dist_bin
        var_dist_bin.jeff(i,j)=var_dist_bin.numer(i,j)/var_dist_bin.denom(i,j);
        var_dist_bin.beta(i,j)=(var_dist_bin.alpha(i,j)^2/var_dist_bin.alpha2(i,j)-1.0)/(var_dist_bin.num(i,j)-1.0);
        var_dist_bin.avgcorr(i,j)=var_dist_bin.corr(i,j)/var_dist_bin.num(i,j);
    end
end

% % compute the mean(.jeffmean) and confidence value(.jeffbeta) of the alpha at each analysis time
% % also compute the mean of (avgcorr) at each analysis time
% for i = 1:len_dist_bin
%     temp(1:totfile)=var_dist_bin.jeff(i,1:totfile);
%     var_dist_bin.jeffmean(i)=mean(temp);
%     var_dist_bin.jeffvar(i)=var(temp);
%     
%     temp1=0.0;
%     temp2=0.0;
%     for j=1:totfile
%         temp1=temp1+var_dist_bin.jeff(i,j);
%         temp2=temp2+var_dist_bin.jeff(i,j)^2;
%     end
%     var_dist_bin.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
%     
%     temp(1:totfile)=var_dist_bin.avgcorr(i,1:totfile);
%     var_dist_bin.avgcorrmean(i)=mean(temp);
% end



%% plot section

%x_dist=[0:dlen*R_cen_lat:(len_dist-1)*dlen*R_cen_lat];
%x_dist_bin=[0:dlen*R_cen_lat*bin:(len_dist_bin-1)*dlen*R_cen_lat*bin];

x_dist = 0:dlen:(len_dist-1)*dlen;
x_dist_bin = 0:dlen*bin:(len_dist_bin-1)*dlen*bin;

% plotting
% x_dist=[0:dlen:len_dist-1];
% x_ml=[1:1:49];

%psfname = sprintf('%s_%s_d0%d.ps',obstypename{idobs},statename{inddim},id);
%disp(sprintf('Removing %s from current directory',psfname))
%system(sprintf('rm %s',psfname));

figure(1);
[ax,h1,h2]=plotyy(x_dist,tot_dist.jeff,x_dist,tot_dist.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
% plot(x_dist,tot_dist.jeff,'b-','LineWidth',2.0);
% hold on;
% plot(x_dist,tot_dist.beta,'r-','LineWidth',2.0);
% hold on;
% plot(x_dist,tot_dist.avgcorr,'k--','LineWidth',2.0);
hold off;
% legend('alpha','beta','abs corr');
hold off;
xlabel('distance (unit: radiance)','FontSize',14);
title(sprintf('%s %s',obstypename{idobs},statename{inddim}),'FontSize',14);
%print('-dpsc',psfname,'-append')
%clf;
% figure(3);
% plot((1:1:100),tot_dist.jeff(1:100));


figure(5);
 [ax,h1,h2]=plotyy(x_dist_bin,tot_dist_bin.jeff,x_dist_bin,tot_dist_bin.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
% plot(x_dist_bin,tot_dist_bin.jeff,'b-','LineWidth',2.0);
% hold on;
% plot(x_dist_bin,tot_dist_bin.beta,'r-','LineWidth',2.0);
% hold on;
% plot(x_dist_bin,tot_dist_bin.avgcorr,'k--','LineWidth',2.0);
hold off;
% legend('alpha','beta','abs corr');
hold off;
xlabel('distance (unit: radiance)','FontSize',14);
title(sprintf('%s %s',obstypename{idobs},statename{inddim}),'FontSize',14);
%print('-dpsc',psfname,'-append')
%clf;
% figure(3);
% plot((1:1:100),tot_dist.jeff(1:100));


for ifile = 1:totfile
    figure(10+ifile);
    [ax,h1,h2]=plotyy(x_dist_bin,var_dist_bin.jeff(:,ifile),x_dist_bin,var_dist_bin.num(:,ifile));
    xlabel('distance (unit: radiance)','FontSize',14);
    title(sprintf('%s %s',obstypename{idobs},statename{inddim}),'FontSize',14);
end



figure(101);
for ifile = 1:totfile
    plot(x_dist_bin,var_dist_bin.jeff(:,ifile));
    hold on;
end
    hold off;


%figure(2);
%plot(x_dist,var_dist.jeffmean,'b-');
%hold on;
%plot(x_dist,var_dist.jeffbeta,'r-.');
%hold on;
%plot(x_dist,var_dist.avgcorrmean,'k--');
%hold off;
%legend('alpha','beta','abs corr');
%hold off;
%xlabel('distance (unit: km)','FontSize',14);
%title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
%print('-dpsc',psfname,'-append')
%clf;

% figure(3);
% [ax,h1,h2]=plotyy(x_ml,tot_ml.jeff,x_ml,tot_ml.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
% % set(h2,'color','b');
% % plot(tot_ml.jeff,'b-');
% hold on;
% plot(tot_ml.beta,'r-','LineWidth',2.0);
% hold on;
% plot(tot_ml.avgcorr,'k--','LineWidth',2.0);
% legend('alpha','beta','abs corr');
% hold off;
% xlabel('model levels','FontSize',14);
% title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
% %print('-dpsc',psfname,'-append')
% %clf;

%figure(4);
%plot(var_ml.jeffmean,'b-');
%hold on;
%plot(var_ml.jeffbeta,'r-.');
%hold on;
%plot(var_ml.avgcorrmean,'k--');
%legend('alpha','beta','abs corr');
%hold off;
%xlabel('model levels','FontSize',14);
%title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
%%print('-dpsc',psfname,'-append')
%%clf;


end   % end inddim
end   % end idobs

