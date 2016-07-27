
totfile = 46;
% inputname = {'U10_U'};
len_dist = 81;
dlen     = 1;
len_ml   = 49;

wrf_ndom = 2;

obstypename = {'SFC_ALTIMETER',...
               'SFC_U_WIND_COMPONENT',...
               'SFC_V_WIND_COMPONENT',...
               'SFC_TEMPERATURE',...
               'SFC_DEWPOINT',...
               'SFC_SPECIFIC_HUMIDITY'...
               'ACARS_U_WIND_COMPONENT',...
               'ACARS_V_WIND_COMPONENT',...
               'ACARS_TEMPERATURE',...
               'ACARS_DEWPOINT',...
               'SAT_U_WIND_COMPONENT',...
               'SAT_V_WIND_COMPONENT',...
               'RADIOSONDE_SURFACE_ALTIMETER',...
               'RADIOSONDE_U_WIND_COMPONENT',...
               'RADIOSONDE_V_WIND_COMPONENT',...
               'RADIOSONDE_TEMPERATURE',...
               'RADIOSONDE_DEWPOINT'
              };
nobstype = size(obstypename,2);

% statename = {'U10','V10','T2','Q2','TH2','MU','PSFC',...   % 2D variables
statename = {'U','V','W','T','PH'...
             'QVAPOR','QCLOUD','QRAIN','QNRAIN',...
             'QICE','QNICE','QSNOW','QGRAUP'};
nstatetype = size(statename,2);

timename={'2010122906','2010122909','2010122912','2010122915','2010122918','2010122921','2010123000',...
          '2010123003','2010123006','2010123009','2010123012','2010123015','2010123018','2010123021','2010123100',...
          '2010123103','2010123106','2010123109','2010123112','2010123115','2010123118','2010123121','2011010100',...
          '2011010103','2011010106','2011010109','2011010112','2011010115','2011010118','2011010121','2011010200',...
          '2011010203','2011010206','2011010209','2011010112','2011010215','2011010218','2011010221','2011010300',...
          '2011010303','2011010306','2011010309','2011010312','2011010315','2011010318','2011010321'};



for id = 2:wrf_ndom
for idobs = 1:1
for inddim = 1:5
    
% load in data
for ifile = 1:totfile
    string = sprintf('load %s_%s_d0%d',obstypename{idobs},statename{inddim},id);
    eval(string);
    string = sprintf('vartmp = %s_%s_d0%d;',obstypename{idobs},statename{inddim},id);
    eval(string);

%     string = sprintf('load %s_%s_d0%d',obstypename{idobs},statename{inddim},id);
%     eval(string);
%     string = sprintf('vartmp = %s_%s_d0%d;',obstypename{idobs},statename{inddim},id);
%     eval(string);

    var_dist.num(1:len_dist,ifile)=vartmp(1,1:len_dist);
    var_dist.numer(1:len_dist,ifile)=vartmp(2,1:len_dist);
    var_dist.denom(1:len_dist,ifile)=vartmp(3,1:len_dist);
    var_dist.alpha(1:len_dist,ifile)=vartmp(4,1:len_dist);
    var_dist.alpha2(1:len_dist,ifile)=vartmp(5,1:len_dist);
    var_dist.corr(1:len_dist,ifile)=vartmp(6,1:len_dist);
    var_dist.abscorr(1:len_dist,ifile)=vartmp(7,1:len_dist);
    
    var_ml.num(1:len_ml,ifile)=vartmp(8,1:len_ml);
    var_ml.numer(1:len_ml,ifile)=vartmp(9,1:len_ml);
    var_ml.denom(1:len_ml,ifile)=vartmp(10,1:len_ml);
    var_ml.alpha(1:len_ml,ifile)=vartmp(11,1:len_ml);
    var_ml.alpha2(1:len_ml,ifile)=vartmp(12,1:len_ml);
    var_ml.corr(1:len_ml,ifile)=vartmp(13,1:len_ml);
    var_ml.abscorr(1:len_ml,ifile)=vartmp(14,1:len_ml);
 
    var100_dist.num(1:len_dist,ifile)=vartmp(15,1:len_dist);
    var100_dist.numer(1:len_dist,ifile)=vartmp(16,1:len_dist);
    var100_dist.denom(1:len_dist,ifile)=vartmp(17,1:len_dist);
    var100_dist.alpha(1:len_dist,ifile)=vartmp(18,1:len_dist);
    var100_dist.alpha2(1:len_dist,ifile)=vartmp(19,1:len_dist);
    var100_dist.corr(1:len_dist,ifile)=vartmp(20,1:len_dist);
    var100_dist.abscorr(1:len_dist,ifile)=vartmp(21,1:len_dist);

    var100_ml.num(1:len_ml,ifile)=vartmp(22,1:len_ml);
    var100_ml.numer(1:len_ml,ifile)=vartmp(23,1:len_ml);
    var100_ml.denom(1:len_ml,ifile)=vartmp(24,1:len_ml);
    var100_ml.alpha(1:len_ml,ifile)=vartmp(25,1:len_ml);
    var100_ml.alpha2(1:len_ml,ifile)=vartmp(26,1:len_ml);
    var100_ml.corr(1:len_ml,ifile)=vartmp(27,1:len_ml);
    var100_ml.abscorr(1:len_ml,ifile)=vartmp(28,1:len_ml);
end

% sum the data in time
tot_dist.num    = sum(var_dist.num,2);
tot_dist.numer  = sum(var_dist.numer,2);
tot_dist.denom  = sum(var_dist.denom,2);
tot_dist.alpha  = sum(var_dist.alpha,2);
tot_dist.alpha2 = sum(var_dist.alpha2,2);
tot_dist.corr   = sum(var_dist.corr,2);
tot_dist.abscorr = sum(var_dist.abscorr,2);

tot100_dist.num    = sum(var100_dist.num,2);
tot100_dist.numer  = sum(var100_dist.numer,2);
tot100_dist.denom  = sum(var100_dist.denom,2);
tot100_dist.alpha  = sum(var100_dist.alpha,2);
tot100_dist.alpha2 = sum(var100_dist.alpha2,2);
tot100_dist.corr   = sum(var100_dist.corr,2);
tot100_dist.abscorr = sum(var100_dist.abscorr,2);

tot_ml.num    = sum(var_ml.num,2);
tot_ml.numer  = sum(var_ml.numer,2);
tot_ml.denom  = sum(var_ml.denom,2);
tot_ml.alpha  = sum(var_ml.alpha,2);
tot_ml.alpha2 = sum(var_ml.alpha2,2);
tot_ml.corr   = sum(var_ml.corr,2);
tot_ml.abscorr = sum(var_ml.abscorr,2);

tot100_ml.num    = sum(var100_ml.num,2);
tot100_ml.numer  = sum(var100_ml.numer,2);
tot100_ml.denom  = sum(var100_ml.denom,2);
tot100_ml.alpha  = sum(var100_ml.alpha,2);
tot100_ml.alpha2 = sum(var100_ml.alpha2,2);
tot100_ml.corr   = sum(var100_ml.corr,2);
tot100_ml.abscorr = sum(var100_ml.abscorr,2);


% .jeff is the alpha defined in Jeff's EOL paper
% .beta is the confidence value of alpha
% .avgcorr is the average of absolute correlation

% for physical distance computations
% compute alpha and beta for total samples (sum in all analysis times)
for i = 1:len_dist
    tot_dist.jeff(i)=tot_dist.numer(i)/tot_dist.denom(i);
    tot_dist.beta(i)=(tot_dist.alpha(i)^2/tot_dist.alpha2(i)-1.0)/(tot_dist.num(i)-1.0);
    tot_dist.avgcorr(i)=tot_dist.corr(i)/tot_dist.num(i);
    tot_dist.avgabscorr(i)=tot_dist.abscorr(i)/tot_dist.num(i);
 
    tot100_dist.jeff(i)=tot100_dist.numer(i)/tot100_dist.denom(i);
    tot100_dist.beta(i)=(tot100_dist.alpha(i)^2/tot100_dist.alpha2(i)-1.0)/(tot100_dist.num(i)-1.0);
    tot100_dist.avgcorr(i)=tot100_dist.corr(i)/tot100_dist.num(i);
    tot100_dist.avgabscorr(i)=tot100_dist.abscorr(i)/tot100_dist.num(i);
end

% compute alpha and beta for each analysis time
for j = 1:totfile
    for i = 1:len_dist
        var_dist.jeff(i,j)=var_dist.numer(i,j)/var_dist.denom(i,j);
        var_dist.beta(i,j)=(var_dist.alpha(i,j)^2/var_dist.alpha2(i,j)-1.0)/(var_dist.num(i,j)-1.0);
        var_dist.avgcorr(i,j)=var_dist.corr(i,j)/var_dist.num(i,j);
        var_dist.avgabscorr(i,j)=var_dist.abscorr(i,j)/var_dist.num(i,j);       
 
        var100_dist.jeff(i,j)=var100_dist.numer(i,j)/var100_dist.denom(i,j);
        var100_dist.beta(i,j)=(var100_dist.alpha(i,j)^2/var100_dist.alpha2(i,j)-1.0)/(var100_dist.num(i,j)-1.0);
        var100_dist.avgcorr(i,j)=var100_dist.corr(i,j)/var100_dist.num(i,j);
        var100_dist.avgabscorr(i,j)=var100_dist.abscorr(i,j)/var100_dist.num(i,j);
    end
end

% compute the mean(.jeffmean) and confidence value(.jeffbeta) of the alpha at each analysis time
% also compute the mean of (avgcorr) at each analysis time
for i = 1:len_dist
    temp(1:totfile)=var_dist.jeff(i,1:totfile);
    var_dist.jeffmean(i)=mean(temp);
    var_dist.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var_dist.jeff(i,j);
        temp2=temp2+var_dist.jeff(i,j)^2;
    end
    var_dist.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var_dist.avgcorr(i,1:totfile);
    var_dist.avgcorrmean(i)=mean(temp);

    temp(1:totfile)=var_dist.avgabscorr(i,1:totfile);
    var_dist.avgabscorrmean(i)=mean(temp);
    
    
    temp(1:totfile)=var100_dist.jeff(i,1:totfile);
    var100_dist.jeffmean(i)=mean(temp);
    var100_dist.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var100_dist.jeff(i,j);
        temp2=temp2+var100_dist.jeff(i,j)^2;
    end
    var100_dist.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var100_dist.avgcorr(i,1:totfile);
    var100_dist.avgcorrmean(i)=mean(temp);

    temp(1:totfile)=var100_dist.avgabscorr(i,1:totfile);
    var100_dist.avgabscorrmean(i)=mean(temp);
end


% for model level computations
% all computations are the same as in physical distance
for i = 1:len_ml
    tot_ml.jeff(i)=tot_ml.numer(i)/tot_ml.denom(i);
    tot_ml.beta(i)=(tot_ml.alpha(i)^2/tot_ml.alpha2(i)-1.0)/(tot_ml.num(i)-1.0);
    tot_ml.avgcorr(i)=tot_ml.corr(i)/tot_ml.num(i);
    tot_ml.avgabscorr(i)=tot_ml.abscorr(i)/tot_ml.num(i);
 
    tot100_ml.jeff(i)=tot100_ml.numer(i)/tot100_ml.denom(i);
    tot100_ml.beta(i)=(tot100_ml.alpha(i)^2/tot100_ml.alpha2(i)-1.0)/(tot100_ml.num(i)-1.0);
    tot100_ml.avgcorr(i)=tot100_ml.corr(i)/tot100_ml.num(i);
    tot100_ml.avgabscorr(i)=tot100_ml.abscorr(i)/tot100_ml.num(i);
end

for j = 1:totfile
    for i = 1:len_ml
        var_ml.jeff(i,j)=var_ml.numer(i,j)/var_ml.denom(i,j);
        var_ml.beta(i,j)=(var_ml.alpha(i,j)^2/var_ml.alpha2(i,j)-1.0)/(var_ml.num(i,j)-1.0);
        var_ml.avgcorr(i,j)=var_ml.corr(i,j)/var_ml.num(i,j);
        var_ml.avgabscorr(i,j)=var_ml.abscorr(i,j)/var_ml.num(i,j);
 
        var100_ml.jeff(i,j)=var100_ml.numer(i,j)/var100_ml.denom(i,j);
        var100_ml.beta(i,j)=(var100_ml.alpha(i,j)^2/var100_ml.alpha2(i,j)-1.0)/(var100_ml.num(i,j)-1.0);
        var100_ml.avgcorr(i,j)=var100_ml.corr(i,j)/var100_ml.num(i,j);
        var100_ml.avgabscorr(i,j)=var100_ml.abscorr(i,j)/var100_ml.num(i,j);
    end
end

for i = 1:len_ml
    temp(1:totfile)=var_ml.jeff(i,1:totfile);
    var_ml.jeffmean(i)=mean(temp);
    var_ml.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var_ml.jeff(i,j);
        temp2=temp2+var_ml.jeff(i,j)^2;
    end
    var_ml.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var_ml.avgcorr(i,1:totfile);
    var_ml.avgcorrmean(i)=mean(temp);

    temp(1:totfile)=var_ml.avgabscorr(i,1:totfile);
    var_ml.avgabscorrmean(i)=mean(temp);
    
    
    temp(1:totfile)=var100_ml.jeff(i,1:totfile);
    var100_ml.jeffmean(i)=mean(temp);
    var100_ml.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var100_ml.jeff(i,j);
        temp2=temp2+var100_ml.jeff(i,j)^2;
    end
    var100_ml.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var100_ml.avgcorr(i,1:totfile);
    var100_ml.avgcorrmean(i)=mean(temp);

    temp(1:totfile)=var100_ml.avgabscorr(i,1:totfile);
    var100_ml.avgabscorrmean(i)=mean(temp);
end

% plotting
x_dist=[-40:dlen:40];
x_ml=[1:1:49];

% psfname = sprintf('%s_%s_d0%d.ps',obstypename{idobs},statename{inddim},id);
% disp(sprintf('Removing %s from current directory',psfname))
% system(sprintf('rm %s',psfname));

figure(1);
[ax,h1,h2]=plotyy(x_dist,tot_dist.jeff,x_dist,tot_dist.num);
set(h1,'color','b');
set(h1,'LineWidth',2.0);
set(h2,'Linestyle','o');
hold on;
plot(x_dist,tot_dist.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist,tot_dist.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
% print('-dpsc',psfname,'-append')
% clf;

figure(2);
plot(x_dist,var_dist.jeffmean,'b-');
hold on;
plot(x_dist,var_dist.jeffbeta,'r-.');
hold on;
plot(x_dist,var_dist.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
% print('-dpsc',psfname,'-append')
% clf;

figure(3);
[ax,h1,h2]=plotyy(x_ml,tot_ml.jeff,x_ml,tot_ml.num);
set(h1,'color','b');
set(h1,'LineWidth',2.0);
set(h2,'Linestyle','o');
% set(h2,'color','b');
% plot(tot_ml.jeff,'b-');
hold on;
plot(tot_ml.beta,'r-','LineWidth',2.0);
hold on;
plot(tot_ml.avgcorr,'k--','LineWidth',2.0);
legend('alpha','beta','abs corr');
hold off;
xlabel('model levels','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
% print('-dpsc',psfname,'-append')
% clf;

figure(4);
plot(var_ml.jeffmean,'b-');
hold on;
plot(var_ml.jeffbeta,'r-.');
hold on;
plot(var_ml.avgcorrmean,'k--');
legend('alpha','beta','abs corr');
hold off;
xlabel('model levels','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14);
% print('-dpsc',psfname,'-append')
% clf;


end   % end inddim
end   % end idobs
end   % end id

