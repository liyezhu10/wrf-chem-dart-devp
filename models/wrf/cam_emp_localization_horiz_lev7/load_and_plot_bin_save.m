
totfile = 46;
% inputname = {'U10_U'};
len_dist = 2000;
dlen     = 0.0002;
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

wrf_ndom = 1;

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

statename = {'U10','V10','T2','Q2','TH2','MU','PSFC',...   % 2D variables
             'U','V','W','T','PH'...
             'QVAPOR','QCLOUD','QRAIN','QNRAIN',...
             'QICE','QNICE','QSNOW','QGRAUP'};
nstatetype = size(statename,2);

timename={'2010122906','2010122909','2010122912','2010122915','2010122918','2010122921','2010123000',...
          '2010123003','2010123006','2010123009','2010123012','2010123015','2010123018','2010123021','2010123100',...
          '2010123103','2010123106','2010123109','2010123112','2010123115','2010123118','2010123121','2011010100',...
          '2011010103','2011010106','2011010109','2011010112','2011010115','2011010118','2011010121','2011010200',...
          '2011010203','2011010206','2011010209','2011010112','2011010215','2011010218','2011010221','2011010300',...
          '2011010303','2011010306','2011010309','2011010312','2011010315','2011010318','2011010321'};



for id = 1:wrf_ndom
for idobs = 2:2
for inddim = 6:20
    
%% load in data
for ifile = 1:totfile
   string = sprintf('load %s_%s_d0%d_%s',obstypename{idobs},statename{inddim},id,timename{ifile});
   eval(string);
   string = sprintf('vartmp = %s_%s_d0%d_%s;',obstypename{idobs},statename{inddim},id,timename{ifile});
   eval(string);

%      string = sprintf('load %s_%s_d0%d',obstypename{idobs},statename{inddim},id);
%      eval(string);
%      string = sprintf('vartmp = %s_%s_d0%d;',obstypename{idobs},statename{inddim},id);
%      eval(string);

    var_dist.num(1:len_dist,ifile)=vartmp(1,1:len_dist);
    var_dist.numer(1:len_dist,ifile)=vartmp(2,1:len_dist);
    var_dist.denom(1:len_dist,ifile)=vartmp(3,1:len_dist);
    var_dist.alpha(1:len_dist,ifile)=vartmp(4,1:len_dist);
    var_dist.alpha2(1:len_dist,ifile)=vartmp(5,1:len_dist);
    var_dist.corr(1:len_dist,ifile)=vartmp(6,1:len_dist);
    
    var100_dist.num(1:len_dist,ifile)=vartmp(7,1:len_dist);
    var100_dist.numer(1:len_dist,ifile)=vartmp(8,1:len_dist);
    var100_dist.denom(1:len_dist,ifile)=vartmp(9,1:len_dist);
    var100_dist.alpha(1:len_dist,ifile)=vartmp(10,1:len_dist);
    var100_dist.alpha2(1:len_dist,ifile)=vartmp(11,1:len_dist);
    var100_dist.corr(1:len_dist,ifile)=vartmp(12,1:len_dist);
    
    % subtract var100 from var => grid points with terrain hihger than 100m
    var100verse_dist.num(1:len_dist,ifile)=var_dist.num(1:len_dist,ifile)-var100_dist.num(1:len_dist,ifile);
    var100verse_dist.numer(1:len_dist,ifile)=var_dist.numer(1:len_dist,ifile)-var100_dist.numer(1:len_dist,ifile);
    var100verse_dist.denom(1:len_dist,ifile)=var_dist.denom(1:len_dist,ifile)-var100_dist.denom(1:len_dist,ifile);
    var100verse_dist.alpha(1:len_dist,ifile)=var_dist.alpha(1:len_dist,ifile)-var100_dist.alpha(1:len_dist,ifile);
    var100verse_dist.alpha2(1:len_dist,ifile)=var_dist.alpha2(1:len_dist,ifile)-var100_dist.alpha2(1:len_dist,ifile);
    var100verse_dist.corr(1:len_dist,ifile)=var_dist.corr(1:len_dist,ifile)-var100_dist.corr(1:len_dist,ifile);
end


%% sum the data in time
tot_dist.num    = sum(var_dist.num,2);
tot_dist.numer  = sum(var_dist.numer,2);
tot_dist.denom  = sum(var_dist.denom,2);
tot_dist.alpha  = sum(var_dist.alpha,2);
tot_dist.alpha2 = sum(var_dist.alpha2,2);
tot_dist.corr   = sum(var_dist.corr,2);

tot100_dist.num    = sum(var100_dist.num,2);
tot100_dist.numer  = sum(var100_dist.numer,2);
tot100_dist.denom  = sum(var100_dist.denom,2);
tot100_dist.alpha  = sum(var100_dist.alpha,2);
tot100_dist.alpha2 = sum(var100_dist.alpha2,2);
tot100_dist.corr   = sum(var100_dist.corr,2);

tot100verse_dist.num    = sum(var100verse_dist.num,2);
tot100verse_dist.numer  = sum(var100verse_dist.numer,2);
tot100verse_dist.denom  = sum(var100verse_dist.denom,2);
tot100verse_dist.alpha  = sum(var100verse_dist.alpha,2);
tot100verse_dist.alpha2 = sum(var100verse_dist.alpha2,2);
tot100verse_dist.corr   = sum(var100verse_dist.corr,2);

%% bin the data
var_dist_bin.num(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.numer(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.denom(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.alpha(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.alpha2(1:len_dist_bin,1:totfile) = 0.0;
var_dist_bin.corr(1:len_dist_bin,1:totfile) = 0.0;

var100_dist_bin.num(1:len_dist_bin,1:totfile) = 0.0;
var100_dist_bin.numer(1:len_dist_bin,1:totfile) = 0.0;
var100_dist_bin.denom(1:len_dist_bin,1:totfile) = 0.0;
var100_dist_bin.alpha(1:len_dist_bin,1:totfile) = 0.0;
var100_dist_bin.alpha2(1:len_dist_bin,1:totfile) = 0.0;
var100_dist_bin.corr(1:len_dist_bin,1:totfile) = 0.0;

var100verse_dist_bin.num(1:len_dist_bin,1:totfile) = 0.0;
var100verse_dist_bin.numer(1:len_dist_bin,1:totfile) = 0.0;
var100verse_dist_bin.denom(1:len_dist_bin,1:totfile) = 0.0;
var100verse_dist_bin.alpha(1:len_dist_bin,1:totfile) = 0.0;
var100verse_dist_bin.alpha2(1:len_dist_bin,1:totfile) = 0.0;
var100verse_dist_bin.corr(1:len_dist_bin,1:totfile) = 0.0;

tot_dist_bin.num(1:len_dist_bin) = 0.0;
tot_dist_bin.numer(1:len_dist_bin) = 0.0;
tot_dist_bin.denom(1:len_dist_bin) = 0.0;
tot_dist_bin.alpha(1:len_dist_bin) = 0.0;
tot_dist_bin.alpha2(1:len_dist_bin) = 0.0;
tot_dist_bin.corr(1:len_dist_bin) = 0.0;

tot100_dist_bin.num(1:len_dist_bin) = 0.0;
tot100_dist_bin.numer(1:len_dist_bin) = 0.0;
tot100_dist_bin.denom(1:len_dist_bin) = 0.0;
tot100_dist_bin.alpha(1:len_dist_bin) = 0.0;
tot100_dist_bin.alpha2(1:len_dist_bin) = 0.0;
tot100_dist_bin.corr(1:len_dist_bin) = 0.0;

tot100verse_dist_bin.num(1:len_dist_bin) = 0.0;
tot100verse_dist_bin.numer(1:len_dist_bin) = 0.0;
tot100verse_dist_bin.denom(1:len_dist_bin) = 0.0;
tot100verse_dist_bin.alpha(1:len_dist_bin) = 0.0;
tot100verse_dist_bin.alpha2(1:len_dist_bin) = 0.0;
tot100verse_dist_bin.corr(1:len_dist_bin) = 0.0;

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
        
        var100_dist_bin.num(i,k)    = var100_dist_bin.num(i,k) + var100_dist.num(j,k);
        var100_dist_bin.numer(i,k)  = var100_dist_bin.numer(i,k) + var100_dist.numer(j,k);
        var100_dist_bin.denom(i,k)  = var100_dist_bin.denom(i,k) + var100_dist.denom(j,k);
        var100_dist_bin.alpha(i,k)  = var100_dist_bin.alpha(i,k) + var100_dist.alpha(j,k);
        var100_dist_bin.alpha2(i,k) = var100_dist_bin.alpha2(i,k) + var100_dist.alpha2(j,k);
        var100_dist_bin.corr(i,k)   = var100_dist_bin.corr(i,k) + var100_dist.corr(j,k);
        
        var100verse_dist_bin.num(i,k)    = var100verse_dist_bin.num(i,k) + var100verse_dist.num(j,k);
        var100verse_dist_bin.numer(i,k)  = var100verse_dist_bin.numer(i,k) + var100verse_dist.numer(j,k);
        var100verse_dist_bin.denom(i,k)  = var100verse_dist_bin.denom(i,k) + var100verse_dist.denom(j,k);
        var100verse_dist_bin.alpha(i,k)  = var100verse_dist_bin.alpha(i,k) + var100verse_dist.alpha(j,k);
        var100verse_dist_bin.alpha2(i,k) = var100verse_dist_bin.alpha2(i,k) + var100verse_dist.alpha2(j,k);
        var100verse_dist_bin.corr(i,k)   = var100verse_dist_bin.corr(i,k) + var100verse_dist.corr(j,k);
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
        
        tot100_dist_bin.num(i)    = tot100_dist_bin.num(i) + tot100_dist.num(j);
        tot100_dist_bin.numer(i)  = tot100_dist_bin.numer(i) + tot100_dist.numer(j);
        tot100_dist_bin.denom(i)  = tot100_dist_bin.denom(i) + tot100_dist.denom(j);
        tot100_dist_bin.alpha(i)  = tot100_dist_bin.alpha(i) + tot100_dist.alpha(j);
        tot100_dist_bin.alpha2(i) = tot100_dist_bin.alpha2(i) + tot100_dist.alpha2(j);
        tot100_dist_bin.corr(i)   = tot100_dist_bin.corr(i) + tot100_dist.corr(j);
        
        tot100verse_dist_bin.num(i)    = tot100verse_dist_bin.num(i) + tot100verse_dist.num(j);
        tot100verse_dist_bin.numer(i)  = tot100verse_dist_bin.numer(i) + tot100verse_dist.numer(j);
        tot100verse_dist_bin.denom(i)  = tot100verse_dist_bin.denom(i) + tot100verse_dist.denom(j);
        tot100verse_dist_bin.alpha(i)  = tot100verse_dist_bin.alpha(i) + tot100verse_dist.alpha(j);
        tot100verse_dist_bin.alpha2(i) = tot100verse_dist_bin.alpha2(i) + tot100verse_dist.alpha2(j);
        tot100verse_dist_bin.corr(i)   = tot100verse_dist_bin.corr(i) + tot100verse_dist.corr(j);
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
    
    tot100_dist.jeff(i)=tot100_dist.numer(i)/tot100_dist.denom(i);
    tot100_dist.beta(i)=(tot100_dist.alpha(i)^2/tot100_dist.alpha2(i)-1.0)/(tot100_dist.num(i)-1.0);
    tot100_dist.avgcorr(i)=tot100_dist.corr(i)/tot100_dist.num(i);
    
    tot100verse_dist.jeff(i)=tot100verse_dist.numer(i)/tot100verse_dist.denom(i);
    tot100verse_dist.beta(i)=(tot100verse_dist.alpha(i)^2/tot100verse_dist.alpha2(i)-1.0)/(tot100verse_dist.num(i)-1.0);
    tot100verse_dist.avgcorr(i)=tot100verse_dist.corr(i)/tot100verse_dist.num(i);
end

% compute alpha and beta for each analysis time
for j = 1:totfile
    for i = 1:len_dist
        var_dist.jeff(i,j)=var_dist.numer(i,j)/var_dist.denom(i,j);
        var_dist.beta(i,j)=(var_dist.alpha(i,j)^2/var_dist.alpha2(i,j)-1.0)/(var_dist.num(i,j)-1.0);
        var_dist.avgcorr(i,j)=var_dist.corr(i,j)/var_dist.num(i,j);
        
        var100_dist.jeff(i,j)=var100_dist.numer(i,j)/var100_dist.denom(i,j);
        var100_dist.beta(i,j)=(var100_dist.alpha(i,j)^2/var100_dist.alpha2(i,j)-1.0)/(var100_dist.num(i,j)-1.0);
        var100_dist.avgcorr(i,j)=var100_dist.corr(i,j)/var100_dist.num(i,j);
        
        var100verse_dist.jeff(i,j)=var100verse_dist.numer(i,j)/var100verse_dist.denom(i,j);
        var100verse_dist.beta(i,j)=(var100verse_dist.alpha(i,j)^2/var100verse_dist.alpha2(i,j)-1.0)/(var100verse_dist.num(i,j)-1.0);
        var100verse_dist.avgcorr(i,j)=var100verse_dist.corr(i,j)/var100verse_dist.num(i,j);
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
    
    
    temp(1:totfile)=var100verse_dist.jeff(i,1:totfile);
    var100verse_dist.jeffmean(i)=mean(temp);
    var100verse_dist.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var100verse_dist.jeff(i,j);
        temp2=temp2+var100verse_dist.jeff(i,j)^2;
    end
    var100verse_dist.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var100verse_dist.avgcorr(i,1:totfile);
    var100verse_dist.avgcorrmean(i)=mean(temp);
end


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
    
    tot100_dist_bin.jeff(i)=tot100_dist_bin.numer(i)/tot100_dist_bin.denom(i);
    tot100_dist_bin.beta(i)=(tot100_dist_bin.alpha(i)^2/tot100_dist_bin.alpha2(i)-1.0)/(tot100_dist_bin.num(i)-1.0);
    tot100_dist_bin.avgcorr(i)=tot100_dist_bin.corr(i)/tot100_dist_bin.num(i);
    
    tot100verse_dist_bin.jeff(i)=tot100verse_dist_bin.numer(i)/tot100verse_dist_bin.denom(i);
    tot100verse_dist_bin.beta(i)=(tot100verse_dist_bin.alpha(i)^2/tot100verse_dist_bin.alpha2(i)-1.0)/(tot100verse_dist_bin.num(i)-1.0);
    tot100verse_dist_bin.avgcorr(i)=tot100verse_dist_bin.corr(i)/tot100verse_dist_bin.num(i);
end

% compute alpha and beta for each analysis time
for j = 1:totfile
    for i = 1:len_dist_bin
        var_dist_bin.jeff(i,j)=var_dist_bin.numer(i,j)/var_dist_bin.denom(i,j);
        var_dist_bin.beta(i,j)=(var_dist_bin.alpha(i,j)^2/var_dist_bin.alpha2(i,j)-1.0)/(var_dist_bin.num(i,j)-1.0);
        var_dist_bin.avgcorr(i,j)=var_dist_bin.corr(i,j)/var_dist_bin.num(i,j);
        
        var100_dist_bin.jeff(i,j)=var100_dist_bin.numer(i,j)/var100_dist_bin.denom(i,j);
        var100_dist_bin.beta(i,j)=(var100_dist_bin.alpha(i,j)^2/var100_dist_bin.alpha2(i,j)-1.0)/(var100_dist_bin.num(i,j)-1.0);
        var100_dist_bin.avgcorr(i,j)=var100_dist_bin.corr(i,j)/var100_dist_bin.num(i,j);
        
        var100verse_dist_bin.jeff(i,j)=var100verse_dist_bin.numer(i,j)/var100verse_dist_bin.denom(i,j);
        var100verse_dist_bin.beta(i,j)=(var100verse_dist_bin.alpha(i,j)^2/var100verse_dist_bin.alpha2(i,j)-1.0)/(var100verse_dist_bin.num(i,j)-1.0);
        var100verse_dist_bin.avgcorr(i,j)=var100verse_dist_bin.corr(i,j)/var100verse_dist_bin.num(i,j);
    end
end

% compute the mean(.jeffmean) and confidence value(.jeffbeta) of the alpha at each analysis time
% also compute the mean of (avgcorr) at each analysis time
for i = 1:len_dist_bin
    temp(1:totfile)=var_dist_bin.jeff(i,1:totfile);
    var_dist_bin.jeffmean(i)=mean(temp);
    var_dist_bin.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var_dist_bin.jeff(i,j);
        temp2=temp2+var_dist_bin.jeff(i,j)^2;
    end
    var_dist_bin.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var_dist_bin.avgcorr(i,1:totfile);
    var_dist_bin.avgcorrmean(i)=mean(temp);
    
    
    temp(1:totfile)=var100_dist_bin.jeff(i,1:totfile);
    var100_dist_bin.jeffmean(i)=mean(temp);
    var100_dist_bin.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var100_dist_bin.jeff(i,j);
        temp2=temp2+var100_dist_bin.jeff(i,j)^2;
    end
    var100_dist_bin.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var100_dist_bin.avgcorr(i,1:totfile);
    var100_dist_bin.avgcorrmean(i)=mean(temp);
    
    
    temp(1:totfile)=var100verse_dist_bin.jeff(i,1:totfile);
    var100verse_dist_bin.jeffmean(i)=mean(temp);
    var100verse_dist_bin.jeffvar(i)=var(temp);
    
    temp1=0.0;
    temp2=0.0;
    for j=1:totfile
        temp1=temp1+var100verse_dist_bin.jeff(i,j);
        temp2=temp2+var100verse_dist_bin.jeff(i,j)^2;
    end
    var100verse_dist_bin.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
    
    temp(1:totfile)=var100verse_dist_bin.avgcorr(i,1:totfile);
    var100verse_dist_bin.avgcorrmean(i)=mean(temp);
end



%% plot section

x_dist=[0:dlen*R_cen_lat:(len_dist-1)*dlen*R_cen_lat];
x_dist_bin=[0:dlen*R_cen_lat*bin:(len_dist_bin-1)*dlen*R_cen_lat*bin];

% plotting
% x_dist=[0:dlen:len_dist-1];
% x_ml=[1:1:49];

psfname = sprintf('%s_%s_d0%d.ps',obstypename{idobs},statename{inddim},id);
disp(sprintf('Removing %s from current directory',psfname))
system(sprintf('rm %s',psfname));

% 0. plot uncontrained subset
% 0.1. unbbined Jeff computation
% figure(1);
% [ax,h1,h2]=plotyy(x_dist,tot_dist.jeff,x_dist,tot_dist.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
plot(x_dist,tot_dist.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist,tot_dist.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist,tot_dist.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',12,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 0.2. bined Jeff computation
% figure(2);
[ax,h1,h2]=plotyy(x_dist_bin,tot_dist_bin.jeff,x_dist_bin,tot_dist_bin.num);
set(h1,'color','b');
set(h1,'LineWidth',2.0);
set(h2,'Linestyle','o');
% plot(x_dist_bin,tot_dist_bin.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot_dist_bin.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot_dist_bin.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',12,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 0.3. unbinned Lili computation (mean of Jeff, and beta of mean)
% figure(3);
plot(x_dist,var_dist.jeffmean,'b-');
hold on;
plot(x_dist,var_dist.jeffbeta,'r-.');
hold on;
plot(x_dist,var_dist.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 0.4. binned Lili computation
% figure(4);
plot(x_dist_bin,var_dist_bin.jeffmean,'b-');
hold on;
plot(x_dist_bin,var_dist_bin.jeffbeta,'r-.');
hold on;
plot(x_dist_bin,var_dist_bin.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('%s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;


% 1. plot subset with terrain height < 100m
% 1.1 unbinned Jeff computation
% figure(11);
% [ax,h1,h2]=plotyy(x_dist,tot100_dist.jeff,x_dist,tot100_dist.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
plot(x_dist,tot100_dist.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist,tot100_dist.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist,tot100_dist.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER<100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 1.2 binned Jeff computation
% figure(12);
[ax,h1,h2]=plotyy(x_dist_bin,tot100_dist_bin.jeff,x_dist_bin,tot100_dist_bin.num);
set(h1,'color','b');
set(h1,'LineWidth',2.0);
set(h2,'Linestyle','o');
% plot(x_dist_bin,tot100_dist_bin.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot100_dist_bin.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot100_dist_bin.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER<100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 1.3 unbbined Lili computation
% figure(13);
plot(x_dist,var100_dist.jeffmean,'b-');
hold on;
plot(x_dist,var100_dist.jeffbeta,'r-.');
hold on;
plot(x_dist,var100_dist.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER<100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% binned Lili computation
% figure(14);
plot(x_dist_bin,var100_dist_bin.jeffmean,'b-');
hold on;
plot(x_dist_bin,var100_dist_bin.jeffbeta,'r-.');
hold on;
plot(x_dist_bin,var100_dist_bin.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER<100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;


% 2. plot subset with terrain > 100m
% 2.1. unbinned Jeff computation
% figure(21);
% [ax,h1,h2]=plotyy(x_dist,tot100_dist.jeff,x_dist,tot100_dist.num);
% set(h1,'color','b');
% set(h1,'LineWidth',2.0);
% set(h2,'Linestyle','o');
plot(x_dist,tot100verse_dist.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist,tot100verse_dist.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist,tot100verse_dist.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER>100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',12,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 2.2. binned Jeff computation
% figure(22);
[ax,h1,h2]=plotyy(x_dist_bin,tot100verse_dist_bin.jeff,x_dist_bin,tot100verse_dist_bin.num);
set(h1,'color','b');
set(h1,'LineWidth',2.0);
set(h2,'Linestyle','o');
% plot(x_dist_bin,tot100_dist_bin.jeff,'b-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot100verse_dist_bin.beta,'r-','LineWidth',2.0);
hold on;
plot(x_dist_bin,tot100verse_dist_bin.avgcorr,'k--','LineWidth',2.0);
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER>100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',12,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 2.3. unbinned Lili computation
% figure(23);
plot(x_dist,var100verse_dist.jeffmean,'b-');
hold on;
plot(x_dist,var100verse_dist.jeffbeta,'r-.');
hold on;
plot(x_dist,var100verse_dist.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER>100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;

% 2.4. binned Lili computation
% figure(24);
plot(x_dist_bin,var100verse_dist_bin.jeffmean,'b-');
hold on;
plot(x_dist_bin,var100verse_dist_bin.jeffbeta,'r-.');
hold on;
plot(x_dist_bin,var100verse_dist_bin.avgcorrmean,'k--');
hold off;
legend('alpha','beta','abs corr');
% ylim([-1 1.5]);
hold off;
xlabel('distance (unit: km)','FontSize',14);
title(sprintf('TER>100m %s %s d0%d',obstypename{idobs},statename{inddim},id),'FontSize',14,'interpreter','none');
print('-dpsc',psfname,'-append')
clf;


end   % end inddim
end   % end idobs
end   % end id

