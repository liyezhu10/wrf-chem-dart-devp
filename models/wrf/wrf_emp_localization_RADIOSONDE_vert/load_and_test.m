
totfile = 1;        %46;
% inputname = {'U10_U'};
group_size  = 2;

%len_dist = 2000;
%dlen     = 0.0002;
len_dist   = 49;

% for plotting use
earth_radius = 6370;              % unit km
cen_lat      = 40.8491;           % center latitude
cen_lat_rad  = cen_lat*pi/180.0;
R_cen_lat    = earth_radius*cos(cen_lat_rad);
% x_ml=[1:1:49];

% process data: bin to larger subsets
%bin = 10;
%len_dist_bin = len_dist/bin+1;

wrf_ndom = 1;

obstypename = {'RADIOSONDE_PRESSURE',...
               'RADIOSONDE_TEMPERATURE',...
               'RADIOSONDE_U_WIND_COMPONENT',...
               'RADIOSONDE_V_WIND_COMPONENT',...
              };
nobstype = size(obstypename,2);

statename = {'U10','V10','T2','Q2','TH2','MU','PSFC',...   % 2D variables
             'U','V','W','T','PH',...
             'QVAPOR','QCLOUD','QNRAIN','QNICE',...
             'QRAIN','QICE','QSNOW','QGRAUP'};
nstatetype = size(statename,2);

obslevelname = {'obslev1'};
varlevelname = {'varlev1'};

timename={'2010122906','2010122909','2010122912','2010122915','2010122918','2010122921','2010123000',...
          '2010123003','2010123006','2010123009','2010123012','2010123015','2010123018','2010123021','2010123100',...
          '2010123103','2010123106','2010123109','2010123112','2010123115','2010123118','2010123121','2011010100',...
          '2011010103','2011010106','2011010109','2011010112','2011010115','2011010118','2011010121','2011010200',...
          '2011010203','2011010206','2011010209','2011010112','2011010215','2011010218','2011010221','2011010300',...
          '2011010303','2011010306','2011010309','2011010312','2011010315','2011010318','2011010321'};

      
% % initialize variable;
% tot_dist.jeff


for id = 1:wrf_ndom
for idobs = 2:2
for inddim = 11:11
%for inddim = 8:20
    
%% load in data
for ifile = 1:totfile
%   vartmp = load(sprintf('%s_%s_d0%d_%s_%s_%s',obstypename{idobs},statename{inddim},id,obslevelname{1},varlevelname{1},timename{ifile}));
    vartmp = load('RADIOSONDE_TEMPERATURE_T_ocean_d01_obslev1');
%    vartmp = load('RADIOSONDE_U_WIND_COMPONENT_U_d01_obslev1_varlev1_2010122906');

    var_dist.num(1:len_dist,ifile)    = vartmp(1,1:len_dist);
    var_dist.beta(1:len_dist,ifile)   = vartmp(2,1:len_dist);
    var_dist.beta2(1:len_dist,ifile)  = vartmp(3,1:len_dist);

    for ig = 1:group_size
       igs = 3+(ig-1)*9+1;   % 9 - number of computations for each group (sumx, sumy, ...)
       var_dist.sumx(1:len_dist,ig,ifile)   = vartmp(igs,1:len_dist);
       var_dist.sumy(1:len_dist,ig,ifile)   = vartmp(igs+1,1:len_dist);
       var_dist.sumx2(1:len_dist,ig,ifile)  = vartmp(igs+2,1:len_dist);
       var_dist.sumy2(1:len_dist,ig,ifile)  = vartmp(igs+3,1:len_dist);
       var_dist.numer(1:len_dist,ig,ifile)  = vartmp(igs+4,1:len_dist);
       var_dist.denom(1:len_dist,ig,ifile)  = vartmp(igs+5,1:len_dist);
       var_dist.alpha(1:len_dist,ig,ifile)  = vartmp(igs+6,1:len_dist);
       var_dist.alpha2(1:len_dist,ig,ifile) = vartmp(igs+7,1:len_dist);
       var_dist.corr(1:len_dist,ig,ifile)   = vartmp(igs+8,1:len_dist);
    end

end


%% sum the data in time
tot_dist.num    = sum(var_dist.num,2);
tot_dist.beta   = sum(var_dist.beta,2);
tot_dist.beta2  = sum(var_dist.beta2,2);

tot_dist.sumx   = sum(var_dist.sumx,3);
tot_dist.sumy   = sum(var_dist.sumy,3);
tot_dist.sumx2  = sum(var_dist.sumx2,3);
tot_dist.sumy2  = sum(var_dist.sumy2,3);
tot_dist.numer  = sum(var_dist.numer,3);
tot_dist.denom  = sum(var_dist.denom,3);
tot_dist.alpha  = sum(var_dist.alpha,3);
tot_dist.alpha2 = sum(var_dist.alpha2,3);
tot_dist.corr   = sum(var_dist.corr,3);


for i = 1:len_dist
    beta(i) = (tot_dist.beta(i)/tot_dist.beta2(i) - 1.0)/(group_size-1);

    for ig = 1:group_size
        alpha_group(i,ig) = tot_dist.numer(i,ig)/tot_dist.denom(i,ig);
    end

end



end   % end inddim
end   % end idobs
end   % end id

