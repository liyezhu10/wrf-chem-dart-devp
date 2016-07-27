clear;
totfile = 1;
inputname = {'SFC_ALTIMETER_T_d01'};
len_dist = 81;
dlen     = 1;
len_ml   = 49;

timename={'2010122906','2010122909','2010122912','2010122915','2010122918','2010122921','2010123000',...
          '2010123003','2010123006','2010123009','2010123012','2010123015','2010123018','2010123021','2010123100',...
          '2010123103','2010123106','2010123109','2010123112','2010123115','2010123118','2010123121','2011010100',...
          '2011010103','2011010106','2011010109','2011010112','2011010115','2011010118','2011010121','2011010200',...
          '2011010203','2011010206','2011010209','2011010112','2011010215','2011010218','2011010221','2011010300',...
          '2011010303','2011010306','2011010309','2011010312','2011010315','2011010318','2011010321'};

% load in data
for ifile = 1:totfile
    string = sprintf('load %s',inputname{1});
    eval(string);
    string = sprintf('vartmp = %s;',inputname{1});
    eval(string);
    
    var_dist.num(1:len_dist,ifile)=vartmp(1,1:len_dist);
    var_dist.numer(1:len_dist,ifile)=vartmp(2,1:len_dist);
    var_dist.denom(1:len_dist,ifile)=vartmp(3,1:len_dist);
    var_dist.alpha(1:len_dist,ifile)=vartmp(4,1:len_dist);
    var_dist.alpha2(1:len_dist,ifile)=vartmp(5,1:len_dist);
    
    var_ml.num(1:len_ml,ifile)=vartmp(6,1:len_ml);
    var_ml.numer(1:len_ml,ifile)=vartmp(7,1:len_ml);
    var_ml.denom(1:len_ml,ifile)=vartmp(8,1:len_ml);
    var_ml.alpha(1:len_ml,ifile)=vartmp(9,1:len_ml);
    var_ml.alpha2(1:len_ml,ifile)=vartmp(10,1:len_ml);
    
    var100_dist.num(1:len_dist,ifile)=vartmp(11,1:len_dist);
    var100_dist.numer(1:len_dist,ifile)=vartmp(12,1:len_dist);
    var100_dist.denom(1:len_dist,ifile)=vartmp(13,1:len_dist);
    var100_dist.alpha(1:len_dist,ifile)=vartmp(14,1:len_dist);
    var100_dist.alpha2(1:len_dist,ifile)=vartmp(15,1:len_dist);
    
    var100_ml.num(1:len_ml,ifile)=vartmp(16,1:len_ml);
    var100_ml.numer(1:len_ml,ifile)=vartmp(17,1:len_ml);
    var100_ml.denom(1:len_ml,ifile)=vartmp(18,1:len_ml);
    var100_ml.alpha(1:len_ml,ifile)=vartmp(19,1:len_ml);
    var100_ml.alpha2(1:len_ml,ifile)=vartmp(20,1:len_ml);
end

% sum the data in time
tot_dist.num(1:len_dist)=0.0;
tot_dist.numer(1:len_dist)=0.0;
tot_dist.denom(1:len_dist)=0.0;
tot_dist.alpha(1:len_dist)=0.0;
tot_dist.alpha2(1:len_dist)=0.0;

tot_ml.num(1:len_ml)=0.0;
tot_ml.numer(1:len_ml)=0.0;
tot_ml.denom(1:len_ml)=0.0;
tot_ml.alpha(1:len_ml)=0.0;
tot_ml.alpha2(1:len_ml)=0.0;

tot100_dist.num(1:len_dist)=0.0;
tot100_dist.numer(1:len_dist)=0.0;
tot100_dist.denom(1:len_dist)=0.0;
tot100_dist.alpha(1:len_dist)=0.0;
tot100_dist.alpha2(1:len_dist)=0.0;

tot100_ml.num(1:len_ml)=0.0;
tot100_ml.numer(1:len_ml)=0.0;
tot100_ml.denom(1:len_ml)=0.0;
tot100_ml.alpha(1:len_ml)=0.0;
tot100_ml.alpha2(1:len_ml)=0.0;

for i = 1:totfile
    temp_dist(1:len_dist)=var_dist.num(1:len_dist,i);
    tot_dist.num(1:len_dist)=tot_dist.num(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var_dist.numer(1:len_dist,i);
    tot_dist.numer(1:len_dist)=tot_dist.numer(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var_dist.denom(1:len_dist,i);
    tot_dist.denom(1:len_dist)=tot_dist.denom(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var_dist.alpha(1:len_dist,i);
    tot_dist.alpha(1:len_dist)=tot_dist.alpha(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var_dist.alpha2(1:len_dist,i);
    tot_dist.alpha2(1:len_dist)=tot_dist.alpha2(1:len_dist)+temp_dist(1:len_dist);
    
    temp_ml(1:len_ml)=var_ml.num(1:len_ml,i);
    tot_ml.num(1:len_ml)=tot_ml.num(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var_ml.numer(1:len_ml,i);
    tot_ml.numer(1:len_ml)=tot_ml.numer(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var_ml.denom(1:len_ml,i);
    tot_ml.denom(1:len_ml)=tot_ml.denom(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var_ml.alpha(1:len_ml,i);
    tot_ml.alpha(1:len_ml)=tot_ml.alpha(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var_ml.alpha2(1:len_ml,i);
    tot_ml.alpha2(1:len_ml)=tot_ml.alpha2(1:len_ml)+temp_ml(1:len_ml);
    
    temp_dist(1:len_dist)=var100_dist.num(1:len_dist,i);
    tot100_dist.num(1:len_dist)=tot100_dist.num(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var100_dist.numer(1:len_dist,i);
    tot100_dist.numer(1:len_dist)=tot100_dist.numer(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var100_dist.denom(1:len_dist,i);
    tot100_dist.denom(1:len_dist)=tot100_dist.denom(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var100_dist.alpha(1:len_dist,i);
    tot100_dist.alpha(1:len_dist)=tot100_dist.alpha(1:len_dist)+temp_dist(1:len_dist);
    temp_dist(1:len_dist)=var100_dist.alpha2(1:len_dist,i);
    tot100_dist.alpha2(1:len_dist)=tot100_dist.alpha2(1:len_dist)+temp_dist(1:len_dist);
    
    temp_ml(1:len_ml)=var100_ml.num(1:len_ml,i);
    tot100_ml.num(1:len_ml)=tot100_ml.num(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var100_ml.numer(1:len_ml,i);
    tot100_ml.numer(1:len_ml)=tot100_ml.numer(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var100_ml.denom(1:len_ml,i);
    tot100_ml.denom(1:len_ml)=tot100_ml.denom(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var100_ml.alpha(1:len_ml,i);
    tot100_ml.alpha(1:len_ml)=tot100_ml.alpha(1:len_ml)+temp_ml(1:len_ml);
    temp_ml(1:len_ml)=var100_ml.alpha2(1:len_ml,i);
    tot100_ml.alpha2(1:len_ml)=tot100_ml.alpha2(1:len_ml)+temp_ml(1:len_ml);
end


% compute alpha and beta
% .jeff is the alpha defined in Jeff's EOL paper
% .beta is the confidence value of alpha
for i = 1:len_dist
    tot_dist.jeff(i)=tot_dist.numer(i)/tot_dist.denom(i);
    tot_dist.beta(i)=(tot_dist.alpha(i)^2/tot_dist.alpha2(i)-1.0)/(tot_dist.num(i)-1.0);
    
    tot100_dist.jeff(i)=tot100_dist.numer(i)/tot100_dist.denom(i);
    tot100_dist.beta(i)=(tot100_dist.alpha(i)^2/tot100_dist.alpha2(i)-1.0)/(tot100_dist.num(i)-1.0);
end

% for j = 1:totfile
%     for i = 1:len_dist
%         var_dist.jeff(i,j)=var_dist.numer(i,j)/var_dist.denom(i,j);
%         var_dist.beta(i,j)=(var_dist.alpha(i,j)^2/var_dist.alpha2(i,j)-1.0)/(var_dist.num(i,j)-1.0);
%         
%         var100_dist.jeff(i,j)=var100_dist.numer(i,j)/var100_dist.denom(i,j);
%         var100_dist.beta(i,j)=(var100_dist.alpha(i,j)^2/var100_dist.alpha2(i,j)-1.0)/(var100_dist.num(i,j)-1.0);
%     end
% end

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
%     
%     temp(1:totfile)=var100_dist.jeff(i,1:totfile);
%     var100_dist.jeffmean(i)=mean(temp);
%     var100_dist.jeffvar(i)=var(temp);
%     
%     temp1=0.0;
%     temp2=0.0;
%     for j=1:totfile
%         temp1=temp1+var100_dist.jeff(i,j);
%         temp2=temp2+var100_dist.jeff(i,j)^2;
%     end
%     var100_dist.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
% end



for i = 1:len_ml
    tot_ml.jeff(i)=tot_ml.numer(i)/tot_ml.denom(i);
    tot_ml.beta(i)=(tot_ml.alpha(i)^2/tot_ml.alpha2(i)-1.0)/(tot_ml.num(i)-1.0);
    
    tot100_ml.jeff(i)=tot100_ml.numer(i)/tot100_ml.denom(i);
    tot100_ml.beta(i)=(tot100_ml.alpha(i)^2/tot100_ml.alpha2(i)-1.0)/(tot100_ml.num(i)-1.0);
end

% for j = 1:totfile
%     for i = 1:len_ml
%         var_ml.jeff(i,j)=var_ml.numer(i,j)/var_ml.denom(i,j);
%         var_ml.beta(i,j)=(var_ml.alpha(i,j)^2/var_ml.alpha2(i,j)-1.0)/(var_ml.num(i,j)-1.0);
%         
%         var100_ml.jeff(i,j)=var100_ml.numer(i,j)/var100_ml.denom(i,j);
%         var100_ml.beta(i,j)=(var100_ml.alpha(i,j)^2/var100_ml.alpha2(i,j)-1.0)/(var100_ml.num(i,j)-1.0);
%     end
% end

% for i = 1:len_ml
%     temp(1:totfile)=var_ml.jeff(i,1:totfile);
%     var_ml.jeffmean(i)=mean(temp);
%     var_ml.jeffvar(i)=var(temp);
%     
%     temp1=0.0;
%     temp2=0.0;
%     for j=1:totfile
%         temp1=temp1+var_ml.jeff(i,j);
%         temp2=temp2+var_ml.jeff(i,j)^2;
%     end
%     var_ml.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
%     
%     
%     temp(1:totfile)=var100_ml.jeff(i,1:totfile);
%     var100_ml.jeffmean(i)=mean(temp);
%     var100_ml.jeffvar(i)=var(temp);
%     
%     temp1=0.0;
%     temp2=0.0;
%     for j=1:totfile
%         temp1=temp1+var100_ml.jeff(i,j);
%         temp2=temp2+var100_ml.jeff(i,j)^2;
%     end
%     var100_ml.jeffbeta(i)=(temp1^2/temp2-1.0)/(totfile-1.0);
% end






