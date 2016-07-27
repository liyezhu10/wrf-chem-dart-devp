zzz_len=107781;

locnum=40;
dloc = 1000;    % unit: meter
locindex(locnum+1)=0;
for i=1:locnum
    locindex(i)=(-1)*(locnum-i+1)*dloc;
    locindex(locnum+1+i)=i*dloc;
end

var.locind(1:2*locnum+1,1:zzz_len)=0;
wrf.obslev(1:zzz_len)=0;
obsnum=1;
var.obsi(obsnum)=zzz(1,1);
var.obsj(obsnum)=zzz(1,2);
var.locind(zzz(1,3),obsnum)=1;
% var.xtruth(zzz(1,3),1,obsnum)=zzz(1,4);
% var.xfmean(zzz(1,3),1,obsnum)=zzz(1,5);
% var.obsinc(zzz(1,3),1,obsnum)=zzz(1,6);
wrf.obslev(obsnum)=1;
wrf.xtruth(1,obsnum)=zzz(1,4);
wrf.xfmean(1,obsnum)=zzz(1,5);
wrf.obsinc(1,obsnum)=zzz(1,6);
wrf.numer(1,obsnum)=zzz(1,7);
wrf.denom(1,obsnum)=zzz(1,8);
wrf.alpha(1,obsnum)=zzz(1,9);
wrf.alpha2(1,obsnum)=zzz(1,10);
var.numer(zzz(1,3),1,obsnum)=zzz(1,7);
var.denom(zzz(1,3),1,obsnum)=zzz(1,8);
var.alpha(zzz(1,3),1,obsnum)=zzz(1,9);
var.alpha2(zzz(1,3),1,obsnum)=zzz(1,10);

for i = 2:zzz_len
%     istatus=0;
%     for j=1:obsnum
%         if ( (var.obsi(j) == zzz(i,1)) && (var.obsj(j) == zzz(i,2)) )
%             istatus=1;
%         end
%     end
%     if ( istatus == 0 ) 
%         obsnum=obsnum+1;
%         var.obsi(obsnum)=zzz(i,1);
%         var.obsj(obsnum)=zzz(i,2);
%     end

    istatus = 0;
    if ( (var.obsi(obsnum) == zzz(i,1)) && (var.obsj(obsnum) == zzz(i,2)))
        istatus = 1;
    end
    if ( istatus == 0 )
        obsnum=obsnum+1;
        var.obsi(obsnum)=zzz(i,1);
        var.obsj(obsnum)=zzz(i,2);
    end

    var.locind(zzz(i,3),obsnum)=var.locind(zzz(i,3),obsnum)+1;
%     var.xtruth(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,4);
%     var.xfmean(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,5);
%     var.obsinc(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,6); 
    wrf.obslev(obsnum)=wrf.obslev(obsnum)+1;
    wrf.xtruth(wrf.obslev(obsnum),obsnum)=zzz(i,4);
    wrf.xfmean(wrf.obslev(obsnum),obsnum)=zzz(i,5);
    wrf.obsinc(wrf.obslev(obsnum),obsnum)=zzz(i,6);
    wrf.numer(wrf.obslev(obsnum),obsnum)=zzz(i,7);
    wrf.denom(wrf.obslev(obsnum),obsnum)=zzz(i,8);
    wrf.alpha(wrf.obslev(obsnum),obsnum)=zzz(i,9);
    wrf.alpha2(wrf.obslev(obsnum),obsnum)=zzz(i,10);
    var.numer(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,7);
    var.denom(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,8);
    var.alpha(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,9);
    var.alpha2(zzz(i,3),var.locind(zzz(i,3),obsnum),obsnum)=zzz(i,10);    
end





% psfname='tutu.ps';

% for iobs = 1:1
%     for ind = 1:2*locnum+1
%         if ( var.locind(ind,iobs) > 0 )
%             for isub = 1:var.locind(ind,iobs)
%                 plot(ind, var.alpha(ind,isub,iobs),'bo');
%                 hold on;
%             end
%         end
%     end
%     
%     string=sprintf('title(''obsi %d, obsj %d'',''FontSize'',14)',var.obsi(iobs),var.obsj(iobs));
%     eval(string);
% 
%     xlim([40,81]);
%     
% %     print('-dpsc',psfname,'-append')
% %     clf;
% end


% temp(1:100)=0.0;
% for iobs = 1:1
%     
% %     for ind = 1:wrf.obslev(iobs)
% %         temp(ind)=(wrf.xtruth(ind,iobs)-wrf.xfmean(ind,iobs))*wrf.obsinc(ind,iobs)/wrf.denominator(ind,iobs);
% %     end
% %     for ind = 1:2*locnum+1
% %         if ( var.locind(ind,iobs) > 0 )
% %             for isub = 1:var.locind(ind,iobs)
% %                 plot(ind, var.xtruth(ind,isub,iobs),'bo');
% %                 hold on;
% %                 plot(ind, var.xfmean(ind,isub,iobs),'b+');
% %                 hold on;
% %             end
% %         end
% %     end
% clf
% %     plot(wrf.xtruth(1:wrf.obslev(iobs),iobs),'bo');
% %     hold on;
% %     plot(wrf.xfmean(1:wrf.obslev(iobs),iobs),'b+');
% %     hold on;
% 
%     plot(wrf.xtruth(1:wrf.obslev(iobs))-wrf.xfmean(1:wrf.obslev(iobs)), 'bo');
%     hold on;
%     plot(wrf.obsinc(1:wrf.obslev(iobs),iobs),'r+');
%     hold on;
% %     plot(wrf.alpha(1:wrf.obslev(iobs),iobs),'k*');
% %     plot(temp(40:var.obslev(iobs)),'r+');
%     string=sprintf('title(''obsi %d, obsj %d'',''FontSize'',14)',var.obsi(iobs),var.obsj(iobs));
%     eval(string);
% 
% %    xlim([40,81]);
%     
% %     print('-dpsc',psfname,'-append')
% %     clf;
% end










% var.loc(1:2*locnum+1)=0;
% for i = 1:zzz_len
%     var.loc(zzz(i,1))=var.loc(zzz(i,1))+1;
%     var.numer(zzz(i,1),var.loc(zzz(i,1)))=zzz(i,2);
%     var.denom(zzz(i,1),var.loc(zzz(i,1)))=zzz(i,3);
%     var.alpha(zzz(i,1),var.loc(zzz(i,1)))=zzz(i,4);
%     var.alpha2(zzz(i,1),var.loc(zzz(i,1)))=zzz(i,4)^2;
% end
% 
% var.loc(locnum+1)=var.loc(locnum)+var.loc(locnum+1)+var.loc(locnum+2);
% var.numer(locnum+1)=var.numer(locnum)+var.numer(locnum+1)+var.numer(locnum+2);
% var.denom(locnum+1)=var.denom(locnum)+var.denom(locnum+1)+var.denom(locnum+2);
% var.alpha(locnum+1)=var.alpha(locnum)+var.alpha(locnum+1)+var.alpha(locnum+2);
% var.alpha2(locnum+1)=var.alpha2(locnum)+var.alpha2(locnum+1)+var.alpha2(locnum+2);
% 
% for i=1:2*locnum+1
%     dimtemp=var.loc(i);
%     var.numermean(i)=mean(var.numer(i,1:dimtemp));
%     var.denommean(i)=mean(var.denom(i,1:dimtemp));
%     var.jeff(i)=var.numermean(i)/var.denommean(i);
%     var.alphasum(i)=sum(var.alpha(i,1:dimtemp));
%     var.alpha2sum(i)=sum(var.alpha2(i,1:dimtemp));
%     var.conf(i)=(var.alphasum(i)^2/var.alpha2sum(i)-1)/(dimtemp-1);
% end


% tempind=205;
% templen=var.loc(tempind);
% tempnumer(1:templen)=var.numer(tempind,1:templen);
% tempdenom(1:templen)=var.denom(tempind,1:templen);
% tempalpha(1:templen)=var.alpha(tempind,1:templen);
% tempalphasq(1:templen)=var.alpha2(tempind,1:templen);
% tempalpha_mean=mean(tempalpha);
% tempalpha_var=(std(tempalpha))^2;
% tempalpha_conf=((sum(tempalpha))^2/(sum(tempalphasq))-1)/(templen-1);
% % for i=1:templen
% 
% clearvars tempalpha2;
% j=0;
% for i=1:templen
%     if ( abs(tempalpha(i)) < 10.0 )
%         j=j+1;
%         tempalpha2(j)=tempalpha(i);
%         tempalpha2sq(j)=tempalphasq(i);
%         tempnumer2(j)=tempnumer(i);
%         tempdenom2(j)=tempdenom(i);
%     end
% end
% tempalpha2_mean=mean(tempalpha2);
% tempalpha2_var=(std(tempalpha2))^2;
% tempalpha2_conf=((sum(tempalpha2))^2/(sum(tempalpha2sq))-1)/(j-1);


% alphasum=sum(var.alpha(202,:));
% alpha2sum=sum(var.alpha2(202,:));

