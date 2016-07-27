zzz_len=107781;

locnum=40;
dloc = 1000;    % unit: meter
locindex(locnum+1)=0;
for i=1:locnum
    locindex(i)=(-1)*(locnum-i+1)*dloc;
    locindex(locnum+1+i)=i*dloc;
end
% 
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

