
% follow process_zzz_modellevel.m

loc.numer(1:50)=0.0;
loc.denom(1:50)=0.0;
loc.alpha(1:50)=0.0;
loc.alpha2(1:50)=0.0;

obsnumtemp=0;

for i = 1:obsnum
    if (wrf.obslev(i) <= 49)
        obsnumtemp=obsnumtemp+1;
    for j = 1:wrf.obslev(i)
        loc.numer(j)=loc.numer(j)+wrf.numer(j,i);
        loc.denom(j)=loc.denom(j)+wrf.denom(j,i);
        loc.alpha(j)=loc.alpha(j)+wrf.alpha(j,i);
        loc.alpha2(j)=loc.alpha2(j)+wrf.alpha2(j,i);
    end
    end
end

for i = 1:49
    loc.jeff(i)=loc.numer(i)/loc.denom(i);
    loc.beta(i)=(loc.alpha(i)^2/loc.alpha2(i)-1)/(obsnumtemp-1);
end














