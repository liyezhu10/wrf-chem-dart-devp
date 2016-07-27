
% follow process_zzz_modellevel.m

locdist.numer(1:2*locnum+1)=0.0;
locdist.denom(1:2*locnum+1)=0.0;
locdist.alpha(1:2*locnum+1)=0.0;
locdist.alpha2(1:2*locnum+1)=0.0;
locdist.samplesize(1:2*locnum+1)=0.0;

obsnumtemp=0;

for i = 1:obsnum
    if (wrf.obslev(i) <= 49)
        obsnumtemp=obsnumtemp+1;
        for ind = 1:2*locnum+1
            if ( var.locind(ind) > 0 )
                for k = 1:var.locind(ind)
                    locdist.samplesize(ind)=locdist.samplesize(ind)+1;
                    locdist.numer(ind)=locdist.numer(ind)+var.numer(ind,k,i);
                    locdist.denom(ind)=locdist.denom(ind)+var.denom(ind,k,i);
                    locdist.alpha(ind)=locdist.alpha(ind)+var.alpha(ind,k,i);
                    locdist.alpha2(ind)=locdist.alpha2(ind)+var.alpha2(ind,k,i);
                end
            end
        end   
    end
end


for i = 1:2*locnum+1
    locdist.jeff(i)=locdist.numer(i)/locdist.denom(i);
    locdist.beta(i)=(locdist.alpha(i)^2/locdist.alpha2(i)-1)/(locdist.samplesize(i)-1);
end



