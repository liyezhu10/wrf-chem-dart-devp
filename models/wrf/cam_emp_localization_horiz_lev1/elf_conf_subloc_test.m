
% nfile = 104;
% for i = 1:nfile
%     string = sprintf('var_%d = load(''RADIOSONDE_TEMPERATURE_T_subloc%d'');',i,i);
%     eval(string);
%     
%     string = sprintf('tmpx_%d = var_%d(:,1);',i,i);
%     eval(string);
%     string = sprintf('tmpy_%d = var_%d(:,2);',i,i);
%     eval(string);
% end
% 
% %% regression toolbox
% for i = 1:nfile
%     string = sprintf('tmpx = tmpx_%d;',i);
%     eval(string);
%     string = sprintf('tmpy = tmpy_%d;',i);
%     eval(string);
% 
%     newy=[ones(size(tmpy)) tmpy];
%     string = sprintf('[b_%d,bint_%d,r_%d,rint_%d,stats_%d]=regress(tmpx,newy);',i,i,i,i,i);
%     eval(string);
%     clearvars tmpx tmpy newy;
% end

alpha(1:nfile) = 0.0;
alpha_conf(1:nfile) = 0.0;
for i = 1:nfile
    string = sprintf('alpha(i) = b_%d(2);',i);
    eval(string);
    string = sprintf('alpha_conf(i) = stats_%d(3);',i);
    eval(string);
end
    
% %% regression analysis
% for i = 1:nfile
%     tmpx=sprintf('tmpx_%d',i);
%     tmpy=sprintf('tmpy_%d',i);
%     
%     sprintf('p_%d',i)=polyfit(tmpy,tmpx,1);
%     
%     ngroup = 8;
%     tmp_nsample = size(tmpx,1);
%     tmp_randsample = floor(tmp_nsample/ngroup);
%     
%     rindex(1:tmp_randsample,1:ngroup) = 0.0;
%     tmpx_group(1:tmp_randsample,1:ngroup) = 0.0;
%     tmpy_group(1:tmp_randsample,1:ngroup) = 0.0;
%     alpha_group(1:ngroup) = 0.0;
%     ssxx(1:ngroup) = 0.0;
%     ssyy(1:ngroup) = 0.0;
%     ssxy(1:ngroup) = 0.0;
%     r2(1:ngroup) = 0.0;
%     se_b(1:ngroup) = 0.0;
%     
%     rindex_tot=randsample(tmp_nsample,tmp_nsample);
%     for j = 1:ngroup
%         for k = 1:tmp_randsample
%             tmpx_group(k,j) = tmpx(rindex(k,j));
%             tmpy_group(k,j) = tmpy(rindex(k,j));
%         end
%     end
%     
%     for j = 1:ngroup
%         p_tmp = polyfit(tmpy_group(:,j),tmpx_group(:,j),1);
%         
%     end
%     
%     
%     clearvars tmpx tmpy rindex tmpx_group tmpy_group alpha_group ssxx ssyy ssxy r2 se_b
% end
% 
% 
% 
% 
% for i = 1:ngroup
%     for j = 1:rand_sample
%         tmpx_group(j,i) = tmpx(rindex(j,i));
%         tmpy_group(j,i) = tmpy(rindex(j,i));
%     end
% end
% for i = 1:ngroup
%     p_tmp = polyfit(tmpy_group(:,i),tmpx_group(:,i),1);
%     alpha_group(i) = p_tmp(1);
%     
%     ssxx(i) = sum((tmpy_group(:,i)-mean(tmpy_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
%     ssyy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpx_group(:,i)-mean(tmpx_group(:,i))));
%     ssxy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
%     r2(i) = ssxy(i)^2/(ssxx(i)*ssyy(i));
%     se_b(i) = sqrt((ssyy(i)-ssxy(i)^2/ssxx(i))/(rand_sample-2))/sqrt(ssxx(i));
% end
% sum_alpha  = 0.0;
% sum_alpha2 = 0.0;
% for i = 1:ngroup
%     sum_alpha  = sum_alpha + alpha_group(i);
%     sum_alpha2 = sum_alpha2 + alpha_group(i)^2;
% end
% beta = (sum_alpha^2/sum_alpha2 - 1)/(ngroup-1);
% 
% 
% 
% 
% 
% 
% 
% % p=polyfit(tmpy,tmpx,1);
% % f=polyval(tmpy,p);
% % 
% % nsample = size(tmpx,1);
% % rcoef = p(1);
% % df = nsample-2;          % degree of freedom
% % var_tmpx=var(tmpx);
% % var_tmpy=var(tmpy);
% % se_b = sqrt(var_tmpx*(nsample-1)/(nsample-2))/sqrt(var_tmpy*(nsample-1));
% % tsig = rcoef/se_b;
% 
% 
% %% pseudo-group alpha
% nsample = size(tmpx,1);
% %rand_sample = 1000;
% %ngroup = floor(nsample/rand_sample);
% ngroup = 8;
% rand_sample = floor(nsample/ngroup);
% 
% rindex(1:rand_sample,1:ngroup) = 0;
% tmpx_group(1:rand_sample,1:ngroup) = 0.0;
% tmpy_group(1:rand_sample,1:ngroup) = 0.0;
% alpha_group(1:ngroup) = 0.0;
% ssxx(1:ngroup) = 0.0;
% ssyy(1:ngroup) = 0.0;
% ssxy(1:ngroup) = 0.0;
% r2(1:ngroup) = 0.0;
% se_b(1:ngroup) = 0.0;
% 
% %for i = 1:ngroup
% %    rindex(1:rand_sample,i) = randi(nsample,rand_sample,1);
% %end
% 
% rindex_tot=randsample(nsample,nsample);
% for i = 1:ngroup
%     rindex(1:rand_sample,i) = rindex_tot((i-1)*rand_sample+1:i*rand_sample);
% end
% 
% for i = 1:ngroup
%     for j = 1:rand_sample
%         tmpx_group(j,i) = tmpx(rindex(j,i));
%         tmpy_group(j,i) = tmpy(rindex(j,i));
%     end
% end
% for i = 1:ngroup
%     p_tmp = polyfit(tmpy_group(:,i),tmpx_group(:,i),1);
%     alpha_group(i) = p_tmp(1);
%     
%     ssxx(i) = sum((tmpy_group(:,i)-mean(tmpy_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
%     ssyy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpx_group(:,i)-mean(tmpx_group(:,i))));
%     ssxy(i) = sum((tmpx_group(:,i)-mean(tmpx_group(:,i))).*(tmpy_group(:,i)-mean(tmpy_group(:,i))));
%     r2(i) = ssxy(i)^2/(ssxx(i)*ssyy(i));
%     se_b(i) = sqrt((ssyy(i)-ssxy(i)^2/ssxx(i))/(rand_sample-2))/sqrt(ssxx(i));
% end
% sum_alpha  = 0.0;
% sum_alpha2 = 0.0;
% for i = 1:ngroup
%     sum_alpha  = sum_alpha + alpha_group(i);
%     sum_alpha2 = sum_alpha2 + alpha_group(i)^2;
% end
% beta = (sum_alpha^2/sum_alpha2 - 1)/(ngroup-1);
% 
% 
% 
% 
% 
% 
% 
% 
