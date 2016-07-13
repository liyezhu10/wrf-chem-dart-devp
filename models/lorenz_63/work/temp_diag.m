n_obs=200;
for i=1:n_obs
    temp = ReadASCIIObsSeq(['./advance_time_' num2str(i,'%.4d') '/obs_seq.final']);
    true_orbit(i,:)=temp.obs(2,:);
    obs_orbit(i,:)=temp.obs(1,:);
    ens_mean(i,:)=temp.obs(4,:);
    
end

plot(obs_orbit(:,1),'g')
hold on;
plot(true_orbit(:,1),'r')
plot(ens_mean(:,1),'b')

temp_all=ReadASCIIObsSeq('obs_seq.final');

for i=1:200
    enkf(i,1)=temp_all.obs(4,1+(i-1)*3)
    enkf(i,2)=temp_all.obs(4,2+(i-1)*3)
    enkf(i,3)=temp_all.obs(4,3+(i-1)*3)
end