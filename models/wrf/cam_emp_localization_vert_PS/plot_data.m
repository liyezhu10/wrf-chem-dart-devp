
x=[-40:dlen:40];

figure(1);
plot(x,tot_dist.jeff,'b-');
hold on;
plot(x,tot_dist.beta,'r-.');
legend('alpha','beta');
hold off;
xlabel('distance (unit: km)','FontSize',14);
title('U10__U, all grid points','FontSize',14);

% figure(2);
% plot(x,var_dist.jeffmean,'b-');
% hold on;
% plot(x,var_dist.jeffbeta,'r-.');
% legend('alpha','beta');
% hold off;
% xlabel('distance (unit: km)','FontSize',14);
% title('U10__U, all grid points','FontSize',14);

figure(3);
plot(tot_ml.jeff,'b-');
hold on;
plot(tot_ml.beta,'r-.');
legend('alpha','beta');
hold off;
xlabel('model levels','FontSize',14);
title('U10__U, all grid points','FontSize',14);

% figure(4);
% plot(var_ml.jeffmean,'b-');
% hold on;
% plot(var_ml.jeffbeta,'r-.');
% legend('alpha','beta');
% hold off;
% xlabel('model levels','FontSize',14);
% title('U10__U, all grid points','FontSize',14);



figure(5);
plot(x,tot100_dist.jeff,'b-');
hold on;
plot(x,tot100_dist.beta,'r-.');
legend('alpha','beta');
hold off;
xlabel('distance (unit: km)','FontSize',14);
title('U10__U, grid points <= 100m','FontSize',14);

% figure(6);
% plot(x,var100_dist.jeffmean,'b-');
% hold on;
% plot(x,var100_dist.jeffbeta,'r-.');
% legend('alpha','beta');
% hold off;
% xlabel('distance (unit: km)','FontSize',14);
% title('U10__U, grid points <= 100m','FontSize',14);

figure(7);
plot(tot100_ml.jeff,'b-');
hold on;
plot(tot100_ml.beta,'r-.');
legend('alpha','beta');
hold off;
xlabel('model levels','FontSize',14);
title('U10__U, grid points <= 100m','FontSize',14);

% figure(8);
% plot(var100_ml.jeffmean,'b-');
% hold on;
% plot(var100_ml.jeffbeta,'r-.');
% legend('alpha','beta');
% hold off;
% xlabel('model levels','FontSize',14);
% title('U10__U, grid points <= 100m','FontSize',14);


% it = 2;
% figure(1);
% plot(x,var_dist.jeff(:,it),'bo');
% hold on;
% plot(x,var_dist.beta(:,it),'r+');
% legend('alpha','beta');
% hold off;
% xlabel('distance (unit: km)','FontSize',14);
% title('U10__U, all grid points','FontSize',14);
% 
% figure(2);
% plot(var_ml.jeff(:,it),'bo');
% hold on;
% plot(var_ml.beta(:,it),'r+');
% legend('alpha','beta');
% hold off;
% xlabel('model levels','FontSize',14);
% title('U10__U, all grid points','FontSize',14);
% 
% 
% figure(3);
% plot(x,var100_dist.jeff(:,it),'bo');
% hold on;
% plot(x,var100_dist.beta(:,it),'r+');
% legend('alpha','beta');
% hold off;
% xlabel('distance (unit: km)','FontSize',14);
% title('U10__U, grid points <= 100m','FontSize',14);
% 
% figure(4);
% plot(var100_ml.jeff(:,it),'bo');
% hold on;
% plot(var100_ml.beta(:,it),'r+');
% legend('alpha','beta');
% hold off;
% xlabel('model levels','FontSize',14);
% title('U10__U, grid points <= 100m','FontSize',14);



