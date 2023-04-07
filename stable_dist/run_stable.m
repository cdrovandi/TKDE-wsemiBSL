mu = 5;
phi = 1;
alpha = 2; % 0 < alpha <= 2
beta = 0.5; % -1 <= beta <= 1
gamma = 1;
delta = 0;
sigma = 0.2;

theta = [mu,phi,alpha,beta,gamma,delta,sigma];

% theta = [5 1 1.95 -0.9 1 0 0.2];
% 
% pd2 = makedist('Stable','alpha',alpha,'beta',beta,'gam',gamma,'delta',delta);
% x = linspace(-200,200,10000);
% figure, plot(x,pdf(pd2,x),'b') % plot the exact pdf

n = 1000;

for i = 1:n
    ssx(i,:) = sim_alphastable_statespace(theta,50);
end
figure
imagesc(abs(corr(ssx)))

figure, hold on

for i = 1:50
    [f,xi,~] = HP_trans_double_KDE(ssx(:,i),linspace(min(ssx(:,i)),max(ssx(:,i)),10000),'double','1');
    plot(xi,f)
end

figure
hold on
for i = 1:50
    ksdensity(ssx(:,i),linspace(min(ssx(:,i)),max(ssx(:,i)),10000));
end
% % 
y = sim_alphastable_statespace(theta,50);
save('50_obs_data_stable_alpha0.7.mat','y','theta')