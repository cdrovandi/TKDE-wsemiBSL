% run standard semiBSL for dataset 1
load('50_obs_data_stable.mat');
M = 50000; % number of MCMC iterations

cov_rw = [0.0846,0.0009;0.0009,0.2421];
param_vec = theta;
para_vec(3:4) = [1.2,0.5];

KernelType = 'normal';

Transformation = 'KDE';

n = 2000; % number of model simulations for estimating synthetic likelihood 
n_cov = 20000;
gamma = 1;

Whitening = 'NoWhitening';
WhiteningType = 3;
Version = 'V1';

tic;
[theta, loglike] = semi_bsl_v2(y,M,n,cov_rw,param_vec,KernelType,...
    Transformation,Version,gamma,Whitening,WhiteningType,n_cov);
time = toc;

save('50_KDE_stable.mat','theta','loglike','time')


