% run standard semiBSL for dataset 2
load('50_obs_data_stable_alpha0.7.mat');
M = 50000; % number of MCMC iterations

cov_rw = [0.0596,-0.0342;-0.0342,0.3760];

param_vec = theta;
% para_vec(3:4) = [1.2,0.5];

KernelType = 'normal';

Transformation = 'KDE';

n = 20000; % number of model simulations for estimating synthetic likelihood 
n_cov = 20000;
gamma = 1;

Whitening = 'NoWhitening';
WhiteningType = 3;
Version = 'V1';

tic;
[theta, loglike] = semi_bsl_v2(y,M,n,cov_rw,param_vec,KernelType,...
    Transformation,Version,gamma,Whitening,WhiteningType,n_cov);
time = toc;

save('50_KDE_stable_alpha0.7.mat','theta','loglike','time')


