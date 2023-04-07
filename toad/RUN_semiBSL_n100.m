% run wsemiBSL with moderate shrinkage

Transformation = 'KDE';
Whitening = 'Whitening';
double = 'single';
Version = 'V1';
WhiteningType = 3;
n_cov = 5000;

load('data_toads_model1.mat')
lag = [1, 2, 4, 8];
ssy = summStat_quantiles3(Y, lag)';
simArgs = struct('ntoads',ntoads,'ndays',ndays,'model',1,'d0',NaN);
sumArgs = struct('lag',lag);

M = 50000; 
start = [1.7, 35, 0.6];

KernelType = 'normal';

n = 100;
gamma = 0.3;

cov_rw = [0.0647,0.0079,-0.0008;0.0079,0.0030,0.0005;-0.0008,0.0005,0.0039];

[theta, loglike] = semi_bsl_v2(ssy',n,M,cov_rw,start,simArgs,sumArgs,KernelType,...
    Transformation, Version, gamma, Whitening, WhiteningType, n_cov, double);

save('50_wsemiBSL_n100_g0.3.mat','theta','loglike')
