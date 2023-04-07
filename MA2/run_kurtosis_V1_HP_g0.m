
% run wsemiBSL TKDE for kurtosis transformed dataset

load('50_data.mat')
M = 50000;
cov_rw = [0.0161,0.0080;0.0080,0.0195];
start = [0.6,0.2];

eps = 0;
delta = 0.6;

KernelType = 'normal';

Transformation = 'HP';

n = 275; n_cov = 20000;
gamma = 0;

Whitening = 'Whitening';
WhiteningType = 3;
Version = 'V1';
double = 'double';

tic;
[theta, loglike] = semi_bsl_v2(y',M,n,cov_rw,start,eps,delta,KernelType,...
    Transformation,Version,gamma,Whitening,WhiteningType,n_cov,double);
time = toc;

filename = strrep(strjoin([M/1000,"_e",eps,"_d",delta,"_semiBSL_",Version,"_",Transformation,"_",Whitening,"_n",n,'.mat']),' ','')
save(filename,'theta','loglike','time')

