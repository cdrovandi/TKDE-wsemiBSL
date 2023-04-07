% running standard semiBSL mg1 example
lower = [0,0,0];

upper = [10,20,0.5];

start = [1,5,0.2];

n = 500; % number of model simulations for estimating synthetic likelihood

M = 100000; % number of MCMC iterations

type = 'KDE'; % using standard KDE not TKDE

load('mg1_dat.mat')

cov_rw = [0.0323,0.0004,-0.0055;...
    0.0004,0.0085,-0.0015;...
    -0.0055,-0.0015,0.0436];
double = 'single'; % no additional transforms

[theta,loglike] = semi_bsl(ssy,n,M,cov_rw,start,lower,upper,type,double);

save('100_mg1_KDE.mat','theta','loglike')
