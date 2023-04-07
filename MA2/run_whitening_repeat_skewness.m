
% perform repeated whitening runs for sknewess dataset

load('50_data.mat')
M = 50000;
cov_rw = [0.0161,0.0080;0.0080,0.0195];
start = [0.6,0.2];

eps = -1;
delta = 1;

KernelType = 'normal';

Transformation = 'KDE';

n = 90; n_cov = 20000;
gamma = 0;

Whitening = 'Whitening';
WhiteningType = 3;
Version = 'V1';
double = 'single';

load('theta_test_whitening_repeat.mat')

theta_all = cell(24,1);

parfor i = 1:24

    [theta, loglike] = semi_bsl_v2(y',M,n,cov_rw,theta_test(i,:),eps,delta,KernelType,...
        Transformation,Version,gamma,Whitening,WhiteningType,n_cov,double);

    theta_all{i} = theta;
end


save('results_whitening_repeat_skewness.mat','theta_all');


%% produce results

load('results_whitening_repeat_skewness.mat')


load('theta_test_whitening_repeat.mat')

load('exact_results.mat')
theta0 = [0.6 0.2];
cov_exact = cov(theta);

len = 0.01;
wid = 0.01;
theta1_grid = 0:len:1;
theta2_grid = -0.2:wid:1;
[X1,X2] = meshgrid(theta1_grid, theta2_grid);
X1 = X1(:);
X2 = X2(:);

[f_true,~] = ksdensity(theta,[X1 X2]);


for i = 1:24
    theta_rep = theta_all{i};
    [f_bsl,~] = ksdensity(theta_rep,[X1 X2]);
    tv_bsl = sum(abs(f_true - f_bsl));
    tv_distance(i)  = 0.5*tv_bsl*len*wid;
    theta_dist(i) = sqrt((theta_test(i,:) - theta0)*inv(cov_exact)*(theta_test(i,:) - theta0)');
end


save('tv_and_dists_whitening_repeat_skewness.mat','tv_distance','theta_dist');





