
% run exact Bayes for MA(2) example

load('50_data.mat')
M = 100000;
cov_rw = [0.0161,0.0080;0.0080,0.0195];
start = [0.6,0.2];


[theta, loglike] = bayes_ma(y',M,cov_rw,start);


save('exact_results.mat','theta','loglike');




