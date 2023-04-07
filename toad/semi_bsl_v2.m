function [theta, loglike] = semi_bsl_v2(ssy,n,M,cov_rw,start,simArgs,sumArgs,KernelType,...
    Transformation,Version,gamma,Whitening,WhiteningType,n_cov,double)
%%
% Generic semiBSL implementation for toad example
%
% M --------------- number of iterations
% n --------------- number of simulated data sets
% cov_rw ---------- random walk covariance matrix
% param_vec ------- initial parameter value
% KernelType ------ kernel for KDE
% Transformation -- 'HP' or 'KDE'
% Version --------- 'V1' or 'V2'
% gamma ----------- Warton shrinkage parameter
% Whitening ------- 'Whitening' or 'NoWhitening'
% WhiteningType --- the type of whitening 3 is PCA whitening
% n_cov ----------- number of model simulations for estimating whitening matrix
% double ---------- whether an additional log transformation will be used and what kind of log transform
%%

ntoads = simArgs.ntoads;
ndays = simArgs.ndays;
model = simArgs.model;
d0 = simArgs.d0;
lag = sumArgs.lag;

loglike = zeros(M,1);
theta = zeros(M,3);

theta_curr = start; 
ns = length(ssy);

% simulate n_cov sets of data
ssx = zeros(n_cov,ns);
parfor k = 1:n_cov
    X = simulate_toads2(theta_curr,ntoads,ndays,model);
    ssx(k,:) = summStat_quantiles3(X,lag);
end

W = whitening(WhiteningType,ssx);

% simulate n sets of data
ssx = zeros(n,ns);
parfor k = 1:n
    X = simulate_toads2(theta_curr,ntoads,ndays,model);
    ssx(k,:) = summStat_quantiles3(X,lag);
end

% semi-parameteric estimate of the log-likelihood
loglike_ind_curr = semipara_kernel_estimate_grc(ssy,ssx,KernelType,Transformation,gamma,Whitening,W,double);

for i = 1:M
    i
    theta_tilde_curr = para_transformation(theta_curr);
    theta_tilde_prop = mvnrnd(theta_tilde_curr,cov_rw);
    theta_prop = para_back_transformation(theta_tilde_prop);
    prob = jacobian_transformation(theta_tilde_prop) / jacobian_transformation(theta_tilde_curr);

    ssx = zeros(n,ns);
    parfor k = 1:n
        X = simulate_toads2(theta_prop,ntoads,ndays,model);
        ssx(k,:) = summStat_quantiles3(X,lag);
    end
    
    loglike_ind_prop = semipara_kernel_estimate_grc(ssy,ssx,KernelType,Transformation,gamma,Whitening,W,double);
        
    if (prob * exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    
    theta(i,:) = theta_curr;
    loglike(i) = loglike_ind_curr;   
    
end

end