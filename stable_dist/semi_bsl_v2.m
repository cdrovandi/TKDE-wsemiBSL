function [theta, loglike] = semi_bsl_v2(y,M,n,cov_rw,param_vec,KernelType,...
    Transformation,Version,gamma,Whitening,WhiteningType,n_cov)
%%
% Generic semiBSL implementation for alpha-stable state space model example
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
%%

loglike=zeros(M,1);
param_vec_curr = param_vec;
theta = zeros(M,2); 
ns = length(y); 

% transform the observed data
ssy = y;

if strcmp('Whitening',Whitening) == 1
    ssx = zeros(n_cov,ns);
    parfor i = 1:n_cov
        ssx(i,:) = sim_alphastable_statespace(param_vec_curr,ns);
    end
    W = whitening_ma_semi(WhiteningType,ssx);
else
    W = [];
end

% simulate data sets and compute summary statistics
ssx = zeros(n,ns);
parfor i = 1:n
    ssx(i,:) = sim_alphastable_statespace(param_vec_curr,ns);
end

if strcmp('V1',Version) == 1
    loglike_ind_curr = semipara_kernel_estimate_grc(ssy,ssx,KernelType,...
        Transformation,gamma,Whitening,W);
elseif strcmp('V2',Version) == 1
    loglike_ind_curr = semipara_kernel_estimate_grc_v2(ssy,ssx,KernelType,...
        Transformation,gamma,Whitening,W);
else
    fprintf('Must use version 1 or version 2\n')
end

param_vec_prop = param_vec_curr;
for i = 1:M
    i
    theta_tilde_curr = para_transformation(param_vec_curr(3:4),[0,-1],[2,1]);
    theta_tilde_prop = mvnrnd(theta_tilde_curr,cov_rw);
    param_vec_prop(3:4) = para_back_transformation(theta_tilde_prop,[0,-1],[2,1]);
     
    prob = jacobian_transformation(theta_tilde_prop,[0,-1],[2,1]) / jacobian_transformation(theta_tilde_curr,[0,-1],[2,1]);

    ssx = zeros(n,ns);
    parfor k = 1:n
        ssx(k,:) = sim_alphastable_statespace(param_vec_prop,ns);
    end

    if strcmp('V1',Version) == 1
        loglike_ind_prop = semipara_kernel_estimate_grc(ssy,ssx,KernelType,...
            Transformation,gamma,Whitening,W);
    elseif strcmp('V2',Version) == 1
        loglike_ind_prop = semipara_kernel_estimate_grc_v2(ssy,ssx,KernelType,...
            Transformation,gamma,Whitening,W);
    else
        fprintf('Must use version 1 or version 2\n')
    end
    
    
    if (prob*exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        param_vec_curr(3:4) = param_vec_prop(3:4);
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = param_vec_curr(3:4);
    loglike(i,1)=loglike_ind_prop;   
end

end


