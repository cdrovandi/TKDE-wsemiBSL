function [theta, loglike] = semi_bsl_v2(y,M,n,cov_rw,start,eps,delta,KernelType,...
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

loglike=zeros(M,1);
theta_curr = start;
theta = zeros(M,2); 
ns = length(y); 

% transform the observed data
ssy = sinh((1/delta)*(asinh(y)+eps));

if strcmp('Whitening',Whitening) == 1
    x = simulate_ma2(theta_curr,ns,n_cov)';
    ssx = sinh((1/delta)*(asinh(x)+eps)); % apply transformation to data  
    W = whitening_ma_semi(WhiteningType,ssx);
else
    W = [];
end

% simulate data sets and compute summary statistics
x = simulate_ma2(theta_curr,ns,n)';
ssx = sinh((1/delta)*(asinh(x)+eps)); % apply transformation to data  

if strcmp('V1',Version) == 1
    loglike_ind_curr = semipara_kernel_estimate_grc(ssy,ssx,KernelType,...
        Transformation,gamma,Whitening,W,double);
elseif strcmp('V2',Version) == 1
    loglike_ind_curr = semipara_kernel_estimate_grc_v2(ssy,ssx,KernelType,...
        Transformation,gamma,Whitening,W,double);
else
    fprintf('Must use version 1 or version 2\n')
end

for i = 1:M
    i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters

    if (sum(theta_prop) < -1 || diff(theta_prop) < -1 || abs(theta_prop(2)) > 1) 
        theta(i,:) = theta_curr;
        loglike(i) = loglike_ind_curr;
        continue;
    end
    
    x = simulate_ma2(theta_prop,ns,n)';
    ssx = sinh((1/delta)*(asinh(x)+eps));
    
    if strcmp('V1',Version) == 1
        loglike_ind_prop = semipara_kernel_estimate_grc(ssy,ssx,KernelType,...
            Transformation,gamma,Whitening,W,double);
    elseif strcmp('V2',Version) == 1
        loglike_ind_prop = semipara_kernel_estimate_grc_v2(ssy,ssx,KernelType,...
            Transformation,gamma,Whitening,W,double);
    else
        fprintf('Must use version 1 or version 2\n')
    end
    
    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i,1)=loglike_ind_prop;   
end

end
