function [theta, loglike] = semi_bsl(ssy,n,M,cov_rw,start,lower,upper,type,double)
%%
% generic semiBSL function for mg1 example
%
% ssy - observed summary
% n - number of model simulations for estimating synthetic likelihood
% M - number of MCMC iterations
% cov_rw - random walk covariance matrix (on transformed space)
% lower - lower limits of prior
% upper - upper limits of prior
% type - whether or not using KDE or TKDE
% double - whether an additional log transform is used and what type
%%

T = 51;
KernelType = 'normal';

loglike = zeros(M,1);
theta = zeros(M,3);

theta_curr = start; 
ns = length(ssy);

% simulate n sets of data
ssx_curr = zeros(n,ns);
parfor k = 1:n
    ssx_curr(k,:) = simulate_mg1(theta_curr,T);
end


% semi-parameteric estimate of the log-likelihood
loglike_ind_curr = semipara_kernel_estimate_grc(ssy,ssx_curr,KernelType,type,1,'NoWhitening',1,double);

i = 1;
while i <= M
    i
       
    % propose new parameter value via normal RW        
    theta_tilde_curr = para_transformation(theta_curr,lower,upper);
    theta_tilde_prop = mvnrnd(theta_tilde_curr,cov_rw);
    theta_prop = para_back_transformation(theta_tilde_prop,lower,upper);
    prob = jacobian_transformation(theta_tilde_prop,lower,upper) /...
        jacobian_transformation(theta_tilde_curr,lower,upper);
    
    if (theta_prop(2) - theta_prop(1) > 10 || theta_prop(2) - theta_prop(1) < 0)
        theta(i,:) = theta_curr; 
        loglike(i) = loglike_ind_curr;
        continue;
    end
      

    ssx_prop = zeros(n,ns);
    parfor k = 1:n
        ssx_prop(k,:) = simulate_mg1(theta_prop,T);
    end
    
    loglike_ind_prop = semipara_kernel_estimate_grc(ssy,ssx_prop,KernelType,type,1,'NoWhitening',1,double);

	
	% accept or reject proposal with probability given by the metropolis-hastings ratio	      
    if (prob * exp(loglike_ind_prop - loglike_ind_curr) > rand)
        fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    
    theta(i,:) = theta_curr;
    loglike(i) = loglike_ind_prop;   
    
    i = i + 1;
        
    
end

end


