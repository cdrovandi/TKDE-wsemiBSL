function [theta, loglike] = bayes_ma(ssy,M,cov_rw,start)

%%
% exact MCMC for MA2 model
%
% ssy - observed time series
% M - number of MCMC iterations is the total number of steps attempted (iterations)
% cov_rw - the covariance matrix used in the random walk (use multivariate random normal for proposals)
% start - starting value for chain
%%

loglike=zeros(M,1);

theta_curr = start;
theta = zeros(M,2);

ns = length(ssy);

the_cov = zeros(ns,ns);
for k = 1:ns
    the_cov(k,k) = 1 + theta_curr(1)^2 + theta_curr(2)^2;
    if (k < ns)
        the_cov(k,k+1) = theta_curr(1) + theta_curr(1)*theta_curr(2);
        the_cov(k+1,k) = the_cov(k,k+1);
    end
    if (k < ns-1)
        the_cov(k,k+2) = theta_curr(2);
        the_cov(k+2,k) = the_cov(k,k+2);
    end
end

loglike_ind_curr = -0.5*log(det(the_cov)) - 0.5*ssy*inv(the_cov)*ssy';

for i = 1:M
    i
    theta_prop = mvnrnd(theta_curr,cov_rw); %Proposing a new pair of parameters
    if (sum(theta_prop) < -1 || diff(theta_prop) < -1 || abs(theta_prop(2)) > 1)
        theta(i,:) = theta_curr;
        loglike(i)=loglike_ind_curr;
        continue;
    end


    the_cov = zeros(ns,ns);
    for k = 1:ns
        the_cov(k,k) = 1 + theta_prop(1)^2 + theta_prop(2)^2;
        if (k < ns)
            the_cov(k,k+1) = theta_prop(1) + theta_prop(1)*theta_prop(2);
            the_cov(k+1,k) = the_cov(k,k+1);
        end
        if (k < ns-1)
            the_cov(k,k+2) = theta_prop(2);
            the_cov(k+2,k) = the_cov(k,k+2);
        end
    end

    loglike_ind_prop = -0.5*log(det(the_cov)) - 0.5*ssy*inv(the_cov)*ssy';

    if (exp(loglike_ind_prop - loglike_ind_curr) > rand)
        %fprintf('*** accept ***\n');
        theta_curr = theta_prop;
        loglike_ind_curr = loglike_ind_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i)=loglike_ind_curr;
end
end

