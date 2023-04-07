function f = semipara_kernel_estimate_grc(ssy,ssx,KernelType,Transformation,gamma,Whitening,W,double)
%%
% generic function for estimating synthetic likelihood semi-parametrically
%
% ssy - observed summary
% ssx - matrix of simulated summaries
% KernelType - which kernel type to use in KDE
% Transformation - indicates if we are using KDE or TKDE
% gamma - amount of shrinkage in Warton estimator if using whitening transformation
% Whitening - indicates if we are using whitening transform
% W - whitening matrix if we are using whitening
% double - indicates if we are using additional log transform with TKDE and what log transform to use
%
%%
if nargin < 3
    KernelType = 'epanechnikov';
end

[n,ns] = size(ssx);

if strcmp('KDE',Transformation)==1
    for j = 1:ns
        [pdfy(j),~,bw] = ksdensity(ssx(:,j),ssy(j),'Kernel',KernelType);
        yu(j) = ksdensity(ssx(:,j),ssy(j),'function','cdf','Kernel',KernelType,'bandwidth',bw);
    end
elseif strcmp('HP',Transformation)==1
    parfor j = 1:ns
        [pdfy(j),~,yu(j)] = HP_trans_double_KDE(ssx(:,j),ssy(j),double,'cdf');
%         [pdfy(j),~,~,yu(j)] = HP_trans_KDE2_CDF(ssx(:,j),ssy(j));
    end
else
    fprintf('Choose a transformation type KDE or HP\n')
end

% Gaussian rank correlation
Rank = [];
for i = 1:ns
    Rank(:,i) = tiedrank(ssx(:,i));
end


if strcmp('Whitening',Whitening) == 1
    
    % compute quantiles
    q = norminv(Rank/(n+1));
    
    % transform quantiles
    q_t = q*W';
    
    [Sigma_warton] = warton(q_t, gamma);
    
    L = chol(Sigma_warton);
    logdetA = 2*sum(log(diag(L)));
    
    Z_y = norminv(yu);
    Z_y_t = Z_y*W';
    
    %loglike_ind_curr
    f = -0.5*logdetA - 0.5*Z_y_t*inv(Sigma_warton)*Z_y_t' + sum(log(pdfy)) - sum(log(normpdf(Z_y, 0, 1)));
    
    if isnan(f) == 1
        f = -Inf;
    end
    
    
elseif strcmp('NoWhitening',Whitening) == 1
    
    RHOHAT=  eye(ns);
    den = sum(norminv((1:n)/(n+1)).^2);
    for i = 2:ns
        for j = 1:(i-1)
            RHOHAT(i,j) = sum(norminv(Rank(:,i)/(n+1)) .* norminv(Rank(:,j)/(n+1))) / den;
            RHOHAT(j,i) = RHOHAT(i,j);
        end
    end
    
    f = log(copulapdf('Gaussian',yu,RHOHAT)) + sum(log(pdfy));
    
else
    fprintf('must specify a valid whitening or not\n');
end

end

