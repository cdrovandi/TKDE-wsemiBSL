function [f] = semipara_kernel_estimate_grc_v2(ssy,ssx,KernelType,Transformation,gamma,Whitening,W,double)


[n,ns] = size(ssx);
if strcmp('HP',Transformation) == 1
   
    for j = 1:ns
        [pdfy(j),~,~] = HP_trans_double_KDE(ssx(:,j),ssy(j),double,'cdf');
    end
    
elseif strcmp('KDE',Transformation) == 1
    
    for j = 1:ns
        [pdfy(j),~,bw] = ksdensity(ssx(:,j),ssy(j),'Kernel',KernelType);
        yu(j) = ksdensity(ssx(:,j),ssy(j),'function','cdf','Kernel',KernelType,'bandwidth',bw);
    end
    
else
    fprintf('Choose a transformation type KDE or HP\n')
end

if strcmp('Whitening',Whitening) == 1
    Rank = [];
    for i = 1:ns
        Rank(:,i) = tiedrank(ssx(:,i));
    end
    
    % compute quantiles
    q = norminv(Rank/(n+1));
    
    % transform quantiles
    q_t = q*W';
    
    Sigma_warton = warton(q_t, gamma);
    
    % rank of observed data in [ssx,ssy]
    Rank2 = [];
    for j = 1:ns
        Rank2(:,j) = tiedrank([ssx(:,j);ssy(j)]);
    end
    Rank_y = Rank2(n+1,:);
    Z_y = norminv(Rank_y/(n+2)); % quantile of observed data rank
    
    % transform quantiles of observed data
    Z_y_t = Z_y*W';
    
    L = chol(Sigma_warton);
    logdetA = 2*sum(log(diag(L)));
    
    %loglike_ind_curr
    f = -0.5*logdetA - 0.5*Z_y_t*inv(Sigma_warton)*Z_y_t' + sum(log(pdfy)) - sum(log(normpdf(Z_y, 0, 1)));
    
elseif strcmp('NoWhitening',Whitening) == 1 % no whitening
    % Gaussian rank correlation
    Rank = [];
    for i = 1:ns
        Rank(:,i) = tiedrank(ssx(:,i));
    end
    RHOHAT =  eye(ns);
    den = sum(norminv((1:n)/(n+1)).^2);
    for i = 2:ns
        for j = 1:(i-1)
            RHOHAT(i,j) = sum(norminv(Rank(:,i)/(n+1)) .* norminv(Rank(:,j)/(n+1))) / den;
            RHOHAT(j,i) = RHOHAT(i,j);
        end
    end
    
    RHOHAT_warton = gamma * RHOHAT + (1 - gamma) * eye(ns);
    
    % rank of observed data in [ssx,ssy]
    Rank2 = [];
    for j = 1:ns
        Rank2(:,j) = tiedrank([ssx(:,j);ssy(j)]);
    end
    Rank_y = Rank2(n+1,:);
    Z_y = norminv(Rank_y/(n+2)); % quantile of observed data rank
    
    f = -0.5*log(det(RHOHAT_warton)) - 0.5*Z_y*(inv(RHOHAT_warton)-eye(ns))*Z_y' + sum(log(pdfy));


else
    fprintf('Please specify a valid whitening type\n')
end

end









