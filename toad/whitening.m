function [W] = whitening(type,ssx)
% Estimate whitening matrix based on simulated summaries ssx
% type 1 - cholesky whitening
% type 2 - ZCA whitenig
% type 3 - PCA whitening
% type 4 - ZCA-corr whitening
% type 5 - PCA-corr whitening

[n,ns] = size(ssx);

% Compute ranks of simulated data sets
Rank = [];
for i = 1:ns
    Rank(:,i) = tiedrank(ssx(:,i));
end

q = [];
q = norminv(Rank/(n+1));

% compute covariance
the_cov = cov(q);

if type == 1 % cholesky
    W = chol(inv(the_cov)); % calculate the whitening matrix
elseif type == 2 % ZCA
    W = the_cov^(-0.5);
elseif type == 3 % PCA
    [U,D] = eig(the_cov);
    W = D^(-0.5)*U';
elseif type == 4 % ZCA-corr
    V = diag(diag(the_cov));
    P =  V^(-0.5) * the_cov * V^(-0.5);
    W = P^(-0.5) * V^(-0.5);
else % type == 5 % PCA-corr
    V = diag(diag(the_cov));
    P =  V^(-0.5) * the_cov * V^(-0.5);
    [G,Theta] = eig(P);
    W = Theta^(-0.5)*G'*V^(-0.5);
end


end