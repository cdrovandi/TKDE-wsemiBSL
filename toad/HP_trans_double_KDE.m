function [pdf,xi,yu] = HP_trans_double_KDE(x,x_eval,double,cdf)

%%
% estimation of marginal density via TKDE with possible additional log
% transform.
%
% x - samples for estimating TKDE
% x_eval - where we want to evaluate estimated density
% double - indicating whether we want additional log transform and what type
% cdf - indicating if we want to also estimate the cdf 
%%


if nargin == 1
    x_eval = linspace(min(x),max(x),100);
end
xi = x_eval;

% transform the observed dataset
x_cent = x - median(x);
pos_inds = find(x_cent>=0);
neg_inds = find(x_cent<0);

if strcmp('double',double)
    x_cent(pos_inds) = log(1+x_cent(pos_inds));
    x_cent(neg_inds) = -log(1-x_cent(neg_inds));
elseif strcmp('pos_skew',double)
    min_xcent = min(x_cent);
    x_eval_min = min(x_eval);
%     x_eval_max = max(x_eval);
    [Y,~] = sort([x;x_eval_min]);
    if Y(1) == x_eval_min
        shift = abs(Y(2)-Y(1))+1;
    else
        shift = 0;
    end
%     shift = abs(max(x)-quantile(x,0.95));
    x_cent = log(1+x_cent-min_xcent+shift);
    med_xcent = median(x_cent);
    x_cent = x_cent - med_xcent;
    pos_inds = find(x_cent>=0);
    neg_inds = find(x_cent<0);
elseif strcmp('neg_skew',double)
    max_xcent = max(x_cent);
    x_eval_max = max(x_eval);
    [Y,~] = sort([x;x_eval_max]);
    if Y(end)==x_eval_max
        shift = abs(Y(end)-Y(end-1))+1;
    else
        shift = 0;
    end
    %     shift = abs(min(x)-quantile(x,0.05));
%     x_cent = -log(1-x_cent+max_xcent+std(x_cent));
    x_cent = -log(1-x_cent+max_xcent+shift);
    med_xcent = median(x_cent);
    x_cent = x_cent - med_xcent;
    pos_inds = find(x_cent>=0);
    neg_inds = find(x_cent<0);
end


% define the hyperbolic power transformation
f = @(x,Alpha,Beta,Lambda) Alpha*sinh(Beta*x).*(sech(Beta*x).^Lambda)/Beta;

%% initialise parameter values
q = 0.95;
x_q = quantile(x_cent,q);
x_p = x_q/2;
p = length(x_cent(x_cent<x_p))/length(x_cent);

z_p = norminv(p,0,1);
z_q = norminv(q,0,1);

if z_q/z_p > x_q/x_p
    % low kurtosis
    Beta_p = acosh(z_q/(2*z_p))/abs(x_p);
%     fprintf('** low kurtosis **\n')
else
    % high kurtosis
    Beta_p = acosh(z_p/(z_q-z_p))/abs(x_p);
%     fprintf('** high kurtosis **\n')
end

% Lambda_p -- if high kurtosis, we expect lambda close to 1, if low
% kurtosis we expect lambda close to 0
Lambda_p = (log(z_p/z_q) + log(sinh(Beta_p*x_q)) - log(sinh(Beta_p*x_p)))/...
    (log(sech(Beta_p*x_p)) - log(sech(Beta_p*x_q)));

q = 0.05;
x_q = quantile(x_cent,q);
x_p = x_q/2;
p = length(x_cent(x_cent<x_p))/length(x_cent);

z_p = norminv(p,0,1); 
z_q = norminv(q,0,1);

if z_q/z_p > x_q/x_p
    % low kurtosis
    Beta_n = acosh(z_q/(2*z_p))/abs(x_p);
%     fprintf('** low kurtosis **\n')
else
    % high kurtosis
    Beta_n = acosh(z_p/(z_q-z_p))/abs(x_p);
%     fprintf('** high kurtosis **\n')
end

% Lambda_n -- if high kurtosis, we expect lambda close to 1, if low
% kurtosis we expect lambda close to 0
Lambda_n = (log(z_p/z_q) + log(sinh(Beta_n*x_q)) - log(sinh(Beta_n*x_p)))/...
    (log(sech(Beta_n*x_p)) - log(sech(Beta_n*x_q)));

x_pos = x_cent(pos_inds);
x_neg = x_cent(neg_inds);

%% perform optimisation
fun_p = @(BL_t) evaluate_loglik(x_pos,BL_t);

% introduce unbounded parameter values for optimisation
Beta_p_t = log(Beta_p); % since Beta > 0
Lambda_p_t = log((Lambda_p+1)/(1-Lambda_p)); % since lambda <= 1

% perform optimisation
opts = optimset('Display','off');
MLE_p = fminsearch(fun_p,[Beta_p_t;Lambda_p_t],opts);

% opts2 = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton');
% MLE_p = fminunc(fun_p,[Beta_p_t;Lambda_p_t],opts2);

% inverse transformation to obtain MLEs
Beta_p_MLE = exp(MLE_p(1));
if MLE_p(2) == Inf
    Lambda_p_MLE = 1;
else
    Lambda_p_MLE = (exp(MLE_p(2))-1)/(1+exp(MLE_p(2)));
end

fun_n = @(BL_t) evaluate_loglik(x_neg,BL_t);

% introduce unbounded parameter values for optimisation
Beta_n_t = log(Beta_n);
Lambda_n_t = log((Lambda_n+1)/(1-Lambda_n));

% perform optimisation
MLE_n = fminsearch(fun_n,[Beta_n_t;Lambda_n_t],opts);
% MLE_n = fminunc(fun_n,[Beta_n_t;Lambda_n_t],opts2);

% inverse transformation to obtain MLEs
Beta_n_MLE = exp(MLE_n(1));
if MLE_n(2) == Inf
    Lambda_n_MLE = 1;
else
    Lambda_n_MLE = (exp(MLE_n(2))-1)/(1+exp(MLE_n(2)));
end


%% transform data
Alpha_MLE = ( (1/length(x_cent)) * ( sum((((sinh(Beta_p_MLE*x_pos).*sech(Beta_p_MLE*x_pos).^Lambda_p_MLE))/Beta_p_MLE).^2) +...
    sum((((sinh(Beta_n_MLE*x_neg).*sech(Beta_n_MLE*x_neg).^Lambda_n_MLE))/Beta_n_MLE).^2) ) )^(-0.5);

pos_trans = f(x_pos, Alpha_MLE, Beta_p_MLE, Lambda_p_MLE);
neg_trans = f(x_neg, Alpha_MLE, Beta_n_MLE, Lambda_n_MLE);

data = [pos_trans; neg_trans];

%% transform evaluation points
x_eval = x_eval - median(x);
pos_inds_eval = find(x_eval>=0);
neg_inds_eval = find(x_eval<0);

jacobian_p = @(x) Alpha_MLE*(1-Lambda_p_MLE*tanh(Beta_p_MLE*x)^2)*sech(Beta_p_MLE*x)^(Lambda_p_MLE-1); % jacobian term for x>median
jacobian_n = @(x) Alpha_MLE*(1-Lambda_n_MLE*tanh(Beta_n_MLE*x)^2)*sech(Beta_n_MLE*x)^(Lambda_n_MLE-1); % jacobian term for x<median

if strcmp('double',double)==1
    x_eval_t = zeros(1,length(x_eval));
    x_eval_t(pos_inds_eval) = log(1+x_eval(pos_inds_eval));
    x_eval_t(neg_inds_eval) = -log(1-x_eval(neg_inds_eval));
    x_eval_tt = zeros(1,length(x_eval));
    x_eval_tt(pos_inds_eval) = f(x_eval_t(pos_inds_eval),Alpha_MLE,Beta_p_MLE,Lambda_p_MLE);
    x_eval_tt(neg_inds_eval) = f(x_eval_t(neg_inds_eval),Alpha_MLE,Beta_n_MLE,Lambda_n_MLE);

    jacobian1 = zeros(1,length(x_eval));
    jacobian2 = zeros(1,length(x_eval));
    
    jacobian1(pos_inds_eval) = Alpha_MLE*(1-Lambda_p_MLE*tanh(Beta_p_MLE*x_eval_t(pos_inds_eval)).^2).*sech(Beta_p_MLE*x_eval_t(pos_inds_eval)).^(Lambda_p_MLE-1);
    jacobian1(neg_inds_eval) = Alpha_MLE*(1-Lambda_n_MLE*tanh(Beta_n_MLE*x_eval_t(neg_inds_eval)).^2).*sech(Beta_n_MLE*x_eval_t(neg_inds_eval)).^(Lambda_n_MLE-1);
    jacobian2(pos_inds_eval) = 1./(1+x_eval(pos_inds_eval));
    jacobian2(neg_inds_eval) = 1./(1-x_eval(neg_inds_eval));
    
    [pdf,~,bw] = ksdensity(data,x_eval_tt,'Kernel','normal');
    pdf = pdf.*abs(jacobian1).*abs(jacobian2);
    
elseif strcmp('pos_skew',double)
    x_eval_t = zeros(1,length(x_eval));
    x_eval_t = log(1+x_eval-min_xcent+shift)-med_xcent;
    
    x_eval_tt = zeros(1,length(x_eval));
    x_eval_tt(pos_inds_eval) = f(x_eval_t(pos_inds_eval),Alpha_MLE,Beta_p_MLE,Lambda_p_MLE);
    x_eval_tt(neg_inds_eval) = f(x_eval_t(neg_inds_eval),Alpha_MLE,Beta_n_MLE,Lambda_n_MLE);
    
    jacobian1 = zeros(1,length(x_eval));
    jacobian2 = zeros(1,length(x_eval));
    
    jacobian1(pos_inds_eval) = Alpha_MLE*(1-Lambda_p_MLE*tanh(Beta_p_MLE*x_eval_t(pos_inds_eval)).^2).*sech(Beta_p_MLE*x_eval_t(pos_inds_eval)).^(Lambda_p_MLE-1);
    jacobian1(neg_inds_eval) = Alpha_MLE*(1-Lambda_n_MLE*tanh(Beta_n_MLE*x_eval_t(neg_inds_eval)).^2).*sech(Beta_n_MLE*x_eval_t(neg_inds_eval)).^(Lambda_n_MLE-1);
    jacobian2 = 1./(1+x_eval-min_xcent+shift);
    
    [pdf,~,bw] = ksdensity(data,x_eval_tt,'Kernel','normal');
    pdf = pdf.*abs(jacobian1).*abs(jacobian2);
    
elseif strcmp('neg_skew',double)
    x_eval_t = zeros(1,length(x_eval));
    x_eval_t = -log(1-x_eval+max_xcent+shift)-med_xcent;
    x_eval_tt = zeros(1,length(x_eval));
    x_eval_tt(pos_inds_eval) = f(x_eval_t(pos_inds_eval),Alpha_MLE,Beta_p_MLE,Lambda_p_MLE);
    x_eval_tt(neg_inds_eval) = f(x_eval_t(neg_inds_eval),Alpha_MLE,Beta_n_MLE,Lambda_n_MLE);
    
    jacobian1 = zeros(1,length(x_eval));
    jacobian2 = zeros(1,length(x_eval));
    
    jacobian1(pos_inds_eval) = Alpha_MLE*(1-Lambda_p_MLE*tanh(Beta_p_MLE*x_eval_t(pos_inds_eval)).^2).*sech(Beta_p_MLE*x_eval_t(pos_inds_eval)).^(Lambda_p_MLE-1);
    jacobian1(neg_inds_eval) = Alpha_MLE*(1-Lambda_n_MLE*tanh(Beta_n_MLE*x_eval_t(neg_inds_eval)).^2).*sech(Beta_n_MLE*x_eval_t(neg_inds_eval)).^(Lambda_n_MLE-1);
    jacobian2 = 1./(1-x_eval+max_xcent+shift);
    
    [pdf,~,bw] = ksdensity(data,x_eval_tt,'Kernel','normal');
    pdf = pdf.*abs(jacobian1).*abs(jacobian2);
else
    x_eval_tt = zeros(1,length(x_eval));
    x_eval_tt(pos_inds_eval) = f(x_eval(pos_inds_eval),Alpha_MLE,Beta_p_MLE,Lambda_p_MLE);
    x_eval_tt(neg_inds_eval) = f(x_eval(neg_inds_eval),Alpha_MLE,Beta_n_MLE,Lambda_n_MLE);
        
    jacobian1 = zeros(1,length(x_eval));
    jacobian1(pos_inds_eval) = Alpha_MLE*(1-Lambda_p_MLE*tanh(Beta_p_MLE*x_eval(pos_inds_eval)).^2).*sech(Beta_p_MLE*x_eval(pos_inds_eval)).^(Lambda_p_MLE-1);
    jacobian1(neg_inds_eval) = Alpha_MLE*(1-Lambda_n_MLE*tanh(Beta_n_MLE*x_eval(neg_inds_eval)).^2).*sech(Beta_n_MLE*x_eval(neg_inds_eval)).^(Lambda_n_MLE-1);
    
    [pdf,~,bw] = ksdensity(data,x_eval_tt,'Kernel','normal');
    pdf = pdf.*abs(jacobian1);
end
    if strcmp('cdf',cdf)==1
        yu = ksdensity(data,x_eval_tt,'function','cdf','Kernel','normal','bandwidth',bw);
    else
        yu = [];
    end
end


function loglik = evaluate_loglik(x,BL_t)
    % performs MLE via minimising the negative loglikelihood
    Beta = exp(BL_t(1));
    
    % transformed lambda 
    if BL_t(2) == Inf
        lambda = 1;
    else
        lambda = (exp(BL_t(2))-1)/(1+exp(BL_t(2)));
    end    
    
    Alpha = ( (1/length(x)) * sum(((sinh(Beta*x).*sech(Beta*x).^lambda)/Beta).^2) )^(-0.5);
    loglik = 0.5*sum((Alpha*sinh(Beta*x).*...
        (sech(Beta*x).^(lambda))/Beta).^2)-...
        length(x)*log(Alpha)-...
        sum(log(1-(lambda)*tanh(Beta*x).^2))-...
        ((lambda)-1)*sum(log(sech(Beta*x)));
end























