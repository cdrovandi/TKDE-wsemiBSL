function y = sim_alphastable_statespace(theta,T)
% code to simulate from alpha-stable state space model

mu = theta(1);
phi = theta(2);
alpha = theta(3);
beta = theta(4);
gamma = theta(5);
delta = theta(6);
sigma = theta(7);
stable = makedist('Stable','alpha',alpha,'beta',beta,'gam',gamma,'delta',delta);
x = 0;
y = zeros(1,T);
v = random(stable,1,T);
for t = 1:T
    x = normrnd(mu + phi*(x-mu), sigma);
    y(t) = exp(0.5*x)*v(t);
end
end
