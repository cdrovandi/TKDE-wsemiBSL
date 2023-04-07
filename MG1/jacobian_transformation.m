function [J] = jacobian_transformation (theta_tilde,a,b)

e_theta_tilde = exp(theta_tilde);

logJ = log(b - a) - log(1 ./ e_theta_tilde + 2 + e_theta_tilde);
J = exp(sum(logJ));

end