function [theta] = para_back_transformation(theta_tilde,a,b)

e_theta_tilde = exp(theta_tilde);

theta = a ./ (1 + e_theta_tilde) + b ./ (1 + 1 ./ e_theta_tilde);

end

