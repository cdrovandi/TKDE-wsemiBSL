function [theta_tilde] = para_transformation(theta,a,b)

theta_tilde = log((theta-a) ./ (b-theta));

end

