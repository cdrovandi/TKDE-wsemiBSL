% compute tv distances required for the plot of results
load('mg1_data3_true_1M.mat')
theta = theta0_post;
theta0 = theta;

% theta 12000:end 1.4 -> 1.9
% theta 22000:end 30 -> 39
% theta 32000:end 0.56 -> 0.66

theta1_grid = linspace(0.5,2,300);
theta2_grid = linspace(3.8,6.3,300);
theta3_grid = linspace(0.15,0.4,300);

len1 = theta1_grid(2) - theta1_grid(1);
len2 = theta2_grid(2) - theta2_grid(1);
len3 = theta3_grid(2) - theta3_grid(1);

a1 = len1*len2;
a2 = len1*len3;
a3 = len2*len3;

[X1,X2] = meshgrid(theta1_grid, theta2_grid);
X1 = X1(:);
X2 = X2(:);

[X3,X4] = meshgrid(theta1_grid, theta3_grid);
X3 = X3(:);
X4 = X4(:);

[X5,X6] = meshgrid(theta2_grid, theta3_grid);
X5 = X5(:);
X6 = X6(:);

[f_true1,~] = ksdensity(theta0(2000:200:end,[1,2]),[X1 X2]);
[f_true2,~] = ksdensity(theta0(2000:200:end,[1,3]),[X3 X4]);
[f_true3,~] = ksdensity(theta0(2000:200:end,[2,3]),[X5 X6]);

tv_distance = zeros(3,3);
%%
load('100_mg1_KDE.mat')
[f_bsl,~] = ksdensity(theta(2000:end,[1,2]),[X1 X2]);
tv_bsl = sum(abs(f_true1 - f_bsl));
tv_distance(1,1)  = 0.5*tv_bsl*a1;

[f_bsl,~] = ksdensity(theta(2000:end,[1,3]),[X3 X4]);
tv_bsl = sum(abs(f_true2 - f_bsl));
tv_distance(1,2)  = 0.5*tv_bsl*a2;

[f_bsl,~] = ksdensity(theta(2000:end,[2,3]),[X5 X6]);
tv_bsl = sum(abs(f_true3 - f_bsl));
tv_distance(1,3)  = 0.5*tv_bsl*a3;


%%
load('100_mg1_HP_single.mat')
[f_bsl,~] = ksdensity(theta(2000:end,[1,2]),[X1 X2]);
tv_bsl = sum(abs(f_true1 - f_bsl));
tv_distance(2,1)  = 0.5*tv_bsl*a1;

[f_bsl,~] = ksdensity(theta(2000:end,[1,3]),[X3 X4]);
tv_bsl = sum(abs(f_true2 - f_bsl));
tv_distance(2,2)  = 0.5*tv_bsl*a2;

[f_bsl,~] = ksdensity(theta(2000:end,[2,3]),[X5 X6]);
tv_bsl = sum(abs(f_true3 - f_bsl));
tv_distance(2,3)  = 0.5*tv_bsl*a3;


%%
load('100_mg1_HP_pos_skew.mat')
[f_bsl,~] = ksdensity(theta(2000:end,[1,2]),[X1 X2]);
tv_bsl = sum(abs(f_true1 - f_bsl));
tv_distance(3,1)  = 0.5*tv_bsl*a1;

[f_bsl,~] = ksdensity(theta(2000:end,[1,3]),[X3 X4]);
tv_bsl = sum(abs(f_true2 - f_bsl));
tv_distance(3,2)  = 0.5*tv_bsl*a2;

[f_bsl,~] = ksdensity(theta(2000:end,[2,3]),[X5 X6]);
tv_bsl = sum(abs(f_true3 - f_bsl));
tv_distance(3,3)  = 0.5*tv_bsl*a3;





save('tv_distance2.mat','tv_distance')