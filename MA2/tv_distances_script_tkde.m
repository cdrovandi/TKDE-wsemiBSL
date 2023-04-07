% compute tv distances for second plot

load('exact_results.mat')
theta0 = theta;

% len = 0.01;
% wid = 0.01;
theta1_grid = linspace(0.2,1.2,500);
theta2_grid = linspace(-0.2,1,500);
len = theta1_grid(2)-theta1_grid(1);
wid = theta2_grid(2)-theta2_grid(1);
[X1,X2] = meshgrid(theta1_grid, theta2_grid);
X1 = X1(:);
X2 = X2(:);

[f_true,~] = ksdensity(theta0,[X1 X2]);

tv_distance = zeros(2,2);
%%
fprintf('1\n')

 load('50_e0_d0.1_semiBSL_V1_KDE_n750.mat')
%load('50_e5_d0.4_semiBSL_V1_KDE_n20000.mat')

[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,1)  = 0.5*tv_bsl*len*wid;

load('50_e0_d0.1_semiBSL_V1_HP_n750.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,2)  = 0.5*tv_bsl*len*wid;

%%
fprintf('2\n')

load('50_e5_d0.4_semiBSL_V1_KDE_n750.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,1)  = 0.5*tv_bsl*len*wid;

load('50_e5_d0.4_semiBSL_V1_HP_n750.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,2)  = 0.5*tv_bsl*len*wid;


save('tvd_tkde.mat','tv_distance')




load('50_e5_d0.4_semiBSL_V1_KDE_n20000.mat')

[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance  = 0.5*tv_bsl*len*wid;

save('tvd_tkde_20000.mat','tv_distance')



