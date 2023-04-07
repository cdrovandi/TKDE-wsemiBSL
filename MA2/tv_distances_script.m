% compute tv distances for first plot

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

tv_distance = zeros(4,4);
%%
fprintf('1\n')

load('50_e0_d1_semiBSL_V1_KDE_NoWhitening_n800.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,1)  = 0.5*tv_bsl*len*wid;

load('50_e0_d1_semiBSL_V1_KDE_Whitening_n80.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,2)  = 0.5*tv_bsl*len*wid;

load('50_e0_d1_semiBSL_V1_HP_NoWhitening_n800.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,3)  = 0.5*tv_bsl*len*wid;

load('50_e0_d1_semiBSL_V1_HP_Whitening_n255.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(1,4)  = 0.5*tv_bsl*len*wid;



%%
fprintf('2\n')

load('50_e-1_d1_semiBSL_V1_KDE_NoWhitening_n700.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,1)  = 0.5*tv_bsl*len*wid;

load('50_e-1_d1_semiBSL_V1_KDE_Whitening_n90.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,2)  = 0.5*tv_bsl*len*wid;

load('50_e-1_d1_semiBSL_V1_HP_NoWhitening_n700.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,3)  = 0.5*tv_bsl*len*wid;

load('50_e-1_d1_semiBSL_V1_HP_Whitening_n240.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(2,4)  = 0.5*tv_bsl*len*wid;


%%
fprintf('3\n')

load('50_e0_d0.6_semiBSL_V1_KDE_NoWhitening_n1000.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(3,1)  = 0.5*tv_bsl*len*wid;

load('50_e0_d0.6_semiBSL_V1_KDE_Whitening_n330.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(3,2)  = 0.5*tv_bsl*len*wid;

load('50_e0_d0.6_semiBSL_V1_HP_NoWhitening_n1000.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(3,3)  = 0.5*tv_bsl*len*wid;

load('50_e0_d0.6_semiBSL_V1_HP_Whitening_n275.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(3,4)  = 0.5*tv_bsl*len*wid;


%%
fprintf('4\n')

load('50_e1_d2_semiBSL_V1_KDE_NoWhitening_n700.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(4,1)  = 0.5*tv_bsl*len*wid;

load('50_e1_d2_semiBSL_V1_KDE_Whitening_n50.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(4,2)  = 0.5*tv_bsl*len*wid;

load('50_e1_d2_semiBSL_V1_HP_NoWhitening_n700.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(4,3)  = 0.5*tv_bsl*len*wid;

load('50_e1_d2_semiBSL_V1_HP_Whitening_n250.mat')
[f_bsl,~] = ksdensity(theta,[X1 X2]);
tv_bsl = sum(abs(f_true - f_bsl));
tv_distance(4,4)  = 0.5*tv_bsl*len*wid;

%%

save('tvd.mat','tv_distance')
