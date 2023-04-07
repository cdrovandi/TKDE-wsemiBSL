% plot the figure in paper

load('tv_distance2.mat')

fig2 = figure(2);
subaxis(3,3,1,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)
rng(5)
load('mg1_data3_true_1M_thinned.mat')
[kx,ky] = ksdensity(theta0_post(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_KDE.mat')
% load('test3.mat')
[kx,ky] = ksdensity(theta(2000:15:end,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

xlabel('\theta_1')
ylabel('\theta_2')
% title('Warton')
xlim([0.5,2])
ylim([3.8,6.3])
txt = {'tv =',num2str(round(tv_distance(1,1),2))};
text(1.7,6,txt)
plot(1,5,'or','MarkerFaceColor','r')
legend('Exact','semiBSL','location','northwest')

% gamma = 0.22
fig2 = figure(2);
subaxis(3,3,2,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)

[kx,ky] = ksdensity(theta0_post(:,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_KDE.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:end,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

plot(1,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_1')
ylabel('\theta_3')
% title('Warton')
xlim([0.5,2])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(1,2),2))};
text(1.7,0.37,txt)


fig2 = figure(2);
subaxis(3,3,3,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)

[kx,ky] = ksdensity(theta0_post(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_KDE.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:end,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

plot(5,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_2')
ylabel('\theta_3')
% title('Warton')
xlim([3.8,6.3])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(1,3),2))};
text(5.6,0.37,txt)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig2 = figure(2);
subaxis(3,3,4,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)
rng(5)

[kx,ky] = ksdensity(theta0_post(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_single.mat')
[kx,ky] = ksdensity(theta(2000:15:end,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.6,0,0.8],'LineStyle','--');
plot(1,5,'or','MarkerFaceColor','r')

xlabel('\theta_1')
ylabel('\theta_2')
legend('Exact','semiBSL TKDE0','location','northwest')
% title('PCA Whitening')
xlim([0.5,2])
ylim([3.8,6.3])
txt = {'tv =',num2str(round(tv_distance(2,1),2))};
text(1.7,6,txt)

fig2 = figure(2);
subaxis(3,3,5,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)

[kx,ky] = ksdensity(theta0_post(:,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_single.mat')
[kx,ky] = ksdensity(theta(2000:15:end,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.6,0,0.8],'LineStyle','--');
plot(1,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_1')
ylabel('\theta_3')
% title('PCA Whitening')
xlim([0.5,2])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(2,2),2))};
text(1.7,0.37,txt)

fig2 = figure(2);
subaxis(3,3,6,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)

[kx,ky] = ksdensity(theta0_post(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_single.mat')
[kx,ky] = ksdensity(theta(2000:15:end,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.6,0,0.8],'LineStyle','--');
plot(5,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_2')
ylabel('\theta_3')
% title('PCA Whitening')
xlim([3.8,6.3])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(2,3),2))};
text(5.6,0.37,txt)
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subaxis(3,3,7,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)
rng(5)

[kx,ky] = ksdensity(theta0_post(:,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_pos_skew.mat')
[kx,ky] = ksdensity(theta(2000:15:end,1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[1,0.6,0],'LineStyle','--');
plot(1,5,'or','MarkerFaceColor','r')

xlabel('\theta_1')
ylabel('\theta_2')
legend('Exact','semiBSL TKDE1','location','northwest')
% title('PCA Whitening')
xlim([0.5,2])
ylim([3.8,6.3])
txt = {'tv =',num2str(round(tv_distance(3,1),2))};
text(1.7,6,txt)



subaxis(3,3,8,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)

[kx,ky] = ksdensity(theta0_post(:,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_pos_skew.mat')
[kx,ky] = ksdensity(theta(2000:15:end,[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[1,0.6,0],'LineStyle','--');
plot(1,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_1')
ylabel('\theta_3')
% title('PCA Whitening')
xlim([0.5,2])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(3,2),2))};
text(1.7,0.37,txt)





subaxis(3,3,9,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)

[kx,ky] = ksdensity(theta0_post(:,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('100_mg1_HP_pos_skew.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:end,2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[1,0.6,0],'LineStyle','--');

plot(5,0.2,'or','MarkerFaceColor','r')

xlabel('\theta_2')
ylabel('\theta_3')
% title('Warton')
xlim([3.8,6.3])
ylim([0.15,0.4])
txt = {'tv =',num2str(round(tv_distance(3,3),2))};
text(5.6,0.37,txt)

set(fig2, 'Position',  [100, 100, 900, 800])
