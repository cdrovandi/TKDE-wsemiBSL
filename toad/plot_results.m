% plot results in paper

load('tv_distance2.mat')
load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,1,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)
[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n50_g0.mat')
% load('test3.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

xlabel('\alpha')
ylabel('\eta')
% title('Warton')
xlim([1.4,1.9])
ylim([30,40])
txt = {'tv =',num2str(round(tv_distance(1,1),2))};
text(1.8,31.5,txt)
plot(1.7,35,'or','MarkerFaceColor','r')
legend('semiBSL, n = 500','wsemiBSL, \gamma = 0, n = 50','location','northwest')

% gamma = 0.22
load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,2,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)

[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n50_g0.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

plot(1.7,0.6,'or','MarkerFaceColor','r')

xlabel('\alpha')
ylabel('p_0')
% title('Warton')
xlim([1.4,1.9])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(1,2),2))};
text(1.8,0.565,txt)


load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,3,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.03,'mb',0.03)

[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n50_g0.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');

plot(35,0.6,'or','MarkerFaceColor','r')

xlabel('\eta')
ylabel('p_0')
% title('Warton')
xlim([30,38])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(1,3),2))};
text(36.5,0.565,txt)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,4,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)

[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n100_g0.3.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor','m','LineStyle','--');
plot(1.7,35,'or','MarkerFaceColor','r')

xlabel('\alpha')
ylabel('\eta')
legend('semiBSL, n = 500','wsemiBSL, \gamma = 0.3, n = 100','location','northwest')
% title('PCA Whitening')
xlim([1.4,1.9])
ylim([30,40])
txt = {'tv =',num2str(round(tv_distance(2,1),2))};
text(1.8,31.5,txt)

load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,5,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)

[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n100_g0.3.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor','m','LineStyle','--');
plot(1.7,0.6,'or','MarkerFaceColor','r')

xlabel('\alpha')
ylabel('p_0')
% title('PCA Whitening')
xlim([1.4,1.9])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(2,2),2))};
text(1.8,0.565,txt)

load('50_wsemiBSL_n500_g1.mat')
fig2 = figure(2);
subaxis(3,3,6,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.03)

[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n100_g0.3.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor','m','LineStyle','--');
plot(35,0.6,'or','MarkerFaceColor','r')

xlabel('\eta')
ylabel('p_0')
% title('PCA Whitening')
xlim([30,38])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(2,3),2))};
text(36.5,0.565,txt)
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subaxis(3,3,7,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)
load('50_wsemiBSL_n500_g1.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n250_g0.7.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),1:2));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
plot(1.7,35,'or','MarkerFaceColor','r')

xlabel('\alpha')
ylabel('\eta')
legend('semiBSL, n = 500','wsemiBSL, \gamma = 0.7, n = 250','location','northwest')
% title('PCA Whitening')
xlim([1.4,1.9])
ylim([30,40])
txt = {'tv =',num2str(round(tv_distance(3,1),2))};
text(1.8,31.5,txt)



subaxis(3,3,8,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)
load('50_wsemiBSL_n500_g1.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n250_g0.7.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),[1,3]));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
plot(1.7,0.6,'or','MarkerFaceColor','r')

xlabel('\alpha')
ylabel('p_0')
% title('PCA Whitening')
xlim([1.4,1.9])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(3,2),2))};
text(1.8,0.565,txt)





subaxis(3,3,9,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0.015,'marginL',0.08,'marginR',0.03,'mt',0.01,'mb',0.05)
load('50_wsemiBSL_n500_g1.mat')
[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_wsemiBSL_n250_g0.7.mat')
% load('test3.mat')

[kx,ky] = ksdensity(theta(2000:15:length(theta),2:3));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,15,'LineWidth',1.5,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');

plot(35,0.6,'or','MarkerFaceColor','r')

xlabel('\eta')
ylabel('p_0')
% title('Warton')
xlim([30,38])
ylim([0.55,0.64])
txt = {'tv =',num2str(round(tv_distance(3,3),2))};
text(36.5,0.565,txt)

set(fig2, 'Position',  [100, 100, 900, 800])
