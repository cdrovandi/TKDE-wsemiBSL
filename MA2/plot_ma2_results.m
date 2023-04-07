% plot results for paper

fig = figure;
thin = 20;

load('tvd.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 0;
delta = 1; % delta < 1 => HIGH kurtosis, delta > 1 => LOW kurtosis 
x_grid = linspace(-4,4,3000);
subaxis(4,5,1,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.05,'marginR',0.08,'mt',0.03,'mb',0.07)
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
% xlabel('s_j')
ylabel('\epsilon = 0, \delta = 1')
% title('\epsilon = 0, \delta = 1')
xlim([-4,4])
ylim([0,0.5])


subaxis(4,5,2,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')

[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_e0_d1_semiBSL_V1_KDE_NoWhitening_n800.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(1,1),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 800')
title('semiBSL')
set(gca,'XTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')



subaxis(4,5,3,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d1_semiBSL_V1_KDE_Whitening_n80.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(1,2),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 80')
title('wsemiBSL, \gamma = 0')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,4,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d1_semiBSL_V1_HP_NoWhitening_n800.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(1,3),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 800')
title('semiBSL TKDE')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,5,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d1_semiBSL_V1_HP_Whitening_n255.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(1,4),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 255')
title('wsemiBSL TKDE, \gamma = 0')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = -1;
delta = 1; % delta < 1 => HIGH kurtosis, delta > 1 => LOW kurtosis 
x_grid = linspace(-10,100,3000);
subaxis(4,5,6,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.05,'marginR',0.08,'mt',0.03,'mb',0.07)
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
% xlabel('s_j')
ylabel('\epsilon = -1, \delta = 1')
xlim([-10,2])
% title('\epsilon = 1.3, \delta = 0.6')


subaxis(4,5,7,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')

[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_e-1_d1_semiBSL_V1_KDE_NoWhitening_n700.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(2,1),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 700')
set(gca,'XTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,8,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e-1_d1_semiBSL_V1_KDE_Whitening_n90.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(2,2),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 90')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,9,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e-1_d1_semiBSL_V1_HP_NoWhitening_n700.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(2,3),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 700')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,10,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e-1_d1_semiBSL_V1_HP_Whitening_n240.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(2,4),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 240')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps = 0;
delta = 0.6; % delta < 1 => HIGH kurtosis, delta > 1 => LOW kurtosis 

x_grid = linspace(-20,20,3000);
subaxis(4,5,11,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.05,'marginR',0.08,'mt',0.03,'mb',0.07)
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
% xlabel('s_j')
ylabel('\epsilon = 0, \delta = 0.6')
% title('\epsilon = 0, \delta = 0.35')
xlim([-20,20])

subaxis(4,5,12,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')

[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_e0_d0.6_semiBSL_V1_KDE_NoWhitening_n1000.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(3,1),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 1000')
set(gca,'XTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,13,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d0.6_semiBSL_V1_KDE_Whitening_n330.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(3,2),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 330')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,14,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d0.6_semiBSL_V1_HP_NoWhitening_n1000.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
text(0.2,0.8,'n = 1000')
% % xlabel('\theta_1')
% % ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(3,3),2))});
text(0.2,0.7,txt)
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,15,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e0_d0.6_semiBSL_V1_HP_Whitening_n275.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(3,4),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 275')
set(gca,'XTickLabel',[],'YTickLabel',[])
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 1;
delta = 2; % delta < 1 => HIGH kurtosis, delta > 1 => LOW kurtosis 

x_grid = linspace(-50,500,3000);
subaxis(4,5,16,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.05,'marginR',0.08,'mt',0.03,'mb',0.07)
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
% xlabel('s_j')
ylabel('\epsilon = 1, \delta = 2')
% title('\epsilon = 5, \delta = 1')
xlim([-0.5,2])
xlabel('s')

subaxis(4,5,17,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')

[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');

hold on
load('50_e1_d2_semiBSL_V1_KDE_NoWhitening_n700.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(4,1),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 700')
% set(gca,'XTickLabel',[])
xlabel('\theta_1')
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,18,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e1_d2_semiBSL_V1_KDE_Whitening_n50.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(4,2),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 50')
set(gca,'YTickLabel',[])
xlabel('\theta_1')
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,19,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e1_d2_semiBSL_V1_HP_NoWhitening_n700.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(4,3),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 700')
set(gca,'YTickLabel',[])
xlabel('\theta_1')
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


subaxis(4,5,20,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on

load('50_e1_d2_semiBSL_V1_HP_Whitening_n250.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
%legend('Exact','TKDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% ylabel('\theta_2')
txt = strjoin({'tv =',num2str(round(tv_distance(4,4),2))});
text(0.2,0.7,txt)
text(0.2,0.8,'n = 250')
set(gca,'YTickLabel',[])
xlabel('\theta_1')
plot(0.6,0.2,'or','linewidth',1.5,'MarkerFaceColor', 'r')


set(fig, 'Position',  [100, 100, 1100, 900])


