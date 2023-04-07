
load('tvd_tkde.mat')
thin = 15;
fig = figure;
subaxis(2,4,1,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.08,'marginR',0.20,'mt',0.03,'mb',0.17)
eps = 0;
delta = 0.1;
x_grid = linspace(-100,100,3000);
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
% xlabel('S(y_i)')
ylabel('\epsilon = 0, \delta = 0.1')


subaxis(2,4,2,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.01,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on
load('50_e0_d0.1_semiBSL_V1_KDE_n750.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
plot(0.6,0.2,'or','linewidth',2,'MarkerFaceColor', 'r')
text(0.2,0.8,'n = 750')
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
xlabel('\theta_1')
ylabel('\theta_2')
set(gca,'XTickLabel',[])
%title('semiBSL')
legend('Exact','semiBSL','location','southeast')
txt = strjoin({'tv =',num2str(round(tv_distance(1,1),2))});
text(0.2,0.7,txt)


subaxis(2,4,3,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.01,'mb',0.07)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on
load('50_e0_d0.1_semiBSL_V1_HP_n750.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0,1,0],'LineStyle','--');
plot(0.6,0.2,'or','linewidth',2,'MarkerFaceColor', 'r')
text(0.2,0.8,'n = 750')
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
% xlabel('\theta_1')
% ylabel('\theta_2')
set(gca,'XTickLabel',[],'YTickLabel',[])
legend('Exact','semiBSL TKDE3','location','southeast')
txt = strjoin({'tv =',num2str(round(tv_distance(1,2),2))});
text(0.2,0.7,txt)





subaxis(2,4,5,'Spacing',0,'SpacingVert',0.03,'SpacingHoriz',0.055,'Padding',0,'marginL',0.08,'marginR',0.20,'mt',0.08,'mb',0.12)
eps = 5;
delta = 0.4;
xmin = 0;
xmax = 1e6;
x_grid = linspace(xmin,xmax,10000);
plot(x_grid,normpdf(sinh(delta*(asinh(x_grid))-eps),0,1) .*abs((delta*cosh(eps - delta*asinh(x_grid))) ./ (sqrt(x_grid.^2+1))),'linewidth',1.5,'Color','k');
xlabel('s')
ylabel('\epsilon = 5, \delta = 0.4')


subaxis(2,4,6,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.09)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on
load('50_e5_d0.4_semiBSL_V1_KDE_n750.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
plot(0.6,0.2,'or','linewidth',2,'MarkerFaceColor', 'r')
text(0.2,0.8,'n = 750')
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
xlabel('\theta_1')
ylabel('\theta_2')
legend('Exact','semiBSL','location','southeast')
txt = strjoin({'tv =',num2str(round(tv_distance(2,1),2))});
text(0.2,0.7,txt)


subaxis(2,4,7,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.09)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on
load('50_e5_d0.4_semiBSL_V1_HP_n750.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[1,0.6,0],'LineStyle','--');
plot(0.6,0.2,'or','linewidth',2,'MarkerFaceColor', 'r')
text(0.2,0.8,'n = 750')
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
xlabel('\theta_1')
% ylabel('\theta_2')
set(gca,'YTickLabel',[])
legend('Exact','semiBSL TKDE1','location','southeast')
txt = strjoin({'tv =',num2str(round(tv_distance(2,2),2))});
text(0.2,0.7,txt)

set(fig, 'Position',  [100, 100, 1000, 600])

% plot(theta(:,1),'Color',[0,0.8,0.2])
% hold on, plot(theta(:,2),'Color',[1,0,0])
% xlabel('Iterations')



subaxis(2,4,8,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0,'Padding',0,'marginL',0.04,'marginR',0.03,'mt',0.03,'mb',0.09)
load('exact_results.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor','k','LineStyle','-');
hold on
load('50_e5_d0.4_semiBSL_V1_KDE_n20000.mat')
[kx,ky] = ksdensity(theta(1:thin:end,:));
Z = reshape(kx,[30,30]);
X = unique(ky(:,1));
Y = unique(ky(:,2));
contour(X,Y,Z,10,'LineWidth',1.3,'ShowText','Off','LineColor',[0.2,0.6,1],'LineStyle','--');
plot(0.6,0.2,'or','linewidth',2,'MarkerFaceColor', 'r')
text(0.2,0.8,'n = 20000')
%legend('Exact','KDE','Location','SouthEast')
xlim([0.15,0.95])
ylim([-0.05,0.85])
xlabel('\theta_1')
set(gca,'YTickLabel',[])
load('tvd_tkde_20000.mat')
legend('Exact','semiBSL','location','southeast')
txt = strjoin({'tv =',num2str(round(tv_distance(1,1),2))});
text(0.2,0.7,txt)


thin = 15;
fig2 = figure;
subaxis(2,3,1,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.07,'mb',0.17)
load('50_e0_d0.1_semiBSL_V1_KDE_n750.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\theta_1','\theta_2')


subaxis(2,3,2,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.07,'mb',0.17)
load('50_e0_d0.1_semiBSL_V1_HP_n750.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\theta_1','\theta_2')


subaxis(2,3,4,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.12,'mb',0.12)
load('50_e5_d0.4_semiBSL_V1_KDE_n750.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\theta_1','\theta_2')


subaxis(2,3,5,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.12,'mb',0.12)
load('50_e5_d0.4_semiBSL_V1_HP_n750.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\theta_1','\theta_2')

subaxis(2,3,6,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.12,'mb',0.12)
load('50_e5_d0.4_semiBSL_V1_KDE_n20000.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\theta_1','\theta_2')





set(fig2, 'Position',  [100, 100, 1400, 500])












