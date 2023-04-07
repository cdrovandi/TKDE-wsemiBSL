% plot the figures in paper 

thin = 15;
fig = figure;


subaxis(2,3,1,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.07,'mb',0.17)
load('ssx_dataset1.mat')
for i = 1:50
%     [f,xi] = ksdensity(ssx(:,i),linspace(-20,20,200));
%     plot(xi,f,'r')
    [pdf,xi,~] = HP_trans_double_KDE(ssx(:,i),linspace(-20,20,200),'double','cdf');
    plot(xi,pdf,'k')
    hold on
end
xlim([-20,20])
xlabel('s')
ylabel('$$(\alpha,\beta)^\top = (1.2,0.5)^\top$$','interpreter','latex')
    


subaxis(2,3,2,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.07,'mb',0.17)
load('100_TKDE_stable.mat')
[f1,xi] = ksdensity(theta(1:thin:end,1));
plot(xi,f1,'linewidth',1.5,'Color',[0,1,0]) 

hold on
load('100_KDE_stable.mat')
[f2,xi] = ksdensity(theta(1:thin:end,1));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.6,1])

load('50_KDE_stable_n20000.mat')
[f2,xi] = ksdensity(theta(1:thin:end,1));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.9,1])

plot(1.2*ones(1,5),linspace(0,4,5),'--r','linewidth',1.5)
xlabel('\alpha')
ylabel('Density')
legend('semiBSL TKDE3, n = 2000','semiBSL, n = 2000','semiBSL, n = 20000','True parameter')

subaxis(2,3,3,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.07,'mb',0.17)
load('100_TKDE_stable.mat')
[f1,xi] = ksdensity(theta(1:thin:end,2));
plot(xi,f1,'linewidth',1.5,'Color',[0,1,0]) 

hold on
load('100_KDE_stable.mat')
[f2,xi] = ksdensity(theta(1:thin:end,2));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.6,1])

load('50_KDE_stable_n20000.mat')
[f2,xi] = ksdensity(theta(1:thin:end,2));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.9,1])


plot(0.5*ones(1,5),linspace(0,2.5,5),'--r','linewidth',1.5)
xlabel('\beta')
ylabel('Density')
% legend('semiBSL TKDE','semiBSL KDE','location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subaxis(2,3,4,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.12,'mb',0.12)
load('ssx_dataset2.mat')
for i = 1:50
%     [f,xi] = ksdensity(ssx(:,i),linspace(-20,20,200));
%     plot(xi,f,'r')
    [pdf,xi,~] = HP_trans_double_KDE(ssx(:,i),linspace(-20,20,200),'double','cdf');
    plot(xi,pdf,'k')
    hold on
end
xlim([-20,20])
xlabel('s')
ylabel('$$(\alpha,\beta)^\top = (0.7,0.5)^\top$$','interpreter','latex')



subaxis(2,3,5,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.12,'mb',0.12)
load('100_TKDE_stable_alpha0.7.mat')
[f1,xi] = ksdensity(theta(1:thin:end,1));
plot(xi,f1,'linewidth',1.5,'Color',[0,1,0]) 

hold on
load('100_KDE_stable_alpha0.7.mat')
[f2,xi] = ksdensity(theta(1:thin:end,1));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.6,1])

plot(0.7*ones(1,5),linspace(0,15,5),'--r','linewidth',1.5)
xlabel('\alpha')
ylabel('Density')
xlabel('\alpha')
ylabel('Density')
% legend('semiBSL TKDE','semiBSL KDE')



subaxis(2,3,6,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.12,'mb',0.12)
load('100_TKDE_stable_alpha0.7.mat')
[f1,xi] = ksdensity(theta(1:thin:end,2));
plot(xi,f1,'linewidth',1.5,'Color',[0,1,0]) 

hold on
load('100_KDE_stable_alpha0.7.mat')
[f2,xi] = ksdensity(theta(1:thin:end,2));
plot(xi,f2,'linewidth',1.5,'Color',[0.2,0.6,1])

plot(0.5*ones(1,5),linspace(0,4,5),'--r','linewidth',1.5)
xlabel('\beta')
ylabel('Density')

xlabel('\beta')
ylabel('Density')
% legend('semiBSL TKDE','semiBSL KDE','location','northwest')




set(fig, 'Position',  [100, 100, 1500, 900])






thin = 15;
fig2 = figure;
subaxis(2,3,1,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.07,'mb',0.17)
load('100_TKDE_stable.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\alpha','\beta')
title('semiBSL TKDE3, n = 2000')

subaxis(2,3,2,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.07,'mb',0.17)
load('100_KDE_stable.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\alpha','\beta')
title('semiBSL, n = 2000')


subaxis(2,3,3,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.07,'mb',0.17)
load('50_KDE_stable_n20000.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\alpha','\beta')
title('semiBSL, n = 20000')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subaxis(2,3,4,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.05,'marginR',0.03,'mt',0.12,'mb',0.12)
load('100_TKDE_stable_alpha0.7.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\alpha','\beta')

subaxis(2,3,5,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.05,'Padding',0,'marginL',0.065,'marginR',0.03,'mt',0.12,'mb',0.12)
load('100_KDE_stable_alpha0.7.mat')
plot(theta(:,1),'Color',[0,0.8,0.2])
hold on, plot(theta(:,2),'Color',[1,0,0])
xlabel('Iterations')
ylabel('Parameter Value')
legend('\alpha','\beta')
set(fig2, 'Position',  [100, 100, 1400, 500])
