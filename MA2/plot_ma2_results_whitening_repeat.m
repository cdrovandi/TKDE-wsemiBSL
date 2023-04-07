

load('exact_results.mat');
thin = 15;


t = tiledlayout(2,4,'TileSpacing','compact','Padding','tight');

nexttile
load('results_whitening_repeat_nothing.mat')

title('\epsilon = 0, \delta = 1','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,1));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,1));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_1','FontSize',14);
ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_skewness.mat')


title('\epsilon = 1, \delta = 1','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,1));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,1));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_1','FontSize',14);
%ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_kurtosis.mat')


title('\epsilon = 0, \delta = 0.6','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,1));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,1));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_1','FontSize',14);
%ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_skewness_and_kurtosis.mat')


title('\epsilon = 1, \delta = 2','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,1));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,1));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_1','FontSize',14);
%ylabel('density','FontSize',14);




nexttile
load('results_whitening_repeat_nothing.mat')


%title('\epsilon = 0, \delta = 1','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,2));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,2));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_2','FontSize',14);
ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_skewness.mat')


%title('\epsilon = 1, \delta = 1','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,2));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,2));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_2','FontSize',14);
%ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_kurtosis.mat')


%title('\epsilon = 0, \delta = 0.6','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,2));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,2));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_2','FontSize',14);
%ylabel('density','FontSize',14);


nexttile
load('results_whitening_repeat_skewness_and_kurtosis.mat')


%title('\epsilon = 1, \delta = 2','FontSize',14);
hold on;
for i = 1:24
    theta_rep = theta_all{i};
    [f,xi] = ksdensity(theta_rep(1:thin:end,1));
    plot(xi,f,'Color',[0.5 0.5 0.5],'LineWidth',1);
end
[f,xi] = ksdensity(theta(1:thin:end,1));

plot(xi,f,'k','LineWidth',2);
xlabel('\theta_2','FontSize',14);
%ylabel('density','FontSize',14);



