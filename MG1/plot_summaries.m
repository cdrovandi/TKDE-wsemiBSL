% plot the distribution of summaries at true parameter value

load('ssx_summaries.mat')
figure
for i = 1:50
    [f,xi] = ksdensity(ssx_curr(:,i))
    plot(xi,f,'k');
    hold on
%     
%     [pdf,xi,~] = HP_trans_double_KDE(ssx_curr(:,i),linspace(-10,70,200),'single','cdf');
%     plot(xi,pdf,'Color',[0.6,0,0.8]);
%     
%     [pdf,xi,~] = HP_trans_double_KDE(ssx_curr(:,i),linspace(-10,70,200),'pos_skew','cdf');
%     plot(xi,pdf,'Color',[1,0.6,0]);
end

xlabel('s')
ylabel('Density')