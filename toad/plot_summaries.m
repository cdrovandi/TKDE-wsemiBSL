% plot distribution of summaries generated at true parameter value

load('ssx_summaries.mat')
fig2 = figure(2);

for i = 1:48
subaxis(6,8,i,'Spacing',0,'SpacingVert',0.06,'SpacingHoriz',0.035,'Padding',0,'marginL',0.08,'marginR',0.03,'mt',0.02,'mb',0.05)
[f,xi] = ksdensity(ssx(:,i));
plot(xi,f,'k','linewidth',1.5)

filename = strrep(strjoin(["s_{",i,"}"]),' ','')
xlabel(filename)

if mod(i,8) == 1
    ylabel('Density')
end

end


set(fig2, 'Position',  [100, 100, 1500, 950])
