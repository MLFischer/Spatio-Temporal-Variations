% find change points and plot
% based on MRES blog Martin H Trauth
% Adapted by Fischer ML
% 2022/11

J = customcolormap([0 0.5 1], {'#ffffff','#ff0000','#000000'});
colorbar; colormap(J); axis off;

%%
maxnumchg = num_chp;
ipoint = findchangepts(y,...
    'Statistic',stat_val,...
    'MaxNumChanges',maxnumchg);
ipoint(end+1) = length(y);
m(1) = mean(y(1:ipoint(1)));
for i = 1 : length(ipoint)-1
    m(i+1) = mean(y(ipoint(i):ipoint(i+1)));
end

m(length(ipoint)+1) = mean(y(ipoint(i):end));
%%
figure('Position',[00 800 1200 300],...
    'Color',[1 1 1]);
axes(...
    'Box','On',...
    'LineWidth',1,...
    'FontSize',14); hold on
line(t,y,...
    'LineWidth',1);
%%
line([0 t(ipoint(1))],[m(1) m(1)],...
    'LineWidth',1,'Color',[0.8510 0.3255 0.0980])

for i = 1 : length(ipoint)-1
    line([t(ipoint(i)) t(ipoint(i+1))],...
    [m(i+1) m(i+1)],...
    'LineWidth',1,'Color',[0.8510 0.3255 0.0980])
end

xline(t(ipoint(1:end-1)),...
    'LineWidth',1,...
    'Color',[0.8510 0.3255 0.0980]);