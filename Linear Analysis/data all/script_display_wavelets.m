% based on MRES http://mres.uni-potsdam.de/index.php/2017/01/31/calculating-the-continuous-1-d-wavelet-transform/


trez = 1/t;
[wt,f,coi] = cwt(series3L,1/agemodelres);
coi_r = 1./coi;

figure('Position',[100 300 800 300],...
   'Color',[1 1 1]);
pcolor(t,1./f,abs(wt)), shading interp
hold on
if testrun  == 1
text(0,2,datastr( xrfindex(i) ) ,...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
end
line(t,coi_r,'Color','w',...
    'LineStyle','--',...
    'LineWidth',2)
xlabel('Time (kyr)')
ylabel('Period (kyr)')
title('Wavelet Power Spectrum')
set(gcf,'Colormap',jet)
set(gca,'XLim',[agemodelmin agemodelmax],...
   'YLim',[10 130])
colorbar
hold off