% Script to display the Chew Bahir and other data vs age
%
% The script allows zooming and paning linked axes.
%
% 21 Aug 2019 - Trauth

MIS =       [0   14  29   57  71   82   87   96   109 123  130 191 243 300 337  374  424  478  533  563];
MIS_label = ["1" "2" "3" "4" "5a" "b" "c" "d" "e" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15"];

MISp =       [0   14   29 57  71    130 191 243 300 337  374  424  478  533  563];
MIS_labelp = ["1" "2" "3" "4" "5"  "6" "7" "8" "9" "10" "11" "12" "13" "14" "15"];


HE = [12 16.8 24 31 38 45 60];
HE_label = ["H0" "H1" "H2" "H3" "H4" "H5" "H6"];

f(1) = figure('Position',[0 0 1200 1000],...
    'Color',[1 1 1]);

a(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'Visible','Off');
tstr(1,:) = strcat(datastr(varselectnum(1))," ",...
    datastr(varselectdem(1)));
t(1) = text(0.99,0.9,tstr(1,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
t(11) = text(0.99,1.30,"KL15 display major XRF results",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');


l(1) = line(data(:,1),data(:,varselectnum(1))./...
    data(:,varselectdem(1)));

if agerefs == 1 
        for i = 1:20
         if MIS(i) < agemodelmax 

               xline(MIS(i),':',{MIS_label(i)}, 'LabelOrientation', 'horizontal');
               
               % old try to just get triangles at the top of the plot °_°
               % YL = get(gca, 'YLim');
               % YLv  = repelem(YL(2),length(MIS));
               % plot(MIS,YLv,'-v','MarkerSize',10,...
               % 'MarkerEdgeColor','red',...
               % 'MarkerFaceColor',[1 .6 .6])
               % Markus and Matlabplot still no friends ;)
         end
        end
end

if agerefs == 2
        for i = 1:7
         if HE(i) < agemodelmax 
               xline(HE(i),':',{HE_label(i)}, 'LabelOrientation', 'horizontal');
         end
        end

end
if agerefs == 3
            for i = 1:20
         if MIS(i) < agemodelmax 
               xline(MIS(i),':',{MIS_label(i)}, 'LabelOrientation', 'horizontal');
         end
        end
            for i = 1:7
         if HE(i) < agemodelmax 
               xline(HE(i),':',{HE_label(i)}, 'LabelOrientation', 'horizontal');
         end
        end

end

a(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'Visible','Off');
tstr(2,:) = strcat(datastr(varselectnum(2))," ",...
    datastr(varselectdem(2)));
t(2) = text(0.99,0.9,tstr(2,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(2) = line(data(:,1),data(:,varselectnum(2))./...
    data(:,varselectdem(2)));

a(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'Visible','Off');
tstr(3,:) = strcat(datastr(varselectnum(3))," ",...
    datastr(varselectdem(3)));
t(3) = text(0.99,0.9,tstr(3,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(3) = line(data(:,1),data(:,varselectnum(3))./...
    data(:,varselectdem(3)));

a(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'Visible','Off');
tstr(4,:) = strcat(datastr(varselectnum(4))," ",...
    datastr(varselectdem(4)));
t(4) = text(0.99,0.9,tstr(4,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(4) = line(data(:,1),data(:,varselectnum(4))./...
    data(:,varselectdem(4)));

a(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'Visible','Off');
tstr(5,:) = strcat(datastr(varselectnum(5))," ",...
    datastr(varselectdem(5)));
t(5) = text(0.99,0.9,tstr(5,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(5) = line(data(:,1),data(:,varselectnum(5))./...
    data(:,varselectdem(5)));

a(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'Visible','Off');
tstr(6,:) = strcat(datastr(varselectnum(6))," ",...
    datastr(varselectdem(6)));
t(6) = text(0.99,0.9,tstr(6,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(6) = line(data(:,1),data(:,varselectnum(6))./...
    data(:,varselectdem(6)));

a(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'Visible','Off');
tstr(7,:) = strcat(datastr(varselectnum(7))," ",...
    datastr(varselectdem(7)));
t(7) = text(0.99,0.9,tstr(7,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(7) = line(data(:,1),data(:,varselectnum(7))./...
    data(:,varselectdem(7)));

a(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'Visible','Off');
tstr(8,:) = strcat(datastr(varselectnum(8))," ",...
    datastr(varselectdem(8)));
t(8) = text(0.99,0.9,tstr(8,:),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'XGrid','On',...
    'Color','None');
l(8) = line(data(:,1),data(:,varselectnum(8))./...
    data(:,varselectdem(8)));
xl(8) = xlabel('Age (kyrs BP)');

linkaxes(b,'x')

z(1) = zoom(f(1));
set(z(1),'Motion','horizontal');
p(1) = pan(f(1));
set(p(1),'Motion','horizontal');


