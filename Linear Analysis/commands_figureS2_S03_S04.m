
% KL 15 Analysis Script
% Figure 02 KL 15 show quality flags
% Markus L Fischer
% 2023/09

clear, clc, close all

% get subfolders
folder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(folder));


ds = get(0,'ScreenSize');


% Read the tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
datakl15xrf     = readtable('data_KL15_XRF.txt');
datakl15agem    = readtable('data_KL15-2_smooth_spline_ages.txt');
datakl15qf      = readtable('data_KL15_qf.txt');

% Rename the ID colum, as this is the depth and the match
datakl15xrf     = renamevars(datakl15xrf, "User_ID","depth");
datakl15agem    = renamevars(datakl15agem, "best","age");

% Do outerjoin, save the variables and get back to array
T               = outerjoin(datakl15agem,datakl15xrf);
datakl15string  = convertCharsToStrings(T.Properties.VariableNames');
datakl15        = table2array(T);
datakl15qf        = table2array(datakl15qf);

% Remove rows without data.
datakl15(isnan(datakl15(:,6))==1,:) = [];

% a to ka
datakl15(:,4)   = datakl15(:,4)/1000;

%only age and counts:
datakl15string;
xrfindex                        =   find(contains(datakl15string,'_Area'))';  % get index numbers of XRF data 
datakl15string                  =   strrep(datakl15string,'_Area','');       % replace underscore with normal spacing


%data = datakl15(:,[4 xrfindex]);
data_string = datakl15string([4 xrfindex],:);


D = table2array(datakl15xrf);
data = D(:,[1 2 4 6 8 10 12 14 16 18 20 22 24 26 28]);



data(:,16)          = sum(data(:,[2:15]),2);
data(1:2118,17)     = diff(data(:,1));
data(2119,17)       = 0;
data(:,18)          = log(data(:,13)./data(:,14));
data(:,19)          = log(data(:,14)./data(:,9));

data_string(16)     = "XRF count sum";
data_string(17)     = "depth diff";
data_string(18)     = "log(Ca/Ti)";
data_string(19)     = "log(Ti/Al)";

% new column with values mean - 2*sd
data(:,20) = data(:,16).* (data(:,16)<(mean(data(:,16))-2*std(data(:,16))));

data((data(:,20)== 0),20) = nan ;

ld=1;



% Plot depth sum and depth differences
fig_s2 = figure('Position',[100 200 1200 400],...
    'Color',[1 1 1]);

a(1) = axes('Position',[0.1 0.7 0.8 0.25],...
    'Visible','Off');
t(1) = text(0.99,1.1,"XRF counts sum",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(1) = axes('Position',[0.1 0.7 0.8 0.25],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');

l(1) = line(data(:,1),data(:,16),"Marker",".");

l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),16), ...
    "LineStyle","none","Marker","diamond","Color","#A2142F",'MarkerSize',12);
l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),16) ...
    ,"LineStyle","none","Marker","square","Color","#A2142F",'MarkerSize',12);


xline(data(find(data(:,17)>1),1),"LineWidth",ld,"LineStyle",":","Color","black")

legend(["all values" "QF 1" "QF 2" "section break"])
ylim([100000 500000])

%yl(1) = ylabel('counts');


a(2) = axes('Position',[0.1 0.4 0.8 0.25],...
    'Visible','Off');
t(2) = text(0.99,1,data_string(18),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(2) = axes('Position',[0.1 0.4 0.8 0.25],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
l(2) = line(data(:,1),data(:,18));
%xline(data(find(data(:,17)>1),1),"LineWidth",ld,"LineStyle",":","Color","black")
l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),18), ...
    "LineStyle","none","Marker","diamond","Color","#A2142F",'MarkerSize',12);
l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),18) ...
    ,"LineStyle","none","Marker","square","Color","#A2142F",'MarkerSize',12);
%yl(2) = ylabel('counts');

a(3) = axes('Position',[0.1 0.1 0.8 0.25 ],...
    'Visible','Off');
t(3) = text(0.99,1,data_string(19),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(3) = axes('Position',[0.1 0.1 0.8 0.25],...
    'XGrid','On',...
    'Color','None');
%xline(data(find(data(:,17)>1),1),"LineWidth",ld,"LineStyle",":","Color","#7E2F8E")
l(3) = line(data(:,1),data(:,19));
%l(3) = line(data(:,1),data(:,18),"LineStyle","none","Marker","none","Color","#D95319","LineWidth",10);
l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==1,1)),19), ...
    "LineStyle","none","Marker","diamond","Color","#A2142F",'MarkerSize',12);
l(1) = line(data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),1) ...
    ,data(ismember( data(:,1),datakl15qf(datakl15qf(:,2)==2,1)),19) ...
    ,"LineStyle","none","Marker","square","Color","#A2142F",'MarkerSize',12);
%yl(3) = ylabel('counts');
xl(2) = xlabel('Depth (cm)');

linkaxes(b,'x')

z(1) = zoom(fig_s2);
set(z(1),'Motion','horizontal');
p(1) = pan(fig_s2);
set(p(1),'Motion','horizontal');
xlim([0 2400])

hold off


%export figure as vector grafic
exportgraphics(fig_s2,'Figure S02.pdf','ContentType','vector')

% Plot all Elements(1):
sel = [2 3 4 5 6 7 8 9];

fig_s3 = figure('Position',[0 00 800 1200],...
    'Color',[1 1 1]);


a(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'Visible','Off');
t(1) = text(0.99,0.9,data_string(sel(1)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")
l(1) = line(data(:,1),data(:,sel(1)));
%l(1) = line(data(:,1),data(:,sel(1)),"LineStyle","none","Marker",".","Color","#D95319");


a(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'Visible','Off');
t(2) = text(0.99,0.9,data_string(sel(2)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")
l(2) = line(data(:,1),data(:,sel(2)));
%l(2) = line(data(:,1),data(:,sel(2)),"LineStyle","none","Marker",".","Color","#D95319");


a(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'Visible','Off');
t(3) = text(0.99,0.9,data_string(sel(3)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(3) = line(data(:,1),data(:,sel(3)));
%l(3) = line(data(:,1),data(:,sel(3)),"LineStyle","none","Marker",".","Color","#D95319");


a(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'Visible','Off');
t(4) = text(0.99,0.9,data_string(sel(4)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(4) = line(data(:,1),data(:,sel(4)));
%l(4) = line(data(:,1),data(:,sel(4)),"LineStyle","none","Marker",".","Color","#D95319");


a(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'Visible','Off');
t(5) = text(0.99,0.9,data_string(sel(5)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(5) = line(data(:,1),data(:,sel(5)));
%l(5) = line(data(:,1),data(:,sel(5)),"LineStyle","none","Marker",".","Color","#D95319");


a(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'Visible','Off');
t(6) = text(0.99,0.9,data_string(sel(6)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(6) = line(data(:,1),data(:,sel(6)));
%l(6) = line(data(:,1),data(:,sel(6)),"LineStyle","none","Marker",".","Color","#D95319");


a(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'Visible','Off');
t(7) = text(0.99,0.9,data_string(sel(7)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(7) = line(data(:,1),data(:,sel(7)));
%l(7) = line(data(:,1),data(:,sel(7)),"LineStyle","none","Marker",".","Color","#D95319");


a(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'Visible','Off');
t(8) = text(0.99,0.9,data_string(sel(8)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(8) = line(data(:,1),data(:,sel(8)));
%l(8) = line(data(:,1),data(:,sel(8)),"LineStyle","none","Marker","o","Color","#D95319");

xl(8) = xlabel('Depth (cm)');


linkaxes(b,'x')

z(1) = zoom(fig_s3);
set(z(1),'Motion','horizontal');
p(1) = pan(fig_s3);
set(p(1),'Motion','horizontal');
xlim([0 2250])

hold off

exportgraphics(fig_s3,'Figure S03.pdf','ContentType','vector')


%
sel = [10 12 13 14 15 16 18 19];

fig_s4 = figure('Position',[0 00 800 1200],...
    'Color',[1 1 1]);


a(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'Visible','Off');
t(1) = text(0.99,0.9,data_string(sel(1)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(1) = axes('Position',[0.1 0.855 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")
l(1) = line(data(:,1),data(:,sel(1)));
%l(1) = line(data(:,1),data(:,sel(1)),"LineStyle","none","Marker",".","Color","#D95319");


a(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'Visible','Off');
t(2) = text(0.99,0.9,data_string(sel(2)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(2) = axes('Position',[0.1 0.745 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")
l(2) = line(data(:,1),data(:,sel(2)));
%l(2) = line(data(:,1),data(:,sel(2)),"LineStyle","none","Marker",".","Color","#D95319");


a(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'Visible','Off');
t(3) = text(0.99,0.9,data_string(sel(3)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(3) = axes('Position',[0.1 0.635 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(3) = line(data(:,1),data(:,sel(3)));
%l(3) = line(data(:,1),data(:,sel(3)),"LineStyle","none","Marker",".","Color","#D95319");


a(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'Visible','Off');
t(4) = text(0.99,0.9,data_string(sel(4)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(4) = axes('Position',[0.1 0.525 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(4) = line(data(:,1),data(:,sel(4)));
%l(4) = line(data(:,1),data(:,sel(4)),"LineStyle","none","Marker",".","Color","#D95319");


a(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'Visible','Off');
t(5) = text(0.99,0.9,data_string(sel(5)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(5) = axes('Position',[0.1 0.415 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(5) = line(data(:,1),data(:,sel(5)));
%l(5) = line(data(:,1),data(:,sel(5)),"LineStyle","none","Marker",".","Color","#D95319");


a(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'Visible','Off');
t(6) = text(0.99,0.9,data_string(sel(6)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(6) = axes('Position',[0.1 0.305 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(6) = line(data(:,1),data(:,sel(6)));
%l(6) = line(data(:,1),data(:,sel(6)),"LineStyle","none","Marker",".","Color","#D95319");


a(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'Visible','Off');
t(7) = text(0.99,0.9,data_string(sel(7)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(7) = axes('Position',[0.1 0.195 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(7) = line(data(:,1),data(:,sel(7)));
%l(7) = line(data(:,1),data(:,sel(7)),"LineStyle","none","Marker",".","Color","#D95319");


a(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'Visible','Off');
t(8) = text(0.99,0.9,data_string(sel(8)),...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(8) = axes('Position',[0.1 0.085 0.8 0.1],...
    'XGrid','On',...
    'Color','None');
xline(data(find(data(:,17)>1),1),"LineWidth",1,"LineStyle",":","Color","#7E2F8E")

l(8) = line(data(:,1),data(:,sel(8)));
%l(8) = line(data(:,1),data(:,sel(8)),"LineStyle","none","Marker","o","Color","#D95319");

xl(8) = xlabel('Depth (cm)');


linkaxes(b,'x')

z(1) = zoom(fig_s4);
set(z(1),'Motion','horizontal');
p(1) = pan(fig_s4);
set(p(1),'Motion','horizontal');
xlim([0 2250])

hold off

exportgraphics(fig_s4,'Figure S04.pdf','ContentType','vector')




%delete local path
rmpath(genpath(folder));