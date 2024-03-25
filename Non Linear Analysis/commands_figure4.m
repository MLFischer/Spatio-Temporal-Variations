%% Based on MATLAB Code of Trauth et al., Quaternary Science Reviews, 2021
% KL15 RQA by Markus L. Fischer, 2023

clear, clc, close all

% get subfolders
folder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(folder));

% Load data
agemodelmin = 5; 
agemodelres = .1;
agemodelmax = 540; %540
inttype     = 'pchip';      

script_read_data_KL15_odp967
script_read_data_KL15_XRF
script_read_data_KL15_CHB_K_RRM_MHT500


data                            =   datakl15age; % all XRF   1:29
datastr                         =   datakl15string([4 7:34],:);       % create a name string of the data
data(:,size(data,2)+1)          = log(data(:,26)./data(:,16));
datastr(size(datastr,1)+1,:)    = "KL15 log(Ti/Al)";
data(:,size(data,2)+1)          = log(data(:,26)./data(:,24));
datastr(size(datastr,1)+1,:)    = "KL15 log(Ti/Ca)";


% Define time and data vector.
t    = 1*data(:,1); t = t';
t_bu = t;
x1   = data(:,31); x1 = x1';
x1_s = "KL15 log(Ti/Ca) Aridity Index";

x2   = data(:,30); x2 = x2';
x2_s = "KL 15 Aridity Index";

x3 = data_chb_age_cb0103_rr_560_rr(:,2)';
x3inv= -1 .* x3;
x3_s = "Chew Bahir Aridity Index";

x4 = dataodp967age(:,5)';
x4_s = dataodp967string(4);

ind2 = [43 156 270 373 440 535];

clear datakl15age datakl15qf datakl15string inttype
clear data_chb_stringv2 data_chb_stringv3 data_chb_stringv4 data_chb_stringv5
clear data_chb_agev2 data_chb_agev3 data_chb_agev4 data_chb_agev5
clear data datastr
clear dataodp967age dataodp967string data_chb_stringfin data_chb_agefin
%

%
% x2 = KL15 Ti/Al
x = x2;
x_s = x2_s;


xmin =  0;
xmax =  550;
ymin = 0;
ymax = 1.5;
ydir = 'normal';      


% Define embedding dimension m and time delay tau.
m   =  5;
tau =  10;
timespan_diff = tau*(m-1);

% Optional filter to remove trend.
filteroption = 0;
filterorder  = 5;
filtercutoff = 1;

theiler = tau;
l_min = 10 ;

% Define window size (w) and step size (ws).
w  = 500;
ws = 10;

% Set correction for border effects on (line_correct = 1) or off
% (line_correct = 0)
line_correct = 1;

threshold_calculation = 'fan';
% Choose norm from ['euc','max']
norm = 'euc';

% Data pretreatment, compute RQ, perform RQA and display results.
script_pretreatment_age
script_noisefilter_age
script_recurrenceplot_age
script_recurrencequantificationanalysis_age


figure1 = figure('Units','normalized',...
    'Position',[0 0 .5 1],...
    'Color',[1 1 1]);

axes1(1) = axes('Units','normalized',...
        'Position',[0.3625 0.7 0.3 0.3],...
    'XLim',[agemodelmin agemodelmax],...
    'YLim',[agemodelmin agemodelmax],...
    'LineWidth',1.5,...
    'Box','On',...
    'Layer','top',...
    'XGrid','On', ...
    'YTickLabel',''); hold on
axis square xy
imagesc(t,t,RR)
title(x_s)
colormap([1 1 1; 0 0 0])


y = [5  5    156  156 43  43 ];
x = [43 156  161  5   5   43 ];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    373  373 270 270];
x = [270 373  373  5   5   270];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    540  540 440 440];
x = [440 540  540  5   5   440];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')


line([0 150], [150 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([150 150], [0 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 

line([0 400], [400 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([400 400], [0 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 

sl = 350;
el = 450;
line([el el], [sl el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#7E2F8E","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#7E2F8E","LineWidth",2.5); % 


sl = 200;
el = 300;
line([el el], [sl el],'Color',"#D95319","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#D95319","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#D95319","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#D95319","LineWidth",2.5); % 

sl = 32;
el = 132;
line([el el], [sl el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#0072BD","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#0072BD","LineWidth",2.5); % 


xlabel('Time (kyr BP)')



b(1) = axes('Position',[0.1 0.6 0.8 0.04],...
    'Visible','Off');

text(0.99,1,"Recurrence Rate",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top','FontSize',12,'FontWeight','bold');


a(1) = axes('Position',[0.1 0.54 0.8 0.06],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","right");

ylabel("KL15")

yt = r_win_e(1,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

line(t,r_win_e(1,:),...
    'LineWidth',0.5, ...
    "color","black")




b(2) = axes('Position',[0.1 0.4 0.8 0.05],...
    'Visible','Off');

t(3) = text(0.99,1,"Determinism",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top','FontSize',12,'FontWeight','bold');

a(2) = axes('Position',[0.1 0.34 0.8 0.06],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","right");

ylabel("KL15")



yt = r_win_e(2,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);


line(t,r_win_e(2,:),...
    'LineWidth',.5, ...
    "color","black")



RQA_KL15 =  r_win_e(1:2,:);

% Chew Bahir Section
% x2 = K 
x = x3;
x_s = x3_s;


xmin =  0;
xmax =  550;
ymin = 0;
ymax = 1.5;
ydir = 'normal';      

script_pretreatment_age
script_noisefilter_age
script_recurrenceplot_age
script_recurrencequantificationanalysis_age


axes1(13) = axes('Units','normalized',...
            'Position',[0.675 0.7 0.3 0.3],...
    'XLim',[agemodelmin agemodelmax],...
    'YLim',[agemodelmin agemodelmax],...
    'LineWidth',1.5,...
    'Box','On',...
    'Layer','top',...
    'XGrid','On', ...
    'YTickLabel',''); hold on
axis square xy
imagesc(t,t,RR)
title(x_s)
colormap([1 1 1; 0 0 0])

y = [5  5    156  156 43  43 ];
x = [43 156  161  5   5   43 ];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    373  373 270 270];
x = [270 373  373  5   5   270];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    540  540 440 440];
x = [440 540  540  5   5   440];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

line([0 150], [150 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([150 150], [0 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 

line([0 400], [400 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([400 400], [0 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 

sl = 350;
el = 450;
line([el el], [sl el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#7E2F8E","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#7E2F8E","LineWidth",2.5); % 


sl = 200;
el = 300;
line([el el], [sl el],'Color',"#D95319","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#D95319","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#D95319","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#D95319","LineWidth",2.5); % 

sl = 32;
el = 132;
line([el el], [sl el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#0072BD","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#0072BD","LineWidth",2.5); % 



xlabel('Time (kyr)')

a(3) = axes('Position',[0.1 0.48 0.8 0.06],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","left");



ylabel("Chew Bahir")


yt = r_win_e(1,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);


line(t,r_win_e(1,:),...
    'LineWidth',.5, ...
    "color","black")

a(4) = axes('Position',[0.1 0.28 0.8 0.07],...
    'XGrid','Off',...
     "YAxisLocation","left", ...
    'YDir', 'normal', ...
    'Color','None');

ylabel("Chew Bahir")
xlabel("Time (kyr)")

yt = r_win_e(2,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt)-0.025;
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

line(t,r_win_e(2,:),...
    'LineWidth',.5, ...
    "color","black")

RQA_CHB =  r_win_e(1:2,:);

%
% ODP 967
%

x = x4;
x_s = x4_s;


xmin =  0;
xmax =  550;
ymin = 0;
ymax = 1.5;
ydir = 'normal';      

script_pretreatment_age
script_noisefilter_age
script_recurrenceplot_age
script_recurrencequantificationanalysis_age

axes1(25) = axes('Units','normalized',...
        'Position',[0.05 0.7 0.3 0.3],...
    'XLim',[agemodelmin agemodelmax],...
    'YLim',[agemodelmin agemodelmax],...
    'LineWidth',1.5,...
    'Box','On',...
    'Layer','top',...
    'XGrid','On'); hold on
axis square xy
imagesc(t,t,RR)
title(x_s)
colormap([1 1 1; 0 0 0])

y = [5  5    156  156 43  43 ];
x = [43 156  161  5   5   43 ];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    373  373 270 270];
x = [270 373  373  5   5   270];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

y = [5   5    540  540 440 440];
x = [440 540  540  5   5   440];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

line([0 150], [150 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([150 150], [0 150],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 

line([0 400], [400 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 
line([400 400], [0 400],'Color',"#A2142F","LineWidth",2,'LineStyle','--'); % 


sl = 350;
el = 450;
line([el el], [sl el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#7E2F8E","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#7E2F8E","LineWidth",2.5); % 


sl = 200;
el = 300;
line([el el], [sl el],'Color',"#D95319","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#D95319","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#D95319","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#D95319","LineWidth",2.5); % 

sl = 32;
el = 132;
line([el el], [sl el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl el], [el el],'Color',"#0072BD","LineWidth",2.5); % 
line([sl sl], [el-10 el+10],'Color',"#0072BD","LineWidth",2.5); % 
line([el-10 el+10], [sl sl],'Color',"#0072BD","LineWidth",2.5); % 

xlabel('Time (kyr)')
ylabel('Time (kyr)')



a(5) = axes('Position',[0.1 0.58 0.8 0.06],...
    'XTickLabels',[],...
    'XGrid','Off',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","left");

ylabel("ODP 967")
yt = r_win_e(1,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);



line(t,r_win_e(1,:),...
    'LineWidth',.5, ...
    "color","black")


a(6) = axes('Position',[0.1 0.39 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
     "YAxisLocation","left", ...
    'YDir', 'normal', ...
    'Color','None');

ylabel("ODP 967")


yt = r_win_e(2,:); % your time series y values
xt = t; % your time series x values
idx = find(~isnan(xt) & ~isnan(yt));

% Create new arrays without NaN values
xt = xt(idx);
yt = yt(idx);

mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);


line(t,r_win_e(2,:),...
    'LineWidth',.5, ...
    "color","black")

RQA_ODP =  r_win_e(1:2,:);


a(7) = axes('Position',[0.1 0.28 0.8 0.39],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","right", ...
    'YColor',"#000000", ...
    'Ylim', [0 1],...
    'Visible','Off');


l(9) = line(0,0);

ind2 = [43 156 270 373 440 535];

x = [ind2(1) ind2(2) ind2(2) ind2(1)];
y = [0 0 1 1];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

x = [ind2(3) ind2(4) ind2(4) ind2(3)];
patch(x,y,'black', 'FaceAlpha',0.05, 'EdgeColor','none')

x = [ind2(5) ind2(6) ind2(6) ind2(5)];
patch(x,y,'black', 'FaceAlpha',0.05, 'EdgeColor','none')

text(30,0.97,"E1",'HorizontalAlignment','Right','FontSize',12);
text(100,0.97,"E2",'HorizontalAlignment','Right','FontSize',12);
text(220,0.97,"E3",'HorizontalAlignment','Right','FontSize',12);
text(320,0.97,"E4",'HorizontalAlignment','Right','FontSize',12);
text(420,0.97,"E5",'HorizontalAlignment','Right','FontSize',12);
text(490,0.97,"E6",'HorizontalAlignment','Right','FontSize',12);

xline(400,'Color',"#A2142F","LineWidth",1.5,'LineStyle','--')
xline(150,'Color',"#A2142F","LineWidth",1.5,'LineStyle','--')

sl = 350;
el = 450;
hl = 0.45;

line([sl el], [hl hl],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl sl], [hl hl],'Color',"#7E2F8E","LineWidth",2.5); % 
line([el el], [hl-0.01 hl+0.01],'Color',"#7E2F8E","LineWidth",2.5); % 
line([sl sl], [hl-0.01 hl+0.01],'Color',"#7E2F8E","LineWidth",2.5); % 

sl = 200;
el = 300;

line([sl el], [hl hl],'Color',"#D95319","LineWidth",2.5); % 
line([sl sl], [hl hl],'Color',"#D95319","LineWidth",2.5); % 
line([el el], [hl-0.01 hl+0.01],'Color',"#D95319","LineWidth",2.5); % 
line([sl sl], [hl-0.01 hl+0.01],'Color',"#D95319","LineWidth",2.5); % 

sl = 32;
el = 132;

line([sl el], [hl hl],'Color',"#0072BD","LineWidth",2.5); % 
line([sl sl], [hl hl],'Color',"#0072BD","LineWidth",2.5); % 
line([el el], [hl-0.01 hl+0.01],'Color',"#0072BD","LineWidth",2.5); % 
line([sl sl], [hl-0.01 hl+0.01],'Color',"#0072BD","LineWidth",2.5); % 


a(8) = axes('Position',[0.1 0.67 0.8 0.0075],...
    'Visible','Off');

x = [0 14 14 0];
y = [0 0 1 1];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')
fsn = 8;
text(-15,0.5,"MIS","FontSize",fsn)
text(5,0.5,"1","FontSize",fsn)
text(20,0.5,"2","FontSize",fsn)
text(40,0.5,"3","FontSize",fsn)
text(62,0.5,"4","FontSize",fsn)
text(100,0.5,"5","FontSize",fsn)
text(162,0.5,"6","FontSize",fsn)
text(215,0.5,"7","FontSize",fsn)
text(270,0.5,"8","FontSize",fsn)
text(318,0.5,"9","FontSize",fsn)
text(350,0.5,"10","FontSize",fsn)
text(395,0.5,"11","FontSize",fsn)
text(450,0.5,"12","FontSize",fsn)
text(505,0.5,"13","FontSize",fsn)



x = [29 57 57 29];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [71 130 130 71];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [191 243 243 191];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [300 337 337 300];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [374 424 424 374];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [478 533 533 478];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')


linkaxes(a,"x")
a(1).XLim = [0 550];

%export figure as vector grafic
exportgraphics(figure1,'Figure 04.pdf','ContentType','vector')

%delete local path
rmpath(genpath(folder));


