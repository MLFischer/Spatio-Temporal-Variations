% KL 15 Analysis Script
% Figure 02 KL 15 XRF dataset exploration
% Markus L Fischer
% 2023/09

clear, clc, close all

% get subfolders
folder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(folder));

% interpolation settings
agemodelmin = 5; 
agemodelres = 0.1;
agemodelmax = 540; %540
inttype     = 'pchip';      

% read KL 15 XRF data
script_read_data_KL15_XRF

% create data array
data                        =   datakl15age;                                % all XRF   1:29

% create name array
datastr                         =   datakl15string([4 7:34],:);       % create a name string of the data
xrfindex                        =   find(contains(datastr,'_Area'));  % get index numbers of XRF data 
datastr                         =   strrep(datastr,'_Area','');       % replace underscore with normal spacing

% create an index guide:
indexguide      = table(datastr,(1:size(datastr,1))');
xrfindex = [4 6 8 16 18 22 26 28];

% get non interpolated xrf data and normalize them:
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
datakl15xrf2     = readtable('data_KL15_XRF.txt');
datakl15agem2    = readtable('data_KL15-2_smooth_spline_ages.txt');
% Rename the ID colum, as this is the depth and the match
datakl15xrf2     = renamevars(datakl15xrf2, "User_ID","depth");
datakl15agem2     = renamevars(datakl15agem2, "best","age");

%read qualitify flags: 
datakl15qf      = readtable('data_KL15_qf.txt');
datakl15qf        = table2array(datakl15qf);

% Do outerjoin, save the variables and get back to array
T               = outerjoin(datakl15agem2,datakl15xrf2);
datakl15string  = convertCharsToStrings(T.Properties.VariableNames');
datakl15        = table2array(T);
datakl15(isnan(datakl15(:,6))==1,:) = []; % deleta NA data
datakl15(:,4) = datakl15(:,4)/1000;       % a to ka

%subset of the variables of interest
data_r     = datakl15(:,[7 9 11 13 15 17 19 21 23 25 27 29 31 33]);
data_r_str = datakl15string([7 9 11 13 15 17 19 21 23 25 27 29 31 33]);
data_r_str =   strrep(data_r_str,'_Area','');       % replace underscore with normal spacing

%second subset of the 10 elements we want to use
data_r     = data_r(:,  [1 2 3 4 8 9 11 12 13 14]);
data_r_str = data_r_str([1 2 3 4 8 9 11 12 13 14]);

data_n = data_r;
data_n = (data_n - mean(data_n))./std(data_n);

% Euclidean distances
Y = pdist(data_n');

fig02 = figure('Position',[0 00 800 800],...
    'Color',[1 1 1]);

% Heatmap of Euclidean distances
ax(1) = axes('Position',[0.51 0.65 0.32 0.32]);
heatmap(flipud(squareform(Y)),...
    'Title','Euclidean Distance',...
    'XDisplayLabels',data_r_str,...
    'YDisplayLabels',flipud(data_r_str))

% Dendrogram plot 
ax(2) = axes('Position',[0.05 0.65 0.4 0.3],'XColor','none');

Z = linkage(Y);
% Dendrogram
[h,nodes,orig] = dendrogram(Z);
ylabel('Euclidian Distance')
set(gca,'XTickLabel',data_r_str(orig))


% Br and LR04
dataLR04      = readtable('data_LR04stack.txt');
dataLR04      = table2array(dataLR04);

script_readco2_age


ax(3) = axes('Position',[0.05 0.475 0.85 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'YDir', 'normal',...
    "YAxisLocation","left");

logBr_Ti = data(:,2);
ylim([0 4000])
l(1) = line(data(:,1),logBr_Ti,"Color","#000000","LineWidth",1);

yyaxis right
ylim([150 350])
l(2) = line(dataepicaco2age(1:5351,1),dataepicaco2age(1:5351,2),"Color","#0072BD","LineWidth",1);
legend(["Br (counts)" "Epica CO_2 (ppm)"])

% Plo9t log(Ca/Ti) and coarse silt
ax(4) = axes('Position',[0.05 0.35 0.85 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'YDir', 'normal',...
    "YAxisLocation","left");
logCa_Ti = log(data(:,24)./data(:,26));

l(1) = line(data(:,1),logCa_Ti,"Color","#000000","LineWidth",1);

yyaxis right
set(gca, 'YDir', 'reverse')

%read Grainsize
datakl15gs      = readtable('data_Kl15_grain_sizes_DF_new.txt');
gs_string  = convertCharsToStrings(datakl15gs.Properties.VariableNames');
gs         = table2array(datakl15gs);

gs2 = readtable('Kl15_grain_sizes_DF_new.csv','NumHeaderLines',1);
gs2 = table2array(gs2);

gs_sel = 12;
line(gs(:,2),gs2(:,gs_sel),"Color","#0072BD","LineWidth",1);


legend(["log(Ca/Ti)" "Coarse Silt"])


% Plot log(Ca/Ti) and fine silt
ax(5) = axes('Position',[0.05 0.225 0.85 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'YDir', 'normal',...
    "YAxisLocation","left");

l(1) = line(data(:,1),logCa_Ti,"Color","#000000","LineWidth",1);

yyaxis right

gs_sel = 10;
line(gs(:,2),gs2(:,gs_sel),"Color","#0072BD","LineWidth",1);
legend(["log(Ca/Ti)" "Fine Silt"])


% Plot log(Ti/Al) and median grain size
ax(6) = axes('Position',[0.05 0.1 0.85 0.1],...
    'XGrid','On',...
    'YDir', 'normal',...
    "YAxisLocation","left");
xlabel('Age (ka)')

log_Ti_Al = log(data(:,26)./data(:,16));
set(gca, 'YDir', 'reverse')
l(1) = line(data(:,1),log_Ti_Al,"Color","#000000","LineWidth",1);
yyaxis right
set(gca, 'YDir', 'reverse')
ylim([0 40])

line(gs(:,2),gs2(:,15),"Color","#0072BD","LineWidth",1);
legend(["log(Ti/Al)" "Grain Size Md."])

%export figure as vector grafic
exportgraphics(fig02,'Figure 02.pdf','ContentType','vector')

%delete local path
rmpath(genpath(folder));

