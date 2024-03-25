% KL 15 Analysis Script
% Figure 03 time series and spearman cross correlation analysis
% Markus L Fischer
% 2023/09

% Importing KL15 data for TSA
clear, clc, close all

% get subfolders
folder = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(folder));

% settings for data interpolation
agemodelmin = 5; 
agemodelres = 0.1;
agemodelmax = 540; %540
inttype     = 'pchip';      

script_read_data_KL15_XRF
script_read_data_KL15_orbitalforcing
script_read_data_KL15_odp967
script_read_data_KL15_odp722
script_read_data_KL15_odp721_722
script_readco2_age
script_read_data_KL15_CHB_K_RRM_MHT500

% create data array
data                        =   datakl15age;                                % all XRF   1:29
data(:,size(data,2)+1)      =   ones(1, length(data(:,1)));
data(:,size(data,2)+1)      =   repelem(1,length(data(:,1)));               % add Ones  31 
data(:,size(data,2)+1 :size(data,2)+3)      =   dataorbitalage(:,[2:4]);    % Ecc Prec Obl.
data(:,size(data,2)+1 :size(data,2)+4)      =   dataodp967age(:,[2:5]);     % ODP967
data(:,size(data,2)+1)      =   dataodp722sstage(:,[2]);                    % ODP722 SST
data(:,size(data,2)+1)      =   dataodp721_722_terr_age(:,2);               % ODP 721 722 terr. influx
data(:,size(data,2)+1)      =   dataepicaco2age(:,2);                       % EPICA CO2
data(:,size(data,2)+1 :size(data,2)+13)      =   ones(length(data(:,1)), 13);
data(:,size(data,2)+1)      =   ones(1, length(data(:,1)));                 % CHB v2
data(:,size(data,2)+1 :size(data,2)+3)      =   datainsage(:,[2:4]);        % Ecc Prec Obl.

% create name array
datastr                         =   datakl15string([4 7:34],:);       % create a name string of the data
xrfindex                        =   find(contains(datastr,'_Area'));  % get index numbers of XRF data 
datastr                         =   strrep(datastr,'_Area','');       % replace underscore with normal spacing

datastr(size(datastr,1)+1,:)    =   'SUS';
datastr(size(datastr,1)+1,:)    =   'ones';
datastr(size(datastr,1)+1:size(datastr,1)+3,:)    =   dataorbitalstring;
datastr(size(datastr,1)+1:size(datastr,1)+4,:)    =   dataodp967string;
datastr(size(datastr,1)+1,:)    =   dataopd722sststring;
datastr(size(datastr,1)+1,:)    =   dataodp721_722_terr_string;
datastr(size(datastr,1)+1,:)    =   dataepicaco2string;
datastr(size(datastr,1)+1:size(datastr,1)+13,:)    =   repelem('ones',13);
datastr(size(datastr,1)+1,:)    =   'ones';
datastr(size(datastr,1)+1:size(datastr,1)+3,:)    =   datainsstring;

% create an index guide:
indexguide      = table(datastr,(1:size(datastr,1))');
xrfindex_chb    =   find(contains(datastr,'CHB'));  % get index numbers of XRF data 

%read grainsize
datakl15gs      = readtable('data_Kl15_grain_sizes_DF_new.txt');
gs_string  = convertCharsToStrings(datakl15gs.Properties.VariableNames');
gs         = table2array(datakl15gs);

%read MIS
dataLR04      = readtable('data_LR04stack.txt');
dataLR04         = table2array(dataLR04);
                
%
data(:,size(data,2)+1)          = log(data(:,26)./data(:,16));
datastr(size(datastr,1)+1,:)    = "KL15 log(Ti/Al)";

data(:,size(data,2)+1)          = (-1).*data(:,55);
datastr(size(datastr,1)+1,:)    = "CHB v2 RRMay2019+MHT500 K(inverted)";
indexguide      = table(datastr,(1:size(datastr,1))');
% Show all data:
varselectnum        = [31 31 60 59 26 16 24 2];       % Numerator
varselectdem        = [31 31 31 31 31 31 31 31];   % Denominator

O = data(:,38);
C = rescale(data_chb_age_cb0103_rr_560_rr(:,2));
K = rescale(data(:,59));


for j = 1:4852

C_OC(j) = corr( O(j:j+499) , C(j:j+499), 'Type', 'Spearman');
C_KO(j) = corr( K(j:j+499) , O(j:j+499), 'Type', 'Spearman');
C_KC(j) = corr( K(j:j+499) , C(j:j+499), 'Type', 'Spearman');

end 


ind2 = [43 156 270 373 440 535];

%run displayscript
script_displayresults_KL15_Figure02_v2

%%
%export figure as vector grafic
exportgraphics(fig03,'Figure 03.pdf','ContentType','vector')

%delete local path
rmpath(genpath(folder));