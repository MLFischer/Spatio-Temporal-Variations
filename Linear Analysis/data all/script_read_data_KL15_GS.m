% Script to read the KL15 SUS data
%
% 11 Nov 2022 - by Fischer

% Read the tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
datakl15gs      = readtable('data_Kl15_grain_sizes_DF_new.txt');

gs_string  = convertCharsToStrings(datakl15gs.Properties.VariableNames');
gs         = table2array(datakl15gs);
