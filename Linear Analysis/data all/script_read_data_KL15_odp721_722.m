% Script to read the ODP721 722 for KL15 XRF data analysis
%
% 11 Nov 2022 - by Fischer

% Read the tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
data_odp721_722_terr     = readtable('data_odp_721_722_terr.txt');
data_odp721_722_terr     = table2array(data_odp721_722_terr);

dataodp721_722_terr_string = ["ODP 721/722 terr"];

% Interpolating data to evenly spaced time axis.
dataodp721_722_terr_age(:,1) = agemodelmin : agemodelres : agemodelmax;


    dataodp721_722_terr_age(:,2) = interp1(data_odp721_722_terr(:,2),...
        data_odp721_722_terr(:,3),dataodp721_722_terr_age(:,1),inttype);
               
clear data_odp721_722_terr 