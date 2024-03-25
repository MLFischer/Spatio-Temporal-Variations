% Script to read the KL15 XRF data
%
% 11 Nov 2022 - by Fischer

% Read the tables
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
datakl15xrf     = readtable('data_KL15_XRF.txt');
datakl15agem    = readtable('data_KL15-2_smooth_spline_ages.txt');
% Rename the ID colum, as this is the depth and the match
datakl15xrf     = renamevars(datakl15xrf, "User_ID","depth");
datakl15agem     = renamevars(datakl15agem, "best","age");

%read qualitify flags: 
datakl15qf      = readtable('data_KL15_qf.txt');
datakl15qf        = table2array(datakl15qf);




% Do outerjoin, save the variables and get back to array
T               = outerjoin(datakl15agem,datakl15xrf);
datakl15string  = convertCharsToStrings(T.Properties.VariableNames');
datakl15        = table2array(T);
% remove QF marked measurements:
datakl15(ismember( datakl15(:,1),datakl15qf(:,1)),:)= [];

% Remove rows without data.
datakl15(isnan(datakl15(:,6))==1,:) = [];

% a to ka
datakl15(:,4) = datakl15(:,4)/1000;


%datakl15string(31,1)
%datakl15string(21,1)
%log(datakl15(1:20,31)./datakl15(1:20,21))

% Interpolating data to evenly spaced time axis.
datakl15age(:,1) = agemodelmin : agemodelres : agemodelmax;

for i = 2:29
    datakl15age(:,i) = interp1(datakl15(:,4),...
        datakl15(:,i+5),datakl15age(:,1),inttype);
end
               
clear datakl15xrf datakl15agem datakl15xrf datakl15 i T
