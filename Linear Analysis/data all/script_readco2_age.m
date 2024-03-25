% Script to read the Antarctica CO2 Data (Bereiter et al. 2015)
%
% 1. Reads data from file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |dataepicaco2age|:
%
%     1 = age (kyrs BP)
%     2 = CO2 (ppm)
%
% 3. Creates string array |dataepicaco2string| with the type of data.
%
% Important:
% There are two different age models in the file, EDC3,N and AICC2012,N
% I use the first one, EDC3,N.
%
% 21 Aug 2019 - Trauth

% Read data from file.
fid = fopen('data_antarctica2015co2.txt');
A1 = textscan(fid,'%f %f %f','Delimiter','\t','Headerlines',15);
fclose(fid);
for i = 1:3
    domccco2(:,i) = A1{i};
end
dataepicaco2(:,1) = domccco2(:,1)/1000;
dataepicaco2(:,2) = domccco2(:,2);

% Sort data.
dataepicaco2 = sortrows(dataepicaco2,1);

% Removing duplicate data points.
dataepicaco22 = dataepicaco2;
for i = 1:size(dataepicaco22,1)-1
    if dataepicaco22(i,1) == dataepicaco22(i+1,1)
        dataepicaco22(i,1) = NaN;
    end
end
dataepicaco22(isnan(dataepicaco22(:,1))==1,:) = [];
dataepicaco2 = dataepicaco22;
clear dataepicaco22

% Replacing NaNs by zeros.
dataepicaco2(isnan(dataepicaco2)==1) = 0;

% Remove rows without data.
dataepicaco2(isnan(dataepicaco2(:,1))==1,:) = [];

% Interpolating data to evenly spaced time axis.
dataepicaco2age(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2
    dataepicaco2age(:,i) = interp1(dataepicaco2(:,1),...
        dataepicaco2(:,i),dataepicaco2age(:,1),inttype);
end
               
% Create string array for graphics.
dataepicaco2string = "EPICA CO2";    

clear A1 ans fid formstrg i dataepicaco2 domccco2




