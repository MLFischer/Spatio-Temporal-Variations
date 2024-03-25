% Script to read the ODP 722 SST Data (Herbert et al. 2010)
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
fid = fopen('data_odp_722_sst.txt');
A1 = textscan(fid,'%f %f %f %f','Delimiter','\t','Headerlines',17);
fclose(fid);
for i = 1:4
    dataherbert(:,i) = A1{i};
end
dataodp722sst(:,1) = dataherbert(:,2);
dataodp722sst(:,2) = dataherbert(:,4);

% Sort data.
dataodp722sst = sortrows(dataodp722sst,1);

% Removing duplicate data points.
dataodp722sst2 = dataodp722sst;
for i = 1:size(dataodp722sst2,1)-1
    if dataodp722sst2(i,1) == dataodp722sst2(i+1,1)
        dataodp722sst2(i,1) = NaN;
    end
end
dataodp722sst2(isnan(dataodp722sst2(:,1))==1,:) = [];
dataodp722sst = dataodp722sst2;
clear dataepicaco22

% Remove rows without data.
dataodp722sst(isnan(dataodp722sst(:,1))==1,:) = [];

% Replacing NaNs by zeros.
dataodp722sst(isnan(dataodp722sst)==1) = 0;

% Interpolating data to evenly spaced time axis.
dataodp722sstage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2
    dataodp722sstage(:,i) = interp1(dataodp722sst(:,1),...
        dataodp722sst(:,i),dataodp722sstage(:,1),inttype);
end
               
% Create string array for graphics.
dataopd722sststring = "ODP 722 SST";    

clear A1 ans fid formstrg i dataodp722sst dataodp722sst2 dataherbert




