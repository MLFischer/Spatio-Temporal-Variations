% Script to read the ODP 967 Data (Grant et al. 2017)
%
% 1. Reads data from file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |dataodp967age|:
%
%     1 = age (kyrs BP)
%     2 = PC1 sapropel
%     3 = PC2 sapropel
%     4 = dust proxy
%     5 = wetness index
%
% 3. Creates string array |dataodp967string| with the type of data.
%
% 21 Aug 2019 - Trauth

% Read data from file.
fid = fopen('data_odp_967_grant.txt');
formstrg = [repmat('%s ',1,4),repmat('%f ',1,23)];
A1 = textscan(fid,formstrg,'Delimiter','\t',...
    'Headerlines',9);

fclose(fid);


for i = 1:23
    datagrant(:,i) = A1{i+4};
end

% Merge into a single array.
dataodp967(:,1) = datagrant(:,22);       % Age
dataodp967(:,2) = datagrant(:,18);       % PC1 sapropel
dataodp967(:,3) = datagrant(:,19);       % PC2 detritus
dataodp967(:,4) = datagrant(:,21);       % Dust proxy
dataodp967(:,5) = datagrant(:,23);       % Wetness index

% Create string array for graphics.
dataodp967string = ["ODP 967 PC1 Detritus";
                    "ODP 967 PC2 Sapropel";
                    "ODP 967 Dust Proxy";
                    "ODP 967 Wetness Index"];             

% Remove rows without data.
dataodp967(isnan(dataodp967(:,1))==1,:) = [];

% Interpolating data to evenly spaced time axis.
dataodp967age(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:5
    dataodp967age(:,i) = interp1(dataodp967(:,1),...
        dataodp967(:,i),dataodp967age(:,1),inttype);
end
               
clear A1 ans datagrant fid formstrg i dataodp967
