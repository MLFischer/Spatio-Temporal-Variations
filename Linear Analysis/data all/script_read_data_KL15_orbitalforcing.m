% Script to read the orbital parameters (Laskar et al. 2004)
%
% 1. Reads data from file.
% 2. Interpolates data to evenly spaced time axis to create a new double
%    array |dataorbitalage|:
%
%     1 = age (kyrs BP)
%     2 = eccentricity
%     3 = precession
%     4 = obliquity
%
% 3. Creates string array |dataorbitalstring| with the type of data.
%
% 22 Aug 2019 - Trauth

% Read data from file.
dataorbital = load('data_orbitalforcing.txt');
dataorbital(:,1) = -dataorbital(:,1);
dataorbital = sortrows(dataorbital,1);

% Create string array for graphics.
dataorbitalstring = ["Laskar 04 Eccentricity";
                     "Laskar 04 Precession";
                     "Laskar 04 Obliquity"];

% Interpolating data to evenly spaced time axis.
dataorbitalage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:4
    dataorbitalage(:,i) = interp1(dataorbital(:,1),...
        dataorbital(:,i),dataorbitalage(:,1),inttype);
end


datains         = readtable('Insolation.txt');
datains        = table2array(datains);
datains(:,1) = datains(:,1)./1000;
datainsstring   = [     "Insolation March Eq.";
                        "Insolation Autumn Eq.";
                        "Insolation Summer 30°N"];
     
datainsage(:,1) = agemodelmin : agemodelres : agemodelmax;
for i = 2:4
    datainsage(:,i) = interp1(datains(:,1),...
        datains(:,i),datainsage(:,1),inttype);
end

clear dataorbital i datains
