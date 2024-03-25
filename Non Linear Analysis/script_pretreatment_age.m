%% Script to pretreat data
%
% Trauth, M.H., Asrat, A., Cohen, A., Duesing, W., Foerster, V.,
% Kaboth-Bahr, S., Kraemer, H.,  Lamb, H., Marwan, N., Maslin, M.,
% Schaebitz, F. (2021) Recurring types of variability and transitions in
% the ~620 kyr record of climate change from the Chew Bahir basin, southern
% Ethiopia, Quaternary Science Reviews.
%
% https://doi.org/10.1016/j.quascirev.2020.106777
%
% Required functions and external scripts:
% (none)
%
% Expected variables:
% - x (vector with proxy values)
% - t (vector with time values)
% - filteroption (scalar with flag whether to apply filter or not)
% - filterorder (scalar with order of the butterworth filter)
% - filtercutoff (scalar with the cutoff frequency of the high-pass filter)

% Delivers:
% - x (vector with filtered proxy values)
%
% Copyright (c) 2021
%
% Martin H. Trauth, University of Potsdam, Potsdam, Germany
% http://www.martinhtrauth.de
%
% Contact: trauth@uni-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

samplingint = abs(mean(diff(t)));
samplingfreq = 1/samplingint;
nyquistfreq = abs(0.5 * mean(diff(t))^(-1));
if filteroption == 1
    [b,a] = butter(filterorder,filtercutoff/nyquistfreq,'high');
    [hfilt,wfilt] = freqz(b,a,1024);
    f = samplingfreq*wfilt/(2*pi);
    x = filtfilt(b,a,x);
end

lenTS = length(x);
