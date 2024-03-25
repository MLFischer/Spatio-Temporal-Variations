%% Script to compute the recurrence plot
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
% - embed.m
% - rp.m
%
% Expected variables:
% - x (vector with proxy values)
% - m (scalar with embedding dimension)
% - tau (scalar with embedding delay)
% - threshold_calculation (string with the recurrence criterion)
% - norm (string with used norm for distance calculation)
% - timespan_diff (scalar with the reduction of the time series due to embedding)
%
% Delivers:
% - RR (recurrence matrix)
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

if strcmp(threshold_calculation,'var')
   e = 0.08;
elseif strcmp(threshold_calculation,'fan')
   e = 0.17;
elseif strcmp(threshold_calculation,'fix')
   e = 0.12*range(x);
end

xVec = embed(x,m,tau);

[R,~,eps] = rp(xVec,e,threshold_calculation,norm);

RR = NaN(size(x,2),size(x,2));
RR(1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2),...
   1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2)) = R;
