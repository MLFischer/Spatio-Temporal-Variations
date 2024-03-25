%% Script to perform the recurrence quantification analysis
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
% - rqa.m
%
% Expected variables:
% - R (recurrence matrix)
% - w (scalar specifying the window size)
% - ws (scalar that sets the window step width)
% - l_min (scalar of minimum line length)
% - theiler (scalar with size of the Theiler correction window)
% - line_correct (scalar with a flag whether to set line correction or not)
% - timespan_diff (scalar with the reduction of the time series due to embedding)
%
% Delivers:
% - r_win_e (vector with the results of RQA)
%
% Copyright (c) 2021
%
% K.H. Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Contact: hkraemer@pik-potsdam.de
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

r_win=zeros(13,ceil((length(R)-w)/ws));
cnt = 1;
for i = 1:ws:(length(R)-w)
     r_win(:,cnt) = rqa(R(i:(i+w),i:(i+w)),l_min,theiler,line_correct);
     cnt = cnt+1;
end

r_win_w(1:size(r_win,1),1:size(R,1)) = NaN;
r_win_w(:,1+w/2:ws:size(R,1)-w/2) = r_win;
for i = 1 : 13
     r_win_w(i,1+w/2:size(R,1)-w/2) = ...
     fillmissing(r_win_w(i,1+w/2:size(R,1)-w/2),'linear');
end

r_win_e(1:size(r_win,1),1:size(x,2)) = NaN;
r_win_e(:,1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2)) = r_win_w;

for i = 1 : 13
   r_win_e(i,isnan(x)==1) = NaN;
end

index_vector_first = [];
index_vector_last = [];
nan_vector = isnan(x);
flag = true;
for i = 1 : size(r_win_e,2)
   if i == 1
     if nan_vector(i) == 0
        flag = false;
     end
   end
   if nan_vector(i) == 0 && flag == true && i == 1
     flag = false;
   elseif nan_vector(i) == 1 && flag == false 
     index_vector_first = [index_vector_first i];
     flag = true;     
   elseif nan_vector(i) == 0 && flag == true 
     index_vector_last = [index_vector_last i-1];
     flag = false;
   end
end

for i = 1 : length(index_vector_first)
     j = index_vector_first(i);
     span = j - (w) - 2;
     gap = (w/2) + mod(span,ws);
     begin = j - gap;    
     for k = 1 : 13
       r_win_e(k,begin:j) = NaN;
     end
end

for i = 1 : length(index_vector_last)
     j = index_vector_last(i);
     span = j + (w);
     gap = (w/2) + mod(span,ws);
     ende = j + gap;    
     for k = 1 : 13
        r_win_e(k,j:ende) = NaN;
     end
end

  
