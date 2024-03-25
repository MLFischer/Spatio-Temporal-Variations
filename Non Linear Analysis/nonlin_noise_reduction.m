function y = nonlin_noise_reduction(x,m,epsilon)
%NONLIN_NOISE_REDUCTION Nonlinear noise reduction.
%
%    Y=NONLIN_NOISE_REDUCTION(X,M,EPSILON) performs a simple nonlinear 
%    noise reduction of values in vector X, embedded in a phase space
%    with dimension M and using a local neighborhood of size EPSILON.
%    Note that by applying this filter, there will be lost M-1
%    data points. We therefore phase-shift each datapoint in the resulting
%    signal by (M;-1)/2. Details about the algorithm can be found
%    in Kantz & Schreiber 2004.
%    
%
% Copyright (c) 2021
%
% K.H. Kraemer, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Contact: hkraemer@pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

if ~isvector(x)
  error("Input data must be a univariate time series.")
end
if rem(m,1) ~= 0 || m < 2
  error("Parameter m must be a positive integer value larger that 1.")
end
if  rem(m,1) ~= 0 || epsilon <= 0
  error("Parameter epsilon must be a positive integer value larger than 0")
end
    
Y = embed(x,m,1);

filtered_signal = zeros(1,length(Y));

Distance_matrix = squareform(pdist(Y,'chebychev'));
[~, sort_indices] = sort(Distance_matrix,2);

for i = 1:length(Distance_matrix)
    neighborhoodsize_index = epsilon;
    dd = 0;
    for neighbors = 2:neighborhoodsize_index
        dd = dd + Y(sort_indices(i,neighbors),ceil(m/2));
    end
    filtered_signal(i) = dd/length(2:neighborhoodsize_index);
end

y = NaN*ones(size(x));
y(1+floor((m-1)/2):end-ceil((m-1)/2)) = filtered_signal;
