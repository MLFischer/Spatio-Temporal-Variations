function y = embed(varargin)
%EMBED   Create embedding vector using time delay embedding
%     Y=EMBED(X,M,T) create the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,17);
%         plot(y(:,1),y(:,2))

% Copyright (c) 2012
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

narginchk(1,3)
nargoutchk(0,1)

try
    t = varargin{3};
catch ME
    t = 1;
end
        
try
    m = varargin{2};
catch ME
    m = 1;
end    

x = varargin{1}(:);

Nx = length(x);
NX = Nx-t*(m-1);

for mi = 1:m;
    jx(1+NX*(mi-1):NX+NX*(mi-1)) = 1+t*(mi-1):NX+t*(mi-1);
end

y = reshape(x(jx),NX,m);
