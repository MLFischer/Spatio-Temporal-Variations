function [y,P,epsilon] = rp(varargin)
% RP computes the recurrence matrix
%
%    [R,DM,EPSILON] = RP(Y,E,THRES_CALC,NORM,TYPE,ALGORITHM) 
%    
%    Calculates a recurrence matrix R from an embedding vector Y and using 
%    the scalar threshold E. Matrix DM is an optional output and is the
%    distancematrix. In case you choose the 'var'-fixed threshold selection
%    method, the optional output vector EPSILON will store the actual calculated
%    value.
%
%    Input:
%
%    Y                      is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
%
%    NORM (optional)        a string that specified the used norm for the 
%                           distance calculation in phasespace; available 
%                           values are 'euc' for euclidean norm or 'max' (default)
%                           for the maximum norm.
%
%    ALGORITHM (optional)   a string that selects different implementations
%                           of calculating the distance matrix; available 
%                           values are 'loops', 'vector' (default), or
%                           'matlabvector'
%
%    THRES_CALC (optional)  a string that specifies how the threshold epsilon
%                           will be applied:
%                           - 'fix' (default) The RP is computed under a fixed
%                             threshold epsilon specified by input value E.
%                           - 'var' The RP is computed under a fixed threshold
%                             epsilon, which corresponds to the lower E% quantile
%                             (specified by input parameter E) of the distance
%                             distribution of all points in phasespace.
%                           - 'fan' The RP is computed under a variable threshold 
%                             ensuring a fixed amount of nearest neighbours in 
%                             phasespace to compute the, corresponding to a
%                             threshold value for each point of the phasespace
%                             trajectory individually.
%
%    TYPE (optional)        a string that specifies the type of the RP:
%                           - 'normal' (default) The RP is computed using the  
%                             definition of Eckmann et al. 1987
%                           - 'diagonal' The RP is computed using the definition 
%                             of Eckmann et al. 1987 and then line corrected 
%                             corresponding to Kraemer and Marwan 2019
%                           - 'shape' The RP is computed using the definition 
%                             of Eckmann et al. 1987 and then shape-converted 
%                             corresponding Donath 2016 (windowshape 3).
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         imagesc(R)

% Copyright (c) 2019
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Modified by Hauke Kraemer,Potsdam Institute for Climate Impact Research, 
% Germany http://www.pik-potsdam.de
%
% Contact: hkraemer@pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%%
narginchk(1,6)
nargoutchk(0,3)

algoLib={'loops','vector','matlabvector'}; % the possible algorithms
try
    algorithm = varargin{6};
    if ~isa(algorithm,'char') || ~ismember(algorithm,algoLib)
        warning(['Specified algorithm should be one of the following possible values:',...
           10,sprintf('''%s'' ',algoLib{:})])
    end
catch
    algorithm = 'vector';
end

typeLib={'normal','diagonal','shape'}; % the possible types
try
    type = varargin{5};
    if ~isa(type,'char') || ~ismember(type,typeLib)
        warning(['Specified RP type should be one of the following possible values:',...
           10,sprintf('''%s'' ',typeLib{:})])
    end
catch
    type = 'normal';
end

methLib={'euc','max'};
try
    meth = varargin{4};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
catch
    meth = 'max';
end

thresLib={'fix','var','fan'};
try
    thres = varargin{3};
    if ~isa(thres,'char') || ~ismember(thres,thresLib)
       warning(['Specified way of calculating threshold should be one of the following possible values:',...
                                10,sprintf('''%s'' ',thresLib{:})])
    end
catch
    thres = 'fix';
end

try
    e = varargin{2};
catch
    e = 1;
end

x = varargin{1};

N = size(x);
if N(1) < N(2)
   warning('Embedding dimension is larger than the length of the vector. Please check!')
end

M=N(1);

switch algorithm 
    case 'loops'
        
         y = zeros(N(1),N(1));
         parfor i = 1:M
               for j = 1:M
                switch lower(meth)
                    case 'euc'
                         d = (x(i,:) - x(j,:)).^2;
                         y(i,j) = sqrt(sum(d));
                    otherwise
                         d = abs(x(i,:) - x(j,:));
                         y(i,j) = max(d);
                end
               end
         end
   case 'vector'
    
        x1 = repmat(x,N(1),1);
        x2 = reshape(repmat(reshape(x,N(1)*N(2),1),1,N(1))',N(1)*N(1),N(2));
        switch lower(meth)
          case 'euc'
              d = (x1 - x2).^2;
              y = sqrt(sum(d,2));
          otherwise
              d = abs(x1 - x2);
              y = max(d,[],2);
        end
        y = reshape(y,N(1), N(1));   
        
    case 'matlabvector'
        
      switch lower(meth)
          case 'euc'
              y = squareform(pdist(x));
              y_dist = squareform(pdist(x));
          otherwise
              y = squareform(pdist(x,'chebychev'));
      end
end
P = y;

if strcmp(thres,'fix')

    y = double(y < e);
    epsilon = e;
    
elseif strcmp(thres,'var')    

    epsilon = quantile(P(:),e);
    y = double(y < epsilon);        

elseif strcmp(thres,'fan')

    q = quantile(P,e);
    thresholds = repmat(q,N(1),1);

    epsilon = e;

    y=double(y<thresholds);
end

if strcmp(type,'diagonal')
    [y, ~] = rp_diagonal(y);
elseif strcmp(type,'shape')
    y = shape3(y);
end

end

function [Y] = shape3(X)

s = floor(size(X,1)/2);
sc = floor(size(X,1)/4);
Y = zeros(size(X,1));
for i = 1:s
    for j = 0:sc-1
        Y(sc-j+i,sc+j+i) = X(sc-j+i,sc+j+i);
        Y(sc-j+i,sc+j+i+1) = X(sc-j+i,sc+j+i+1);
        Y(sc+j+i,sc-j+i) = X(sc+j+i,sc-j+i);
        Y(sc+j+i,sc-j+i+1) = X(sc+j+i,sc-j+i+1);
    end
end

end

function [X_new,dl_new] = rp_diagonal(varargin)

narginchk(1,3)
nargoutchk(1,2)

X = varargin{1};

[N_org,M_org] = size(X);
if N_org~=M_org
    error('Input needs to be a squared, binary recurrence matrix')
end
if sum(sum(X>1))~=0 || sum(sum(X<0))~=0 || sum(sum(rem(X,1)))~=0
    error('Input needs to be a squared, binary recurrence matrix')
end

if issymmetric(X)
    symm = true;

    X2 = tril(X);
    
    X_cl = convertRP(X2);
 
    [~, lines_1] = dl_h(X_cl);
        
    lines_1_copy = lines_1;
    
    Nlines = size(lines_1,2);
    
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
      
    
else
    symm = false;

    X2 = tril(X);
    X3 = triu(X)';
    
    X_cl = convertRP(X2);
    X_cl2 = convertRP(X3);
    
    [~, lines_1] = dl_h(X_cl);
    [~, lines_2] = dl_h(X_cl2);
    
    lines_1_copy = lines_1;
    lines_2_copy = lines_2;
    
    Nlines = size(lines_1,2);
    Nlines2 = size(lines_2,2);
    
    X_hori = zeros(size(X_cl));
    for i = 1:size(lines_1,2)
        line_ind = lines_1(2,i);
        column_ind = lines_1(3,i);
        for j = 0:lines_1(1,i)-1
            X_hori(line_ind,column_ind+j) = lines_1(1,i);
        end
    end
    X_hori2 = zeros(size(X_cl2));
    for i = 1:size(lines_2,2)
        line_ind = lines_2(2,i);
        column_ind = lines_2(3,i);
        for j = 0:lines_2(1,i)-1
            X_hori2(line_ind,column_ind+j) = lines_2(1,i);
        end
    end
end

line_matrix_final = zeros(3,1);

[N,M] = size(X_hori);
for l_ind = 1:Nlines

    if ~ismember(lines_1(:,l_ind)',lines_1_copy','rows')
        continue
    end
    
    linei = lines_1(2,l_ind); 
    columni = lines_1(3,l_ind);
    
    line_matrix_final = horzcat(line_matrix_final,lines_1(:,l_ind));
    
    X_hori = delete_line_from_RP(X_hori,lines_1(:,l_ind));
 
    l_max = lines_1(1,l_ind);
    for l = 1:l_max
        
        for index = -1:2:1
            
            if linei+index > N | linei+index == 0
                break
            end
                           
            if X_hori(linei+index,columni+l-1)

                [X_hori,lines_1_copy]=scan_lines(X_hori,lines_1_copy,...
                    linei+index,columni+l-1);

            end

        end
        
    end

end

if ~symm

    line_matrix_final2 = zeros(3,1);
    for l_ind = 1:Nlines2

        if ~ismember(lines_2(:,l_ind)',lines_2_copy','rows')
            continue
        end

        linei = lines_2(2,l_ind); 
        columni = lines_2(3,l_ind);

        line_matrix_final2 = horzcat(line_matrix_final2,lines_2(:,l_ind));

        X_hori2 = delete_line_from_RP(X_hori2,lines_2(:,l_ind));

        l_max = lines_2(1,l_ind);
        for l = 1:l_max

            for scan = 1:2

                if scan == 1
                    index = 1;

                    if linei+index > N
                        break
                    end
                else
                    index = -1;

                    if linei+index == 0
                        break
                    end 
                end

                if X_hori2(linei+index,columni+l-1)

                    [X_hori2,lines_2_copy]=scan_lines(X_hori2,lines_2_copy,...
                        linei+index,columni+l-1);

                end

            end

        end      

    end
end

X_cl_new = zeros(N,M);
if ~symm
   X_cl2_new = zeros(N,M); 
end

for i = 1:size(line_matrix_final,2)
    l_max = line_matrix_final(1,i);
    linei = line_matrix_final(2,i);
    columni = line_matrix_final(3,i);
    for j = 1:l_max
        X_cl_new(linei,columni+j-1) = 1;
    end
end

if symm

    XX = revertRP(X_cl_new);    
    X_new = XX + (XX-eye(size(XX)))';
else 

    for i = 1:size(line_matrix_final2,2)
        l_max = line_matrix_final2(1,i);
        linei = line_matrix_final2(2,i);
        columni = line_matrix_final2(3,i);
        for j = 1:l_max
            X_cl2_new(linei,columni+j-1) = 1;
        end
    end

    XX = revertRP(X_cl_new);
    XXX= revertRP(X_cl2_new);
    X_new = XX + (XXX-eye(size(XXX)))';
    
end

if nargout > 1
    [~, lines3] = dl_e(X_new);
    dl_new = sortrows(lines3','descend')';
end

end

function [a_out, b_out] = dl_e(X)

[Y,~] = size(X);
lines(:,1) = getLinesOnDiag(X,-Y+1);
for j=-Y+2:Y-1
    lines = horzcat(lines,getLinesOnDiag(X,j)); 
end

zero_lines = lines(1,:)==0;
lines(:,zero_lines) = []; 

b_out= sortrows(lines','descend')';
a_out = mean(b_out(1,:));
end

function lines = getLinesOnDiag(M,j)

    d = diag(M,j);
    if ~any(d)
        lines = [0;0;0];
        return
    end
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);

    lines = zeros(3,numel(starts));
    for n=1:numel(starts)
        ll = get_indices(starts(n),j);
        for k = 1:length(ll)
            lll = ll{k};
            lines(2,n) = lll(1);
            lines(3,n) = lll(2);
            lines(1,n) = ends(n) - starts(n) +1;
        end
    end
    
end

function tuples = get_indices(indices,position_from_LOI)

narginchk(2,2)
nargoutchk(1,1)

if size(indices,1)<size(indices,2)
    indices = indices';
end
if rem(position_from_LOI,1)~=0
    error('position_from_LOI needs to be a integer')
end

tuples = cell(1,length(indices));

for i = 1:length(indices)
    
    if position_from_LOI < 0
        start_line = indices(i) + abs(position_from_LOI);
        start_column = indices(i);
    elseif position_from_LOI > 0
        start_line = indices(i);
        start_column = indices(i) + abs(position_from_LOI); 
    elseif position_from_LOI == 0
        start_line = indices(i);
        start_column = indices(i);  
    end    
   
    tuples{i}=[start_line start_column];
end


end

function [a_out, b_out]=dl_h(x)

narginchk(1,1)
nargoutchk(0,2)

[N,~] = size(x);

liness = zeros(3,1);
for j = 1:N
    d = x(j,:)';
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);
    
    if ~isempty(starts)
        lines = zeros(3,numel(starts));
        for n=1:numel(starts)
            lines(2,n) = j;
            lines(3,n) = starts(n);
            lines(1,n) = ends(n) - starts(n) +1;        
        end
    else
        lines = zeros(3,1);
    end
    liness = horzcat(liness,lines);
end

zero_lines = liness(1,:)==0;
liness(:,zero_lines) = []; 

b_out= sortrows(liness','descend')';
a_out = mean(b_out(1,:));
end

function Y = convertRP(X)

N = size(X);

Y = zeros(2*N(1)+1,N(1));

for i = 0:N(1)-1
   Y(N(1)+i+1,(1:(N(1)-i))) = diag(X,i);
end

for i = 0:N(1)-1
   Y(N(1)-i+1,(1:(N(1)-i))+i) = diag(X,-i);
end

end

function X = revertRP(Y)

N = size(Y);

X = zeros(N(2),N(2));

Z = [Y zeros(N(1),N(2)+1)];
Z = flipud(Z);

for i = 1:N(2)
    di = diag(Z,-i);
    X(:,N(2)-i+1) = di(1:N(2));
end
end

function RP = delete_line_from_RP(RP,l_vec)

    RP(l_vec(2),l_vec(3)+(1:l_vec(1))-1) = 0;

end

function [XX,YY]= scan_lines(XX,l_vec,line,column)
    
    index = 0;
    while true        

        loc_line = find(line==l_vec(2,:));
        del_ind = loc_line(column+index==l_vec(3,loc_line));
        
        if del_ind
            break
        else
            index = index - 1;
        end
    end   
    
    XX(l_vec(2,del_ind),l_vec(3,del_ind)+(1:l_vec(1,del_ind))-1) = 0;
    
    len = l_vec(1,del_ind);
    li = l_vec(2,del_ind);
    co = l_vec(3,del_ind);
    
    l_vec(:,del_ind) = [];
    
    [N,M] = size(XX);
    
    if li-1 < 1
        flag1 = false;
    else
        flag1 = true;
    end
    if li+1 > N
        flag2 = false;
    else
        flag2 = true;
    end
    
    for i = 1:len
        
        if li-1 < 1 || co+i-2 == 0
            flag1b = false;
        else
            flag1b = true;
        end
        if li-1 < 1 || co+i > M
            flag1c = false;
        else
            flag1c = true;
        end
        if li+1 > N || co+i-2 == 0
            flag2b = false;
        else
            flag2b = true;
        end
        if li+1 > N || co+i > M
            flag2c = false;
        else
            flag2c = true;
        end
        
        if flag1b && XX(li-1,co+i-2)

            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-2);
            
        elseif flag1 && XX(li-1,co+i-1)

            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i-1);
            
        elseif flag1c && XX(li-1,co+i)

            [XX,l_vec] = scan_lines(XX,l_vec,li-1,co+i);
            
        elseif flag2b && XX(li+1,co+i-2)

            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-2);
            
        elseif flag2 && XX(li+1,co+i-1)

            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i-1);
             
        elseif flag2c && XX(li+1,co+i)

            [XX,l_vec] = scan_lines(XX,l_vec,li+1,co+i);
        end
    end
    
    YY = l_vec;

end

function y = embed(varargin)

narginchk(1,3)
nargoutchk(0,1)

try
    t = varargin{3};
catch
        t = 1;
end
        
try
    m = varargin{2};
catch
        m = 1;
end    

x = varargin{1}(:);

Nx = length(x);
NX = Nx-t*(m-1);
if t*(m-1) > Nx
   warning('embedding timespan exceeding length of time series')
end

for mi = 1:m
    jx(1+NX*(mi-1):NX+NX*(mi-1)) = 1+t*(mi-1):NX+t*(mi-1);
end

y = reshape(x(jx),NX,m);
end
