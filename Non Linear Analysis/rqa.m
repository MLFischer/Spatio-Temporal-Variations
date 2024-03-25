function y = rqa(varargin)
% RQA performs recurrence quantification analysis
%
%    Y=RQA(X,L,T,LINE_CORRECT) calculates measures of recurrence 
%    quantification analysis from a recurrence matrix X using minimal line 
%    length L and a Theiler window T. Optional input LINE_CORRECT can
%    be used in order to correct for border effects while counting the
%    diagonal line length histograms (Kraemer & Marwan 2019).
%    LINE_CORRECT = 0 (default) specifies the conventional method,
%    LINE_CORRECT = 1 forces to use the 'kelo'-correction.
%
%    Output:
%      Y(1)  = RR     (recurrence rate)
%      Y(2)  = DET    (determinism)
%      Y(3)  = <L>    (mean diagonal line length)
%      Y(4)  = Lmax   (maximal diagonal line length)
%      Y(5)  = ENTR   (normalized entropy of the diagonal line lengths)
%      Y(6)  = LAM    (laminarity)
%      Y(7)  = TT     (trapping time)
%      Y(8)  = Vmax   (maximal vertical line length)
%      Y(9)  = RTmax  (maximal white vertical line length)
%      Y(10) = T2     (recurrence time of 2nd type)
%      Y(11) = RTE    (recurrence time entropy, i.e., RPDE)
%      Y(12) = Clust  (clustering coefficient)
%      Y(13) = Trans  (transitivity)
% 
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         Y = rqa(R);
%
% Copyright (c) 2016
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Modified by Hauke Krämer,Potsdam Institute for Climate Impact Research, 
% Germany http://www.pik-potsdam.de
%
% Contact: hkraemer@pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

narginchk(1,4)
nargoutchk(0,1)

try 
    line_correct = varargin{4};
    if line_correct ~= 1 && line_correct ~= 0
        warning('input "line_correct" must be set to 0 (conventional) or 1 (kelo correction). Now set to 0.')
        line_correct = 0;
    end
catch
    line_correct = 0;
end

try
    theiler_window = varargin{3};
catch
    theiler_window = 1;
end

try
    l_min = varargin{2};
catch
    l_min = 2;
end
x = varargin{1};

N = size(x);

y = zeros(13,1);

if theiler_window
   x_theiler = double(triu(x,theiler_window) + tril(x,-theiler_window));
else
   x_theiler = double(x);
end

N_all = N(1)*N(2);
N_all = N_all - N(1) - 2*((theiler_window-1)*N(1) - sum(1:(theiler_window-1)));

N_recpoints = sum(x_theiler(:));
y(1) = N_recpoints/N_all; 

if line_correct == 0
    [~,l_hist] = dl_conventional(x_theiler);
else
    [~,l_hist] = dl_kelo(x_theiler);
end

l_hist = hist(l_hist,1:N(1));

y(2) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

if isnan(sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1))))
    y(3) = 0;
else
    y(3) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1)));
end

if any(l_hist)
   y(4) = find(l_hist,1,'last');
end

l_hist = l_hist(l_min:end);

l_classes = sum(l_hist~=0); 

l_prob = l_hist/sum(l_hist);
ent_Sum = (l_prob .* log(l_prob));

if l_classes > 1
    y(5) = -nansum(ent_Sum)/log(length(l_min:N(1)));
else
    y(5) = -nansum(ent_Sum);
end

v_hist = zeros(1,N(1));
for i = 1:N(1)
   cnt = 0;
   for j = 1:N(2)
      if x_theiler(j,i)
         cnt = cnt+1;
      else
         if cnt
             v_hist(cnt) = v_hist(cnt) + 1;
         end
         cnt = 0;
      end
   end
   if cnt
       v_hist(cnt) = v_hist(cnt) + 1;
   end
end

y(6) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

if isnan(sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1))))
    y(7) = 0;
else
    y(7) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1)));
end

if any(v_hist)
   y(8) = find(v_hist,1,'last');
end

rt_hist = zeros(1,N(1));
for i = 1:N(1)
   cnt = 0;
   
   first_flag = false;
   
   for j = 1:N(2)
      if ~x(j,i)
         if first_flag
            cnt = cnt + 1;
         end
      else
         first_flag = true;
         if cnt
             rt_hist(cnt) = rt_hist(cnt) + 1;
         end
         cnt = 0;
      end
   end
end

if any(rt_hist)
    y(9) = find(rt_hist,1,'last');
end

y(10) = sum(rt_hist .* (1:N(1))) / sum(rt_hist);
if isnan(y(10))
    y(10)=0;
end

rt_classes = sum(rt_hist~=0);
rt_prob = rt_hist/sum(rt_hist);
ent_Sum = (rt_prob .* log(rt_prob));
if rt_classes > 1
    y(11) = -nansum(ent_Sum)/log(N(1));
else
    y(11) = -nansum(ent_Sum);
end

kv = sum(x_theiler,1);
y(12) = nanmean(diag(x_theiler*x_theiler*x_theiler)' ./ (kv .* (kv-1)));

denom = sum(sum(x_theiler * x_theiler));
y(13) = trace(x_theiler*x_theiler*x_theiler)/denom;

end

function [a_out, b_out]=dl_kelo(varargin)

X = varargin{1};
styleLib={'normal','semi'};
try
    type = varargin{2};
    if ~isa(type,'char') || ~ismember(type,styleLib)
        warning(['Specified RP type should be one of the following possible values:',...
           10,sprintf('''%s'' ',styleLib{:})])
    end
catch
    type = 'normal';
end

[Y,~] = size(X);
if issymmetric(X)
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type);
    for j=-Y+2:-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end

    lines = horzcat(lines,lines);

    lines = horzcat(lines,max(borderlines),max(borderlines));
else
    [lines(1),borderlines(1)] = getLinesOnDiag(X,-Y+1,type);
    for j=-Y+2:Y-1
        [ll,bl] = getLinesOnDiag(X,j,type);
        lines = horzcat(lines,ll);
        borderlines = horzcat(borderlines,bl);
    end
    borderlines = sort(borderlines,'descend');

    lines = horzcat(lines,borderlines(2:3));
end

zero_lines = lines(:)==0;
lines(zero_lines) = []; 

b_out= sort(lines,'descend')';
a_out = mean(b_out);
end

function [lines, borderline] = getLinesOnDiag(M,j,type)
    d = diag(M,j);
    border_line_length = length(d);
    if ~any(d)
        lines = 0;
        borderline = 0;
        return
    end
    starts = find(diff([0; d],1)==1);
    ends = find(diff([d; 0],1)==-1);

    lines = zeros(1,numel(starts));
    borderline = zeros(1,numel(starts));
    
    if strcmp(type,'normal')
        for n=1:numel(starts)
            if ends(n) - starts(n) + 1 < border_line_length
                lines(n) = ends(n) - starts(n) +1;
            elseif ends(n) - starts(n) + 1 == border_line_length
                borderline = ends(n) - starts(n) +1;
            end
        end
    elseif strcmp(type,'semi')
        for n=1:numel(starts)
            if ends(n) ~= border_line_length && starts(n) ~=1               
                lines(n) = ends(n) - starts(n) +1;
            else
                borderline(n) = ends(n) - starts(n) +1;
            end
        end    
    end
    
end

function [a_out, b_out]=dl_conventional(x)

narginchk(1,1)
nargoutchk(0,2)

warning off
if any(x(:))

  if min(size(x))>100000

    x2=uint8(x);
    N=size(x2);
    x3=zeros(2*N(2)+N(1),N(2));
    x3(N(2)+1:N(2)+N(1),1:N(2))=x2;
    N3=size(x3);
    
    i2=repmat(((1:1+N(2))+N(1)+N(2))',1,N(2));
    i4=i2+repmat((2*N(2)+N(1)+1)*[0:N(2)-1],size(i2,1),1);
    i4(:,end)=[];
    i4=reshape(i4,size(i4,1)*size(i4,2),1);
    x3(i4)=[];
    x3(end)=[];
    x2=(reshape(x3,N(1)+N(2),N(2)))';
  
    x2(end+1,:)=0;
    x=reshape(x2,size(x2,1)*size(x2,2),1);
    x2=x(2:end);x(end)=[];
    z0=find(x==0&x2==1);
    z1=find(x2==0&x==1);
  
  else
  
    N=size(x);

    x1=spdiags(double(x));
    z=reshape(x1,size(x1,1)*size(x1,2),1);
    z2(2:length(z)+1)=z;z2(1)=0;z2(end+1)=0;
    z=diff(z2);
    z0=find(z==1);
    z1=find(z==-1);
  
  end
  if length(z0)>length(z1), z0(end)=[]; end
  if length(z1)>length(z0), z1(end)=[]; end
  
  if isempty(z0), z0=0; end
  if isempty(z1), z1=0; end
  
  if z0(1)>z1(1)
    z0(2:end+1)=z0(1:end);z0(1)=0; 
    if length(z0)>length(z1) 
       z0(end)=[];
    end
  end

  l=sort(z1-z0);
  l1=l(find(l-1));
  
  if nargout==2
     b_out=zeros(length(l),1);
     b_out=l';
  end
  
  if nargout>0
     a_out=mean(l1);
  else
     mean(l1)
  end
  
else

  if nargout==2
     b_out=NaN;
  end

  if nargout>0
     a_out=NaN;
  else
     NaN
  end

end

warning on
end
