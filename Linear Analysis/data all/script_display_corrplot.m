% Produce the input lower triangular matrix data
% https://de.mathworks.com/matlabcentral/answers/699755-fancy-correlation-plots-in-matlab
% adapted by Markus L Fischer
% 2022/11

if symask == 1
    C = tril(C,-1);
end

C(logical(eye(size(C)))) = 1;

% Set [min,max] value of C to scale colors
clrLim = [-1,1];
 load('CorrColormap.mat') % Uncomment for custom CorrColormap
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];
myLabel = C_lab;
% Compute center of each circle
% This assumes the x and y values were not entered in imagesc()
x = 1 : 1 : size(C,2); % x edges
y = 1 : 1 : size(C,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0)=nan; % eliminate cordinates for zero correlations
% Set color of each rectangle
% Set color scale

cmap = CorrColormap; % Uncomment for CorrColormap
%cmap= customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
Cscaled = (C - clrLim(1))/range(clrLim); % always [0:1]
colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size between [0 1]
Cscaled = (abs(C) - 0)/1;
diamSize = Cscaled * range(diamLim) + diamLim(1);
% Create figure

%fh = figure('Position',[00 500 600 600])

colormap(ax(1),cmap);
% colormap(CorrColormap) %Uncomment for CorrColormap
tickvalues = 1:length(C);
x = zeros(size(tickvalues));
text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right');
x(:) = length(C)+1;
text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90);
% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax(1),'equal')
axis(ax(1),'tight')
set(ax(1),'YDir','Reverse')
colorbar()
caxis(clrLim);
axis off
hold(ax(1), 'off')

clear Cscaled cmap diamSize tickvalues h ax Cscaled Cscaled y myLabel diamLim theta x xAll yAll clrLim colIdx CorrColormap fh