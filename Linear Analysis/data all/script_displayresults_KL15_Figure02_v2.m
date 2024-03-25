fig03 = figure('Position',[0 00 800 1200],...
    'Color',[1 1 1]);


b(1) = axes('Position',[0.1 0.95 0.8 0.006],...
    'Visible','Off');

x = [0 14 14 0];
y = [0 0 1 1];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')
fsn = 8;
text(-15,0.5,"MIS","FontSize",fsn)
text(5,0.5,"1","FontSize",fsn)
text(20,0.5,"2","FontSize",fsn)
text(40,0.5,"3","FontSize",fsn)
text(62,0.5,"4","FontSize",fsn)
text(100,0.5,"5","FontSize",fsn)
text(162,0.5,"6","FontSize",fsn)
text(215,0.5,"7","FontSize",fsn)
text(270,0.5,"8","FontSize",fsn)
text(318,0.5,"9","FontSize",fsn)
text(350,0.5,"10","FontSize",fsn)
text(395,0.5,"11","FontSize",fsn)
text(450,0.5,"12","FontSize",fsn)
text(505,0.5,"13","FontSize",fsn)


x = [29 57 57 29];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [71 130 130 71];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [191 243 243 191];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [300 337 337 300];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [374 424 424 374];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

x = [478 533 533 478];
patch(x,y,'black', ...
    'FaceAlpha',0.15, 'EdgeColor','none')

a(1) = axes('Position',[0.1 0.90 0.8 0.05],...
    'Visible','Off');
t(1) = text(0.99,1.0,"LR04",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');

b(2) = axes('Position',[0.1 0.90 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse',...
    "YAxisLocation","right");

l(1) = line(dataLR04(1:540,1),dataLR04(1:540,2),"Color","#000000","LineWidth",0.5);



yt = dataLR04(1:540,2); % your time series y values
xt = dataLR04(1:540,1); % your time series x values
mean_yt = min(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt<mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#0072BD', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.1, 'EdgeAlpha', 0);





a(12) = axes('Position',[0.1 0.87 0.8 0.04],...
    'Visible','Off');
t(12) = text(0.99,1,"Eccentricity",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');

b(3) = axes('Position',[0.1 0.87 0.8 0.04],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse');

line(data(:,1),-data(:,32),"LineStyle","-","Color","#000000");


line(data(-50+10*ind2,1),-data(-50+10*ind2,32),"LineStyle","none", ...
    "Color","#000000","Marker","diamond");

text(30,0.005,"E1",'HorizontalAlignment','Right');
text(100,0.005,"E2",'HorizontalAlignment','Right');
text(220,0.005,"E3",'HorizontalAlignment','Right');
text(320,0.005,"E4",'HorizontalAlignment','Right');
text(420,0.005,"E5",'HorizontalAlignment','Right');
text(490,0.005,"E6",'HorizontalAlignment','Right');

a(2) = axes('Position',[0.1 0.81 0.8 0.04],...
    'Visible','Off');


t(2) = text(0.99,1.15,"Precession and Insolation NHS 30°N",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(4) = axes('Position',[0.1 0.81 0.8 0.04],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse', ...
    "YAxisLocation","right");

line(data(:,1),data(:,33),"LineStyle","-","Color","#000000");


text(data(50,1),data(50,33)-0.02,"E1.1","FontSize",8)
text(data(280,1),data(280,33)-0.02,"E1.1","FontSize",8)
text(data(520,1),data(520,33)-0.02,"E2.1","FontSize",8)
text(data(770,1),data(770,33)-0.02,"E2.2","FontSize",8)
text(data(990,1),data(990,33)-0.02,"E2.3","FontSize",8)
text(data(1220,1),data(1220,33)-0.02,"E2.4","FontSize",8)
text(data(1450,1),data(1450,33)-0.02,"E2.5","FontSize",8)


b(5) = axes('Position',[0.1 0.78 0.8 0.04],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'normal', ...
    "YAxisLocation","left");

line(data(:,1),data(:,58),"LineStyle","-","Color","#000000");
ylabel("W/m^2")



a(3) = axes('Position',[0.1 0.67 0.8 0.1],...
    'Visible','Off');

t(3) = text(0.99,1,"ODP 967 Wetness Index",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(6) = axes('Position',[0.1 0.67 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'normal', ...
    'YColor',"#000000");

yt = data(:,38); % your time series y values
xt = data(:,1); % your time series x values
mean_yt = 1;
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt<mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#0072BD', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt>mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);



l(3) = line(data(:,1),data(:,38),"Color","#000000","LineWidth",0.5);




a(4) = axes('Position',[0.1 0.57 0.8 0.1],...
    'Visible','Off');
t(4) = text(0.99,1,"KL15 Aridity Index",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(7) = axes('Position',[0.1 0.57 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse', ...
    "YAxisLocation","right", ...
    'YColor',"#000000");

ylabel("rescaled log(Ti/Al)")


yt = K; % your time series y values
xt = data(:,1); % your time series x values
mean_yt = mean(yt) -0.05;
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#0072BD', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

l(4) = line(data(:,1),K,"Color","#000000","LineWidth",0.5);



a(5) = axes('Position',[0.1 0.45 0.8 0.1],...
    'Visible','Off');

t(5) = text(0.99,1.2,"Chew Bahir Aridity Index",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(8) = axes('Position',[0.1 0.45 0.8 0.1],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    'YDir', 'reverse', ...
    "YAxisLocation","left",...
    'YColor',"#000000");

ylabel("rescaled K counts")
ylim([0 1])

yt = C; % your time series y values
xt = data(:,1); % your time series x values
mean_yt = mean(yt)-0.075;
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#0072BD', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);


x = [108 115.4 115.4 108]; % x-coordinates of the vertices of the complex area
y = [4 4 0 0]; % y-coordinates of the vertices of the complex area
fill(x, y, 'white','FaceAlpha',1, 'EdgeAlpha', 0); % shade the complex area in red

x = [116.5 121.75 121.75 116.5]; % x-coordinates of the vertices of the complex area
y = [4 4 0 0]; % y-coordinates of the vertices of the complex area
fill(x, y, 'white','FaceAlpha',1, 'EdgeAlpha', 0); % shade the complex area in red


l(5) = line(data(:,1),C,"Color","#000000","LineWidth",0.5);



a(6) = axes('Position',[0.1 0.39 0.8 0.05],...
    'Visible','Off');

t(6) = text(0.99,0.96,"Insolation 0°N",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(9) = axes('Position',[0.1 0.39 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","right", ...
    'YColor',"#000000");

ylabel("W/m^2")


l(6) = line(data(:,1),data(:,56),"Color","#000000","LineWidth",0.5);
text(data(4090,1),500,"Spring")
line([data(4090,1) data(4090,1)],[data(4090,56)+10 485],'LineStyle','-',"Color","black")

l(6) = line(data(:,1),data(:,57),"Color","#000000","LineWidth",0.5);
text(data(4010,1),500,"Autumn",'HorizontalAlignment','Right')
line([data(4000,1) data(4000,1)],[data(4000,57)+10 485],'LineStyle','-',"Color","black")


a(7) = axes('Position',[0.1 0.33 0.8 0.05],...
    'Visible','Off');
t(7) = text(0.99,1.1,"ODP 722 SST",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(10) = axes('Position',[0.1 0.33 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","left", ...
    'YColor',"#000000");

ylabel("°C")


l(7) = line(data(20:5351,1),data(20:5351,39),"Color","#000000","LineWidth",0.5);



a(8) = axes('Position',[0.1 0.28 0.8 0.05],...
    'Visible','Off');

t(8) = text(0.99,1,"Epica CO_2",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(11) = axes('Position',[0.1 0.28 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
        'YDir', 'normal', ...
    "YAxisLocation","right", ...
    'YColor',"#000000");
ylabel("ppm")


l(8) = line(dataepicaco2age(:,1),dataepicaco2age(:,2),"Color","#000000","LineWidth",0.5);





a(9) = axes('Position',[0.1 0.21 0.8 0.05],...
    'Visible','Off');
t(9) = text(0.99,1.3,"ODP 967 ⋆ Chew Bahir",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(12) = axes('Position',[0.1 0.21 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
        'YDir', 'reverse', ...
    "YAxisLocation","left", ...
    'YColor',"#000000");
ylabel("Xcorr")


tsub = data(250:5101,1);


yt = C_OC'; % your time series y values
xt = data(250:5101,1); % your time series x values
mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

l(9) = line(data(250:5101,1),C_OC,"Color","#000000","LineWidth",0.5);





a(10) = axes('Position',[0.1 0.16 0.8 0.05],...
    'Visible','Off');

t(10) = text(0.99,1.1,"KL15 ⋆ ODP 967",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(13) = axes('Position',[0.1 0.16 0.8 0.05],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
        'YDir', 'reverse', ...
    "YAxisLocation","right", ...
    'YColor',"#000000");
ylabel("Xcorr")

yt = C_KO'; % your time series y values
xt = data(250:5101,1); % your time series x values
mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

l(10) = line(data(250:5101,1),C_KO,"Color","#000000","LineWidth",0.5);




a(13) = axes('Position',[0.1 0.085 0.8 0.05],...
    'Visible','Off');

t(13) = text(0.99,1.3,"KL15 ⋆ Chew Bahir",...
    'HorizontalAlignment','Right',...
    'VerticalAlignment','Top');
b(14) = axes('Position',[0.1 0.085 0.8 0.05],...
    'XGrid','Off',...
     "YAxisLocation","left", ...
    'YDir', 'normal', ...
    'Color','None', ...
    'Ylim',[-1 0.5]);
ylabel("Xcorr")


yt = C_KC'; % your time series y values
xt = data(250:5101,1); % your time series x values
mean_yt = mean(yt);
hold on
fill_x = [xt, fliplr(xt)];


fill_y = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y(yt>mean_yt) = mean_yt;
fill_y(1,1) = mean_yt;
fill_y(length(fill_y(:,1)),1) = mean_yt;
color_rgb = sscanf('#D95319', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y,color_rgb,'FaceAlpha', 0.25, 'EdgeAlpha', 0);

fill_y2 = [yt, fliplr(repmat(mean_yt, size(yt)))];
fill_y2(yt<mean_yt) = mean_yt;
fill_y2(1,1) = mean_yt;
fill_y2(length(fill_y2(:,1)),1) = mean_yt;
color_rgb = sscanf('#77AC30', '#%2x%2x%2x', [1 3])/255;
fill(fill_x, fill_y2, color_rgb, 'FaceAlpha', 0.25, 'EdgeAlpha', 0);

x = [100 130 130 100]; % x-coordinates of the vertices of the complex area
y = [-1 -1 1 1]; % y-coordinates of the vertices of the complex area
fill(x, y, 'white','FaceAlpha',.9, 'EdgeAlpha', 0); % shade the complex area in red


l(13) = line(data(250:5101,1),C_KC,"Color","#000000","LineWidth",0.5);

xl(10) = xlabel('Time (kyrs BP)');

b(15) = axes('Position',[0.1 0.085 0.8 0.865],...
    'XTickLabels',[],...
    'XGrid','On',...
    'Color','None', ...
    'XColor', 'none', 'XTick', [], ...
    "YAxisLocation","right", ...
    'YColor',"#000000", ...
    'Ylim', [0 1],...
    'Visible','Off');
l(9) = line(0,0);


x = [ind2(1) ind2(2) ind2(2) ind2(1)];
y = [0 0 1 1];
patch(x,y,'black', ...
    'FaceAlpha',0.05, 'EdgeColor','none')

x = [ind2(3) ind2(4) ind2(4) ind2(3)];
patch(x,y,'black', 'FaceAlpha',0.05, 'EdgeColor','none')

x = [ind2(5) ind2(6) ind2(6) ind2(5)];
patch(x,y,'black', 'FaceAlpha',0.05, 'EdgeColor','none')


xline(400,'Color',"#A2142F","LineWidth",.5,'LineStyle','--')
xline(150,'Color',"#A2142F","LineWidth",.5,'LineStyle','--')

sl = 350;
el = 450;
hl = .21;
hd = .005;
wd = .5;
line([sl el], [hl hl],'Color',"#7E2F8E","LineWidth",wd); % 
line([sl sl], [hl hl],'Color',"#7E2F8E","LineWidth",wd); % 
line([el el], [hl-hd hl+hd],'Color',"#7E2F8E","LineWidth",wd); % 
line([sl sl], [hl-hd hl+hd],'Color',"#7E2F8E","LineWidth",wd); % 

sl = 200;
el = 300;


line([sl el], [hl hl],'Color',"#D95319","LineWidth",wd); % 
line([sl sl], [hl hl],'Color',"#D95319","LineWidth",wd); % 
line([el el], [hl-hd hl+hd],'Color',"#D95319","LineWidth",wd); % 
line([sl sl], [hl-hd hl+hd],'Color',"#D95319","LineWidth",wd); % 

sl = 32;
el = 132;

line([sl el], [hl hl],'Color',"#0072BD","LineWidth",wd); % 
line([sl sl], [hl hl],'Color',"#0072BD","LineWidth",wd); % 
line([el el], [hl-hd hl+hd],'Color',"#0072BD","LineWidth",wd); % 
line([sl sl], [hl-hd hl+hd],'Color',"#0072BD","LineWidth",wd); % 


linkaxes(b,"x")
b(1).XLim = [0 550];

z(1) = zoom(fig03);
set(z(1),'Motion','horizontal');

p(1) = pan(fig03);
set(p(1),'Motion','horizontal');

hold off
