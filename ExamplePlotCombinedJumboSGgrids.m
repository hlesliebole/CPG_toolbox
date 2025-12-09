% Example code to make a plot of the 1m gridded data from a Jumbo
%  survey for a combined range of mop SG mat file in a new CG struct array.

%% Uses the m-script function CG=CombineSGdata(mapth,MopStart,MopEnd)

%% settings
%mpath='/volumes/group/MOPS/'; % path to cpg mop files
mpath='/Users/William/Desktop/MOPS/';
MopStart=666%666;%582; % start mop number
MopEnd=683%683;%590; % ending mop number

%% load the combined SG mat file data.  Return the combined
%  data to a struct array SG instead of CG to use the normal 
%  single mop SG gridding and plotting code with the combined data.
SG=CombineSGdata(mpath,MopStart,MopEnd);

%% identify jumbos by checking the original survey file names for the 
%  word jumbo
jumbo=find(contains({SG.File},'umbo'));

fprintf('The SG struct array has %i Jumbo Surveys.\n',numel(jumbo))

fprintf('Plotting the most recent Jumbo Survey in the database.\n',numel(jumbo))



%% Need to turn the saved SG struct array grid points into a 2d grid array


% make grid area encompassing all survey data
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));
%[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
[X,Y]=meshgrid(minx:maxx,miny:maxy);

figure('position',[ 21          56        1385         728]);

%--- a make regular grid for a specific survey

%% 2016 

% first survey

didx=find([SG(jumbo).Datenum] ==  datenum(2015,11,12));

SurvNum=jumbo(didx); 
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices

% nextjumbo survey
didx=find([SG(jumbo).Datenum] ==  datenum(2016,2,9));
SurvNum2=jumbo(didx);

x=SG(SurvNum2).X;
y=SG(SurvNum2).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices


% plot the xutm,yutm,z 3d surface

ax1=axes('position',[0.05    0.0500    0.1504    0.95]);
p1=surf(X,Y,Z2-Z1,'linestyle','none');

set(gca,'xlim',[473300 474100]);hold on;
set(gca,'dataaspectratio',[1 1 1])

grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);

y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');
zlabel('Elevation (m.NAVD88)');

view(2)
%axis equal

d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];
text(473350,3653120,str,'fontweight','bold','fontsize',18);

set(gca,'clim',[-3 3]);
polarmap;

title({
    datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
    'fontsize',16)

%% -------------------------------------
%% 2017 

% first survey

didx=find([SG(jumbo).Datenum] ==  datenum(2016,9,29));

SurvNum=jumbo(didx); 
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices

% nextjumbo survey
didx=find([SG(jumbo).Datenum] ==  datenum(2017,1,26));
SurvNum2=jumbo(didx);

x=SG(SurvNum2).X;
y=SG(SurvNum2).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices

% plot the xutm,yutm,z 3d surface

ax2=axes('position',[0.35    0.0500    0.1504    0.95]);
p2=surf(X,Y,Z2-Z1,'linestyle','none');

set(gca,'xlim',[473300 474100]);hold on;
set(gca,'dataaspectratio',[1 1 1])

grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);

y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');
zlabel('Elevation (m.NAVD88)');

view(2)
%axis equal

d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];
text(473350,3653120,str,'fontweight','bold','fontsize',18);

set(gca,'clim',[-3 3]);
polarmap;

title({
    datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
    'fontsize',16)

%% -------------------------------------

%% -------------------------------------
%% 2021 

% first survey

didx=find([SG(jumbo).Datenum] ==  datenum(2020,10,1));

SurvNum=jumbo(didx); 
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices

% nextjumbo survey
didx=find([SG(jumbo).Datenum] ==  datenum(2021,1,28));
SurvNum2=jumbo(didx);

x=SG(SurvNum2).X;
y=SG(SurvNum2).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices

% plot the xutm,yutm,z 3d surface

ax3=axes('position',[0.35    0.0500    0.1504    0.95]);
p2=surf(X,Y,Z2-Z1,'linestyle','none');

set(gca,'xlim',[473300 474100]);hold on;
set(gca,'dataaspectratio',[1 1 1])

grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);

y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');
zlabel('Elevation (m.NAVD88)');

view(2)
%axis equal

d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];
text(473350,3653120,str,'fontweight','bold','fontsize',18);



set(gca,'clim',[-3 3]);
polarmap;

title({
    datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
    'fontsize',16)

%% -------------------------------------
%% -------------------------------------
%% 2023 

% first survey

didx=find([SG(jumbo).Datenum] ==  datenum(2022,10,12));

SurvNum=jumbo(didx); 
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices

% nextjumbo survey
didx=find([SG(jumbo).Datenum] ==  datenum(2023,1,23));
SurvNum2=jumbo(didx);

x=SG(SurvNum2).X;
y=SG(SurvNum2).Y;

% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);

Z2=X*NaN; % initialize the 2d elevation Z array as NaNs
Z2(idx)=SG(SurvNum2).Z; % overlay the valid data points using the 1d indices

% plot the xutm,yutm,z 3d surface

ax4=axes('position',[0.35    0.0500    0.1504    0.95]);
p2=surf(X,Y,Z2-Z1,'linestyle','none');

set(gca,'xlim',[473300 474100]);hold on;
set(gca,'dataaspectratio',[1 1 1])

grid on;
set(gca,'color',[.7 .7 .7],'fontsize',14);

y_labels = get(gca, 'YTick');
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');
zlabel('Elevation (m.NAVD88)');

view(2)
%axis equal

d=Z2-Z1;erode=round(sum(d(d<0))/1000);dep=round(sum(d(d>0))/1000);
str=['-V= ' num2str(erode) 'K  ; +V= ' num2str(dep) 'K'];
text(473350,3653120,str,'fontweight','bold','fontsize',18);

set(gca,'clim',[-3 3]);
polarmap;

% title({
%     datestr(SG(SurvNum).Datenum), ' to ', datestr(SG(SurvNum2).Datenum)},...
%     'fontsize',16)

title({
    datestr(SG(SurvNum).Datenum), ' to ','23-Jan-2023'},...
    'fontsize',16)

%t=text(ax4,473550,3652000,10,{'PLEASE','STAND BY'},'fontweight','bold','fontsize',20)
%% -------------------------------------

cb=colorbar;
cb.Label.String='Elevation Change (m)';
cb.FontSize=14;

set(ax1,'position',[0.075    0.0500    0.1504    0.95]);
set(ax2,'position',[0.3    0.0500    0.1504    0.95]);
set(ax3,'position',[0.525    0.0500    0.1504    0.95]);
set(ax4,'position',[0.75    0.0500    0.1504    0.95]);
colormap(flipud(colormap))

makepng('CardiffSevereWinterChanges.png')