% Example code to make a plot of the most recent 
% 1m gridded data from a truck survey in a Mops SG mat file.
%
% The code needs access to the reefbreak group/MOPS and
%  group/MOPS/toolbox directories so you need to add paths
%  before calling this script
%
% eg. for a mac
% addpath '/volumes/group/MOPS/'; % path to cpg mop mat files
% addpath '/volumes/group/MOPS/toolbox'; % path to cpg mop mat files

% pick a mop number
MopNumber=553;  % blacks north (wide terrace and back beach area)

%% load the SG mat file
matfile=['M' num2str(MopNumber,'%5.5i') 'SG.mat'];
fprintf('\nLoading %s ...\n',matfile)
load(matfile,'SG');

%% identify jumbos by checking the original survey file names for the 
%  word jumbo
%jumbo=find(contains({SG.File},'umbo'));

%% identify truck surveys by matching .Source field name
trk=find(strcmpi({SG.Source},'Trk'));

fprintf('The SG struct array has %i Truck Surveys.\n',numel(trk))
fprintf('Starting  %s , and ending %s\n',datestr(SG(trk(1)).Datenum),datestr(SG(trk(end)).Datenum))

fprintf('Plotting the most recent Truck Survey in the database.\n',numel(trk))

% Just look at the most recent truck survey
SurvNum=trk(end); 

%% The SG struct array only save 1m spatial resolution x,y grid points that
%   have valid elevation data for this mop transect area (+/-50m in the
%   alongshore on each side of the transect line).

%  Need to turn the saved SG struct array grid points into a 2d grid array

% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;

% Make 2d x,y utm arrays encompassing the valid points
[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
idx=sub2ind(size(X),y-min(y)+1,x-min(x)+1); % 1d indices of 2d utm grid

Z=X*NaN; % initialize the 2d elevation Z array as NaNs
Z(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices

% X,Y, and Z are now 2D arrays. X,Y points with no valid data have Z=NaN

%% --------------------------------------------------------------
%% plot the xutm,yutm,z 3d surface
figure('position',[32 143 1010  652]);
p1=surf(X,Y,Z,'linestyle','none');
set(gca,'view',[ -10.3232   26.1279]);
BeachColorbar; % use designer beach profile color scale    
grid on;
set(gca,'color',[.75 .75 .75],'fontsize',14);

y_labels = get(gca, 'YTick');
set(gca,'YTick',y_labels);
set(gca, 'YTickLabel', num2str(y_labels'));ytickangle(45);
ylabel('northings (m)');
x_labels = get(gca, 'XTick');
set(gca,'XTick',x_labels);
set(gca, 'XTickLabel', num2str(x_labels'));xtickangle(45);
xlabel('eastings (m)');
zlabel('Elevation (m.NAVD88)');

title([' Mop ' num2str(MopNumber) ' Truck Survey: ' datestr(SG(SurvNum).Datenum)],'fontsize',16)
