% Example code to display the historical spatial coverage of
% 1m gridded data from a truck survey in a Mops SG mat file.
%
% The code needs access to the reefbreak group/MOPS and
%  group/MOPS/toolbox directories, so you need to add paths
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

%% identify truck surveys indices by matching .Source field name
trk=find(strcmpi({SG.Source},'Trk'));

fprintf('The SG struct array has %i Truck Surveys.\n',numel(trk))
fprintf('Starting  %s , and ending %s\n',datestr(SG(trk(1)).Datenum),datestr(SG(trk(end)).Datenum))

fprintf('Plotting the most recent Truck Survey in the database.\n',numel(trk))

SurvNum=trk(1); % consider all the truck surveys
%% The SG struct array only save 1m spatial resolution x,y grid points that
%   have valid elevation data for this mop transect area (+/-50m in the
%   alongshore on each side of the transect line).

%  Need to turn the saved SG struct array grid points into a 2d grid array

% Mop area 1m grid points with valid data for all truck surveys
x=vertcat(SG(trk).X);
y=vertcat(SG(trk).Y);

% universal grid limits
xmin=min(x);xmax=max(x);ymin=min(y);ymax=max(y);

% Make 2d x,y utm arrays encompassing the valid points of all the truck
% surveys
[X,Y]=meshgrid(min(x):max(x),min(y):max(y));

%% now build 3D Z array where the 3rd dim is the truck survey number

% initialize the 3d elevation Z array as NaNs
Z=NaN(size(X,1),size(X,2),numel(trk)); 

nt=0; % number of truck survey counter
for SurvNum=trk

    Z2d=NaN(size(X)); % initialize temp 2d z array
    % get valid x,y gid point for this survey
    x=SG(SurvNum).X;
    y=SG(SurvNum).Y;

  idx=sub2ind(size(X),y-ymin+1,x-xmin+1); % 1d indices of universal 2d utm grid

  Z2d(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
  % add to 3d matrix
  nt=nt+1;
  Z(:,:,nt)=Z2d;
end

% get truck gridded survey stats
NumSurv2d=sumnd(~isnan(Z),3); % number surveys at each grid point
GlobalMeanZ2d=meannd(Z,3); % mean elevation at each grid point
GlobalMedianZ2d=mediannd(Z,3); % median elevation at each grid point
GlobalStdZ2d=stdnd(Z,3); % elevation standard deviation at each grid point

figure('position',[48   154   804   570]);
surf(NumSurv2d);view(2);colormap(jet);colorbar;
xlabel('East (m)');ylabel('North (m)');set(gca,'fontsize',14);
title('Number of Truck Surveys with Valid Data')

figure('position',[68   154   804   570]);
surf(GlobalMeanZ2d);view(2);colormap(jet);colorbar;
xlabel('East (m)');ylabel('North (m)');set(gca,'fontsize',14);
title('Mean Elevations (m, navd88) for All Truck Surveys with Valid Data')

figure('position',[88   154   804   570]);
surf(GlobalMedianZ2d);view(2);colormap(jet);colorbar;
xlabel('East (m)');ylabel('North (m)');set(gca,'fontsize',14);
title('Median Elevations (m, navd88) for All Truck Surveys with Valid Data')

figure('position',[108   154   804   570]);
surf(GlobalStdZ2d);view(2);colormap(jet);colorbar;
xlabel('East (m)');ylabel('North (m)');set(gca,'fontsize',14);
title('Elevation Standard Deviations (m) for All Truck Surveys with Valid Data')


function M=sumnd(M,dim)
s=size(M);
M=permute(M,[setdiff(1:ndims(M),dim),dim]);
M=reshape(M,[],s(dim));
M=sum(M,2,'omitnan');
s(dim)=1;
M=reshape(M,s);
end

function M=meannd(M,dim)
s=size(M);
M=permute(M,[setdiff(1:ndims(M),dim),dim]);
M=reshape(M,[],s(dim));
M=mean(M,2,'omitnan');
s(dim)=1;
M=reshape(M,s);
end

function M=mediannd(M,dim)
s=size(M);
M=permute(M,[setdiff(1:ndims(M),dim),dim]);
M=reshape(M,[],s(dim));
M=median(M,2,'omitnan');
s(dim)=1;
M=reshape(M,s);
end

function M=stdnd(M,dim)
s=size(M);
M=permute(M,[setdiff(1:ndims(M),dim),dim]);
M=reshape(M,[],s(dim));
M=std(M,0,2,'omitnan');
s(dim)=1;
M=reshape(M,s);
end
