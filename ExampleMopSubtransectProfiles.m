% Example code to get subtransect profiles from CPG MOP
%   gridded SG.mat files

% set paths to cpg mop files and toolbox m-scripts
addpath /volumes/group/MOPS
addpath /volumes/group/MOPS/toolbox

% load the 1m gridded (MOPS/*SG.mat) survey data for mop 276 along point loma
MopNumber=276;

matfile=['M' num2str(MopNumber,'%5.5i') 'SG.mat' ];
load(matfile,'SG'); % load the SG struct array
% SG contains grids for all the historical surveys in the CPG MOP
% database for that Mop.  Some gridded data is better than others depending
% on the quality and spatial density of the underlying survey data.
fprintf('\nSG struct array contains the following surveys:\n\n')
% list survey types and dates 
for n=1:size(SG,2)
 fprintf('%s %s\n',datestr(SG(n).Datenum),SG(n).Source)
end

%% plot the first survey in the gridded data set
figure('position',[59   202   682   523]);
SurvNum=1;
PlotCMstruct(SG,SurvNum); % use toolbox function to view grid
set(gca,'fontsize',12);
title([ 'Mop: ' num2str(MopNumber) ' '...
    SG(SurvNum).Source ' ' datestr(SG(SurvNum).Datenum)],...
    'fontsize',16);

%% Subdivide the mop area into 20 cross-shore transects that don't overlap
load('MopTableUTM.mat','Mop'); % load "Mop" table variable of Mop transect information
% get 20 subtransets for this MopNumber. 
% Extend the transect line -100m landward of the official Mop back beach
%  point and 0m seward of the Mop offshore point to include bluffs
ExtendLine=[-100 0];
NumSub=20; % number of subtransects
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,NumSub,ExtendLine);
% All the returned subtransects are the same length, n, of 1m spaced points
% x1d(n) = the xshore distance from the Mop back beach of the n transect points
% xt(n) = UTM eastings of the original Mop transect line points 
% yt(1:n) = UTM northings of the original Mop transect line points 
% xst(NumSub,n) = UTM eastings matrix of the NumSub subtransect line points
% yst(NumSub,n) = UTM northings matrix of the NumSub subtransect line points

%% plot the subtransects
figure('position',[89   170   682   523]);
PlotCMstruct(SG,SurvNum);hold on;
for n=1:NumSub
    plot(xst(n,:),yst(n,:),'g-')
end
set(gca,'fontsize',12);
title([ 'Mop: ' num2str(MopNumber) ' with ' num2str(NumSub) ' Subtransects'],...
    'fontsize',16);
%% This section interpolates gridded data onto the subtransects

if ~isempty(SG(SurvNum).X) % check to make sure it isn't an empty grid
    
    % First, reconstruct a full 2d grid from the valid grid point data
    %  kept in the SG struct array
    xg=SG(SurvNum).X;yg=SG(SurvNum).Y;zg=SG(SurvNum).Z; %valid grid points
    nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
    tg=nan(ny,nx); % temp grid of Nans
    idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
    tg(idx)=zg; % add valid data to temp grid
    [X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1); % X,Y grid indices

    % now get 2d interpolated z values for all the subtransect points
    zst=xst*NaN; % initialize transect elevation points 
    zst(:) = interp2(X,Y,tg,xst(:),yst(:));  
    % zst(NumSub,:) now contains elevations for all the subtransects
    %  with UTM coords xst and yst. Points with no data = NaNs
    
end

%% make a plot with all the subtransect profiles on it. For the cross-shore
%   coordinate use the distance from the Mop back beach line which is in
%   in the vector x1d returned by GetTransectLines

figure('position',[129   140   682   523]);
hold on;
for n=1:NumSub
    plot(x1d,zst(n,:),'-','DisplayName',['Transect ' num2str(n)]);
end
set(gca,'xdir','reverse');grid on;
legend('location','eastoutside');
xlabel('Cross-shore Distance from Mop Back Beach Line (m)');
ylabel('Elevation (m, NAVD88)');
set(gca,'fontsize',12);
title([ 'Mop: ' num2str(MopNumber) ' '...
    SG(SurvNum).Source ' ' datestr(SG(SurvNum).Datenum)],...
    'fontsize',16);


