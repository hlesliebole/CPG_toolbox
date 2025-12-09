% Example code to take input Coastal Frontiers (CF) profile Xutm, Yutm data points and
%
% 1. Find the nearest Mop transect(s) to the CF data. CF profile data can
%  potentially be split onto two adjacent Mop transects if the CF profile
%  orientation is significantly different than the Mops.
%
% 2. Project all the CF profile data onto the single Mop transect with the 
%     most CF points nearest to it.
%
% 3. Plot the CF data as a Mop profile on top of all the historical
%     Mop profiles at that transect.
%

% add paths to MOPS and MOPS/toolbox to use the function
% eg. for a mac
addpath /Volumes/group/MOPS  % folder with MOP mat files
addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

%%  You want to load the CF Xutm Yutm Znavd88 data for a single profile survey date here
%     For now, using the most recent wheelie data at Mop 870 saved as Xutm,Yutm,Znavd88
%     saved in TestUTM.mat to test the code
load('TestUTM.mat','Xutm','Yutm','Znavd88')

%% Find the nearest Mop(s) to the data
Nmops=FindNearestMopTransectsUTM(Xutm,Yutm);

%% Select the Mop number with the most nearest points using the mode function
Nmop=mode(Nmops);
fprintf('Nearest Mop is: %i\n',Nmop)

%% project the CF points onto the Nmop transect to make a profile

% set nearest point profile tolerances 
Ytol=50; % max dist (m) a survey point can be from the transect
Xtol=5; % max alongtransect gap (m) that will filled with linear interpolation

% GetNearestPointsProfile returns:
%   Xmop = vector of 1m xshore grid points (fixed range starting -100m back
%         point and extending out to the 10m depth offshore Mop point) 
%   Zmop = vector of the nearest point profile elevations
%       
%   for input CF Xutm, Yutm, Znavd88 survey data points

[Xmop,Zmop]=GetNearestPointsProfile(Nmop,Ytol,Xtol,Xutm,Yutm,Znavd88);

%% Now load all the exist data for Nmop in CPG SA file
load(['M' num2str(Nmop,'%5.5i') 'SA.mat'],'SA');

% get nearest point profiles for all the sur ey sin the SA file
[X1D,Z1D]=GetAllNearestPointsProfiles(SA,Ytol,Xtol);

%% Make example plot

idx=1:size(SA,2); % plot all surveys

col=jet(numel(idx));
m=0;
figure('position',[102         181        1249         531])
% plot survey profiles in database
for n=idx
 m=m+1;
 pl(m)=plot(X1D,Z1D(n,:),'-','color',col(m,:),'linewidth',2,'DisplayName',...
        datestr(SA(n).Datenum,'mm/dd/yy'));hold on;
end
% now overlay new CF profile
m=m+1;
pl(m)=plot(X1D,Z1D(n,:),'-','color','k','linewidth',3,'DisplayName',...
        'CF Data Profile');hold on;

set(gca,'ylim',[-12 6],'fontsize',14);
xl=get(gca,'xlim');
plot(xl,xl*0+1.344,'k--');text(xl(2)-5,1.65,'MHW','fontsize',12);
set(gca,'xdir','reverse');
lg=legend(pl,'location','eastoutside','numcolumns',4);
title(lg,'Survey Dates')
box on;grid on;
title(['Mop ' num2str(Nmop) ' Profiles']);
xlabel('Xshore Distance (m)');
ylabel('Elev (m, NAVD88)');


function [X1D,Z1D]=GetNearestPointsProfile(Nmop,Ytol,Xtol,Xutm,Yutm,Znavd88)

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');
 
%  Get the Nmop main transect line as UTM x,y points ~1m apart (xt,yt)
%   This also sets the fixed x1d transect point starting
%   -100m behind the mop back beach point out to the offshore point
[X1D,xt,yt,xst,yst]=GetTransectLines(Mop,Nmop,1,[-100 0]);

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xst),NearIdx); 

% xshore distance (1m xshore resolution) along the Mop transect for each survey point
X=X1D(col);

% For the xshore distances (1m xshore resolution) with data (Xuniq), 
%  find the nearest survey point elevation for each (Znear) and
%  how far away from the line it was (Zdist).

Xuniq=unique(X);

n=0;
for x=Xuniq
    n=n+1; % xshore point counter
    idx=find(X == x);
    [Zdist(n),imin]=min(dp(idx));
    Znear(n)=Znavd88(idx(imin));   
end

% make an interpolated profile based on distance from 
% mop line (dy) and xshore spatial gap (dx) tolerances

% set nearest points with Zdist > Ytol to NaNs
Znear(Zdist > Ytol)=NaN;

% seed the regular spaced profile vectors (x1d and z1d) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(X1D,Xuniq));
z1d=X1D*NaN;z1d(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(z1d);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= dx are removed. NaN
%  no data points in the gaps larger than dx are kept.

x1dv=X1D;x1dv(isnan(z1d) & sz <= Xtol)=[];
z1dv=z1d;z1dv(isnan(z1d) & sz <= Xtol)=[];

% interpolate to fill the small gaps
%z1di=interp1(x1dv,z1dv,x1d);
Z1D=interp1(x1dv,z1dv,X1D);

end


