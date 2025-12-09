% Example code to plot SA survey data points projected onto
%  the Mop line and color coded by distance from the mop line

% Rather than interpolating a profile from gridded survey points,
% this code 
%
%  1) loads the 1m spatial averaged survey point data and
%  2) finds the distance of each from the mop line defined 
%      as a line of xshore points with 1m spacing.
%  3) for each xshore point that has one or more survey points
%      that are closer to it than the other mop line points,
%      find the closest of those survey points 
%
% The resulting profile variables are
%
%     Xuniq = the unique xshore distances with survey data (1m integers)
%     Znear = the elevation of the survey point nearest to each Xuniq point (m,
%             navd88)
%     Zdist = the distance of each Znear point from the mop line (m)
%
% It then makes mop line profile vectors x1d and z1di that use a 
%   distance from the mop line threshhold, dy, to define valid profile 
%   Znear data, and a maximum allowable spacing between Xuniq points, dx, 
%   that can be filled using 1d interpolation of the remaining Znear points,
%   if interpolated/regular spaced xshore profile info is needed.
 

% set path to CPG MOP files from your machine
% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% set the mop number of interest
MopNumber=743; % Beacons
MopNumber=900; % Oside

% load the SA mat file
matfile=sprintf('M%5.5iSA.mat',MopNumber);
load(matfile,'SA');

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');
 
%  Get thethe main transect line as UTM x,y points ~1m apart
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,1,[-100 0]);

% select a survey to analyize
%  Use USACE SHOALS data from 26 Oct 2009 as a complex example
nsurv=19;%22;
nsurv=size(SA,2);%last survey;

% find nearest subtransect line point to the input location
[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(SA(nsurv).Y),double(SA(nsurv).X)],'euclidean','smallest',1);

% define X based on the nearest transect line point
%   row=nearest subtransect number; col = xshore distance indice on
%   the nearest subtransect
[row,col] = ind2sub(size(xst),NearIdx); 

% xshore distance (1m xshore resolution) along the Mop transect for each survey point
X=x1d(col);

% For the xshore distances (1m xshore resolution) with data (Xuniq), 
%  find the nearest survey point elevation for each (Znear) and
%  how far away from the line it was (Zdist).

Xuniq=unique(X);

n=0;
for x=Xuniq
    n=n+1; % xshore point counter
    idx=find(X == x);
    [Zdist(n),imin]=min(dp(idx));
    Znear(n)=SA(nsurv).Z(idx(imin));   
end

% example to make an interpolated profile based on distance from 
% mop line (dy) and xshore spatial gap (dx) tolerances
dy=10; % alongshore distance from mop line tolerance
dx=5; % xshore date gap tolerance

% set nearest points with Zdist > dy to NaNs
Znear(Zdist > dy)=NaN;

% seed the regular spaced profile vectors (x1d and z1d) with the 
%  valid nearest point elevations (both valid and NaNs)
ndx=find(ismember(x1d,Xuniq));
z1d=x1d*NaN;z1d(ndx)=Znear;

% identify the size of the gap each xshore point is in (0 = not in a gap)
sz=gapsize(z1d);

% make vectors of valid data points to be used in interpolation
%  where the x,z profile points in gaps <= dx are removed. NaN
%  no data points in the gaps larger than dx are kept.

x1dv=x1d;x1dv(isnan(z1d) & sz <= dx)=[];
z1dv=z1d;z1dv(isnan(z1d) & sz <= dx)=[];

% interpolate to fill the small gaps
z1di=interp1(x1dv,z1dv,x1d);




close all
figure('position',[440    97   629   700]);
subplot(3,1,1)
[ScatterPlot,ColorBarPlot]=ColorScatterProfile(X,SA(nsurv).Z,dp);
xlabel('Mop Line Xshore Distance (m)');
ylabel('Elevation (m, NAVD88)');
grid on;set(gca,'xdir','reverse','fontsize',12);
title({['Mop Area' num2str(MopNumber)...
    '; All 1m Spatial Avg Survey Points Projected on Mop Line']}...
    ,{[SA(nsurv).Source ' ' datestr(SA(nsurv).Datenum)]});
ColorBarPlot.Label.String='Distance From Mop Line (m)';
hold on;
pn=plot(Xuniq,Znear,'m.','markersize',8);
legend(pn,'Nearest Survey Points to Line w/ 1m xshore resolution','location','northwest');
ColorBarPlot.Location='south';
xl=get(gca,'xlim');

subplot(3,1,2)
plot(Xuniq,Zdist,'m.','markersize',10);
grid on;set(gca,'xdir','reverse','fontsize',12,'xlim',xl);
xlabel('Mop Line Xshore Distance (m)');
ylabel('Distance from Mop Line (m)');
title({'Nearest 1m Spatial Avg Survey Point Distances '},...
    {'from the Mop Line defined with 1m xshore resolution points'})


subplot(3,1,3)
plot(x1d,z1di,'k-','linewidth',2);
grid on;set(gca,'xdir','reverse','fontsize',12,'xlim',xl);
xlabel('Mop Line Xshore Distance (m)');
ylabel('Elevation (m, NAVD88)');
title({'Interpolated Profile using Nearest Points '},...
    {['Max distance from Mop Line Threshold =' num2str(dy)...
    'm ; Max xshore filled gap threshold =' num2str(dx) 'm']})


function [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,NumTrans,ExtendLine)

%  Used by BuildSMmatfiles.m

%  Calculates the interpolated 1m xshore resolution transect line, and NumTrans 
%   subtransect lines, for/within the MopNumber survey area bounds. Extends
%   the lines -/+ ExtendLine(1)/ ExtendLine(2) meters from the 
%   back beach (-) and offshore (+) mop line, using the Mop table for
%   reference.

% eg. for 20 subtransects that extend -100m back from back beach
%      mop line:
%
%  [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,580,20,[-100 0])

%  returns transects all of the same length = longest of the
%  subtransects


%-------------------------------------------------------------
% first, the subtransect lines
%-------------------------------------------------------------
% tweak offshore mop point locations to be at .999 of the distance of the
 % actual mop transect to prevent neighboring mops sharing the exact same 
 % offshore mop location.
 
 for m=MopNumber-1:MopNumber+1
     
     % check if back and offshore points are same (in a river channel)
     if Mop.BackXutm(m) == Mop.OffXutm(m)
         Mop.OffXutm(m)=Mop.OffXutm(m)-1; % tweak x offshore value if so
     end
         
  xa=Mop.BackXutm(m)+0.999*(Mop.OffXutm(m)-Mop.BackXutm(m));
  xyratio=(xa-Mop.BackXutm(m))/(Mop.OffXutm(m)-Mop.BackXutm(m));
  ya=Mop.BackYutm(m)+xyratio*(Mop.OffYutm(m)-Mop.BackYutm(m));
  Mop.OffXutm(m)=xa;
  Mop.OffYutm(m)=ya;
 end

 % back beach neighboring mop midpoints
 xb1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
 yb1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
 xb2=Mop.BackXutm(MopNumber);
 yb2=Mop.BackYutm(MopNumber);
 xb3=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
 yb3=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
 % offshore neighboring mop midpoints
 xo1=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
 yo1=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);
 xo2=Mop.OffXutm(MopNumber);
 yo2=Mop.OffYutm(MopNumber);
 xo3=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
 yo3=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
 
 nst=NumTrans; % number of subtransects 
 % increase nst by 2 when calculating x,y step size to keep
 %  first and last subtransects inside mop area boundaries
 %  for 2d interpolation. 
%  xbst=linspace(xb1,xb2,nst+2);xbst=xbst(2:end-1);
%  ybst=linspace(yb1,yb2,nst+2);ybst=ybst(2:end-1);
%  xost=linspace(xo1,xo2,nst+2);xost=xost(2:end-1);
%  yost=linspace(yo1,yo2,nst+2);yost=yost(2:end-1);
 
 [xbst,ybst]=EqualSpacedPoints([xb1 xb2 xb3],[yb1 yb2 yb3],NumTrans);
 
 [xost,yost]=EqualSpacedPoints([xo1 xo2 xo3],[yo1 yo2 yo3],NumTrans);
 
%  270-atan2d(yb2-yo2,xb2-xo2)
%  270-atan2d(yb2-yo2+2,xb2-xo2)
%  270-atan2d(yb2-yo2-2,xb2-xo2)
%  270-atan2d(ybst-yost,xbst-xost)
%  mean(270-atan2d(ybst-yost,xbst-xost))
 
 % get max subtransect length
 dmax=0;
 for n=1:nst
     dist=pdist([xbst(n),ybst(n);xost(n),yost(n)]);
     dmax=max([dmax dist]);
 end
 
 % include the actual Mop transect distance 
 x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
 x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
 dist=pdist([x1,y1;x2,y2]);
 
 dmax=max([dmax dist]); % max of the transect & subtransects lengths
 dmax=ceil(dmax); % round up to whole meter
 

 % now loop through subtransect endpoints and make 1m xshore interpolated
 %  transects
 for n=1:nst
   x1=xbst(n);y1=ybst(n);x2=xost(n);y2=yost(n);
   ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
   dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
   % adjust offshore extension so all transects have same xshore length
   ExtendOff=ExtendLine(2)+(dmax-dist);
 
   % 1m spaced points along transect 
   % x coords of 1m spaced points along transect
   xst(n,:)=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
   % y coords of 1m spaced points along transect
   yst(n,:)=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
 
 end
 
 % xshore distance m from back beach point
 x1d=ExtendLine(1):dmax+ExtendLine(2);
 
%-------------------------------------------------------------
% now the single Mop transect line
%-------------------------------------------------------------
    
 % Mop transect 1m spaced xshore points
 x1=Mop.BackXutm(MopNumber);y1=Mop.BackYutm(MopNumber);
 x2=Mop.OffXutm(MopNumber);y2=Mop.OffYutm(MopNumber);
 ang=atan2(y2-y1,x2-x1); % radian angle between back and offshore point
 dist=pdist([x1,y1;x2,y2]); % distance between back and offshore point
 % adjust offshore extension so all transects have same xshore length
 ExtendOff=ExtendLine(2)+(dmax-dist);
 
 % x coords of 1m spaced points along transect
 xt=x1+ExtendLine(1)*cos(ang):cos(ang):x2+ExtendOff*cos(ang);
 % y coords of 1m spaced points along transect
 yt=y1+ExtendLine(1)*sin(ang):sin(ang):y2+ExtendOff*sin(ang);  
 
end

function sz=gapsize(x)
% Calculates the gap size of a vector. 
%
% A gap is defined as the number of consequtive nans.
%
% USAGE: sz=gapsize(x)
% sz is same length as input. 
%
% example:
% x=rand(20,1);
% x(x>.5)=nan; 
% [x gapsize(x)]
%
% aslak grinsted 2010
x=~isnan(x);
hasdata=[0;find(x(:)); length(x)+1];
sz=zeros(size(x));
for ii=1:length(hasdata)-1
    ix=hasdata(ii)+1:hasdata(ii+1)-1;
    sz(ix)=length(ix);
end

end

function [ScatterPlot,ColorBarPlot]=ColorScatterProfile(x,y,z)

[z,i]=sort(z,'descend');
x=x(i);
y=y(i);

zmin=quantile(z,.05);zmax=quantile(z,.95);

zrange=zmax-zmin; % set max slope for coloring
zscaled = 1+64*(z-zmin)/zrange;
zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
idx=find(~isnan(zscaled)); % non NaN points 
                                     
cm =jet(64);

scp=scatter(x(idx), y(idx), 7, cm(ceil(zscaled(idx)),:), 'filled');
scp.MarkerFaceAlpha = .9;
scp.MarkerEdgeAlpha = .9;
%view(-10,80)
colormap(jet(64))

cb=colorbar;cb.Label.String='Dist From Mop Line (m)';
set(gca,'clim',[zmin zmax]);
set(gca,'xlim',[min(x) max(x)]);
%set(gca,'zlim',[min(z) max(z)]);
set(gca,'color',[.7 .7 .7]);

ScatterPlot=scp;
ColorBarPlot=cb;

end

