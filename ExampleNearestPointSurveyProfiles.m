% example code to plot all the iG8wheel profiles for a Mop line using
%  the nearest point profile method

load('MopTableUTM.mat','Mop');

MopNumber=510;
    
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])

idx=find(strcmp({SA.Source},'iG8wheel'));

figure('position',[ 116         159        1215         531]);
m=0;
for n=idx
    m=m+1;

     [x1d,z1di]=GetNonGriddedProfile(Mop,SA,MopNumber,n);
    
     p(m)=plot(x1d,z1di,'-',...
        'DisplayName',datestr(SA(n).Datenum),...
        'linewidth',1);hold on;

end

grid on;
set(gca,'fontsize',12);xlabel('Distance From Mop Back Beach Point (m)');ylabel('Elevation (m, NAVD88)')
set(gca,'ylim',[0 3.5],'xdir','reverse','fontsize',14);

title([{['Mop ' num2str(MopNumber)]},{'iG8wheel profiles' }],...
   'fontsize',16);
xl=get(gca,'xlim');
plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);  
legend(p,'location','eastoutside','fontsize',12,'numcolumns',3);

function [x1d,z1di]=GetNonGriddedProfile(Mop,SA,MopNumber,nsurv)

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
 
%  Get thethe main transect line as UTM x,y points ~1m apart
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,1,[-100 0]);

% select a survey to analyize
%  Use USACE SHOALS data from 26 Oct 2009 as a complex example
% nsurv=19;%22;
% nsurv=size(SA,2);%last survey;

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
dy=25; % alongshore distance from mop line tolerance
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

end

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
 %if MopNumber == 1;dm=0;else;dm=1;end
 
 dm=1;
 if MopNumber == 1;MopNumber =2 ;end
 for m=MopNumber-dm:MopNumber+1
     
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
 xb1=mean([Mop.BackXutm(MopNumber-dm) Mop.BackXutm(MopNumber)]);
 yb1=mean([Mop.BackYutm(MopNumber-dm) Mop.BackYutm(MopNumber)]);
 xb2=Mop.BackXutm(MopNumber);
 yb2=Mop.BackYutm(MopNumber);
 xb3=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
 yb3=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
 % offshore neighboring mop midpoints
 xo1=mean([Mop.OffXutm(MopNumber-dm) Mop.OffXutm(MopNumber)]);
 yo1=mean([Mop.OffYutm(MopNumber-dm) Mop.OffYutm(MopNumber)]);
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

function [marker_x,marker_y]=EqualSpacedPoints(x,y,n)

% used by GetTransectLines.m , which is used by BuildSMmatfiles.m
   
%    x = [0 5 10];
%    y = [2 4 1];
%    n=10;
   dist_from_start = cumsum( [0, sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2)] );
   marker_dist = dist_from_start(end)/n;
   marker_locs = marker_dist/2 : marker_dist : dist_from_start(end);   %replace with specific distances if desired
   marker_indices = interp1( dist_from_start, 1 : length(dist_from_start), marker_locs);
   marker_base_pos = floor(marker_indices);
   weight_second = marker_indices - marker_base_pos;
   marker_x = x(marker_base_pos) .* (1-weight_second) + x(marker_base_pos+1) .* weight_second;
   marker_y = y(marker_base_pos) .* (1-weight_second) + y(marker_base_pos+1) .* weight_second;
%    figure;
%    plot(x, y);
%    hold on;
%    plot(marker_x, marker_y, 'r+');
%    hold off
   
   end



 