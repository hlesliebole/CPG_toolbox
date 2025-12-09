function CP=SM2CP(SM,ShorelineElev,UpperShorefaceElev,SmoothWindow)
% clear all
% load M00582SM.mat
% ShorelineElev=0.774;
% UpperShorefaceElev=1.566;
% SmoothWindow=1;

% eg. CP=SM2CP(SM,0.774,1.566,2) ; with msl=0.744, mhhw=1.566

% Returns a struct array of "Characteristic Xshore Profiles" (CP) for a MOP
% line, where the profiles are defined as a function of the MSL shoreline 
% location.
  

% input:  
%   1. SM = 1 x Nsurvey survey morpho parmater struct array
%   2. ShorelineElev = navd88 elevation that defines the shoreline location
%   3. UpperShorefaceElev = navd88 elevation used with the shoreline elev
%                          to define the shoreface slope for runup calcs
%   4. SmoothWindow = MSL X location dimension smoothing window (meters)


% returns: 
%   CP 1 x Mprofiles characteristic profile struct array
%     where M is the number of characteristic profiles over the 
%     observed range of MSL shoreline locations with 1m xshore resolution.
%   CP(m).Xshoreline = mop profile shoreline xshore location (whole meters);
%                      Same as beach width base on Mop backbeach point
%                      location.
%   CP(m).ShorefaceSlope = mop profile slope;
%   CP(m).X1D = profile x locations (meters, relative to back beach point);
%   CP(m).Z1Dmean = profile elevations (meters);

% dependencies:

% Uses intersections.m to find intersection of the profiles with different
%  elevation levels

%----------------------------------------------------------------

%  find longest xshore profile
  xmax=max([SM.X1D]);  % max crossshore location

% loop thru survey mean profiles. Remove any with insufficient data
%  to identify the shoreline location or calculate a shoreface slope.
%  Temporarily append the shoreline x location, shoreface slope , 
%  and subarial profile area values to each survey in SM.
  
for ns=1:size(SM,2) % loop thru SM surveys
    
    % clip off -x points and buffer seaward mean profiles with NaNs so they 
    % all have the same xshore range between x=0 and xmax
    SM(ns).Z1Dmean(SM(ns).X1D < 0)=[];
    SM(ns).X1D(SM(ns).X1D < 0)=[];
    if SM(ns).X1D(end) < xmax
     for n=SM(ns).X1D(end)+1:xmax
         SM(ns).X1D(n)=n;
         SM(ns).Z1Dmean(n)=NaN;
     end
    end
      
    % find shoreline elevation x location using intersections.m
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)]; % ends of profile
    sel=[ShorelineElev ShorelineElev]; % shoreline elevation  
    xz=intersections(xl,sel,SM(ns).X1D,SM(ns).Z1Dmean); % x intersections
    
    if isempty(xz)
        %  no intersections found
        SM(ns).Xshoreline=NaN;
        SM(ns).BeachArea=NaN;
    else
        % shoreline x location
        SM(ns).Xshoreline=max(xz);
        subair=find(SM(ns).X1D > 0 & SM(ns).X1D <= max(xz));
        % estimate subarial "volume" (subarial profile area)
        ixd=find(~isnan(SM(ns).Z1Dmean));
        zi=interp1(SM(ns).X1D(ixd),SM(ns).Z1Dmean(ixd),SM(ns).X1D,'linear','extrap');
        SM(ns).BeachArea=nansum(zi(subair));
    end
    
    % find ShorelineElev-to-UpperShorefaceElev shore slope
    
    xl=[SM(ns).X1D(1) SM(ns).X1D(end)]; 
    sel=[ShorelineElev ShorelineElev];
    x1=intersections(xl,sel,SM(ns).X1D,SM(ns).Z1Dmean);
    if ~isempty(x1); x1=nanmin(x1); end
    sel=[UpperShorefaceElev UpperShorefaceElev];
    x2=intersections(xl,sel,SM(ns).X1D,SM(ns).Z1Dmean);
    if ~isempty(x2); x2=nanmin(x2); end
    
    if ~isempty(x1) && ~isempty(x2)
        SM(ns).ShorefaceSlope=(UpperShorefaceElev-ShorelineElev)/abs(x2-x1);
    else
        SM(ns).ShorefaceSlope=NaN;
        %  insufficient data for the shoreface slope calc
        SM(ns).Xshoreline=NaN;
        SM(ns).BeachArea=NaN;
    end

end

%--------------------------------------------------------
% sort the appended SM struct by Xshoreline location

  T=struct2table(SM); % sort by date before saving
  sortedT = sortrows(T, 'Xshoreline');
  %sortedT = sortrows(T, 'BeachArea'); % option to sort by area instead
  SM=table2struct(sortedT)';
  
  % find max nonNaN profile
  nsmax=find(~isnan([SM.Xshoreline]), 1, 'last' );

  %--------------------------------------------------------
  % make basic profile elevations 2d matrix from SM profiles
   Z=nan([nsmax xmax+1]); 
   for n=1:nsmax
       Z(n,:)=SM(n).Z1Dmean;
   end
   
   % round shoreline locations to nearest 1m in x
   xsl1m=round([SM.Xshoreline]);
   xslmin=round(nanmin([SM.Xshoreline]));
   xslmax=round(nanmax([SM.Xshoreline]));

   % Use Z to make a version of the SM struct array, SM1m, which averages
   % profile together in 1m xshore shoreline location increments

   m=0;
   for n=xslmin:xslmax % loop through 1m increment x shoreline location range
       m=m+1;
       SM1m(m).Xshoreline=n;
       SM1m(m).X1D=0:xmax;
       idx=find(xsl1m == n); % find profiles with same shorelin location 
       if ~isempty(idx)
        SM1m(m).Z1Dmean=nanmean(Z(idx,:),1); % avg them together
       else
        SM1m(m).Z1Dmean=nan*Z(1,:); % no profiles for this shoreline location
       end
   end

% new max number of profiles
nsmax=size(SM1m,2);

% grid the 1m shoreline increment profiles. This fills in any
%  data gaps owing to short xshore profiles and/or 1m shoreline location
%  increments with no observed profiles
  
[X,Y]=meshgrid(0:xmax,1:nsmax);
N=[];Xg=[];Zg=[];
for ns=1:nsmax
      X1=SM1m(ns).X1D;Z1=SM1m(ns).Z1Dmean;...
      X1(isnan(Z1))=[];Z1(isnan(Z1))=[];...
      Xg=[Xg X1];Zg=[Zg Z1];N=[N X1*0+ns];
end

Z=griddata(Xg,N,Zg,X,Y);


%  option to view the Z matrix
% figure;surf(X,Y,Z);shading flat;set(gca,'xdir','reverse');
% colormap(BeachColorMap);set(gca,'clim',[-8 5]);view(2)

% calculate shoreface slopes for each gridded profile
xl=[0 xmax];
sel=[ShorelineElev ShorelineElev];
usf=[UpperShorefaceElev UpperShorefaceElev]; 

for i=1:size(Z,1)
       xshore(i)=max(intersections(xl,sel,0:xmax,Z(i,:)));
       xusf(i)=max(intersections(xl,usf,0:xmax,Z(i,:)));
       xslope(i)=(UpperShorefaceElev-ShorelineElev)/abs(xusf(i)-xshore(i));
end


% smooth the profiles/Z matrix in the shoreline x location
%  direction base don the input SmoothWindow length scale
win=round(SmoothWindow);
if win < 1;win=1;end
C=movmean(Z,win,1); % running mean smoothing window 

% build characteristic profile array
xl=[0 xmax];
for m=1:size(C,1)  
  selx=max(intersections(xl,sel,0:xmax,C(m,:)));
  usfx=max(intersections(xl,usf,0:xmax,C(m,:)));
  CP(m).Xshoreline=selx;
  CP(m).ShorefaceSlope=(UpperShorefaceElev-ShorelineElev)/abs(usfx-selx);
  CP(m).X1D=0:xmax;
  CP(m).Z1Dmean=C(m,:);  
end


% make figure of data in charactersitic profile CP struct array

figure('position',[440 214 774 584]);

ax2=axes('position',[.25 .1 .5 .8]);
imagesc(0:xmax,min([CP.Xshoreline]):max([CP.Xshoreline]),C,'AlphaData',~isnan(C));
%colormap(BeachColorMap);
set(gca,'clim',[-8 5],'ydir','normal','xdir','reverse');
grid on;set(gca,'fontsize',13);
xlabel('Xshore Distance (m)');
title({'Characterstic Mean Xshore Profiles vs. Shoreline Location',...
    ['MOP: ' num2str(SM(1).Mopnum)], ['Shoreline Location (Y dimension) Smoothing Window: '...
    num2str(round(SmoothWindow)) ' m']});

hold on;
for ns=1:size(SM,2)
    idx=find(~isnan([SM(ns).Z1Dmean]));
    if ~isempty(idx)
    plot(0,SM(ns).Xshoreline,...
        'k.','markersize',10);
    scp=plot(SM(ns).X1D(idx),SM(ns).Xshoreline*ones(size(idx)),...
        'k-','linewidth',2);scp.Color(4)=.1;
    end
end
yl=get(ax2,'ylim');p1=plot([0 0],yl,'.-','color','k','markersize',10);
legend(p1,'surveys','location','southwest');

BeachBarColorbar

ax1=axes('position',[.06 .1 .15 .8]);
plot([CP.ShorefaceSlope],[CP.Xshoreline]);
ylabel('Beach Width = Mop MSL Shoreline X Location (m)');
xlabel({'MSL-MHHW';'Shoreface Slope'});
yl=get(ax2,'ylim');
set(gca,'ylim',yl,'fontsize',13);grid on;

end

