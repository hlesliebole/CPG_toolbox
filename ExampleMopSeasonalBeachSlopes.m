% Example script to calculate seasonal beach face slopes (MSL to MHW)
% for an input Mop number.  
%
% - slopes are estimated from mean monthly, quarterly, and global profile
%    shapes derived from all the historical survey data
% - slopes include a mean slope between MSL and MHW and lower/upper beach
%    face slopes between MSL and the MSL-MHW midpoint, and the MSL-MHW midpoint
%    and MHW, respectively.
%
%  - Reads survey data from the CPG MOP *SA.mat and generates profiles
%     directly from this 1m spatial averge survey data rather than using
%     the *SM.mat file profiles.
%  - Uses the function GetMeanNearestProfiles(SA) at the end of this file
%      to calculate the profiles using the nearest to transect method rather
%      than gridding.

clearvars
close all
%%-------------------------------
% Set CPG MOP path and Mop number
%%--------------------------------

addpath '/Volumes/group/MOPS'
addpath '/Volumes/group/MOPS/toolbox'

MopNum=507;

%% load the SA struct array for this mop
load(['M' num2str(MopNum,'%5.5i') 'SA.mat'],'SA');

%% Example option to reduce to dates beginning when there are truck surveys
% idx=find(contains({SA.Source},'Trk'));
% SA=SA(idx(1):end);

%% Get mean profile info for this Mop from its SA struct array using 
%    matlab function GetMeanNearestProfiles(SA)
%    
%   profiles are derived from the nearest survey points to the transect
%   projected on to the transect line.
%
%   X0BeachOnly = Mop cross-shore back beach point x offset of the returned 
%                 X1Dt locations which are based on the Truck LiDAR beach-only
%                 location on the mop transect.
%   X1dt = xshore distance relative to the truck beach-only boundary if
%          it exists for this profile, otherwise the Mop back beach point
%   Zf = size(SA,2) profiles of individual surveys 
%   Zm = 12 mean month profiles
%   Zq = 4 mean quarter profiles
%   Zs = 2 mean seasonal profiles
%   Zg = 1 global mean profile
[X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanNearestProfiles(SA);

%% Calculate global mean profile slopes
% get the MSL shoreline position from the global mean profile
gxMSL=intersections([X1Dt(1) X1Dt(end)],[0.774 0.774],X1Dt,Zg);
if ~isempty(gxMSL);gxMSL=gxMSL(end);end
% get the MHW shoreline position from the global mean profile
gxMHW=intersections([X1Dt(1) X1Dt(end)],[1.344 1.344],X1Dt,Zg);
if ~isempty(gxMHW);gxMHW=gxMHW(end);end
% halfway between msl and mhw
z=mean([0.774 1.344]);
gxHalf=intersections([X1Dt(1) X1Dt(end)],[z z],X1Dt,Zg);
if ~isempty(gxHalf);gxHalf=gxHalf(end);end
% mean MSL-MHW slope
gBeta=(1.344-0.744)/(gxMSL-gxMHW);
% upper beach face slope
gBetaUpper=(1.344-z)/(gxHalf-gxMHW);
% lower beach face slope
gBetaLower=(z-0.774)/(gxMSL-gxHalf);

%%  Calculate slopes from the 3-month running mean monthly profiles
for n=1:size(Zm3,1) % loop through 12 months
 
    % MSL location
    xint=intersections([X1Dt(1) X1Dt(end)],[0.774 0.774],X1Dt,Zm3(n,:));
    if ~isnan(xint);xMSL(n)=xint(end);else;xMSL(n)=NaN;end
     
    % MHW location
    xint=intersections([X1Dt(1) X1Dt(end)],[1.344 1.344],X1Dt,Zm3(n,:));
    if ~isnan(xint);xMHW(n)=xint(end);else;xMHW(n)=NaN;end

    % (MSL+MHW)/2 location
    z=mean([0.774 1.344]);
    xint=intersections([X1Dt(1) X1Dt(end)],[z z],X1Dt,Zm3(n,:));
    if ~isnan(xint);xHalf(n)=xint(end);else;xHalf(n)=NaN;end
   
    Beta(n)=(1.344-0.774)/(xMSL(n)-xMHW(n)); % mean beach face slope
    BetaUpper(n)=(1.344-z)/(xHalf(n)-xMHW(n)); % mean upper beach face slope
    BetaLower(n)=(z-0.774)/(xMSL(n)-xHalf(n)); % mean lower beach face slope

end

%%  Calculate quarterly profile slopes
for n=1:size(Zq,1) % loop through 12 months
 
    % MSL location
    xint=intersections([X1Dt(1) X1Dt(end)],[0.774 0.774],X1Dt,Zq(n,:));
    if ~isnan(xint);xMSL(n)=xint(end);else;xMSL(n)=NaN;end
     
    % MHW location
    xint=intersections([X1Dt(1) X1Dt(end)],[1.344 1.344],X1Dt,Zq(n,:));
    if ~isnan(xint);xMHW(n)=xint(end);else;xMHW(n)=NaN;end

    % (MSL+MHW)/2 location
    z=mean([0.774 1.344]);
    xint=intersections([X1Dt(1) X1Dt(end)],[z z],X1Dt,Zq(n,:));
    if ~isnan(xint);xHalf(n)=xint(end);else;xHalf(n)=NaN;end
   
    BetaQ(n)=(1.344-0.774)/(xMSL(n)-xMHW(n)); % mean beach face slope
    BetaUpperQ(n)=(1.344-z)/(xHalf(n)-xMHW(n)); % mean upper beach face slope
    BetaLowerQ(n)=(z-0.774)/(xMSL(n)-xHalf(n)); % mean lower beach face slope

end

%% Make illustrative figure of global beach face profile and slopes
figure('position',[6   310   984   475]);
p(1)=plot(X1Dt,Zg,'-','color',[.8 .8 .8],'linewidth',4,'DisplayName','Global Mean Beach Face');hold on;grid on
set(gca,'xdir','reverse','ylim',[.4 1.8],'fontsize',14,'linewidth',2)
xlabel('Cross-shore Distance (m)');ylabel('Elevation (m, navd88)');
title(['Mop ' num2str(MopNum) ' Global Mean Beach Face Profile'],'fontsize',18);
xl=get(gca,'xlim');
plot(xl,[0.774 0.774],'k--','linewidth',2);text(xl(2)*.98,0.85,'MSL','fontsize',14)
plot(xl,[1.344 1.344],'k--','linewidth',2);text(xl(2)*.98,1.4,'MHW','fontsize',14)
plot(xl,[z z],'k:','linewidth',2);text(xl(2)*.98,z+.06,'(MSL+MHW) / 2','fontsize',14)

p(3)=plot([gxMSL gxHalf],[0.774 z],'-','color',[0 .8 0],'linewidth',2,'DisplayName',['Mean Upper Beach Face Slope {\beta_{Lower}}=' num2str(gBetaLower,'%5.3f')]);
p(4)=plot([gxHalf gxMHW],[z 1.344],'r-','linewidth',3,'DisplayName',['Mean Upper Beach Face Slope {\beta_{Upper}}=' num2str(gBetaUpper,'%5.3f')]);
p(2)=plot([gxMSL gxMHW],[0.774 1.344],'b-','linewidth',2,'DisplayName',['Mean Beach Face Slope {\beta}=' num2str(gBeta,'%5.3f')]);
plot([gxMSL gxHalf gxMHW],[0.774 z 1.344],'k.','markersize',20)
legend(p,'location','north');
makepng('GlobalBeachFaceProfile.png')

%% Make figure of monthly mean profiles

col=jet(12);
charmon = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}';
figure('position',[244   237   830   514]);
for n=1:12
    pm(n)=plot(X1Dt,Zm(n,:),'-','color',col(n,:),'linewidth',2,'DisplayName',charmon{n});hold on
end
grid on;
set(gca,'xdir','reverse','ylim',[.4 1.8],'fontsize',14,'linewidth',2)
xlabel('Cross-shore Distance (m)');ylabel('Elevation (m, navd88)');
xl=get(gca,'xlim');
plot(xl,[0.774 0.774],'k--','linewidth',2);text(xl(2)*.98,0.85,'MSL','fontsize',14)
plot(xl,[1.344 1.344],'k--','linewidth',2);text(xl(2)*.98,1.4,'MHW','fontsize',14)
plot(xl,[z z],'k:','linewidth',2);text(xl(2)*.98,z+.06,'(MSL+MHW) / 2','fontsize',14)
legend(pm,'location','eastoutside')
title(['Mop ' num2str(MopNum) ' Monthly Mean Beach Face Profile'],'fontsize',18);
makepng('MonthlyBeachFaceProfiles.png')


%% Make figure of 3-month running mean monthly profiles

col=jet(12);
char3mon = {'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ'}';
figure('position',[321   174   830   514]);
for n=1:12
    pm(n)=plot(X1Dt,Zm3(n,:),'-','color',col(n,:),'linewidth',2,'DisplayName',char3mon{n});hold on
end
grid on;
set(gca,'xdir','reverse','ylim',[.4 1.8],'fontsize',14,'linewidth',2)
xlabel('Cross-shore Distance (m)');ylabel('Elevation (m, navd88)');
xl=get(gca,'xlim');
plot(xl,[0.774 0.774],'k--','linewidth',2);text(xl(2)*.98,0.85,'MSL','fontsize',14)
plot(xl,[1.344 1.344],'k--','linewidth',2);text(xl(2)*.98,1.4,'MHW','fontsize',14)
plot(xl,[z z],'k:','linewidth',2);text(xl(2)*.98,z+.06,'(MSL+MHW) / 2','fontsize',14)
legend(pm,'location','eastoutside')
title(['Mop ' num2str(MopNum) ' Running 3-month Monthly Mean Beach Face Profiles'],'fontsize',18);
makepng('Running3MonthBeachFaceProfiles.png')

%% Make figure of quarterly mean profiles

col=jet(4);
charqtr = {'Jan-Mar','Apr-Jun','Jul-Sep','Oct-Dec'}';
figure('position',[385   129   830   514]);
for n=1:4
    pq(n)=plot(X1Dt,Zq(n,:),'-','color',col(n,:),'linewidth',2,'DisplayName',charqtr{n});hold on
end
grid on;
set(gca,'xdir','reverse','ylim',[.4 1.8],'fontsize',14,'linewidth',2)
xlabel('Cross-shore Distance (m)');ylabel('Elevation (m, navd88)');
xl=get(gca,'xlim');
plot(xl,[0.774 0.774],'k--','linewidth',2);text(xl(2)*.98,0.85,'MSL','fontsize',14)
plot(xl,[1.344 1.344],'k--','linewidth',2);text(xl(2)*.98,1.4,'MHW','fontsize',14)
plot(xl,[z z],'k:','linewidth',2);text(xl(2)*.98,z+.06,'(MSL+MHW) / 2','fontsize',14)
legend(pq,'location','eastoutside')
title(['Mop ' num2str(MopNum) ' Quarterly Mean Beach Face Profiles'],'fontsize',18);
makepng('QuarterlyBeachProfiles.png')

%%  Final Figure of the various slopes
figure('position',[ 460    25   862   546]);hold on;
set(gca,'xlim',[0.5 12.5],'fontsize',14,'linewidth',2);grid on
xlabel('12 Overlapping 3-month Periods');ylabel('Slope');box on
xl=get(gca,'xlim');
plot(xl,[gBeta gBeta],'b--','linewidth',2,'DisplayName','Global {\beta}');
plot(xl,[gBetaLower gBetaLower],'--','color',[0 .8 0],'linewidth',2,'DisplayName','Global {\beta_{Lower}}');
plot(xl,[gBetaUpper gBetaUpper],'r--','linewidth',2,'DisplayName','Global {\beta_{Upper}}');

plot(2.:3:11.,BetaQ,'bo','markersize',10,'linewidth',2,'DisplayName','Quarterly {\beta}');
plot(2.:3:11.,BetaLowerQ,'o','color',[0 .8 0],'markersize',10,'linewidth',2,'DisplayName','Quarterly {\beta_{Lower}}');
plot(2.:3:11.,BetaUpperQ,'ro','markersize',10,'linewidth',2,'DisplayName','Quarterly {\beta_{Upper}}');

plot(1:12,Beta,'b.-','markersize',20,'linewidth',2,'DisplayName','Monthly {\beta}');
plot(1:12,BetaLower,'.-','color',[0 .8 0],'markersize',20,'linewidth',2,'DisplayName','Monthly {\beta_{Lower}}');
plot(1:12,BetaUpper,'r.-','markersize',20,'linewidth',2,'DisplayName','Monthly {\beta_{Upper}}');

% Fit a sine wave to Beta 
SineParams=sineFit(1:12,Beta,0);
% Syntax:
%       [SineParams]=sineFit(x,y,optional)
%       Input: x and y values, y=offs+amp+sin(2*pi*f*x+phi)+noise
%              optional: plot graphics if ommited. Do not plot if 0.
%       Output: SineParams(1): offset (offs)
%               SineParams(2): amplitude (amp)
%               SineParams(3): frequency (f)
%               SineParams(4): phaseshift (phi)
%               SineParams(5): MSE , if negative then SineParams are from FFT 
%       yOut=offs+amp*sin(2*pi*f*x+phi)
%
BetaMean=SineParams(1);
BetaAmp=SineParams(2);
BetaFreq=SineParams(3);
BetaPhase=SineParams(4);
BetaMonthly=BetaMean+BetaAmp*sin(2*pi*BetaFreq*(1:12)+BetaPhase);
plot(1:12,BetaMonthly,'b.:','markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta}');

% Fit a sine wave to BetaLower 
SP=sineFit(1:12,BetaLower,0);
BetaLowerMonthly=SP(1)+SP(2)*sin(2*pi*SP(3)*(1:12)+SP(4));
plot(1:12,BetaLowerMonthly,'.:','color',[0 .8 0],'markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta_{Lower}}');

% Fit a sine wave to BetaUpper
SP=sineFit(1:12,BetaUpper,0);
BetaUpperMonthly=SP(1)+SP(2)*sin(2*pi*SP(3)*(1:12)+SP(4));
plot(1:12,BetaUpperMonthly,'r.:','markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta_{Upper}}');

set(gca,'xtick',1:12,'xticklabels',char3mon)
title(['Mop ' num2str(MopNum) ' Mean Profile Beach Face Slopes'],'fontsize',18);
legend('location','eastoutside')
makepng('SeasonalBeachFaceSlopes.png')


%% ---------------------------------------------------------------------
function [X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanNearestProfiles(SA)

%% returns mean profiles for the input SA struct array
%
%  X0BeachOnly = Mop xshore location (offset) of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to X0BeachOnly
%  Zf(M,N) = M survey profiles relative to X0BeachOnly
%  Zm(12,N) = monthly mean profiles (calculated by year first and then across all years)
%  Zm3(12,N) = running 3-month monthly mean profiles (calculated by year first and then across all years)             
%  Zq(4,N) = quarterly mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  Zs(2,N) = seasonal mean profiles (Jan-Jun,Jul-Dec)
%  Zg(1,N) = global mean profile

%% get nearest point profile
Ytol=50; % 50m alongcoast nearest point tolerance
Xtol=10; % 5m cross shore gap interpolation tolerance

[X1D,Z1D]=GetNearestPointsProfiles(SA,Ytol,Xtol);

%X2D=repmat(X1D,size(Z1D,1),1);
%% find the truck back beach boundary along the mop transect
trk=find(strcmp({SA.Source},'Trk'));
if ~isempty(trk)
for n=1:numel(trk)
    xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
end
else
    xback=0;
end
X0BeachOnly=median(xback);
% fprintf('Truck Back Beach x= %6.1f \n',X0BeachOnly)

%% reduce profiles to common back beach x location based on truck beach-only data
idx=find(X1D >= X0BeachOnly);
Z1Dt=Z1D(:,idx);
X1Dt=X1D(idx);

%% remove any elevations seaward of the min profile elevation
%  (lidar swash filter)
Zf=Z1Dt;
[zmin,imin]=min(Zf');
for n=1:size(Zf,1)
    if imin(n) > 0
        Zf(n,imin(n)+1:end)=NaN;
    end
end

% remove any additional outliers
TF=isoutlier(Zf,"mean");
Zf(TF == 1)=NaN;

% survey datetimes
Zdt=datetime([SA.Datenum],'convertfrom','datenum');

% year-month means
n=0;
for y=year(Zdt(1)):year(Zdt(end))
   for m=1:12
       n=n+1;
       Zyear(n)=y;
       Zmonth(n)=m;
       Zym(n,:)=Zf(1,:)*NaN;
       idx=find(year(Zdt) == y & month(Zdt) == m);
       if numel(idx) == 1
           Zym(n,:)=Zf(idx,:);
       elseif numel(idx) > 1
           Zym(n,:)=mean(Zf(idx,:),'omitnan');
       end
   end
end

% monthly means across all years
n=0;
for m=1:12
       n=n+1;
       Zm(n,:)=Zym(1,:)*NaN;
       idx=find(Zmonth == m);
       if numel(idx) == 1
           Zm(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
           Zm(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% 3-month running means across all years
n=0;
for m3=1:12
       n=n+1;
       Zm3(n,:)=Zym(1,:)*NaN;
       if m3 == 1
           m=[1 2 12];
       elseif m3 == 12
           m=[1 11 12];
       else
           m=m3-1:m3+1;
       end

       idx=find(Zmonth == m(1) | Zmonth == m(2) | Zmonth == m(3));
       if numel(idx) == 1
           Zm3(n,:)=Zym(idx,:);
       elseif numel(idx) > 1
           Zm3(n,:)=mean(Zym(idx,:),'omitnan');
       end
end


% quarterly means from monthly means
n=0;
for q=1:4
       n=n+1;
       Zq(n,:)=Zm(1,:)*NaN;      
       idx=(q-1)*3+1:q*3;
       Zq(n,:)=mean(Zm(idx,:),'omitnan');
end

% seasonal means from quarterly means
Zs(1,:)=mean(Zq(1:2,:),'omitnan');
Zs(2,:)=mean(Zq(3:4,:),'omitnan');
   
% global mean
Zg=mean(Zs,'omitnan');

end

