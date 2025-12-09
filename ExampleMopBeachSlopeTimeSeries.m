% Example script to calculate Mopbeach face slopes (MSL to MHW)
% for an input Mop range.  
%
% -  Only plots the MSL-MHW slope but code includes the mean slope 
%    between MSL and MHW and lower/upper beach face slopes 
%     between MSL and the MSL-MHW midpoint, and the MSL-MHW midpoint
%     and MHW, respectively.  
%
%  - Reads survey data from the CPG MOP *SA.mat and generates profiles
%     directly from this 1m spatial averge survey data rather than using
%     the *SM.mat file profiles.
%
%  - Uses the function GetMeanNearestProfiles(SA) at the end of this file
%      to calculate the profiles using the nearest to transect method rather
%      than gridding.

clearvars
close all
%%-------------------------------
% Set CPG MOP path and Mop numbers
%%--------------------------------

addpath '/Volumes/group/MOPS'
addpath '/Volumes/group/MOPS/toolbox'

MopStart=580;
MopEnd=585;

%% ----------------------------------------------------
figure('position',[ 71         377        1092         346]);hold on;
hold on;

mm=0;
col=jet(numel(MopStart:MopEnd));
for MopNum=MopStart:MopEnd
    fprintf('Getting Mop  %i\n',MopNum)
    mm=mm+1;
    Beta=[];

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

%%  Calculate slopes from individual profiles
for n=1:size(Zf,1) % loop through all surveys
 
    % MSL location
    xint=intersections([X1Dt(1) X1Dt(end)],[0.774 0.774],X1Dt,Zf(n,:));
    if ~isnan(xint);xMSL(n)=xint(end);else;xMSL(n)=NaN;end
     
    % MHW location
    xint=intersections([X1Dt(1) X1Dt(end)],[1.344 1.344],X1Dt,Zf(n,:));
    if ~isnan(xint);xMHW(n)=xint(end);else;xMHW(n)=NaN;end

    % (MSL+MHW)/2 location
    z=mean([0.774 1.344]);
    xint=intersections([X1Dt(1) X1Dt(end)],[z z],X1Dt,Zf(n,:));
    if ~isnan(xint);xHalf(n)=xint(end);else;xHalf(n)=NaN;end
   
    Beta(n)=(1.344-0.774)/(xMSL(n)-xMHW(n)); % mean beach face slope
    BetaUpper(n)=(1.344-z)/(xHalf(n)-xMHW(n)); % mean upper beach face slope
    BetaLower(n)=(z-0.774)/(xMSL(n)-xHalf(n)); % mean lower beach face slope

end


%%  Final Figure of the various slopes
x=datetime([SA.Datenum],'convertfrom','datenum');
%plot(x,BetaUpper,'.','markersize',20,'color',col(mm,:),'linewidth',2,'DisplayName',['Mop ' num2str(MopNum)]);
plot(x,Beta,'ko','markersize',4,'markerfacecolor',col(mm,:),'linewidth',1,'DisplayName',['Mop ' num2str(MopNum)]);

end
%set(gca,'xlim',[0.5 12.5],'fontsize',14,'linewidth',2);grid on
title('SIO MOPs')
xlabel('Date');ylabel('MSL-MHW Slope');box on
xl=get(gca,'xlim');

%set(gca,'fontsize',14,'xlim',[datetime(2022,10,15) datetime(2023,6,15)])
grid on;
legend('location','eastoutside')


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

