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
fid=fopen('SlopeSineFits.dat','w');
fprintf('              Beta                      Beta Lower                  Beta Upper\n')
fprintf('Mop  Mean    Amp   Freq    Phase  Mean    Amp    Freq   Phase  Mean    Amp    Freq   Phase   \n')
for MopNum=520:915

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


%% Fit a sine wave to non-NaN Betas 
SP=zeros(1,5);
SPL=zeros(1,5);
SPH=zeros(1,5);

x=1:12;
if sum(~isnan(Beta)) > 0

SP=sineFit(x(~isnan(Beta)),Beta(~isnan(Beta)),0);
BetaMean=SP(1);
BetaAmp=SP(2);
BetaFreq=SP(3);
BetaPhase=SP(4);
BetaMonthly=BetaMean+BetaAmp*sin(2*pi*BetaFreq*(1:12)+BetaPhase);
%plot(1:12,BetaMonthly,'b.:','markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta}');

% Fit a sine wave to BetaLower 
SPL=sineFit(x(~isnan(BetaLower)),BetaLower(~isnan(BetaLower)),0);
BetaLowerMonthly=SPL(1)+SPL(2)*sin(2*pi*SPL(3)*(1:12)+SPL(4));
%plot(1:12,BetaLowerMonthly,'.:','color',[0 .8 0],'markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta_{Lower}}');

% Fit a sine wave to BetaUpper
SPH=sineFit(x(~isnan(BetaUpper)),BetaUpper(~isnan(BetaUpper)),0);
BetaUpperMonthly=SPH(1)+SPH(2)*sin(2*pi*SPH(3)*(1:12)+SPH(4));
%plot(1:12,BetaUpperMonthly,'r.:','markersize',15,'linewidth',2,'DisplayName','Sine Fit {\beta_{Upper}}');

end

fprintf('%i %5.3f  %5.4f  %5.3f  %5.2f | %4.3f  %4.4f  %4.3f  %4.2f | %4.3f  %4.4f  %4.3f  %4.2f \n',...
    MopNum,SP(1:4),SPL(1:4),SPH(1:4))
fprintf(fid,'%i %5.3f  %5.4f  %5.3f  %5.2f %4.3f  %4.4f  %4.3f  %4.2f %4.3f  %4.4f  %4.3f  %4.2f \n',...
    MopNum,SP(1:4),SPL(1:4),SPH(1:4));
end

fclose(fid);
% set(gca,'xtick',1:12,'xticklabels',char3mon)
% title(['Mop ' num2str(MopNum) ' Mean Profile Beach Face Slopes'],'fontsize',18);
% legend('location','eastoutside')
% makepng('SeasonalBeachFaceSlopes.png')


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
    mdx=find(~isnan(Z1D(trk(n),:)), 1 );
    if ~isempty(mdx)
    xback(n)=X1D(find(~isnan(Z1D(trk(n),:)), 1 ));
    else
    xback(n)=0;
    end
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

