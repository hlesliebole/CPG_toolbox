close all
clearvars

% want to maximize the number of months that have
%  at least one "complete subaerial profile" 
%  while minimizing the amount of profile 
%  interpolation/extrapolation required to
%  complete partial profiles.
%
%  Ideally complete profiles start at the truck
%   beach only boundary but a more seaward point
%   may result in sigficantly more viable profiles
%   (eg. blacks beach)
%
%  complete profile qc/construction variables in order of implementation:
%
%  1. min profile elevation (eg. msl or mhw)
%  2. complete profiles origin (relative to truck beach only boundary)
%  3. max profile gap interpolation/extrapolation tolerance (m)
%  4. max percent missing data tolerance
%  5. max gap to the back beach origin tolerance

% first find all the complete profiles down to the min profile
%  elevation.

% loop though incomplete profiles and complete them while
%  keeping track of percent missing data, any gap in
%  meters to the truck origin, and the largest gap filled
%  in meters.
clearvars
close all

load M00654SA.mat;
nsurv=size(SA,2);
fprintf('%i surveys.\n',nsurv)

[X0BeachOnly,X1Dt,Zf,Zm,Zm3,Zq,Zs,Zg]=GetMeanNearestProfiles(SA);

%% Idetify differnt types of surveys
utair=find(contains({SA.Source},'UTAir'));
trk=find(contains({SA.Source},'Trk'));
usace=find(contains({SA.Source},'USACE'));
usgs=find(contains({SA.Source},'USGS'));
ccc=find(contains({SA.Source},'CCC'));
jumbo=find(contains({SA.File},'umbo'));
wheel=find(contains({SA.Source},'iG8wheel'));

Zdt=datetime([SA.Datenum],'convertfrom','datenum');
Nmonths=numel(unique(year(Zdt)*100+month(Zdt)));


zc=0.774; %msl
dztol=0.1;
[zmin,imin]=min(abs(Zf-zc)');
imin(isnan(zmin))=NaN;
imin(zmin > dztol)=NaN;

figure
subplot(2,1,1);plot(zmin);hold on;plot([0 numel(zmin)],[dztol dztol],'k--');
set(gca,'yscale','log');
xlabel('Survey Number');ylabel('Dist from contour (m)');
title(['Min Profile Point Distance from ' num2str(zc) ' m contour']);
subplot(2,1,2);plot(imin(~isnan(imin)));hold on;plot(X1Dt(imin(~isnan(imin))),'s')
fprintf('%i surveys have a valid %6.3f m contour location data.\n',numel(find(~isnan(imin))),zc)
xlabel('Survey Number');ylabel('Dist from back beach (m)');
title(['Contour Distance from (Truck) back beach (m)' num2str(zc) ' m contour']);

% TF=isoutlier(Z);
% Zf=Z;Zf(TF == 1)=NaN;
% Zmin=min(Zf);

% figure('position',[136         156        1227         582]);
% plot(X2D(:),Z1D(:),'k.',X2Df(:),Zf(:),'r.')
% hold on;plot(X1D,Z1D(end,:),'g.');
% zmed=median(Zf,'omitnan');
% plot(X1Dt,zmed,'y-','linewidth',2);
% title(SA(end).File)
% plot([X0BeachOnly X0BeachOnly],[0 3],'k-');
% set(gca,'xdir','reverse')

figure;
hold on;
MonthColormap;
for m=1:12
    plot(X1Dt,Zm(m,:),'-','color',col(m,:)/255,'linewidth',2)
    hold on;
end
legend

figure;
hold on;
for m=1:4
    plot(X1Dt,Zq(m,:),'-','linewidth',2)
    hold on;
end
legend

% seasonal means from quarterly means
Zs(1,:)=mean(Zq(1:2,:),'omitnan');
Zs(2,:)=mean(Zq(3:4,:),'omitnan');

figure;
hold on;
    plot(X1Dt,Zs(1,:),'r-','linewidth',2)
    plot(X1Dt,Zs(2,:),'g-','linewidth',2)
    plot(X1Dt,Zg,'k-','linewidth',2)
    hold on;
legend
    
zmin=min(Zf');zmax=max(Zf');
zmin(zmin < -2)=-2;
figure;hold on;
for i=1:numel(zmin)
    plot([Zdt(i) Zdt(i)],[zmin(i) zmax(i)],'k.-')
    plot([Zdt(i)],[zmax(i)],'r.')
end
plot([Zdt(1) Zdt(end)],[1.344 1.344],'-','color',[.7 .7 .7])
plot([Zdt(1) Zdt(end)],[0.774 0.774],'-','color',[.7 .7 .7])
plot([Zdt(1) Zdt(end)],[0.218 0.218],'-','color',[.7 .7 .7])
set(gca,'ylim',[-1 5])

% use mhw as min elev
%  loop through transects find number that are complete down
%  to min elev
Zmin=0.774;%1.344;
for Zmin=[0.218 0.774 1.344]
fprintf(' ------ Zmin = %6.3f\n',Zmin)

m=0;
for n=1:size(Zf,1)
    % locate min elevation crossing if it exists
    idx=find( Zf(n,:) <= Zmin ,1,'first');
     if ~isempty(idx)
      GoodFrac(n)=sum(~isnan(Zf(n,1:idx)))/idx;
     else
      
      idx=find( Zf(n,:) > Zmin ,1,'last');
        if ~isempty(idx)
          GoodFrac(n)=-sum(~isnan(Zf(n,1:idx)))/idx;
        else
          GoodFrac(n)=0;  
        end
     end
end

NCmonths=numel(unique(year(Zdt(GoodFrac == 1))*100+month(Zdt(GoodFrac == 1))));
N90months=numel(unique(year(Zdt(abs(GoodFrac) >= .9))*100+month(Zdt(abs(GoodFrac) >= .9))));
N80months=numel(unique(year(Zdt(abs(GoodFrac) >= .8))*100+month(Zdt(abs(GoodFrac) >= .8))));
N70months=numel(unique(year(Zdt(abs(GoodFrac) >= .7))*100+month(Zdt(abs(GoodFrac) >= .7))));
N60months=numel(unique(year(Zdt(abs(GoodFrac) >= .6))*100+month(Zdt(abs(GoodFrac) >= .6))));
fprintf('Num Surveys: %i  - Num unique months: %i\n',size(Zf,1),Nmonths);
fprintf('Complete Profiles: %i - Num unique months: %i\n',numel(find(GoodFrac == 1)),NCmonths);
fprintf('90%%+ Complete Profiles: %i - Num unique months: %i\n',numel(find(abs(GoodFrac) >= .9)),N90months);
fprintf('80%%+ Complete Profiles: %i - Num unique months: %i\n',numel(find(abs(GoodFrac) >= .8)),N80months);
fprintf('70%%+ Complete Profiles: %i - Num unique months: %i\n',numel(find(abs(GoodFrac) >= .7)),N70months);
fprintf('60%%+ Complete Profiles: %i - Num unique months: %i\n',numel(find(abs(GoodFrac) >= .6)),N60months);
end
figure;plot(GoodFrac,'*')

Zt=Zf;Zt(Zt < 1.344)=NaN;
figure('position',[266    29   552   720]);
imagesc(~isnan(Zt))
set(gca,'xlim',[1 100])

 idx=find(~isnan(Zt(2,:))); yq=interp1(idx,Zt(2,idx),1:100,'linear','extrap');

