% Code to view survey SA data behind shoreline position (beach width) 
% estimates.


addpath /Users/William/Desktop/MOPS/toolbox/profiles

%dt=ncread('MOP581_complete.nc','time');
dt=ncread('MOP581_complete.nc','time');
xmhw=ncread('MOP581_complete.nc','xmhw');
dt=dt+datenum(1970,1,1);
t=datetime(dt,'convertfrom','datenum');
%%
% load the desired SA struct array
load M00581SA.mat

% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=10; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=5; % patch any cross-shore profile gaps < 5m wide
YdistTol=25; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

% get some beach widths
Znavd88=1.344; % MHW
[BWmin,BWmax]=GetCpgProfileBeachWidths(Znavd88,X1Dmop,Z1Dmedian); 

%%
figure
plot(Zdatetime,BWmin-mean(BWmin,'omitnan'),'.-');hold on
plot(t,xmhw,'.-');

%%
% figure
% XbackgapTol=10;
% [BVmin,BVmax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dmop,Z1Dmedian);
% plot(Zdatetime,BVmin,'.-',Zdatetime,BVmax,'o')

%%

load M00581SA.mat
% first find all the official jumbo survey data indexes
jdx=find(contains({SA.File},'umbo'));

% now loop through them to see if there are any Trk or AtvMR surveys within 2 %  days
for Surv=jdx 
 ldx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
 ([SA.Datenum] >= SA(Surv).Datenum-2 & [SA.Datenum] <= SA(Surv).Datenum+2));

%  if yes, add any LiDAR x,y,z data to the jumbo data fields
 if ~isempty(ldx)
  for SameSurv=ldx % loop through same time lidar surveys
      SA(Surv).X=vertcat(SA(SameSurv).X,SA(Surv).X);
      SA(Surv).Y=vertcat(SA(SameSurv).Y,SA(Surv).Y);
      SA(Surv).Z=vertcat(SA(SameSurv).Z,SA(Surv).Z);
  end
 end
end

% now reduce to just the (LiDAR enhanced) jumbos
SA=SA(jdx);

% now pass to the profile function using a modest YdistTol that should
%  pull in most jumbo data that is near the transect
NumSubTrans=1; % 1 will just return the main transect info; 100 would make
                % statistical profiles with complete use of any 1m res LiDAR
XgapTol=25; % patch any cross-shore profile gaps < 25m wide for jumbos
YdistTol=50; % can go as large as 50m up- downcoast for a single Mop area

% get nearest point profiles

[X1Dmop,X1Dcpg,Zdatetime,Z1Dtrans,Z1Dmean,Z1Dmedian,Z1Dmin,Z1Dmax,Z1Dstd]=...
  GetCpgNearestPointProfiles(SA,NumSubTrans,XgapTol,YdistTol);
% (Note: It can take a few minutes to run depending on the number of 
%    NumSubTrans subtransects being considered.)

figure
Znavd88=-7;
XbackgapTol=10;
[BVmin,BVmax]=GetCpgProfileBeachVolumes(Znavd88,XbackgapTol,X1Dmop,Z1Dtrans);
plot(Zdatetime,BVmin,'.-',Zdatetime,BVmax,'o')

%
imean=find(Zdatetime > datetime(2007,1,1));
figure
Znavd88=-7;
XbackgapTol=10;
[BVmin,BVmax]=GetCpgProfileMobileVolumes(Znavd88,XbackgapTol,X1Dmop,Z1Dtrans);
plot(Zdatetime,BVmin-mean(BVmin(imean),'omitnan'),'.-')

hold on

Znavd88=-4;
XbackgapTol=10;
[BVmin,BVmax]=GetCpgProfileMobileVolumes(Znavd88,XbackgapTol,X1Dmop,Z1Dtrans);
plot(Zdatetime,BVmin-mean(BVmin(imean),'omitnan'),'.-')

Znavd88=.774;
XbackgapTol=10;
[BVmin,BVmax]=GetCpgProfileMobileVolumes(Znavd88,XbackgapTol,X1Dmop,Z1Dtrans);
plot(Zdatetime,BVmin-mean(BVmin(imean),'omitnan'),'.-')

%%
zc=Z1Dtrans(:,1:522);
xc=X1Dmop(1:522);

figure;
Znavd88=-6;
XbackgapTol=10;
[BVmin,BVmax]=GetCpgProfileMobileVolumes(Znavd88,XbackgapTol,xc,zc);
plot(Zdatetime,BVmin-mean(BVmin(imean),'omitnan'),'.-')
