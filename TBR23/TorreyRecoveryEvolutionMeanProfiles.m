clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

% get global mean profile info for this reach
GMP=GetGMPmopRange(580,589);

% figure out how many surveys to process starting 6 Apr
MopNumber=580;
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));
nsurveys=numel(jumbo)-jstart;
col=jet(nsurveys);

% loop throu TP mops 578-595
%close all
figure('position',[1          59        1020         738]);
hold on
m=0;

for mm=1:nsurveys+1%nsurveys+1
m=m+1;
z1t=[];
z2t=[]; 
zt=[];
xlag=0;
j=0; % mop counter
for MopNumber=580:589
 j=j+1;

% find matching global mop profile index for this mop
gdx=find([GMP.Mop] == MopNumber);
% progressively adjust xshore reference frame relative to starting Mop
xlag=xlag+GMP(gdx).Xlag;
    
load(['M00' num2str(MopNumber,'%3.3i') 'SA.mat'])
jumbo=find(contains({SA.File},'umbo') | contains({SA.File},'etski'));
% everything relative to 6 Apr jumbo
jstart=find([SA(jumbo).Datenum] == datenum(2023,4,6));

% get survey index realtive to 6 Apr
n=jumbo(jstart+m-1);

% use nearest points to transect profile estimation method
[x1d,z1]=GetNonGriddedProfile(MopNumber,n);

% shift xshore reference frame relative to first mop
x1d=x1d+xlag;

% find overlap with first mop reference frame and add to profile
%  matrix
 if j == 1
    xt=x1d;
    zt(j,:)=z1;
 else
    zl=xt*NaN;
    [idx,ndx]=ismember(x1d,xt);
    zl(ndx(ndx > 0))=z1(idx == 1);
    zt(j,:)=zl;
 end

end

% plot mean profile across mops
plot(xt,mean(zt,'omitnan'),'-','linewidth',2,'DisplayName',datestr(SA(n).Datenum));
hold on;

end

grid on
box on
set(gca,'fontsize',14,'linewidth',2)
set(gca,'xdir','reverse')
legend('location','northwest')
xlabel('Cross-shore Distance (m)')
ylabel('Elevation (m, NAVD88)')
title('Alongshore Averaged Xshore Profiles: Mops 580-589')
makepng('TorreyAlongshoreAvgProfiles.png')
