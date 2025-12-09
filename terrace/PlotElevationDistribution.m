%  Calculates the most recent changes to an alongshore
%  range mop profile elevations as a function of their current 
%  elevations.  A form of alongshore averaging when the number
%  of mops is > 1.

clearvars
close all
addpath /Users/William/Desktop/Mops
addpath /Users/William/Desktop/Mops/toolbox
%DefineMopPath

%% settings
MopStart=507; % startting mop
MopEnd=512; % ending mop
MaxDistFrmTransect=10; % max distance (m) from transect tolerance for profile estimation
MaxXshoreGapTolerance=5; % max xshore gap in profile that is interpolated, otherwise =NaN

load Mops507to514GBP.mat

all_marks = {'o','+','*','x','s','d','^','v','>','<','p','h'};
figure('position',[100 31 1076 766]);

%nx=numel(-20:90);
%Z=NaN(nx,8);

% initialize vectors that will hold all the profile elevations
z1da=[];
z1dap=[];

m=0; % mop counter
% loop throu SIO mops 507-511
for MopNumber=507:512
m=m+1;

% load the 1m averge point data SA struct array for this mop
load(['M00' num2str(MopNumber,'%2.2i') 'SA.mat'])

% most recent (last) profile index
n=size(SA,2);

% get all the nearest point profiles
%[X1d,Z1d]=GetNearestProfiles(SA,MaxDistFrmTransect,MaxXshoreGapTolerance);
% get most recent nongridded (nearest points) profile

[x1d,z1di]=GetNonGriddedProfile(MopNumber,n);
% get most recent nongridded (nearest points) profile
[x1p,z1dp]=GetNonGriddedProfile(MopNumber,n-1);
% combine elevations of all the mop profiles
z1da=[z1da z1di];
z1dap=[z1dap z1dp];


% x1d=x1d+GBP(m).Xlag;
% nidx=x1d+21;
% x1d=x1d(nidx > 0 & nidx <= nx);
% z1di=z1di(nidx > 0 & nidx <= nx);
% nidx=nidx(nidx > 0 & nidx <= nx);
% 
% Z(nidx,m)=z1di;

%  p(m)=plot(x1d,z1di-z1dp,'LineStyle','-',...
%         'Marker',all_marks{mod(1+round(11*rand),13)},...
%         'DisplayName',num2str(GBP(m).Mop),...
%         'linewidth',2);hold on;
    p(m)=plot(z1dp,z1di-z1dp,...
         all_marks{mod(1+round(10*rand),12)},...
        'linewidth',2,...
        'DisplayName',num2str(GBP(m).Mop));hold on;
 hold on
end

% hold on;pm=plot(x1d,mean(Z(:,1:5)','omitnan'),'k-','linewidth',4)
%yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
 
grid on;
set(gca,'fontsize',12);xlabel('Previous profile elevation');ylabel('Elevation Change (m, NAVD88)')
%set(gca,'xdir','reverse','fontsize',14);

title(['Surveys on ' datestr(SA(n).Datenum)],...
   'fontsize',16);

yl=get(gca,'ylim');
%set(gca,'xlim',[-20 90]);xl=get(gca,'xlim');
plot([2.119 2.119],yl,'k--');text(2.26,.95*yl(2),' HAT','fontsize',14,'horizontalalign','right');
plot([1.566 1.566],yl,'k--');text(1.7,.95*yl(2),' MHHW','fontsize',14,'horizontalalign','right');
plot([1.344 1.344],yl,'k--');text(1.44,.95*yl(2),' MHW','fontsize',14,'horizontalalign','right');
plot([.774 .774],yl,'k--');text(.9,.95*yl(2),' MSL','fontsize',14,'horizontalalign','right');
plot([0.218 0.218],yl,'k--');text(.34,.95*yl(2),' MLW','fontsize',14,'horizontalalign','right');
%plot([-0.058 -0.058],yl,'k--');text(0.05,yl(2),' MLLW','fontsize',14,'horizontalalign','right');  
legend(p,'location','eastoutside','fontsize',14)
makepng(['SioAllMopProfiles.png'])

% vertical resolution for summed changes
Res=0.1;
% round survey x,y to desired resolution
zr=Res*round(z1da/Res); % round to Res meter spatial resolution

dz=z1da-z1dap;

% bin and average rounded survey data by placing in unique
%  x,y data array
[uz, ~, zidz] = unique(zr);
%array of counts of the number of points at each unique x/y combination
zcount = accumarray(zidz(:), 1); 
dz(isnan(dz))=0;
% array of total change in z that fall into each unique current z combination
dztot = accumarray(zidz(:), dz.');

figure('position',[ 108         215        1230         631]);
subplot(2,1,2)
bar(uz,dztot/5);set(gca,'xlim',[.5 3.5],'xtick',0.5:.1:3.5,'fontsize',14);grid on
idx=uz > .7 & dztot' < 0;
ev=sum(dztot(idx))/5;
idx=find(uz > .7 & dztot' > 0);
av=sum(dztot(idx))/5;
dv=sum(dztot(uz > .7)/5);

title(['Vol Changes Above MSL ( erosion | accretion | net ): ' ....
    num2str(ev,'%4.1f') ' | +' num2str(av,'%4.1f') ' | ' ...
    num2str(dv,'%4.1f') 'm^{3}/m'],'fontsize',16);  
xlabel('Current Profile Elevation (m, Navd88)');ylabel('Volume Change (m^{3}/m-shoreline)');
yl=get(gca,'ylim');set(gca,'ylim',yl);
hold on;
%set(gca,'xlim',[-20 90]);xl=get(gca,'xlim');
plot([2.119 2.119],yl,'k--');text(2.26,.95*yl(2),' HAT','fontsize',14,'horizontalalign','right');
plot([1.566 1.566],yl,'k--');text(1.7,.95*yl(2),' MHHW','fontsize',14,'horizontalalign','right');
plot([1.344 1.344],yl,'k--');text(1.44,.95*yl(2),' MHW','fontsize',14,'horizontalalign','right');
plot([.774 .774],yl,'k--');text(.9,.95*yl(2),' MSL','fontsize',14,'horizontalalign','right');

subplot(2,1,1)
bar(uz,zcount/5);set(gca,'xlim',[.5 3.5],'xtick',0.5:.1:3.5,'fontsize',14);grid on
xlabel('Current Profile Elevation (m, Navd88)');ylabel('Profile Area (m^{2}/m-shoreline)');
yl=get(gca,'ylim');
hold on;
%set(gca,'xlim',[-20 90]);xl=get(gca,'xlim');
plot([2.119 2.119],yl,'k--');text(2.26,.95*yl(2),' HAT','fontsize',14,'horizontalalign','right');
plot([1.566 1.566],yl,'k--');text(1.7,.95*yl(2),' MHHW','fontsize',14,'horizontalalign','right');
plot([1.344 1.344],yl,'k--');text(1.44,.95*yl(2),' MHW','fontsize',14,'horizontalalign','right');
plot([.774 .774],yl,'k--');text(.9,.95*yl(2),' MSL','fontsize',14,'horizontalalign','right');


