clear all
%close all
%addpath /volumes/group/MOPS
addpath /Users/William/Desktop/MOPS
for MopNumber=676%722%638:667%550:550;%582:582
%nr=10; % number of most recent profiles to plot
nr=16; % number of most recent profiles to plot
nr=10;
nr=50;
nr=8;

figure('position',[124         205        1140         480]);
 %----- add cobble sightings
 
matfile=['/Users/William/Desktop/MOPS/M' num2str(MopNumber,'%5.5i') 'SA.mat'];
load(matfile,'SA');
 %  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% divide mop area into 20 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,20,[-100 0]);

idx=find(vertcat(SA.Class) > 1);
Xutm=vertcat(SA.X);Yutm=vertcat(SA.Y);
Xutm=Xutm(idx);Yutm=Yutm(idx);
Z=vertcat(SA.Z);Z=Z(idx);
dt=[];
for n=1:size(SA,2)
    dt=[dt' SA(n).Datenum*ones(size(SA(n).Z))']';
end
dt=dt(idx);

[dp,NearIdx]=...
    pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);

[row,col] = ind2sub(size(xst),NearIdx);

hold on;pc=plot(x1d(col),Z,'m.','DisplayName','Past ATV Cobble Sightings');

%--------------

load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');

figure('position',[ 124         205        1140         480]);
nn=0;
% for n=1:1 
%     m=find(strcmp({SM.Source},'USACE') & [SM.Datenum] == datenum(2014,9,4));
%     %m=size(SM,2)-nr+n;
% %     if n==1;SM(m).Datenum=datenum(2021,10,11);end
% %     if n==2;SM(m).Datenum=datenum(2021,10,12);end
% %     if n==3;SM(m).Datenum=datenum(2021,10,13);end
%     if ~isnan(min(SM(m).Z1Dtransect))
%     nn=nn+1;
%     z=SM(m).Z1Dtransect;%z(z <-1.0)=NaN;
%     xMSL=intersections([SM(m).X1D(1) SM(m).X1D(end)],[0.774 0.774],SM(m).X1D,SM(m).Z1Dtransect);
%     if isempty(xMSL);xMSL=999.9;end
%     p(nn)=plot(SM(m).X1D,z,'-','linewidth',2,'DisplayName',...
%         [datestr(SM(m).Datenum,'mm/dd/yy') ' xMSL=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
%     end
% end
dt=datetime([SM.Datenum],'convertfrom','Datenum');
mm=find( ( strcmp({SM.Source},'Trk') | strcmp({SM.Source},'AtvMR') ) & month(dt) == 9 );
col=jet(nr);
for n=1:nr    
    m=mm(numel(mm)-nr+n);
%     if n==1;SM(m).Datenum=datenum(2021,10,11);end
%     if n==2;SM(m).Datenum=datenum(2021,10,12);end
%     if n==3;SM(m).Datenum=datenum(2021,10,13);end
    if ~isnan(min(SM(m).Z1Dtransect))
    nn=nn+1;
    z=SM(m).Z1Dtransect;%z(z <-1.0)=NaN;
    if n == nr
        z=z+0.08; %AtvMR correction
    end
    %xMSL=intersections([SM(m).X1D(1) SM(m).X1D(end)],[0.774 0.774],SM(m).X1D,SM(m).Z1Dtransect);
    xMSL=intersections([SM(m).X1D(1) SM(m).X1D(end)],[1.344 1.344],SM(m).X1D,SM(m).Z1Dtransect);
    if isempty(xMSL);xMSL=999.9;end
    p(nn)=plot(SM(m).X1D,z,'-','color',col(nn,:),'linewidth',2,'DisplayName',...
        [datestr(SM(m).Datenum,'mm/dd/yy') ' xMHW=' num2str(xMSL(end),'%4.1f') 'm ' SM(m).Source]);hold on;
    end
end

% best historical fit index
% bfidx=279;
% nn=nn+1;
% z=SM(bfidx).Z1Dmean;z(z < -0.7)=NaN;
% p(nn)=plot(SM(bfidx).X1D,z,'k:','linewidth',2,'DisplayName',...
%         [datestr(SM(bfidx).Datenum,'mm/dd/yy') ' Best Past Fit z=0.5-1.0m']);
    
xl=get(gca,'xlim');
set(gca,'xlim',[-30 xl(2)]);
%yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);

plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);
%ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
set(gca,'xdir','reverse','fontsize',14);grid on;
legend([p],'location','eastoutside','numcolumns',1);
title(['MOP ' num2str(MopNumber) ' Transect Profiles']);
xlabel('Cross-shore Distance (m)');
ylabel('Elevation (m, NAVD88)');

if nr > 3
makepng(['MOP' num2str(MopNumber) 'profiles.png'])
else
  makepng(['MOP' num2str(MopNumber) 'SeptProfiles.png'])  
end
end  