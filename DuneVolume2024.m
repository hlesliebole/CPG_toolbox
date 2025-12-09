clearvars
% estimate the volume of the cardiff dune
DefineMopPath
% combine Mop SG files for 668-679
Mop1=670;Mop2=678;
CG=SGcombineMops(Mop1,Mop2);
SurvNum=find([CG.Datenum] == datenum(2015,10,6));
%SurvNum=find([CG.Datenum] == datenum(2023,1,9));
SurvNum=find([CG.Datenum] == datenum(2019,9,12));
[X,Y,Z]=SG2grid(CG,SurvNum);
Z(Z > 4.7)=NaN; % clip back of beach at approx base of road riprap
Z(Z < 3.5)=NaN;
%figure;surf(X,Y,Z);shading flat;BeachBarColorbar;
SurvNum=find([CG.Datenum] == datenum(2022,1,19));
%SurvNum=find([CG.Datenum] == datenum(2019,9,12));
SurvNum=find([CG.Datenum] == datenum(2023,1,9));
SurvNum=find([CG.Datenum] == datenum(2024,1,9));
SurvNum=SurvNum(2);

[X2,Y2,Z2]=SG2grid(CG,SurvNum);
%figure;surf(X2,Y2,Z2);shading flat;BeachBarColorbar;
% step thru northing grid lines of truck survey and find
%  last dune point.  If above 4.7m NAVD88, interpolate
%  down to 4.7m in 3 1m grid point step to complete back of
%  dune based on design shape.
for n=1:size(Z2,1)
    ixmax=find(~isnan(Z2(n,:)), 1, 'last' );
    if ~isempty(ixmax)
        if Z2(n,ixmax) > 4.7
            z4=Z2(n,ixmax)-[1 2 3 4]*(Z2(n,ixmax)-4.7)/4;
            Z2(n,ixmax+1:ixmax+4)=z4;
        end
    end
end
% find any missing grid points in the baseline grid where
%  the recent truck survey has an elevation over 4.7m. Set
%  the elevation to 4.7m in the baseline grid.
imiss=find( Z2(:) > 4.7 & isnan(Z(:)) );Z(imiss)=4.7;

Zd=Z2-Z;Zd(Z < 3.2)=NaN;

figure;surf(Zd);shading flat;colormap(jet);view(2);
%Zd(Zd < 0)=NaN;

DuneVol=round(sum(Zd(:),'omitnan'));
DuneVolperShorelength=round(sum(Zd(:),'omitnan'))/((Mop2-Mop1+1)*100);
DuneVolperCardiffShorelength=round(sum(Zd(:),'omitnan'))/((682-666+1)*100);

idx=find(~isnan(Z(:)));
[y,x]=utm2deg([min(X(idx)) max(X(idx))],[min(Y(idx))-400 max(Y(idx))+400],...
    repmat('11 S',[2 1]));
figure('position',[ 215          93        1060         700]);
ax1=axes;set(ax1,'xlim',x,'ylim',y);hold on;
% load cddem.mat
% imAlpha=ones(size(cddem));imAlpha(isnan(cddem))=0;
plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)
hold on;
[ylat,xlon]=utm2deg(X(:),Y(:),repmat('11 S',[numel(X(:)) 1]));
Xl=X;Yl=Y;Xl(:)=xlon;Yl(:)=ylat;
surf(Xl,Yl,Zd);shading flat;colormap(jet);view(2);set(gca,'clim',[0 3])
view(0,90);set(gca,'dataaspectratio',[1.0000    2.8255   50.0000]);
cb=colorbar;cb.Label.String='Dune Elevation above Oct 2015 Baseline Beach (m)';
cb.FontSize=14;
set(gca,'fontsize',12);
% zoom in on dune
set(gca,'xlim',[-117.2821 -117.2767]);
set(gca,'ylim',[33.0004   33.0157]);


dim = [.25 .05 .3 .3];
str = [{'                 \bf ERODED DUNE VOLUME \rm'}, {[num2str(DuneVol) ' m^{3}']},{},...
    {[num2str(DuneVolperShorelength,'%4.1f') ' m^{3}/m-Dune Shorelength']}];%,...
    %{[num2str(DuneVolperCardiffShorelength,'%4.1f') ' m^{3}/m-Cardiff SB Shorelength']}];
if exist('a1','var');delete(a1);end
a1=annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'backgroundcolor','w','fontsize',14);
% title('Cardiff SB Living Shoreline Dune : AS BUILT (SIO 1/19/2022 Truck LiDAR Survey)',...
%     'fontsize',16)
% print(gcf,'-dpng','-r100','-loose','CardiffDuneAsBuilt.png')

title('Cardiff SB Living Shoreline Dune : Post-Storm (SIO 1/9/2024 ATV LiDAR Survey)',...
    'fontsize',16)
print(gcf,'-dpng','-r100','-loose','CardiffDune9Jan2024.png')


% [y,x]=utm2deg([min(X(idx)) max(X(idx))],[min(Y(idx)) max(Y(idx))],...
%     repmat('11 S',[2 1]));
% figure('position',[501 503 713/3 855]);
% ax1=axes;set(ax1,'xlim',[x],'ylim',[y]);hold on;
% % load cddem.mat
% % imAlpha=ones(size(cddem));imAlpha(isnan(cddem))=0;
% plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)
% hold on;
% surf(Xl,Yl,Zd);shading flat;colormap(jet);view(2);set(gca,'clim',[0 3])

% CF profiles = mops 673 677

% SurvNum=find(strcmp({CG.Source},'USACE') & year(datetime([CG.Datenum],'convertfrom','datenum')) == 2014);
% [X,Y,ZS]=SG2grid(CG,SurvNum);
% figure;surf(X,Y,ZS);shading flat;BeachBarColorbar;

% % relevant surveys
% SurveyDatenum3=datenum(2018,8,9);SurveySource3='Trk';% pre-dune but some beach nourishing?
% SurveyDatenum2=datenum(2016,3,3);SurveySource2='Gps';% post el nino 
% %SurveyDatenum3=datenum(2016,3,22);SurveySource3='KMair';% post el nino 
% SurveyDatenum3=datenum(1998,4,8);SurveySource3='UTAir';% post el nino
% SurveyDatenum1=datenum(2015,10,6);SurveySource1='KMair';% Melville airborne lidar prior to el nino
% %SurveyDatenum2=datenum(2016,9,29);SurveySource2='Gps';% post-el nino fall jumbo
% SurveyDatenum4=datenum(2019,9,12);SurveySource4='Trk';% post-dune 
% SurveyDatenum5=datenum(2010,9,23);SurveySource5='CCC';% post- el nino
% %SurveyDatenum5=datenum(2021,10,11);SurveySource5='Trk';% last fall 
% SurveyDatenum6=datenum(2022,1,19);SurveySource6='Trk';% most recent
% 
% % plot profiles for those 4 dates at 676
% 
% figure('position',[440   304   591   493]);
% MopNumber=677;
% matfile=['M' num2str(MopNumber,'%5.5i') 'SM.mat' ];
% fprintf('Loading %s\n',matfile);
% load(matfile,'SM');
% 
% % plot road
% SurvNum=find([SM.Datenum] == SurveyDatenum5 & strcmp({SM.Source},SurveySource5));
% n=find(SM(SurvNum).X1D < 0);
% p5=plot(SM(SurvNum).X1D(n),SM(SurvNum).Z1Dtransect(n),'-','linewidth',2,'color',[.8 .8 .8],...
%     'DisplayName',[datestr(SurveyDatenum5) ' ' SurveySource5]);
% hold on;
% text(-17,5.2,'ROAD','color',[.4 .4 .4],'fontsize',12,'fontweight','bold')
% 
% SurvNum=find([SM.Datenum] == SurveyDatenum1 & strcmp({SM.Source},SurveySource1));
% p1=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',3,...
%     'DisplayName',[datestr(SurveyDatenum1) ' ' SurveySource1]);
% hold on;
% SurvNum=find([SM.Datenum] == SurveyDatenum2 & strcmp({SM.Source},SurveySource2));
% p2=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
%     'DisplayName',[datestr(SurveyDatenum2) ' ' SurveySource2]);
% % 
% SurvNum=find([SM.Datenum] == SurveyDatenum3 & strcmp({SM.Source},SurveySource3));
% p3=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
%     'DisplayName',[datestr(SurveyDatenum3) ' ' SurveySource3]);
% % 
% % 
% SurvNum=find([SM.Datenum] == SurveyDatenum4 & strcmp({SM.Source},SurveySource4));
% p4=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
%     'DisplayName',[datestr(SurveyDatenum4) ' ' SurveySource4]);
% % profile offset for location of dune design revetment and toe.
% % revetment estimted to be from -1.5 to +6m from the peak elev location 
% %  in the post construction truck survey. cobble toe +10.5m peak
% [zmax,imax]=max(SM(SurvNum).Z1Dtransect);x0=SM(SurvNum).X1D(imax);
% 
% % 
% % 
% SurvNum=find([SM.Datenum] == SurveyDatenum6 & strcmp({SM.Source},SurveySource6));
% p6=plot(SM(SurvNum).X1D,SM(SurvNum).Z1Dtransect,'-','linewidth',2,...
%     'DisplayName',[datestr(SurveyDatenum6) ' ' SurveySource6]);
% 
% plot([-60 80],[.774 .774],'k--');text(-50,.65,'MSL','fontweight','bold');
% 
% p7=plot([x0+6 x0-1.5],[1.22 3.96],'k:','linewidth',5,'DisplayName','Buried Revetment 2T Stone');
% p8=plot(x0+10.5,4,'m.','markersize',20,'DisplayName','Dune Design Toe');
% 
% set(gca,'xdir','reverse');grid on;
% legend([p3 p1 p2 p4 p6 p7 p8],' 8 Apr 1998 Airborne LiDAR; post-El Nino eroded beach',...
%     ' 6 Oct 2015  Airborne LiDAR; pre-El Nino BASELINE RECOVERED BEACH',...
%     ' 3 Mar 2016  ATV-Dolly Survey; post-El Nino eroded beach',...
%     '12 Sep 2019  Truck LiDAR; post-Dune Completion',...
%     '19 Jan 2022  Truck LiDAR; most recent survey',...
%     'Dune Design Buried Revetment 2T Stone','Dune Design Toe Location',...
%     'location','northwest','fontsize',12);
% set(gca,'xlim',[-60 80],'ylim',[0 9],'fontsize',12,'fontweight','demi')
% xlabel('Mop Cross-Shore Profile Distance (m)')
% ylabel('Elevation (m, NAVD88)')
% title(['Mop ' num2str(MopNumber) ' Profiles'],'fontsize',16,'fontweight','bold')
% 
% print(gcf,'-dpng','-r100','-loose','DuneProfiles.png');
% % 
% % 
% % 
% % 
% % 
% % [X,Y,Z]=CombineSGgrids(669,679,SG(176).Datenum,SG(176).Source);
% % figure;
% % surf(X,Y,Z);shading flat;BeachColorbar;set(gca,'dataaspectratio',[1 1 .1]);view(2);
% % for MopNumber=MopStart:MopEnd
% % PlotLabelMopTransectUTM(MopNumber,'3d','k','ShadeOff');