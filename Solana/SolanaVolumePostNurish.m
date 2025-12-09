%  Builds grid from combined truck lidar and gps data
%  
% 1. the combined 1m average surveys are first reduced to nearest point
%    profiles.
%
% 2. the profile elevations from -12 to MSL are converted to xg,yg space and gridded with a 250m
%    tolerance for added boundary NaNs
%
% 3. the ATVmr in xg,yg space is gridded separately with a 5m NaN boundary tolerance
%     and overlaid on top of the gridded profiles


%% 1. & 2.

load SolanaShoreboxMap.mat

%% Build a shorebox baseline elev grid 
% make 2d X,Y grid arrays
[X,Y]=meshgrid(xgimg,ygimg);
Z=nan(size(X));% initialize composite elev grid as no data NaNs 

%%  Combine jumbo and atv mr SA data files
figure;hold on;
sdate(1)=datenum(2024,4,10);
sdate(2)=datenum(2024,4,18);
dx=[];
dy=[];
% loop through even numbered mops that have dolly-jetski surveys  
for m=str2num(MopSB.Name{1}(2:5)):str2num(MopSB.Name{end}(2:5))
    matfile=['M00' num2str(m,'%3.3i') 'SA.mat'];
    load(matfile,'SA')
    % identify Gps (dolly-jetski) 
    ndx=find(strcmp({SA.Source},'Gps') & [SA.Datenum] == sdate(1));
    adx=find(strcmp({SA.Source},'Gps') & [SA.Datenum] == sdate(2));
    %adx=find(strcmp({SA.Source},'AtvMR') & [SA.Datenum] == sdate(1));
    %adx=[];
    % only use ATV mr data above 0.72 until processed level 2 file is cleaned up
    %idx=find([SA(adx).Z] > 0.72);
    % SA(adx).X=SA(adx).X(idx);
    % SA(adx).Y=SA(adx).Y(idx);
    % SA(adx).Z=SA(adx).Z(idx);
    % combine surveys under SA(1)
    SA=SA([ndx adx]);
    SA(1).X=vertcat(SA.X);
    SA(1).Y=vertcat(SA.Y);
    SA(1).Z=vertcat(SA.Z);
    % just keep SA(1)
    SA=SA(1);
    % reduce to a nearest point profiles, interpolate up to 100m in xshore
    %   to close gaps between dolly and jetski
    if numel([SA(1).X]) > 0
    %[X1D,Z1D]=GetAllNearestPointsProfiles(SA,25,100);
    plot(SA(1).X,SA(1).Y,'.')
    % convert profiles to utm coords
    %[Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(m,X1D);
    %  convert utm corrds to shorebox coords
    [Xsb,Ysb]=utm2ShoreBox(SA(1).X,SA(1).Y,str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));
    dx=[dx Xsb'];
    dy=[dy Ysb'];
    % add points to shorebox grid indices
    gdx=sub2ind(size(X),round(Ysb)-ygimg(1)+1,round(Xsb)-xgimg(1)+1);% 1d grid indices with data 
    Z(gdx)=SA(1).Z; % overlay on grid
    end
end

%%

% Grid the combined elevation surveys using Delaunay tesselation 
maxdist=350; %set distance threshold of Nan boundary points for 200m spaced profiles 
idx=find(~isnan(Z(:)));
[xadd,yadd,zadd]=addNoDataAreaPoints(X(idx),Y(idx),Z(idx),maxdist);
    Z=griddata(double(xadd),double(yadd),double(zadd),...
        double(xgimg),double(ygimg));

figure
%imagesc(xgimg(1050:3000),ygimg,Z(:,1050:3000));
imagesc(xgimg,ygimg,Z);
colorbar
set(gca,'clim',[-12 6],'ydir','normal')
colormap(jet);

%%  3.  truck in xg,yg space is gridded separately with a 5m NaN boundary tolerance
%     and overlaid on top of the gridded profiles


%% load combined SA grids for the reach
CS=SAcombineMops(str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));

sdate(1)=datenum(2024,4,16);

Z2=Z*nan;
for ns=1:1   
ndx=find(strcmp({CS.Source},'Trk')  & ...
             [CS.Datenum] == sdate(ns));
for idx=ndx   
 zdx=find(CS(idx).Z > -2);
 if numel(zdx) > 0
  [Xsb,Ysb]=utm2ShoreBox(CS(idx).X(zdx),CS(idx).Y(zdx),str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));

  gdx=sub2ind(size(X),round(Ysb)-ygimg(1)+1,round(Xsb)-xgimg(1)+1);% 1d grid indices with data 
    Z2(gdx)=CS(idx).Z(zdx); % overlay on grid
 end
end

end

% ndx=find(strcmp({CS.Source},'Multibeam'));
% for idx=ndx   
%  zdx=find(CS(idx).Z > -15);
%  if numel(zdx) > 0
%   [Xsb,Ysb]=utm2ShoreBox(CS(idx).X(zdx),CS(idx).Y(zdx),str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));
% 
%   gdx=sub2ind(size(X),round(Ysb)-ygimg(1)+1,round(Xsb)-xgimg(1)+1);% 1d grid indices with data 
%     Z2(gdx)=CS(idx).Z(zdx); % overlay on grid
%  end
% end


% 
% 
% zmax=5;
% 
% for i=1:3577    
%     j=min([find(~isnan(Z2(:,i)))]);
%     if isempty(j)
%         jmin=1;
%     else
%         jmin=j-1;
%     end
%     Z2(jmin,i)=10000;
% end

% Grid the combined elevation surveys using Delaunay tesselation 
%   tight boundary constraint for lidar data
maxdist=5;
idx=find(~isnan(Z2(:)));
[xadd,yadd,zadd]=addNoDataAreaPoints(X(idx),Y(idx),Z2(idx),maxdist);
    Z2=griddata(double(xadd),double(yadd),double(zadd),...
        double(xgimg),double(ygimg));
%Z2(Z2 > 5)=NaN;
%Z2(Z2 < -1)=NaN;

%%
figure
%imagesc(xgimg(1050:3000),ygimg,Z(:,1050:3000));
imagesc(xgimg,ygimg,Z2);
colorbar
set(gca,'clim',[-2 4],'ydir','normal')
colormap(jet);
%%
Z1=Z;
Z0=Z;%Z0(Z0 > 1.5)=NaN;
idx=find(~isnan(Z2(:)));
Z0(idx)=Z2(idx);
figure
%imagesc(xgimg(1050:3000),ygimg,Z(:,1050:3000));
imagesc(xgimg,ygimg,Z0);
colorbar
set(gca,'clim',[-12 5],'ydir','normal')
colormap(jet);
hold on;
plot(dx,dy,'k.')
%%
Z=Z0;
save SolanaPostnourishmentGrid.mat X Y Z

%%
% % second layer is multibeam
% 
% 
% % third layer is truck or with z >= 0.5m navd88
% 
% 
% CS=SGcombineMops(str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));
% ndx=find(strcmp({CS.Source},'Multibeam'));
% ndx=find(strcmp({CS.Source},'USACE'));
% % load tmpSG.mat
% % CS=SG;
% %%
% ndx=find(strcmp({CS.Source},'Gps'));
% ndx=ndx(end-3):ndx(end);
% for idx=ndx 
%     if(strcmp(CS(idx).Source,'Gps'))
%         zdx=find(CS(idx).Z < 0.5);
%     else
%         zdx=find(CS(idx).Z < 10);
%     end
%     zoff=0;if(CS(idx).Datenum == 739312);zoff=0.7;end
% [Xsb,Ysb]=utm2ShoreBox(CS(idx).X(zdx),CS(idx).Y(zdx),str2num(MopSB.Name{1}(2:5)),str2num(MopSB.Name{end}(2:5)));
% [ScatterPlot,ColorBarPlot]=ColorScatter2dTopoBathy(Xsb,Ysb,CS(idx).Z(zdx)+zoff);
% %pl=plot(Xsb,Ysb,'.','markersize',1);
% end

%set(gca,'xlim',[-200 600],'ylim',[100 2100])
% 
% iii=0;
% for nn=1:3
%     for mm=1:6
%         iii=iii+1;
%         axs(iii)=axes('position',[.06+(nn-1)*.23 .055+(mm-1)*.152 .21 .145],'box','on');
%         set(axs(iii),'xlim',[-20 90],'xdir','reverse','ylim',[-.5 3.5],'fontsize',12);
%         grid on;
%         hold on;
%         MopNum=515-(nn-1)*6-(7-mm);%499+iii;
%         %PlotCanyonProfile;   
%         if iii == 3% || iii == 9
% 
%             ylabel('Elev (m, NAVD88)');
% 
% 
%         end
%     if iii == 7
%         xlabel('Cross-Shore Distance (m)')
%         else
%             set(axs(iii),'xticklabels',' ')
%     end
%         % if iii == 8;title('SOUTH CANYON HEAD MOPs');end
%         % if iii == 16;title('NORTH CANYON HEAD MOPs');end
%         % if MopNum < 518
%         % legend([p1 p2],'2003-2007 only OBS','8 Feb 2023 JETSKI',...
%         %         'location','northwest','fontweight','bold')
%         % %text(550,-8,['MOP ' num2str(MopNum)],'fontsize',14,'fontweight','bold')
%         % elseif MopNum > 517 && MopNum < 520
%         %     legend([p1 p2],'2003-2016 only OBS','8 Feb 2023 JETSKI',...
%         %         'location','northwest','fontweight','bold')
%         % else
%         %     legend([p1 p2],'2003-2022 OBS','8 Feb 2023 JETSKI',...
%         %         'location','northwest','fontweight','bold')
%         % end
%         if iii > 20
%             text(550,-14,['MOP ' num2str(MopNum)],'fontsize',14,'fontweight','bold')
%         else
%             text(75,2.6,['MOP ' num2str(MopNum)],'fontsize',16,'fontweight','bold')
%         end
% 
%         if iii == 3
% xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
% %plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');text(xl(2),1.8,' MHHW','fontsize',14);
% %plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
% plot(xl,[.774 .774],'k-');text(xl(2),1.,' MSL','fontsize',14);
% plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.15,' MLLW','fontsize',14); 
%         else
%             xl=get(gca,'xlim');%xl(2)=90;set(gca,'xlim',xl);
% %plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
% plot(xl,[1.566 1.566],'k--');%text(xl(2),1.8,' MHHW','fontsize',14);
% %plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
% plot(xl,[.774 .774],'k-');%text(xl(2),1.,' MSL','fontsize',14);
% plot(xl,[-0.058 -0.058],'k--');%text(xl(2),0.15,' MLLW','fontsize',14);
%         end
% 
% load(['M00' num2str(MopNum,'%3.3i') 'SA.mat'])
% [x1d,z1di]=GetNonGriddedProfile(MopNum,2);
% eln=plot(x1d,z1di,'k-','linewidth',4,'displayname',[datestr(SA(2).Datenum) ' NASA ATM']);
% hold on;
% idx=find(strcmp({SA.Source},'Trk'));
% n=idx(end-1);
% [x1d,z1di]=GetNonGriddedProfile(MopNum,n);
% pl1=plot(x1d,z1di,'c-',...
%         'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
%         'linewidth',4);
% 
% 
% n=size(SA,2);
% [x1d,z1di]=GetNonGriddedProfile(MopNum,n);
% pl2=plot(x1d,z1di,'r-',...
%         'DisplayName',[datestr(SA(n).Datenum) ' ' SA(n).Source],...
%         'linewidth',4);
% if MopNum == 509
%     lg=legend([eln pl1 pl2],'location','northoutside','numcolumns',3,'fontsize',16);
%     lg.Position=[0.2697 0.9692 0.2469 0.0277];
% end
% 
%     end
% end
% 
% makepng('LaJollaShores11Jan24.png')


% % sl5=text(-170,-435,{'Torrey Pines',' State Beach'},'color','k','fontsize',14,'fontweight','bold','rotation',90);
% % sl5.Extent;x1=sl5.Extent(1);x2=x1+sl5.Extent(3);y1=sl5.Extent(2);y2=y1+sl5.Extent(4);
% % fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[.8 .8 .8],'EdgeColor',...
% %     [.8 .8 .8],'FaceColor',[.8 .8 .8],'FaceAlpha','.7');
% % sl5=text(-170,-435,{'Torrey Pines',' State Beach'},'color','k','fontsize',14,'fontweight','bold','rotation',90);
% % 
% % sl5=text(4200,-435,{'    San','Dieguito','   River'},'color','k','fontsize',14,'fontweight','bold','rotation',90);
% % sl5.Extent;x1=sl5.Extent(1);x2=x1+sl5.Extent(3);y1=sl5.Extent(2);y2=y1+sl5.Extent(4);
% % fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[.8 .8 .8],'EdgeColor',...
% %     [.8 .8 .8],'FaceColor',[.8 .8 .8],'FaceAlpha','.6');
% % sl5=text(4200,-435,{'    San','Dieguito','   River'},'color','k','fontsize',14,'fontweight','bold','rotation',90);
% % sl5=text(2600,-435,{'Powerhouse','     Park'},'color','k','fontsize',14,...
% %     'fontweight','bold','rotation',90);
% % sl5.Extent;x1=sl5.Extent(1);x2=x1+sl5.Extent(3);y1=sl5.Extent(2);y2=y1+sl5.Extent(4);
% % fill([x1 x2 x2 x1 x1],[y1 y1 y2 y2 y1],[.8 .8 .8],'EdgeColor',...
% %     [.8 .8 .8],'FaceColor',[.8 .8 .8],'FaceAlpha','.6');
% % sl5=text(2600,-435,{'Powerhouse','     Park'},'color','k','fontsize',14,...
% %     'fontweight','bold','rotation',90);
% % 
% % %a1=text(2300,400,'\leftarrow','color','red','fontsize',40,'fontweight','bold');
% % a1=text(1900,400,'\leftarrow','color','red','fontsize',25,'fontweight','bold');
% % a1=text(1500,400,'\leftarrow','color','red','fontsize',30,'fontweight','bold');
% % a1=text(1100,400,'\leftarrow','color','red','fontsize',40,'fontweight','bold');
% % a1=text(700,400,'\leftarrow','color','red','fontsize',50,'fontweight','bold');
% % a1=text(400,400,'\leftarrow','color','red','fontsize',40,'fontweight','bold');
% % a1=text(100,400,'\leftarrow','color','red','fontsize',25,'fontweight','bold');
% % 
% % % a1=text(1900,400,'\leftarrow','color','red','fontsize',20,'fontweight','bold');
% % % a1=text(1500,400,'\leftarrow','color','red','fontsize',30,'fontweight','bold');
% % a1=text(4700,400,'\leftarrow','color','red','fontsize',40,'fontweight','bold');
% % a1=text(4300,400,'\leftarrow','color','red','fontsize',50,'fontweight','bold');
% % a1=text(3950,400,'\leftarrow','color','red','fontsize',35,'fontweight','bold');
% % a1=text(3700,400,'\leftarrow','color','red','fontsize',25,'fontweight','bold');
% % a1=text(3450,400,'\leftarrow','color','red','fontsize',20,'fontweight','bold');
% % 
% % a1=text(3100,400,'\leftarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(3100,400,'\rightarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(2700,400,'\leftarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(2700,400,'\rightarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(2200,400,'\leftarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(2200,400,'\rightarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(5100,400,'\leftarrow','color','green','fontsize',30,'fontweight','bold');
% % a1=text(5100,400,'\rightarrow','color','green','fontsize',30,'fontweight','bold');
% % %a1=text(4000,400,'\leftarrow','color','red','fontsize',25,'fontweight','bold');
% % 
% 
% % ar = annotation("arrow",'x',[.2 .23],'y',[.2 .23],'linewidth',20);
% % c = ar.Color;
% % ar.Color = "red";
% % 
% % % Sxy axes
% % ax2=axes('position',[pos(1) .55 pos(3) .4],'fontsize',14,'linewidth',2);
% % set(ax2,'xlim',xl);hold on;box on
% % 
% % load DelMar594to646netTransport.mat
% % net(net(:)==0)=NaN;
% % % col=jet(size(net,1));
% % % figure;for n=1:size(net,1);p(n)=plot(MopYears,net(n,:),'.-','color',col(n,:),'DisplayName',num2str(MopNumbers(n)));hold on;end
% % % legend(p(1:2:end),'location','eastoutside')
% % 
% % col=jet(size(net,2));
% % for n=1:size(net,2)
% %     if MopYears(n) == 2016
% %             pl(n)=plot(MopSB.BackLon(:),net(:,n),'.-','color','m','linewidth',3,'markersize',15,'DisplayName',[num2str(MopYears(n))]);hold on;
% %     elseif MopYears(n) == 2010
% %         pl(n)=plot(MopSB.BackLon(:),net(:,n),'.:','color','m','linewidth',3,'markersize',15,'DisplayName',[num2str(MopYears(n))]);hold on;
% %     elseif MopYears(n) == 2023
% %         pl(n)=plot(MopSB.BackLon(:),net(:,n),'.-','color','k','linewidth',4,'markersize',15,'DisplayName',[num2str(MopYears(n))]);hold on;
% %     else
% %     pl(n)=plot(MopSB.BackLon(:),net(:,n),'.-','color',col(n,:),'linewidth',1,'markersize',15,'DisplayName',num2str(MopYears(n)));hold on;
% %    end
% % end
% % plot(xl,[0 0],'k:','linewidth',2)
% % set(gca,'fontsize',14,'linewidth',2,'ylim',[-175 125]);grid on;
% % lg=legend(pl,'location','east');
% % pos2=get(lg,'position');set(lg,'position',[.92 0.48 0.0526 0.4993])
% % ylabel({'Net Annual Relative Longshore Transport','K\cdot\Sigma{H_s}^{1/2} S_{xy} (m^{3}/yr)' });
% % sl=text(-250,-120,{'SOUTH   <<---','TRANSPORT'},'fontsize',16,'fontweight','bold','rotation',90);
% % sl2=text(-250,5,{'--->  NORTH',' TRANSPORT'},'fontsize',16,'fontweight','bold','rotation',90);
% % sl3=text(3300,-130,'El Nino Winters -->','color','m','fontsize',16,'fontweight','bold');
% % sl4=text(2500,80,{'\color[rgb]{0, .5, 0}LOW TRANSPORT','\color{black}NORTH DEL MAR'},'color','k','fontsize',16,'fontweight','bold');
% % sl5=text(600,80,{'\color{red}HIGH SOUTH TRANSPORT','      \color{black}SOUTH DEL MAR'},'color','k','fontsize',16,'fontweight','bold');
% % sl6=text(4000,80,{'\color{red}HIGH SOUTH TRANSPORT','      \color{black}SOUTH SOLANA BCH'},'color','k','fontsize',16,'fontweight','bold');
% % 
% addpath /Users/William/Desktop/SANDAG
% load Mops512to929GMP.mat
% % EGMP=GetEnsembleMopProfile(512,683,0.1);
% % save TP2CardiffEnsembleMopProfiles.mat EGMP
% %load TP2CardiffEnsembleMopProfiles.mat 
% load TP2OsideEnsembleMopProfiles.mat 
% 
% Nmops=size(EGMP,2); % number of mops in reach
% MopNumbers=[EGMP.MopBase]; % vector of mop numbers
% 
% % turn EGMP into a mop number vs xshore distance relative to the global
% % profiles' navd88 = 0 locations.
% 
% ZofXorigin=0.774; % global profile elevation to use as the xshore xorigin
% 
% % loop through mop global profile and adjust their xshore axes to the
% %  alongshore global contour elevation-based reference frame and only
% %  keep the xshore ramge with valid global profile data at each mop.
% %
% xmin=[];
% xmax=[];
% Vxs=0*MopNumbers';
% for n=1:Nmops
%     xdata=find(~isnan(EGMP(n).Zglobal)); % xshore points with data 
%     EGMP(n).X=EGMP(n).X(xdata); % reduce to locations with data
%     EGMP(n).Zglobal=EGMP(n).Zglobal(xdata); % reduce to elevations with data
%     %EGMP(n).AnnualAnomalies=EGMP(n).AnnualAnomalies(:,xdata); % reduce to anoms with data
%     EGMP(n).QuarterlyAnomalies=EGMP(n).QuarterlyAnomalies(:,xdata); % reduce to anoms with data
%     x0=find(EGMP(n).Zglobal > ZofXorigin , 1, 'last' ); % find xshore indice of the origin
%     if numel(x0) > 0
%     EGMP(n).X0=EGMP(n).X(x0);
%     EGMP(n).X=EGMP(n).X-EGMP(n).X(x0); % adjust xshore locations relative to the new origin   
% 
%     end
%     xmin=min([xmin min(EGMP(n).X)]);
%     xmax=max([xmax max(EGMP(n).X)]);
% end
% 
% Z2D=NaN(Nmops,numel(xmin:xmax));
% 
% %figure('position',[214 450 1045 350]);
% for n=1:Nmops
%     ix=find(ismember(xmin:xmax,EGMP(n).X));
%     Z2D(n,ix)=EGMP(n).Zglobal;
% end
% 
% QA2D=NaN(Nmops,numel(xmin:xmax));
% for n=1:Nmops
% ndx=find(month(datetime([EGMP(n).QuarterlyDatenums],'convertfrom','datenum')) == 8);
% ix=find(ismember(xmin:xmax,EGMP(n).X));
% QA2D(n,ix)=mean(EGMP(n).QuarterlyAnomalies(ndx,:),'omitnan');
% end
% 
% QA2DS=NaN(Nmops,numel(xmin:xmax));
% for n=1:Nmops
% ndx=find(month(datetime([EGMP(n).QuarterlyDatenums],'convertfrom','datenum')) == 5);
% ix=find(ismember(xmin:xmax,EGMP(n).X));
% QA2DS(n,ix)=mean(EGMP(n).QuarterlyAnomalies(ndx,:),'omitnan');
% end
% 
% % %contour(Z2D',[-8 -6 -4 -2 0 2],'linecolor','k');hold on;
% % imagesc(MopNumbers',xmin:xmax,Z2D','AlphaData',~isnan(Z2D'));
% % hold on;
% % [C,h]=contour(MopNumbers',xmin:xmax,Z2D',[-8 -6 -4 -2 0 2],'linecolor','w');
% % clabel(C,h,'color','w','labelspacing',1000);
% % shading flat;BeachColorbar
% % %set(gca,'ydir','normal','xtick',520:5:595);
% % set(gca,'ydir','normal','xtick',512:10:929);
% % 
% % xlabel('Mop Number');ylabel('Xshore Distance (m)')
% % title('Torrey Pines SB : Global Mean Jumbo Survey TopoBathy (m, NAVD88)')
% % %pause
% 
% %figure('position',[214 50 1045 850]);
% MinBeachYear=min([EGMP.AnnualBeachYears]);
% MaxBeachYear=max([EGMP.AnnualBeachYears]);
% 
% for by=2023%MinBeachYear:MaxBeachYear
%     fprintf('%6i: ',by)
% APdiff=NaN(Nmops,numel(xmin:xmax));    
%     %if by > MinBeachYear
%         %AP=A2D;%AP(AP < 0)=0;AP(isnan(AP))=0;    
%     %end
%     A2D=NaN(Nmops,numel(xmin:xmax));
% %by=2016;
% for n=1:Nmops
%     % 2d array xshore indices for this mop's profile
%     %ix=find(ismember(EGMP(n).X,xmin:xmax));
%     %iy=find(EGMP(n).AnnualBeachYears == by); 
%     iy=find(EGMP(n).QuarterlyDatenums == datenum(by,8,15));
%     if ~isempty(iy)
%     ix=find(ismember(xmin:xmax,EGMP(n).X));
%     %A2D(n,ix)=EGMP(n).AnnualAnomalies(iy,:);
%     A2D(n,ix)=EGMP(n).QuarterlyAnomalies(iy,:)-QA2D(n,ix);
%     AP(n,ix)=Z2D(n,ix);
% 
%     end
% end
% 
% mdx=find(MopNumbers > 593 & MopNumbers < 647);
% 
% addpath /Users/William/Desktop/MOPS/toolbox
% 
% hold on;
% for mm=mdx
%  sidx=find(strcmp({MopSB.Name{:}},['D0' num2str(MopNumbers(mm))]));
%  xs=MopSB.BackLon(sidx); % shorebox mop back beach x
%  ys=MopSB.BackLat(sidx); % shorebox mop back beach y
%  %plot(xs,ys,'r+')
%  iy=find(~isnan(A2D(mm,:)));
%  Y=xmin:xmax;Y=ys+Y(iy)+EGMP(mm).X0;% mop xshore location 
%  AX=A2D(mm,iy); % mop xshore anom
%  ColorScatterPolarmap(Y*0+xs,Y,AX)
% 
%  %[Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(mm,X);
%  % [sbx,sby]=LatLon2shorebox(Lat,Lon);
%  % polarscatterplot(sbx,sby);
% end
% 
% % 
% % 
% % if by > 2001 & by < 2024
% % % clf;
% % % subplot(2,1,1)
% % ax1=axes('position',[.121 .052 .492 .94]);
% % imagesc(xmin:xmax,MopNumbers',A2D,'AlphaData',~isnan(A2D));hold on;
% % [C,h]=contour(xmin:xmax,MopNumbers',Z2D-.774,[-10 -4 0],'linecolor','w','linestyle','-','linewidth',4);
% % [C,h]=contour(xmin:xmax,MopNumbers',Z2D-.774,[-10 -4 0],'linecolor','k','linestyle','-','linewidth',3);
% % % clabel(C,h,'color','w','labelspacing',1400,'fontsize',24,'fontweight','bold');
% % % clabel(C,h,'color','k','labelspacing',1400,'fontsize',22,'fontweight','bold');
% % text(600,700,'-10m','FontSize',24,'fontweight','bold','Color','k')
% % text(320,700,'-4m','FontSize',24,'fontweight','bold','Color','k')
% % text(30,700,'MSL','FontSize',24,'fontweight','bold','Color','k','rotation',90)
% % shading flat;
% colormap(flipud(polarmap));set(gca,'clim',[-1 1]);
% % set(gca,'ydir','normal','fontsize',12);
% % set(gca,'ydir','normal','xdir','reverse','ytick',[]);
% % set(gca,'xlim',[-100 1000])
% % grid on;
% % set(gca,'color','none')
% % %xlabel('Mop Number');
% % xlabel('Xshore Distance (m)');
% % %xlabel('                   South Torrey                           North Torrey       S Del Mar      N Del Mar          Solana Bch        Cardiff  ')
% %cb=colorbar('position',[.20 .15 .015 .65],'fontsize',12,'fontweight','bold');
% cb=colorbar('horizontal');
% set(cb,'position',[.33 .8 .4 .03]);
% cb.Label.String='ELEVATION ANOMALY (m)';
% cb.FontSize=16;
% % set(gca,'fontweight','bold','fontsize',14)
% % % title(sprintf('Summer %i Nearshore Elevation Anomalies',by),'fontsize',20);
% % % t1=text(575,600,'Long-Term Mean Depth Contours (m,NAVD88)','fontsize',12,'fontweight','bold')
% %  set(gcf,'inverthardcopy','off')
% end
% 
% set(gca,'xlim',[min(xgimg) max(xgimg)],'ylim',[min(ygimg) max(ygimg)])
% makepng('DelMar2023Anomalies.png')
% %makepng('FigureS2.png')
% 
% function [ScatterPlot]=ColorScatterPolarmap(x,y,d)
% 
% [ad,ix]=sort(abs(d));d=d(ix);x=x(ix);y=y(ix);
% 
% %zmin=quantile(z,.05);zmax=quantile(z,.95);
% zmin=-1;zmax=1;
% zrange=zmax-zmin; % set max slope for coloring
% zscaled = 1+64*(d-zmin)/zrange;
% zscaled(zscaled < 1)=1;zscaled(zscaled > 64)=64;
% idx=find(~isnan(zscaled)); % non NaN points 
% 
% load PolarMap.mat
% cm = flipud(PolarMap);
% 
% scp=scatter(x(idx), y(idx), 75, cm(ceil(zscaled(idx)),:), 'filled');
% scp.MarkerFaceAlpha = .9;
% scp.MarkerEdgeAlpha = .9;
% end