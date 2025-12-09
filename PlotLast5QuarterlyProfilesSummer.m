
% Makes quarterly prifle plots for an defined beach year

close all
clearvars

Byear=2024; % Survey beach year to plot quarterly means

% get MOP transect location info
load MOPTableUTM.mat

% open mop prolile definition file
fid=fopen('MopLatLonEutnNutmOrientation.dat','w');

% define 9 state beach names for files
bname{1}='DelMar';
% bname{2}='BlacksNorth';
% bname{3}='Cardiff';
% bname{2}='SilverStrand';
% bname{3}='TorreyPines';
% bname{4}='Cardiff';
% bname{5}='SanElijo';
% bname{6}='Moonlight';
% bname{7}='Leucadia';
% bname{8}='SouthCarlsbad';
% bname{9}='Carlsbad';

% define 9 state beach names for titles
btitle{1}='Del Mar';
% btitle{2}='Blacks North';
% btitle{3}='Cardiff SB';
% btitle{2}='Silver Strand SB';
% btitle{3}='Torrey Pines SB';
% btitle{4}='Cardiff SB';
% btitle{5}='San Elijo SB';
% btitle{6}='Moonlight SB';
% btitle{7}='Leucadia SB';
% btitle{8}='South Carlsbad SB';
% btitle{9}='Carlsbad SB';

% define 9 mop ranges for the state beaches

MopRng=[624 624];

% MopRng=[624 624
%         552 552
%         674 674];

% MopRng=[2 28 
%         85 157
%         536 607
%         664 683
%         684 706
%         717 725
%         740 757
%         762 819
%         825 854];

% legend names for 4 profiles
DispName(1,:)=[num2str(Byear-1) ' Summer Quarter (JAS)' ];
DispName(2,:)=[num2str(Byear-1) '   Fall Quarter (OND)' ];
DispName(3,:)=[num2str(Byear) ' Winter Quarter (JFM)' ];
DispName(4,:)=[num2str(Byear) ' Spring Quarter (AMJ)' ];
DispName(5,:)=[num2str(Byear) ' Summer Quarter (JAS)' ];


% line colors for 4 profiles
col(1,:)=[0.7 0 0];
col(2,:)=[0 0 0.7];
col(3,:)=[0.7 0.6 0];
col(4,:)=[0 0.7 0];
col(5,:)=[0.7 0 0];


%% loop thru state beaches

for sb=1:3%1:9

 % loop through mop range for this beach
 for MopNum=MopRng(sb,1):MopRng(sb,2)

    % load SA mat file
load(['M' num2str(MopNum,'%5.5i') 'SA.mat'],'SA');

%% get quarterly mean profiles

%  X0BeachOnly = Mop xshore location of truck back beach boundary
%  X1Dt(N) = N profile xshore x values (1m res) relative to Mop back beach point
%  QY(M) = M Quartery profile years
%  Q(M) = M Quartery profile quarters (1-4)
%  Zyq(M,N) = M QY(M)/Q(M) year/quarter mean profiles (Jan-Mar,Apr-Jun,Jul-Sep,Oct-Dec)
%  NS(M) = M QY(M)/Q(M) year/quarter number of surveys in quarter mean 

[X0BeachOnly,X1Dt,QY,Q,Zyq,NS]=GetMeanQuarterlyNearestProfiles(SA,Byear);

%  Get utm and latlon coords of xshore profile points
[Lat,Lon,Xutm,Yutm]=MopxshoreX2LatLonUTM(MopNum,X1Dt);

%figure;
for yq=1:numel(QY)
    %plot3(X1Dt,(QY(yq)+Q(yq)/4)*ones(size(X1Dt)),Zyq(yq,:),'-');hold on;
    ofile=[bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(QY(yq)) '_' num2str(Q(yq)) '.dat'];
    fprintf('%s\n',ofile)
    % find last valid data point in the profile
    ilast=find(~isnan(Zyq(yq,:)),1,'last');
    fid=fopen(ofile,'w');
    for m=1:ilast
        %fprintf(fid,'%i %8.3f\n',n-1,Zyq(yq,n));
        %for m=1:numel(X1Dt)
        fprintf(fid,'%i %10.1f %10.1f %13.7f %13.7f %8.2f\n',...
            X1Dt(m),Xutm(m),Yutm(m),Lon(m),Lat(m),Zyq(yq,m));
        %end
    end
    fclose(fid);
end

[X0BeachOnly,X1Dt,QY,Q,Zyq,NS]=GetMeanQuarterlyNearestProfiles(SA,Byear+1);
for yq=1%:numel(QY)
    %plot3(X1Dt,(QY(yq)+Q(yq)/4)*ones(size(X1Dt)),Zyq(yq,:),'-');hold on;
    ofile=[bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(QY(yq)) '_' num2str(Q(yq)) '.dat'];
    fprintf('%s\n',ofile)
    % find last valid data point in the profile
    ilast=find(~isnan(Zyq(yq,:)),1,'last');
    fid=fopen(ofile,'w');
    for m=1:ilast
        %fprintf(fid,'%i %8.3f\n',n-1,Zyq(yq,n));
        %for m=1:numel(X1Dt)
        fprintf(fid,'%i %10.1f %10.1f %13.7f %13.7f %8.2f\n',...
            X1Dt(m),Xutm(m),Yutm(m),Lon(m),Lat(m),Zyq(yq,m));
        %end
    end
    fclose(fid);
end




    % make a new figure
    close all
    p=[];
    figure('position',[ 86         102        1238         628]);
    % profile axes
    ax1=axes('position',[0.05 .1 .45 .8]);
    % variables to keep track of the min-max cross-range of the profiles
    %  where one of more quarter has data
    dmax=0;dmin=Inf;

    
  % loop through 5 quarters starting with 5th quarter of previous year
  nq=0;np=0; % quarter and valid profile counters  
  for q5=[3 4 1 2 3]

      q=4;if q5 < 5;q=q5;end
      nq=nq+1;

      % figure out the survey year based on the beach year and quarter
      if q5 == 3 || q5 == 4; Syear = Byear -1; else Syear = Byear; end

      % name of the quarterly data file create by BuildQuarterlyProfiles.m
        pfile=[bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(Syear) '_' num2str(q) '.dat'];
      % plot if it exists  
        if exist(pfile,'file')
          fprintf('Plotting %s\n',pfile);
          d=load(pfile); 
          np=np+1;
          % plot the profile
          if nq < 5
          p(np)=plot(d(:,1),d(:,6),'-','color',col(nq,:),'linewidth',3,...
              'DisplayName',DispName(nq,:));
          else
          p(np)=plot(d(:,1),d(:,6),'m:','color',col(nq,:),'linewidth',3,...
             'DisplayName',DispName(nq,:));
          end
          hold on;
          %  check for possible new min or max xshore data limits
          if d(1,1) < dmin;dmin=d(1,1);x1=d(1,4);y1=d(1,5);end
          if d(end,1) > dmax;dmax=d(end,1);x2=d(end,4);y2=d(end,5);end
        end

     end
 

  % set x axes limits base don valida data range
  set(ax1,'xlim',[dmin-5 dmax+5])
  % add some tide levels
  xl=get(gca,'xlim');
  %plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
  plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
  %plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
  plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
  plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14); 

  set(gca,'xdir','reverse','fontsize',16,'linewidth',2);
  xlabel('Cross-Shore Distance from MOP Back Beach Point (m)');
  ylabel('Elevation (m, NAVD88)');
  legend(p,'location','northwest','fontsize',16)
  title('Quarterly Mean Profiles','fontsize',20)

  %% On a second axes, plot portion of Mop transect with valida quarterly
  %   profile data on a google map

  ax2=axes('position',[0.58 .1 .38 .8]);
  
  % plot the transect range with valid quarterly data
  p1=plot([x1 x2],[y1 y2],'-','color','m','linewidth',5,'DisplayName',...
      'Subaerial Transect with Data');
  hold on;
  % mark the location of MOP back beach point
  p2=plot(Mop.BackLon(MopNum),Mop.BackLat(MopNum),'o','markersize',15,...
      'markeredgecolor','m','markerfacecolor','k','linewidth',2,...
      'DisplayName','MOP Back Beach Point'); 
  % text(mean(d(:,4)),mean(d(:,5))+0.00005,num2str(d(end,1)-d(1,1),'%im'),...
  %     'horizontalalign','center','fontsize',20,'fontweight','bold','color','m');
  set(gca,'xlim',[mean([x1 x2])-.001 mean([x1 x2])+.001],...
      'ylim',[mean([y1 y2])-.001 mean([y1 y2])+.001])
  set(gca,'fontsize',16,'linewidth',2);
  xlabel('Longitude');
  ylabel('Latitude');
  
  xl=get(gca,'xlim');yl=get(gca,'ylim');
  % ds=0.25;
  % set(gca,'xlim',[xl(1)-ds*diff(xl) xl(2)+ds*diff(xl)]);
  % set(gca,'ylim',[yl(1)-ds*diff(yl) yl(2)+ds*diff(yl)]);
  plot_google_map('MapType', 'satellite');
  legend([p1 p2],'location','northwest','fontsize',16)
   title(['MOP D' num2str(MopNum,'%4.4i')],'fontsize',20)

   ax3=axes('position',[.52 .92 .01 .01],'xtick',[],'ytick',[]);
   title(btitle{sb},'fontsize',24)
   ax3.YAxis.Visible = 'off';   % remove y-axis
   ax3.XAxis.Visible = 'off';   % remove x-axis

   % make a jpeg file (much smaller than png equivalent)

    jpgfile=[bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(Syear) '.jpg'];

   print(gcf,'-djpeg',jpgfile,'-r72');

 end
   
end
