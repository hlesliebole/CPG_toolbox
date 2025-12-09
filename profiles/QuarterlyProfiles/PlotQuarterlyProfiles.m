
% Makes quarterly prifle plots for an defined beach year

close all
clearvars

Byear=2024; % Survey beach year to plot quarterly means

% get MOP transect location info
load MOPTableUTM.mat

% open mop prolile definition file
fid=fopen('MopLatLonEutnNutmOrientation.dat','w');

% define 9 state beach names for files
bname{1}='BorderField';
bname{2}='SilverStrand';
bname{3}='TorreyPines';
bname{4}='Cardiff';
bname{5}='SanElijo';
bname{6}='Moonlight';
bname{7}='Leucadia';
bname{8}='SouthCarlsbad';
bname{9}='Carlsbad';

% define 9 state beach names for titles
btitle{1}='Border Field SP';
btitle{2}='Silver Strand SB';
btitle{3}='Torrey Pines SB';
btitle{4}='Cardiff SB';
btitle{5}='San Elijo SB';
btitle{6}='Moonlight SB';
btitle{7}='Leucadia SB';
btitle{8}='South Carlsbad SB';
btitle{9}='Carlsbad SB';

% define 9 mop ranges for the state beaches
MopRng=[2 28 
        85 157
        536 607
        664 683
        684 706
        717 725
        740 757
        762 819
        825 854];

% legend names for 4 profiles
DispName(1,:)=[num2str(Byear-1) ' 4th Quarter (OND)' ];
DispName(2,:)=[num2str(Byear) ' 1st Quarter (JFM)' ];
DispName(3,:)=[num2str(Byear) ' 2nd Quarter (AMJ)' ];
DispName(4,:)=[num2str(Byear) ' 3rd Quarter (JAS)' ];

% line colors for 4 profiles
col(1,:)=[0.7 0 0];
col(2,:)=[0 0 0.7];
col(3,:)=[0.7 0.6 0];
col(4,:)=[0 0.7 0];


%% loop thru state beaches

for sb=2:9

 % loop through mop range for this beach
 for MopNum=MopRng(sb,1):MopRng(sb,2)

    % make a new figure
    close all
    p=[];
    figure('position',[ 86         102        1238         628]);
    % profile axes
    ax1=axes('position',[0.05 .1 .45 .8]);
    % variables to keep track of the min-max cross-range of the profiles
    %  where one of more quarter has data
    dmax=0;dmin=Inf;

    
  % loop through 4 quarters starting with 4th quarter of previous year
  nq=0;np=0; % quarter and valid profile counters  
  for q=[4 1 2 3]
      nq=nq+1;

      % figure out the survey year based on the beach year and quarter
      if q ==4; Syear = Byear -1; else Syear = Byear; end

      % name of the quarterly data file create by BuildQuarterlyProfiles.m
        pfile=[bname{sb} '/' bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(Syear) '_' num2str(q) '.dat'];
      % plot if it exists  
        if exist(pfile,'file')
          fprintf('Plotting %s\n',pfile);
          d=load(pfile); 
          np=np+1;
          % plot the profile
          p(np)=plot(d(:,1),d(:,6),'-','color',col(nq,:),'linewidth',3,...
              'DisplayName',DispName(nq,:));
          hold on;
          %  check for possible new min or max xshore data limits
          if d(1,1) < dmin;dmin=d(1,1);x1=d(1,4);y1=d(1,5);end
          if d(end,1) > dmax;dmax=d(end,1);x2=d(end,4);y2=d(end,5);end
        end

     end
 

  % set x axes limits based on valid data range
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

    jpgfile=['Plots/' bname{sb} '_D' num2str(MopNum,'%4.4i') '_'...
        num2str(Syear) '.jpg'];

    print(gcf,'-djpeg',jpgfile,'-r72');

 end
   
end
