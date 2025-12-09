%  make a nearshore-beach-bluff colormap


%zc = [-8:0.25:-1,-0.8,-0.631,0.218,0.774,1.344,1.566,2.119,2.4:.2:5,5.25:.25:24]; % 127 elevation levels
zc = [-8:0.25:-1.25,-0.932,-0.631,-0.058,0.218,0.774,1.344,1.566,2.119,2.4:.2:5]; % navd88
zlabel = [-8:1:-1,-0.63,-0.06,0.22,.774,1.34,1.57,2.12,3:1:5]; % ticks to label

tz=1:length(zc);
%zc=zc-0.774;

close all
clf;pcolor(1:2,zc,[zc ; zc]');shading flat
set(gca,'ytick',zlabel,'yticklabel',num2str(zlabel'))
colormap(demcmap(zc-.774,length(zc)-1))
cols=colormap;land=find(cols(:,2) == 0.4 & cols(:,3) == 0.2);
cols(end-3,:)=1.1*cols(end-2,:);
cols(end-2,:)=mean([cols(end-1,:)' cols(end-2,:)']');
cols(end,:)=.9*cols(end,:);
colormap([cols(1:land-1,:)' flipud(cols(land:end,:))']');


seamap=parula(14+land-1);%sea=sea(end-land+2:end,:);
upper=gray(18);
upper=flipud(parula(58));
colormap([seamap(1:end-14,:)' flipud(cols(land+11:end,:))' flipud(upper(2:13,:))']');




yl=get(gca,'ylim');

yyaxis right;set(gca,'ylim',yl,'ytick',...
    [-8+.774,-7+.774,-6+.774,-5+.774,-4+.774,-3+.774,-2+.774,...
    -0.63,-0.06,0.22,0.774,1.34,1.63,2.11,2+.774,3+.774,4+.774],'yticklabel',...
    str2mat('-8','-7','-6','-5','-4','-3','-2','LAT','MLLW','MLW',...
    'MSL','MHW','MHHW','HAT','+2','+3','+4'))

title({'Elev. (m) ';'NAVD88 |  MSL ';''},'fontweight','normal','fontsize',10)
set(gca,'xtick',[]);
%xlabel({'Elev. (m) ';'NAVD88 |  MSL '},'fontweight','normal','fontsize',10)

%figure('position',[440    96   316   702]); 
figure

seamap=parula(14+land-1);%sea=sea(end-land+2:end,:);
%upper=gray(18);
upper=flipud(parula(84));
% upper(1:15,2)=upper(1:15,2)-.5;
% upper(1:15,2)=upper(1:15,2)-.2;
%upper(1:15,3)=upper(1:15,3)-0.05;

colormap([seamap(1:end-15,:)' 1.4*flipud(cols(land+13:end,:))' flipud(upper(1:15,:))']');

izx = interp1(zc,1:length(zc),-8:0.01:5,'linear'); % here we transform the data
BeachColorMap=colormap;
BeachColorMap=BeachColorMap(floor(izx),:);
ic=find(round(10000*BeachColorMap(:,1)) == 5307);
BeachColorMap(ic,:)=BeachColorMap(ic,:)*.65;
ic=find(round(10000*BeachColorMap(:,1)) == 4507);
BeachColorMap(ic,:)=BeachColorMap(ic,:)*.75;
ic=find(round(10000*BeachColorMap(:,1)) == 3765);
BeachColorMap(ic,:)=BeachColorMap(ic,:)*.85;

BeachColorBarLims=[-8 5];
set(gca,'clim',BeachColorBarLims); BeachColorBar=colorbar; colormap(BeachColorMap)
BeachColorBarTitleString={'Elev. (m)';'MSL | NAVD88'};
BeachColorBar.Title.String=BeachColorBarTitleString;
BeachColorBarTicks=[-8+.774,-7+.774,-6+.774,-5+.774,-4+.774,-3+.774,-2+.774,...
    -0.63,-0.06,0.22,0.774,1.34,1.566,2.119,2+.774,3+.774,4+.774];
BeachColorBar.Ticks=BeachColorBarTicks;
BeachColorBarTickLabels=...
    str2mat('-8','-7','-6','-5','-4','-3','-2','LAT','MLLW','MLW',...
    'MSL','MHW','MHHW','HAT','+2','+3','+4');
BeachColorBar.TickLabels=BeachColorBarTickLabels;
BeachColorBar.TickLength=0;
NavdAxes=axes('position',BeachColorBar.Position,'color','none','xtick',[],...
    'ytick',[-8:5],'ylim',[-8 5],'TickLength',[0 .01]);

save BeachColorMap.mat BeachColorMap BeachColorBarLims...
    BeachColorBarTitleString BeachColorBarTicks BeachColorBarTickLabels

% d = 1.25*peaks(123);
% dtr = d;
% dtr(:) = interp1(zc,tz,d(:),'pchip'); % here we transform the data - only for displaying...
% subplot(1,2,1)
% imagesc(d)
% colorbar
% subplot(1,2,2)
% imagesc(dtr)         
% cblh = colorbar;        % and here we have to pay the price of being cunning/lazy
% set(cblh,'ytick',1:length(zc))  % and set ticks and ticklabels manually
% set(cblh,'ytick',1:length(zc),'yticklabel',num2str(zc'))