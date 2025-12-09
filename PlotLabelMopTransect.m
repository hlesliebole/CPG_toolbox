function PlotLabelMopTransect(MopID,PlotMethod,LineColor,ShadeOnOff)

% used by PlotCMstruct.m
hold on;
load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% Get numeric mop number
if isnumeric(MopID);MopNumber=MopID;else;...
        MopNumber=find(contains(Mop.Name,MopID));end

% plot transect
  n=MopNumber;

 if contains(PlotMethod,'3d','IgnoreCase',true) % 3D plot
     
     zd=.774; % plot 3D line and shading at msl elevation
     
   % option to shade Mop area being viewed on 3D plot
   %   makes shaded area 50% transparent

  if strcmpi(ShadeOnOff,'ShadeOn') && MopNumber > 2
    x1=mean([Mop.BackLon(MopNumber-1) Mop.BackLon(MopNumber)]);
    y1=mean([Mop.BackLat(MopNumber-1) Mop.BackLat(MopNumber)]);
    x2=mean([Mop.BackLon(MopNumber) Mop.BackLon(MopNumber+1)]);
    y2=mean([Mop.BackLat(MopNumber) Mop.BackLat(MopNumber+1)]);
    x3=mean([Mop.OffLon(MopNumber) Mop.OffLon(MopNumber+1)]);
    y3=mean([Mop.OffLat(MopNumber) Mop.OffLat(MopNumber+1)]);
    x4=mean([Mop.OffLon(MopNumber-1) Mop.OffLon(MopNumber)]);
    y4=mean([Mop.OffLat(MopNumber-1) Mop.OffLat(MopNumber)]);
    fill3([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[zd zd zd zd zd],[.95 .95 .95],...
        'edgecolor','none','facealpha',.2)
    hold on;
  end
       
    % 3D line plot
    plot3([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],...
        [zd zd],'-','Color',LineColor);
    hold on;
    plot3([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n) ],...
        [zd zd],'.','Color',LineColor,'markersize',25);
    zlabel('NAVD88 (m)');
    hold on;
    text(Mop.OffLon(n),Mop.OffLat(n),zd,[num2str(MopNumber) '  '],...
        'horizontalalign','right','fontsize',14,'fontweight','bold',...
        'color',LineColor)
  
 else  % 2D plot
     
% option to shade Mop area being viewed on 2D plot

  if strcmpi(ShadeOnOff,'ShadeOn') && MopNumber > 2
    x1=mean([Mop.BackLon(MopNumber-1) Mop.BackLon(MopNumber)]);
    y1=mean([Mop.BackLat(MopNumber-1) Mop.BackLat(MopNumber)]);
    x2=mean([Mop.BackLon(MopNumber) Mop.BackLon(MopNumber+1)]);
    y2=mean([Mop.BackLat(MopNumber) Mop.BackLat(MopNumber+1)]);
    x3=mean([Mop.OffLon(MopNumber) Mop.OffLon(MopNumber+1)]);
    y3=mean([Mop.OffLat(MopNumber) Mop.OffLat(MopNumber+1)]);
    x4=mean([Mop.OffLon(MopNumber-1) Mop.OffLon(MopNumber)]);
    y4=mean([Mop.OffLat(MopNumber-1) Mop.OffLat(MopNumber)]);
    fill([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[.95 .95 .95],'edgecolor','none')
    hold on;
  end
    
    % 2D line plot
    plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],...
        '-','Color',LineColor);
    hold on;
    plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n) ],...
        '.','Color',LineColor,'markersize',25);
    text(Mop.OffLon(n),Mop.OffLat(n),[num2str(MopNumber) '  '],...
        'horizontalalign','right','fontsize',14,'fontweight','bold',...
        'color',LineColor)
 end
 
end

 
