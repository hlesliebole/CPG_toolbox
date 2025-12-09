
x=[];y=[];z=[];
load M00578SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00579SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00580SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00581SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00582SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00583SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
load M00584SG.mat
n=size(SG,2);
% n=find([SG.Datenum]==datenum(2021,10,28) &...
%     strcmpi({SG.Source},'iG8wheel') == 1);
x=[x' SG(n).X']';
y=[y' SG(n).Y']';
z=[z' SG(n).Z']';
figure('position',[191   261   846   491]);
ScatterPlotBeachUTM(vertcat(x),vertcat(y),vertcat(z),'3d')
grid on;
%set(gca,'position',[0.1202    -0.05    0.7164    0.8093]);
view(-90,90)
%zlabel('Elevation (m, NAVD88)')
BeachColorbar
set(gca,'fontsize',14)
xl=get(gca,'xlim');yl=get(gca,'ylim');
for mop=578:584
PlotMopTransectUTM(mop,'3D','k','ShadeOff');
end
set(gca,'xlim',[xl(1) xl(2)+20],'ylim',[yl(1)-20 yl(2)+20]);

title([ datestr(SG(end).Datenum,'mm/dd/yyyy') " Gridded Wheelie Survey"])
makepng('LatestWheelieGrid.png')
% title([ datestr(SG(n).Datenum,'mm/dd/yyyy') " Gridded Wheelie Survey"])
% makepng('28oct2021WheelieGrid.png')

%---------------------------------------------------------------------

function PlotMopTransectUTM(MopID,PlotMethod,LineColor,ShadeOnOff)

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
     zd=10.;
     
   % option to shade Mop area being viewed on 3D plot
   %   makes shaded area 50% transparent

  if strcmpi(ShadeOnOff,'ShadeOn') && MopNumber > 2
    x1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
    y1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
    x2=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
    y2=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
    x3=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
    y3=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
    x4=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
    y4=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);
    fill3([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[zd zd zd zd zd],[.95 .95 .95],...
        'edgecolor','none','facealpha',.2)
    hold on;
  end
       
    % 3D line plot
    plot3([Mop.BackXutm(n) Mop.OffXutm(n)],[Mop.BackYutm(n) Mop.OffYutm(n)],...
        [zd zd],'-','Color',LineColor);
    hold on;
    plot3([Mop.BackXutm(n) Mop.OffXutm(n)],[Mop.BackYutm(n) Mop.OffYutm(n) ],...
        [zd zd],'.','Color',LineColor,'markersize',25);
    zlabel('NAVD88 (m)');
    hold on;
    text(Mop.BackXutm(n)+10 ,Mop.BackYutm(n) ,zd,...
        num2str(MopID))
  
 else  % 2D plot
     
% option to shade Mop area being viewed on 2D plot

  if strcmpi(ShadeOnOff,'ShadeOn') && MopNumber > 2
    x1=mean([Mop.BackXutm(MopNumber-1) Mop.BackXutm(MopNumber)]);
    y1=mean([Mop.BackYutm(MopNumber-1) Mop.BackYutm(MopNumber)]);
    x2=mean([Mop.BackXutm(MopNumber) Mop.BackXutm(MopNumber+1)]);
    y2=mean([Mop.BackYutm(MopNumber) Mop.BackYutm(MopNumber+1)]);
    x3=mean([Mop.OffXutm(MopNumber) Mop.OffXutm(MopNumber+1)]);
    y3=mean([Mop.OffYutm(MopNumber) Mop.OffYutm(MopNumber+1)]);
    x4=mean([Mop.OffXutm(MopNumber-1) Mop.OffXutm(MopNumber)]);
    y4=mean([Mop.OffYutm(MopNumber-1) Mop.OffYutm(MopNumber)]);
    fill([x1 x2 x3 x4 x1],[y1 y2 y3 y4 y1],[.95 .95 .95],'edgecolor','none')
    hold on;
  end
    
    % 2D line plot
    plot([Mop.BackXutm(n) Mop.OffXutm(n)],[Mop.BackYutm(n) Mop.OffYutm(n)],...
        '-','Color',LineColor);
    hold on;
    plot([Mop.BackXutm(n) Mop.OffXutm(n)],[Mop.BackYutm(n) Mop.OffYutm(n) ],...
        '.','Color',LineColor,'markersize',25);
    
 end
 
end
