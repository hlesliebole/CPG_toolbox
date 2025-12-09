
addpath altmany-export_fig-3

load MopTableUTM.mat

MopStart=find(strcmp([Mop.Name],'D0645'));
MopEnd=find(strcmp([Mop.Name],'D0657'));
MopStart=find(strcmp([Mop.Name],'D0632'));
MopEnd=find(strcmp([Mop.Name],'D0647'));
MopStart=find(strcmp([Mop.Name],'D0636'));
MopEnd=find(strcmp([Mop.Name],'D0665'));

% san dieguito + solana
MopStart=find(strcmp([Mop.Name],'D0630'));
%MopEnd=find(strcmp([Mop.Name],'D0646'));
MopEnd=find(strcmp([Mop.Name],'D0667'));

DateStart=datenum(2023,12,27);
DateEnd=datenum(2026,1,1);

CS=SAcombineMops(MopStart,MopEnd);

% solana beach lidar dates
load M00654SA.mat
sdx=find((strcmp({SA.Source},'AtvMR') | strcmp({SA.Source},'Trk')) &...
    [SA.Datenum] >= DateStart & [SA.Datenum] <= DateEnd);

% Reduce to complete surveys between 8 mar and 29 may
gdate(1)=datenum(2024,3,8);
gdate(2)=datenum(2024,3,25);
gdate(3)=datenum(2024,4,4);
gdate(4)=datenum(2024,4,10);
gdate(5)=datenum(2024,4,26);
gdate(6)=datenum(2024,5,1);
gdate(7)=datenum(2024,5,10);
gdate(8)=datenum(2024,5,14);
gdate(9)=datenum(2024,5,24);
gdate(10)=datenum(2024,5,29);

sdx2=find(ismember([SA(sdx).Datenum],gdate));
%%
figure('position',[292   442   594   355])

PrevMon=0;
nn=1;
ns=0;
for surv=sdx(sdx2)
    
    clf

Mon=month(datetime(SA(surv).Datenum,'convertfrom','datenum'));
% use same global monthly mean for all last N survey frames
if Mon ~= PrevMon & nn == 1

%% build Global month mean grid for mop range
PrevMon=Mon;
CGM=GMcombineMopsMonth(MopStart,MopEnd,Mon);

xmin=min(vertcat(CGM.X2D));xmax=max(vertcat(CGM.X2D));
ymin=min(vertcat(CGM.Y2D));ymax=max(vertcat(CGM.Y2D));
xmin=floor(xmin);xmax=ceil(xmax);
ymin=floor(ymin);ymax=ceil(ymax);
x=CGM.X2D;
y=CGM.Y2D;
% x,y grid arrays
[Xg,Yg]=meshgrid(xmin:xmax,ymin:ymax);
idx=sub2ind(size(Xg),y-ymin+1,x-xmin+1);
% Global Month Mean z grid array
Zmean=Xg*NaN; % initialize as NaNs
Zmean(idx)=CGM.Z2Dmean; % assign valid z grid data points
end

Zmon=Xg*NaN;
%mm=mm+1;

ndx=find((strcmp({CS.Source},'AtvMR') | strcmp({CS.Source},'Trk')) &...
    [CS.Datenum] >= SA(surv).Datenum-2 & [CS.Datenum] <= SA(surv).Datenum+2);

for idx = ndx
x=CS(idx).X;
y=CS(idx).Y;
z=CS(idx).Z;
ibad=find(y-ymin+1 < 1 | y-ymin+1 > size(Xg,1));
x(ibad)=[];y(ibad)=[];z(ibad)=[];
ibad=find(x-xmin+1 < 1 | x-xmin+1 > size(Xg,2));
x(ibad)=[];y(ibad)=[];z(ibad)=[];
gdx=sub2ind(size(Xg),y-ymin+1,x-xmin+1);
% overlay gridded recent survey data on grid
Zmon(gdx)=z; % assign valid z grid data points
end

Anom=Zmon-Zmean;
mdx=find(~isnan(Anom(:)));

%figure;surf(Xg,Yg,Zjuly-Zmean);colormap(flipud(polarmap));shading flat;view(2)
%
[lat,lon]=utm2deg(Xg(mdx),Yg(mdx),repmat('11 S',[length(Xg(mdx)) 1]));

ax0=axes;
[ScatterPlot,ColorBarPlot]=ColorScatterPolarmap(lon,lat,Anom(mdx));
delete(ColorBarPlot)
hold on



%
% ndx=find((strcmp({CS.Source},'AtvMR') | strcmp({CS.Source},'Trk')) &...
%     [CS.Datenum] >= DateStart & [CS.Datenum] <= DateEnd);
% 
% if numel(ndx) > 0
% 
% x=[];y=[];z=[];
% for idx = ndx 
% x=[x' vertcat(CS(idx).X)']';
% y=[y' vertcat(CS(idx).Y)']';
% z=[z' vertcat(CS(idx).Z)']';
% end
% % trim higher elevations
% zidx=find(z < 10);
% x=x(zidx);y=y(zidx);z=z(zidx);
% % convert to lat lon
% [lat,lon]=utm2deg(x,y,repmat('11 S',[length(x) 1]));

%hold on;

% ScatterPlotBeachLatLon(lat,lon,z,'2d');
% hold on
% load MopTableUTM.mat
% for n=MopStart:MopEnd
% plot([Mop.BackLon(n) Mop.OffLon(n)],[Mop.BackLat(n) Mop.OffLat(n)],'m-')
% text(Mop.OffLon(n),Mop.OffLat(n),num2str(n),'horizontalalign','right','color','w')
% end
%plot(lon2(idx2),lat2(idx2),'m.')

% set(gca,'clim',[-13 -4])
% set(gca,'clim',[-1 6])
set(gca,'fontsize',16)

% title({[datestr(CS(idx).Datenum) ' |  Multibeam '],...
%     'Mops 645 to 655'},...
%     'fontsize',16);
sdatenum=CS(idx).Datenum;
t1=text(-117.27,32.9777,datestr(CS(idx).Datenum),...
    'fontsize',35,'fontweight','bold','color','k','backgroundcolor',[.8 .8 .8]);
tl=text(-117.275,mean([32.9729   32.9787]),{'Time of Year','  Elevation','Anomaly (m)'},'color','w','fontsize',16,'horizontalalign','left');
% t2=text(-117.2738,32.9915,'Fletcher Cove','color','k','fontsize',18,'fontweight','bold','backgroundcolor',[.8 .8 .8]);
% t3=text(-117.2685,32.975,'San Dieguito Inlet','color','k','fontsize',18,'fontweight','bold','backgroundcolor',[.8 .8 .8]);
% t4=text(-117.2775,33.001,'South Cardiff SB','color','k','fontsize',18,'fontweight','bold','backgroundcolor',[.8 .8 .8]);
set(gca,'xtick',[],'ytick',[]);box on;
set(gca,'xlim',[-117.2762 -117.2651]);
set(gca,'ylim',[32.9729   32.9787]);
dx=.1;set(gca,'position',[0.1300-dx    0.1100-dx    0.7750*1.2    0.8150*1.2])

plot_google_map('MapType', 'satellite')
% if sdatenum == datenum(2024,1,10)
%     Zpre=Zmon;
% end
% 
% if sdatenum > datenum(2024,1,10)
%     Znew=Zmon-Zpre;
%     ns=ns+1;
%     Sdatetime(ns)=datetime(sdatenum,'convertfrom','datenum');
%     np=find(Znew(:) > 0);
%     sum(Znew(np))
%     SposAnom(ns)=sum(Znew(np));
% end


% if sdatenum == datenum(2024,1,29)
%     t5=text(-117.2738,32.9915,{'Nourishment Started','       17-Jan-2024'},...
%         'color','g','fontsize',22,'fontweight','bold','backgroundcolor',[0 0 0]);
% end
% if round(sdatenum) == round(datenum(2024,2,21))
%     t6=text(-117.276,32.9985,{'\Leftarrow Sand Spreads North','        from Pad Edge'},...
%         'color','g','fontsize',22,'fontweight','bold','backgroundcolor',[0 0 0]);
% end
% if sdatenum == datenum(2024,3,8)
%     t7=text(-117.272,32.9827,{'Nourishment',' Completed',' 7-Mar-2024'},...
%         'color','g','fontsize',22,'fontweight','bold','backgroundcolor',[0 0 0]);
% end
% if sdatenum == datenum(2024,4,10)
%     % t8=text(-117.270,32.9775,{'Pad Sand Reaches','Inlet 1 Month Later'},...
%     %     'color','g','fontsize',22,'fontweight','bold','backgroundcolor',[0 0 0]);
% end


%%
ax1=axes('position',[.04 .075 .01 .85]);
colormap(flipud(polarmap(64)))
cb=colorbar;
set(gca,'clim',[-2 2]);
%cb.Label.String='Time of Year Elevation Anomaly (m)';
set(cb,'position',[0.05  0.0747    0.015    0.8505])
%set(NavdAxes,'position',[0.05  0.0747    0.01    0.8505],'fontsize',12)
set(cb,'fontsize',14,'color','w')
%title(BeachColorBarTitleString)
axes(ax1)
axis off
%text(2.75,.5,{'Time of Year',' Elevation','Anomaly (m)'},'color','w','fontsize',16,'horizontalalign','left')
%%
% pos=get(gca,'position');
% set(gca,'position',[pos(1)-0.085 0.05 pos(3)-.12 0.85])
% pos=get(gca,'position');
% fac=0.15;
% dx=pos(3)*fac;dy=pos(4)*fac;set(gca,'position',[pos(1)-dx pos(2)-dy pos(3)+dx pos(4)+dy])

if nn == 1
    gif('SolanaDogBeachNourishmentArrival.gif','overwrite',true,'DelayTime',1.)%,'resolution',72,'nodither')
    nn=nn+1;
% elseif sdatenum == datenum(2024,1,29)
% 
     for n=1:2
       gif('frame',gcf)
     end
%      %delete(t5)
% elseif round(sdatenum) == round(datenum(2024,2,21))
% 
%      for n=1:4
%       gif('frame',gcf)
%      end
%      %delete(t6)
% elseif sdatenum == datenum(2024,3,8)
% 
%      for n=1:4
%       gif('frame',gcf)
%      end
%      %delete(t7)
% elseif sdatenum == datenum(2024,4,10)
% 
%      for n=1:4
%       gif('frame',gcf)
%      end
%      %delete(t7);delete(t8)
else
    gif('frame',gcf)
    nn=nn+1;
    % if mn == 11 || (ismember(yr,[2004 2005 2013 2016]) && mn == 8)
      %   for n=1:3
      % gif('frame',gcf)
      %  end
    %end
end


end

for n=1:2
    gif('frame',gcf)
end

%%

%
%end
%end
% title({'Solana Beach, CA',datestr(CS(idx).Datenum)},'fontsize',18)
% set(gca,'xtick',[],'ytick',[]);box on;
% %%
% pos=get(gca,'position');
% set(gca,'position',[pos(1)-0.085 0.05 pos(3)-.12 0.85])
% BeachColorbar
% %
% ax1=axes('position',[.04 .075 .01 .85]);
% colormap(flipud(polarmap(64)))
% cb=colorbar;
% set(gca,'clim',[-1 1]);
% %cb.Label.String='Time of Year Elevation Anomaly (m)';
% set(cb,'position',[0.05  0.0747    0.01    0.8505])
% %set(NavdAxes,'position',[0.05  0.0747    0.01    0.8505],'fontsize',12)
% set(cb,'fontsize',14)
% %title(BeachColorBarTitleString)
% axes(ax1)
% axis off
% text(1.75,1.05,{'Time of Year','Elevation Anomaly (m)'},'color','k','fontsize',16,'horizontalalign','center')
% %%
% makepng('DogBeachAnomalyEvolutionJuly2024v2.png')
% 
% if nn == 1
%     gif('JumboBlueMovieFast.gif','overwrite',true,'DelayTime',.5)
%     nn=nn+1
%    else
%     gif('frame',gcf)
%     nn=nn+1;
%     if mn == 11 || (ismember(yr,[2004 2005 2013 2016]) && mn == 8)
%         for n=1:3
%         gif('frame',gcf)
%         end
%     end
%    end
% end
% 
% end
% 
% %end
%     for n=1:10
%     gif('frame',gcf)
%     end


 function gif(varargin)
% gif is the simplest way to make gifs. Simply call
% 
%   gif('myfile.gif') 
% 
% to write the first frame, and then call 
% 
%   gif
% 
% to write each subsequent frame. That's it. 
% 
%% Syntax
% 
%  gif('filename.gif') 
%  gif(...,'DelayTime',DelayTimeValue,...) 
%  gif(...,'LoopCount',LoopCountValue,...) 
%  gif(...,'frame',handle,...) 
%  gif(...,'resolution',res)
%  gif(...,'nodither') 
%  gif(...,'overwrite',true)
%  gif 
%  gif('clear') 
% 
%% Description 
% 
% gif('filename.gif') writes the first frame of a new gif file by the name filename.gif. 
% 
% gif(...,'DelayTime',DelayTimeValue,...) specifies a the delay time in seconds between
% frames. Default delay time is 1/15. 
% 
% gif(...,'LoopCount',LoopCountValue,...) specifies the number of times the gif animation 
% will play. Default loop count is Inf. 
% 
% gif(...,'frame',handle,...) uses the frame of the given figure or set of axes. The default 
% frame handle is gcf, meaning the current figure. To turn just one set of axes into a gif, 
% use 'frame',gca. This behavior changed in Jan 2021, as the default option changed from
% gca to gcf.
% 
% gif(...,'resolution',res) specifies the resolution (in dpi) of each frame. This option
% requires export_fig (https://www.mathworks.com/matlabcentral/fileexchange/23629).
%
% gif(...,'nodither') maps each color in the original image to the closest color in the new 
% without dithering. Dithering is performed by default to achieve better color resolution, 
% albeit at the expense of spatial resolution.
% 
% gif(...,'overwrite',true) bypasses a dialoge box that would otherwise verify 
% that you want to overwrite an existing file by the specified name. 
%
% gif adds a frame to the current gif file. 
% 
% gif('clear') clears the persistent variables associated with the most recent gif. 
% 
%% Example 
% For examples, type 
% 
%   cdt gif
% 
%% Author Information 
% This function was written by Chad A. Greene of the University of Texas 
% Institute for Geophysics (UTIG), June 2017. 
% 
% See also: imwrite, getframe, and rgb2ind. 
% Define persistent variables: 
persistent gif_filename firstframe DelayTime DitherOption LoopCount frame resolution
%% Parse Inputs
if nargin>0 
   
   % The user may want to clear things and start over: 
   if any(strcmpi(varargin,'clear'))
            
      % Clear persistent variables associated with this function: 
      clear gif_filename firstframe DelayTime DitherOption LoopCount frame resolution
   end
   
   % If the first input ends in .gif, assume this is the first frame:
   if strcmpi(varargin{1}(end-3:end),'.gif')
      
      % This is what the user wants to call the new .gif file: 
      gif_filename = varargin{1}; 
      
      % Check for an existing .gif file by the same name: 
      if exist(gif_filename,'file')==2
         OverWrite = false; % By default, do NOT overwrite an existing file by the input name. 
         if nargin>1
            tmp = strncmpi(varargin,'overwrite',4); 
            if any(tmp)
               OverWrite = varargin{find(tmp)+1}; 
               assert(islogical(OverWrite),'Error: Overwrite input must be either true or false.')
            end
         end
         
         if ~OverWrite
         
            % Ask the user if (s)he wants to overwrite the existing file: 
            choice = questdlg(['The file  ',gif_filename,' already exists. Overwrite it?'], ...
               'The file already exists.','Overwrite','Cancel','Cancel');
            if strcmp(choice,'Overwrite')
               OverWrite = true; 
            end
         end
         
         % Overwriting basically means deleting and starting from scratch: 
         if OverWrite
            delete(gif_filename) 
         else 
            clear gif_filename firstframe DelayTime DitherOption LoopCount frame
            error('The giffing has been canceled.') 
         end
         
      end
      
      firstframe = true; 
      
      % Set defaults: 
      DelayTime = 1/15; 
      DitherOption = 'dither'; 
      LoopCount = Inf; 
      frame = gcf; 
      resolution = 0; % When 0, it's used as a boolean to say "don't use export_fig". If greater than zero, the boolean says "use export_fig and use the specified resolution."  
   end
   
   tmp = strcmpi(varargin,'DelayTime'); 
   if any(tmp) 
      DelayTime = varargin{find(tmp)+1}; 
      assert(isscalar(DelayTime),'Error: DelayTime must be a scalar value.') 
   end
   
   if any(strcmpi(varargin,'nodither'))
      DitherOption = 'nodither'; 
   end
   
   tmp = strcmpi(varargin,'LoopCount'); 
   if any(tmp) 
      LoopCount = varargin{find(tmp)+1}; 
      assert(isscalar(LoopCount),'Error: LoopCount must be a scalar value.') 
   end
   
   tmp = strncmpi(varargin,'resolution',3); 
   if any(tmp) 
      resolution = varargin{find(tmp)+1}; 
      assert(isscalar(resolution),'Error: resolution must be a scalar value.') 
      assert(exist('export_fig.m','file')==2,'export_fig not found. If you wish to specify the image resolution, get export_fig here :https://www.mathworks.com/matlabcentral/fileexchange/23629. Otherwise remove the resolution from the gif inputs to use the default (lower quality) built-in getframe functionality.')  
      warning off export_fig:exportgraphics
   end
   
   tmp = strcmpi(varargin,'frame'); 
   if any(tmp) 
      frame = varargin{find(tmp)+1}; 
      assert(ishandle(frame)==1,'Error: frame must be a figure handle or axis handle.') 
   end
   
else
   assert(isempty(gif_filename)==0,'Error: The first call of the gif function requires a filename ending in .gif.') 
end
%% Perform work: 
if resolution % If resolution is >0, it means use export_fig
   
   if isgraphics(frame,'figure')
      f = export_fig('-nocrop',['-r',num2str(resolution)]);
   else
      % If the frame is a set of axes instead of a figure, use default cropping: 
      f = export_fig(['-r',num2str(resolution)]);
   end
      
else
   % Get frame: 
   fr = getframe(frame); 
   f =  fr.cdata; 
end
% Convert the frame to a colormap and corresponding indices: 
[imind,cmap] = rgb2ind(f,256,DitherOption);    
% Write the file:     
if firstframe
   imwrite(imind,cmap,gif_filename,'gif','LoopCount',LoopCount,'DelayTime',DelayTime)
   firstframe = false;
else
   imwrite(imind,cmap,gif_filename,'gif','WriteMode','append','DelayTime',DelayTime)
end
 end