clearvars
close all
DefineMopPath

cmon={'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct'};
figure('position',[178  145  1109  429]);

MopNumber=553;
%--------------

load(['M' num2str(MopNumber,'%5.5i') 'SM.mat' ],'SM');

ndx=find(strcmpi({SM.Source},'ig8wheel'));
ndx(4:7)=[];
nn=0;
for n=ndx
    nn=nn+1;
    Z(nn,:)=SM(n).Z1Dtransect;
end

[Xq,Yq] = meshgrid(SM(1).X1D,SM(ndx(1)).Datenum:SM(ndx(end)).Datenum);
Zq=interp2(SM(1).X1D,[SM(ndx).Datenum],Z,Xq,Yq);
T=SM(ndx(1)).Datenum:SM(ndx(end)).Datenum;


set(gca,'xlim',[20 120],'color',[.8 .8 .8]);
xl=get(gca,'xlim');
set(gca,'ylim',[.7 4]);
%yl=get(gca,'ylim');set(gca,'ylim',[0 yl(2)]);
hold on
plot(xl,[2.119 2.119],'k--');text(xl(2),2.26,' HAT','fontsize',14);
plot(xl,[1.566 1.566],'k--');text(xl(2),1.7,' MHHW','fontsize',14);
plot(xl,[1.344 1.344],'k--');text(xl(2),1.44,' MHW','fontsize',14);
plot(xl,[.774 .774],'k--');text(xl(2),.9,' MSL','fontsize',14);
plot(xl,[-0.058 -0.058],'k--');text(xl(2),0.05,' MLLW','fontsize',14);
%ps=plot(77,-0.31,'k.','markersize',20,'DisplayName','Paros');
set(gca,'xdir','reverse','fontsize',14);grid on;box on;
%legend([p pc],'location','eastoutside');
title(['MOP ' num2str(MopNumber) ' Transect Profiles']);
xlabel('Cross-shore Distance (m)');
ylabel('Elevation (m, NAVD88)');

hold on
dt=datetime([SM(ndx(1:end-3)).Datenum],'convertfrom','datenum');
monthcolor=jet(10);
dtmon=month(dt);
nn=0;
nm=0;
mprev=0;

for m=ndx(1:end-3)
    dx=0;
    nn=nn+1;
    %dx=SM(m).Datenum-SM(ndx(1)).Datenum;
    %if ~isnan(min(SM(m).Z1Dmean)) 
    if dtmon(nn) ~= mprev
        mprev=dtmon(nn);
        nm=nm+1;
        lp(nm)=plot([130 140],[0 0],'color',[monthcolor(dtmon(nn),:) .5],...
            'linewidth',3,'displayname',cmon{nm});     
    end
    
    z=SM(m).Z1Dtransect;z(z < 0.7)=NaN;z(z > 3.25)=NaN;
%    z=Zq(nn,:);
    
    p(nn)=plot(SM(1).X1D+dx,Z(nn,:),'-','color','m','linewidth',3);
    %p(nn)=plot(SM(1).X1D+dx,z,'-','color',[0 0 0],'linewidth',2);
    if(nn > 1);set(p(nn-1),'color',[0.3 0.3 0.3],'linewidth',1);end
    if(nn > 2);set(p(nn-2),'color',[0.6 0.6 0.6],'linewidth',1);end
    if(nn > 3);set(p(nn-3),'color',[0.9 0.9 0.9],'linewidth',1);end
    if(nn > 4);set(p(nn-4),'color',[monthcolor(dtmon(nn),:) .5],'linewidth',2);end
    %lg=legend(lp,'fontsize',16,'position',[ 0.2989    0.4709    0.0649    0.4172]);
    lg=legend(lp,'fontsize',20,'location','northwest','numcolumns',9);
     
%     if(nn > 4);set(p(nn-4),'color',[0 0.6 0],'linewidth',1);end
%     if(nn > 5);set(p(nn-5),'color',[0 0  .3],'linewidth',1);end
%     if(nn > 6);set(p(nn-6),'color',[0 0  .6],'linewidth',1);end
%     if(nn > 7);set(p(nn-7),'color',[0.7 0.7  .7],'linewidth',1);end
%     if(nn > 8);set(p(nn-8),'color',[.9 .9 .9],'linewidth',1);end
    
    
    
    set(gca,'xlim',[20 120]);
    set(gca,'ylim',[.7 3.6]);
    title(datestr(SM(m).Datenum),'fontsize',30)
    
   if nn == 1
    gif('BNevolutionColor.gif','overwrite',true,'DelayTime',.1)
   else
    gif('frame',gcf)
   end
    end
%end
    for n=1:20
    gif('frame',gcf)
    end
 

%set(gca,'xlim',[0 xl(2)]);


% 
% if nr > 3
% makepng(['MOP' num2str(MopNumber) 'GBR22profiles.png'])
% else
% makepng(['MOP' num2str(MopNumber) 'GBR22profiles.png'])  
% end
% end  

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


