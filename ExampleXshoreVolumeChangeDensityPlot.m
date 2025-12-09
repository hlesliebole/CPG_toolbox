% Example code to plot the jumbo survey cross-shore volume evolution for
% a specified reach (Mop range) and date range (changes plotted 
% relative to the first survey date in the range).

% Jumbo surveys have been gridded into 1m x,y (UTM) spatial resolution
% elevation points (m, NAVD88) and stored in the CPGMOP M*SG.mat files 
% for each Mop number.  So the elevation difference at a grid point from 
% two surveys is equal to the volume difference in m^3.

% Uses the m-script function CG=CombineSGdata(mapth,MopStart,MopEnd)
%  in group/MOPS/toolbox

%% -------------------------------------------------------------
%% Script Settings

%% 1. Paths to CPGMOP data files and the MOPS/toolbox
mpath='/volumes/group/MOPS/'; % reefbreak on a mac
addpath '/volumes/group/MOPS/toolbox/';

%mpath='/Users/William/Desktop/MOPS/';

%% 2. Mop range to use in volume change calculations

%MopStart=668;MopEnd=682;MopWaves=675; % mop range for Cardiff 
 MopStart=568;MopEnd=598;MopWaves=584; % mop range for north TP

%% 3. Date range to consider. The first and last dates to not
%    have to match survey dates precisely

%StartDate=datenum(2010,11,22);  % Cardiff time with frequent surveys
%EndDate=datenum(2011,3,10);

StartDate=datenum(2021,10,24); % TP time with many frequent (SCARP)
EndDate=datenum(2021,12,18);

%% 4. Elevation bin sizes (m) to use when calculating volume change as 
%    a function of the starting survey grid elevations. Larger
%    values give smoother results with less cross-shore detail.  
zRes=0.5; % 50 cm

%% ---------------------------------------------------------------

% Alongshore reach length (m). 
% Used to normalize first moment of the seaward deposition distribution
L=100*(MopEnd-MopStart+1);

%% load the combined SG gridded mat file data.  Return the combined
%  data to a struct array SG 
SG=CombineSGdata(mpath,MopStart,MopEnd);

%% identify jumbos with jetskis by checking the original survey file names for the 
%  word jumbo
jumbo=find(contains({SG.File},'umbo'));
% Sometimes a survey file labeled jumbo does not have jetski data, so
%  find jumbo surveys with data below -3m navd88 as official jetski jumbos
m=0; % jetski survey counter
jetski=[];
for j=1:length(jumbo)
    if( min(SG(jumbo(j)).Z) < -3 ) % check if min depth below -3m
        m=m+1; jetski(m)=jumbo(j); % if yes, add to SG index list
    end       
end

fprintf('The SG struct array has %i Jumbo-Jetski Surveys.\n',numel(jetski))

% find the jetski surveys that fall in the specified date range
idx=find( [SG(jetski).Datenum] >= StartDate & [SG(jetski).Datenum] <= EndDate);

fprintf('%i Jumbo-Jetski Surveys found in the date range.\n',numel(idx))
% reduce jetski indices to just these survey indices
jetski=jetski(idx);
% print a listing of jetski dates in the date range
datestr([SG(jetski).Datenum])

%% The 1m spatial res gridded survey data saved in the SG struct array 
%  is just the valid SG.X,SG.Y,SG.Z (x,y in UTM coords) data grid points 
%  to save storage space/memory.  So you have to turn the saved SG struct 
%  array grid points back into into actual 2d grid arrays with NaNs for
%  the no data grid points

% Figure out a universal grid area encompassing all gridded survey data
%  in the SG struct array, so everything will be in the same grid ref
%  frame.
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));

% use meshgrid to create universal 2d UTM grid indice matrices
[X,Y]=meshgrid(minx:maxx,miny:maxy);

%% Z0 is a grid of the first survey in the date range. All grid
%  changes will be calculated relative to this reference grid.
SurvNum=jetski(1); % SG index number of first survey in date range
% Mop area 1m gridded points with valid data for this survey
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% get 1d indice values of the valid x,y grid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
% initialize the 2d elevation Z0 array as NaNs
Z0=X*NaN; 
% overlay the valid data point elevations using the 1d indices
Z0(idx)=SG(SurvNum).Z; 

%% Define Z0 elevation bins for the volume change densities using the 
%  specified bin elevation resolut zRes. Make th eZ0 elevations relative 
%  to MSL for easier interpretation of the results.  
z0bin=round((Z0-0.774)/zRes); 

%% create a figure
figure('position',[ 21          56        1385         728]);
set(gcf,'color','w'); % generate images with white backgrounds

%% Get wave Hs time series 

% use MopWaves mop number for wave info

stn=['D0' num2str(MopWaves)];

urlbase=...  % nearshore station
'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';

% read hindcast file dates and wave height info
urlend = '_hindcast.nc'; 
dsurl=strcat(urlbase,stn,urlend); % full url path and filename
hctime=ncread(dsurl,'waveTime'); % read wave times
% convert mop unix time to UTC
hctime=datetime(hctime,'ConvertFrom','posixTime');
hchs=ncread(dsurl,'waveHs'); % read wave heights

% now read nowcast file
urlend = '_nowcast.nc';
dsurl=strcat(urlbase,stn,urlend);% full url path and filename

% read in netcdf wave forecast times 
nctime=ncread(dsurl,'waveTime');
% convert mop unix time to UTC
nctime=datetime(nctime,'ConvertFrom','posixTime');
% combine nowcast with hindcast times
wavetime=vertcat(hctime,nctime); 

% read nowcast Hs data
nchs=ncread(dsurl,'waveHs'); 
% combine with hindcast file
wavehs=vertcat(hchs,nchs); 

%% First plot wave time series as a lower 2nd axes

ax2=axes('position',[0.07500    0.1100    0.750 0.2150]);
set(ax2,'color',[.9 .9 .9]);box on;hold on;
% plot time series
plot(wavetime,wavehs,'k-','linewidth',2); hold on;

% just show date range of jetski surveys
set(ax2,'xlim',[datetime(StartDate-7,'convertfrom','datenum') ...
datetime(EndDate+7,'convertfrom','datenum')]);

ax2ylims=get(ax2,'ylim');

ylabel(['Hs (m), Mop ' num2str(MopWaves)])
xlabel('Date')
set(ax2,'fontsize',14);
grid on;
datetick('x','mm/dd')


%% Loop through surveys in the date range and compare their grid
% elevations changes (aka volume changes with the 1m spatial res grid) 
% to Z0 as a function of the Z0 elevation bands with zRes bin resolution 

ax1=axes('position',[0.07500    0.4100    0.750   0.5150]);
set(ax1,'color',[.9 .9 .9])
set(ax1,'xtick',min(z0bin(:)):max(z0bin(:))*zRes,'xlim',[min(z0bin(:)) max(z0bin(:))]*zRes);
xlabel([ datestr(SG(jetski(1)).Datenum) ' Nearshore Elevation, Z_o (m, MSL)']);
ylabel('Net Volume Change Density (m^{3} x 1000 / m of Z_o )');
set(ax1,'fontsize',16);
title(['Net Cross-shore Volume Change relative to ' datestr(SG(jetski(1)).Datenum)...
    ' : Mops ' num2str(MopStart) ' to ' num2str(MopEnd)],'fontsize',20)
grid on;box on;hold on;

set(ax1,'ylim',[-50 50]);

% show starting survey as grey dashed line
x=datetime(SG(jetski(1)).Datenum,'convertfrom','datenum');
        plot(ax2,[x x],ax2ylims,'--','color',[.6 .6 .6],'linewidth',2); hold on;
plot(ax1,[min(z0bin(:)) max(z0bin(:))]*zRes,[0 0],...
    '--','color',[.6 .6 .6],'linewidth',2,'DisplayName',...
    datestr(SG(jetski(1)).Datenum)); hold on;
        
col=jet(length(jetski)-1);
for n=2:length(jetski)
    
    SurvNum=jetski(n); % survey in SG struct array to process
    
    % Turn stored valid survey grid points into an actual grid
    %  with NaNs at the no data grid points.
    
    % Mop area 1m gridded points with valid data for this survey
    x=SG(SurvNum).X;
    y=SG(SurvNum).Y;
    % get 1d indice values of the valid x,y grid points
    idx=sub2ind(size(X),y-miny+1,x-minx+1);
    % initialize the 2d elevation Z0 array as NaNs
    Z=X*NaN; 
    % overlay the valid data point elevations onto the 2d grid 
    %  of NaNs using the 1d indices
    Z(idx)=SG(SurvNum).Z; 
    
    % get the change grid relative to initial grid
    dZ=Z-Z0;

    % initialize the volume change vector for this survey
    dv(n-1,1:length(min(z0bin(:)):max(z0bin(:))))=0;
    % initialize first moment parameter for this survey
    mu(n-1)=0;
    
    m=0; % Z0 elevation bin counter
    
    % loop through the Z0 reference grid elevation bins
    for iz=min(z0bin(:)):max(z0bin(:))
      m=m+1;
      % net volume change "density" in the elevation bin
      dv(n-1,m)=sum(dZ(z0bin==iz),'omitnan')/zRes;
      % add to first moment mu
      if(iz < 0 && dv(n-1,m) > 0) mu(n-1)=mu(n-1)+dv(n-1,m)*iz*zRes;end
    end
    
    % normalize by mop range alongshore distance, so moment results for
    % different coastal reaches can be compared
    mu(n-1)=round(-mu(n-1)/L); 
    % print result to the command window
    fprintf('%s\n',[datestr(SG(jetski(n)).Datenum) '  \mu_{SSD} = ' num2str(mu(n-1))])
      
end

%% loop through the number of surveys again to make animation

    % set fixed y-axis limits that contain all the surveys
    ymax=10*ceil(0.1*max(dv(:))/1000);
    set(ax1,'ylim',[-ymax ymax]);
    
   for n=2:length(jetski)
     
     % reset previous survey line plot to normal linewidth 
     if n > 2
        set(pl(n-1),'markersize',15,'linewidth',2)
     end
    
    % plot current survey as bold line
    pl(n)=plot(ax1,[min(z0bin(:)):max(z0bin(:))]*zRes,dv(n-1,:)/1000,'.-','color',col(n-1,:),...
    'markersize',25,'DisplayName',...
    [datestr(SG(jetski(n)).Datenum) '  \mu_{SSD} = ' num2str(mu(n-1))] ,'linewidth',4);hold on

    % mark date of survey on wave time series 
        x=datetime(SG(jetski(n)).Datenum,'convertfrom','datenum');
        plot(ax2,[x x],ax2ylims,'-','color',col(n-1,:),'linewidth',2); hold on;

     % redo the legend with new survey date   
     legend('location','eastoutside','numcolumns',1)
     % reset axes location after legend messes it up
     set(ax1,'position',[0.07500    0.4100    0.750   0.5150]);
 
    % create or add to gif animation
    if n == 2
     % long output image filenames with Mop and Date range
     str=['XshoreVolChangeMops' num2str(MopStart) 'to' num2str(MopEnd)...
        'from' datestr(SG(jetski(1)).Datenum,'yyyymmdd')...
        'to' datestr(SG(jetski(end)).Datenum,'yyyymmdd')];
     gif([ str '.gif'],'overwrite',true,'DelayTime',1.)
    else
     gif('frame',gcf)
    end
         
   end

% make pause at end of gif loop by repeating last frame a few times
for n=1:3
    gif('frame',gcf)
end

% make an static image of the final plot
set(pl(length(jetski)),'markersize',15,'linewidth',2); % unbold last survey first
% make static png image 
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpng','-r300','-loose',[ str '.png']);

fprintf('\n%s\n%s\n%s\n','Output image and gif movie written to:',...
    [ str '.png'],...
    [ str '.gif'])

%-----------------------------------------------------------------------
%-----------------------------------------------------------------------

%% function to make gifs

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


