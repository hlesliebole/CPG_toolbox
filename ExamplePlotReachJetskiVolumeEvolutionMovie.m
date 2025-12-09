% Example code plot the jumbo survey cross-shore volume evolution for
% a specified reach (Mop range) and date range (changes plotted 
% relative to the first survey date in the range).

% Jumbo surveys have been gridded into 1m x,y (UTM) spatial resolution
% elevation points (m, NAVD88) and stored in the CPGMOP M*SG.mat files 
% for each Mop number.  So he elevation difference at a grid point from 
% two surveys is equal to the volume difference in m^3.


%% Uses the m-script function CG=CombineSGdata(mapth,MopStart,MopEnd)

%% settings
% 1. Path to CPGMOP data files
%mpath='/volumes/group/MOPS/'; % reefbreak on a mac
mpath='/Users/William/Desktop/MOPS/';

% 2. Mop range to use in volume change calculations
MopStart=668;MopEnd=682; % start and mop numbers for north TP 666-683
%MopStart=568;MopEnd=598; % start and mop numbers for north TP

% 3. Date range to consider. The first and last dates to not
%    have to match survey dates precisely
StartDate=datenum(2010,11,22);
EndDate=datenum(2011,3,10);

% 4. Elevation bin sizes (m) to use when calculating volume change as 
%    a function of the starting survey grid elevations. Larger
%    values give smoother results with less cross-shore detail.  
zRes=0.5; % 10 cm

% Alongshore reach length (m). 
% Used to normalize first moment of the seaward deposition distribution
L=100*(MopEnd-MopStart+1);

%% load the combined SG gridded mat file data.  Return the combined
%  data to a struct array SG instead of CG to use the normal 
%  single mop SG gridding and plotting code with the combined data.
%SG=CombineSGdata(mpath,MopStart,MopEnd);

%% identify jumbos with jetskis by checking the original survey file names for the 
%  word jumbo
jumbo=find(contains({SG.File},'umbo'));
% find jumbo surveys with jetski data
m=0;
jetski=[];
for j=1:length(jumbo)
    %fprintf('%s %5.1f\n',datestr(SM(jumbo(j)).Datenum),min(SM(jumbo(j)).Z1Dmean));
    if( min(SG(jumbo(j)).Z) < -3 )
        m=m+1;
        jetski(m)=jumbo(j);
    end       
end

fprintf('The SG struct array has %i Jumbo-Jetski Surveys.\n',numel(jetski))

% find the jetski surveys that fall in the date range
idx=find( [SG(jetski).Datenum] >= StartDate & [SG(jetski).Datenum] <= EndDate);

fprintf('%i Jumbo-Jetski Surveys found in the date range.\n',numel(idx))
% reduce indices to just these survey indices
jetski=jetski(idx);
datestr([SG(jetski).Datenum])

%% Need to turn the saved SG struct array grid points into a 2d grid arrays

% make base grid area encompassing all gridded survey data
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));

% 2d UTM grid indice matrices
[X,Y]=meshgrid(minx:maxx,miny:maxy);

%% Z0 is a grid of the first survey in the date range. All grid
%  changes will be calculated relative to this reference grid.
SurvNum=jetski(1); 
% Mop area 1m gridded points with valid data for this survey
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% get 1d indice values of the valid x,y grid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
% initialize the 2d elevation Z0 array as NaNs
Z0=X*NaN; 
% overlay the valid data point elevations using the 1d indices
Z0(idx)=SG(SurvNum).Z; 

% rounded Z0 elevations relative to MSL and scaled into elevation bins  
z0bin=round((Z0-0.774)/zRes); 



figure('position',[ 21          56        1385         728]);

%% add wave Hs time series as 2nd axes

% use Mop in middle of Cardiff for wave info
MopNumber=675;
stn=['D0' num2str(MopNumber)];

urlbase=...  % nearshore station
'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/model/MOP_alongshore/';

% read hindcast file dates and wave height info
urlend = '_hindcast.nc'; 
dsurl=strcat(urlbase,stn,urlend); % full url path and filename
wavetime=ncread(dsurl,'waveTime'); % read wave times
% convert mop unix time to UTC
hctime=datetime(wavetime,'ConvertFrom','posixTime');
hchs=ncread(dsurl,'waveHs'); % read wave heights

ax2=axes('position',[0.1300    0.1100    0.7750 0.2150]);

% plot time series
plot(hctime,hchs,'k-','linewidth',2); hold on;

% just show date range of jetski surveys
set(ax2,'xlim',[datetime(StartDate-7,'convertfrom','datenum') ...
datetime(EndDate+7,'convertfrom','datenum')]);

ax2ylims=get(gca,'ylim');



% replot time series
%plot(hctime,hchs,'k-','linewidth',2); hold on;

ylabel('Hs (m)')
xlabel('Date')
set(gca,'fontsize',14);
grid on;
datetick('x','mm/dd')


%% Loop through other surveys in the date range and compare their grid
% elevations changes (volume changes) to Z0 as a function of the Z0 

ax1=axes('position',[0.1300    0.4100    0.7750 0.5150]);
set(gca,'xtick',min(z0bin(:)):max(z0bin(:))*zRes,'xlim',[min(z0bin(:)) max(z0bin(:))]*zRes);
xlabel([ datestr(SG(jetski(1)).Datenum) ' Nearshore Elevation, Z_o (m, MSL)']);
ylabel('Net Volume Change Density (m^{3} x 1000 / dZ_o )');
set(gca,'fontsize',14);
title(['Net Cross-shore Volume Change since ' datestr(SG(jetski(1)).Datenum)...
    ' : Mops ' num2str(MopStart) ' to ' num2str(MopEnd)],'fontsize',18)
grid on;box on;hold on;
set(ax1,'ylim',[-50 50]);


col=jet(length(jetski)-1);
for n=2:length(jetski)
    
   SurvNum=jetski(n); 
    % Mop area 1m gridded points with valid data for this survey
    x=SG(SurvNum).X;
    y=SG(SurvNum).Y;
    % get 1d indice values of the valid x,y grid points
    idx=sub2ind(size(X),y-miny+1,x-minx+1);
    % initialize the 2d elevation Z0 array as NaNs
    Z=X*NaN; 
    % overlay the valid data point elevations using the 1d indices
    Z(idx)=SG(SurvNum).Z; 
    
    % change grid relative to initial grid
    dz=Z-Z0;

    dv=[];
    m=0;
    mu=0; % initialize first moment parameter
    % loop through the Z0 reference grid elevation bins
    for iz=min(z0bin(:)):max(z0bin(:))
      m=m+1;
      % net volume change "density" in the elevation bin
      dv(m)=sum(dz(z0bin==iz),'omitnan')/zRes;
      % add to first moment
      if(iz < 0 && dv(m) > 0) mu=mu+dv(m)*iz*zRes;end
    end
    % normalize by mop range alongshore distance, so results for
    % different coastal reaches can be compared
    mu=round(-mu/L); 
    % plot results
    fprintf('%s\n',[datestr(SG(jetski(n)).Datenum) '  \mu_{SSD} = ' num2str(mu)])
    if n > 2
        set(pl(n-1),'markersize',15,'linewidth',2)
    end
    pl(n)=plot(ax1,[min(z0bin(:)):max(z0bin(:))]*zRes,dv/1000,'.-','color',col(n-1,:),...
    'markersize',25,'DisplayName',...
    [datestr(SG(jetski(n)).Datenum) '  \mu_{SSD} = ' num2str(mu)] ,'linewidth',4);hold on

    % mark dates of surveys 
    %for n=2:length(jetski)
        x=datetime(SG(jetski(n)).Datenum,'convertfrom','datenum');
        plot(ax2,[x x],ax2ylims,'-','color',col(n-1,:),'linewidth',2); hold on;
    %end
%     if n < 10
%         legend('numcolumns',1)
%     else
        legend('numcolumns',2)
 %   end
    
    if n == 2
       
    gif('CardiffVolumeChange.gif','overwrite',true,'DelayTime',1.)
   else
    gif('frame',gcf)
   end
    
end

% pause end of loop
for n=1:2
    gif('frame',gcf)
end

% make an image of the final plot
makepng('CardiffFrequentJetskiEvolution.png')

fprintf('%s\n% s\n %s\n','Output image and gif movie written to:',...
    'CardiffFrequentJetskiEvolution.png',...
    'CardiffVolumeChange.gif')

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


