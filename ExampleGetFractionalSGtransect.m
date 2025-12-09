% Example code to get fractional Mop transect lines
% to 1 decimal place (eg. for 583.X) from the CPG MOP 
% SG gridded survey data files

% set path to CPG MOP files from your machine
% addpath /Volumes/group/MOPS  % folder with MOP mat files
% addpath /Volumes/group/MOPS/toolbox  % folder with MOP m-scripts

% set fractional mop number
FracMopNumber=583.3;
% make sure it is to one decimal place
FracMopNumber=round(FracMopNumber,1);

% get the Mop line number that is closest by rounding
% to whole mop number
MopNumber=round(FracMopNumber);

% The mop area is going to be divided into 11 subtransects, 0.1 "mops"
% apart, centered on the Mop transect line. So, eg for MopNumber=583, 
% subtransect 1 is 582.5, subtransect 6 = 583.0, and subtransect 11=583.5

% figure out the desired subtransect number
Frac=FracMopNumber-MopNumber;
SubTransectNumber=6+round(10*Frac);
  
% load the SG mat file
matfile=sprintf('M%5.5iSG.mat',MopNumber);
load(matfile,'SG');

%  load Mop Transect Info
load('MopTableUTM.mat','Mop');

% divide mop area into 11 mop subtransects at 1m xshore resolution,
%  with an extra 100m of back beach for each. x1d is the xshore distance
%  from the backbeach line (line connecting Mop backbeach points).
%  xt,yt are the 2d grid locations of the main Mop transect line, and
%  xst,yst are the 2d grid locations for all the subtransects (including
%  the main line)
[x1d,xt,yt,xst,yst]=GetTransectLines(Mop,MopNumber,11,[-100 0]);
    
% choose a survey SG struct array index number to extract the transect from 

% wheelie survey SG struct arrray index numbers
idx=find(strcmp({SG.Source},'iG8wheel'));
% use the most recent wheelie survey index number
sn=idx(end);

% The SG struct array only retains info for grid points that
%  have valid data, so reconstruct the full x,y grid "tg" with the
%  no data grid points as NaN's

xg=SG(sn).X;yg=SG(sn).Y;zg=SG(sn).Z;
nx=max(xg)-min(xg)+3;ny=max(yg)-min(yg)+3; % temp grid dimensions
tg=nan(ny,nx); % temp grid of Nans
idx=sub2ind([ny nx],yg-min(yg)+2,xg-min(xg)+2); % data point 1d indices
tg(idx)=zg; % add data to temp grid of NaNs
[X,Y]=meshgrid(min(xg)-1:max(xg)+1,min(yg)-1:max(yg)+1); % define grid x,y points

% now 2d interpolate z values for all the subtransect points from the
zst=xst*NaN; % initialize transect elevation points as NaNs 
zst(:) = interp2(X,Y,tg,xst(:),yst(:)); % interpolate onto subtransects

% The subtransect info you want is
%  x1d,zst(SubTransectNumber,:)
xFracTransect=x1d;
zFracTransect=zst(SubTransectNumber,:);

% plot the section of the subtransect with data
figure;
igood=find(~isnan(zFracTransect));
plot(x1d(igood),zFracTransect(igood),'.-');
set(gca,'xdir','reverse');
grid on;
xlabel('Xshore Distance from Mop Back Beach Line (m)');
ylabel('Elevation (m, NAVD88)')
title([{matfile};{datestr(SG(sn).Datenum)};...
{['Fractional Mop Transect:' num2str(FracMopNumber,'%6.1f')]} ]);

