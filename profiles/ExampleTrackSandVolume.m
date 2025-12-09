% Example code to make a plot of the 1m gridded data from a Jumbo
%  survey for a combined range of mop SG mat file in a new CG struct array.

%% Uses the m-script function CG=CombineSGdata(mapth,MopStart,MopEnd)

%% settings
%mpath='/volumes/group/MOPS/'; % path to cpg mop files
mpath='/Users/William/Desktop/MOPS/';
MopStart=624%666;%582; % start mop number
MopEnd=627%683;%590; % ending mop number

%% load the combined SG mat file data.  Return the combined
%  data to a struct array SG instead of CG to use the normal 
%  single mop SG gridding and plotting code with the combined data.
SG=CombineSGdata(mpath,MopStart,MopEnd);

%% identify jumbos by checking the original survey file names for the 
%  word jumbo
utair=find(contains({SG.Source},'UTAir'));
trk=find(contains({SG.Source},'Trk'));
usace=find(contains({SG.Source},'USACE'));
usgs=find(contains({SG.Source},'USGS'));
ccc=find(contains({SG.Source},'CCC'));
utrk=unique([utair trk usace usgs ccc]);

% datestr([SG(jumbo).Datenum])
% fprintf('The SG struct array has %i UTAir Surveys.\n',numel(jumbo))
% 
% fprintf('Plotting the most recent airborne Survey in the database.\n',numel(jumbo))

%% Need to turn the saved SG struct array grid points into a 2d grid array

% make grid area encompassing all survey data
minx=min(vertcat(SG.X));
maxx=max(vertcat(SG.X));
miny=min(vertcat(SG.Y));
maxy=max(vertcat(SG.Y));
%[X,Y]=meshgrid(min(x):max(x),min(y):max(y));
[X,Y]=meshgrid(minx:maxx,miny:maxy);

% inlet location
[e,n,zn]=deg2utm(32.9338,-117.26);

figure('position',[ 21          56        1385         728]);

%SurvNum=trk(end-1);
SurvNum=size(SG,2);
Zmn=X*0+100; % initialize high surface
Zcount=X*0; % count num surveys with grid data at each point
Zmean=X*0;
for SurvNum=1:size(SG,2)
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
Zcount(idx)=Zcount(idx)+1;
Zmean(idx)=Zmean(idx)+Z1(idx);
Zmn=min(Zmn,Z1);
end
Zmn(Zmn == 100)=NaN; % if still high surface, no data
Zmean=Zmean./Zcount;

% now find subaerial volumes
mm=0;
v=[];
nv=[];
for SurvNum=1:size(SG,2)%utrk(2:end)
    mm=mm+1;
% Mop area 1m grid points with valid data
x=SG(SurvNum).X;
y=SG(SurvNum).Y;
% Make 2d x,y utm arrays encompassing the valid points
idx=sub2ind(size(X),y-miny+1,x-minx+1);
Z1=X*NaN; % initialize the 2d elevation Z array as NaNs
Z1(idx)=SG(SurvNum).Z; % overlay the valid data points using the 1d indices
Z1(Z1 < 0.774)=NaN;%Z1(Z1 < 1.344)=NaN;%
Z1(Z1 > 3.0)=NaN;

dz=Z1-Zmn;
v(mm)=sum(dz(:),'omitnan');
if v(mm) < 0
    SurvNum
end
nv(mm)=numel(find(~isnan(dz(:))));
end

%v(nv < 100)=NaN;

plot([SG(1:size(SG,2)).Datenum],v/(100*(MopEnd-MopStart+1)),...
    'k.-','markersize',35);datetick
grid on;
set(gca,'fontsize',18);
xlabel('Date');
ylabel('Mobile Sand Volume (m^{3}/m-shoreline)');
title({'Measured Volume of "Known to be Mobile" Beach Sand',...
    '(Sand Shoreward of MSL)'},...
    'fontsize',24);
box on
set(gca,'linewidth',2)
%makepng('SIOsandVolume.png')