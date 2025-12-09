%  Example code to find the truck tif grid point global minimas for a input
%  Mop area.

MopNumber=552; % Set mop number

%% --------------------------------

load('MopTableUTM.mat','Mop'); % load Mop transect info

%  With reefbreak mounted on a mac, make a structure array listing of all truck
%  beach only tif files. Modify path for windows machine and/or 
%  cliff tif files.
fprintf('Making a list of truck tif files...\n')
TifFiles=dir('/Volumes/group/LiDAR/VMZ2000_Truck/LiDAR_Processed_Level2/*/Beach_Only/2*ground.tif');

% now loop through files names and add the survey file datetime and
%  the starting and ending mop numbers from the file name as new 
%   fields in the TifFiles structure array

for n=1:size(TifFiles,1)
    % parse out date and mop range from the file name
    info=regexprep(TifFiles(n).name,'_', ' '); % separate file name parts with spaces
    info=regexp(info,'\S*','match'); % parse into a cell array
    TifFiles(n).Datetime=datetime(info{1},'inputformat','yyyyMMyy');
    TifFiles(n).MopStart=str2double(info{2});
    TifFiles(n).MopEnd=str2double(info{3});
end

TifFiles % display struct array

%% -----------------------------------

% find tif files that include the mop number
idx=find([TifFiles.MopStart] <= MopNumber & [TifFiles.MopEnd] >= MopNumber);
fprintf('%i tif files found for Mop %i\n',numel(idx),MopNumber)

% make a new struct array Tdata that contains the tif data closest to the
%  mop number for all the truck surveys

Tdata=[];
n=0;
for i=idx
    n=n+1; % survey counter
    
    Tdata(n).Datetime=TifFiles(i).Datetime;
    
    % load tif data
    % suppress matlab warning about undefined geokey in tif file    
    warning('off','map:geotiff:undefinedGTModelTypeGeoKey');
        
    % load x y z data, make double precision
    fprintf('Reading: %s\n',TifFiles(i).name)

    % read in survey data and bound info
    [surv,R]=geotiffread([TifFiles(i).folder '/' TifFiles(i).name]);
    surv=double(surv);
    % extract data array indices with elevations
    [i,j]=find(surv > -9999);
    z=surv(surv > -9999);
    % assign utm coordinates to each elevation
    % (there may be a small offset error here, have
    %  Adam check if this is correct conversion)
    xutm=R.XWorldLimits(1)+j;
    yutm=R.YWorldLimits(2)-i;
    
    % find points closest to the mop number
    N=XY2MopNumsV2(xutm,yutm,Mop);
    
    ndx=find(N == MopNumber);
    if ~isempty(ndx)
    % save the xutm, yutm, z values in the struct array
     Tdata(n).X=xutm(ndx);
     Tdata(n).Y=yutm(ndx);
     Tdata(n).Z=z(ndx);
    else
     Tdata(n).X=[];
     Tdata(n).Y=[];
     Tdata(n).Z=[];
    end
        
end

% define a utm universal grid area that encapsulates all the 
%   tif points for the mop area
minx=min(vertcat(Tdata.X));
maxx=max(vertcat(Tdata.X));
miny=min(vertcat(Tdata.Y));
maxy=max(vertcat(Tdata.Y));
% 2d universal grid indices
[X,Y]=meshgrid(minx:maxx,miny:maxy);

% make a 2d matrix Zt(j,k) of all the tif grid survey points,
%  where j is the survey number, and k is the 1d indices of
%  the 2d universal grid points

j=0; % survey counter
Zt=NaN(size(Tdata,2),size(X(:),1)); % initialize Zt as no data NaNs

    for n=1:size(Tdata,2)
      j=j+1;
      x=Tdata(n).X;
      y=Tdata(n).Y;
      % find 1d indices of valid x,y tif grid points for this survey
      ndx=sub2ind(size(X),y-miny+1,x-minx+1);
      % add then to the Zt matrix
      Zt(j,ndx)=Tdata(n).Z;
    end

% Count the number of valid surveys at each grid point
Zcount=X*NaN;Zcount(:)=sum(~isnan(Zt)); % sum of valid grid point values
Zcount(Zcount == 0)=NaN; % grid points with no values set to no data
figure('position',[101   261   701   536]);
surf(X,Y,Zcount);colormap(jet);colorbar;view(2)
xlabel('E UTM');ylabel('N UTM');
title(['Mop ' num2str(MopNumber) ' Number of Values at each tif grid point']);
    
% Get the min z at each tif grid point that has at least
%  one valid survey point.
Zmin=X*NaN;Zmin(:)=min(Zt);

figure('position',[176   362   859   406]);
surf(Zmin);shading flat;colormap(jet);colorbar
view(-20,20);xlabel('E UTM');ylabel('N UTM');zlabel('Elev (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ' Minimum Tif Elevations']);

% Get the max z at each tif grid point that has at least
%  one valid survey point.
Zmax=X*NaN;Zmax(:)=max(Zt);

figure('position',[257   306   859   406]);
surf(Zmax);shading flat;colormap(jet);colorbar
view(-20,20);xlabel('E UTM');ylabel('N UTM');zlabel('Elev (m, NAVD88)');
title(['Mop ' num2str(MopNumber) ' Maximum Tif Elevations']);

% use the matlab isoutlier function to flag suspect elevations at
% each tif grid point
TF=isoutlier(Zt);

% plot all of the tif elevations as a grey point cloud
figure('position',[347    61   992   573]);
hold on;
for n=1:size(TF,1)
    plot3(X(TF(n,:)),Y(TF(n,:)),Zt(n,TF(n,:)),'m.','markersize',15)
end
scp=scatter3(vertcat(Tdata.X),vertcat(Tdata.Y),vertcat(Tdata.Z),'.',...
    'MarkerFaceColor',[.8 .8 .8],'MarkerEdgeColor',[.8 .8 .8]);
scp.MarkerFaceAlpha = .1;
scp.MarkerEdgeAlpha = .1;
scp.SizeData=1;
grid on;
view(0.01,0)
title(['Mop ' num2str(MopNumber) ' isoutliers'])

%------------------------------------------------------------
function  N=XY2MopNumsV2(xutm,yutm,Mop)
%------------------------------------------------------------

sflag=0;
if(numel(xutm) == 1);xutm=[xutm;xutm];yutm=[yutm;yutm];sflag=1;end
% Returns a vector of the nearest mop transect numbers, N, based on the mop transect
% definition table, Mop, for the input x,y data point vectors.

%------------------------------------------------------------
%fprintf('Making approximate Mop match to start...\n')
% To narrow down the possible range of Mops in the survey area but avoid 
% excluding Mops at theedges of the survey region, set a distance tolerance, 
% tol, in meters around the min-max box of survey points.
tol=1000; % use 1000m to be safe
minx=min(xutm)-tol;maxx=max(xutm)+tol;
miny=min(yutm)-tol;maxy=max(yutm)+tol;
% find the Mops whose back beach or offshore point falls within
% the survey min-max + tolerance box.
mopidx=find( ((Mop.BackXutm >= minx & Mop.BackXutm <= maxx) | ...
        (Mop.OffXutm >= minx &  Mop.OffXutm <= maxx)) & ...
        ((Mop.BackYutm >= miny & Mop.BackYutm <= maxy) | ...
        (Mop.OffYutm >= miny &  Mop.OffYutm <= maxy)));
% reduce to possible start and end Mop numbers
MopStart=min(mopidx);MopEnd=max(mopidx);
% As first approximate matching of points to Mop transects, find
% closest of thse Mop transect back beach points (not the Mop transect line) 
% to each survey data point.
mtx=[];
mty=[];
mnum=[];
mtx=Mop.BackXutm(MopStart:MopEnd)'; 
mty=Mop.BackYutm(MopStart:MopEnd)'; 
mopnum=MopStart:MopEnd;
% match survey utm points to the nearest back beach points
[dp,NearIdx]=pdist2([mty',mtx'],[double(yutm),double(xutm)],'euclidean','smallest',1);
% find min and max of actual matched mop numbers
ApproxMop=mopnum(NearIdx);
%fprintf('Mop Match Range: %i to %i\n',min(ApproxMop),max(ApproxMop))
%  Get max distance of the data points from their nearest back beach point
%  Use this when defining Mop transect lines as series of closely spaced
%  points
MaxDist=max(ceil(dp)); 
%fprintf('Max Distance from Back Beach Points: %i\n',MaxDist)
%-----------------------------------------------------------------
% Now loop through approximately matched Mop numbers and do a more
% exact matching of data points to the Mop transect lines

%fprintf('Making more precise data match to each Mop transect line...\n')

N=ApproxMop;  % initialize N as the approximate Mop numbers

% Step up the coast for each Mop number. Get the survey points with appox Mop
%  numbers equal to this Mop and the next upcoast Mop and find which data
%  points are actually closest to the Mop's transect line rather than just
%  the back beach point.

for mn=min(ApproxMop)-1:max(ApproxMop)+1
    
    % mn is the target Mop number for exact nearest data matching
    if mn == min(ApproxMop)-1
%     ll=fprintf('Finding points closest to Mop %i of %i ',mn,max(ApproxMop)+1);
    else
%     fprintf(repmat('\b',1,ll))
%     ll=fprintf('Finding points closest to Mop %i of %i ',mn,max(ApproxMop)+1);
    end
    
    if mn > 0 && mn < 11594 % can't go below Mop #1
        
    % get indexes of all nearby data points with approx nearest mop number
    % equal to mn or the next upcoast mop.   
    idx=find(N >= mn & N <= mn+1); 
    
    pt=[xutm(idx) yutm(idx) xutm(idx)*0]; % make 3d [x y 0] nearby point array 

    % find distance of points to Mop mn transect line
    Mopnum=mn;
    v1=[Mop.BackXutm(Mopnum),Mop.BackYutm(Mopnum),0];
    v2=[Mop.OffXutm(Mopnum),Mop.OffYutm(Mopnum),0];
    d1 = point_to_line(pt, v1, v2); % distance to mn line
    
    % find distance of points to Mop mn+1 transect line
    Mopnum=mn+1;  
    v1=[Mop.BackXutm(Mopnum),Mop.BackYutm(Mopnum),0];
    v2=[Mop.OffXutm(Mopnum),Mop.OffYutm(Mopnum),0];
    d2 = point_to_line(pt, v1, v2); % distance to mn+1 line
    
    % Correct the overall nearest mop number vector N
    N(idx(d1 <= d2))=mn; % points closest to target mop
    N(idx(d1 > d2))=mn+1; % points closer to target mop+1
    end     
end

if sflag ==1;N=N(1);end
fprintf('Nearest Mop range: %i %i\n',min(N),max(N))
%fprintf('\n');
%fprintf('\n%i Data Pts Corrected from approx to actual nearest Mop transect.\n',length(find((ApproxMop - N) > 0)))
end

%------------------------------------------------------------
function d = point_to_line(pt, v1, v2)
%------------------------------------------------------------

% pt should be nx3
% v1 and v2 are vertices on the line (each 1x3)
% d is a nx1 vector with the orthogonal distances
% eg.
% v1 = [0,0,0];
% v2 = [3,0,0];
% pt = [0,5,0];
% distance = point_to_line(pt,v1,v2)


v1 = repmat(v1,size(pt,1),1);
v2 = repmat(v2,size(pt,1),1);
a = v1 - v2;
b = pt - v2;
d = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
end