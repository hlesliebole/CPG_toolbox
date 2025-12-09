%function buildNOAAmopDataFiles(NOAAdirPath,SurveySource,MaxElev,MaxDist)
%
% Code to read custom ascii x,y,z *.txt files downloaded via the 
% coast.noaa.gov digital coast data viewer
%
% Custom download "cart" options are lat/lon, elev in meters, points,
%  "ground" and "bathymetry"/"green laser returns", and "last return"
%
%  the resulting downloaded zip file contains one or more Job*.txt files
%
%  this code reads through those txt files and creates a series of
%  x,y,z,t files of the form YYYYMMDD_FirstMop#_LastMop#_SurveySource.llz
%
%  MaxDist is the maximum distance (meters) a survey point can be from a
%  mop point to be saved in a processed data file (eg. some lidar surveys
%  go way inlnad)

MaxElev=5;
MaxDist=200;
% load mop info
load Moptable Mop

NOAAdirPath='/Volumes/LaCie/ShoreBox/USGS*1998';
SurveySource='USGSlidar';

%NOAAdirPath='/Volumes/LaCie/ShoreBox/USGS*2016';
%SurveySource='USGSlidar';

%NOAAdirPath='/Volumes/LaCie/ShoreBox/CCC*20092011';
%SurveySource='CCClidar';

%NOAAdirPath='/Volumes/LaCie/ShoreBox/USACEshoals*2014';
%SurveySource='USACEShoals';


filesd=dir([NOAAdirPath '/Job*.las']); %get list of txt files

% loop thorugh files
for n=1:size(filesd,1)
%     d=lasdata('Job545912_33117_16_35.las','loadall');
%     t=d.gps_time(1)
%     td=(10^9+t)/(60*60*24) + datenum('1980-01-06 00:00');
    %read in the data
    dfile=[filesd(n).folder '/' filesd(n).name];
    fprintf('%g of %g : %s\n',n,size(filesd,1),dfile)
    d=lasdata(dfile,'loadall');
    fprintf('%12i points read in.\n',size(d.z,1))
    
    i=find(d.z <= MaxElev); % make an initial beach/estuary-only elevation cutoff
    if(~isempty(i))
     fprintf('Keeping %12i points below %g m elevation.\n',length(i),MaxElev)
     
        
    x=d.x(i); % lon
    y=d.y(i); % lat
    z=d.z(i); % elev navd88
    t=d.gps_time(i); % gps time
    clear d
    
    % find mop range covered by data
    [xmin,i]=min(x);[xmax,j]=max(x);
    [ymin,k]=min(y);[ymax,l]=max(y);
    % add an edge to include bounding mop points for nearest search
    edge=0.0001*MaxDist;
    xmin=xmin-edge;xmax=xmax+edge;
    ymin=ymin-edge;ymax=ymax+edge;
    
    imp=find((Mop.BackLon >= xmin & Mop.BackLon <= xmax) | ...
        (Mop.BackLat >= ymin & Mop.BackLat <= ymax));
    if(~isempty(imp))
    fprintf('Possible Mop Bounds are: %g to %g\n',min(imp),max(imp))
    
    % find nearest mop back beach point and distance to each lidar point
    % assign survey points to the nearest mop line based on closest
    %  mop backbeach point, use survey  y value with cosine latitude
    %  for consistent approx local scaling
    fprintf('Finding Nearest Mop Points to Survey Points.\n')
    [dp,jmop]=pdist2([Mop.BackLon(imp)*cosd(y(1)),Mop.BackLat(imp)],...
       [x*cosd(y(1)),y],'euclidean','Smallest',1);
    dp=dp*111139; % approx distances from mop points in meters
    imop=imp(jmop)';
       
    % make second data cutoff if a land elevation is more than MaxDist
    %  from a mop point
    i=find(dp' > MaxDist & z > 1);
    if(~isempty(i))
    fprintf('Removing %12i land points farther than %g m from Mop points.\n',...
        length(i),MaxDist)
    x(i)=[];y(i)=[];z(i)=[];t(i)=[];imop(i)=[];dp(i)=[];
    end
    
    if(~isempty(z))
%     
%     z=z-0.778;
%     zscaled = 1+(z-min(z))*10;
%     cn = ceil(max(zscaled));
%     cm = colormap(demcmap(z*10,cn));
%     figure
%     scatter(x,y,5, cm(ceil(zscaled),:));
%     colorbar
%     set(gca,'dataaspectratio',[1 1 1]);
%     set(gca,'clim',[min(z) max(z)])

    % survey point timestamps as datenums
    %  shoals data has unusual time stamp
%    if(strcmp(SurveySource,'USACEShoals'))
%     td=(2*10^9+t)/(60*60*24) + datenum('1980-01-06 00:00');
%    else
     td=(10^9+t)/(60*60*24) + datenum('1980-01-06 00:00');
%    end
    td=floor(td); % truncate to day
    ud=unique(td); % unique survey day datenums
    fprintf('Data is from %g survey days\n',length(ud));
    
    % loop through ranges of mop numbers in 'finc' mop increments
    % and write data to llz files
    finc=50; % number of mops per output file
    
    % open output mop range files
    mops=min(imop);mope=max(imop);
    mopf1=1+finc*(floor((mops-1)/finc)); % starting mop file range
    mopf2=1+finc*(floor((mope-1)/finc)); % ending mop file range
    nfiles=1+floor((mopf2-mopf1)/finc);
    
    % figure out date for filename based on gps time stamp of first
      %  data point to go in the file
     tdn = (10^9+t(1))/(60*60*24) + datenum('1980-01-06 00:00');
     
    % make vector of output file numbers to go with the matched mop data
    ifile=1+floor((imop-mopf1)/finc);
    for k=1:length(ud) % loop thru survey days
        tdn=ud(k); % survey date
     for i=1:nfiles % loop through mop range output files
          % make file name
        strmop1=num2str(mopf1+(i-1)*finc,'%3.3i');
        strmop2=num2str((finc-1)+mopf1+(i-1)*finc,'%3.3i');
        ofile=[datestr(tdn,'yyyymmdd') '_' strmop1 '_' strmop2 '_'...
            SurveySource '.llz'];
        % find data points associated with this file
        j=find(ifile == i & td' == tdn);
        if(~isempty(j))
        fprintf('-----------------\n')
        fprintf('Opening and appending %12i points to\n %s\n',length(j),ofile)
        fprintf('-----------------\n')
         fid=fopen(ofile,'a');
         fprintf(fid,'%f %f %f\n',[y(j)';x(j)';z(j)']);
         fclose(fid);  
        end
     end
    end
    
    else
        fprintf('No data remaining to write to an output file.\n')
    end
    
    else
        fprintf('No data remaining to write to an output file.\n')
    end
    
    else
        fprintf('No data remaining to write to an output file.\n')
    end
    clear x y z t imop dp
 end

%end