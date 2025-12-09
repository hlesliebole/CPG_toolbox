function [xutm,yutm,z,c,utmzone]=readSurveyFileUTM2(dataSource,dataFilename)

xutm=[];yutm=[];z=[];c=[];utmzone=[];


if exist(dataFilename,'file')
%-----------------------------------------
% GPS survey data loading
%-----------------------------------------

if strcmpi(dataSource,'gps') || strcmpi(dataSource,'iG8wheel') ...
        || strcmpi(dataSource,'OCSwheel') % GPS llz file

% load x y z data, make double precision
fprintf('Reading: %s\n',dataFilename)
llz=load(dataFilename);
%  GPS survey *.navd88 data files can have 3,4,5 or 7 columns
%   depending on how old they are.
%
% The oldest 3 col files do not have sample date/time
% check if data file has a time column or not

    y=double(llz(:,1)); % Lat as y
    x=double(llz(:,2)); % Lon as x
    if(size(llz,2) == 3) % original 3 column files of lat, lon, z
     z=double(llz(:,3)); % elev m, navd88 
     c=z*0; % unknown substrate classification
     %t=zeros(size(z));
    elseif(size(llz,2) == 4) % original 4 column files of lat, lon, z, time
     z=double(llz(:,3)); % elev m, navd88
     c=z*0; % unknown substrate classification
     %t=double(llz(:,4)); % time stamp
    elseif(size(llz,2) == 5) % 5 column files with substrate code in col 5
     z=double(llz(:,3));
     c=double(llz(:,5)); % GPS substrate classification
     %t=double(llz(:,4));
    elseif(size(llz,2) == 6) % 5 column files with substrate code in col 5
     z=double(llz(:,5));
     c=z*0; % unknown substrate classification
     %t=double(llz(:,6));
    else % 7 column files with UTM coords in cols 3 4, and substrate in 7
     z=double(llz(:,5));
     c=double(llz(:,7)); % GPS substrate classification
     %t=double(llz(:,6));    
    end

    % change any gps substrate code > 3 (eg. jetski =7) to 0
    %  0=unknown, 1=sand; 2=cobble; 3=bedrock
    c(c > 3)=0;

    % convert to UTM
    [xutm,yutm,utmzone]=deg2utm(y,x);
    
    
%-----------------------------------------
% Lidar survey data loading
%-----------------------------------------    
   
else % otherwise, source is some form of lidar data
    
    ext=dataFilename(end-2:end); % get filename extension

  %-----------------------------------------       
  % Truck/ATV MiniRanger/Drone MiniRanger tif files
  %-----------------------------------------   

    if(strcmpi(ext,'tif')) % use tifread if tif file
    
    % suppress matlab warning about undefined geokey in tif file    
    warning('off','map:geotiff:undefinedGTModelTypeGeoKey');
        
    % load x y z data, make double precision
    fprintf('Reading: %s\n',dataFilename)

    % read in survey data and bound info
    [surv,R]=geotiffread(dataFilename);
    surv=double(surv);
    % extract data array indices with elevations
    [i,j]=find(surv > -9999);
    z=surv(surv > -9999);
    % tif has utm coordinates
   if strcmp(R.CoordinateSystemType,'planar')
    % assign utm coordinates to each elevation
    % (there may be a small offset error here, have
    %  Adam check if this is correct conversion)
    xutm=R.XWorldLimits(1)+j;
    yutm=R.YWorldLimits(2)-i;
    utmzone='11 S';
    % convert utm coords to lat lons
    %[y,x]=utm2deg(Eutm,Nutm,repmat('11 S',[length(Eutm) 1]));
    c=z*0; % unknown substrate classification
   else
    % tif has lat lon coordinates 
    xlon=R.LongitudeLimits(1)+j*R.CellExtentInLongitude;
    ylat=R.LatitudeLimits(2)-i*R.CellExtentInLatitude;
    
    c=z*0; % unknown substrate classification
    % convert to UTM
    [xutm,yutm,utmzone]=deg2utm(ylat,xlon);
    utmzone=utmzone(1,:);
   end
    

%-----------------------------------------       
% RTK Drone las files
%-----------------------------------------   


    elseif(strcmpi(ext,'las')) % las file eg. RTK Drone
    fprintf('Reading: %s\n',dataFilename)
    d=lasdata(dataFilename,'loadall');
    xutm=d.x;
    yutm=d.y;
    z=d.z;
    c=z*0; % unknown substrate classification
    utmzone='11 S';

  %-----------------------------------------       
  % UTexas and 97/98 el nino airborne lidar llz files
  %-----------------------------------------   

    else % otherwise load columns of x,y,z lidar data (txt, or lidar llz file)
        if(strcmpi(dataSource,'utair')) % UT airborne lidar saved as lon, lat, z
        
        % load x y z data, make double precision
        fprintf('Reading: %s\n',dataFilename)
        llz=load(dataFilename);
        
        x=double(llz(:,1)); % Lon as x
        y=double(llz(:,2)); % Lat as y
        z=double(llz(:,3));
        c=z*0; % unknown substrate classification
        
        % convert to UTM
        [xutm,yutm,utmzone]=deg2utm(y,x);
       
  %-----------------------------------------       
  % Other airborne lidar llz files processed from NOAA coast site
  %-----------------------------------------   
     
      
        else % other airborne lidar saved at lat, lon, z 
            
        % load x y z data, make double precision
        fprintf('Reading: %s\n',dataFilename)
        llz=load(dataFilename);
        
        x=double(llz(:,2)); % Lon as x
        y=double(llz(:,1)); % Lat as y
        z=double(llz(:,3));
        c=z*0; % unknown substrate classification
        
        % convert to UTM
        [xutm,yutm,utmzone]=deg2utm(y,x);
        
        end
    end
end


% remove duplcate points before returning data


%----------------------------
%  Now do standard QC of duplicate x,y,z points
%------------------------------

% qc for duplicate survey locations (same lat,lon) and duplicate
%  survey data (same lat, lon, z)
npts=size(xutm,1);
% check for duplicate x, y, z
[uxy,ia,ic]=unique([xutm, yutm, z],'rows');uniqxyz=size(uxy,1);
dupxyz=npts-uniqxyz;
% remove those points
xutm=xutm(ia);yutm=yutm(ia);z=z(ia);c=c(ia);%t=t(ia);
nxyzpts=size(xutm,1);

% now check for duplicate x, y with different z
[uxy,ia,ic]=unique([xutm, yutm],'rows');uniqxy=size(uxy,1);
dupxy=nxyzpts-uniqxy;
% remove those points too
xutm=xutm(ia);yutm=yutm(ia);z=z(ia);c=c(ia);%t=t(ia);
nxyzpts=size(xutm,1);

fprintf(1,'%s\n','Num_Points  Unique_Pts  Dup_Points  Dup_Point_diff_z')
fprintf(1,'%8d %8d %10d %10d\n',npts,nxyzpts,dupxyz,dupxy)

end

end

