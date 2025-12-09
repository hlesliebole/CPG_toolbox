
clearvars
% s=webread('http://coastsat.wrl.unsw.edu.au/time-series/usa_CA_0011-0015/');
% T=readgeotable('/Users/William/desktop/CoastSat_transect_layer.geojson');
% load Mop point info table
load MopTableUTM.mat

tislands=[23 24 29 34 36 41 46 48 50 79 80 81 83 86 87 91];
%g=readgeota('CoastSat_transect_layer.geojson');
s=fileread('CoastSat_transect_layer.geojson');
ss=strsplit(s,'TransectId');
n=0;
for nn=2:size(ss,2)
    st=strsplit(ss{nn});
    VosName=regexprep(st{2}, {'"',','}, '');
    if contains(VosName,'usa_CA')
        %fprintf('%s\n',VosName)
        n=n+1;
    VosTCA(n).Name=regexprep(st{2}, {'"',','}, '');
    VosTCA(n).Orientation=str2double(regexprep(st{6}, {'"',','}, ''));
    VosTCA(n).Slope=str2double(regexprep(st{8}, {'"',','}, ''));
    VosTCA(n).Trend=str2double(regexprep(st{10}, {'"','},'}, ''));
    VosTCA(n).BackLon=str2double(regexprep(st{15}, {'[[',','}, ''));
    VosTCA(n).BackLat=str2double(regexprep(st{16}, {']',','}, ''));
    VosTCA(n).OffLon=str2double(regexprep(st{17}, {'[',','}, ''));
    VosTCA(n).OffLat=str2double(regexprep(st{18}, {']]','}},'}, ''));
    [XutmV,YutmV,UTMzone]=deg2utm(VosTCA(n).BackLat,VosTCA(n).BackLon);
     VosTCA(n).BackXutm=XutmV;
     VosTCA(n).BackYutm=YutmV; 
     [XutmV,YutmV,UTMzone]=deg2utm(VosTCA(n).OffLat,VosTCA(n).OffLon);
     VosTCA(n).OffXutm=XutmV;
     VosTCA(n).OffYutm=YutmV; 
     VosTCA(n).UTMzone=UTMzone;
     orient=90-atan2d(VosTCA(n).OffYutm-VosTCA(n).BackYutm,...
                      VosTCA(n).OffXutm-VosTCA(n).BackXutm);
     if orient < 0;orient=orient+360;end
     VosTCA(n).OrientationUTM=orient;
     VosTCA(n).BackMop=NaN; % default is no nearest Mop

    % if contains(VosTCA(n).Name,'usa_CA')
    %     tnum=str2num(VosTCA(n).Name(8:11));
    %     %fprintf('%s %d\n',VosT(n).Name,tnum)
    % else
    %     tnum=0;
    % end
    % 
    % if contains(VosTCA(n).Name,'usa_CA') & ~ismember(tnum,tislands)
    % 
    %     %[MopFrac,X]=LatLon2MopFractionXshoreX(VosT(n).BackLat,VosT(n).BackLon);
    % [XutmV,YutmV,UTMzone]=deg2utm(VosTCA(n).BackLat,VosTCA(n).BackLon);
    %  VosTCA(n).BackXutm=XutmV;
    %  VosTCA(n).BackYutm=YutmV; 
    %  VosTCA(n).UTMzone=UTMzone; 
    % 
    % % find nearest mop to the Vos transect back beach point
    % [NmopV,dpV]=FindNearestMopTransectsUTM(XutmV,YutmV);
    % % shift Vos back beach location on Vos transect line by
    % %  the distance dp to be one or seaward of the mop back beach
    % %  line
    % Xutm=XutmV+dpV*cosd(90-VosTCA(n).Orientation);
    % Yutm=YutmV+dpV*sind(90-VosTCA(n).Orientation);
    % 
    % % find the nearest Mop again
    % [Nmop,dp]=FindNearestMopTransectsUTM(Xutm,Yutm);
    % if Nmop ~= NmopV
    %     fprintf('*** Different Mop with Vos seaward shift: %i %i %d %d\n',NmopV,Nmop,dpV,dp)
    % end

    % if Nmop == 1
    %     MopFrac=1;
    %     X=NaN;
    % elseif Nmop == 11594
    %     MopFrac=11594;
    %     X=NaN;
    
    %else 
    
    % % divide mop area into 101 mop subtransects (~1m apart alongshore) 
    % % with 1m xshore resolution, with an extra 100m of back beach for each
    % % in case the point falls landward of the back beach line
    % Ntransects=101;
    % [x1d,xt,yt,xst,yst]=GetTransectLines(Mop,Nmop,Ntransects,[-1000 0]);
    % 
    % % find nearest subtransect line point to the input location
    % [dp,NearIdx]=...
    %     pdist2([yst(:),xst(:)],[double(Yutm),double(Xutm)],'euclidean','smallest',1);
    % 
    % % define X based on the nearest transect line point
    % %   row=nearest subtransect number; col = xshore distance indice on
    % %   the nearest subtransect
    %   [row,col] = ind2sub(size(xst),NearIdx); 
    % 
    % % fractional mop number based on nearest subtransect number
    % MopFrac=Nmop+(row-(Ntransects+1)/2)/(Ntransects-1);
    % % xshore distance along the subtransect
    % X=x1d(col);
    % end
% 
%         VosTCA(n).BackMop=Nmop;
%         %VosT(n).BackMopFrac=MopFrac;
%         %VosT(n).BackMopXoffset=X;
%     else
%         VosTCA(n).BackXutm=NaN;
%         VosTCA(n).BackYutm=NaN; 
%         VosTCA(n).BackMop=NaN;
%         VosTCA(n).UTMzone=NaN;
%         %VosT(n).BackMopFrac=NaN;
%         %VosT(n).BackMopXoffset=NaN;
end

end

% ibad=find(isnan([VosTCA.BackMop]));
% VosTCA(ibad)=[];



% now loop through Mops and find closest Vos transect

load('MopTableUTM.mat','Mop');  % Load "Mop" table array
% 
for n=1:size(Mop,1)
%     Xutm=Mop.BackXutm(n);
%     Yutm=Mop.BackYutm(n);
[Xutm,Yutm,UTMzone]=deg2utm(Mop.BackLat(n),Mop.BackLon(n));
% find indices of all CoastSat transects in same UTM zone to
%  narrow the nearest search
idx=find(strcmp({VosTCA.UTMzone},UTMzone)); 
% use pdist function with the CoastSat transect midpoints to find 
%  the closest midpoint to the Mop back beach point
    [dp,lmop]=pdist2([Yutm,Xutm],...
        [([VosTCA(idx).BackYutm]'+[VosTCA(idx).OffYutm]')/2,...
        ([VosTCA(idx).BackXutm]'+[VosTCA(idx).OffXutm]')/2],...
          'euclidean','smallest',1);

   [dpV,imin]=min(dp); % closest CoastSat midpoint is dPV meters away 
   vdx=idx(imin); % nearest CoastSat transect struct array indice is vdx

% if the distance is greater than 100m, there is no close 
%  Vos transect
if dpV < 150
    VosTCA(idx(imin)).BackMop=n;
    fprintf('Mop %i Vos %s\n',n,VosTCA(idx(imin)).Name)
end

end


    % Xutm=XutmV+dpV*cosd(90-VosTCA(n).Orientation);
    % Yutm=YutmV+dpV*sind(90-VosTCA(n).Orientation);
% end
save VosTCA.mat VosTCA

function [Nmop,dp]=FindNearestMopTransectsUTM(Xutm,Yutm)

% Used by BuildSGmatfiles.m

% Assigns Xutm (easting)  Yutm (northing) points to the nearest
%  Mop transect line number Nmop

load('MopTableUTM.mat','Mop');  % Load "Mop" table array

% Find range of possible mops in play by find mop back beach
% points closest to data bounding box corners.
Xmin=min(Xutm);Xmax=max(Xutm);Ymin=min(Yutm);Ymax=max(Yutm);
Xbound=[Xmin Xmax Xmax Xmin];Ybound=[Ymin Ymin Ymax Ymax];   
[dp,lmop]=pdist2([Mop.BackYutm,Mop.BackXutm],[Ybound',Xbound'],...
    'euclidean','smallest',1);
NearMops=min(lmop):max(lmop);

%---------------------- 
%  make interpolated transect line points with approx 1m spatial resolution
mtx=[];
mty=[];
mnum=[];
extend=100; % extend transect off and onshore, in meters


for MopNum=NearMops(1):NearMops(end)
 % figure out x step size for approx 1m alongshore transect resolution
 xstep=min([abs(Mop.BackXutm(MopNum)-Mop.OffXutm(MopNum)) abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum))])/...
        abs(Mop.BackYutm(MopNum)-Mop.OffYutm(MopNum));

 xi=min([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)])-extend:xstep:max([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)]+extend);
 if Mop.BackXutm(MopNum) == Mop.OffXutm(MopNum) && Mop.BackYutm(MopNum) == Mop.OffYutm(MopNum)
     yi=0;
 else
 yi=interp1([Mop.BackXutm(MopNum) Mop.OffXutm(MopNum)],[Mop.BackYutm(MopNum) Mop.OffYutm(MopNum)],...
    xi,'linear','extrap');
 end

mtx=[mtx xi]; 
mty=[mty yi];
mnum=[mnum MopNum*ones(size(xi))];
end

%---------------------------------------------------------------------

 % match survey utm points to the nearest mop transect line

[dp,lmop]=pdist2([mty',mtx'],[Yutm,Xutm],...
    'euclidean','smallest',1);
Nmop=mnum(lmop);
if dp(1) > 400
    fprintf('Warning: Nearest Mop is Far away %d\n',dp(1))
    %pause
end

end


