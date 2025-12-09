
clearvars
% s=webread('http://coastsat.wrl.unsw.edu.au/time-series/usa_CA_0011-0015/');
% T=readgeotable('/Users/William/desktop/CoastSat_transect_layer.geojson');
% load Mop point info table
load MopTableUTM.mat

tislands=[23 24 29 34 36 41 46 48 50 79 80 81 83 86 87 91];
%g=readgeota('CoastSat_transect_layer.geojson');
s=fileread('CoastSat_transect_layer.geojson');
ss=strsplit(s,'TransectId');
for n=2:size(ss,2)
    st=strsplit(ss{n});
    VosT(n).Name=regexprep(st{2}, {'"',','}, '');
    VosT(n).Orientation=str2double(regexprep(st{6}, {'"',','}, ''));
    VosT(n).Slope=str2double(regexprep(st{8}, {'"',','}, ''));
    VosT(n).Trend=str2double(regexprep(st{10}, {'"','},'}, ''));
    VosT(n).BackLon=str2double(regexprep(st{15}, {'[[',','}, ''));
    VosT(n).BackLat=str2double(regexprep(st{16}, {']',','}, ''));
    VosT(n).OffLon=str2double(regexprep(st{17}, {'[',','}, ''));
    VosT(n).OffLat=str2double(regexprep(st{18}, {']]','}},'}, ''));

    if contains(VosT(n).Name,'usa_CA')
        tnum=str2num(VosT(n).Name(8:11));
        %fprintf('%s %d\n',VosT(n).Name,tnum)
    else
        tnum=0;
    end

    if contains(VosT(n).Name,'usa_CA') & ~ismember(tnum,tislands)
        
        %[MopFrac,X]=LatLon2MopFractionXshoreX(VosT(n).BackLat,VosT(n).BackLon);
    [XutmV,YutmV,UTMzone]=deg2utm(VosT(n).BackLat,VosT(n).BackLon);
     VosT(n).BackXutm=XutmV;
     VosT(n).BackYutm=YutmV; 
     VosT(n).UTMzone=UTMzone; 
    
    % find nearest mop to the Vos transect back beach point
    [NmopV,dpV]=FindNearestMopTransectsUTM(XutmV,YutmV);
    % shift Vos back beach location on Vos transect line by
    %  the distance dp to be one or seaward of the mop back beach
    %  line
    Xutm=XutmV+dpV*cosd(90-VosT(n).Orientation);
    Yutm=YutmV+dpV*sind(90-VosT(n).Orientation);
    
    % find the nearest Mop again
    [Nmop,dp]=FindNearestMopTransectsUTM(Xutm,Yutm);
    if Nmop ~= NmopV
        fprintf('*** Different Mop with Vos seaward shift: %i %i %d %d\n',NmopV,Nmop,dpV,dp)
    end

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

        VosT(n).BackMop=Nmop;
        %VosT(n).BackMopFrac=MopFrac;
        %VosT(n).BackMopXoffset=X;
    else
        VosT(n).BackXutm=NaN;
        VosT(n).BackYutm=NaN; 
        VosT(n).BackMop=NaN;
        VosT(n).UTMzone=NaN;
        %VosT(n).BackMopFrac=NaN;
        %VosT(n).BackMopXoffset=NaN;
end

end

ibad=find(isnan([VosT.BackMop]));
VosT(ibad)=[];

save VosTransects4.mat VosT

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


