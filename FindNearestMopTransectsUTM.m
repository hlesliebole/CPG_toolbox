function Nmop=FindNearestMopTransectsUTM(Xutm,Yutm)

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
% size(mtx)
% size(mty)
% size(Yutm)
% size(Xutm)
% size([mty',mtx'])
if size([Yutm',Xutm'],2) ~= 2
    Yutm=Yutm';Xutm=Xutm';
end


[dp,lmop]=pdist2([mty',mtx'],[Yutm',Xutm'],...
    'euclidean','smallest',1);
Nmop=mnum(lmop);

end

