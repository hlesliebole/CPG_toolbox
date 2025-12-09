
load MopTableUTM
Mop1=495;%507;
Mop2=517;

%% set lat,lonrange based on mop locations
edge=0.001;
ymin=min([Mop.BackLat(Mop1:Mop2)' Mop.OffLat(Mop1:Mop2)']-edge);
ymax=max([Mop.BackLat(Mop1:Mop2)' Mop.OffLat(Mop1:Mop2)']+edge);
xmin=min([Mop.BackLon(Mop1:Mop2)' Mop.OffLon(Mop1:Mop2)']-edge);
xmax=max([Mop.BackLon(Mop1:Mop2)' Mop.OffLon(Mop1:Mop2)']+edge);

figure('position',[210 403 713 855]);
ax1=axes;set(ax1,'xlim',[xmin xmax],'ylim',[ymin ymax]);hold on;
% plot google map
plot_google_map('MapType', 'satellite','Alpha', 1,'axis',ax1)

hold on;
%SG=SGcombineMops(Mop1,Mop2);
trk=find(strcmp({SG.Source},'Trk'));
shoals=find(strcmp({SG.Source},'USACE'));
n=shoals(end-3:end);
[y,x]=utm2deg(vertcat(SG(n).X),vertcat(SG(n).Y),repmat('11 S',[length(vertcat(SG(n).X)) 1]));
z=vertcat(SG(n).Z);
x(z > 1)=[];
y(z > 1)=[];
z(z > 1)=[];

load BeachColorMap

zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    (BeachColorBarLims(2)-BeachColorBarLims(1));
zscaled(zscaled < 1)=1;
zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);
zscaled(isnan(z))=1;

%cn = ceil(max(zscaled));                                       
cm = BeachColorMap;

scp=scatter(x, y, 10, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;

n=trk(end);
[y,x]=utm2deg(vertcat(SG(n).X),vertcat(SG(n).Y),repmat('11 S',[length(vertcat(SG(n).X)) 1]));
z=vertcat(SG(n).Z);
zscaled = 1+size(BeachColorMap,1)*(z-BeachColorBarLims(1))/...
    (BeachColorBarLims(2)-BeachColorBarLims(1));
zscaled(zscaled < 1)=1;
zscaled(zscaled > size(BeachColorMap,1))=size(BeachColorMap,1);
zscaled(isnan(z))=1;

%cn = ceil(max(zscaled));                                       
cm = BeachColorMap;

scp=scatter(x, y, 10, cm(ceil(zscaled),:), 'filled');
scp.MarkerFaceAlpha = .6;
scp.MarkerEdgeAlpha = .6;

BeachColorbar
