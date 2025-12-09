%a1 = [];
fid=fopen('20240117_waypointtracks.txt','r');
eutm=[];
nutm=[];
while feof(fid) == 0
    t=textscan(fid,'%s %f %f %f','Delimiter',',','HeaderLines',3);
    %a2 = textscan(fidi,format,'delimiter',',','HeaderLines',1,'CollectOutput',1);
    %a1 = [a1;a2{1}];
    %fprinf('%i %i\n',numel(t{2}).numel(t{3}))
    if numel(t{2}) == numel(t{3})
     eutm=[eutm;t{2}];
     nutm=[nutm;t{3}];
    end
end
fclose(fid);
[mLat,mLon]=utm2deg(eutm,nutm,repmat('11 S',[length(eutm) 1]));
hold on;
%plot(mLon,mLat,'g.','markersize',1)
% p=0;
% for n=1:5000:20*5000
% p=p+1;
% text(mLon(n),mLat(n),num2str(p),'color','k')
% end

mLat=mLat(1:15:end);
mLon=mLon(1:15:end);

idx=find(diff(mLat) < 0);
plot(mLon(idx),mLat(idx),'k.','markersize',5)

idx=find(diff(mLat) > 0);
%plot(mLon(idx),mLat(idx),'c.','markersize',1)
%plot(mLon(idx),mLat(idx),'m.','markersize',5)
%%

tx=text(-117.285,32.9813,'17 Jan 2024 Tracklines: Solid Black = Southbound ; Dotted Black= Northbound','color','y','fontsize',16);
makejpg('MultibeamSandWavesTracklines.jpg')